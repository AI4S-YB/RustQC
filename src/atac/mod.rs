//! ATAC-Seq quality control and Tn5 preprocessing.
//!
//! Implements bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, and library
//! complexity, plus optional +4/-5 Tn5 shift and fixed-interval NFR/mono/di/tri
//! BAM split. Numerical fidelity targets ATACseqQC 1.36.0.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::bam_io;
use crate::cli::AtacArgs;
use crate::config::AtacConfig;

use bam_qc::{BamQcAccum, PbcChromAccum};
use frag_size::FragSizeAccum;
use lib_complexity::DupFreqAccum;
use tss_cov::TssCov;

/// Small helper for the mate-pair buffer.
struct FirstMateInfo {
    pos1: i64,
    tlen: i64,
}

/// Helper: median of a Vec<f64>.
fn median(mut xs: Vec<f64>) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = xs.len() / 2;
    if xs.len() % 2 == 0 {
        (xs[mid - 1] + xs[mid]) / 2.0
    } else {
        xs[mid]
    }
}

/// Resolve the output path for a metric file, either flat or in a subdirectory.
fn metric_path(outdir: &str, flat: bool, subdir: &str, filename: &str) -> PathBuf {
    if flat {
        Path::new(outdir).join(filename)
    } else {
        Path::new(outdir).join(subdir).join(filename)
    }
}

/// Entry point for the `rustqc atac` subcommand.
///
/// Loads the merged YAML config (XDG system → user → env → -c flag), extracts
/// the `atac:` section, and merges it with CLI flags via [`resolve`].
pub fn run(args: AtacArgs) -> Result<()> {
    // Load merged config from all sources (same pattern as run_rna in main.rs).
    let (full_cfg, _config_sources) =
        crate::config::load_merged_config(args.config.as_deref())?;
    let atac_cfg = full_cfg.atac;

    let cfg = resolve(&args, &atac_cfg);

    // Validate: Phase 13 supports exactly one BAM input.
    if cfg.inputs.len() != 1 {
        anyhow::bail!(
            "rustqc atac expects exactly one BAM input, got {}",
            cfg.inputs.len()
        );
    }
    let input = &cfg.inputs[0];

    // PE-check the BAM.
    pe_check::assert_paired_end(Path::new(input))
        .with_context(|| format!("paired-end check failed for {}", input))?;

    // Derive sample name.
    let sample = cfg.sample_name.clone().unwrap_or_else(|| {
        Path::new(input)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string()
    });

    eprintln!("[rustqc atac] sample: {}", sample);

    // Load TSS list from GTF.
    let tss_list = crate::gtf::extract_tss(Path::new(&cfg.gtf))
        .with_context(|| format!("failed to parse GTF: {}", cfg.gtf))?;
    if tss_list.is_empty() {
        eprintln!("[rustqc atac] WARNING: no TSS entries extracted from GTF — TSS metrics will be empty");
    } else {
        eprintln!("[rustqc atac] loaded {} TSS entries from GTF", tss_list.len());
    }

    // Resolve flank.
    let flank = resolve_flank(cfg.tsse_flank);
    let mut tss_cov = TssCov::new(tss_list, flank);

    // Open BAM.
    let (mut reader, header) = bam_io::open(Path::new(input))
        .with_context(|| format!("failed to open BAM: {}", input))?;

    // Get chromosome names.
    let seq_names: Vec<String> = bam_io::reference_sequences(&header)
        .into_iter()
        .map(|(n, _)| n)
        .collect();

    // Detect mitochondrial chromosome.
    let mito = cfg.mito_chrom.clone().or_else(|| {
        mito::detect_mito(&seq_names).map(String::from)
    });
    if let Some(ref m) = mito {
        eprintln!("[rustqc atac] mito chromosome: {}", m);
    }

    // Allocate per-chromosome PBC accumulators.
    let n_chroms = seq_names.len();
    let mut pbc: Vec<PbcChromAccum> = (0..n_chroms).map(|_| PbcChromAccum::default()).collect();

    // Initialize other accumulators.
    let mut bq = BamQcAccum::new();
    let mut frag = FragSizeAccum::new();
    let mut dup = DupFreqAccum::default();

    // Mate buffer for PBC pair reconstruction: qname → first-seen mate info.
    let mut mate_buffer: HashMap<Vec<u8>, FirstMateInfo> = HashMap::new();

    // ── Single-pass BAM scan ──────────────────────────────────────────────────
    let mut n_records: u64 = 0;
    for record_result in reader.records() {
        let record = record_result.context("failed to read BAM record")?;
        let flags = u16::from(record.flags());

        // Skip secondary / supplementary.
        if flags & (0x100 | 0x800) != 0 {
            continue;
        }

        let tid_i = bam_io::tid(&record);
        if tid_i < 0 {
            // Unmapped reads with no tid: still feed into bamQC (unmapped flag set).
            let q = bam_io::qname(&record);
            let qname_str = std::str::from_utf8(q).unwrap_or("").to_string();
            let mapq = bam_io::mapq(&record);
            bq.update_flags(flags, mapq, false, &qname_str);
            n_records += 1;
            continue;
        }

        let tid = tid_i as usize;
        let chrom_name = &seq_names[tid];
        let is_mito = mito.as_deref().is_some_and(|m| chrom_name == m);

        let q = bam_io::qname(&record);
        let qname_str = std::str::from_utf8(q).unwrap_or("").to_string();
        let mapq = bam_io::mapq(&record);
        let tlen = record.template_length() as i64;

        // bamQC flag/MAPQ counters.
        bq.update_flags(flags, mapq, is_mito, &qname_str);
        n_records += 1;

        // fragSize: skip unmapped (0x4) and QC-fail (0x200).
        if flags & (0x4 | 0x200) == 0 {
            frag.update(tlen);
        }

        // 5'-end position for TssCov.
        let pos0 = bam_io::pos_0based(&record);
        if pos0 >= 0 && (flags & 0x4) == 0 {
            let is_reverse = flags & 0x10 != 0;
            let pos5p = if is_reverse {
                let endp = bam_io::end_pos_0based(&record);
                if endp >= 0 {
                    endp as u64
                } else {
                    (pos0 + 1) as u64
                }
            } else {
                (pos0 + 1) as u64
            };
            tss_cov.add_5prime(chrom_name, pos5p);
        }

        // PBC fingerprint: only for proper pairs where both mates are mapped and on same chrom.
        if (flags & 0x1) != 0 && (flags & 0x4) == 0 && (flags & 0x8) == 0 {
            let is_first = flags & 0x40 != 0;
            let mate_tid = bam_io::mtid(&record);
                if mate_tid as usize == tid {
                // Same-chromosome pairs only.
                if let Some(prev) = mate_buffer.remove(q) {
                    let (p1, t1, p2, t2) = if is_first {
                        (pos0 + 1, tlen, prev.pos1, prev.tlen)
                    } else {
                        (prev.pos1, prev.tlen, pos0 + 1, tlen)
                    };
                    pbc[tid].add_pe(p1, t1, p2, t2);
                    // DupFreq fingerprint: (chrom_id, leftmost_pos, |isize|).
                    let leftpos = p1.min(p2);
                    let abs_isize = t1.abs().max(t2.abs());
                    dup.add_pe(tid as u32, leftpos, abs_isize);
                } else {
                    mate_buffer.insert(q.to_vec(), FirstMateInfo { pos1: pos0 + 1, tlen });
                }
            }
        }
    }

    eprintln!("[rustqc atac] processed {} primary records", n_records);

    // ── Finalize metrics ──────────────────────────────────────────────────────
    let bq_report = bam_qc::finalize(&bq, &pbc);
    let tsse_result = tsse::compute(&tss_cov);
    let nfr_rows = nfr_score::compute(&tss_cov);
    let pt_rows = pt_score::compute(&tss_cov);
    let frag_rows = frag.finalize();
    let dup_hist = dup.histogram();
    let lib_rows = lib_complexity::estimate(&dup_hist, 100)
        .context("library complexity estimation failed")?;

    // Compute median NFR/PT scores for JSON summary.
    let nfr_median = median(nfr_rows.iter().map(|r| r.nfr_score).collect());
    let pt_median = median(pt_rows.iter().map(|r| r.pt_score).collect());

    // ── Create output directories ─────────────────────────────────────────────
    let outdir = &cfg.outdir;
    let flat = cfg.flat_output;

    let mk = |sub: &str| -> Result<()> {
        if !flat {
            fs::create_dir_all(Path::new(outdir).join(sub))
                .with_context(|| format!("create_dir_all {}/{}", outdir, sub))?;
        } else {
            fs::create_dir_all(outdir)
                .with_context(|| format!("create_dir_all {}", outdir))?;
        }
        Ok(())
    };

    fs::create_dir_all(outdir)
        .with_context(|| format!("create_dir_all {}", outdir))?;
    mk("bamqc")?;
    mk("fragsize")?;
    mk("tsse")?;
    mk("nfr")?;
    mk("pt")?;
    mk("lib_complexity")?;

    // ── Write bamQC TSV ───────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "bamqc", &format!("{}.bamqc.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "sample\ttotal_qnames\tduplicate_rate\tmitochondria_rate\tproper_pair_rate\tunmapped_rate\thas_unmapped_mate_rate\tnot_passing_qc_rate\tnrf\tpbc1\tpbc2")?;
        writeln!(
            w,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            sample,
            bq_report.total_qnames,
            bq_report.duplicate_rate,
            bq_report.mitochondria_rate,
            bq_report.proper_pair_rate,
            bq_report.unmapped_rate,
            bq_report.has_unmapped_mate_rate,
            bq_report.not_passing_qc_rate,
            bq_report.nrf,
            bq_report.pbc1,
            bq_report.pbc2,
        )?;
    }

    // ── Write MAPQ histogram TSV ──────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "bamqc", &format!("{}.mapq.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "mapq\tcount")?;
        for (q, c) in &bq_report.mapq_hist {
            writeln!(w, "{}\t{}", q, c)?;
        }
    }

    // ── Write fragSize TSV ────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "fragsize", &format!("{}.fragsize.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        frag_size::write_tsv(&mut w, &frag_rows)?;
    }

    // ── Write TSSE TSV ────────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "tsse", &format!("{}.tsse.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "window_idx\tnorm_signal")?;
        for (i, v) in tsse_result.values.iter().enumerate() {
            writeln!(w, "{}\t{:.8}", i + 1, v)?;
        }
    }

    // ── Write NFR TSV ─────────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "nfr", &format!("{}.nfr.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "tss_idx\tn1\tnf\tn2\tnfr_score\tlog2meancov")?;
        for r in &nfr_rows {
            writeln!(w, "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                r.tss_idx, r.n1, r.nf, r.n2, r.nfr_score, r.log2_mean_cov)?;
        }
    }

    // ── Write PT TSV ──────────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "pt", &format!("{}.pt.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "tss_idx\tpromoter\tbody\tpt_score\tlog2meancov")?;
        for r in &pt_rows {
            writeln!(w, "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                r.tss_idx, r.promoter, r.body, r.pt_score, r.log2_mean_cov)?;
        }
    }

    // ── Write lib_complexity TSV ──────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "lib_complexity", &format!("{}.libcomplexity.tsv", sample));
        let mut w = BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "relative_size\tdistinct_fragments\tputative_reads")?;
        for r in &lib_rows {
            writeln!(w, "{:.2}\t{:.2}\t{:.2}", r.relative_size, r.distinct_fragments, r.putative_reads)?;
        }
    }

    // ── Write SVG plots ───────────────────────────────────────────────────────
    {
        let p = metric_path(outdir, flat, "fragsize", &format!("{}.fragsize.svg", sample));
        plots::fragsize_svg(&frag_rows, &p, &sample)
            .with_context(|| format!("fragsize SVG: {}", p.display()))?;
    }
    {
        let p = metric_path(outdir, flat, "tsse", &format!("{}.tsse.svg", sample));
        plots::tsse_svg(&tsse_result.values, &p, &sample)
            .with_context(|| format!("TSSE SVG: {}", p.display()))?;
    }
    {
        let p = metric_path(outdir, flat, "lib_complexity", &format!("{}.libcomplexity.svg", sample));
        plots::lib_complexity_svg(&lib_rows, &p, &sample)
            .with_context(|| format!("lib_complexity SVG: {}", p.display()))?;
    }

    // ── Build and write JSON summary ──────────────────────────────────────────
    // Derive relative TSV paths (relative to outdir).
    let fragsize_tsv_path = if flat {
        format!("{}.fragsize.tsv", sample)
    } else {
        format!("fragsize/{}.fragsize.tsv", sample)
    };
    let tsse_tsv_path = if flat {
        format!("{}.tsse.tsv", sample)
    } else {
        format!("tsse/{}.tsse.tsv", sample)
    };
    let nfr_tsv_path = if flat {
        format!("{}.nfr.tsv", sample)
    } else {
        format!("nfr/{}.nfr.tsv", sample)
    };
    let pt_tsv_path = if flat {
        format!("{}.pt.tsv", sample)
    } else {
        format!("pt/{}.pt.tsv", sample)
    };
    let libcomplexity_tsv_path = if flat {
        format!("{}.libcomplexity.tsv", sample)
    } else {
        format!("lib_complexity/{}.libcomplexity.tsv", sample)
    };

    // extrapolated_total: row where relative_size == 1.0; None when NaN.
    let extrapolated_total: Option<f64> = lib_rows
        .iter()
        .find(|r| (r.relative_size - 1.0).abs() < 1e-9)
        .and_then(|r| if r.distinct_fragments.is_nan() { None } else { Some(r.distinct_fragments) });

    let atac_summary = summary::AtacSummary {
        schema_version: "1.0".to_string(),
        sample: sample.clone(),
        tool_versions: summary::ToolVersions {
            rustqc: env!("CARGO_PKG_VERSION").to_string(),
            atacseqqc_replicates: "1.36.0".to_string(),
        },
        split_method: "fixed_intervals_v1",
        bamqc: {
            let mut mapq_histogram = serde_json::Map::new();
            for (k, v) in &bq_report.mapq_hist {
                mapq_histogram.insert(k.to_string(), serde_json::Value::Number((*v).into()));
            }
            summary::BamqcSection {
                total_qnames: bq_report.total_qnames,
                duplicate_rate: bq_report.duplicate_rate,
                mitochondria_rate: bq_report.mitochondria_rate,
                proper_pair_rate: bq_report.proper_pair_rate,
                unmapped_rate: bq_report.unmapped_rate,
                has_unmapped_mate_rate: bq_report.has_unmapped_mate_rate,
                not_passing_qc_rate: bq_report.not_passing_qc_rate,
                nrf: bq_report.nrf,
                pbc1: bq_report.pbc1,
                pbc2: bq_report.pbc2,
                mapq_histogram,
            }
        },
        fragsize: summary::FragsizeSection {
            total_pairs: frag_rows.iter().map(|(_, c, _)| c).sum(),
            tsv_path: fragsize_tsv_path,
        },
        tsse: summary::TsseSection {
            score: tsse_result.tsse_score,
            n_windows: tsse_result.values.len() as u32,
            values: tsse_result.values.clone(),
            tsv_path: tsse_tsv_path,
        },
        nfr: summary::ScoreSection {
            n_tss: nfr_rows.len() as u32,
            median_score: nfr_median,
            tsv_path: nfr_tsv_path,
        },
        pt: summary::ScoreSection {
            n_tss: pt_rows.len() as u32,
            median_score: pt_median,
            tsv_path: pt_tsv_path,
        },
        lib_complexity: summary::LibComplexitySection {
            n_rows: lib_rows.len() as u32,
            extrapolated_total,
            tsv_path: libcomplexity_tsv_path,
        },
    };

    // Determine JSON output path.
    let json_path_str = cfg.json_summary.clone().unwrap_or_else(|| {
        Path::new(outdir)
            .join(format!("{}.atac.summary.json", sample))
            .to_string_lossy()
            .to_string()
    });

    if json_path_str == "-" {
        serde_json::to_writer_pretty(std::io::stdout(), &atac_summary)
            .context("write JSON summary to stdout")?;
        println!();
    } else {
        let p = Path::new(&json_path_str);
        let f = File::create(p).with_context(|| format!("create JSON {}", p.display()))?;
        serde_json::to_writer_pretty(BufWriter::new(f), &atac_summary)
            .context("write JSON summary")?;
        eprintln!("[rustqc atac] wrote JSON summary: {}", p.display());
    }

    // Wire EmitWriters scaffold (Phase 14 implements actual file writes).
    let _ = (cfg.emit_shifted_bam, cfg.emit_split_bams);

    eprintln!("[rustqc atac] done");
    Ok(())
}

#[derive(Debug, Clone)]
pub struct ResolvedAtacConfig {
    pub inputs: Vec<String>,
    pub gtf: String,
    /// Reserved for CRAM decoding (Phase N+1).
    #[allow(dead_code)]
    pub reference: Option<String>,
    pub outdir: String,
    pub sample_name: Option<String>,
    pub flat_output: bool,
    pub json_summary: Option<String>,
    pub mito_chrom: Option<String>, // None ⇒ auto-detect at runtime
    pub tsse_flank: u32,
    pub emit_shifted_bam: bool,
    pub emit_split_bams: bool,
    /// Reserved: rayon worker count for the streaming driver (Phase N+1).
    #[allow(dead_code)]
    pub threads: usize,
    /// Reserved: MAPQ threshold; informational only today (matches ATACseqQC,
    /// which does not filter by MAPQ globally). Wired in if/when filtering is added.
    #[allow(dead_code)]
    pub mapq_cut: u8,
    /// Reserved for future logging verbosity control.
    #[allow(dead_code)]
    pub quiet: bool,
    /// Reserved for future logging verbosity control.
    #[allow(dead_code)]
    pub verbose: bool,
}

const DEFAULT_TSSE_FLANK: u32 = 1000;

/// Merge CLI flags with YAML `AtacConfig`, applying defaults as a last resort.
///
/// Precedence (highest → lowest):
/// - CLI flag present → use CLI value
/// - YAML field set   → use YAML value
/// - Neither          → built-in default
///
/// For `Option<X>` fields: `args.X.clone().or(atac_cfg.X.clone()).unwrap_or(default)`
/// For `bool` fields:      `args.X || atac_cfg.X`  (either source enabling → enabled)
pub fn resolve(args: &AtacArgs, atac_cfg: &AtacConfig) -> ResolvedAtacConfig {
    ResolvedAtacConfig {
        inputs: args.input.clone(),
        gtf: args.gtf.clone(),
        reference: args.reference.clone(),
        outdir: args.outdir.clone(),
        sample_name: args.sample_name.clone(),
        flat_output: args.flat_output,
        json_summary: args.json_summary.clone(),
        mito_chrom: args.mito_chrom.clone().or_else(|| atac_cfg.mito_chrom.clone()),
        tsse_flank: args
            .tsse_flank
            .or(atac_cfg.tsse_flank)
            .unwrap_or(DEFAULT_TSSE_FLANK),
        emit_shifted_bam: args.emit_shifted_bam || atac_cfg.emit_shifted_bam,
        emit_split_bams: args.emit_split_bams || atac_cfg.emit_split_bams,
        threads: args.threads,
        mapq_cut: args.mapq_cut,
        quiet: args.quiet,
        verbose: args.verbose,
    }
}

pub mod bam_qc;
pub mod bam_writer;
pub mod frag_size;
pub mod lib_complexity;
pub mod loess;
pub mod mito;
pub mod nfr_score;
pub mod pe_check;
pub mod plots;
pub mod pt_score;
pub mod shift;
pub mod split;
pub mod summary;
pub mod tss_cov;
pub mod tsse;

/// Return the effective TSS coverage flank, enforcing the PTscore body requirement.
///
/// PTscore needs `[TSS - 2000, TSS + 500 + body]` — the promoter + gene-body window
/// comfortably fits within 3000 bp of each TSS. The user-supplied `tsse_flank` is
/// raised to 3000 if lower.
pub fn resolve_flank(tsse_flank: u32) -> u32 {
    const PT_REQUIREMENT: u32 = 3000;
    tsse_flank.max(PT_REQUIREMENT)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Commands};
    use clap::Parser;

    fn parse(args: &[&str]) -> AtacArgs {
        match Cli::parse_from(args).command {
            Commands::Atac(a) => a,
            _ => panic!("expected Atac"),
        }
    }

    #[test]
    fn flank_floor_at_3000() {
        assert_eq!(resolve_flank(1000), 3000);
        assert_eq!(resolve_flank(5000), 5000);
    }

    #[test]
    fn resolve_applies_defaults() {
        let r = resolve(
            &parse(&["rustqc", "atac", "x.bam", "--gtf", "g.gtf"]),
            &AtacConfig::default(),
        );
        assert_eq!(r.tsse_flank, DEFAULT_TSSE_FLANK);
        assert_eq!(r.threads, 1);
        assert_eq!(r.mapq_cut, 30);
        assert!(!r.emit_shifted_bam);
        assert!(r.mito_chrom.is_none());
    }

    #[test]
    fn resolve_passes_through_overrides() {
        let r = resolve(
            &parse(&[
                "rustqc",
                "atac",
                "x.bam",
                "--gtf",
                "g.gtf",
                "--mito-chrom",
                "MT",
                "--tsse-flank",
                "500",
            ]),
            &AtacConfig::default(),
        );
        assert_eq!(r.tsse_flank, 500);
        assert_eq!(r.mito_chrom.as_deref(), Some("MT"));
    }

    #[test]
    fn resolve_yaml_config_overrides_default_when_cli_absent() {
        let atac_cfg = AtacConfig {
            mito_chrom: Some("XYZ".into()),
            tsse_flank: Some(2500),
            emit_shifted_bam: true,
            emit_split_bams: false,
        };
        let r = resolve(
            &parse(&["rustqc", "atac", "x.bam", "--gtf", "g.gtf"]),
            &atac_cfg,
        );
        assert_eq!(r.mito_chrom.as_deref(), Some("XYZ"));
        assert_eq!(r.tsse_flank, 2500);
        assert!(r.emit_shifted_bam);
    }

    #[test]
    fn resolve_cli_args_override_yaml() {
        let atac_cfg = AtacConfig {
            mito_chrom: Some("XYZ".into()),
            tsse_flank: Some(2500),
            emit_shifted_bam: true,
            emit_split_bams: false,
        };
        let r = resolve(
            &parse(&[
                "rustqc",
                "atac",
                "x.bam",
                "--gtf",
                "g.gtf",
                "--mito-chrom",
                "MT",
                "--tsse-flank",
                "500",
            ]),
            &atac_cfg,
        );
        assert_eq!(r.mito_chrom.as_deref(), Some("MT"));
        assert_eq!(r.tsse_flank, 500);
        // emit_shifted_bam is true in YAML and CLI doesn't set it, so still true
        assert!(r.emit_shifted_bam);
    }
}
