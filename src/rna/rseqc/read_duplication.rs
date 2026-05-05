//! Read duplication rate analysis (reimplementation of RSeQC's read_duplication.py).
//!
//! Computes position-based and sequence-based duplication histograms from
//! a BAM/SAM/CRAM file in a single pass.

use anyhow::{Context, Result};
use log::debug;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

// ===================================================================
// Result structures
// ===================================================================

/// Histogram of duplication rates.
///
/// Maps occurrence count → number of distinct sequences/positions
/// observed at that count.
pub type DupHistogram = BTreeMap<u64, u64>;

/// Results from read duplication analysis.
#[derive(Debug)]
pub struct ReadDuplicationResult {
    /// Position-based duplication histogram.
    pub pos_histogram: DupHistogram,
    /// Sequence-based duplication histogram.
    pub seq_histogram: DupHistogram,
}

// ===================================================================
// Output
// ===================================================================

/// Write a duplication histogram to a TSV file in RSeQC format.
///
/// Format: `Occurrence\tUniqReadNumber\n`
///
/// # Arguments
/// * `histogram` - The duplication histogram to write
/// * `path` - Output file path
fn write_histogram(histogram: &DupHistogram, path: &Path) -> Result<()> {
    let mut file = std::fs::File::create(path)
        .with_context(|| format!("Failed to create output file: {}", path.display()))?;

    writeln!(file, "Occurrence\tUniqReadNumber")?;
    for (occurrence, uniq_count) in histogram {
        writeln!(file, "{}\t{}", occurrence, uniq_count)?;
    }

    Ok(())
}

/// Write read duplication results to output files.
///
/// Creates two files:
/// - `{stem}.pos.DupRate.xls` — position-based duplication histogram
/// - `{stem}.seq.DupRate.xls` — sequence-based duplication histogram
///
/// # Arguments
/// * `result` - The duplication analysis results
/// * `outdir` - Output directory
/// * `stem` - File name stem (typically BAM file name without extension)
pub fn write_read_duplication(
    result: &ReadDuplicationResult,
    outdir: &Path,
    stem: &str,
) -> Result<()> {
    let pos_path = outdir.join(format!("{}.pos.DupRate.xls", stem));
    let seq_path = outdir.join(format!("{}.seq.DupRate.xls", stem));

    write_histogram(&result.pos_histogram, &pos_path).with_context(|| {
        format!(
            "Failed to write position duplication file: {}",
            pos_path.display()
        )
    })?;
    debug!(
        "Wrote position-based duplication rates to {}",
        pos_path.display()
    );

    write_histogram(&result.seq_histogram, &seq_path).with_context(|| {
        format!(
            "Failed to write sequence duplication file: {}",
            seq_path.display()
        )
    })?;
    debug!(
        "Wrote sequence-based duplication rates to {}",
        seq_path.display()
    );

    // Write the R plotting script (matching upstream read_duplication.py output)
    let r_path = outdir.join(format!("{}.DupRate_plot.r", stem));
    write_r_script(&result.pos_histogram, &result.seq_histogram, stem, &r_path)
        .with_context(|| format!("Failed to write duplication R script: {}", r_path.display()))?;
    debug!(
        "Wrote duplication R plotting script to {}",
        r_path.display()
    );

    Ok(())
}

/// Write the R plotting script matching upstream RSeQC read_duplication.py.
fn write_r_script(
    pos_hist: &BTreeMap<u64, u64>,
    seq_hist: &BTreeMap<u64, u64>,
    stem: &str,
    path: &Path,
) -> Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(path)?;

    writeln!(f, "pdf('{stem}.DupRate_plot.pdf')")?;
    writeln!(f, "par(mar=c(5,4,4,5),las=0)")?;

    // seq_occ and seq_uniqRead vectors
    let seq_occ: Vec<String> = seq_hist.keys().map(|k| k.to_string()).collect();
    let seq_vals: Vec<String> = seq_hist.values().map(|v| v.to_string()).collect();
    writeln!(f, "seq_occ=c({})", seq_occ.join(","))?;
    writeln!(f, "seq_uniqRead=c({})", seq_vals.join(","))?;

    // pos_occ and pos_uniqRead vectors
    let pos_occ: Vec<String> = pos_hist.keys().map(|k| k.to_string()).collect();
    let pos_vals: Vec<String> = pos_hist.values().map(|v| v.to_string()).collect();
    writeln!(f, "pos_occ=c({})", pos_occ.join(","))?;
    writeln!(f, "pos_uniqRead=c({})", pos_vals.join(","))?;

    // Plot commands matching upstream exactly
    writeln!(f, "plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Occurrence of read',pch=4,cex=0.8,col='blue',xlim=c(1,500),yaxt='n')")?;
    writeln!(
        f,
        "points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')"
    )?;
    writeln!(f, "ym=floor(max(log10(pos_uniqRead)))")?;
    writeln!(
        f,
        "legend(300,ym,legend=c('Mapping-based','Sequence-based'),col=c('blue','red'),pch=c(4,20))"
    )?;
    writeln!(f, "axis(side=2,at=0:ym,labels=0:ym)")?;

    // Right-axis annotation: percentage labels at first 4 positions
    writeln!(f, "axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))")?;
    writeln!(f, "mtext(4, text = \"Reads %\", line = 2)")?;
    writeln!(f, "dev.off()")?;

    Ok(())
}

// ===================================================================
// Tests
// ===================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_io::{self as bam, CigarKind};
    use log::debug;
    use std::collections::HashMap;
    use std::time::Instant;

    /// Build the position key for a read from its CIGAR alignment.
    ///
    /// Constructs `{chrom}:{start}:{exon1_start}-{exon1_end}:{exon2_start}-{exon2_end}:...`
    /// matching RSeQC's `fetch_exon` + position key logic.
    fn build_position_key(chrom: &str, pos: i64, cigar_ops: &[(CigarKind, u32)]) -> String {
        use CigarKind as K;

        let mut key = format!("{}:{}:", chrom, pos);
        let mut ref_pos = pos;

        for (kind, len) in cigar_ops {
            let len = *len as i64;
            match kind {
                K::Match | K::SequenceMatch | K::SequenceMismatch => {
                    let end = ref_pos + len;
                    key.push_str(&format!("{}-{}:", ref_pos, end));
                    ref_pos = end;
                }
                K::Deletion | K::Skip => {
                    ref_pos += len;
                }
                K::SoftClip => {
                    ref_pos += len;
                }
                K::Insertion | K::HardClip | K::Pad => {}
            }
        }

        key
    }

    /// Standalone read duplication analysis (retained for tests only).
    fn read_duplication(bam_path: &str, mapq_cut: u8) -> Result<ReadDuplicationResult> {
        let start = Instant::now();

        let (mut reader, header) = bam::open(Path::new(bam_path))
            .with_context(|| format!("Failed to open BAM file: {}", bam_path))?;

        let target_names: Vec<String> = bam::reference_sequences(&header)
            .into_iter()
            .map(|(name, _)| name)
            .collect();

        let mut seq_dup: HashMap<Vec<u8>, u64> = HashMap::new();
        let mut pos_dup: HashMap<String, u64> = HashMap::new();
        let mut total_processed = 0u64;

        for result in reader.records() {
            let record = result.context("Failed to read BAM record")?;
            let flags = record.flags();

            if flags.is_unmapped() || flags.is_qc_fail() || bam::mapq(&record) < mapq_cut {
                continue;
            }

            total_processed += 1;

            // Decode the sequence to ASCII uppercase bases for the seq dedup map.
            let seq_len = record.sequence().len();
            let mut seq_upper = Vec::with_capacity(seq_len);
            for i in 0..seq_len {
                let nibble = bam::encoded_base(&record, i);
                let ch = match nibble {
                    1 => b'A',
                    2 => b'C',
                    4 => b'G',
                    8 => b'T',
                    15 => b'N',
                    _ => b'N',
                };
                seq_upper.push(ch);
            }
            *seq_dup.entry(seq_upper).or_insert(0) += 1;

            let tid = bam::tid(&record);
            if tid >= 0 && (tid as usize) < target_names.len() {
                let chrom = &target_names[tid as usize];
                let pos = bam::pos_0based(&record);
                let cigar_ops = bam::cigar_ops(&record).context("malformed CIGAR")?;
                let key = build_position_key(chrom, pos, &cigar_ops);
                *pos_dup.entry(key).or_insert(0) += 1;
            }
        }

        debug!(
            "Read duplication: processed {} reads, {} distinct sequences, {} distinct positions",
            total_processed,
            seq_dup.len(),
            pos_dup.len()
        );

        let mut pos_histogram: DupHistogram = BTreeMap::new();
        for count in pos_dup.values() {
            *pos_histogram.entry(*count).or_insert(0) += 1;
        }

        let mut seq_histogram: DupHistogram = BTreeMap::new();
        for count in seq_dup.values() {
            *seq_histogram.entry(*count).or_insert(0) += 1;
        }

        let elapsed = start.elapsed();
        debug!(
            "Read duplication analysis complete in {:.1}s ({} reads processed)",
            elapsed.as_secs_f64(),
            total_processed
        );

        Ok(ReadDuplicationResult {
            pos_histogram,
            seq_histogram,
        })
    }

    #[test]
    fn test_build_position_key_simple() {
        let cigar_ops = vec![(CigarKind::Match, 50u32)];
        let key = build_position_key("chr1", 1000, &cigar_ops);
        assert_eq!(key, "chr1:1000:1000-1050:");
    }

    #[test]
    fn test_build_position_key_spliced() {
        let cigar_ops = vec![
            (CigarKind::Match, 10u32),
            (CigarKind::Skip, 500u32),
            (CigarKind::Match, 20u32),
        ];
        let key = build_position_key("chr1", 1000, &cigar_ops);
        assert_eq!(key, "chr1:1000:1000-1010:1510-1530:");
    }

    #[test]
    fn test_read_duplication_small() {
        let bam_path = "tests/data/test.bam";
        if !Path::new(bam_path).exists() {
            return; // Skip if test data not available
        }

        let result = read_duplication(bam_path, 30).unwrap();

        assert!(
            !result.pos_histogram.is_empty(),
            "Position histogram should not be empty"
        );
        assert!(
            !result.seq_histogram.is_empty(),
            "Sequence histogram should not be empty"
        );

        let max_pos_occ = *result.pos_histogram.keys().last().unwrap();
        assert!(max_pos_occ >= 1, "Should have at least occurrence=1");

        let max_seq_occ = *result.seq_histogram.keys().last().unwrap();
        assert!(max_seq_occ >= 1, "Should have at least occurrence=1");
    }
}
