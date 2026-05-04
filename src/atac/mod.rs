//! ATAC-Seq quality control and Tn5 preprocessing.
//!
//! Implements bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, and library
//! complexity, plus optional +4/-5 Tn5 shift and fixed-interval NFR/mono/di/tri
//! BAM split. Numerical fidelity targets ATACseqQC 1.36.0.

use anyhow::{Context, Result};

use crate::cli::AtacArgs;
use crate::config::AtacConfig;

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
    for input in &cfg.inputs {
        pe_check::assert_paired_end(std::path::Path::new(input))
            .with_context(|| format!("paired-end check failed for {}", input))?;
    }
    anyhow::bail!("rustqc atac is not yet implemented (PE check passed; metrics pending)");
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ResolvedAtacConfig {
    pub inputs: Vec<String>,
    pub gtf: String,
    pub reference: Option<String>,
    pub outdir: String,
    pub sample_name: Option<String>,
    pub flat_output: bool,
    pub json_summary: Option<String>,
    pub mito_chrom: Option<String>, // None ⇒ auto-detect at runtime
    pub tsse_flank: u32,
    pub emit_shifted_bam: bool,
    pub emit_split_bams: bool,
    pub threads: usize,
    pub mapq_cut: u8,
    pub quiet: bool,
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
pub mod frag_size;
pub mod lib_complexity;
pub mod loess;
pub mod mito;
pub mod nfr_score;
pub mod pe_check;
pub mod pt_score;
pub mod tss_cov;
pub mod tsse;

#[allow(dead_code)]
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
