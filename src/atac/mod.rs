//! ATAC-Seq quality control and Tn5 preprocessing.
//!
//! Implements bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, and library
//! complexity, plus optional +4/-5 Tn5 shift and fixed-interval NFR/mono/di/tri
//! BAM split. Numerical fidelity targets ATACseqQC 1.36.0.

use anyhow::{Context, Result};

use crate::cli::AtacArgs;

/// Entry point for the `rustqc atac` subcommand.
pub fn run(args: AtacArgs) -> Result<()> {
    let cfg = resolve(&args);
    for input in &cfg.inputs {
        pe_check::assert_paired_end(std::path::Path::new(input))
            .with_context(|| format!("paired-end check failed for {}", input))?;
    }
    anyhow::bail!("rustqc atac is not yet implemented (PE check passed; metrics pending)");
}

#[derive(Debug, Clone)]
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

pub fn resolve(args: &AtacArgs) -> ResolvedAtacConfig {
    ResolvedAtacConfig {
        inputs: args.input.clone(),
        gtf: args.gtf.clone(),
        reference: args.reference.clone(),
        outdir: args.outdir.clone(),
        sample_name: args.sample_name.clone(),
        flat_output: args.flat_output,
        json_summary: args.json_summary.clone(),
        mito_chrom: args.mito_chrom.clone(),
        tsse_flank: args.tsse_flank.unwrap_or(DEFAULT_TSSE_FLANK),
        emit_shifted_bam: args.emit_shifted_bam,
        emit_split_bams: args.emit_split_bams,
        threads: args.threads,
        mapq_cut: args.mapq_cut,
        quiet: args.quiet,
        verbose: args.verbose,
    }
}

pub mod mito;
pub mod pe_check;

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
    fn resolve_applies_defaults() {
        let r = resolve(&parse(&["rustqc", "atac", "x.bam", "--gtf", "g.gtf"]));
        assert_eq!(r.tsse_flank, DEFAULT_TSSE_FLANK);
        assert_eq!(r.threads, 1);
        assert_eq!(r.mapq_cut, 30);
        assert!(!r.emit_shifted_bam);
        assert!(r.mito_chrom.is_none());
    }

    #[test]
    fn resolve_passes_through_overrides() {
        let r = resolve(&parse(&[
            "rustqc",
            "atac",
            "x.bam",
            "--gtf",
            "g.gtf",
            "--mito-chrom",
            "MT",
            "--tsse-flank",
            "500",
        ]));
        assert_eq!(r.tsse_flank, 500);
        assert_eq!(r.mito_chrom.as_deref(), Some("MT"));
    }
}
