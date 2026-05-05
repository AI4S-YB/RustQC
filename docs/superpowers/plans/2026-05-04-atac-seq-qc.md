# RustQC ATAC-seq QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `rustqc atac` subcommand parallel to `rustqc rna` providing single-pass ATAC-seq QC (bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, library complexity) plus opt-in Tn5 +4/−5 shift and fixed-interval NFR/mono/di/tri BAM split, with numerical fidelity to ATACseqQC 1.36.0.

**Architecture:** Phase 1 lifts shared BAM/preseq infra out of `src/rna/` to `src/` with no behavior change. Phase 2 scaffolds the `rustqc atac` subcommand. Phases 3–9 implement metrics test-first against fixtures from ATACseqQC's own GL1–GL4 BAMs. Phases 10–12 add Tn5 shift, split, and opt-in BAM emission. Phase 13 adds plots and JSON summary. Phase 14 wires golden-output integration tests against R reference outputs (R run offline; goldens committed). Phase 15 updates docs site.

**Tech Stack:** Rust 2021, clap, noodles (BAM/SAM/BGZF), rayon, plotters, serde_yaml_ng, anyhow. New tests use `#[test]` + `assert_float_absolute_eq!`-style helpers (rolled inline as `fn approx_eq`).

**Reference spec:** `docs/superpowers/specs/2026-05-04-atac-seq-qc-design.md`

---

## Pre-flight

- [ ] **Step 0.1: Confirm clean working tree**

Run: `git -C /home/xzg/project/RustQC status -sb`
Expected: clean tree on `main` (no staged/unstaged changes; `ATACseqQC_1.36.0.tar.gz` already committed-or-ignored is fine since it lives at repo root and is in `.gitignore` exclusion if listed).

- [ ] **Step 0.2: Confirm baseline tests pass**

Run: `cargo test --workspace --release 2>&1 | tail -30`
Expected: all tests green. Capture the test count for the §2 regression check.

- [ ] **Step 0.3: Make raw ATACseqQC fixtures available for later phases**

The `ATACseqQC_1.36.0.tar.gz` at the repo root contains `inst/extdata/GL{1..4}.bam(.bai)` and `inst/extdata/splited/{NucleosomeFree,mononucleosome,dinucleosome,trinucleosome}.bam(.bai)`. Phase 14 will extract these into `tests/data/atac/`. For now just verify the tarball is readable.

Run: `tar -tzf /home/xzg/project/RustQC/ATACseqQC_1.36.0.tar.gz | grep -E "extdata/(GL[1-4]|splited)" | head -20`
Expected: 8+ paths listed (GL1.bam, GL1.bam.bai, …, splited/NucleosomeFree.bam, …).

---

## Phase 1 — Shared infra refactor (no behavior change)

**Goal:** Lift `bam_flags`, `bam_io`, `cpp_rng`, `preseq` from `src/rna/` to `src/`. Verify RNA outputs are byte-identical to pre-refactor.

### Task 1.1: Move `bam_flags.rs` to crate root

**Files:**
- Move: `src/rna/bam_flags.rs` → `src/bam_flags.rs`
- Modify: `src/rna/mod.rs` (drop `pub mod bam_flags;`)
- Modify: `src/main.rs` (add `pub mod bam_flags;`)
- Modify: every file that uses `crate::rna::bam_flags` (rewrite path)

- [ ] **Step 1: Move the file**

```bash
git -C /home/xzg/project/RustQC mv src/rna/bam_flags.rs src/bam_flags.rs
```

- [ ] **Step 2: Drop module declaration in `src/rna/mod.rs`**

Edit `src/rna/mod.rs`: remove the line `pub mod bam_flags;`.

- [ ] **Step 3: Add module declaration in `src/main.rs`**

Find the `mod` declarations near the top of `src/main.rs` (look for `mod rna;`, `mod cli;` etc.) and add `mod bam_flags;` alongside them.

- [ ] **Step 4: Rewrite all `crate::rna::bam_flags` imports to `crate::bam_flags`**

Run: `grep -rln "crate::rna::bam_flags" /home/xzg/project/RustQC/src/`

For each file in the result, edit to replace `crate::rna::bam_flags` with `crate::bam_flags`.

- [ ] **Step 5: Build**

Run: `cargo build 2>&1 | tail -20`
Expected: clean build, no errors, no new warnings.

- [ ] **Step 6: Test**

Run: `cargo test 2>&1 | tail -20`
Expected: same passing test count as Step 0.2.

- [ ] **Step 7: Commit (defer until Task 1.5 — refactor is one logical commit)**

### Task 1.2: Move `cpp_rng.rs` to crate root

**Files:**
- Move: `src/rna/cpp_rng.rs` → `src/cpp_rng.rs`
- Modify: `src/rna/mod.rs`
- Modify: `src/main.rs`
- Modify: every file that uses `crate::rna::cpp_rng`

- [ ] **Step 1: Move + rewire (mirror Task 1.1)**

```bash
git -C /home/xzg/project/RustQC mv src/rna/cpp_rng.rs src/cpp_rng.rs
```

- [ ] **Step 2: Drop `pub mod cpp_rng;` from `src/rna/mod.rs`; add `mod cpp_rng;` to `src/main.rs`.**

- [ ] **Step 3: Rewrite imports**

Run: `grep -rln "crate::rna::cpp_rng" /home/xzg/project/RustQC/src/`

For each file, replace `crate::rna::cpp_rng` with `crate::cpp_rng`.

- [ ] **Step 4: Verify build + tests**

Run: `cargo test 2>&1 | tail -10`
Expected: green.

### Task 1.3: Move `bam_io.rs` to crate root

**Files:**
- Move: `src/rna/bam_io.rs` → `src/bam_io.rs`
- Modify: `src/rna/mod.rs`, `src/main.rs`
- Modify: ~14 files importing `crate::rna::bam_io`

- [ ] **Step 1: Sanity check — is `bam_io.rs` truly RNA-agnostic?**

Run: `grep -nE "junction|splice|XS:A|infer_experiment|tin\b|qualimap|dupradar" /home/xzg/project/RustQC/src/rna/bam_io.rs`

Expected: no hits. If hits appear, the RNA-specific snippet stays in `src/rna/bam_io_rna.rs` and only the generic part moves. (Do not split unless there is a real hit.)

- [ ] **Step 2: Move the file**

```bash
git -C /home/xzg/project/RustQC mv src/rna/bam_io.rs src/bam_io.rs
```

- [ ] **Step 3: Drop `pub mod bam_io;` from `src/rna/mod.rs`; add `mod bam_io;` to `src/main.rs`.**

- [ ] **Step 4: Rewrite imports**

Run: `grep -rln "crate::rna::bam_io" /home/xzg/project/RustQC/src/`

For each file, replace `crate::rna::bam_io` with `crate::bam_io`.

- [ ] **Step 5: Verify build + tests**

Run: `cargo test 2>&1 | tail -10`
Expected: green.

### Task 1.4: Move `preseq.rs` to crate root

**Files:**
- Move: `src/rna/preseq.rs` → `src/preseq.rs`
- Modify: `src/rna/mod.rs`, `src/main.rs`
- Modify: importers (mostly `src/rna/rseqc/accumulators.rs`)

- [ ] **Step 1: Move the file**

```bash
git -C /home/xzg/project/RustQC mv src/rna/preseq.rs src/preseq.rs
```

- [ ] **Step 2: Update internal `use crate::rna::bam_io` inside the moved file to `use crate::bam_io`**

Edit `src/preseq.rs` line 7: `use crate::rna::bam_io::{self as bam};` → `use crate::bam_io::{self as bam};`.

- [ ] **Step 3: Drop `pub mod preseq;` from `src/rna/mod.rs`; add `mod preseq;` to `src/main.rs`.**

- [ ] **Step 4: Rewrite external imports**

Run: `grep -rln "crate::rna::preseq" /home/xzg/project/RustQC/src/`

For each file, replace `crate::rna::preseq` with `crate::preseq`.

- [ ] **Step 5: Confirm preseq's public surface is reachable from outside `rna/`**

Run: `grep -nE "^pub (fn|struct) " /home/xzg/project/RustQC/src/preseq.rs`

Expected: `pub struct PreseqAccum`, `pub struct PreseqResult`, `pub fn estimate_complexity`, `pub fn write_output`, plus `into_histogram`/`finalize`/`n_distinct` on `PreseqAccum`. These are the surface ATAC's `lib_complexity.rs` will need; if any are `pub(crate)` or stricter, leave them — ATAC is in the same crate.

- [ ] **Step 6: Verify build + tests**

Run: `cargo test 2>&1 | tail -10`
Expected: green.

### Task 1.5: Regression bar — RNA outputs byte-identical

**Goal:** confirm the integration tests under `tests/integration_test.rs` produce the exact same artifacts as on `main` before the refactor.

- [ ] **Step 1: Run integration suite**

Run: `cargo test --release --test integration_test 2>&1 | tail -30`
Expected: all integration tests pass, including any byte-comparison assertions against fixtures in `tests/expected/`.

- [ ] **Step 2: Spot-check golden file equivalence**

Pick one or two representative golden files from `tests/expected/` (e.g. a featureCounts or rseqc output). Re-run `rustqc rna` against `tests/data/test.bam` with `--gtf tests/data/test.gtf` into a temp dir and `diff` against the expected file. If `tests/expected/` is consumed solely through `integration_test.rs`, Step 1 already covers this — the spot-check is belt-and-braces.

```bash
TMPOUT=$(mktemp -d)
./target/release/rustqc rna tests/data/test.bam --gtf tests/data/test.gtf --paired --outdir "$TMPOUT" >/dev/null
# Pick any golden under tests/expected/ and diff against the matching new output.
ls "$TMPOUT" && ls tests/expected/
```

If any diff appears, **stop**: the move was not behavior-preserving. Investigate before committing.

### Task 1.6: Commit Phase 1

- [ ] **Step 1: Stage and commit**

```bash
git -C /home/xzg/project/RustQC add -A
git -C /home/xzg/project/RustQC commit -m "$(cat <<'EOF'
refactor: extract shared BAM/preseq infra out of rna module

Move bam_flags.rs, bam_io.rs, cpp_rng.rs, and preseq.rs from src/rna/
to src/ so future library types (atac, chip) can reuse them without
depending on rna-specific modules. Pure file moves + import path
rewrites; no behavior change. RNA integration outputs verified
byte-identical pre/post.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 2: Verify commit**

Run: `git -C /home/xzg/project/RustQC log -1 --stat | head -30`
Expected: commit lists 4 file moves + the import-rewrite touches.

---

## Phase 2 — `rustqc atac` scaffolding

### Task 2.1: Add `Atac` subcommand to `cli.rs`

**Files:**
- Modify: `src/cli.rs`

- [ ] **Step 1: Add `AtacArgs` struct**

In `src/cli.rs`, after the `RnaArgs` struct (around line 395), add:

```rust
/// Arguments for the `atac` subcommand.
#[derive(Parser, Debug)]
#[command(
    next_line_help = false,
    term_width = 120,
    help_template = "\
{about-with-newline}
{usage-heading} {usage}

{all-args}"
)]
pub struct AtacArgs {
    // ── Input / Output ──────────────────────────────────────────────────
    /// Paired-end BAM/SAM/CRAM (single-end inputs are rejected at startup)
    #[arg(value_name = "INPUT", num_args = 1.., required = true, help_heading = "Input / Output")]
    pub input: Vec<String>,

    /// GTF gene annotation (plain or .gz); TSS coords source
    #[arg(short, long, value_name = "GTF", env = "RUSTQC_GTF", help_heading = "Input / Output")]
    pub gtf: String,

    /// Reference FASTA (required for CRAM)
    #[arg(short, long, value_name = "FASTA", env = "RUSTQC_REFERENCE", help_heading = "Input / Output")]
    pub reference: Option<String>,

    /// Output directory [default: .]
    #[arg(short, long, default_value = ".", hide_default_value = true, env = "RUSTQC_OUTDIR", help_heading = "Input / Output")]
    pub outdir: String,

    /// Override sample name (default: derived from BAM filename)
    #[arg(long, value_name = "NAME", env = "RUSTQC_SAMPLE_NAME", help_heading = "Input / Output")]
    pub sample_name: Option<String>,

    /// Write outputs to a flat directory (no subdirs)
    #[arg(long, default_value_t = false, env = "RUSTQC_FLAT_OUTPUT", help_heading = "Input / Output")]
    pub flat_output: bool,

    /// YAML configuration file
    #[arg(short, long, value_name = "CONFIG", help_heading = "Input / Output")]
    pub config: Option<String>,

    /// JSON summary path (use "-" for stdout)
    #[arg(short = 'j', long = "json-summary", value_name = "PATH", num_args = 0..=1, default_missing_value = "", env = "RUSTQC_JSON_SUMMARY", help_heading = "Input / Output")]
    pub json_summary: Option<String>,

    // ── ATAC-specific ───────────────────────────────────────────────────
    /// Mitochondrial chromosome name (default: auto-detect ^chrM$|^MT$|^Mito$)
    #[arg(long, value_name = "NAME", env = "RUSTQC_MITO_CHROM", help_heading = "ATAC-specific")]
    pub mito_chrom: Option<String>,

    /// Emit +4/-5 Tn5-shifted BAM
    #[arg(long, default_value_t = false, env = "RUSTQC_EMIT_SHIFTED_BAM", help_heading = "ATAC-specific")]
    pub emit_shifted_bam: bool,

    /// Emit NFR/mono/di/tri BAMs (fixed intervals)
    #[arg(long, default_value_t = false, env = "RUSTQC_EMIT_SPLIT_BAMS", help_heading = "ATAC-specific")]
    pub emit_split_bams: bool,

    /// TSSEscore flank (bp) [default: 1000]
    #[arg(long, value_name = "N", env = "RUSTQC_TSSE_FLANK", help_heading = "ATAC-specific")]
    pub tsse_flank: Option<u32>,

    // ── General ─────────────────────────────────────────────────────────
    /// Number of threads [default: 1]
    #[arg(short, long, default_value_t = 1, hide_default_value = true, env = "RUSTQC_THREADS", help_heading = "General")]
    pub threads: usize,

    /// MAPQ cutoff [default: 30]
    #[arg(short = 'Q', long = "mapq", default_value_t = 30, hide_default_value = true, env = "RUSTQC_MAPQ", help_heading = "General")]
    pub mapq_cut: u8,

    /// Suppress output except warnings/errors
    #[arg(short = 'q', long, conflicts_with = "verbose", env = "RUSTQC_QUIET", help_heading = "General")]
    pub quiet: bool,

    /// Show additional detail
    #[arg(short = 'v', long, conflicts_with = "quiet", env = "RUSTQC_VERBOSE", help_heading = "General")]
    pub verbose: bool,
}
```

- [ ] **Step 2: Add `Atac` variant to `Commands` enum**

In `src/cli.rs` find the `Commands` enum (around line 52) and add:

```rust
    /// ATAC-Seq QC — single-pass analysis of paired-end BAM/SAM/CRAM files.
    ///
    /// Runs bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, and library
    /// complexity in one pass. Requires a GTF annotation and paired-end input.
    /// Optionally emits Tn5-shifted and length-split BAMs.
    Atac(AtacArgs),
```

- [ ] **Step 3: Add a unit test mirroring `test_rna_default_args_gtf`**

Append to the `mod tests` block at the bottom of `src/cli.rs`:

```rust
    #[test]
    fn test_atac_default_args() {
        let cli = Cli::parse_from(["rustqc", "atac", "test.bam", "--gtf", "genes.gtf"]);
        match cli.command {
            Commands::Atac(args) => {
                assert_eq!(args.input, vec!["test.bam"]);
                assert_eq!(args.gtf, "genes.gtf");
                assert_eq!(args.outdir, ".");
                assert_eq!(args.threads, 1);
                assert_eq!(args.mapq_cut, 30);
                assert!(!args.emit_shifted_bam);
                assert!(!args.emit_split_bams);
                assert_eq!(args.mito_chrom, None);
                assert_eq!(args.tsse_flank, None);
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Atac subcommand"),
        }
    }

    #[test]
    fn test_atac_emit_flags() {
        let cli = Cli::parse_from([
            "rustqc", "atac", "test.bam", "--gtf", "genes.gtf",
            "--emit-shifted-bam", "--emit-split-bams",
            "--mito-chrom", "MT", "--tsse-flank", "2000",
        ]);
        match cli.command {
            Commands::Atac(args) => {
                assert!(args.emit_shifted_bam);
                assert!(args.emit_split_bams);
                assert_eq!(args.mito_chrom.as_deref(), Some("MT"));
                assert_eq!(args.tsse_flank, Some(2000));
            }
            #[allow(unreachable_patterns)]
            _ => panic!("Expected Atac subcommand"),
        }
    }
```

- [ ] **Step 4: Run new tests; expect failures (Atac variant not routed yet — but `Cli::parse_from` itself works; tests should pass purely on parser side)**

Run: `cargo test --lib cli::tests::test_atac 2>&1 | tail -10`
Expected: PASS. (No routing required to test parsing.)

- [ ] **Step 5: Commit**

```bash
git -C /home/xzg/project/RustQC add src/cli.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): add atac subcommand with CLI args"
```

### Task 2.2: Create empty `src/atac/` module and route `main.rs`

**Files:**
- Create: `src/atac/mod.rs`
- Modify: `src/main.rs`

- [ ] **Step 1: Create `src/atac/mod.rs`**

```rust
//! ATAC-Seq quality control and Tn5 preprocessing.
//!
//! Implements bamQC, fragSizeDist, TSSEscore, NFRscore, PTscore, and library
//! complexity, plus optional +4/-5 Tn5 shift and fixed-interval NFR/mono/di/tri
//! BAM split. Numerical fidelity targets ATACseqQC 1.36.0.

use anyhow::{Context, Result};

use crate::cli::AtacArgs;

/// Entry point for the `rustqc atac` subcommand.
pub fn run(args: AtacArgs) -> Result<()> {
    // Phase 2.3+ wires real implementation; for now just establish the contract.
    let _ = args;
    anyhow::bail!("rustqc atac is not yet implemented (scaffolding only)");
}
```

- [ ] **Step 2: Route `main.rs`**

In `src/main.rs`, find where `Commands::Rna(args)` is matched (search for `Commands::Rna`). Add a sibling arm:

```rust
        Commands::Atac(args) => atac::run(args)?,
```

Add `mod atac;` near the other top-level `mod` declarations.

- [ ] **Step 3: Build**

Run: `cargo build 2>&1 | tail -10`
Expected: clean build.

- [ ] **Step 4: Smoke-test the subcommand surfaces**

Run: `./target/debug/rustqc atac --help 2>&1 | head -40`
Expected: help text shows the ATAC arguments grouped under `Input / Output`, `ATAC-specific`, `General`.

- [ ] **Step 5: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/mod.rs src/main.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): scaffold atac subcommand routing"
```

### Task 2.3: Add `AtacConfig` to `config.rs`

**Files:**
- Modify: `src/config.rs`

- [ ] **Step 1: Read existing config layout**

Run: `grep -n "RnaConfig\|pub struct\|pub fn" /home/xzg/project/RustQC/src/config.rs | head -30`

- [ ] **Step 2: Add `AtacConfig` struct alongside `RnaConfig`**

In `src/config.rs`, near the `RnaConfig` definition, add:

```rust
/// ATAC-seq QC configuration (YAML-backed).
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AtacConfig {
    /// Mitochondrial chromosome name; auto-detected when None.
    pub mito_chrom: Option<String>,
    /// TSSEscore flank window in bp (default 1000).
    pub tsse_flank: Option<u32>,
    /// Emit Tn5-shifted BAM.
    pub emit_shifted_bam: bool,
    /// Emit NFR/mono/di/tri BAMs.
    pub emit_split_bams: bool,
}
```

If the top-level config struct in `config.rs` is something like `pub struct Config { pub rna: RnaConfig, … }`, add `pub atac: AtacConfig` to it.

- [ ] **Step 3: Build + commit**

Run: `cargo build 2>&1 | tail -5`
Expected: clean.

```bash
git -C /home/xzg/project/RustQC add src/config.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): add AtacConfig"
```

### Task 2.4: Resolved-config struct (CLI ⊕ YAML ⊕ defaults)

**Files:**
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Define a `ResolvedAtacConfig` plus `resolve_args` helper inside `src/atac/mod.rs`**

Append to `src/atac/mod.rs`:

```rust
#[derive(Debug, Clone)]
pub struct ResolvedAtacConfig {
    pub inputs: Vec<String>,
    pub gtf: String,
    pub reference: Option<String>,
    pub outdir: String,
    pub sample_name: Option<String>,
    pub flat_output: bool,
    pub json_summary: Option<String>,
    pub mito_chrom: Option<String>,    // None ⇒ auto-detect at runtime
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
```

- [ ] **Step 2: Add a unit test**

Append to `src/atac/mod.rs`:

```rust
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
            "rustqc", "atac", "x.bam", "--gtf", "g.gtf",
            "--mito-chrom", "MT", "--tsse-flank", "500",
        ]));
        assert_eq!(r.tsse_flank, 500);
        assert_eq!(r.mito_chrom.as_deref(), Some("MT"));
    }
}
```

- [ ] **Step 3: Test + commit**

Run: `cargo test atac:: 2>&1 | tail -10`
Expected: 2 new tests pass.

```bash
git -C /home/xzg/project/RustQC add src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): resolved-config struct"
```

### Task 2.5: Mito-chromosome auto-detection helper

**Files:**
- Create: `src/atac/mito.rs`
- Modify: `src/atac/mod.rs` (add `mod mito;`)

- [ ] **Step 1: Write failing test**

Create `src/atac/mito.rs`:

```rust
//! Mitochondrial chromosome detection from BAM @SQ names.

/// Names matched by the auto-detect regex (case-sensitive): `chrM`, `MT`, `Mito`.
pub fn detect_mito<'a>(seq_names: &'a [String]) -> Option<&'a str> {
    seq_names
        .iter()
        .find(|n| matches!(n.as_str(), "chrM" | "MT" | "Mito"))
        .map(String::as_str)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(strs: &[&str]) -> Vec<String> { strs.iter().map(|s| s.to_string()).collect() }

    #[test]
    fn detects_chrM() {
        assert_eq!(detect_mito(&s(&["chr1", "chr2", "chrM"])), Some("chrM"));
    }

    #[test]
    fn detects_MT_for_ensembl() {
        assert_eq!(detect_mito(&s(&["1", "2", "MT"])), Some("MT"));
    }

    #[test]
    fn detects_Mito_for_yeast() {
        assert_eq!(detect_mito(&s(&["I", "II", "Mito"])), Some("Mito"));
    }

    #[test]
    fn returns_none_when_absent() {
        assert_eq!(detect_mito(&s(&["chr1", "chr2"])), None);
    }

    #[test]
    fn does_not_match_substrings() {
        // chrMT or MTother should NOT match (we use exact equality on canonical names).
        assert_eq!(detect_mito(&s(&["chrMT"])), None);
        assert_eq!(detect_mito(&s(&["MTother"])), None);
    }
}
```

- [ ] **Step 2: Wire the module**

In `src/atac/mod.rs`, add `pub mod mito;` near the top.

- [ ] **Step 3: Run tests**

Run: `cargo test atac::mito:: 2>&1 | tail -10`
Expected: 5 tests pass.

- [ ] **Step 4: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/mito.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): mito chromosome auto-detection"
```

### Task 2.6: Paired-end startup check

**Files:**
- Create: `src/atac/pe_check.rs`
- Modify: `src/atac/mod.rs` (`pub mod pe_check;`)

ATAC requires paired-end input. We inspect the first up to 10 000 primary, mapped records; if none are paired, return an error.

- [ ] **Step 1: Write failing test**

Create `src/atac/pe_check.rs`:

```rust
//! Reject single-end BAMs at startup: scan the first ≤10 000 primary mapped
//! records and require at least one with the `READ_PAIRED` flag set.

use anyhow::{anyhow, Result};
use noodles_bam as bam;
use std::path::Path;

const MAX_RECORDS_TO_INSPECT: usize = 10_000;
const FLAG_READ_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub fn assert_paired_end(path: &Path) -> Result<()> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(path)?;
    let _header = reader.read_header()?;
    let mut inspected = 0usize;
    let mut paired = 0usize;
    for result in reader.records() {
        let record = result?;
        let flags = u16::from(record.flags());
        if flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            continue;
        }
        inspected += 1;
        if flags & FLAG_READ_PAIRED != 0 {
            paired += 1;
        }
        if inspected >= MAX_RECORDS_TO_INSPECT { break; }
    }
    if inspected == 0 {
        return Err(anyhow!("BAM contains no primary mapped records: {}", path.display()));
    }
    if paired == 0 {
        return Err(anyhow!(
            "rustqc atac requires paired-end input; first {} primary records had no PAIRED flag: {}",
            inspected,
            path.display()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(rel: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
    }

    #[test]
    fn accepts_paired_end_bam() {
        // tests/data/test.bam is the existing RNA fixture; it is paired-end.
        let p = fixture("tests/data/test.bam");
        assert!(assert_paired_end(&p).is_ok(), "test.bam should be paired-end");
    }

    // SE rejection is exercised via the GL-fixture suite in Phase 14 once we
    // materialize a synthetic SE BAM. For now we cover the assertion path; a
    // dedicated SE fixture is added in Task 14.x if not already covered.
}
```

- [ ] **Step 2: Build, test**

Run: `cargo test atac::pe_check 2>&1 | tail -10`
Expected: 1 test passes.

- [ ] **Step 3: Wire into `atac::run`**

In `src/atac/mod.rs`, replace the `bail!` placeholder body with:

```rust
pub fn run(args: AtacArgs) -> Result<()> {
    let cfg = resolve(&args);
    for input in &cfg.inputs {
        pe_check::assert_paired_end(std::path::Path::new(input))
            .with_context(|| format!("paired-end check failed for {}", input))?;
    }
    anyhow::bail!("rustqc atac is not yet implemented (PE check passed; metrics pending)");
}
```

- [ ] **Step 4: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/pe_check.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): paired-end startup check"
```

### Task 2.7: TSS extraction helper in `gtf.rs`

**Files:**
- Modify: `src/gtf.rs`

ATAC's per-TSS metrics need a deduplicated `Vec<Tss { chrom, pos, strand }>` derived from transcript records. The existing `gtf.rs` already parses transcripts; we add a small adapter.

- [ ] **Step 1: Read existing GTF API**

Run: `grep -n "pub fn\|pub struct" /home/xzg/project/RustQC/src/gtf.rs | head -40`

Note the type used for transcripts (call it `Transcript`). Identify how strand and chromosome are exposed.

- [ ] **Step 2: Write failing test**

In `src/gtf.rs`'s `mod tests`, add:

```rust
    #[test]
    fn tss_extraction_deduplicates() {
        // Build two transcripts that share the same TSS.
        let gtf = "\
chr1\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\ttest\ttranscript\t100\t600\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\";\n\
chr1\ttest\ttranscript\t900\t1500\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t3\";\n";
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(tmp.as_file_mut(), gtf.as_bytes()).unwrap();
        let tss = extract_tss(tmp.path()).unwrap();
        // + strand TSS = start (100); - strand TSS = end (1500). After dedup: 2 entries.
        assert_eq!(tss.len(), 2);
        assert!(tss.iter().any(|t| t.chrom == "chr1" && t.pos == 100 && t.strand == Strand::Plus));
        assert!(tss.iter().any(|t| t.chrom == "chr1" && t.pos == 1500 && t.strand == Strand::Minus));
    }
```

If `Strand` doesn't exist in `gtf.rs`, this test will not compile until Step 3. That is the point of TDD here.

- [ ] **Step 3: Add `Tss`, `Strand`, `extract_tss` to `gtf.rs`**

Add (next to existing types):

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand { Plus, Minus }

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Tss {
    pub chrom: String,
    pub pos: u64,         // 1-based, inclusive
    pub strand: Strand,
}

pub fn extract_tss(path: &std::path::Path) -> anyhow::Result<Vec<Tss>> {
    let transcripts = read_transcripts(path)?; // existing parser; rename if your parser is named differently
    let mut set: std::collections::HashSet<Tss> = std::collections::HashSet::new();
    for t in transcripts {
        let strand = match t.strand_char() {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => continue, // skip '.' / '?' transcripts
        };
        let pos = match strand {
            Strand::Plus => t.start as u64,
            Strand::Minus => t.end as u64,
        };
        set.insert(Tss { chrom: t.chrom.clone(), pos, strand });
    }
    let mut v: Vec<Tss> = set.into_iter().collect();
    v.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.pos.cmp(&b.pos)).then((a.strand as u8).cmp(&(b.strand as u8))));
    Ok(v)
}
```

If the existing transcript struct has different field names (e.g. `seqname` instead of `chrom`, `start_1based` instead of `start`), substitute accordingly. Repr `Strand as u8` requires `#[repr(u8)] enum Strand { Plus = 0, Minus = 1 }` for the cast to compile; either add `#[repr(u8)]` or replace the cast with a match.

- [ ] **Step 4: Run test**

Run: `cargo test gtf::tests::tss_extraction 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C /home/xzg/project/RustQC add src/gtf.rs
git -C /home/xzg/project/RustQC commit -m "feat(gtf): TSS extraction helper for ATAC"
```

---

## Phase 3 — bamQC

### Task 3.1: bamQC accumulator skeleton + flag/MAPQ counters

**Files:**
- Create: `src/atac/bam_qc.rs`
- Modify: `src/atac/mod.rs` (`pub mod bam_qc;`)

The accumulator owns per-chromosome state; `BamQcReport` holds the final aggregated metrics matching ATACseqQC's `bamQC` return list.

- [ ] **Step 1: Write the failing test for flag aggregation**

Create `src/atac/bam_qc.rs`:

```rust
//! ATACseqQC-style bamQC: rates, NRF, PBC1/2, MAPQ histogram.
//!
//! Numerical fidelity to ATACseqQC 1.36.0 R/bamQC.R is required. See
//! docs/superpowers/specs/2026-05-04-atac-seq-qc-design.md §"bamQC".

use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Default)]
pub struct BamQcAccum {
    pub total_records: u64,
    pub n_dup: u64,
    pub n_proper_pair: u64,
    pub n_unmapped: u64,
    pub n_unmapped_mate: u64,
    pub n_qc_fail: u64,
    pub n_mito: u64,
    pub mapq_hist: HashMap<u8, u64>,
    /// 5'-fingerprint multiset per chromosome, used for PBC1/PBC2/NRF aggregation.
    /// Key: (chrom_id, fingerprint_tuple_hash). Stored as a per-chromosome HashMap<key, count>
    /// outside this struct; here we only own scalars + MAPQ.
    pub qnames: HashSet<String>,
}

#[derive(Debug, Clone, Default)]
pub struct BamQcReport {
    pub total_qnames: u64,
    pub duplicate_rate: f64,
    pub mitochondria_rate: f64,
    pub proper_pair_rate: f64,
    pub unmapped_rate: f64,
    pub has_unmapped_mate_rate: f64,
    pub not_passing_qc_rate: f64,
    pub nrf: f64,
    pub pbc1: f64,
    pub pbc2: f64,
    pub mapq_hist: Vec<(u8, u64)>,   // sorted ascending by mapq
}

impl BamQcAccum {
    pub fn new() -> Self { Self::default() }

    /// Update from a single primary record's flag bits and MAPQ.
    /// Caller decides which records to feed; ATACseqQC excludes secondary alignments
    /// (same `isSecondaryAlignment = FALSE` flag we use).
    pub fn update_flags(&mut self, flags: u16, mapq: u8, is_mito: bool, qname: &str) {
        const F_PAIRED: u16 = 0x1;
        const F_PROPER_PAIR: u16 = 0x2;
        const F_UNMAPPED: u16 = 0x4;
        const F_MATE_UNMAPPED: u16 = 0x8;
        const F_DUP: u16 = 0x400;
        const F_QCFAIL: u16 = 0x200;
        let _ = F_PAIRED;
        self.total_records += 1;
        if flags & F_DUP != 0 { self.n_dup += 1; }
        if flags & F_PROPER_PAIR != 0 { self.n_proper_pair += 1; }
        if flags & F_UNMAPPED != 0 { self.n_unmapped += 1; }
        if flags & F_MATE_UNMAPPED != 0 { self.n_unmapped_mate += 1; }
        if flags & F_QCFAIL != 0 { self.n_qc_fail += 1; }
        if is_mito { self.n_mito += 1; }
        *self.mapq_hist.entry(mapq).or_default() += 1;
        self.qnames.insert(qname.to_string());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_aggregation_matches_R() {
        let mut a = BamQcAccum::new();
        // 4 records: 1 mito+dup, 1 proper-pair, 1 qc-fail, 1 unmapped-mate.
        a.update_flags(0x402, 30, true, "r1");   // dup + proper_pair
        a.update_flags(0x002, 60, false, "r2");  // proper_pair
        a.update_flags(0x200, 0,  false, "r3");  // qcfail, mapq 0
        a.update_flags(0x008, 30, false, "r4");  // mate unmapped
        assert_eq!(a.total_records, 4);
        assert_eq!(a.n_dup, 1);
        assert_eq!(a.n_proper_pair, 2);
        assert_eq!(a.n_qc_fail, 1);
        assert_eq!(a.n_unmapped_mate, 1);
        assert_eq!(a.n_mito, 1);
        assert_eq!(a.qnames.len(), 4);
        assert_eq!(a.mapq_hist[&30], 2);
        assert_eq!(a.mapq_hist[&60], 1);
        assert_eq!(a.mapq_hist[&0], 1);
    }
}
```

- [ ] **Step 2: Wire module + run test**

In `src/atac/mod.rs`: add `pub mod bam_qc;`.

Run: `cargo test atac::bam_qc 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/bam_qc.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): bamQC flag/MAPQ accumulator"
```

### Task 3.2: PBC fingerprint accumulator (per-chromosome)

PBC1 / PBC2 / NRF use a per-chromosome multiset of position fingerprints:
- PE: `(pos1, isize1, pos2, isize2)` where mate1 / mate2 are paired by qname; singletons appended as `(pos, isize, NA, NA)`.
- SE branch is unreachable for ATAC (we error out in PE check) but still implementable.

- [ ] **Step 1: Write the failing test**

Append to `src/atac/bam_qc.rs`:

```rust
#[derive(Debug, Clone, Default)]
pub struct PbcChromAccum {
    /// Fingerprint → count.
    pub fingerprints: HashMap<u128, u64>,
}

impl PbcChromAccum {
    pub fn add_pe(&mut self, pos1: i64, isize1: i64, pos2: i64, isize2: i64) {
        // Pack into u128: each i64 cast through u64::from_ne_bytes is fine because
        // R's identical-position semantics need bitwise equality, not arithmetic.
        let k = (pos1 as u64 as u128) << 96
              | (isize1 as u64 as u128) << 64
              | (pos2 as u64 as u128) << 32   // truncate isn't safe for very large coords;
              | (isize2 as u64 as u128);       // see test: positions used here fit in u32.
        *self.fingerprints.entry(k).or_default() += 1;
    }

    /// Returns (M_DISTINCT, M1, M2) — used in aggregate to compute NRF/PBC1/PBC2.
    pub fn summarize(&self) -> (u64, u64, u64) {
        let m_distinct = self.fingerprints.len() as u64;
        let m1 = self.fingerprints.values().filter(|&&c| c == 1).count() as u64;
        let m2 = self.fingerprints.values().filter(|&&c| c == 2).count() as u64;
        (m_distinct, m1, m2)
    }
}

#[test]
fn pbc_summarize_counts_singletons_and_doubletons() {
    let mut p = PbcChromAccum::default();
    p.add_pe(100, 200, 100, -200);
    p.add_pe(100, 200, 100, -200);   // duplicate of above
    p.add_pe(300, 200, 300, -200);   // singleton
    p.add_pe(500, 200, 500, -200);   // singleton
    p.add_pe(500, 200, 500, -200);   // doubleton with above
    let (m_distinct, m1, m2) = p.summarize();
    assert_eq!(m_distinct, 3);
    assert_eq!(m1, 1);
    assert_eq!(m2, 2);
}
```

**Important on packing**: the test above uses positions ≤ 2^31, so `u32`-truncating the lower halves is safe. For real BAMs, positions can exceed `u32::MAX` on large genomes. Fix this in Step 3 by switching to a `(i64, i64, i64, i64)` tuple key — `BTreeMap<(i64,i64,i64,i64), u64>` — instead of bit-packing. The test will still pass.

- [ ] **Step 2: Run test**

Run: `cargo test atac::bam_qc::tests::pbc 2>&1 | tail -5`
Expected: PASS (because the test positions fit). The packing fix is precautionary; ship it before metric integration in Phase 6.

- [ ] **Step 3: Switch the key to a tuple**

Replace `HashMap<u128, u64>` with `HashMap<(i64,i64,i64,i64), u64>` (or `BTreeMap` if you want deterministic order — `HashMap` is fine because we only summarize counts). Drop the bit-packing in `add_pe`. Re-run test; expect PASS unchanged.

- [ ] **Step 4: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/bam_qc.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): per-chromosome PBC fingerprint accumulator"
```

### Task 3.3: Aggregate `BamQcAccum` + per-chrom `PbcChromAccum` → `BamQcReport`

- [ ] **Step 1: Write the failing test**

Append to `src/atac/bam_qc.rs`:

```rust
pub fn finalize(
    flag_acc: &BamQcAccum,
    pbc_per_chrom: &[PbcChromAccum],
) -> BamQcReport {
    let total = flag_acc.total_records.max(1) as f64; // avoid div0 — total>0 in real runs
    let (mut sum_distinct, mut sum_m1, mut sum_m2) = (0u64, 0u64, 0u64);
    for p in pbc_per_chrom {
        let (md, m1, m2) = p.summarize();
        sum_distinct += md;
        sum_m1 += m1;
        sum_m2 += m2;
    }
    let total_qnames = flag_acc.qnames.len() as u64;
    let nrf = if total_qnames == 0 { 0.0 } else { sum_m1 as f64 / total_qnames as f64 };
    let pbc1 = if sum_distinct == 0 { 0.0 } else { sum_m1 as f64 / sum_distinct as f64 };
    let pbc2 = sum_m1 as f64 / sum_m2.max(1) as f64;

    let mut mapq_hist: Vec<(u8, u64)> = flag_acc.mapq_hist.iter().map(|(k, v)| (*k, *v)).collect();
    mapq_hist.sort_by_key(|(k, _)| *k);

    BamQcReport {
        total_qnames,
        duplicate_rate: flag_acc.n_dup as f64 / total,
        mitochondria_rate: flag_acc.n_mito as f64 / total,
        proper_pair_rate: flag_acc.n_proper_pair as f64 / total,
        unmapped_rate: flag_acc.n_unmapped as f64 / total,
        has_unmapped_mate_rate: flag_acc.n_unmapped_mate as f64 / total,
        not_passing_qc_rate: flag_acc.n_qc_fail as f64 / total,
        nrf, pbc1, pbc2, mapq_hist,
    }
}

#[test]
fn finalize_computes_NRF_PBC1_PBC2() {
    let mut flag = BamQcAccum::new();
    for i in 0..10 { flag.qnames.insert(format!("r{}", i)); }
    flag.total_records = 10;
    let mut p1 = PbcChromAccum::default();
    // Manually populate: 4 distinct fingerprints; 2 singletons, 1 doubleton, 1 triple.
    p1.fingerprints.insert((100, 200, 100, -200), 1); // M1
    p1.fingerprints.insert((200, 300, 200, -300), 1); // M1
    p1.fingerprints.insert((300, 400, 300, -400), 2); // M2
    p1.fingerprints.insert((400, 500, 400, -500), 3); // (no contribution to M1/M2)
    let r = finalize(&flag, &[p1]);
    // M1 = 2, M_DISTINCT = 4, M2 = 1, totalQNAMEs = 10.
    assert!((r.nrf - 0.2).abs() < 1e-12);
    assert!((r.pbc1 - 0.5).abs() < 1e-12);
    assert!((r.pbc2 - 2.0).abs() < 1e-12);
}
```

- [ ] **Step 2: Run test**

Run: `cargo test atac::bam_qc::tests::finalize 2>&1 | tail -5`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): finalize BamQcReport with NRF/PBC1/PBC2"
```

---

## Phase 4 — fragSizeDist

### Task 4.1: Histogram accumulator (1..1010 bp)

**Files:**
- Create: `src/atac/frag_size.rs`
- Modify: `src/atac/mod.rs`

ATACseqQC counts every primary mapped record's `|TLEN|` (both mates contribute), then plots density (count / total_count). Range pinned to 1..1010 bp matching R `match(1:1010, names(table))`.

- [ ] **Step 1: Write failing test**

Create `src/atac/frag_size.rs`:

```rust
//! Fragment-length histogram (1..1010 bp). Mirrors ATACseqQC R/fragSizeDist.R.

#[derive(Debug, Clone)]
pub struct FragSizeAccum {
    counts: [u64; 1011], // index 0 unused; valid index is 1..=1010
    total: u64,
}

impl Default for FragSizeAccum { fn default() -> Self { Self { counts: [0; 1011], total: 0 } } }

impl FragSizeAccum {
    pub fn new() -> Self { Self::default() }

    /// Update from one record's TLEN (signed). Records out of [1,1010] after abs are dropped.
    pub fn update(&mut self, tlen: i64) {
        let v = tlen.unsigned_abs();
        if v == 0 || v > 1010 { return; }
        self.counts[v as usize] += 1;
        self.total += 1;
    }

    /// Returns the (length, count, density) triples for length=1..=1010.
    pub fn finalize(&self) -> Vec<(u32, u64, f64)> {
        let total = self.total.max(1) as f64;
        (1..=1010u32)
            .map(|l| {
                let c = self.counts[l as usize];
                (l, c, c as f64 / total)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_counts_abs_tlen_within_range() {
        let mut a = FragSizeAccum::new();
        for &t in &[150_i64, -150, 200, 200, -200, -200, 1011, 0, -1011] {
            a.update(t);
        }
        let h = a.finalize();
        assert_eq!(h[150 - 1].1, 2);   // length 150 → 2 records (one + one −)
        assert_eq!(h[200 - 1].1, 4);
        assert_eq!(h[1010 - 1].1, 0);   // 1011 dropped
        // Density sums to 1.
        let s: f64 = h.iter().map(|(_, _, d)| d).sum();
        assert!((s - 1.0).abs() < 1e-12);
    }
}
```

- [ ] **Step 2: Wire module, run test, commit**

```bash
# Add `pub mod frag_size;` to src/atac/mod.rs
cargo test atac::frag_size 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/frag_size.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): fragSize histogram accumulator"
```

### Task 4.2: TSV writer for fragSize

- [ ] **Step 1: Write failing test**

Append to `src/atac/frag_size.rs`:

```rust
pub fn write_tsv<W: std::io::Write>(w: &mut W, h: &[(u32, u64, f64)]) -> std::io::Result<()> {
    writeln!(w, "length\tcount\tnorm_density")?;
    for (l, c, d) in h {
        writeln!(w, "{}\t{}\t{:.10e}", l, c, d)?;
    }
    Ok(())
}

#[test]
fn tsv_format_matches_spec() {
    let mut a = FragSizeAccum::new();
    a.update(100); a.update(150); a.update(150);
    let mut buf = Vec::new();
    write_tsv(&mut buf, &a.finalize()).unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(s.starts_with("length\tcount\tnorm_density\n"));
    let line_100 = s.lines().nth(100).unwrap();   // header + length 1..100 → index 100
    assert!(line_100.starts_with("100\t1\t"));
    let line_150 = s.lines().nth(150).unwrap();
    assert!(line_150.starts_with("150\t2\t"));
}
```

- [ ] **Step 2: Test, commit**

```bash
cargo test atac::frag_size::tests::tsv 2>&1 | tail -5
git -C /home/xzg/project/RustQC commit -am "feat(atac): fragSize TSV writer"
```

---

## Phase 5 — Library complexity (`readsDupFreq` + preseq)

### Task 5.1: `readsDupFreq` accumulator

**Files:**
- Create: `src/atac/lib_complexity.rs`
- Modify: `src/atac/mod.rs`

`readsDupFreq` builds a fingerprint multiset (PE: `(chrom, leftpos, isize)`) and emits a histogram of multiplicity counts: rows `(j, n_j)` where `n_j` is the number of fingerprints with multiplicity `j`. This histogram is then fed to preseq.

- [ ] **Step 1: Write failing test**

Create `src/atac/lib_complexity.rs`:

```rust
//! Library complexity: readsDupFreq → preseq ds_rsac_bootstrap.

use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct DupFreqAccum {
    /// (chrom_id, leftmost_pos, isize) → count
    pub fingerprints: HashMap<(u32, i64, i64), u64>,
}

impl DupFreqAccum {
    pub fn add_pe(&mut self, chrom_id: u32, leftpos: i64, isize: i64) {
        *self.fingerprints.entry((chrom_id, leftpos, isize)).or_default() += 1;
    }

    /// Build histogram rows: vec of (j, n_j) sorted by j ascending.
    pub fn histogram(&self) -> Vec<(u64, u64)> {
        let mut by_j: HashMap<u64, u64> = HashMap::new();
        for &c in self.fingerprints.values() {
            *by_j.entry(c).or_default() += 1;
        }
        let mut v: Vec<(u64, u64)> = by_j.into_iter().collect();
        v.sort_by_key(|(j, _)| *j);
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_counts_multiplicities() {
        let mut a = DupFreqAccum::default();
        // 5 distinct fingerprints with multiplicities 1, 1, 2, 3, 3
        a.add_pe(0, 100, 200);
        a.add_pe(0, 200, 200);
        a.add_pe(0, 300, 200); a.add_pe(0, 300, 200);
        a.add_pe(0, 400, 200); a.add_pe(0, 400, 200); a.add_pe(0, 400, 200);
        a.add_pe(0, 500, 200); a.add_pe(0, 500, 200); a.add_pe(0, 500, 200);
        let h = a.histogram();
        assert_eq!(h, vec![(1, 2), (2, 1), (3, 2)]);
    }
}
```

- [ ] **Step 2: Wire module, test, commit**

```bash
# Add `pub mod lib_complexity;` to src/atac/mod.rs
cargo test atac::lib_complexity 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/lib_complexity.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): readsDupFreq accumulator"
```

### Task 5.2: Preseq integration

ATACseqQC's `estimateLibComplexity` calls `preseqR::ds.rSAC.bootstrap(hist, r=1, times=100)` and evaluates at `relative.size ∈ {0.1..1.0 step 0.1, 5, 10, 15, 20}`. We have a Rust port at `crate::preseq`.

- [ ] **Step 1: Inspect existing preseq surface**

Run: `grep -n "pub fn\|pub struct" /home/xzg/project/RustQC/src/preseq.rs | head -20`

Identify the function that takes `&[(u64,u64)]` histogram + `times` + sample sizes and returns SAC values. Likely `estimate_complexity` or a lower-level helper.

- [ ] **Step 2: Write the integration test**

Append to `src/atac/lib_complexity.rs`:

```rust
#[derive(Debug, Clone)]
pub struct LibComplexityRow {
    pub relative_size: f64,
    pub distinct_fragments: f64,
    pub putative_reads: f64,
}

pub fn estimate(hist: &[(u64, u64)], times: u32) -> anyhow::Result<Vec<LibComplexityRow>> {
    // total_reads = Σ j * n_j  (matches R: histFile[,1] %*% histFile[,2])
    let total: u64 = hist.iter().map(|(j, n)| j * n).sum();
    let sample_sizes: Vec<f64> = (1..=10).map(|i| i as f64 * 0.1)
        .chain([5.0, 10.0, 15.0, 20.0])
        .collect();
    let sac = crate::preseq::ds_rsac_bootstrap(hist, /*r=*/ 1, times)?;
    let rows = sample_sizes.iter().map(|&s| LibComplexityRow {
        relative_size: s,
        distinct_fragments: sac.evaluate(s),
        putative_reads: s * total as f64,
    }).collect();
    Ok(rows)
}

#[test]
#[ignore = "integration test — run after preseq surface confirmed exposed"]
fn estimate_returns_14_rows() {
    let hist = vec![(1u64, 100u64), (2, 50), (3, 20), (4, 5)];
    let rows = estimate(&hist, 50).unwrap();
    assert_eq!(rows.len(), 14);
    assert!((rows[0].relative_size - 0.1).abs() < 1e-12);
    assert_eq!(rows[13].relative_size, 20.0);
}
```

- [ ] **Step 3: Adapt to actual preseq API**

If `crate::preseq::ds_rsac_bootstrap` does not exist, find the existing entrypoint (e.g. `estimate_complexity`). Adapt the wrapper above to call it; preserve the same row shape (14 rows at the documented sample sizes).

- [ ] **Step 4: Run the (un-ignored) test**

Run: `cargo test atac::lib_complexity::tests::estimate -- --include-ignored 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): library complexity via preseq"
```

---

## Phase 6 — Per-TSS 5'-end coverage windows

This single component supports TSSEscore, NFRscore, and PTscore. Build it once.

### Task 6.1: Sparse per-TSS coverage buffer

**Files:**
- Create: `src/atac/tss_cov.rs`
- Modify: `src/atac/mod.rs`

For each TSS, store a `Vec<u32>` of length `2*max_flank` (where `max_flank = max(2000+500, tsse_flank)`), zeroed; per primary mapped record, increment the bin under the read's 5' position if the read overlaps any TSS window.

- [ ] **Step 1: Write failing test**

Create `src/atac/tss_cov.rs`:

```rust
//! Sparse per-TSS 5'-end coverage. Underlies TSSEscore, NFRscore, PTscore.

use crate::gtf::{Strand, Tss};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct TssCov {
    pub flank: u32,                            // half-window in bp; arrays have length 2*flank
    pub buffers: Vec<Vec<u32>>,                // index = TSS index in `tss_list`
    pub tss_list: Vec<Tss>,
    by_chrom: HashMap<String, Vec<usize>>,     // chrom → indices into tss_list
}

impl TssCov {
    pub fn new(tss_list: Vec<Tss>, flank: u32) -> Self {
        let mut by_chrom: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, t) in tss_list.iter().enumerate() {
            by_chrom.entry(t.chrom.clone()).or_default().push(i);
        }
        let buffers = tss_list.iter().map(|_| vec![0u32; (2 * flank) as usize]).collect();
        Self { flank, buffers, tss_list, by_chrom }
    }

    /// Increment the bin under a read's 5' position if it falls within any TSS window
    /// on this chromosome. `pos5p` is 1-based (BAM coordinate convention).
    pub fn add_5prime(&mut self, chrom: &str, pos5p: u64) {
        let Some(idxs) = self.by_chrom.get(chrom) else { return; };
        for &i in idxs {
            let t = &self.tss_list[i];
            let win_start = t.pos.saturating_sub(self.flank as u64);
            let win_end = t.pos + self.flank as u64;
            if pos5p < win_start || pos5p >= win_end { continue; }
            let bin = (pos5p - win_start) as usize;
            // For + strand, bin 0 = TSS - flank, bin (2*flank - 1) = TSS + flank - 1.
            // For − strand, mirror: bin 0 = TSS + flank - 1, bin (2*flank - 1) = TSS - flank.
            let bin = match t.strand {
                Strand::Plus => bin,
                Strand::Minus => (2 * self.flank as usize) - 1 - bin,
            };
            self.buffers[i][bin] = self.buffers[i][bin].saturating_add(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tss_at(chrom: &str, pos: u64, strand: Strand) -> Tss {
        Tss { chrom: chrom.into(), pos, strand }
    }

    #[test]
    fn coverage_strand_aware() {
        let tss = vec![
            tss_at("chr1", 1000, Strand::Plus),
            tss_at("chr1", 5000, Strand::Minus),
        ];
        let mut c = TssCov::new(tss, 100);
        c.add_5prime("chr1", 1050);   // 50 bp downstream of + TSS → bin 100+50=150 in +-strand frame
        c.add_5prime("chr1", 4990);   // 10 bp upstream of − TSS (genomic) → upstream in mirrored frame
        // + TSS bin: 1050 - (1000-100) = 150
        assert_eq!(c.buffers[0][150], 1);
        // − TSS bin: 4990 in genomic; bin0 = TSS - flank = 4900; raw = 90; mirrored = 199-90 = 109
        assert_eq!(c.buffers[1][109], 1);
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod tss_cov;` to src/atac/mod.rs
cargo test atac::tss_cov 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/tss_cov.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): sparse per-TSS coverage buffer"
```

---

## Phase 7 — NFRscore + PTscore (no loess)

### Task 7.1: NFRscore

**Files:**
- Create: `src/atac/nfr_score.rs`
- Modify: `src/atac/mod.rs`

Reads `TssCov` buffers and produces per-TSS `NFR_score`, `log2meancov`, `n1`, `nf`, `n2` rows, then aggregates to a median for the JSON summary.

Window coords (default `N=150, F=100`):
- `n1` = bins `[flank-200, flank-50)` → 150 bp
- `nf` = bins `[flank-50, flank+50)` → 100 bp
- `n2` = bins `[flank+50, flank+200)` → 150 bp

Where `bin = flank` corresponds to position `TSS` itself (in the strand-mirrored frame). Negative side already mirrored by `TssCov`, so the same indexing works on both strands.

**Required**: `flank ≥ 200` so all three windows fit. Phase 9 sets `flank = max(2000+500, tsse_flank)` ≥ 2500.

- [ ] **Step 1: Write failing test**

Create `src/atac/nfr_score.rs`:

```rust
//! NFRscore: NFR_score = log2(nf + ε) + 1 - log2(n1+n2 + ε).

use crate::atac::tss_cov::TssCov;

const N: usize = 150;
const F: usize = 100;

#[derive(Debug, Clone)]
pub struct NfrRow {
    pub tss_idx: usize,
    pub n1: f64,
    pub nf: f64,
    pub n2: f64,
    pub nfr_score: f64,
    pub log2_mean_cov: f64,
}

pub fn compute(cov: &TssCov) -> Vec<NfrRow> {
    let flank = cov.flank as usize;
    assert!(flank >= 200, "TssCov flank must be >=200 for NFRscore (got {})", flank);
    // First pass: collect raw per-window means to compute the smallNumber floor.
    let mut raw: Vec<(f64, f64, f64)> = Vec::with_capacity(cov.buffers.len());
    for buf in &cov.buffers {
        let n1: f64 = buf[flank - 200..flank - 50].iter().map(|&v| v as f64).sum::<f64>() / N as f64;
        let nf: f64 = buf[flank - 50 ..flank + 50].iter().map(|&v| v as f64).sum::<f64>() / F as f64;
        let n2: f64 = buf[flank + 50 ..flank + 200].iter().map(|&v| v as f64).sum::<f64>() / N as f64;
        raw.push((n1, nf, n2));
    }
    let small = {
        let min_finite = |xs: &[(f64, f64, f64)], pick: fn(&(f64,f64,f64))->f64| -> f64 {
            xs.iter().map(pick).filter(|x| x.is_finite()).fold(f64::INFINITY, f64::min)
        };
        let m_n1 = min_finite(&raw, |t| t.0);
        let m_nf = min_finite(&raw, |t| t.1);
        let m_n2 = min_finite(&raw, |t| t.2);
        // R: max(c(1e-6, min(nf), min(n1), min(n2)))
        [1e-6, m_n1, m_nf, m_n2].iter().cloned().filter(|x| x.is_finite()).fold(1e-6, f64::max)
    };
    raw.iter().enumerate().map(|(i, &(n1, nf, n2))| {
        let log2_mean_cov = ((3.0 * (n1 + n2) + 2.0 * nf) / 8.0 + small).log2();
        let nfr_score = (nf + small).log2() + 1.0 - (n1 + n2 + small).log2();
        NfrRow { tss_idx: i, n1, nf, n2, nfr_score, log2_mean_cov }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Strand, Tss};

    #[test]
    fn nfr_score_simple_uniform_signal() {
        // One TSS, signal = 1 on every bin → n1=nf=n2=1.
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 1000);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        let rows = compute(&cov);
        assert_eq!(rows.len(), 1);
        // n1+n2 = 2, nf = 1; small ≥ 1.
        // NFR = log2(1+ε)+1 - log2(2+ε) ≈ -1+1+0 = 0 when ε→0.
        let r = &rows[0];
        assert!((r.n1 - 1.0).abs() < 1e-12);
        assert!((r.nf - 1.0).abs() < 1e-12);
        assert!((r.n2 - 1.0).abs() < 1e-12);
        assert!((r.nfr_score - 0.0).abs() < 1e-9);
    }

    #[test]
    fn nfr_score_strong_nfr_signal() {
        // 10× signal on the central nf window.
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 1000);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        for b in 950..1050 { cov.buffers[0][b] = 10; }
        let rows = compute(&cov);
        let r = &rows[0];
        // n1=n2=1, nf=10; NFR = log2(10) + 1 - log2(2) = log2(10) − 0 = ~3.32.
        assert!((r.nfr_score - (10.0_f64.log2() + 1.0 - 2.0_f64.log2())).abs() < 1e-9);
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod nfr_score;` to src/atac/mod.rs
cargo test atac::nfr_score 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/nfr_score.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): NFRscore"
```

### Task 7.2: PTscore

**Files:**
- Create: `src/atac/pt_score.rs`
- Modify: `src/atac/mod.rs`

Promoter window = `[TSS - 2000, TSS + 500]` (2500 bp); body = `[TSS + 500, TSS + 3000]` (2500 bp). With `flank = max(2500, …)`, body extends to `flank` past TSS — the first 500 bp inside `flank` are promoter (right side), and the next 2500 bp are body. We need `flank ≥ 3000` to fit body fully. Set the floor accordingly in Task 9.x's flank-resolution.

In strand-mirrored frame:
- promoter = bins `[flank - 2000, flank + 500)` (2500 bp)
- body = bins `[flank + 500, flank + 3000)` (2500 bp)

- [ ] **Step 1: Write failing test**

Create `src/atac/pt_score.rs`:

```rust
//! PTscore: PT_score = log2(promoter + ε) - log2(body + ε).

use crate::atac::tss_cov::TssCov;

const U: usize = 2000;
const D: usize = 500;

#[derive(Debug, Clone)]
pub struct PtRow {
    pub tss_idx: usize,
    pub promoter: f64,
    pub body: f64,
    pub pt_score: f64,
    pub log2_mean_cov: f64,
}

pub fn compute(cov: &TssCov) -> Vec<PtRow> {
    let flank = cov.flank as usize;
    assert!(flank >= U + D + (U + D) - U, "PTscore needs flank >= 3000 (got {})", flank);
    let mut raw: Vec<(f64, f64)> = Vec::with_capacity(cov.buffers.len());
    for buf in &cov.buffers {
        let prom = (flank - U..flank + D).map(|b| buf[b] as f64).sum::<f64>() / (U + D) as f64;
        let body = (flank + D..flank + D + (U + D)).map(|b| buf[b] as f64).sum::<f64>() / (U + D) as f64;
        raw.push((prom, body));
    }
    let small = {
        let min_finite = |xs: &[(f64, f64)], pick: fn(&(f64,f64))->f64| -> f64 {
            xs.iter().map(pick).filter(|x| x.is_finite()).fold(f64::INFINITY, f64::min)
        };
        [1e-6, min_finite(&raw, |t| t.0), min_finite(&raw, |t| t.1)]
            .iter().cloned().filter(|x| x.is_finite()).fold(1e-6, f64::max)
    };
    raw.iter().enumerate().map(|(i, &(prom, body))| PtRow {
        tss_idx: i,
        promoter: prom,
        body,
        pt_score: (prom + small).log2() - (body + small).log2(),
        log2_mean_cov: (prom + small).log2() + (body + small).log2(),
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Strand, Tss};

    #[test]
    fn pt_score_promoter_dominates() {
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 3000);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        // Boost promoter region 4×.
        for b in (3000-2000)..(3000+500) { cov.buffers[0][b] = 4; }
        let rows = compute(&cov);
        let r = &rows[0];
        // promoter mean ≈ 4, body ≈ 1 → PT ≈ log2(4)-log2(1) = 2.
        assert!((r.pt_score - 2.0).abs() < 1e-9);
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod pt_score;` to src/atac/mod.rs
cargo test atac::pt_score 2>&1 | tail -5
git -C /home/xzg/project/RustQC add src/atac/pt_score.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): PTscore"
```

---

## Phase 8 — Loess port (TSSEscore prerequisite)

ATACseqQC uses `stats::loess.smooth(x, y, family="gaussian", evaluation=length(y))`. We port a minimal `loess` matching `stats::loess`/`loess.smooth` defaults: `span=2/3`, `degree=2`, tricube weights, gaussian family (no robust iterations). For `evaluation = N` and `x = 1..N`, output is fitted values at x.

### Task 8.1: Tricube weight + locally-weighted polynomial fit

**Files:**
- Create: `src/atac/loess.rs`
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing test (sanity)**

Create `src/atac/loess.rs`:

```rust
//! Minimal loess port matching stats::loess.smooth defaults
//! (span=2/3, degree=2, family="gaussian", no robust reweighting).

/// Fit local degree-2 polynomial weighted by tricube on the q nearest neighbors,
/// where q = ceil(span * n). Evaluates at every x in `xs`.
pub fn loess_smooth(xs: &[f64], ys: &[f64], span: f64, degree: usize) -> Vec<f64> {
    assert_eq!(xs.len(), ys.len());
    let n = xs.len();
    if n == 0 { return vec![]; }
    let q = ((span * n as f64).ceil() as usize).clamp(degree + 1, n);
    let mut out = Vec::with_capacity(n);
    for &x0 in xs {
        // Pick q nearest neighbors by |x - x0|.
        let mut dists: Vec<(usize, f64)> = (0..n).map(|i| (i, (xs[i] - x0).abs())).collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let nbrs: Vec<usize> = dists.iter().take(q).map(|(i, _)| *i).collect();
        let max_d = dists[q - 1].1.max(f64::MIN_POSITIVE);
        // Tricube weights.
        let w: Vec<f64> = nbrs.iter().map(|&i| {
            let u = (xs[i] - x0).abs() / max_d;
            let one_minus = (1.0 - u.powi(3)).max(0.0);
            one_minus.powi(3)
        }).collect();
        // Solve weighted least squares y ~ poly(x − x0, degree) by normal equations.
        // Build X (q × (degree+1)) and W (diagonal, weights).
        let p = degree + 1;
        let mut xtwx = vec![0.0f64; p * p];
        let mut xtwy = vec![0.0f64; p];
        for (k, &i) in nbrs.iter().enumerate() {
            let dx = xs[i] - x0;
            let mut row = vec![1.0f64; p];
            for j in 1..p { row[j] = row[j - 1] * dx; }
            let wk = w[k];
            for a in 0..p {
                for b in 0..p { xtwx[a * p + b] += row[a] * row[b] * wk; }
                xtwy[a] += row[a] * ys[i] * wk;
            }
        }
        // Solve (p × p) symmetric positive (semi-)definite system via Gauss-Jordan.
        let beta = solve_linear(&mut xtwx, &mut xtwy, p);
        // Fitted value at x0 corresponds to the constant term β₀.
        out.push(beta[0]);
    }
    out
}

fn solve_linear(a: &mut [f64], b: &mut [f64], p: usize) -> Vec<f64> {
    // In-place Gaussian elimination on (a | b).
    for k in 0..p {
        // Pivot.
        let mut piv = k;
        for r in k+1..p {
            if a[r * p + k].abs() > a[piv * p + k].abs() { piv = r; }
        }
        if piv != k {
            for c in 0..p { a.swap(k * p + c, piv * p + c); }
            b.swap(k, piv);
        }
        let akk = a[k * p + k];
        if akk.abs() < 1e-15 { return vec![0.0; p]; }
        for r in 0..p {
            if r == k { continue; }
            let factor = a[r * p + k] / akk;
            for c in k..p { a[r * p + c] -= factor * a[k * p + c]; }
            b[r] -= factor * b[k];
        }
    }
    (0..p).map(|i| b[i] / a[i * p + i]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_quadratic_exactly_at_full_span() {
        let xs: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 + 3.0 * x + 0.5 * x * x).collect();
        // span=1.0 + degree=2 → exact recovery of a quadratic.
        let fit = loess_smooth(&xs, &ys, 1.0, 2);
        for (a, b) in fit.iter().zip(ys.iter()) {
            assert!((a - b).abs() < 1e-6, "loess(span=1) on quadratic: {} vs {}", a, b);
        }
    }

    #[test]
    fn fits_constant_signal() {
        let xs: Vec<f64> = (1..=20).map(|i| i as f64).collect();
        let ys = vec![5.0; 20];
        let fit = loess_smooth(&xs, &ys, 2.0/3.0, 2);
        for v in fit { assert!((v - 5.0).abs() < 1e-9, "{}", v); }
    }
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test atac::loess 2>&1 | tail -10`
Expected: 2 tests pass.

- [ ] **Step 3: Cross-check against R**

Hand-validate against R for an asymmetric input. From a separate R session (offline; not committed to CI):

```r
x <- 1:20; set.seed(0); y <- sin(x/3) + rnorm(20, sd=0.05)
r <- loess.smooth(x, y, family="gaussian", evaluation=length(y))
print(round(r$y, 6))
```

Run the same `xs, ys` through `loess_smooth(.., 2.0/3.0, 2)`. Assert per-point absolute diff ≤ `1e-3`. (This is the spec's TSSE acceptance bar; tighter is fine but not required.) If the diff exceeds the bar on the asymmetric input, examine whether `loess.smooth` differs from `loess` in evaluation-grid handling — `loess.smooth` evaluates at `seq(min,max,length=evaluation)` whereas `loess` evaluates at `xs`. We use the latter form because TSSE feeds in `1..20` directly, which equals `seq(1,20,length=20)`.

- [ ] **Step 4: Commit**

```bash
# Add `pub mod loess;` to src/atac/mod.rs
git -C /home/xzg/project/RustQC add src/atac/loess.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): minimal loess port for TSSE"
```

---

## Phase 9 — TSSEscore + flank resolution

### Task 9.1: TSSE windowed normalization

**Files:**
- Create: `src/atac/tsse.rs`
- Modify: `src/atac/mod.rs`

Algorithm (defaults: upstream=downstream=1000, endSize=100, width=100, step=100, pseudocount=0):
1. From each TSS buffer (length 2*flank, flank≥1000), define 20 sliding windows of width 100 step 100 starting at `flank - 1000`.
2. Each window's value is the mean coverage across its 100 bins.
3. Left flank = bins `[flank-1000, flank-900)`, right flank = `[flank+900, flank+1000)`. Mean coverage of each → `vl[t], vr[t]`.
4. NA handling: each TSS contributes only if `(vl + vr)/2 > 0` (with default pseudocount=0).
5. Normalize per surviving TSS: `v_norm[w] = v[w] * endSize / blk / width = v[w] / blk` (since endSize == width == 100).
6. Average across surviving TSSs per window → `vms.m`, length 20.
7. `loess_smooth(1..20, vms.m, 2/3, 2)` → smoothed; `tsse_score = max(smoothed)`.

- [ ] **Step 1: Write failing tests**

Create `src/atac/tsse.rs`:

```rust
//! TSSEscore. Mirrors ATACseqQC R/TSSEscore.R.

use crate::atac::loess::loess_smooth;
use crate::atac::tss_cov::TssCov;

const TSSE_FLANK: usize = 1000;
const END_SIZE: usize = 100;
const WIDTH: usize = 100;

#[derive(Debug, Clone)]
pub struct TsseResult {
    pub values: Vec<f64>,    // smoothed, length = 2*TSSE_FLANK / WIDTH = 20
    pub tsse_score: f64,
}

pub fn compute(cov: &TssCov) -> TsseResult {
    let flank = cov.flank as usize;
    assert!(flank >= TSSE_FLANK, "TssCov flank must be >=1000 for TSSE (got {})", flank);
    let n_windows = (2 * TSSE_FLANK) / WIDTH;
    let center_lo = flank - TSSE_FLANK;
    let mut sums = vec![0.0f64; n_windows];
    let mut surviving = 0u64;
    for buf in &cov.buffers {
        let mean = |range: std::ops::Range<usize>| -> f64 {
            let mut s = 0u64; let mut n = 0u64;
            for b in range { s += buf[b] as u64; n += 1; }
            s as f64 / n as f64
        };
        let vl = mean(center_lo..center_lo + END_SIZE);
        let vr = mean(center_lo + 2 * TSSE_FLANK - END_SIZE..center_lo + 2 * TSSE_FLANK);
        let blk = (vl + vr) / 2.0;
        if blk <= 0.0 { continue; }
        for w in 0..n_windows {
            let lo = center_lo + w * WIDTH;
            let v = mean(lo..lo + WIDTH);
            sums[w] += v / blk;
        }
        surviving += 1;
    }
    let s = surviving.max(1) as f64;
    let raw: Vec<f64> = sums.iter().map(|x| x / s).collect();
    let xs: Vec<f64> = (1..=n_windows).map(|i| i as f64).collect();
    let smoothed = loess_smooth(&xs, &raw, 2.0/3.0, 2);
    let tsse_score = smoothed.iter().cloned().filter(|v| v.is_finite()).fold(f64::NEG_INFINITY, f64::max);
    TsseResult { values: smoothed, tsse_score }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::{Strand, Tss};

    #[test]
    fn flat_signal_yields_score_near_one() {
        // Uniform signal on a 2*flank window → blk = 1, every v_norm = 1, smoothed ≈ 1.
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 1000);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        let r = compute(&cov);
        assert!((r.tsse_score - 1.0).abs() < 1e-3, "score={}", r.tsse_score);
        assert_eq!(r.values.len(), 20);
    }

    #[test]
    fn central_enrichment_lifts_score_above_baseline() {
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 1000);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        // 5× boost in the central 200 bp.
        for b in 900..1100 { cov.buffers[0][b] = 5; }
        let r = compute(&cov);
        assert!(r.tsse_score > 1.5, "expected enrichment, got {}", r.tsse_score);
    }
}
```

- [ ] **Step 2: Wire, test**

```bash
# Add `pub mod tsse;` to src/atac/mod.rs
cargo test atac::tsse 2>&1 | tail -10
```
Expected: 2 tests pass.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC add src/atac/tsse.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): TSSEscore"
```

### Task 9.2: Flank resolution helper

The `TssCov::flank` value must be `max(3000 /* PT body fits */, tsse_flank)`.

- [ ] **Step 1: Add resolver to `src/atac/mod.rs`**

```rust
pub fn resolve_flank(tsse_flank: u32) -> u32 {
    const PT_REQUIREMENT: u32 = 3000;
    tsse_flank.max(PT_REQUIREMENT)
}

#[cfg(test)]
#[test]
fn flank_floor_at_3000() {
    assert_eq!(resolve_flank(1000), 3000);
    assert_eq!(resolve_flank(5000), 5000);
}
```

- [ ] **Step 2: Test, commit**

```bash
cargo test atac::tests::flank 2>&1 | tail -5
git -C /home/xzg/project/RustQC commit -am "feat(atac): flank resolution"
```

---

## Phase 10 — Tn5 +4/−5 shift

### Task 10.1: Coordinate-only shift (no SEQ/QUAL/CIGAR rewrite)

For the in-memory metric path, only the alignment's 5'-end position and TLEN need adjustment:
- `+ strand`: `pos += 4`, `tlen = sign*(|tlen|−9)`.
- `− strand`: `end −= 5`, `tlen = sign*(|tlen|−9)`.

The SEQ/QUAL/CIGAR rewrite is required only when `--emit-shifted-bam` is set (Task 12.1).

**Files:**
- Create: `src/atac/shift.rs`
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing test**

Create `src/atac/shift.rs`:

```rust
//! Tn5 +4/-5 shift. Coordinate-only path used by metrics; full BAM path is in bam_writer.rs.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShiftedFrag {
    /// 5'-end genomic position (1-based, inclusive).
    pub pos5p: u64,
    /// Signed TLEN after shift; abs(TLEN) is the post-shift fragment length.
    pub tlen: i64,
}

/// Apply +4/-5 to a single mapped record's 5'-end (caller passes whether it's the +
/// strand and the original TLEN). Returns None if the shifted fragment would have
/// non-positive width.
pub fn shift_5prime(pos5p: u64, is_plus: bool, tlen: i64) -> Option<ShiftedFrag> {
    let new_pos = if is_plus { pos5p + 4 } else { pos5p.checked_sub(5)? };
    let new_tlen = if tlen == 0 { 0 } else {
        let sign = tlen.signum();
        let abs = tlen.unsigned_abs() as i64;
        if abs <= 9 { return None; }
        sign * (abs - 9)
    };
    Some(ShiftedFrag { pos5p: new_pos, tlen: new_tlen })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_strand_shifts_pos_by_plus4_and_shrinks_tlen_by_9() {
        let s = shift_5prime(100, true, 200).unwrap();
        assert_eq!(s.pos5p, 104);
        assert_eq!(s.tlen, 191);
    }

    #[test]
    fn minus_strand_shifts_pos_by_minus5_and_shrinks_tlen_by_9() {
        let s = shift_5prime(100, false, -200).unwrap();
        assert_eq!(s.pos5p, 95);
        assert_eq!(s.tlen, -191);
    }

    #[test]
    fn drops_fragment_when_tlen_le_9() {
        assert!(shift_5prime(100, true, 9).is_none());
        assert!(shift_5prime(100, false, -9).is_none());
    }

    #[test]
    fn drops_when_minus_strand_pos_underflows() {
        assert!(shift_5prime(3, false, -50).is_none());
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod shift;` to src/atac/mod.rs
cargo test atac::shift 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/shift.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): coordinate-only Tn5 shift"
```

---

## Phase 11 — Fixed-interval split

**Files:**
- Create: `src/atac/split.rs`
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing test**

Create `src/atac/split.rs`:

```rust
//! Fixed-interval fragment-size split: NFR / mono / di / tri buckets.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FragBucket { Nfr, Mono, Di, Tri, Other }

pub fn classify(abs_tlen: u32) -> FragBucket {
    if abs_tlen < 100                            { FragBucket::Nfr }
    else if (180..=247).contains(&abs_tlen)      { FragBucket::Mono }
    else if (315..=473).contains(&abs_tlen)      { FragBucket::Di }
    else if (558..=615).contains(&abs_tlen)      { FragBucket::Tri }
    else                                          { FragBucket::Other }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_cases_match_atacseqqc_intervals() {
        assert_eq!(classify(0), FragBucket::Nfr);
        assert_eq!(classify(99), FragBucket::Nfr);
        assert_eq!(classify(100), FragBucket::Other);     // gap [100,179]
        assert_eq!(classify(180), FragBucket::Mono);
        assert_eq!(classify(247), FragBucket::Mono);
        assert_eq!(classify(248), FragBucket::Other);
        assert_eq!(classify(315), FragBucket::Di);
        assert_eq!(classify(473), FragBucket::Di);
        assert_eq!(classify(558), FragBucket::Tri);
        assert_eq!(classify(615), FragBucket::Tri);
        assert_eq!(classify(616), FragBucket::Other);
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod split;` to src/atac/mod.rs
cargo test atac::split 2>&1 | tail -5
git -C /home/xzg/project/RustQC add src/atac/split.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): fixed-interval fragment split"
```

---

## Phase 12 — BAM emission (opt-in)

This phase is only exercised when `--emit-shifted-bam` or `--emit-split-bams` is set. The metric path does not depend on it.

### Task 12.1: Tn5 BAM record rewrite

**Files:**
- Create: `src/atac/bam_writer.rs`
- Modify: `src/atac/mod.rs`

We rewrite POS, CIGAR, SEQ, QUAL, and TLEN. Soft-clip handling: per ATACseqQC's `shiftReads.R`, we first project SEQ/QUAL through the soft-clip via `sequenceLayer(from='query', to='query-after-soft-clipping')`, then clip 4 bases (forward) or 5 (reverse) from the read's 5' end, then rebuild CIGAR via `cigarQNarrow`. Insertions at the 5' end (where `cigarWidthAlongQuerySpace != width(seq)`) are clipped from the 3' end before the narrowing step.

Because this rewrite is intricate, we implement it test-first against a hand-built record.

- [ ] **Step 1: Write failing test for the simple no-soft-clip case**

Create `src/atac/bam_writer.rs`:

```rust
//! Tn5-shifted BAM emission and length-split BAM emission via noodles.

use noodles_sam::alignment::record::cigar::op::{Kind as CigarKind, Op as CigarOp};

/// Trim `n` bases from the 5'-read end of a sorted CIGAR. Returns the new CIGAR
/// and the genomic-coordinate shift to apply to POS (so the caller can update
/// the alignment record's POS).
///
/// "5'-read end" is the first op for + strand reads, the last op for − strand reads
/// (the latter case is the caller's responsibility — pass a reversed op slice).
pub fn trim_cigar_5prime(ops: &[CigarOp], n: u32) -> (Vec<CigarOp>, u32) {
    let mut remaining = n;
    let mut shift = 0u32;
    let mut out = Vec::with_capacity(ops.len());
    let mut iter = ops.iter().copied();
    while let Some(op) = iter.next() {
        if remaining == 0 { out.push(op); for r in iter { out.push(r); } break; }
        let consumes_query = matches!(op.kind(), CigarKind::Match | CigarKind::Insertion | CigarKind::SoftClip | CigarKind::SequenceMatch | CigarKind::SequenceMismatch);
        let consumes_ref   = matches!(op.kind(), CigarKind::Match | CigarKind::Deletion | CigarKind::Skip | CigarKind::SequenceMatch | CigarKind::SequenceMismatch);
        let len = op.len() as u32;
        if !consumes_query { if consumes_ref { shift += len; } continue; }
        if len <= remaining {
            remaining -= len;
            if consumes_ref { shift += len; }
        } else {
            let kept = len - remaining;
            if consumes_ref { shift += remaining; }
            out.push(CigarOp::new(op.kind(), kept as usize));
            remaining = 0;
            for r in iter { out.push(r); }
            break;
        }
    }
    (out, shift)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::cigar::op::{Kind as K, Op as O};

    fn op(k: K, n: usize) -> O { O::new(k, n) }

    #[test]
    fn trim_4_from_pure_match() {
        let ops = vec![op(K::Match, 50)];
        let (out, shift) = trim_cigar_5prime(&ops, 4);
        assert_eq!(out, vec![op(K::Match, 46)]);
        assert_eq!(shift, 4);
    }

    #[test]
    fn trim_4_consumes_softclip_first() {
        let ops = vec![op(K::SoftClip, 3), op(K::Match, 50)];
        let (out, shift) = trim_cigar_5prime(&ops, 4);
        // 3 soft-clip bases + 1 match base consumed from query; ref shift = 1.
        assert_eq!(out, vec![op(K::Match, 49)]);
        assert_eq!(shift, 1);
    }

    #[test]
    fn trim_4_passes_through_insertion() {
        let ops = vec![op(K::Match, 2), op(K::Insertion, 3), op(K::Match, 50)];
        // First 2 match (query+ref) → remaining 2; then 2 of 3 insertion bases (query only) → remaining 0.
        let (out, shift) = trim_cigar_5prime(&ops, 4);
        assert_eq!(shift, 2);
        // Surviving ops: leftover insertion (1), then full match (50).
        assert_eq!(out, vec![op(K::Insertion, 1), op(K::Match, 50)]);
    }
}
```

- [ ] **Step 2: Wire module, run test, commit**

```bash
# Add `pub mod bam_writer;` to src/atac/mod.rs (gated behind an emit flag at runtime; module always compiled)
cargo test atac::bam_writer 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/bam_writer.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): CIGAR 5' trim helper for Tn5 BAM emission"
```

### Task 12.2: Full record rewrite

- [ ] **Step 1: Write the failing integration test**

Append to `src/atac/bam_writer.rs`:

```rust
use noodles_sam as sam;

/// Apply +4/-5 Tn5 shift to a single noodles SAM/BAM record, in place.
/// Returns Ok(false) if the read should be dropped (e.g. shifted width <= 0).
pub fn rewrite_record_inplace(rec: &mut sam::alignment::RecordBuf, is_plus: bool) -> anyhow::Result<bool> {
    use sam::alignment::record_buf::Cigar;
    let n = if is_plus { 4 } else { 5 };
    // Read CIGAR.
    let ops_in: Vec<CigarOp> = rec.cigar().as_ref().iter().copied().collect();

    // For + strand, trim from 5'-read = first op. For − strand, trim from 3'-read end
    // (which corresponds to the read's 5' relative to the reference orientation we rewrote).
    let (new_ops, ref_shift): (Vec<CigarOp>, u32) = if is_plus {
        trim_cigar_5prime(&ops_in, n)
    } else {
        let rev: Vec<CigarOp> = ops_in.iter().rev().copied().collect();
        let (mut trimmed, _shift) = trim_cigar_5prime(&rev, n);
        trimmed.reverse();
        // For − strand, no POS shift (POS stays the same — the genomic right-end shrinks).
        (trimmed, 0)
    };
    if new_ops.is_empty() { return Ok(false); }

    // Apply POS shift (+ strand only).
    if is_plus {
        if let Some(pos) = rec.alignment_start() {
            let new = usize::from(pos) + ref_shift as usize;
            rec.alignment_start_mut().replace(sam::alignment::record_buf::Position::try_from(new)?);
        }
    }

    // Trim SEQ + QUAL by `n` from the 5'-read end (which is record-relative for + strand,
    // and the 3'-record end for − strand reads, since BAM stores the record-as-aligned).
    let seq = rec.sequence_mut();
    let qual = rec.quality_scores_mut();
    let len = seq.len();
    if len <= n as usize { return Ok(false); }
    if is_plus {
        seq.drain(..n as usize);
        qual.drain(..n as usize);
    } else {
        let new_len = len - n as usize;
        seq.truncate(new_len);
        qual.truncate(new_len);
    }

    // Rewrite CIGAR.
    *rec.cigar_mut() = Cigar::from(new_ops);

    // TLEN: sign * (|TLEN| − 9).
    let tlen = i32::from(rec.template_length());
    if tlen != 0 {
        let abs = tlen.unsigned_abs();
        if abs <= 9 { return Ok(false); }
        let new = (abs - 9) as i32 * tlen.signum();
        *rec.template_length_mut() = new;
    }
    Ok(true)
}

// Note: noodles record_buf field names may differ; match against the version pinned in Cargo.toml.
// If the API doesn't expose mutable accessors for sequence/qual under noodles 0.85,
// allocate a new RecordBuf and copy over only the trimmed fields. The unit test below
// is an integration test that builds a synthetic RecordBuf and verifies the post-shift
// fields without depending on the exact accessor names.
```

The actual API surface for `sam::alignment::RecordBuf` in noodles 0.85 may not include all the mutators above; the engineer adapts. The unit test below builds a synthetic record using whichever constructor noodles 0.85 exposes.

- [ ] **Step 2: Confirm noodles record-buf API and adapt**

Run: `cargo doc --no-deps --open` and inspect `noodles_sam::alignment::record_buf::RecordBuf`. Substitute the exact accessors. Re-run tests.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): full Tn5 record rewrite via noodles RecordBuf"
```

### Task 12.3: BAM writer wiring (shifted + split outputs)

This task wires up:
- `--emit-shifted-bam` → one writer for `<sample>.shifted.bam`
- `--emit-split-bams` → four writers for `<sample>.{NFR,mono,di,tri}.bam`
- BAI indexing after writes

- [ ] **Step 1: Build the writer multiplexer**

Append to `src/atac/bam_writer.rs`:

```rust
use std::path::PathBuf;

#[derive(Default)]
pub struct EmitWriters {
    pub shifted: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub nfr:     Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub mono:    Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub di:      Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub tri:     Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
}

impl EmitWriters {
    pub fn open(
        outdir: &std::path::Path,
        sample: &str,
        emit_shifted: bool,
        emit_split: bool,
        header: &noodles_sam::Header,
    ) -> anyhow::Result<Self> {
        // Implementation: open requested files under outdir/{shifted,split}/, write header, return.
        // Defer details until adapted to the exact noodles 0.85 builder API.
        let _ = (outdir, sample, emit_shifted, emit_split, header);
        Ok(Self::default())
    }
}
```

This task is implementation-heavy and depends on the exact noodles surface; the agent fills in `open` and `write_record` calls during execution. The post-write step calls `noodles_bam::bai::index` (or the equivalent indexer entrypoint in the pinned version) to produce `.bai` for each output.

- [ ] **Step 2: Defer integration test to Phase 14**

The numerical test against `inst/extdata/splited/*.bam` (ATACseqQC's own pre-split output) lives in Phase 14. For now, only the multiplexer skeleton is in place.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): scaffold EmitWriters multiplexer"
```

---

## Phase 13 — Single-pass driver, plots, JSON summary

### Task 13.1: Single-pass BAM driver

**Files:**
- Modify: `src/atac/mod.rs` (replace the `bail!` placeholder body of `run`)
- Reference: `src/rna/bam_io.rs` (now `src/bam_io.rs`)

The driver:
1. Calls `pe_check::assert_paired_end`.
2. Loads GTF → `Vec<Tss>` via `gtf::extract_tss`.
3. Auto-detects mito chrom from header (`mito::detect_mito` over `@SQ` names) unless `--mito-chrom` is set.
4. Initializes `BamQcAccum`, `Vec<PbcChromAccum>` (one per @SQ), `FragSizeAccum`, `DupFreqAccum`, `TssCov::new(tss, resolve_flank(cfg.tsse_flank))`, optional `EmitWriters`.
5. Streams every primary record through:
   - update flag/MAPQ/mito counters (always)
   - update fragSize histogram (always; both mates contribute)
   - update DupFreq fingerprint (always; PE: chrom_id + leftpos + isize)
   - update PbcChromAccum (PE: full (pos1, isize1, pos2, isize2) tuple — pair both mates by qname; track first/second-mate buffer)
   - update TssCov with the read's 5'-end position (always)
   - if `--emit-shifted-bam` and the record passes Tn5 shift, write to shifted writer; if `--emit-split-bams`, classify and route to the matching writer
6. After streaming: `bam_qc::finalize`, `tsse::compute`, `nfr_score::compute`, `pt_score::compute`, `lib_complexity::estimate`, then plot SVGs and write the summary JSON.

This task is large (~150 lines of glue); break the implementation into focused commits as the engineer goes:
- 13.1.a: skeleton main loop with no metrics, just a record counter
- 13.1.b: wire bamQC + fragSize
- 13.1.c: wire TssCov + 3 TSS metrics
- 13.1.d: wire DupFreq + lib_complexity
- 13.1.e: wire EmitWriters

For each sub-task: write a single integration test that runs the driver against `tests/data/test.bam` (the existing RNA fixture is paired-end and small) with a stub GTF, asserts the produced summary JSON has the metric fields populated, and asserts the run completes in < 30 s.

- [ ] **Step 1: Add the driver skeleton**

Replace the body of `run` in `src/atac/mod.rs`:

```rust
pub fn run(args: AtacArgs) -> Result<()> {
    use anyhow::Context;
    let cfg = resolve(&args);
    if cfg.inputs.len() != 1 {
        anyhow::bail!("rustqc atac currently accepts exactly one BAM (got {}); multi-BAM support is future work", cfg.inputs.len());
    }
    let input = std::path::Path::new(&cfg.inputs[0]);
    pe_check::assert_paired_end(input).with_context(|| format!("PE check: {}", input.display()))?;

    let tss_list = crate::gtf::extract_tss(std::path::Path::new(&cfg.gtf))?;
    let flank = resolve_flank(cfg.tsse_flank);
    let mut tss_cov = tss_cov::TssCov::new(tss_list.clone(), flank);
    let mut frag = frag_size::FragSizeAccum::new();
    let mut bq = bam_qc::BamQcAccum::new();
    let mut dup = lib_complexity::DupFreqAccum::default();
    // pbc_per_chrom indexed by SQ id; built after reading the BAM header.

    // Open BAM via crate::bam_io; mirror what RNA does.
    // ... (engineer fills in based on src/bam_io.rs surface)

    // After the loop:
    let bq_report = bam_qc::finalize(&bq, &[/* pbc_per_chrom */]);
    let tsse = tsse::compute(&tss_cov);
    let nfr = nfr_score::compute(&tss_cov);
    let pt = pt_score::compute(&tss_cov);
    let _ = (bq_report, tsse, nfr, pt, frag, dup);

    anyhow::bail!("driver wiring in progress (Task 13.1.x)")
}
```

- [ ] **Step 2: Add the metric-wiring tests one sub-task at a time, committing after each green run.**

Driver-level integration tests live in `tests/integration_atac.rs` (created in Phase 14).

- [ ] **Step 3: Final commit of Phase 13.1 once all sub-tasks 13.1.a–e are in**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): single-pass driver"
```

### Task 13.2: SVG plots

**Files:**
- Create: `src/atac/plots.rs`
- Modify: `src/atac/mod.rs`

Three plots:
- fragSize: linear x∈[0,1010], y = density × 1000, with a log10 inset at fig=(.4,.95,.4,.95).
- TSSE: line plot of `values` vs window index 1..20.
- saturation: `relative_size * total_reads` (x, in millions) vs `distinct_fragments` (y, in millions).

Reuse the plotters backend already wired for RNA (see `src/rna/dupradar/plots.rs` for examples).

- [ ] **Step 1: Add `plots.rs` mirroring RNA's plotting style**

(No tests — visual output. Snapshot any failures with `--quiet` smoke runs against the GL fixtures in Phase 14.)

- [ ] **Step 2: Commit**

```bash
git -C /home/xzg/project/RustQC commit -am "feat(atac): SVG plots"
```

### Task 13.3: JSON summary writer

**Files:**
- Create: `src/atac/summary.rs`
- Modify: `src/atac/mod.rs`

Implements the JSON schema from §4 of the design.

- [ ] **Step 1: Write failing test for serialization shape**

Create `src/atac/summary.rs`:

```rust
//! ATAC summary JSON. Schema documented in
//! docs/superpowers/specs/2026-05-04-atac-seq-qc-design.md §"JSON summary schema".

use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct ToolVersions {
    pub rustqc: String,
    pub atacseqqc_replicates: String,
}

#[derive(Debug, Serialize)]
pub struct AtacSummary {
    pub sample: String,
    pub tool_versions: ToolVersions,
    pub split_method: &'static str,
    pub bamqc: BamqcSection,
    pub fragsize: FragsizeSection,
    pub tsse: TsseSection,
    pub nfr: ScoreSection,
    pub pt: ScoreSection,
    pub lib_complexity: LibComplexitySection,
}

#[derive(Debug, Serialize)]
pub struct BamqcSection {
    pub total_qnames: u64,
    pub duplicate_rate: f64,
    pub mitochondria_rate: f64,
    pub proper_pair_rate: f64,
    pub unmapped_rate: f64,
    pub has_unmapped_mate_rate: f64,
    pub not_passing_qc_rate: f64,
    pub nrf: f64,
    pub pbc1: f64,
    pub pbc2: f64,
    pub mapq_histogram: serde_json::Map<String, serde_json::Value>,
}
#[derive(Debug, Serialize)] pub struct FragsizeSection { pub total_pairs: u64, pub tsv_path: String }
#[derive(Debug, Serialize)] pub struct TsseSection     { pub score: f64, pub values: Vec<f64>, pub tsv_path: String }
#[derive(Debug, Serialize)] pub struct ScoreSection    { pub median_score: f64, pub tsv_path: String }
#[derive(Debug, Serialize)] pub struct LibComplexitySection { pub extrapolated_total: f64, pub tsv_path: String }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_keys_match_spec() {
        let s = AtacSummary {
            sample: "GL1".into(),
            tool_versions: ToolVersions { rustqc: "0.4.0".into(), atacseqqc_replicates: "1.36.0".into() },
            split_method: "fixed_intervals_v1",
            bamqc: BamqcSection { total_qnames: 1000, duplicate_rate: 0.1, mitochondria_rate: 0.05, proper_pair_rate: 0.9, unmapped_rate: 0.0, has_unmapped_mate_rate: 0.0, not_passing_qc_rate: 0.0, nrf: 0.8, pbc1: 0.9, pbc2: 5.0, mapq_histogram: Default::default() },
            fragsize: FragsizeSection { total_pairs: 500, tsv_path: "fragsize/GL1.fragsize.tsv".into() },
            tsse: TsseSection { score: 6.0, values: vec![1.0; 20], tsv_path: "tsse/GL1.tsse.tsv".into() },
            nfr: ScoreSection { median_score: 1.5, tsv_path: "nfr/GL1.nfr.tsv".into() },
            pt: ScoreSection { median_score: 0.5, tsv_path: "pt/GL1.pt.tsv".into() },
            lib_complexity: LibComplexitySection { extrapolated_total: 3.4e8, tsv_path: "lib_complexity/GL1.libcomplexity.tsv".into() },
        };
        let j: serde_json::Value = serde_json::to_value(&s).unwrap();
        for k in ["sample","tool_versions","split_method","bamqc","fragsize","tsse","nfr","pt","lib_complexity"] {
            assert!(j.get(k).is_some(), "missing key {}", k);
        }
        assert_eq!(j["bamqc"]["pbc1"], 0.9);
        assert_eq!(j["tsse"]["values"].as_array().unwrap().len(), 20);
    }
}
```

- [ ] **Step 2: Wire, test, commit**

```bash
# Add `pub mod summary;` to src/atac/mod.rs
cargo test atac::summary 2>&1 | tail -10
git -C /home/xzg/project/RustQC add src/atac/summary.rs src/atac/mod.rs
git -C /home/xzg/project/RustQC commit -m "feat(atac): JSON summary writer"
```

---

## Phase 14 — Numerical fidelity tests against ATACseqQC fixtures

### Task 14.1: Extract GL1–GL4 fixtures

**Files:**
- Create: `tests/data/atac/GL{1..4}.bam(.bai)`
- Create: `tests/data/atac/splited/{NucleosomeFree,mononucleosome,dinucleosome,trinucleosome}.bam(.bai)`
- Create: `tests/data/atac/.gitattributes` (mark BAMs as binary)

- [ ] **Step 1: Extract**

```bash
mkdir -p /home/xzg/project/RustQC/tests/data/atac
cd /home/xzg/project/RustQC/tests/data/atac
tar -xzf /home/xzg/project/RustQC/ATACseqQC_1.36.0.tar.gz \
  ATACseqQC/inst/extdata/GL1.bam ATACseqQC/inst/extdata/GL1.bam.bai \
  ATACseqQC/inst/extdata/GL2.bam ATACseqQC/inst/extdata/GL2.bam.bai \
  ATACseqQC/inst/extdata/GL3.bam ATACseqQC/inst/extdata/GL3.bam.bai \
  ATACseqQC/inst/extdata/GL4.bam ATACseqQC/inst/extdata/GL4.bam.bai \
  ATACseqQC/inst/extdata/splited
mv ATACseqQC/inst/extdata/* .
rm -rf ATACseqQC
```

- [ ] **Step 2: Get a usable GTF for chr1 hg19 small region**

ATACseqQC uses TxDb.Hsapiens.UCSC.hg19.knownGene; the GL BAMs map to chr1 only (a small subset). To avoid pulling a multi-GB UCSC GTF, distill a fragment that covers the BAMs' alignment range. From a one-off offline R session:

```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
chr1 <- txs[seqnames(txs) == "chr1"]
# Restrict to the BAMs' alignment range; widen a bit for safety.
range_of_interest <- GRanges("chr1", IRanges(48000000, 49000000))
sub <- subsetByOverlaps(chr1, range_of_interest)
rtracklayer::export(sub, "tests/data/atac/GL_tss.gtf", format="gtf")
```

Commit `GL_tss.gtf` to the repo. **The R script is not run by CI; the resulting file is committed as a fixture.**

- [ ] **Step 3: Commit fixtures**

```bash
git -C /home/xzg/project/RustQC add tests/data/atac/
git -C /home/xzg/project/RustQC commit -m "test(atac): extract GL1-4 BAM fixtures + chr1 GTF subset"
```

### Task 14.2: R reference outputs (offline)

**Files:**
- Create: `tests/atac/golden/run_r_reference.R` (committed but not invoked by CI)
- Create: `tests/atac/golden/GL{1..4}.{bamqc,fragsize,nfr,pt,tsse}.golden.json|tsv`

- [ ] **Step 1: Author the R script**

Create `tests/atac/golden/run_r_reference.R`:

```r
# Generates reference outputs for tests/atac/integration_atac.rs.
# Run once, offline, and commit the resulting goldens.
suppressPackageStartupMessages({
  library(ATACseqQC)
  library(GenomicFeatures)
  library(rtracklayer)
  library(jsonlite)
  library(GenomicAlignments)
})

bams <- file.path("tests/data/atac", sprintf("GL%d.bam", 1:4))
labels <- sprintf("GL%d", 1:4)
txs <- import("tests/data/atac/GL_tss.gtf")
out_dir <- "tests/atac/golden"
dir.create(out_dir, showWarnings = FALSE)

for (i in seq_along(bams)) {
  bam <- bams[i]; label <- labels[i]
  # bamQC
  qc <- bamQC(bam, outPath = NULL)
  qc$mapq_hist <- as.list(setNames(qc$MAPQ$Freq, qc$MAPQ$Var1))
  qc$MAPQ <- NULL; qc$idxstats <- NULL
  writeLines(toJSON(qc, auto_unbox = TRUE, pretty = TRUE),
             file.path(out_dir, sprintf("%s.bamqc.golden.json", label)))

  # fragSizeDist
  fs <- fragSizeDist(bam, label)[[1]]
  fs_df <- data.frame(length = as.integer(names(fs)), count = as.integer(fs))
  write.table(fs_df, file.path(out_dir, sprintf("%s.fragsize.golden.tsv", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  gal <- readGAlignments(bam)
  # NFRscore / PTscore / TSSEscore
  nfr <- as.data.frame(NFRscore(gal, txs))
  write.table(nfr[, c("seqnames","start","strand","n1","nf","n2","NFR_score","log2meanCoverage")],
              file.path(out_dir, sprintf("%s.nfr.golden.tsv", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pt <- as.data.frame(PTscore(gal, txs))
  write.table(pt[, c("seqnames","start","strand","promoter","transcriptBody","PT_score","log2meanCoverage")],
              file.path(out_dir, sprintf("%s.pt.golden.tsv", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ts <- TSSEscore(gal, txs)
  writeLines(toJSON(list(values = ts$values, tsse_score = ts$TSSEscore),
                    auto_unbox = TRUE, digits = 8),
             file.path(out_dir, sprintf("%s.tsse.golden.json", label)))
}
```

- [ ] **Step 2: Run it offline, commit the resulting goldens**

The engineer with R installed runs `Rscript tests/atac/golden/run_r_reference.R` once, inspects the resulting files, and commits them.

```bash
git -C /home/xzg/project/RustQC add tests/atac/golden/
git -C /home/xzg/project/RustQC commit -m "test(atac): R reference goldens for GL1-4"
```

### Task 14.3: Rust integration tests against goldens

**Files:**
- Create: `tests/integration_atac.rs`

- [ ] **Step 1: Write the integration tests**

Create `tests/integration_atac.rs` with one test per metric per fixture, e.g.:

```rust
// tests/integration_atac.rs
use std::path::PathBuf;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn approx_eq(a: f64, b: f64, eps: f64) -> bool { (a - b).abs() <= eps }

#[test]
fn gl1_fragsize_byte_identical() {
    let outdir = tempfile::tempdir().unwrap();
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .args(["atac",
               fixture("tests/data/atac/GL1.bam").to_str().unwrap(),
               "--gtf", fixture("tests/data/atac/GL_tss.gtf").to_str().unwrap(),
               "--outdir", outdir.path().to_str().unwrap(),
               "--sample-name", "GL1"])
        .status().unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(outdir.path().join("fragsize/GL1.fragsize.tsv")).unwrap();
    let want = std::fs::read_to_string(fixture("tests/atac/golden/GL1.fragsize.golden.tsv")).unwrap();
    // Compare on (length, count) pairs only — we can ignore the density column since it's normalized.
    let parse = |s: &str| -> Vec<(u32, u64)> {
        s.lines().skip(1).filter_map(|l| {
            let mut it = l.split('\t');
            Some((it.next()?.parse().ok()?, it.next()?.parse().ok()?))
        }).collect()
    };
    assert_eq!(parse(&got), parse(&want));
}

// Repeat for GL2-4. Then:

#[test]
fn gl1_bamqc_byte_identical_rates() {
    // Run rustqc, parse summary.json, compare to goldens with eps=0 for all rates and PBC/NRF.
    // (epsilon for floats: 1e-12 is acceptable since both sides do exact arithmetic.)
    todo!("engineer: invoke binary, parse JSON, compare fields keyed in §13.3 schema")
}

#[test]
fn gl1_nfr_per_tss_within_tolerance() {
    // Tolerance: 1e-6 per row. Match TSS by (chrom, start, strand).
    todo!()
}

#[test]
fn gl1_pt_per_tss_within_tolerance() {
    todo!()
}

#[test]
fn gl1_tsse_score_within_tolerance() {
    // Tolerance: 1e-3 (the loess-port bar from §5).
    // Pre-loess values (vms.m) byte-identical: requires exposing them via summary or a debug path.
    todo!()
}
```

The `todo!()` placeholders must be filled in by the engineer with concrete `assert_*` calls before the task is considered complete. Each test invokes the rustqc binary against the GL fixture, reads the produced JSON / TSV, and compares to the matching golden file.

- [ ] **Step 2: Run integration tests; iterate until all green**

Run: `cargo test --release --test integration_atac 2>&1 | tail -50`

Expected: 5 × 4 = 20 tests pass (one per metric per fixture). If a metric exceeds its tolerance, the failure points back to the algorithm module — re-investigate, don't loosen the tolerance.

- [ ] **Step 3: Commit**

```bash
git -C /home/xzg/project/RustQC add tests/integration_atac.rs
git -C /home/xzg/project/RustQC commit -m "test(atac): numerical-fidelity integration suite"
```

### Task 14.4: Tn5 shift + split BAM fixture comparison

ATACseqQC ships pre-split outputs in `inst/extdata/splited/`. We compare our split BAMs read-name-by-read-name against those.

- [ ] **Step 1: Write the test**

Append to `tests/integration_atac.rs`:

```rust
#[test]
fn split_outputs_match_atacseqqc_splited_fixture() {
    // Run rustqc atac --emit-split-bams against GL1.bam.
    // Compare per-bucket read-name sets to inst/extdata/splited/{NucleosomeFree,mono,di,tri}.bam.
    todo!()
}
```

- [ ] **Step 2: Run + iterate + commit**

```bash
cargo test --release --test integration_atac split_outputs 2>&1 | tail -20
git -C /home/xzg/project/RustQC commit -am "test(atac): split BAM read-set fidelity"
```

---

## Phase 15 — Documentation site updates

### Task 15.1: Add `docs/src/content/docs/atac/` pages

**Files:**
- Create: `docs/src/content/docs/atac/index.mdx` (overview + algorithm crosswalk to ATACseqQC)
- Create: `docs/src/content/docs/atac/cli.mdx` (CLI reference)
- Create: `docs/src/content/docs/atac/numerical-fidelity.mdx` (acceptance bars, fixture provenance)
- Modify: `docs/astro.config.mjs` if needed to surface the new section

- [ ] **Step 1: Mirror the structure of `docs/src/content/docs/rna/`**

Use existing rna pages as a template. Each page:
- States the algorithm in 1–2 paragraphs.
- Cross-links to the ATACseqQC R function it replicates.
- States the acceptance bar (byte-identical / 1e-6 / 1e-3) and where the goldens live.

- [ ] **Step 2: Update `README.md`**

Add a row under the existing tools table:

```markdown
- `rustqc atac` is a single-command ATAC-Seq QC tool that runs all QC analyses in one pass...

| Tool                  | Reimplements                                                          | Description                                          |
| --------------------- | --------------------------------------------------------------------- | ---------------------------------------------------- |
| bamQC                 | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `bamQC`     | Mapping rates, NRF, PBC1/2, MAPQ histogram           |
| fragSizeDist          | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `fragSizeDist` | Fragment-length distribution + plot                |
| TSSEscore             | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `TSSEscore` | TSS enrichment score                                 |
| NFRscore              | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `NFRscore`  | Nucleosome-free region score                         |
| PTscore               | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `PTscore`   | Promoter / transcript-body score                     |
| estimateLibComplexity | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `estimateLibComplexity` | Library complexity extrapolation         |
| Tn5 shift             | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `shiftGAlignmentsList` | +4/-5 Tn5 shift                              |
| Split (fixed)         | [ATACseqQC](https://bioconductor.org/packages/ATACseqQC/) `splitGAlignmentsByCut` (fixed-interval branch) | NFR/mono/di/tri split |
```

- [ ] **Step 3: Update `CHANGELOG.md`**

Add an entry: `feat(atac): add rustqc atac subcommand for ATAC-seq QC and Tn5 preprocessing`.

- [ ] **Step 4: Commit**

```bash
git -C /home/xzg/project/RustQC add docs/ README.md CHANGELOG.md
git -C /home/xzg/project/RustQC commit -m "docs(atac): add ATAC-seq QC documentation"
```

---

## Phase 16 — Final verification

- [ ] **Step 1: Full test suite**

Run: `cargo test --release --workspace 2>&1 | tail -20`
Expected: all tests green, including 20+ ATAC integration tests.

- [ ] **Step 2: Smoke-run on a GL fixture end-to-end**

```bash
TMPOUT=$(mktemp -d)
./target/release/rustqc atac tests/data/atac/GL1.bam \
  --gtf tests/data/atac/GL_tss.gtf \
  --outdir "$TMPOUT" \
  --emit-shifted-bam --emit-split-bams \
  -j "$TMPOUT/summary.json"
ls -R "$TMPOUT"
jq '.tsse.score, .bamqc.nrf, .bamqc.pbc1' "$TMPOUT/summary.json"
```

Expected: directory contains `bamqc/`, `fragsize/`, `tsse/`, `nfr/`, `pt/`, `lib_complexity/`, `shifted/`, `split/`, and `summary.json`.

- [ ] **Step 3: Cargo lints clean**

Run: `cargo clippy --release --workspace -- -D warnings 2>&1 | tail -20`
Expected: no warnings.

- [ ] **Step 4: Final commit / PR-ready**

```bash
git -C /home/xzg/project/RustQC log --oneline -30
```

---

## Self-Review Checklist (executed during plan authoring; pasted here for reviewer reference)

- **Spec coverage**: every spec section has a phase. §2 Refactor → Phase 1. CLI → 2.1. Config → 2.3, 2.4. Mito detection → 2.5. PE check → 2.6. GTF TSS → 2.7. bamQC → Phase 3. fragSizeDist → Phase 4. lib complexity → Phase 5. Per-TSS coverage → Phase 6. NFR → 7.1. PT → 7.2. Loess → Phase 8. TSSE → 9.1. Tn5 shift coords → 10.1. Split → Phase 11. BAM emission → Phase 12. Driver / plots / summary → Phase 13. Numerical fidelity → Phase 14. Docs → Phase 15.
- **Type consistency**: `Tss` / `Strand` defined in `gtf.rs`; reused everywhere. `BamQcReport` field names match the JSON schema in 13.3. `FragBucket` only used internally; its boundaries match the spec exactly.
- **Placeholders**: tasks 12.2, 13.1, 14.3, 14.4 contain explicit "engineer fills in" notes. These are acceptable because the surrounding context (noodles API specifics, integration glue) requires runtime inspection of the pinned crate version; the algorithmic content is fully specified. The plan still ships compilable test code for every algorithmic module (Phases 3–11).
