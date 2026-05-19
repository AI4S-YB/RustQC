# ATAC Tn5 Shift CLI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add explicit `rustqc atac` controls for Tn5 shift semantics, including `--tn5-shift yes|no`, `--input-is-shifted`, safe conflict validation, basic-only QC output when shift is disabled for unshifted input, and JSON metadata.

**Architecture:** Model shift intent as a small `Tn5Shift` enum in `src/cli.rs`, resolve it with YAML-backed `AtacConfig` in `src/atac/mod.rs`, and expose derived behavior through `ResolvedAtacConfig` helper methods. The ATAC driver gates only TSS-dependent work and outputs; bamQC, fragment size, and library complexity remain unchanged. The summary schema keeps stable top-level keys by serializing skipped `tsse`, `nfr`, and `pt` sections as JSON `null`.

**Tech Stack:** Rust 2021, clap derive `ValueEnum`, serde/serde_yaml_ng/serde_json, anyhow, existing ATAC integration fixtures under `tests/data/atac`.

---

## File Structure

- Modify `src/cli.rs`
  - Add `Tn5Shift` enum.
  - Add `AtacArgs.tn5_shift: Option<Tn5Shift>`.
  - Add `AtacArgs.input_is_shifted: Option<bool>` with bare-flag and explicit-false support.
  - Extend CLI parser tests.
- Modify `src/config.rs`
  - Add YAML fields `AtacConfig.tn5_shift: Option<Tn5Shift>` and `AtacConfig.input_is_shifted: Option<bool>`.
  - Add config parsing tests.
- Modify `src/atac/mod.rs`
  - Add resolved fields and validation helpers.
  - Use the shift-state helpers during TSS coverage accumulation.
  - Skip TSSE/NFR/PT computation, directories, files, and plots when TSS-dependent metrics are disabled.
- Modify `src/atac/summary.rs`
  - Add `Tn5ShiftSection`.
  - Change `tsse`, `nfr`, and `pt` to `Option<...>`.
  - Add summary schema tests for both full and skipped metric states.
- Modify `tests/integration_atac.rs`
  - Add one integration test for basic-only output.
  - Add one integration test for the contradictory flag combination.
- Modify `docs/src/content/docs/atac/cli.mdx`
  - Document the new flags, YAML semantics, behavior matrix, and skipped-output behavior.

---

### Task 1: CLI and YAML Parsing

**Files:**
- Modify: `src/cli.rs`
- Modify: `src/config.rs`

- [ ] **Step 1: Write failing CLI parser tests**

Add these tests inside `#[cfg(test)] mod tests` in `src/cli.rs`, after `test_atac_emit_flags`:

```rust
#[test]
fn test_atac_tn5_shift_flags() {
    let cli = Cli::parse_from([
        "rustqc",
        "atac",
        "test.bam",
        "--gtf",
        "genes.gtf",
        "--tn5-shift",
        "no",
        "--input-is-shifted",
    ]);
    match cli.command {
        Commands::Atac(args) => {
            assert_eq!(args.tn5_shift, Some(Tn5Shift::No));
            assert_eq!(args.input_is_shifted, Some(true));
        }
        #[allow(unreachable_patterns)]
        _ => panic!("Expected Atac subcommand"),
    }
}

#[test]
fn test_atac_input_is_shifted_accepts_explicit_false() {
    let cli = Cli::parse_from([
        "rustqc",
        "atac",
        "test.bam",
        "--gtf",
        "genes.gtf",
        "--input-is-shifted=false",
    ]);
    match cli.command {
        Commands::Atac(args) => {
            assert_eq!(args.input_is_shifted, Some(false));
        }
        #[allow(unreachable_patterns)]
        _ => panic!("Expected Atac subcommand"),
    }
}
```

- [ ] **Step 2: Run CLI tests to verify they fail**

Run:

```bash
cargo test cli::tests::test_atac
```

Expected: compilation fails because `Tn5Shift`, `AtacArgs.tn5_shift`, and `AtacArgs.input_is_shifted` do not exist.

- [ ] **Step 3: Add `Tn5Shift` enum and CLI fields**

In `src/cli.rs`, after the `Strandedness` enum `impl std::fmt::Display for Strandedness`, add:

```rust
/// Whether RustQC should apply Tn5 +4/-5 shift for ATAC insertion-site metrics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, ValueEnum, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Tn5Shift {
    /// Apply Tn5 shift when the input is not already shifted.
    #[default]
    Yes,
    /// Do not apply Tn5 shift.
    No,
}

impl std::fmt::Display for Tn5Shift {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Tn5Shift::Yes => write!(f, "yes"),
            Tn5Shift::No => write!(f, "no"),
        }
    }
}
```

In `src/cli.rs`, inside `pub struct AtacArgs`, after `pub mito_chrom: Option<String>,`, add:

```rust
    /// Apply Tn5 +4/-5 shift for insertion-site metrics: yes or no [default: yes]
    #[arg(
        long = "tn5-shift",
        value_enum,
        value_name = "yes|no",
        env = "RUSTQC_TN5_SHIFT",
        help_heading = "ATAC-specific"
    )]
    pub tn5_shift: Option<Tn5Shift>,

    /// Input BAM coordinates are already Tn5 shifted
    #[arg(
        long,
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "true",
        value_parser = clap::value_parser!(bool),
        env = "RUSTQC_INPUT_IS_SHIFTED",
        help_heading = "ATAC-specific"
    )]
    pub input_is_shifted: Option<bool>,
```

In `test_atac_default_args`, add:

```rust
                assert_eq!(args.tn5_shift, None);
                assert_eq!(args.input_is_shifted, None);
```

- [ ] **Step 4: Write failing YAML parsing tests**

Add this test in `src/config.rs` inside `#[cfg(test)] mod tests`, after `test_empty_top_level_config`:

```rust
#[test]
fn test_atac_shift_config() {
    let yaml = r#"
atac:
  tn5_shift: no
  input_is_shifted: true
"#;
    let config: Config = serde_yaml_ng::from_str(yaml).unwrap();
    assert_eq!(config.atac.tn5_shift, Some(crate::cli::Tn5Shift::No));
    assert_eq!(config.atac.input_is_shifted, Some(true));
}
```

- [ ] **Step 5: Run config test to verify it fails**

Run:

```bash
cargo test config::tests::test_atac_shift_config
```

Expected: compilation fails because `AtacConfig` has no `tn5_shift` or `input_is_shifted` fields.

- [ ] **Step 6: Add YAML fields**

In `src/config.rs`, change the import at the top from:

```rust
use crate::cli::Strandedness;
```

to:

```rust
use crate::cli::{Strandedness, Tn5Shift};
```

In `src/config.rs`, update `AtacConfig` to:

```rust
/// ATAC-seq QC configuration (YAML-backed).
#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct AtacConfig {
    /// Mitochondrial chromosome name; auto-detected when None.
    pub mito_chrom: Option<String>,
    /// TSSEscore flank window in bp (default 1000).
    pub tsse_flank: Option<u32>,
    /// Whether RustQC should apply Tn5 +4/-5 shift for insertion-site metrics.
    pub tn5_shift: Option<Tn5Shift>,
    /// Whether input BAM coordinates are already Tn5 shifted.
    pub input_is_shifted: Option<bool>,
    /// Emit Tn5-shifted BAM.
    pub emit_shifted_bam: bool,
    /// Emit NFR/mono/di/tri BAMs.
    pub emit_split_bams: bool,
}
```

- [ ] **Step 7: Run parser tests**

Run:

```bash
cargo test cli::tests::test_atac
cargo test config::tests::test_atac_shift_config
```

Expected: all three tests pass.

- [ ] **Step 8: Commit**

Run:

```bash
git add src/cli.rs src/config.rs
git commit -m "feat(atac): parse Tn5 shift inputs"
```

---

### Task 2: Resolve and Validate Effective Shift State

**Files:**
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing unit tests for defaults, precedence, and conflict validation**

In `src/atac/mod.rs`, inside `#[cfg(test)] mod tests`, update existing `AtacConfig` literals to include the new fields:

```rust
            tn5_shift: None,
            input_is_shifted: None,
```

Add these tests after `resolve_cli_args_override_yaml`:

```rust
#[test]
fn resolve_shift_defaults_to_yes_unshifted() {
    let r = resolve(
        &parse(&["rustqc", "atac", "x.bam", "--gtf", "g.gtf"]),
        &AtacConfig::default(),
    );
    assert_eq!(r.tn5_shift, crate::cli::Tn5Shift::Yes);
    assert!(!r.input_is_shifted);
    assert!(r.apply_tn5_shift());
    assert!(r.tss_dependent_metrics_enabled());
    assert!(r.validate_shift_state().is_ok());
}

#[test]
fn resolve_shift_no_unshifted_disables_tss_metrics() {
    let r = resolve(
        &parse(&[
            "rustqc",
            "atac",
            "x.bam",
            "--gtf",
            "g.gtf",
            "--tn5-shift",
            "no",
        ]),
        &AtacConfig::default(),
    );
    assert_eq!(r.tn5_shift, crate::cli::Tn5Shift::No);
    assert!(!r.input_is_shifted);
    assert!(!r.apply_tn5_shift());
    assert!(!r.tss_dependent_metrics_enabled());
    assert!(r.validate_shift_state().is_ok());
}

#[test]
fn resolve_shift_no_input_shifted_keeps_tss_metrics() {
    let r = resolve(
        &parse(&[
            "rustqc",
            "atac",
            "x.bam",
            "--gtf",
            "g.gtf",
            "--tn5-shift",
            "no",
            "--input-is-shifted",
        ]),
        &AtacConfig::default(),
    );
    assert_eq!(r.tn5_shift, crate::cli::Tn5Shift::No);
    assert!(r.input_is_shifted);
    assert!(!r.apply_tn5_shift());
    assert!(r.tss_dependent_metrics_enabled());
    assert!(r.validate_shift_state().is_ok());
}

#[test]
fn validate_rejects_shift_yes_with_input_shifted() {
    let r = resolve(
        &parse(&[
            "rustqc",
            "atac",
            "x.bam",
            "--gtf",
            "g.gtf",
            "--tn5-shift",
            "yes",
            "--input-is-shifted",
        ]),
        &AtacConfig::default(),
    );
    let err = r.validate_shift_state().unwrap_err().to_string();
    assert!(err.contains("--tn5-shift yes"), "error was: {}", err);
    assert!(err.contains("--input-is-shifted"), "error was: {}", err);
}

#[test]
fn resolve_cli_shift_overrides_yaml() {
    let atac_cfg = AtacConfig {
        mito_chrom: None,
        tsse_flank: None,
        tn5_shift: Some(crate::cli::Tn5Shift::No),
        input_is_shifted: Some(false),
        emit_shifted_bam: false,
        emit_split_bams: false,
    };
    let r = resolve(
        &parse(&[
            "rustqc",
            "atac",
            "x.bam",
            "--gtf",
            "g.gtf",
            "--tn5-shift",
            "yes",
        ]),
        &atac_cfg,
    );
    assert_eq!(r.tn5_shift, crate::cli::Tn5Shift::Yes);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
cargo test atac::tests::resolve_shift
cargo test atac::tests::validate_rejects_shift_yes_with_input_shifted
```

Expected: compilation fails because resolved shift fields and helper methods do not exist.

- [ ] **Step 3: Add resolved fields and helper methods**

In `src/atac/mod.rs`, add this import near the existing crate imports:

```rust
use crate::cli::Tn5Shift;
```

In `ResolvedAtacConfig`, after `pub tsse_flank: u32,`, add:

```rust
    pub tn5_shift: Tn5Shift,
    pub input_is_shifted: bool,
```

After the `ResolvedAtacConfig` struct, add:

```rust
impl ResolvedAtacConfig {
    /// Return true when RustQC should apply Tn5 +4/-5 shift in memory.
    pub fn apply_tn5_shift(&self) -> bool {
        self.tn5_shift == Tn5Shift::Yes && !self.input_is_shifted
    }

    /// Return true when TSSE/NFR/PT can be computed from insertion-site coordinates.
    pub fn tss_dependent_metrics_enabled(&self) -> bool {
        self.apply_tn5_shift() || self.input_is_shifted
    }

    /// Reject contradictory Tn5 shift settings before any BAM work starts.
    pub fn validate_shift_state(&self) -> Result<()> {
        if self.tn5_shift == Tn5Shift::Yes && self.input_is_shifted {
            anyhow::bail!(
                "--tn5-shift yes cannot be combined with --input-is-shifted; \
                 use --tn5-shift no --input-is-shifted for already shifted input, \
                 or omit --input-is-shifted for ordinary unshifted ATAC BAMs"
            );
        }
        Ok(())
    }
}
```

In `resolve`, after `tsse_flank`, add:

```rust
        tn5_shift: args
            .tn5_shift
            .or(atac_cfg.tn5_shift)
            .unwrap_or(Tn5Shift::Yes),
        input_is_shifted: args
            .input_is_shifted
            .or(atac_cfg.input_is_shifted)
            .unwrap_or(false),
```

At the start of `run`, immediately after `let cfg = resolve(&args, &atac_cfg);`, add:

```rust
    cfg.validate_shift_state()?;
```

- [ ] **Step 4: Run resolve tests**

Run:

```bash
cargo test atac::tests::resolve_shift
cargo test atac::tests::validate_rejects_shift_yes_with_input_shifted
```

Expected: all listed tests pass.

- [ ] **Step 5: Commit**

Run:

```bash
git add src/atac/mod.rs
git commit -m "feat(atac): resolve Tn5 shift behavior"
```

---

### Task 3: JSON Summary Schema for Shift State and Skipped Metrics

**Files:**
- Modify: `src/atac/summary.rs`
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing summary schema tests**

In `src/atac/summary.rs`, inside `schema_keys_match_spec`, update the synthetic `AtacSummary` initializer to include:

```rust
            tn5_shift: Tn5ShiftSection {
                requested: true,
                input_is_shifted: false,
                applied: true,
                tss_dependent_metrics_enabled: true,
            },
```

and wrap the existing TSS sections:

```rust
            tsse: Some(TsseSection {
                score: 7.5,
                n_windows: 20,
                values: vec![1.0; 20],
                tsv_path: "tsse/test_sample.tsse.tsv".to_string(),
            }),
            nfr: Some(ScoreSection {
                n_tss: 10,
                median_score: 0.5,
                tsv_path: "nfr/test_sample.nfr.tsv".to_string(),
            }),
            pt: Some(ScoreSection {
                n_tss: 10,
                median_score: 1.2,
                tsv_path: "pt/test_sample.pt.tsv".to_string(),
            }),
```

Add `"tn5_shift"` to the top-level key loop, and change TSS assertions to:

```rust
        assert_eq!(j["tn5_shift"]["requested"].as_bool(), Some(true));
        assert_eq!(j["tn5_shift"]["input_is_shifted"].as_bool(), Some(false));
        assert_eq!(j["tn5_shift"]["applied"].as_bool(), Some(true));
        assert_eq!(
            j["tn5_shift"]["tss_dependent_metrics_enabled"].as_bool(),
            Some(true)
        );

        assert_eq!(j["tsse"]["values"].as_array().unwrap().len(), 20);
        assert!(j["tsse"].get("tsv_path").is_some(), "tsse missing tsv_path");
        assert!(j["nfr"].get("tsv_path").is_some(), "nfr missing tsv_path");
        assert!(j["pt"].get("tsv_path").is_some(), "pt missing tsv_path");
```

Add this second test after `schema_keys_match_spec`:

```rust
#[test]
fn schema_serializes_skipped_tss_metrics_as_null() {
    let mut mapq_hist = serde_json::Map::new();
    mapq_hist.insert("30".to_string(), serde_json::Value::Number(500u64.into()));

    let s = AtacSummary {
        schema_version: "1.0".to_string(),
        sample: "test_sample".to_string(),
        tool_versions: ToolVersions {
            rustqc: "0.3.0".to_string(),
            atacseqqc_replicates: "1.36.0".to_string(),
        },
        split_method: "fixed_intervals_v1",
        tn5_shift: Tn5ShiftSection {
            requested: false,
            input_is_shifted: false,
            applied: false,
            tss_dependent_metrics_enabled: false,
        },
        bamqc: BamqcSection {
            total_qnames: 1000,
            duplicate_rate: 0.05,
            mitochondria_rate: 0.02,
            proper_pair_rate: 0.95,
            unmapped_rate: 0.01,
            has_unmapped_mate_rate: 0.02,
            not_passing_qc_rate: 0.005,
            nrf: 0.8,
            pbc1: 0.9,
            pbc2: 3.5,
            mapq_histogram: mapq_hist,
        },
        fragsize: FragsizeSection {
            total_pairs: 800,
            tsv_path: "fragsize/test_sample.fragsize.tsv".to_string(),
        },
        tsse: None,
        nfr: None,
        pt: None,
        lib_complexity: LibComplexitySection {
            n_rows: 14,
            extrapolated_total: Some(750.0),
            tsv_path: "lib_complexity/test_sample.libcomplexity.tsv".to_string(),
        },
    };

    let json = serde_json::to_string(&s).unwrap();
    let v: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert!(v["tsse"].is_null());
    assert!(v["nfr"].is_null());
    assert!(v["pt"].is_null());
    assert_eq!(
        v["tn5_shift"]["tss_dependent_metrics_enabled"].as_bool(),
        Some(false)
    );
}
```

- [ ] **Step 2: Run summary tests to verify they fail**

Run:

```bash
cargo test atac::summary::tests::schema
```

Expected: compilation fails because `Tn5ShiftSection` does not exist and TSS fields are not optional.

- [ ] **Step 3: Update summary schema**

In `src/atac/summary.rs`, update `AtacSummary`:

```rust
/// Complete JSON summary for one ATAC-seq sample.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtacSummary {
    pub schema_version: String,
    pub sample: String,
    pub tool_versions: ToolVersions,
    /// Always `"fixed_intervals_v1"` per spec.
    pub split_method: &'static str,
    pub tn5_shift: Tn5ShiftSection,
    pub bamqc: BamqcSection,
    pub fragsize: FragsizeSection,
    pub tsse: Option<TsseSection>,
    pub nfr: Option<ScoreSection>,
    pub pt: Option<ScoreSection>,
    pub lib_complexity: LibComplexitySection,
}
```

After `ToolVersions`, add:

```rust
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tn5ShiftSection {
    pub requested: bool,
    pub input_is_shifted: bool,
    pub applied: bool,
    pub tss_dependent_metrics_enabled: bool,
}
```

- [ ] **Step 4: Update driver summary construction to compile**

In `src/atac/mod.rs`, in `summary::AtacSummary { ... }`, add after `split_method`:

```rust
        tn5_shift: summary::Tn5ShiftSection {
            requested: cfg.tn5_shift == Tn5Shift::Yes,
            input_is_shifted: cfg.input_is_shifted,
            applied: cfg.apply_tn5_shift(),
            tss_dependent_metrics_enabled: cfg.tss_dependent_metrics_enabled(),
        },
```

Wrap current TSS sections in `Some(...)`:

```rust
        tsse: Some(summary::TsseSection {
            score: tsse_result.tsse_score,
            n_windows: tsse_result.values.len() as u32,
            values: tsse_result.values.clone(),
            tsv_path: tsse_tsv_path,
        }),
        nfr: Some(summary::ScoreSection {
            n_tss: nfr_rows.len() as u32,
            median_score: nfr_median,
            tsv_path: nfr_tsv_path,
        }),
        pt: Some(summary::ScoreSection {
            n_tss: pt_rows.len() as u32,
            median_score: pt_median,
            tsv_path: pt_tsv_path,
        }),
```

- [ ] **Step 5: Run summary tests**

Run:

```bash
cargo test atac::summary::tests::schema
```

Expected: both tests pass.

- [ ] **Step 6: Commit**

Run:

```bash
git add src/atac/summary.rs src/atac/mod.rs
git commit -m "feat(atac): report Tn5 shift state in summary"
```

---

### Task 4: Gate TSS-Dependent Metric Computation and Output

**Files:**
- Modify: `src/atac/mod.rs`

- [ ] **Step 1: Write failing unit tests for coordinate helper**

In `src/atac/mod.rs`, before `pub fn run(args: AtacArgs) -> Result<()>`, add this test-only target by first adding tests below in the existing test module:

```rust
#[test]
fn tss_position_shift_helper_applies_plus4_minus5() {
    assert_eq!(tss_position_for_coverage(100, false, true), Some(104));
    assert_eq!(tss_position_for_coverage(100, true, true), Some(95));
}

#[test]
fn tss_position_shift_helper_can_leave_input_unchanged() {
    assert_eq!(tss_position_for_coverage(100, false, false), Some(100));
    assert_eq!(tss_position_for_coverage(100, true, false), Some(100));
}

#[test]
fn tss_position_shift_helper_drops_underflowing_reverse_shift() {
    assert_eq!(tss_position_for_coverage(3, true, true), None);
}
```

- [ ] **Step 2: Run helper tests to verify they fail**

Run:

```bash
cargo test atac::tests::tss_position_shift_helper
```

Expected: compilation fails because `tss_position_for_coverage` does not exist.

- [ ] **Step 3: Add coordinate helper**

In `src/atac/mod.rs`, after `metric_path`, add:

```rust
/// Return the 1-based 5' position to accumulate for TSS-dependent metrics.
fn tss_position_for_coverage(pos5p: u64, is_reverse: bool, apply_tn5_shift: bool) -> Option<u64> {
    if !apply_tn5_shift {
        return Some(pos5p);
    }
    if is_reverse {
        pos5p.checked_sub(5)
    } else {
        Some(pos5p + 4)
    }
}
```

- [ ] **Step 4: Refactor TSS loading and coverage accumulation**

In `run`, after sample logging, replace the unconditional TSS load block with:

```rust
    let tss_metrics_enabled = cfg.tss_dependent_metrics_enabled();
    let mut tss_cov = if tss_metrics_enabled {
        let tss_list = crate::gtf::extract_tss(Path::new(&cfg.gtf))
            .with_context(|| format!("failed to parse GTF: {}", cfg.gtf))?;
        if tss_list.is_empty() {
            eprintln!(
                "[rustqc atac] WARNING: no TSS entries extracted from GTF — TSS metrics will be empty"
            );
        } else {
            eprintln!(
                "[rustqc atac] loaded {} TSS entries from GTF",
                tss_list.len()
            );
        }
        let flank = resolve_flank(cfg.tsse_flank);
        Some(TssCov::new(tss_list, flank))
    } else {
        eprintln!(
            "[rustqc atac] TSS-dependent metrics disabled because --tn5-shift no was used without --input-is-shifted"
        );
        None
    };
```

Remove the old lines that resolved `flank` and initialized `let mut tss_cov = TssCov::new(...)`.

In the single-pass scan, replace the block:

```rust
            let pos5p_shifted = if is_reverse {
                pos5p.checked_sub(5)
            } else {
                Some(pos5p + 4)
            };
            if let Some(p) = pos5p_shifted {
                tss_cov.add_5prime(chrom_name, p);
            }
```

with:

```rust
            if let Some(tss_cov) = tss_cov.as_mut() {
                if let Some(p) = tss_position_for_coverage(pos5p, is_reverse, cfg.apply_tn5_shift())
                {
                    tss_cov.add_5prime(chrom_name, p);
                }
            }
```

- [ ] **Step 5: Refactor finalization and output gating**

Replace the current TSS finalization block:

```rust
    let tsse_result = tsse::compute(&tss_cov);
    let nfr_rows = nfr_score::compute(&tss_cov);
    let pt_rows = pt_score::compute(&tss_cov);
```

with:

```rust
    let tss_metric_results = tss_cov.as_ref().map(|cov| {
        let tsse_result = tsse::compute(cov);
        let nfr_rows = nfr_score::compute(cov);
        let pt_rows = pt_score::compute(cov);
        let nfr_median = median(nfr_rows.iter().map(|r| r.nfr_score).collect());
        let pt_median = median(pt_rows.iter().map(|r| r.pt_score).collect());
        (tsse_result, nfr_rows, pt_rows, nfr_median, pt_median)
    });
```

Remove the old standalone `nfr_median` and `pt_median` computation.

In output directory creation, replace:

```rust
    mk("tsse")?;
    mk("nfr")?;
    mk("pt")?;
```

with:

```rust
    if tss_metric_results.is_some() {
        mk("tsse")?;
        mk("nfr")?;
        mk("pt")?;
    }
```

Wrap the TSSE TSV writer in:

```rust
    if let Some((tsse_result, _, _, _, _)) = &tss_metric_results {
        let p = metric_path(outdir, flat, "tsse", &format!("{}.tsse.tsv", sample));
        let mut w =
            BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "window_idx\tnorm_signal")?;
        for (i, v) in tsse_result.values.iter().enumerate() {
            writeln!(w, "{}\t{:.8}", i + 1, v)?;
        }
    }
```

Wrap the NFR TSV writer in:

```rust
    if let Some((_, nfr_rows, _, _, _)) = &tss_metric_results {
        let p = metric_path(outdir, flat, "nfr", &format!("{}.nfr.tsv", sample));
        let mut w =
            BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "tss_idx\tn1\tnf\tn2\tnfr_score\tlog2meancov")?;
        for r in nfr_rows {
            writeln!(
                w,
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                r.tss_idx, r.n1, r.nf, r.n2, r.nfr_score, r.log2_mean_cov
            )?;
        }
    }
```

Wrap the PT TSV writer in:

```rust
    if let Some((_, _, pt_rows, _, _)) = &tss_metric_results {
        let p = metric_path(outdir, flat, "pt", &format!("{}.pt.tsv", sample));
        let mut w =
            BufWriter::new(File::create(&p).with_context(|| format!("create {}", p.display()))?);
        use std::io::Write as _;
        writeln!(w, "tss_idx\tpromoter\tbody\tpt_score\tlog2meancov")?;
        for r in pt_rows {
            writeln!(
                w,
                "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                r.tss_idx, r.promoter, r.body, r.pt_score, r.log2_mean_cov
            )?;
        }
    }
```

Wrap the TSSE SVG writer in:

```rust
    if let Some((tsse_result, _, _, _, _)) = &tss_metric_results {
        let p = metric_path(outdir, flat, "tsse", &format!("{}.tsse.svg", sample));
        plots::tsse_svg(&tsse_result.values, &p, &sample)
            .with_context(|| format!("TSSE SVG: {}", p.display()))?;
    }
```

When building summary paths, replace the unconditional TSS path variables with optional path variables:

```rust
    let tsse_tsv_path = tss_metric_results.as_ref().map(|_| {
        if flat {
            format!("{}.tsse.tsv", sample)
        } else {
            format!("tsse/{}.tsse.tsv", sample)
        }
    });
    let nfr_tsv_path = tss_metric_results.as_ref().map(|_| {
        if flat {
            format!("{}.nfr.tsv", sample)
        } else {
            format!("nfr/{}.nfr.tsv", sample)
        }
    });
    let pt_tsv_path = tss_metric_results.as_ref().map(|_| {
        if flat {
            format!("{}.pt.tsv", sample)
        } else {
            format!("pt/{}.pt.tsv", sample)
        }
    });
```

In `summary::AtacSummary`, replace the TSS section construction with:

```rust
        tsse: match (&tss_metric_results, &tsse_tsv_path) {
            (Some((tsse_result, _, _, _, _)), Some(tsv_path)) => Some(summary::TsseSection {
                score: tsse_result.tsse_score,
                n_windows: tsse_result.values.len() as u32,
                values: tsse_result.values.clone(),
                tsv_path: tsv_path.clone(),
            }),
            _ => None,
        },
        nfr: match (&tss_metric_results, &nfr_tsv_path) {
            (Some((_, nfr_rows, _, nfr_median, _)), Some(tsv_path)) => Some(summary::ScoreSection {
                n_tss: nfr_rows.len() as u32,
                median_score: *nfr_median,
                tsv_path: tsv_path.clone(),
            }),
            _ => None,
        },
        pt: match (&tss_metric_results, &pt_tsv_path) {
            (Some((_, _, pt_rows, _, pt_median)), Some(tsv_path)) => Some(summary::ScoreSection {
                n_tss: pt_rows.len() as u32,
                median_score: *pt_median,
                tsv_path: tsv_path.clone(),
            }),
            _ => None,
        },
```

This avoids `expect` in production code while still keeping the path and result state tied together.

- [ ] **Step 6: Run ATAC unit and smoke tests**

Run:

```bash
cargo test atac::tests
cargo test atac::summary::tests
cargo test --test integration_atac rustqc_atac_runs_on_gl1_fixture
cargo test --test integration_atac gl1_metrics_smoke
```

Expected: all listed tests pass, and default behavior still writes TSSE/NFR/PT outputs.

- [ ] **Step 7: Commit**

Run:

```bash
git add src/atac/mod.rs
git commit -m "feat(atac): gate TSS metrics on Tn5 shift state"
```

---

### Task 5: Integration Tests and Documentation

**Files:**
- Modify: `tests/integration_atac.rs`
- Modify: `docs/src/content/docs/atac/cli.mdx`

- [ ] **Step 1: Add integration tests for basic-only and conflict behavior**

In `tests/integration_atac.rs`, after `rustqc_atac_runs_on_gl1_fixture`, add:

```rust
#[test]
fn tn5_shift_no_unshifted_input_writes_basic_qc_only() {
    let outdir = tempfile::tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .args([
            "atac",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL1.bam"),
            "--gtf",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL_tss.gtf"),
            "--outdir",
            outdir.path().to_str().unwrap(),
            "--sample-name",
            "GL1",
            "--tn5-shift",
            "no",
        ])
        .status()
        .unwrap();
    assert!(
        status.success(),
        "rustqc atac --tn5-shift no exited non-zero: {:?}",
        status
    );

    for sub in [
        "bamqc/GL1.bamqc.tsv",
        "bamqc/GL1.mapq.tsv",
        "fragsize/GL1.fragsize.tsv",
        "fragsize/GL1.fragsize.svg",
        "lib_complexity/GL1.libcomplexity.tsv",
        "lib_complexity/GL1.libcomplexity.svg",
        "GL1.atac.summary.json",
    ] {
        assert!(
            outdir.path().join(sub).exists(),
            "missing basic output file: {}",
            sub
        );
    }

    for sub in [
        "tsse/GL1.tsse.tsv",
        "tsse/GL1.tsse.svg",
        "nfr/GL1.nfr.tsv",
        "pt/GL1.pt.tsv",
    ] {
        assert!(
            !outdir.path().join(sub).exists(),
            "TSS-dependent output should be skipped: {}",
            sub
        );
    }

    let summary = read_summary(&outdir, "GL1");
    assert_eq!(summary["tn5_shift"]["requested"].as_bool(), Some(false));
    assert_eq!(summary["tn5_shift"]["input_is_shifted"].as_bool(), Some(false));
    assert_eq!(summary["tn5_shift"]["applied"].as_bool(), Some(false));
    assert_eq!(
        summary["tn5_shift"]["tss_dependent_metrics_enabled"].as_bool(),
        Some(false)
    );
    assert!(summary["tsse"].is_null());
    assert!(summary["nfr"].is_null());
    assert!(summary["pt"].is_null());
}

#[test]
fn tn5_shift_yes_with_input_is_shifted_errors() {
    let outdir = tempfile::tempdir().unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .args([
            "atac",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL1.bam"),
            "--gtf",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL_tss.gtf"),
            "--outdir",
            outdir.path().to_str().unwrap(),
            "--sample-name",
            "GL1",
            "--tn5-shift",
            "yes",
            "--input-is-shifted",
        ])
        .output()
        .unwrap();
    assert!(
        !output.status.success(),
        "contradictory Tn5 flags should fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--tn5-shift yes") && stderr.contains("--input-is-shifted"),
        "stderr did not explain contradictory flags: {}",
        stderr
    );
}
```

- [ ] **Step 2: Run integration tests**

Run:

```bash
cargo test --test integration_atac tn5_shift
```

Expected: both tests pass.

- [ ] **Step 3: Update CLI docs**

In `docs/src/content/docs/atac/cli.mdx`, in the ATAC-specific options table, add rows after `--mito-chrom`:

```md
| `--tn5-shift <yes\|no>` | `RUSTQC_TN5_SHIFT` | `yes` | Whether RustQC should apply +4/−5 Tn5 shift for insertion-site metrics |
| `--input-is-shifted[=true\|false]` | `RUSTQC_INPUT_IS_SHIFTED` | false | Declare that input BAM coordinates are already Tn5 shifted |
```

After the examples block, add:

```md
## Tn5 shift semantics

TSS-dependent metrics (`TSSEscore`, `NFRscore`, and `PTscore`) use Tn5 insertion-site coordinates.
By default, `rustqc atac` assumes ordinary unshifted ATAC BAM input and applies +4/−5 shift in
memory for those metrics.

| `--tn5-shift` | `--input-is-shifted` | Behavior |
|---------------|----------------------|----------|
| `yes` | absent or `false` | Apply in-memory Tn5 shift and write all QC outputs |
| `yes` | `true` | Exit with an error to avoid double-shifting |
| `no` | `true` | Treat input coordinates as already shifted and write all QC outputs |
| `no` | absent or `false` | Write only basic QC outputs; `tsse`, `nfr`, and `pt` are `null` in the JSON summary |

Basic QC outputs include bamQC, MAPQ histogram, fragment size, library complexity, and the JSON
summary. `--emit-shifted-bam` remains an output-emission flag and does not control metric
semantics.
```

- [ ] **Step 4: Run full relevant verification**

Run:

```bash
cargo fmt --check
cargo test cli::tests::test_atac
cargo test config::tests::test_atac_shift_config
cargo test atac::tests
cargo test atac::summary::tests
cargo test --test integration_atac rustqc_atac_runs_on_gl1_fixture
cargo test --test integration_atac gl1_metrics_smoke
cargo test --test integration_atac tn5_shift
```

Expected: all commands pass.

- [ ] **Step 5: Commit**

Run:

```bash
git add tests/integration_atac.rs docs/src/content/docs/atac/cli.mdx
git commit -m "test(atac): cover Tn5 shift CLI behavior"
```

---

## Final Verification

- [ ] Run formatting:

```bash
cargo fmt --check
```

Expected: success.

- [ ] Run ATAC-focused tests:

```bash
cargo test cli::tests::test_atac
cargo test config::tests::test_atac_shift_config
cargo test atac::tests
cargo test atac::summary::tests
cargo test --test integration_atac
```

Expected: success.

- [ ] Run broader release-gate tests if time permits:

```bash
cargo test
cargo clippy -- -D warnings
```

Expected: success.

- [ ] Inspect final diff:

```bash
git status --short
git log --oneline -5
```

Expected: worktree is clean after the planned commits; recent commits show the parsing, resolution, driver gating, and integration-test changes.
