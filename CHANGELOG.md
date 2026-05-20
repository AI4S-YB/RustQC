# RustQC Changelog

## [Unreleased]

## [Version 0.4.1](https://github.com/AI4S-YB/RustQC/releases/tag/v0.4.1) - 2026-05-20

### Added

- Explicit ATAC Tn5 shift controls: `--tn5-shift <yes|no>` and
  `--input-is-shifted`, with matching YAML/env configuration and validation
  for contradictory shifted-input settings.

### Changed

- `rustqc atac --tn5-shift no` on unshifted input now emits shift-independent
  QC only, keeps skipped TSS-dependent summary sections as JSON `null`, and
  reports the resolved `tn5_shift` state in the summary.

## [Version 0.4.0](https://github.com/AI4S-YB/RustQC/releases/tag/v0.4.0) - 2026-05-05

### Added

- `rustqc atac` subcommand for single-pass ATAC-seq QC, replicating
  ATACseqQC 1.36.0's `bamQC`, `fragSizeDist`, `TSSEscore`, `NFRscore`,
  `PTscore`, and `estimateLibComplexity`. Numerical fidelity targets the
  upstream R package (byte-identical / 1e-6 / 1e-3 per metric — see the
  documentation site).
- Tn5 +4/−5 shift and fixed-interval (NFR / mono / di / tri) length split
  helpers — opt-in via `--emit-shifted-bam` / `--emit-split-bams` (file
  writing reserved for a follow-up release).
- Shared BAM/preseq infrastructure lifted out of `src/rna/` to crate root
  (`src/bam_flags.rs`, `src/bam_io.rs`, `src/cpp_rng.rs`, `src/preseq.rs`).
  No behavior change to existing `rustqc rna` outputs.

### Known gaps

- Factor-footprinting metrics (`factorFootprints`, `vPlot`, etc.) —
  requires BSgenome / PWM motif DBs.
- Random-forest split branch (only the fixed-interval split lands here).
- `EmitWriters::open` actual file writing for `--emit-shifted-bam` /
  `--emit-split-bams`.

## [Version 0.3.0](https://github.com/AI4S-YB/RustQC/releases/tag/v0.3.0) - 2026-04-22

First release of the AI4S-YB fork. Focus: Windows support and a pure-Rust
alignment-file backend. This is a **hard fork** of `seqeralabs/RustQC` v0.2.1
and is not published to crates.io.

### Breaking changes

- **CRAM input is no longer supported.** The `rust-htslib` (HTSlib C library)
  backend was replaced with [`noodles`](https://github.com/zaeleus/noodles).
  `noodles-cram` 0.88 uses Rust 1.88+ syntax (let-chains), which conflicts
  with the project's MSRV of 1.87, so CRAM was dropped rather than bumping
  the toolchain. BAM and SAM remain fully supported. CRAM can return by
  raising MSRV to ≥ 1.88 and re-adding `noodles-cram`.

### New features

- **Native Windows builds.** `x86_64-pc-windows-msvc` is now a first-class
  target. Release artifacts include `rustqc-windows-x86_64.zip` alongside
  the existing `.tar.gz` archives for Linux and macOS.

### Internal changes

- New `src/rna/bam_io.rs` facade module centralizes the semantic differences
  between `rust-htslib` and `noodles` (MAPQ=255 sentinel handling, 1-based
  ↔ 0-based position conversion, aux-tag integer extraction, CIGAR op
  collection, QNAME `*` fallback, 4-bit encoded-base access for sequence
  hashing). `open()` auto-detects SAM vs BAM and transcodes SAM through
  noodles at load time so downstream code only sees `bam::Record`.
- Build pipeline simplified: the htslib-era C dependency chain
  (`hts-sys`, `openssl-sys`, `curl-sys`, `libz-sys`, `bzip2-sys`,
  `lzma-sys`, `libclang`) is gone. Linux CI deps shrank from eight
  system packages to two (`libfontconfig1-dev`, `pkg-config`).
- `build.rs` no longer shells out to `date`; `qualimap/report.rs`
  uses `chrono::Local` for local-time formatting. These were the
  last POSIX-only hold-outs blocking Windows.
- `plotters`'s `fontconfig-dlopen` feature is now gated behind
  `cfg(not(windows))`.

### Known regressions

- Multithreaded BGZF decode (the old `bam.set_threads(n)` path) is
  disabled. Large BAM throughput is lower than upstream v0.2.1 until
  `noodles_bgzf::io::MultithreadedReader` is wired in; correctness is
  unaffected. See `TODO(noodles-threading)` markers in
  `src/rna/dupradar/counting.rs`.

### Tests

- All 234 unit and integration tests pass on Linux, macOS and Windows.
- `preseq lc_extrap` output is byte-identical to the pre-migration
  reference (`tests/data/test.preseq_lc_extrap.txt`).

## [Version 0.2.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.1) - 2026-04-09

### Bug fixes

- Fix SIMD builds by scoping rustflags to target triple with explicit `--target` (#90, #91)

### Other changes

- Trigger releases via workflow dispatch instead of tag push (#92)

## [Version 0.2.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.2.0) - 2026-04-09

### Features

- Ship SIMD-optimized binaries with CPU detection and upgrade hints (#81)
- Write `CITATIONS.md` with upstream tool versions (#87)
- Add XDG config discovery and deep-merge support (#88)

### Bug fixes

- Replace header-based duplicate check with flag-based detection (#84)
- Use `.log` extension for junction_annotation output (#80)

### Other changes

- Bump docker/login-action from 4.0.0 to 4.1.0 (#78)
- Bump vite from 7.3.1 to 7.3.2 in docs (#77)
- Bump defu from 6.1.4 to 6.1.6 in docs (#74)

## [Version 0.1.1](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.1) - 2026-04-02

### Bug fixes

- Fix featureCounts summary to use gene-level stats; add biotype summary (#66)
- Fix inner_distance histogram to include overflow bucket in bulk cutoff loop (#67)

### Other changes

- Add crates.io publishing to release workflow (#62)
- Documentation fixes (#70)

## [Version 0.1.0](https://github.com/seqeralabs/RustQC/releases/tag/v0.1.0) - 2026-04-01

Initial release of RustQC -- fast quality control tools for sequencing data, written in Rust.

A single `rustqc rna` command runs 15 QC analyses in one pass over the BAM file, producing output that is format- and numerically identical to the upstream tools and fully compatible with [MultiQC](https://multiqc.info/).

### Tools

- **dupRadar** -- PCR duplicate rate vs. expression analysis with density scatter plots, boxplots, and expression histograms. 14-column duplication matrix with logistic-regression model fitting matching the [R dupRadar](https://github.com/ssayols/dupRadar) package.
- **featureCounts** -- Gene-level read counting with assignment summary, compatible with [Subread featureCounts](http://subread.sourceforge.net/). Includes per-biotype read counting and MultiQC integration.
- **RSeQC** (8 tools) -- [RSeQC](https://rseqc.sourceforge.net/)-compatible implementations of bam_stat, infer_experiment, read_duplication, read_distribution, junction_annotation, junction_saturation, inner_distance, and TIN (Transcript Integrity Number). Includes native plot generation (PNG + SVG) with no R dependency.
- **preseq** -- Library complexity extrapolation (`lc_extrap`) matching [preseq](http://smithlabresearch.org/software/preseq/) v3.2.0, including C++ RNG FFI for reproducible bootstrap sampling.
- **Qualimap rnaseq** -- Gene body coverage profiling, read genomic origin, junction analysis, and transcript coverage bias matching [Qualimap](http://qualimap.conesalab.org/).
- **samtools** -- flagstat, idxstats, and full stats output (SN section + all histogram sections) matching [samtools](http://www.htslib.org/).

### Features

- Single static binary with no runtime dependencies
- SAM, BAM, and CRAM input support (auto-detected)
- Multi-threaded alignment processing with htslib thread pool
- GTF annotation support (gzip-compressed files accepted)
- YAML configuration for output control, chromosome name mapping, and tool parameters
- Multiple BAM file support via positional arguments
- `--sample-name` flag to override BAM-derived sample name in output filenames
- Per-tool seed flags (`--tin-seed`, `--junction-saturation-seed`, `--preseq-seed`) for reproducible results
- Docker container at `ghcr.io/seqeralabs/rustqc`
- Cross-platform builds (Linux x86_64/aarch64, macOS x86_64/aarch64)
