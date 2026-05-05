# RustQC ATAC-seq QC Subcommand — Design

- **Date**: 2026-05-04
- **Author**: brainstormed with user
- **Status**: design approved, awaiting plan
- **Reference**: ATACseqQC 1.36.0 (Bioconductor 3.23, R package, GPL ≥ 2)

## Goal

Add a new `rustqc atac` subcommand parallel to the existing `rustqc rna`, providing single-pass ATAC-seq QC and Tn5 preprocessing as a fast, dependency-free Rust tool. Numerical outputs match the upstream ATACseqQC R package; pipeline-level outputs (file naming, MultiQC compatibility) align with nf-core/atacseq.

## Scope

### v1 — In scope

- **Core QC metrics** (numerical fidelity to ATACseqQC):
  - `bamQC` — total reads, duplicate / mito / proper-pair / unmapped / QC-fail rates, NRF, PBC1, PBC2, MAPQ histogram, idxstats
  - `fragSizeDist` — fragment-length histogram (1–1010 bp) plus dual-window plot (linear + log10 inset)
  - `TSSEscore` — TSS enrichment score with sliding-window normalization
  - `NFRscore` — nucleosome-free region score per TSS
  - `PTscore` — promoter / transcript-body ratio per TSS
  - `estimateLibComplexity` / `readsDupFreq` / `saturationPlot` — library complexity extrapolation reusing the existing Rust preseq port
- **Tn5 preprocessing** (opt-in BAM emission):
  - `shiftReads` — +4 / −5 strand-aware shift; in-memory by default
  - `splitGAlignmentsByCut` — fixed-interval split into NFR / mono / di / tri BAMs

### v1 — Out of scope (TODO)

- Factor-footprinting suite (`factorFootprints`, `footprintsScanner`, `vPlot`, `enrichedFragments`, `distanceDyad`, `plotFootprints`) — depends on BSgenome / PWM motif DBs, conflicts with the static-binary distribution model.
- Random-forest model for `splitGAlignmentsByCut` — only the fixed-interval branch is implemented; the RF model is a future refinement.
- ChIP-seq subcommand — out of scope, but the §2 refactor extracts shared infra so it could be added later without disturbing `rna/` or `atac/`.

## Inputs and CLI

```
rustqc atac <INPUT.bam>... --gtf <GTF> [options]

Input / Output:
  -g, --gtf <GTF>              GTF annotation, plain or .gz; TSS coords source (required)
  -r, --reference <FASTA>      Reference FASTA (required for CRAM)
  -o, --outdir <DIR>           Output directory [default: .]
      --sample-name <NAME>     Override sample name (default: derived from BAM filename)
      --flat-output            Write outputs to a flat directory (no subdirs)
  -c, --config <YAML>          YAML configuration file
  -j, --json-summary <PATH>    JSON summary path (use "-" for stdout)

ATAC-specific:
      --mito-chrom <NAME>      Mitochondrial chromosome name (default: auto-detect ^chrM$|^MT$|^Mito$)
      --emit-shifted-bam       Emit +4/−5 Tn5-shifted BAM
      --emit-split-bams        Emit NFR/mono/di/tri BAMs (fixed intervals)
      --tsse-flank <N>         TSSEscore flank window (default 1000)

General (shared with rna):
  -t, --threads, -Q --mapq, -q --quiet, -v --verbose
```

**Paired-end requirement**: `rustqc atac` rejects single-end BAMs at startup. The first N records (e.g. 10 000) are inspected for the `READ_PAIRED` flag; if none are paired, the tool exits with a clear error message rather than producing meaningless fragSize / NFR / TSSE numbers.

## Architecture

### Module layout

```
src/
  bam_flags.rs        # moved from rna/ (§2 refactor)
  bam_io.rs           # moved from rna/; PE pair reconstruction, noodles wrappers
  cpp_rng.rs          # moved from rna/
  preseq.rs           # moved from rna/
  cli.rs              # add Atac(AtacArgs) to Commands
  config.rs           # add atac: AtacConfig
  main.rs             # route atac → atac::run
  gtf.rs              # shared, with a small TSS-extraction helper
  rna/
    mod.rs
    dupradar/
    featurecounts/
    qualimap/
    rseqc/
  atac/
    mod.rs
    bam_qc.rs         # bamQC: rates, NRF, PBC1/2, MAPQ histogram
    frag_size.rs      # fragSizeDist
    tsse.rs           # TSSEscore
    nfr_score.rs      # NFRscore
    pt_score.rs       # PTscore
    lib_complexity.rs # readsDupFreq + saturationPlot via crate::preseq
    shift.rs          # Tn5 +4/−5 shift, in-memory
    split.rs          # fixed-interval NFR/mono/di/tri split
    bam_writer.rs     # opt-in BAM + .bai output via noodles
    plots.rs          # plotters backend for fragSize / TSSE / saturation SVGs
    summary.rs        # aggregate JSON
```

### §2 — Shared-infra refactor (prerequisite commit)

Before touching ATAC code, lift four files out of `src/rna/` to the crate root, **with no behavior change**:

| Current | New |
|---|---|
| `src/rna/bam_flags.rs` | `src/bam_flags.rs` |
| `src/rna/bam_io.rs` | `src/bam_io.rs` |
| `src/rna/cpp_rng.rs` | `src/cpp_rng.rs` |
| `src/rna/preseq.rs` | `src/preseq.rs` |

Done as a single commit: `refactor: extract shared BAM/preseq infra out of rna module`. Visibility of `preseq::ds_rsac_bootstrap` (and any other entrypoints `lib_complexity.rs` needs) is widened to `pub(crate)` so `atac/` can call it. Regression bar: `cargo test --workspace` green and existing RNA integration outputs byte-identical to pre-refactor.

If during the move a part of `bam_io.rs` proves RNA-specific (e.g. splice-junction-aware CIGAR helpers, `XS:A:` strand tag), the RNA-specific subset stays in `src/rna/bam_io_rna.rs` and only the truly generic part moves to `src/bam_io.rs`. Goal: keep RNA byte-identical first, share aggressively second.

### Single-pass BAM scan

All ATAC metrics are computed in one streaming pass over the BAM:

1. Per record: update fragSize histogram (every primary record contributes `|TLEN|`, so each fragment is counted twice — matches the R behavior, where the downstream density normalization cancels the 2× factor), per-chromosome MAPQ histogram, mito counter, duplicate / proper-pair / QC-fail / unmapped flag counters, idxstats counter.
2. Track 5'-end positions for PBC1/PBC2/NRF and for `readsDupFreq` (duplicate-frequency histogram fed to preseq).
3. For records overlapping any TSS ± `max(promoter_window, tsse_flank + nfr_pad, pt_pad)` window, accumulate per-TSS 5'-end coverage into a sparse buffer keyed by transcript id.
4. If `--emit-shifted-bam` / `--emit-split-bams`: stream Tn5-shifted records to the corresponding writer(s) inline.

After the streaming phase: aggregate per-chromosome stats into bamQC summary; finalize TSSE / NFR / PT from the per-TSS coverage windows; run preseq bootstrap; emit plots and JSON.

## Algorithms — numerical fidelity

Every formula below is read directly from ATACseqQC 1.36.0 source (`R/*.R`).

### bamQC (`R/bamQC.R`)

For each `@SQ` (chromosome) in the header, scan all primary records:

- `lenQ` = number of records.
- `isMitochondria[i]` = `seqname == mito_chrom`.
- Flags from BAM bitfield: `isDuplicate`, `isProperPair`, `isUnmappedQuery`, `hasUnmappedMate`, `isNotPassingQualityControls`.
- MAPQ → per-chromosome histogram, summed across chromosomes at the end.
- Position fingerprint per chromosome:
  - PE: `(pos1, isize1, pos2, isize2)` tuple per read pair (singletons get NA second half).
  - SE: `(pos, qwidth)` per read.
- `M_DISTINCT` = number of distinct fingerprints.
- `M1` = fingerprints occurring exactly once (`!duplicated && !duplicatedFromLast`).
- `M2` = number of fingerprints occurring exactly twice.

Final aggregation across chromosomes:

```
totalQNAMEs           = |unique(qname)|
duplicateRate         = Σ isDup       / Σ lenQ
mitochondriaRate      = Σ isMito      / Σ lenQ
properPairRate        = Σ isProperPair/ Σ lenQ
unmappedRate          = Σ isUnmappedQuery / Σ lenQ
hasUnmappedMateRate   = Σ hasUnmappedMate / Σ lenQ
notPassingQCRate      = Σ isNotPassingQualityControls / Σ lenQ
NRF                   = ΣM1 / totalQNAMEs
PBC1                  = ΣM1 / ΣM_DISTINCT
PBC2                  = ΣM1 / max(1, ΣM2)
```

### fragSizeDist (`R/fragSizeDist.R`)

- Records: secondary excluded, unmapped excluded, QC-fail excluded.
- Fragment length: `|TLEN|` (R uses `abs(isize)`).
- Histogram: integer keys 1..1010 (R uses `match(1:1010, names(table))`), zero-fill missing.
- Plot:
  - Linear: x ∈ [0,1010], y = density × 1000, line plot.
  - log10 inset: x ∈ [0,1010], y = log10(density), inset at `fig=c(.4,.95,.4,.95)`, custom log-axis tick formatter.
- Output TSV columns: `length\tcount\tnorm_density`.

### TSSEscore (`R/TSSEscore.R`)

- `obj` = 5'-end coverage of all reads (per-base counts on each chromosome).
- `txs` = transcripts; deduplicated; restricted to chromosomes present in `obj`.
- `sel.center` = `promoters(txs, upstream=1000, downstream=1000)` (default).
- `sliding_windows(width=100, step=100)` → 20 windows per TSS.
- `vms.center[t,w]` = `viewMeans(coverage)` for window w of TSS t.
- `sel.left.flank` = leftmost `endSize=100bp` of `sel.center`; `sel.right.flank` = rightmost 100bp.
- `vms.left[t]`, `vms.right[t]` = mean coverage in those flanks.
- For each TSS:
  - Treat NA on one side as the other side; remaining NAs filled with `pseudocount=0`.
  - `blk = (vl + vr) / 2`; if `pseudocount ≤ 0`, drop TSS where `blk == 0`.
  - For surviving windows: `v_norm = v * endSize / blk / width` ≡ `v * 100 / blk / 100` ≡ `v / blk`.
- Mean across surviving TSS per window → 20-element vector `vms.m`.
- `loess.smooth(x=1..20, y=vms.m, family="gaussian", evaluation=20)`.
- `TSSE = max(y_loess[!is.infinite])`.
- Output: `values` (20 floats post-loess), `tsse_score` (scalar).

**Loess strategy**: port a minimal `loess` (locally-weighted polynomial regression, default tricube weights, `span=2/3`, degree=2) tuned to match `stats::loess.smooth`. Acceptance bar: TSSE numerical diff ≤ 1e-3 vs R on the GL1–GL4 fixtures. The pre-loess `vms.m` vector is byte-identical to R.

### NFRscore (`R/NFRscore.R`)

For each TSS, build three sub-windows (strand-aware), following ATACseqQC's `promoters()` / `shift()` semantics exactly:

```
nucleosomeSize    N = 150
nucleosomeFreeSize F = 100

sel = promoters(tss, upstream = N + floor(F/2), downstream = N + ceiling(F/2))
n1  = promoters(sel, upstream = 0, downstream = N)        # 150 bp upstream nucleosome window
n2  = shift(n1, +(N+F))     for + strand                  # 150 bp downstream nucleosome window
n2  = shift(n1, -(N+F))     for - strand
nf  = shift(n1, +N), then width <- F      for + strand    # 100 bp central NFR window
nf  = shift(n1, -N), then start <- end-F+1 for - strand
```

For default `(N, F) = (150, 100)` and a `+` strand TSS this yields `n1 = [TSS-200, TSS-51]`, `nf = [TSS-50, TSS+49]`, `n2 = [TSS+50, TSS+199]` (1-based, inclusive). Unit tests pin these exact boundaries so the implementation cannot drift off-by-one.

Mean coverage per window from the 5'-end coverage track. Then:

```
smallNumber  = max(1e-6, min(n1), min(n2), min(nf))
log2meanCov  = log2((3*(n1+n2) + 2*nf) / 8 + smallNumber)
NFR_score    = log2(nf + smallNumber) + 1 - log2(n1 + n2 + smallNumber)
```

Output one row per transcript: `gene_id\tn1\tnf\tn2\tnfr_score\tlog2meancov`. Acceptance bar: per-row diff ≤ 1e-6 vs R.

### PTscore (`R/PTscore.R`)

```
upstream   U = 2000
downstream D =  500

promoter   = [TSS-U, TSS+D]              (strand-aware)
body       = [TSS+D, TSS+D+U+D]          (shift forward by U+D, − strand reflected)

smallNumber  = max(1e-6, min(promoter), min(body))
log2meanCov  = log2(promoter + smallNumber) + log2(body + smallNumber)
PT_score     = log2(promoter + smallNumber) - log2(body + smallNumber)
```

Output: `gene_id\tpromoter\tbody\tpt_score\tlog2meancov`. Acceptance bar: per-row diff ≤ 1e-6 vs R.

### estimateLibComplexity (`R/estimateLibComplexity.R`)

- `readsDupFreq`: per BAM record, build fingerprint key (PE: `(chr, leftpos, isize)`, SE: `(chr, pos, qwidth)`); count occurrences `j`; emit histogram `(j, #fingerprints with multiplicity j)`.
- Feed histogram to `preseqR::ds.rSAC.bootstrap(hist, r=1, times=100)` → fitted SAC function `f`.
- Evaluate at `relative.size ∈ {0.1, 0.2, …, 1.0, 5, 10, 15, 20}`.
- Output rows: `relative_size\tdistinct_fragments\tputative_reads`.

We reuse `crate::preseq::ds_rsac_bootstrap` from the moved-up preseq module; the §2 refactor exposes whatever entrypoints `lib_complexity.rs` needs (`pub(crate)` is sufficient).

### shiftReads (`R/shiftReads.R`)

For each alignment:

- `+ strand`: shift POS by +4; trim 4 bases from the 5' end of SEQ/QUAL (after first resolving any soft-clip via CIGAR).
- `− strand`: shift end by −5 (i.e. reduce alignment length from the 3'-genomic / 5'-read end); trim 5 bases from the 5' end of SEQ/QUAL.
- TLEN: `sign(TLEN) * (|TLEN| - 9)` when TLEN is set.
- CIGAR: equivalent to `cigarQNarrow(start = strand=='-' ? 1 : 5, end = strand=='-' ? -6 : -1)`; soft-clipped 5' bases are first folded into the read via `sequenceLayer(from='query', to='query-after-soft-clipping')`.
- Edge cases:
  - Insertions at the 5' end such that `cigarWidthAlongQuerySpace != width(seq)`: clip from 3' end.
  - `qwidth` after shift must be > 0; reads that would have non-positive width are dropped (with a counter logged).
- Output BAM (when `--emit-shifted-bam`): same header as input, records reordered by Tn5-shifted POS within each chromosome (re-sort phase), `.bai` written.

### splitGAlignmentsByCut — fixed-interval branch

Bucket each fragment by `|TLEN|` (after Tn5 shift):

```
NFR  : [0, 100)
mono : [180, 247]
di   : [315, 473]
tri  : [558, 615]
(other: dropped — not written to any output BAM)
```

Random-forest model branch from R is **TODO**; CLI help and the JSON summary include `"split_method": "fixed_intervals_v1"` so downstream tools can tell the two apart later.

## Outputs and pipeline integration

### Directory layout (default)

```
outdir/
  bamqc/             <sample>.bamqc.json, <sample>.bamqc.tsv, <sample>.mapq.tsv
  fragsize/          <sample>.fragsize.tsv, <sample>.fragsize.svg
  tsse/              <sample>.tsse.tsv, <sample>.tsse.svg
  nfr/               <sample>.nfr.tsv
  pt/                <sample>.pt.tsv
  lib_complexity/    <sample>.libcomplexity.tsv, <sample>.libcomplexity.svg
  shifted/           <sample>.shifted.bam(.bai)        # only with --emit-shifted-bam
  split/             <sample>.NFR.bam(.bai), .mono.bam(.bai), .di.bam(.bai), .tri.bam(.bai)
                                                       # only with --emit-split-bams
  <sample>.atac.summary.json
```

`--flat-output` collapses everything into `outdir/`. `-j <path>` redirects only the summary JSON.

### nf-core/atacseq + MultiQC alignment

- Per-file naming (`<sample>.fragsize.tsv`, `<sample>.tsse.tsv`, etc.) follows the nf-core/atacseq convention of `<sample>.<tool>` prefixes so the existing MultiQC `atacseqqc` / custom_content modules pick them up unchanged.
- TSV column order matches the R package's `write.csv` / `print.data.frame` conventions where applicable, so consumers parsing nf-core outputs see the same shape.

### JSON summary schema

```json
{
  "sample": "<derived or --sample-name>",
  "tool_versions": {
    "rustqc": "<crate version>",
    "atacseqqc_replicates": "1.36.0"
  },
  "split_method": "fixed_intervals_v1",
  "bamqc": {
    "total_qnames": <u64>,
    "duplicate_rate": <f64>,
    "mitochondria_rate": <f64>,
    "proper_pair_rate": <f64>,
    "unmapped_rate": <f64>,
    "has_unmapped_mate_rate": <f64>,
    "not_passing_qc_rate": <f64>,
    "nrf": <f64>, "pbc1": <f64>, "pbc2": <f64>,
    "mapq_histogram": { "<mapq>": <count>, ... }
  },
  "fragsize": { "total_pairs": <u64>, "tsv_path": "fragsize/<sample>.fragsize.tsv" },
  "tsse":     { "score": <f64>, "values": [<20 floats>], "tsv_path": "..." },
  "nfr":      { "median_score": <f64>, "tsv_path": "..." },
  "pt":       { "median_score": <f64>, "tsv_path": "..." },
  "lib_complexity": {
    "extrapolated_total": <f64>,
    "tsv_path": "lib_complexity/<sample>.libcomplexity.tsv"
  }
}
```

## Testing strategy

### Unit tests (in-tree, `cargo test`)

- `bam_qc.rs`: PBC1/PBC2/NRF formulas, mito-detection regex, MAPQ histogram aggregation, edge case `M2 == 0`.
- `frag_size.rs`: histogram bounds (1..1010), zero-fill, density normalization.
- `tsse.rs`: per-TSS normalization with `pseudocount = 0` (TSS dropped), `pseudocount > 0` (TSS kept).
- `nfr_score.rs` / `pt_score.rs`: strand-aware coordinate construction, `smallNumber` floor activation, log2 sign on degenerate inputs.
- `shift.rs`: + / − strand shift, soft-clip absorption, TLEN halving, qwidth-zero drop.
- `split.rs`: bucket boundaries 100 / 180 / 247 / 315 / 473 / 558 / 615.

### Integration / numerical-fidelity tests

`tests/atac/` driven by ATACseqQC's own fixtures `inst/extdata/GL{1..4}.bam` (extracted from the tarball into `tests/data/atac/`). Reference outputs are produced by a one-off R script (`tests/atac/golden/run_r_reference.R`); the resulting TSV/JSON fixtures are committed to the repo. CI does **not** invoke R — Rust tests read the fixtures directly.

| Metric | Acceptance bar |
|---|---|
| `fragSizeDist` | byte-identical histogram |
| `bamQC` rates / NRF / PBC1 / PBC2 / MAPQ histogram | byte-identical |
| `NFRscore` (per TSS) | abs diff ≤ 1e-6 |
| `PTscore` (per TSS) | abs diff ≤ 1e-6 |
| `TSSEscore` (post-loess scalar) | abs diff ≤ 1e-3 |
| `TSSEscore` (pre-loess `vms.m`) | byte-identical |
| `shiftReads` BAM | `samtools view` POS / CIGAR / SEQ / QUAL byte-identical on GL1.bam |
| `splitGAlignmentsByCut` (fixed) | each output BAM's read-name set byte-identical |

### RNA regression after §2 refactor

After lifting `bam_flags`, `bam_io`, `cpp_rng`, `preseq` to crate root, the existing RNA integration tests under `tests/` must produce byte-identical output to pre-refactor. This is the gate before any ATAC code is written.

### Performance bar

On nf-core/atacseq's standard test sample, single-pass `rustqc atac` is at least ~10× faster than the equivalent ATACseqQC R workflow (which scans the BAM multiple times). Exact numbers go in the README after benchmarking; this design only commits to the ≥10× target as a sanity check, not a contractual SLA.

## Open items / future work

- **Random-forest split branch**: port ATACseqQC's RF model for `splitGAlignmentsByCut`; would need to ship serialized model weights with the binary or as a side-loaded asset.
- **Factor footprinting (scope C)**: `factorFootprints`, `footprintsScanner`, `vPlot`, `enrichedFragments`. Requires reference FASTA + PWM motif input. Probably best as a separate `rustqc atac-footprint` subcommand to keep the QC binary lean.
- **CRAM input**: covered today by the same `--reference` flag pattern as `rustqc rna`; ATAC will inherit whatever CRAM support `rna` adds.
- **ChIP-seq subcommand**: not in this project; the §2 refactor enables it without further disruption.
