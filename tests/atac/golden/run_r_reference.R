#!/usr/bin/env Rscript
# =============================================================================
# RustQC Phase 14 — R Reference Golden Generator
# =============================================================================
#
# PURPOSE
# -------
# This script generates JSON/TSV "golden" output files for numerical fidelity
# tests in Phase 14. The Rust integration tests in tests/integration_atac.rs
# compare rustqc atac metric outputs against these R-generated values.
#
# This script is OFFLINE-ONLY and is NOT run by CI. It must be run manually
# by a developer who has R and the required packages installed.
#
# After running this script, commit the resulting golden files so that the
# integration tests auto-activate and perform numerical comparisons:
#   tests/atac/golden/GL1.bamqc.golden.json
#   tests/atac/golden/GL1.fragsize.golden.tsv
#   tests/atac/golden/GL1.nfr.golden.tsv
#   tests/atac/golden/GL1.pt.golden.tsv
#   tests/atac/golden/GL1.tsse.golden.json
#   (same for GL2, GL3, GL4)
#
# CHROMOSOME ASSUMPTIONS
# ----------------------
# All GL1-4 BAMs align exclusively to chr1 (hg19/GRCh37). GL2/GL3/GL4 reads
# are concentrated in chr1:565,550-996,878; GL1 reads span
# chr1:565,608-249,202,181. The GTF used here (GL_tss.gtf) covers 14 synthetic
# transcripts across chr1 matching this distribution.
#
# If BAM chromosomes change (they are fixed ATACseqQC 1.36.0 fixtures and
# should not change), regenerate tests/data/atac/GL_tss.gtf accordingly.
#
# REQUIRED R PACKAGES
# -------------------
# - ATACseqQC (>= 1.36.0)   via BiocManager::install("ATACseqQC")
# - GenomicFeatures          via BiocManager::install("GenomicFeatures")
# - GenomicAlignments        via BiocManager::install("GenomicAlignments")
# - jsonlite                 via install.packages("jsonlite")
# - rtracklayer              via BiocManager::install("rtracklayer")
# - BSgenome.Hsapiens.UCSC.hg19  (for shiftReads; optional but recommended)
#     via BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#
# HOW TO RUN
# ----------
# From the repository root:
#   cd /path/to/RustQC
#   Rscript tests/atac/golden/run_r_reference.R
#
# Or from within R:
#   source("tests/atac/golden/run_r_reference.R")
#
# GOLDEN FILE FORMAT CONVENTIONS
# --------------------------------
# .bamqc.golden.json  — JSON object with same fields as rustqc JSON "bamqc" section
# .fragsize.golden.tsv — TSV: columns frag_size, count (one row per size 1..2000)
# .nfr.golden.tsv     — TSV: columns tss_idx, n1, nf, n2, nfr_score, log2meancov
# .pt.golden.tsv      — TSV: columns tss_idx, promoter, body, pt_score, log2meancov
# .tsse.golden.json   — JSON object: {"tsse_score": float, "values": [float x 20]}
#
# =============================================================================

suppressPackageStartupMessages({
  library(ATACseqQC)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(jsonlite)
  library(rtracklayer)
  if (requireNamespace("txdbmaker", quietly = TRUE)) library(txdbmaker)
})

# Bioconductor moved makeTxDbFromGFF from GenomicFeatures to txdbmaker.
# Use the txdbmaker version when available, otherwise fall back to the
# GenomicFeatures version (older Bioc releases).
make_txdb_from_gff <- function(path) {
  # txdbmaker uses lowercase "gtf"; GenomicFeatures historically accepted "GTF"
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(path, format = "gtf")
  } else {
    GenomicFeatures::makeTxDbFromGFF(path, format = "GTF")
  }
}

# ---- Paths ------------------------------------------------------------------
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else if (!is.null(sys.frame(1)$ofile)) {
  sys.frame(1)$ofile
} else {
  "."
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root  <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
data_dir    <- file.path(repo_root, "tests", "data", "atac")
golden_dir  <- file.path(repo_root, "tests", "atac", "golden")
gtf_path    <- file.path(data_dir, "GL_tss.gtf")
bam_dir     <- data_dir

# Validate paths
stopifnot(
  "data/atac directory not found"    = dir.exists(data_dir),
  "GL_tss.gtf not found"             = file.exists(gtf_path),
  "golden directory not found"       = dir.exists(golden_dir)
)

# ---- Load TSS from synthetic GTF -------------------------------------------
message("Loading TSS from GL_tss.gtf ...")
txdb_local <- make_txdb_from_gff(gtf_path)
HG19_CHR1_LEN <- 249250621L
txs <- transcripts(txdb_local)
seqlevels(txs, pruning.mode = "coarse") <- "chr1"
suppressWarnings(seqlengths(txs) <- c(chr1 = HG19_CHR1_LEN))
# GRanges spanning chr1 — used for the `which=` arg of readBamFile.
chr1_full <- GRanges(seqnames = "chr1",
                     ranges   = IRanges(1, HG19_CHR1_LEN))
message("  Loaded ", length(txs), " transcripts on chr1")

# ---- Helper: write JSON golden for bamQC ------------------------------------
write_bamqc_golden <- function(sample, bam_path) {
  message("bamQC for ", sample, " ...")
  qc <- bamQC(bam_path, outPath = NULL)
  # ATACseqQC bamQC returns a list; key fields:
  #   duplicateRate, mitochondriaRate, properPairRate, NRF, PBC1, PBC2,
  #   totalQnameSorted
  # ATACseqQC 1.34+ uses these field names (older versions had `totalQnameSorted`,
  # `NRF`, `PBC1`, `PBC2` — keep both via %||% fall-through for portability).
  out <- list(
    total_qnames        = as.integer(qc$totalQNAMEs                  %||% qc$totalQnameSorted %||% NA),
    duplicate_rate      = as.numeric(qc$duplicateRate                %||% NA),
    mitochondria_rate   = as.numeric(qc$mitochondriaRate             %||% NA),
    proper_pair_rate    = as.numeric(qc$properPairRate               %||% NA),
    unmapped_rate       = as.numeric(qc$unmappedRate                 %||% NA),
    has_unmapped_mate_rate = as.numeric(qc$hasUnmappedMateRate       %||% NA),
    not_passing_qc_rate    = as.numeric(qc$notPassingQualityControlsRate %||% NA),
    nrf                 = as.numeric(qc$nonRedundantFraction         %||% qc$NRF  %||% NA),
    pbc1                = as.numeric(qc$PCRbottleneckCoefficient_1   %||% qc$PBC1 %||% NA),
    pbc2                = as.numeric(qc$PCRbottleneckCoefficient_2   %||% qc$PBC2 %||% NA)
  )
  out_path <- file.path(golden_dir, paste0(sample, ".bamqc.golden.json"))
  write_json(out, out_path, auto_unbox = TRUE, digits = 15)
  message("  Written: ", out_path)
}

# ---- Helper: write TSV golden for fragSizeDist ------------------------------
write_fragsize_golden <- function(sample, bam_path) {
  message("fragSizeDist for ", sample, " ...")
  # fragSizeDist returns a named integer vector of counts per fragment size
  # ATACseqQC 1.34+ removed maxFragmentLength arg; pass only required args.
  # Return shape changed: now a named list whose element is a `table` keyed
  # by fragment size strings.
  fsd <- fragSizeDist(bam_path, bamFiles.labels = sample,
                      index = paste0(bam_path, ".bai"))
  tab <- if (is.list(fsd)) fsd[[1]] else fsd[, 1]
  sizes  <- as.integer(if (!is.null(names(tab))) names(tab) else rownames(tab))
  counts <- as.integer(tab)
  # Mirror rustqc's [1, 2000] cap (matches ATACseqQC's older `maxFragmentLength`
  # default; longer TLENs are mostly chimeric noise).
  keep <- sizes >= 1 & sizes <= 2000
  sizes  <- sizes[keep]
  counts <- counts[keep]
  out_path <- file.path(golden_dir, paste0(sample, ".fragsize.golden.tsv"))
  write.table(
    data.frame(frag_size = sizes, count = counts),
    file = out_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("  Written: ", out_path)
}

# ---- Helper: load shifted reads (needed for TSSE/NFR/PT) -------------------
load_shifted <- function(sample, bam_path) {
  message("Loading + shifting reads for ", sample, " ...")
  tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YZ")
  gal <- readBamFile(bam_path, tag = tags, which = chr1_full,
                     asMates = TRUE, bigFile = TRUE)
  tmp_shifted <- tempfile(fileext = ".bam")
  gal1 <- shiftGAlignmentsList(gal, outbam = tmp_shifted)
  gal1
}

# ---- Helper: write NFR golden -----------------------------------------------
write_nfr_golden <- function(sample, gal1) {
  message("NFRscore for ", sample, " ...")
  nfr <- NFRscore(gal1, txs)
  # ATACseqQC sorts the result by NFR_score descending; for byte-stable
  # comparison against rustqc (which emits TSS in GTF input order), we
  # re-sort by tx_name ascending so row N matches rustqc's tss_idx N.
  nfr_df <- as.data.frame(nfr)
  if ("tx_name" %in% colnames(nfr_df)) {
    nfr_df <- nfr_df[order(as.character(nfr_df$tx_name)), , drop = FALSE]
  }
  out_path <- file.path(golden_dir, paste0(sample, ".nfr.golden.tsv"))
  write.table(
    nfr_df,
    file = out_path, sep = "\t", quote = FALSE, row.names = TRUE
  )
  message("  Written: ", out_path)
}

# ---- Helper: write PT golden ------------------------------------------------
write_pt_golden <- function(sample, gal1) {
  message("PTscore for ", sample, " ...")
  pt <- PTscore(gal1, txs)
  pt_df <- as.data.frame(pt)
  if ("tx_name" %in% colnames(pt_df)) {
    pt_df <- pt_df[order(as.character(pt_df$tx_name)), , drop = FALSE]
  }
  out_path <- file.path(golden_dir, paste0(sample, ".pt.golden.tsv"))
  write.table(
    pt_df,
    file = out_path, sep = "\t", quote = FALSE, row.names = TRUE
  )
  message("  Written: ", out_path)
}

# ---- Helper: write TSSE golden ----------------------------------------------
write_tsse_golden <- function(sample, gal1) {
  message("TSSEscore for ", sample, " ...")
  tsse <- TSSEscore(gal1, txs)
  out <- list(
    tsse_score = as.numeric(tsse$TSSEscore),
    values     = as.numeric(tsse$values)
  )
  out_path <- file.path(golden_dir, paste0(sample, ".tsse.golden.json"))
  write_json(out, out_path, auto_unbox = TRUE, digits = 15)
  message("  Written: ", out_path)
}

# ---- Null-coalescing helper -------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---- Main loop over GL1-4 ---------------------------------------------------
samples <- c("GL1", "GL2", "GL3", "GL4")

for (sample in samples) {
  bam_path <- file.path(bam_dir, paste0(sample, ".bam"))
  if (!file.exists(bam_path)) {
    warning("BAM not found, skipping: ", bam_path)
    next
  }
  message("\n=== Processing ", sample, " ===")

  # bamQC (does not require shifted reads)
  tryCatch(write_bamqc_golden(sample, bam_path),  error = function(e) warning(e))

  # Fragment size distribution (does not require shifted reads)
  tryCatch(write_fragsize_golden(sample, bam_path), error = function(e) warning(e))

  # Metrics requiring shifted reads
  gal1 <- tryCatch(load_shifted(sample, bam_path), error = function(e) {
    warning("Could not load shifted reads for ", sample, ": ", e)
    NULL
  })

  if (!is.null(gal1)) {
    tryCatch(write_nfr_golden(sample, gal1),  error = function(e) warning(e))
    tryCatch(write_pt_golden(sample, gal1),   error = function(e) warning(e))
    tryCatch(write_tsse_golden(sample, gal1), error = function(e) warning(e))
  }
}

message("\nDone! Golden files written to: ", golden_dir)
message("Review output, then commit the resulting *.golden.{json,tsv} files.")
message("The Rust integration tests will auto-activate when goldens are present.")
