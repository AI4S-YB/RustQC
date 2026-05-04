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
})

# ---- Paths ------------------------------------------------------------------
script_dir  <- dirname(normalizePath(if (interactive()) "." else commandArgs(trailingOnly = FALSE)[4]))
repo_root   <- file.path(script_dir, "..", "..", "..")
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
txdb_local <- makeTxDbFromGFF(gtf_path, format = "GTF")
txs        <- transcripts(txdb_local)
seqlevels(txs, pruning.mode = "coarse") <- "chr1"
message("  Loaded ", length(txs), " transcripts on chr1")

# ---- Helper: write JSON golden for bamQC ------------------------------------
write_bamqc_golden <- function(sample, bam_path) {
  message("bamQC for ", sample, " ...")
  qc <- bamQC(bam_path, outPath = NULL)
  # ATACseqQC bamQC returns a list; key fields:
  #   duplicateRate, mitochondriaRate, properPairRate, NRF, PBC1, PBC2,
  #   totalQnameSorted
  out <- list(
    total_qnames      = as.integer(qc$totalQnameSorted %||% NA),
    duplicate_rate    = as.numeric(qc$duplicateRate    %||% NA),
    mitochondria_rate = as.numeric(qc$mitochondriaRate  %||% NA),
    proper_pair_rate  = as.numeric(qc$properPairRate    %||% NA),
    nrf               = as.numeric(qc$NRF               %||% NA),
    pbc1              = as.numeric(qc$PBC1              %||% NA),
    pbc2              = as.numeric(qc$PBC2              %||% NA)
  )
  out_path <- file.path(golden_dir, paste0(sample, ".bamqc.golden.json"))
  write_json(out, out_path, auto_unbox = TRUE, digits = 15)
  message("  Written: ", out_path)
}

# ---- Helper: write TSV golden for fragSizeDist ------------------------------
write_fragsize_golden <- function(sample, bam_path) {
  message("fragSizeDist for ", sample, " ...")
  # fragSizeDist returns a named integer vector of counts per fragment size
  fsd <- fragSizeDist(bam_path, bamFiles.labels = sample,
                      maxFragmentLength = 2000, index = paste0(bam_path, ".bai"))
  # fsd is a matrix with one column per sample; rows = fragment sizes
  sizes  <- as.integer(rownames(fsd))
  counts <- as.integer(fsd[, 1])
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
  gal <- readBamFile(bam_path, tag = tags, which = as(seqinfo(txdb_local)["chr1"], "GRanges"),
                     asMates = TRUE, bigFile = TRUE)
  tmp_shifted <- tempfile(fileext = ".bam")
  gal1 <- shiftGAlignmentsList(gal, outbam = tmp_shifted)
  gal1
}

# ---- Helper: write NFR golden -----------------------------------------------
write_nfr_golden <- function(sample, gal1) {
  message("NFRscore for ", sample, " ...")
  nfr <- NFRscore(gal1, txs)
  out_path <- file.path(golden_dir, paste0(sample, ".nfr.golden.tsv"))
  write.table(
    as.data.frame(nfr),
    file = out_path, sep = "\t", quote = FALSE, row.names = TRUE
  )
  message("  Written: ", out_path)
}

# ---- Helper: write PT golden ------------------------------------------------
write_pt_golden <- function(sample, gal1) {
  message("PTscore for ", sample, " ...")
  pt <- PTscore(gal1, txs)
  out_path <- file.path(golden_dir, paste0(sample, ".pt.golden.tsv"))
  write.table(
    as.data.frame(pt),
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
