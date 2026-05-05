//! TSSEscore. Mirrors ATACseqQC R/TSSEscore.R.

use crate::atac::loess::loess_smooth;
use crate::atac::tss_cov::TssCov;
use crate::gtf::Strand;

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
    assert!(
        flank >= TSSE_FLANK + END_SIZE,
        "TssCov flank must be >= {} for TSSE (got {})",
        TSSE_FLANK + END_SIZE,
        flank
    );
    let n_windows = (2 * TSSE_FLANK) / WIDTH;
    let center_lo = flank - TSSE_FLANK;
    let center_hi = flank + TSSE_FLANK;
    let mut sums = vec![0.0f64; n_windows];
    let mut surviving = 0u64;
    for (i, buf) in cov.buffers.iter().enumerate() {
        let strand = cov.tss_list[i].strand;
        let mean = |range: std::ops::Range<usize>| -> f64 {
            let mut s = 0u64; let mut n = 0u64;
            for b in range { s += buf[b] as u64; n += 1; }
            s as f64 / n as f64
        };
        // Background flanks per ATACseqQC: 100 bp OUTSIDE the [TSS-1000, TSS+999]
        // sliding-window region (the R `flank()` operator on `sel.center`).
        let vl = mean(center_lo - END_SIZE..center_lo);
        let vr = mean(center_hi..center_hi + END_SIZE);
        let blk = (vl + vr) / 2.0;
        if blk <= 0.0 { continue; }
        // R's TSSEscore averages bins in GENOMIC order (slidingWindows on the
        // promoters() range, which is not strand-aware after the slide). Our
        // TssCov buffer is strand-mirrored so that bin index `flank` is the TSS
        // for both strands; for + strand, transcription-upstream → window 0;
        // for − strand the same mirroring puts transcription-upstream at
        // window 0 too — but R's bin 0 for − strand is genomic-leftmost
        // (= transcription-downstream). Reverse the accumulation index for
        // − strand to match R's strand-mixed colMeans.
        for w in 0..n_windows {
            let lo = center_lo + w * WIDTH;
            let v = mean(lo..lo + WIDTH);
            let target_w = match strand {
                Strand::Plus => w,
                Strand::Minus => n_windows - 1 - w,
            };
            sums[target_w] += v / blk;
        }
        surviving += 1;
    }
    let s = surviving.max(1) as f64;
    let raw: Vec<f64> = sums.iter().map(|x| x / s).collect();
    let xs: Vec<f64> = (1..=n_windows).map(|i| i as f64).collect();
    // R's loess.smooth() defaults to degree=1; ATACseqQC's TSSEscore uses
    // those defaults (no `degree` override), so match degree=1 here.
    let smoothed = loess_smooth(&xs, &raw, 2.0/3.0, 1);
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
        let mut cov = TssCov::new(tss, 1100);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        let r = compute(&cov);
        assert!((r.tsse_score - 1.0).abs() < 1e-3, "score={}", r.tsse_score);
        assert_eq!(r.values.len(), 20);
    }

    #[test]
    fn central_enrichment_lifts_score_above_baseline() {
        let tss = vec![Tss { chrom: "chr1".into(), pos: 1_000_000, strand: Strand::Plus }];
        let mut cov = TssCov::new(tss, 1100);
        for b in 0..cov.buffers[0].len() { cov.buffers[0][b] = 1; }
        // 5× boost in the central 200 bp (TSS at bin = flank = 1100).
        for b in 1000..1200 { cov.buffers[0][b] = 5; }
        let r = compute(&cov);
        assert!(r.tsse_score > 1.5, "expected enrichment, got {}", r.tsse_score);
    }
}
