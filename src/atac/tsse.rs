//! TSSEscore. Mirrors ATACseqQC R/TSSEscore.R.

use crate::atac::loess::loess_smooth;
use crate::atac::tss_cov::TssCov;

const TSSE_FLANK: usize = 1000;
const END_SIZE: usize = 100;
const WIDTH: usize = 100;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TsseResult {
    pub values: Vec<f64>,    // smoothed, length = 2*TSSE_FLANK / WIDTH = 20
    pub tsse_score: f64,
}

#[allow(dead_code)]
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
