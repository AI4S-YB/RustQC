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
        // R's formula: max(c(1e-6, min(nf), min(n1), min(n2)))
        [1e-6, m_n1, m_nf, m_n2]
            .iter()
            .cloned()
            .filter(|x| x.is_finite())
            .fold(1e-6, f64::max)
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
        // small = max(1e-6, 1, 1, 1) = 1; NFR = log2(1+1)+1-log2(1+1+1) = 2-log2(3)
        let expected = 2.0_f64 - 3.0_f64.log2();
        assert!((r.nfr_score - expected).abs() < 1e-9);
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
        // n1=n2=1, nf=10; small = max(1e-6, min(n1)=1, min(nf)=10, min(n2)=1) = 10
        // NFR = log2(10+10)+1-log2(1+1+10) = log2(20)+1-log2(12)
        let expected = 20.0_f64.log2() + 1.0 - 12.0_f64.log2();
        assert!((r.nfr_score - expected).abs() < 1e-9);
    }
}
