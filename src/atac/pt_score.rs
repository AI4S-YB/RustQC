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
    assert!(flank >= 3000, "PTscore needs flank >= 3000 (got {})", flank);
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
        // R's formula: max(c(1e-6, min(promoter), min(body)))
        [1e-6, min_finite(&raw, |t| t.0), min_finite(&raw, |t| t.1)]
            .iter()
            .cloned()
            .filter(|x| x.is_finite())
            .fold(1e-6, f64::max)
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
        // promoter=4, body=1; small = max(1e-6, min(prom)=4, min(body)=1) = 4
        // PT = log2(4+4)-log2(1+4) = log2(8)-log2(5)
        let expected = 8.0_f64.log2() - 5.0_f64.log2();
        assert!((r.pt_score - expected).abs() < 1e-9);
    }
}
