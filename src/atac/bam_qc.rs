//! ATACseqQC-style bamQC: rates, NRF, PBC1/2, MAPQ histogram.
//!
//! Numerical fidelity to ATACseqQC 1.36.0 R/bamQC.R is required. See
//! docs/superpowers/specs/2026-05-04-atac-seq-qc-design.md §"bamQC".

use std::collections::{HashMap, HashSet};

#[allow(dead_code)]
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
    /// Unique QNAME set; `len()` yields ATACseqQC's `totalQNAMEs` value used as the
    /// NRF denominator. Grows proportional to distinct reads — see Phase 13 handoff
    /// notes for the eventual streaming-driver memory tradeoff.
    pub qnames: HashSet<String>,
}

#[allow(dead_code)]
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
    pub mapq_hist: Vec<(u8, u64)>, // sorted ascending by mapq
}

#[allow(dead_code)]
impl BamQcAccum {
    pub fn new() -> Self {
        Self::default()
    }

    /// Update from a single primary record's flag bits and MAPQ.
    /// Caller decides which records to feed; ATACseqQC excludes secondary alignments
    /// (same `isSecondaryAlignment = FALSE` flag we use).
    pub fn update_flags(&mut self, flags: u16, mapq: u8, is_mito: bool, qname: &str) {
        const F_PROPER_PAIR: u16 = 0x2;
        const F_UNMAPPED: u16 = 0x4;
        const F_MATE_UNMAPPED: u16 = 0x8;
        const F_DUP: u16 = 0x400;
        const F_QCFAIL: u16 = 0x200;
        self.total_records += 1;
        if flags & F_DUP != 0 {
            self.n_dup += 1;
        }
        if flags & F_PROPER_PAIR != 0 {
            self.n_proper_pair += 1;
        }
        if flags & F_UNMAPPED != 0 {
            self.n_unmapped += 1;
        }
        if flags & F_MATE_UNMAPPED != 0 {
            self.n_unmapped_mate += 1;
        }
        if flags & F_QCFAIL != 0 {
            self.n_qc_fail += 1;
        }
        if is_mito {
            self.n_mito += 1;
        }
        *self.mapq_hist.entry(mapq).or_default() += 1;
        self.qnames.insert(qname.to_string());
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct PbcChromAccum {
    /// Fingerprint → count.
    /// Key is a `(pos1, isize1, pos2, isize2)` tuple (all i64) to safely handle
    /// genomes whose coordinates exceed u32::MAX without truncation.
    pub fingerprints: HashMap<(i64, i64, i64, i64), u64>,
}

#[allow(dead_code)]
impl PbcChromAccum {
    /// Record one PE fragment fingerprint.
    /// `pos1`/`isize1` come from mate 1; `pos2`/`isize2` from mate 2.
    /// For singletons, pass sentinel values (e.g. i64::MIN) for the missing mate.
    pub fn add_pe(&mut self, pos1: i64, isize1: i64, pos2: i64, isize2: i64) {
        *self.fingerprints.entry((pos1, isize1, pos2, isize2)).or_default() += 1;
    }

    /// Returns `(M_DISTINCT, M1, M2)` — used in the aggregate to compute NRF/PBC1/PBC2.
    ///
    /// - `M_DISTINCT`: number of unique fingerprints.
    /// - `M1`: fingerprints occurring exactly once.
    /// - `M2`: fingerprints occurring exactly twice.
    pub fn summarize(&self) -> (u64, u64, u64) {
        let m_distinct = self.fingerprints.len() as u64;
        let m1 = self.fingerprints.values().filter(|&&c| c == 1).count() as u64;
        let m2 = self.fingerprints.values().filter(|&&c| c == 2).count() as u64;
        (m_distinct, m1, m2)
    }
}

/// Aggregate flag counters and per-chromosome PBC fingerprints into a `BamQcReport`.
///
/// Formulas match ATACseqQC 1.36.0 `R/bamQC.R`:
/// - `NRF  = ΣM1 / totalQNAMEs`  (0.0 when totalQNAMEs == 0)
/// - `PBC1 = ΣM1 / ΣM_DISTINCT`  (0.0 when ΣM_DISTINCT == 0)
/// - `PBC2 = ΣM1 / max(1, ΣM2)`
/// - `mapq_hist` sorted ascending by MAPQ value.
#[allow(dead_code)]
pub fn finalize(flag_acc: &BamQcAccum, pbc_per_chrom: &[PbcChromAccum]) -> BamQcReport {
    let total = flag_acc.total_records.max(1) as f64; // avoid div0 — total > 0 in real runs
    let (mut sum_distinct, mut sum_m1, mut sum_m2) = (0u64, 0u64, 0u64);
    for p in pbc_per_chrom {
        let (md, m1, m2) = p.summarize();
        sum_distinct += md;
        sum_m1 += m1;
        sum_m2 += m2;
    }
    let total_qnames = flag_acc.qnames.len() as u64;
    let nrf = if total_qnames == 0 {
        0.0
    } else {
        sum_m1 as f64 / total_qnames as f64
    };
    let pbc1 = if sum_distinct == 0 {
        0.0
    } else {
        sum_m1 as f64 / sum_distinct as f64
    };
    let pbc2 = sum_m1 as f64 / sum_m2.max(1) as f64;

    let mut mapq_hist: Vec<(u8, u64)> = flag_acc
        .mapq_hist
        .iter()
        .map(|(k, v)| (*k, *v))
        .collect();
    mapq_hist.sort_by_key(|(k, _)| *k);

    BamQcReport {
        total_qnames,
        duplicate_rate: flag_acc.n_dup as f64 / total,
        mitochondria_rate: flag_acc.n_mito as f64 / total,
        proper_pair_rate: flag_acc.n_proper_pair as f64 / total,
        unmapped_rate: flag_acc.n_unmapped as f64 / total,
        has_unmapped_mate_rate: flag_acc.n_unmapped_mate as f64 / total,
        not_passing_qc_rate: flag_acc.n_qc_fail as f64 / total,
        nrf,
        pbc1,
        pbc2,
        mapq_hist,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flag_aggregation_matches_r() {
        let mut a = BamQcAccum::new();
        // 4 records: 1 mito+dup, 1 proper-pair, 1 qc-fail, 1 unmapped-mate.
        a.update_flags(0x402, 30, true, "r1"); // dup + proper_pair
        a.update_flags(0x002, 60, false, "r2"); // proper_pair
        a.update_flags(0x200, 0, false, "r3"); // qcfail, mapq 0
        a.update_flags(0x008, 30, false, "r4"); // mate unmapped
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

    #[test]
    fn pbc_summarize_counts_singletons_and_doubletons() {
        let mut p = PbcChromAccum::default();
        p.add_pe(100, 200, 100, -200);
        p.add_pe(100, 200, 100, -200); // duplicate of above
        p.add_pe(300, 200, 300, -200); // singleton
        p.add_pe(500, 200, 500, -200); // singleton
        p.add_pe(500, 200, 500, -200); // doubleton with above
        let (m_distinct, m1, m2) = p.summarize();
        assert_eq!(m_distinct, 3);
        assert_eq!(m1, 1);
        assert_eq!(m2, 2);
    }

    #[test]
    #[allow(non_snake_case)]
    fn finalize_computes_NRF_PBC1_PBC2() {
        let mut flag = BamQcAccum::new();
        for i in 0..10 {
            flag.qnames.insert(format!("r{}", i));
        }
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

    #[test]
    fn finalize_pbc2_does_not_divide_by_zero_when_m2_is_zero() {
        // Spec edge case: when no fingerprint occurs exactly twice (M2 == 0),
        // PBC2 = M1 / max(1, M2) must yield M1 / 1, not panic.
        let mut flag = BamQcAccum::new();
        flag.qnames.insert("r1".into());
        flag.total_records = 1;
        let mut p = PbcChromAccum::default();
        p.fingerprints.insert((1, 2, 3, 4), 1); // singleton only — M1=1, M2=0
        let r = finalize(&flag, &[p]);
        assert!((r.pbc2 - 1.0).abs() < 1e-12);
    }
}
