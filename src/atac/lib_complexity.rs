//! Library complexity: readsDupFreq → preseq bootstrap SAC curve.
//!
//! `DupFreqAccum` tracks PE fragment fingerprints `(chrom_id, leftpos, isize)`
//! and emits a histogram of multiplicity counts used by the preseq adapter.

use std::collections::HashMap;

/// Accumulates paired-end fragment fingerprints and produces a
/// frequency-of-frequencies histogram for library complexity estimation.
///
/// Fingerprint key: `(chrom_id, leftmost_pos, isize)` — matches ATACseqQC's
/// `readsDupFreq` which keys on chromosome + leftmost alignment position +
/// insert size.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct DupFreqAccum {
    /// (chrom_id, leftmost_pos, isize) → observation count
    pub fingerprints: HashMap<(u32, i64, i64), u64>,
}

#[allow(dead_code)]
impl DupFreqAccum {
    /// Record one observation of the PE fragment defined by
    /// `(chrom_id, leftpos, isize)`.
    pub fn add_pe(&mut self, chrom_id: u32, leftpos: i64, isize: i64) {
        *self.fingerprints.entry((chrom_id, leftpos, isize)).or_default() += 1;
    }

    /// Build histogram rows: `Vec<(j, n_j)>` sorted by `j` ascending,
    /// where `n_j` is the number of distinct fingerprints observed exactly `j`
    /// times.
    pub fn histogram(&self) -> Vec<(u64, u64)> {
        let mut by_j: HashMap<u64, u64> = HashMap::new();
        for &c in self.fingerprints.values() {
            *by_j.entry(c).or_default() += 1;
        }
        let mut v: Vec<(u64, u64)> = by_j.into_iter().collect();
        v.sort_by_key(|(j, _)| *j);
        v
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_counts_multiplicities() {
        let mut a = DupFreqAccum::default();
        // 5 distinct fingerprints with multiplicities 1, 1, 2, 3, 3
        a.add_pe(0, 100, 200);
        a.add_pe(0, 200, 200);
        a.add_pe(0, 300, 200);
        a.add_pe(0, 300, 200);
        a.add_pe(0, 400, 200);
        a.add_pe(0, 400, 200);
        a.add_pe(0, 400, 200);
        a.add_pe(0, 500, 200);
        a.add_pe(0, 500, 200);
        a.add_pe(0, 500, 200);
        let h = a.histogram();
        assert_eq!(h, vec![(1, 2), (2, 1), (3, 2)]);
    }
}
