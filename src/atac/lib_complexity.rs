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
pub struct DupFreqAccum {
    /// (chrom_id, leftmost_pos, isize) → observation count
    pub fingerprints: HashMap<(u32, i64, i64), u64>,
}

impl DupFreqAccum {
    /// Record one observation of the PE fragment defined by
    /// `(chrom_id, leftpos, isize)`.
    pub fn add_pe(&mut self, chrom_id: u32, leftpos: i64, isize: i64) {
        *self
            .fingerprints
            .entry((chrom_id, leftpos, isize))
            .or_default() += 1;
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
// Library complexity estimation
// ============================================================================

/// One row of the library complexity curve, corresponding to one sample size.
#[derive(Debug, Clone)]
pub struct LibComplexityRow {
    /// Fraction of observed total reads (e.g. 0.1 = 10 %, 5.0 = 500 %).
    pub relative_size: f64,
    /// Expected number of distinct fragments at this sample size.
    pub distinct_fragments: f64,
    /// Putative number of reads at this sample size (`relative_size × total`).
    pub putative_reads: f64,
}

/// Evaluate the library complexity curve at 14 standard ATACseqQC sample sizes.
///
/// Sample sizes (relative_size): `{0.1, 0.2, …, 1.0, 5, 10, 15, 20}`.
///
/// Matches `ATACseqQC::estimateLibComplexity` which calls
/// `preseqR::ds.rSAC.bootstrap(hist, r=1, times=<times>)` then evaluates at
/// those 14 relative sizes.
///
/// # Arguments
/// * `hist`  – Frequency-of-frequencies `(j, n_j)` from `DupFreqAccum::histogram`.
/// * `times` – Number of bootstrap replicates (100 in production; use ≥50).
///
/// # Returns
/// 14 rows in sample-size order.  `distinct_fragments` may be `NaN` for rows
/// where bootstrap convergence failed (very small or degenerate histograms).
pub fn estimate(hist: &[(u64, u64)], times: u32) -> anyhow::Result<Vec<LibComplexityRow>> {
    // total = Σ j·n_j  (matches R: histFile[,1] %*% histFile[,2])
    let total: u64 = hist.iter().map(|(j, n)| j * n).sum();
    let n_distinct: u64 = hist.iter().map(|(_, n)| n).sum();

    // 14 standard relative sizes
    let relative_sizes: Vec<f64> = (1u32..=10)
        .map(|i| i as f64 * 0.1)
        .chain([5.0, 10.0, 15.0, 20.0])
        .collect();

    // Convert relative sizes → absolute read counts for preseq
    let targets: Vec<f64> = relative_sizes.iter().map(|&s| s * total as f64).collect();

    // Bootstrap SAC curve at our explicit targets.
    // We use seed 408 (matching upstream preseq default) for reproducibility.
    let estimates =
        crate::preseq::estimate_at_targets(hist, total, n_distinct, &targets, times, 408)?;

    let rows = relative_sizes
        .into_iter()
        .zip(targets)
        .zip(estimates)
        .map(|((rel, abs_reads), distinct)| LibComplexityRow {
            relative_size: rel,
            distinct_fragments: distinct,
            putative_reads: abs_reads,
        })
        .collect();

    Ok(rows)
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

    #[test]
    fn estimate_returns_14_rows() {
        let hist = vec![(1u64, 100u64), (2, 50), (3, 20), (4, 5)];
        let rows = estimate(&hist, 50).unwrap();
        assert_eq!(rows.len(), 14);
        assert!(
            (rows[0].relative_size - 0.1).abs() < 1e-12,
            "first row relative_size should be 0.1, got {}",
            rows[0].relative_size
        );
        assert_eq!(
            rows[13].relative_size, 20.0,
            "last row relative_size should be 20.0"
        );
        // Sanity: distinct_fragments should be finite and non-negative for all rows
        for row in &rows {
            assert!(
                row.distinct_fragments.is_finite() && row.distinct_fragments >= 0.0,
                "distinct_fragments={} at relative_size={} is not finite/non-negative",
                row.distinct_fragments,
                row.relative_size
            );
        }
    }
}
