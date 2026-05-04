//! Fragment-length histogram (1..1010 bp). Mirrors ATACseqQC R/fragSizeDist.R.

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FragSizeAccum {
    counts: [u64; 1011], // index 0 unused; valid index is 1..=1010
    total: u64,
}

impl Default for FragSizeAccum { fn default() -> Self { Self { counts: [0; 1011], total: 0 } } }

impl FragSizeAccum {
    #[allow(dead_code)]
    pub fn new() -> Self { Self::default() }

    /// Update from one record's TLEN (signed). Records out of [1,1010] after abs are dropped.
    #[allow(dead_code)]
    pub fn update(&mut self, tlen: i64) {
        let v = tlen.unsigned_abs();
        if v == 0 || v > 1010 { return; }
        self.counts[v as usize] += 1;
        self.total += 1;
    }

    /// Returns the (length, count, density) triples for length=1..=1010.
    #[allow(dead_code)]
    pub fn finalize(&self) -> Vec<(u32, u64, f64)> {
        let total = self.total.max(1) as f64;
        (1..=1010u32)
            .map(|l| {
                let c = self.counts[l as usize];
                (l, c, c as f64 / total)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_counts_abs_tlen_within_range() {
        let mut a = FragSizeAccum::new();
        for &t in &[150_i64, -150, 200, 200, -200, -200, 1011, 0, -1011] {
            a.update(t);
        }
        let h = a.finalize();
        assert_eq!(h[150 - 1].1, 2);   // length 150 → 2 records (one + one −)
        assert_eq!(h[200 - 1].1, 4);
        assert_eq!(h[1010 - 1].1, 0);   // 1011 dropped
        // Density sums to 1.
        let s: f64 = h.iter().map(|(_, _, d)| d).sum();
        assert!((s - 1.0).abs() < 1e-12);
    }
}
