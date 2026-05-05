//! Fragment-length histogram (1..2000 bp). Mirrors ATACseqQC R/fragSizeDist.R
//! (which historically used `maxFragmentLength=2000` as the default cap).

const FRAG_MAX: usize = 2000;

#[derive(Debug, Clone)]
pub struct FragSizeAccum {
    counts: [u64; FRAG_MAX + 1], // index 0 unused; valid index is 1..=FRAG_MAX
    total: u64,
}

impl Default for FragSizeAccum {
    fn default() -> Self { Self { counts: [0; FRAG_MAX + 1], total: 0 } }
}

impl FragSizeAccum {
    pub fn new() -> Self { Self::default() }

    /// Update from one record's TLEN (signed). Records out of [1, FRAG_MAX]
    /// after abs are dropped (matches ATACseqQC's fragment-length cap).
    pub fn update(&mut self, tlen: i64) {
        let v = tlen.unsigned_abs() as usize;
        if v == 0 || v > FRAG_MAX { return; }
        self.counts[v] += 1;
        self.total += 1;
    }

    /// Returns the (length, count, density) triples for length=1..=FRAG_MAX.
    pub fn finalize(&self) -> Vec<(u32, u64, f64)> {
        let total = self.total.max(1) as f64;
        (1..=FRAG_MAX as u32)
            .map(|l| {
                let c = self.counts[l as usize];
                (l, c, c as f64 / total)
            })
            .collect()
    }
}

pub fn write_tsv<W: std::io::Write>(w: &mut W, h: &[(u32, u64, f64)]) -> std::io::Result<()> {
    writeln!(w, "length\tcount\tnorm_density")?;
    for (l, c, d) in h {
        writeln!(w, "{}\t{}\t{:.10e}", l, c, d)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_counts_abs_tlen_within_range() {
        let mut a = FragSizeAccum::new();
        for &t in &[150_i64, -150, 200, 200, -200, -200, 2001, 0, -2001] {
            a.update(t);
        }
        let h = a.finalize();
        assert_eq!(h[150 - 1].1, 2);   // length 150 → 2 records (one + one −)
        assert_eq!(h[200 - 1].1, 4);
        assert_eq!(h[2000 - 1].1, 0);  // length 2001 dropped (above cap)
        // Density sums to 1.
        let s: f64 = h.iter().map(|(_, _, d)| d).sum();
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn tsv_format_matches_spec() {
        let mut a = FragSizeAccum::new();
        a.update(100); a.update(150); a.update(150);
        let mut buf = Vec::new();
        write_tsv(&mut buf, &a.finalize()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("length\tcount\tnorm_density\n"));
        let line_100 = s.lines().nth(100).unwrap();   // header + length 1..100 → index 100
        assert!(line_100.starts_with("100\t1\t"));
        let line_150 = s.lines().nth(150).unwrap();
        assert!(line_150.starts_with("150\t2\t"));
    }
}
