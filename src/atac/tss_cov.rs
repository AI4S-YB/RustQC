//! Sparse per-TSS 5'-end coverage. Underlies TSSEscore, NFRscore, PTscore.

use crate::gtf::{Strand, Tss};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct TssCov {
    pub flank: u32,                            // half-window in bp; arrays have length 2*flank
    pub buffers: Vec<Vec<u32>>,                // index = TSS index in `tss_list`
    pub tss_list: Vec<Tss>,
    by_chrom: HashMap<String, Vec<usize>>,     // chrom → indices into tss_list
}

impl TssCov {
    pub fn new(tss_list: Vec<Tss>, flank: u32) -> Self {
        let mut by_chrom: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, t) in tss_list.iter().enumerate() {
            by_chrom.entry(t.chrom.clone()).or_default().push(i);
        }
        // Buffer length 2*flank + 1: bin index `flank` is the TSS itself (offset 0),
        // bin 0 is TSS-flank, bin 2*flank is TSS+flank — for both strands.
        let buf_len = (2 * flank + 1) as usize;
        let buffers = tss_list.iter().map(|_| vec![0u32; buf_len]).collect();
        Self { flank, buffers, tss_list, by_chrom }
    }

    /// Increment the bin under a read's 5' position if it falls within any TSS window
    /// on this chromosome. `pos5p` is 1-based (BAM coordinate convention).
    ///
    /// Bin indexing places the TSS at bin index `flank` for both strands.
    ///
    /// - `+` strand: `bin = pos5p - (TSS - flank)`, `∈ [0, 2*flank]`
    /// - `−` strand: `bin = (TSS + flank) - pos5p`, `∈ [0, 2*flank]`
    ///
    /// (mirror = `2*flank − bin_raw`, no −1; matches ATACseqQC's symmetric
    /// `promoters()` window definitions for both strands.)
    pub fn add_5prime(&mut self, chrom: &str, pos5p: u64) {
        let Some(idxs) = self.by_chrom.get(chrom) else { return; };
        for &i in idxs {
            let t = &self.tss_list[i];
            let win_start = t.pos.saturating_sub(self.flank as u64);
            let win_end = t.pos + self.flank as u64;
            if pos5p < win_start || pos5p > win_end { continue; }
            let bin_raw = (pos5p - win_start) as usize;
            let bin = match t.strand {
                Strand::Plus => bin_raw,
                Strand::Minus => (2 * self.flank as usize) - bin_raw,
            };
            self.buffers[i][bin] = self.buffers[i][bin].saturating_add(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tss_at(chrom: &str, pos: u64, strand: Strand) -> Tss {
        Tss { chrom: chrom.into(), pos, strand }
    }

    #[test]
    fn coverage_strand_aware() {
        let tss = vec![
            tss_at("chr1", 1000, Strand::Plus),
            tss_at("chr1", 5000, Strand::Minus),
        ];
        let mut c = TssCov::new(tss, 100);
        c.add_5prime("chr1", 1050);   // 50 bp downstream of + TSS in genomic = bin flank+50=150
        c.add_5prime("chr1", 4990);   // 10 bp upstream of − TSS (genomic) = bin 2*flank-90=110
        // + TSS bin: 1050 - (1000-100) = 150
        assert_eq!(c.buffers[0][150], 1);
        // − TSS: pos=5000, flank=100. bin_raw = 4990 - (5000-100) = 90.
        // mirror = 2*flank - bin_raw = 200 - 90 = 110.
        assert_eq!(c.buffers[1][110], 1);
    }
}
