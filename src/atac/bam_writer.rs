//! BAM emission helpers for opt-in Tn5 shift output.
//!
//! Phase 12 delivers:
//!   - [`trim_cigar_5prime`]      — CIGAR 5'-trim helper (Task 12.1)
//!   - [`rewrite_record_inplace`] — Tn5 +4/−5 shift applied to a `RecordBuf` (Task 12.2)
//!   - [`EmitWriters`]            — scaffold multiplexer stub (Task 12.3; Phase 14 fills it)

#[allow(unused_imports)]
pub(crate) use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles_sam::alignment::record::cigar::Op as CigarOp;

// ─────────────────────────────────────────────────────────────────────────────
// Task 12.1 — trim_cigar_5prime
// ─────────────────────────────────────────────────────────────────────────────

/// Trim `n` query-consuming bases from the 5'-end of a CIGAR operation slice.
///
/// Returns `(new_ops, ref_shift)` where `ref_shift` is the number of reference-
/// genome bases consumed by the trimmed portion — i.e. the amount by which the
/// alignment start (POS) must be advanced on the + strand.
///
/// Algorithm:
/// - Walk left-to-right.
/// - Non-query-consuming ops (D/N/H/P) encountered while `remaining > 0` are
///   **dropped** and their reference length is added to `shift`; this matches
///   the `cigarQNarrow` convention in R's ATACseqQC.
/// - Query-consuming ops are split at the trim boundary.
#[allow(dead_code)]
pub fn trim_cigar_5prime(ops: &[CigarOp], n: u32) -> (Vec<CigarOp>, u32) {
    let mut remaining = n as usize;
    let mut shift: u32 = 0;
    let mut out: Vec<CigarOp> = Vec::with_capacity(ops.len());

    let mut iter = ops.iter().copied();
    while remaining > 0 {
        let op = match iter.next() {
            Some(o) => o,
            None => break, // exhausted — empty result
        };
        let kind = op.kind();
        let len = op.len();

        let consumes_query = kind.consumes_read();
        let consumes_ref = kind.consumes_reference();

        if !consumes_query {
            // D/N/H/P: no query consumed; drop and advance reference shift.
            if consumes_ref {
                shift += len as u32;
            }
            continue;
        }

        // Query-consuming op.
        if len <= remaining {
            // Fully consumed by trim.
            remaining -= len;
            if consumes_ref {
                shift += len as u32;
            }
        } else {
            // Partially consumed: emit the tail, stop trimming.
            let leftover = len - remaining;
            if consumes_ref {
                shift += remaining as u32;
            }
            remaining = 0;
            out.push(CigarOp::new(kind, leftover));
        }
    }

    // Emit any remaining ops after the trim boundary.
    out.extend(iter);

    (out, shift)
}

// ─────────────────────────────────────────────────────────────────────────────
// Task 12.2 — rewrite_record_inplace
// ─────────────────────────────────────────────────────────────────────────────

/// Apply the Tn5 +4/−5 read-shift to a [`noodles_sam::alignment::RecordBuf`] in-place.
///
/// - `is_plus`: `true` → trim 4 bases from the 5'-end (left end) and advance POS by the
///   reference shift; `false` → trim 5 bases from the 3'-end (right end in read
///   coordinates, i.e. the tail of the stored sequence).
///
/// Returns `Ok(true)` if the record was modified successfully, `Ok(false)` if the record
/// should be dropped (read too short after shift, or TLEN would become degenerate ≤ 0).
///
/// # noodles 0.81 API notes
///
/// All required mutators are available:
/// - `rec.alignment_start_mut()` → `&mut Option<Position>`
/// - `rec.cigar_mut()`           → `&mut record_buf::Cigar`  (wraps `Vec<Op>`)
/// - `rec.sequence_mut()`        → `&mut Sequence`           (`AsMut<Vec<u8>>` via `as_mut()`)
/// - `rec.quality_scores_mut()`  → `&mut QualityScores`      (`AsMut<Vec<u8>>` via `as_mut()`)
/// - `rec.template_length_mut()` → `&mut i32`
#[allow(dead_code)]
pub fn rewrite_record_inplace(
    rec: &mut noodles_sam::alignment::RecordBuf,
    is_plus: bool,
) -> anyhow::Result<bool> {
    use noodles_core::Position;
    use noodles_sam::alignment::record_buf::Cigar as RecordCigar;

    let n: usize = if is_plus { 4 } else { 5 };

    // ── 1. Collect CIGAR ops ──────────────────────────────────────────────────
    let ops_in: Vec<CigarOp> = rec.cigar().as_ref().to_vec();

    // ── 2. Trim CIGAR ─────────────────────────────────────────────────────────
    // + strand: trim from the left (5' of the read = left of CIGAR).
    // − strand: reverse ops, trim 5 from the left, reverse back.
    //           ref_shift is 0 for − strand (POS does not move; the *tail* shrinks).
    let (new_ops, ref_shift) = if is_plus {
        trim_cigar_5prime(&ops_in, n as u32)
    } else {
        let rev: Vec<CigarOp> = ops_in.iter().rev().copied().collect();
        let (mut trimmed, _) = trim_cigar_5prime(&rev, n as u32);
        trimmed.reverse();
        (trimmed, 0u32)
    };

    if new_ops.is_empty() {
        return Ok(false);
    }

    // ── 3. SEQ + QUAL sanity check ────────────────────────────────────────────
    let seq_len = rec.sequence().len();
    if seq_len <= n {
        return Ok(false);
    }
    let new_seq_len = seq_len - n;

    // ── 4. POS shift (+ strand only) ─────────────────────────────────────────
    if is_plus && ref_shift > 0 {
        if let Some(start) = rec.alignment_start() {
            let new_pos_1based = usize::from(start) + ref_shift as usize;
            // Position is 1-based; usize::from(Position) yields the 1-based value.
            let new_pos = Position::new(new_pos_1based).ok_or_else(|| {
                anyhow::anyhow!(
                    "alignment start overflow after Tn5 shift: {}",
                    new_pos_1based
                )
            })?;
            *rec.alignment_start_mut() = Some(new_pos);
        }
    }

    // ── 5. SEQ trim ───────────────────────────────────────────────────────────
    {
        let seq_vec: &mut Vec<u8> = rec.sequence_mut().as_mut();
        if is_plus {
            seq_vec.drain(..n);
        } else {
            seq_vec.truncate(new_seq_len);
        }
    }

    // ── 6. QUAL trim ─────────────────────────────────────────────────────────
    {
        let qual_vec: &mut Vec<u8> = rec.quality_scores_mut().as_mut();
        if !qual_vec.is_empty() {
            // QUAL may legitimately be absent (empty) in some BAM records.
            if qual_vec.len() >= n {
                if is_plus {
                    qual_vec.drain(..n);
                } else {
                    qual_vec.truncate(new_seq_len);
                }
            } else {
                qual_vec.clear();
            }
        }
    }

    // ── 7. Replace CIGAR ─────────────────────────────────────────────────────
    *rec.cigar_mut() = RecordCigar::from(new_ops);

    // ── 8. TLEN adjustment ───────────────────────────────────────────────────
    // ATACseqQC adjusts |TLEN| by 9 = 4 + 5 (total shift across both mates).
    let tlen = rec.template_length();
    if tlen != 0 {
        let abs = tlen.unsigned_abs();
        if abs <= 9 {
            return Ok(false);
        }
        let new_tlen = (abs - 9) as i32 * tlen.signum();
        *rec.template_length_mut() = new_tlen;
    }

    Ok(true)
}

// ─────────────────────────────────────────────────────────────────────────────
// Task 12.3 — EmitWriters scaffold
// ─────────────────────────────────────────────────────────────────────────────

use std::path::Path;

/// Multiplexer for optional shifted/split BAM output streams.
///
/// Each field is `None` until Phase 14's `open()` implementation creates the
/// underlying file and BAI-indexed writer.  For now the scaffold compiles cleanly
/// and lets Phase 13's driver wiring reference the type.
#[allow(dead_code)]
#[derive(Default)]
pub struct EmitWriters {
    pub shifted: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub nfr: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub mono: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub di: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
    pub tri: Option<noodles_bam::io::Writer<std::io::BufWriter<std::fs::File>>>,
}

impl EmitWriters {
    /// Open output BAM writers for the requested emission modes.
    ///
    /// **Phase 14 placeholder** — currently returns an all-`None` default so that
    /// Phase 13's driver wiring compiles without requiring real file I/O.
    /// Phase 14 will add file creation, header writing, and BAI-index scaffolding.
    #[allow(dead_code)]
    pub fn open(
        outdir: &Path,
        sample: &str,
        emit_shifted: bool,
        emit_split: bool,
        header: &noodles_sam::Header,
    ) -> anyhow::Result<Self> {
        // TODO(phase-14): create files, write SAM headers, set up BAI indexing.
        let _ = (outdir, sample, emit_shifted, emit_split, header);
        Ok(Self::default())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use CigarKind::{Insertion, Match, SoftClip};

    fn op(kind: CigarKind, len: usize) -> CigarOp {
        CigarOp::new(kind, len)
    }

    // ── Task 12.1 tests ───────────────────────────────────────────────────────

    /// [M50] trim 4 → [M46], shift=4.
    #[test]
    fn trim_4_from_pure_match() {
        let ops = vec![op(Match, 50)];
        let (new_ops, shift) = trim_cigar_5prime(&ops, 4);
        assert_eq!(new_ops, vec![op(Match, 46)]);
        assert_eq!(shift, 4);
    }

    /// [S3, M50] trim 4 → [M49], shift=1.
    ///
    /// - S3 consumes 3 query bases (not reference): remaining 4→1, shift stays 0.
    /// - M50 is query+ref-consuming: trim 1 base → emit M49, shift += 1 → shift=1.
    #[test]
    fn trim_4_consumes_softclip_first() {
        let ops = vec![op(SoftClip, 3), op(Match, 50)];
        let (new_ops, shift) = trim_cigar_5prime(&ops, 4);
        assert_eq!(new_ops, vec![op(Match, 49)]);
        assert_eq!(shift, 1);
    }

    /// [M2, I3, M50] trim 4 → [I1, M50], shift=2.
    ///
    /// Step-by-step:
    /// - M2: query+ref, len=2 ≤ remaining=4 → drop; remaining=2, shift=2.
    /// - I3: query-only (no ref), len=3 > remaining=2 → emit I(3-2)=I1; shift unchanged; remaining=0.
    /// - M50: emitted intact.
    #[test]
    fn trim_4_passes_through_insertion() {
        let ops = vec![op(Match, 2), op(Insertion, 3), op(Match, 50)];
        let (new_ops, shift) = trim_cigar_5prime(&ops, 4);
        assert_eq!(new_ops, vec![op(Insertion, 1), op(Match, 50)]);
        assert_eq!(shift, 2);
    }

    // ── Task 12.2 test ────────────────────────────────────────────────────────

    /// Calling rewrite_record_inplace on a record with TLEN ≤ 9 returns Ok(false).
    #[test]
    fn rewrite_drops_small_tlen() {
        use noodles_sam::alignment::{
            record::cigar::{op::Kind, Op},
            RecordBuf,
        };
        let cigar: noodles_sam::alignment::record_buf::Cigar =
            [Op::new(Kind::Match, 50)].into_iter().collect();
        let seq = noodles_sam::alignment::record_buf::Sequence::from(b"ACGT".as_slice());
        let mut rec = RecordBuf::builder()
            .set_cigar(cigar)
            .set_sequence(seq)
            .set_template_length(5) // ≤ 9 → should be dropped
            .build();
        let result = rewrite_record_inplace(&mut rec, true).unwrap();
        assert!(!result, "record with TLEN=5 (≤9) must be dropped");
    }
}
