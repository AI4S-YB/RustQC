//! Tn5 +4/-5 coordinate-only shift for ATAC-seq reads.
//!
//! This module implements the coordinate-only path used by the metric pipeline
//! (Phase 13). Full record rewrite (CIGAR/SEQ/QUAL) is Phase 12.

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShiftedFrag {
    pub pos5p: u64,
    pub tlen: i64,
}

#[allow(dead_code)]
pub fn shift_5prime(pos5p: u64, is_plus: bool, tlen: i64) -> Option<ShiftedFrag> {
    let new_pos = if is_plus {
        pos5p + 4
    } else {
        pos5p.checked_sub(5)?
    };
    let new_tlen = if tlen == 0 {
        0
    } else {
        let sign = tlen.signum();
        let abs = tlen.unsigned_abs() as i64;
        if abs <= 9 {
            return None;
        }
        sign * (abs - 9)
    };
    Some(ShiftedFrag {
        pos5p: new_pos,
        tlen: new_tlen,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_strand_shifts_pos_by_plus4_and_shrinks_tlen_by_9() {
        let result = shift_5prime(100, true, 200).unwrap();
        assert_eq!(result.pos5p, 104);
        assert_eq!(result.tlen, 191);
    }

    #[test]
    fn minus_strand_shifts_pos_by_minus5_and_shrinks_tlen_by_9() {
        let result = shift_5prime(100, false, -200).unwrap();
        assert_eq!(result.pos5p, 95);
        assert_eq!(result.tlen, -191);
    }

    #[test]
    fn drops_fragment_when_tlen_le_9() {
        assert!(shift_5prime(100, true, 9).is_none());
        assert!(shift_5prime(100, false, -9).is_none());
    }

    #[test]
    fn drops_when_minus_strand_pos_underflows() {
        assert!(shift_5prime(3, false, -50).is_none());
    }
}
