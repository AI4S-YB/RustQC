//! Mitochondrial chromosome detection from BAM @SQ names.

/// Names matched by the auto-detect logic (case-sensitive): `chrM`, `MT`, `Mito`.
pub fn detect_mito<'a>(seq_names: &'a [String]) -> Option<&'a str> {
    seq_names
        .iter()
        .find(|n| matches!(n.as_str(), "chrM" | "MT" | "Mito"))
        .map(String::as_str)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(strs: &[&str]) -> Vec<String> {
        strs.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn detects_chrm() {
        assert_eq!(detect_mito(&s(&["chr1", "chr2", "chrM"])), Some("chrM"));
    }

    #[test]
    fn detects_mt_for_ensembl() {
        assert_eq!(detect_mito(&s(&["1", "2", "MT"])), Some("MT"));
    }

    #[test]
    fn detects_mito_for_yeast() {
        assert_eq!(detect_mito(&s(&["I", "II", "Mito"])), Some("Mito"));
    }

    #[test]
    fn returns_none_when_absent() {
        assert_eq!(detect_mito(&s(&["chr1", "chr2"])), None);
    }

    #[test]
    fn does_not_match_substrings() {
        // chrMT or MTother should NOT match (we use exact equality on canonical names).
        assert_eq!(detect_mito(&s(&["chrMT"])), None);
        assert_eq!(detect_mito(&s(&["MTother"])), None);
    }
}
