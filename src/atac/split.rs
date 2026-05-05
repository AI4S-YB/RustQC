//! Fixed-interval fragment-size split: NFR / mono / di / tri buckets.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(dead_code)]
pub enum FragBucket {
    Nfr,
    Mono,
    Di,
    Tri,
    Other,
}

#[allow(dead_code)]
pub fn classify(abs_tlen: u32) -> FragBucket {
    if abs_tlen < 100 {
        FragBucket::Nfr
    } else if (180..=247).contains(&abs_tlen) {
        FragBucket::Mono
    } else if (315..=473).contains(&abs_tlen) {
        FragBucket::Di
    } else if (558..=615).contains(&abs_tlen) {
        FragBucket::Tri
    } else {
        FragBucket::Other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_cases_match_atacseqqc_intervals() {
        assert_eq!(classify(0), FragBucket::Nfr);
        assert_eq!(classify(99), FragBucket::Nfr);
        assert_eq!(classify(100), FragBucket::Other);
        assert_eq!(classify(180), FragBucket::Mono);
        assert_eq!(classify(247), FragBucket::Mono);
        assert_eq!(classify(248), FragBucket::Other);
        assert_eq!(classify(315), FragBucket::Di);
        assert_eq!(classify(473), FragBucket::Di);
        assert_eq!(classify(558), FragBucket::Tri);
        assert_eq!(classify(615), FragBucket::Tri);
        assert_eq!(classify(616), FragBucket::Other);
    }
}
