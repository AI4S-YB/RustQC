//! JSON summary schema for `rustqc atac`.
//!
//! The `AtacSummary` struct is serialised to JSON by the driver.  The spec
//! schema lives in `docs/superpowers/specs/2026-05-04-atac-seq-qc-design.md`.

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Top-level summary
// ---------------------------------------------------------------------------

/// Complete JSON summary for one ATAC-seq sample.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtacSummary {
    pub schema_version: String,
    pub sample: String,
    pub tool_versions: ToolVersions,
    pub bamqc: BamqcSection,
    pub fragsize: FragsizeSection,
    pub tsse: TsseSection,
    pub nfr: ScoreSection,
    pub pt: ScoreSection,
    pub lib_complexity: LibComplexitySection,
}

// ---------------------------------------------------------------------------
// Subsections
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolVersions {
    pub rustqc: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BamqcSection {
    pub total_qnames: u64,
    pub duplicate_rate: f64,
    pub mitochondria_rate: f64,
    pub proper_pair_rate: f64,
    pub nrf: f64,
    pub pbc1: f64,
    pub pbc2: f64,
    /// MAPQ value (as string) → read count.
    pub mapq_histogram: serde_json::Map<String, serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragsizeSection {
    /// Total fragment pairs counted (sum of counts across all lengths 1..=1010).
    pub total_pairs: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TsseSection {
    /// Peak TSSE score (max of the smoothed profile).
    pub score: f64,
    /// Number of windows in the profile (typically 20).
    pub n_windows: u32,
}

/// Shared section for NFR and PT scores.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoreSection {
    pub n_tss: u32,
    pub median_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibComplexitySection {
    /// Number of rows in the saturation curve (typically 14).
    pub n_rows: u32,
    /// Expected distinct fragments at 1× sequencing depth (relative_size == 1.0).
    pub distinct_at_1x: f64,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_keys_match_spec() {
        // Build a synthetic AtacSummary and round-trip through JSON.
        let mut mapq_hist = serde_json::Map::new();
        mapq_hist.insert("30".to_string(), serde_json::Value::Number(500u64.into()));
        mapq_hist.insert("60".to_string(), serde_json::Value::Number(200u64.into()));

        let s = AtacSummary {
            schema_version: "1.0".to_string(),
            sample: "test_sample".to_string(),
            tool_versions: ToolVersions {
                rustqc: "0.3.0".to_string(),
            },
            bamqc: BamqcSection {
                total_qnames: 1000,
                duplicate_rate: 0.05,
                mitochondria_rate: 0.02,
                proper_pair_rate: 0.95,
                nrf: 0.8,
                pbc1: 0.9,
                pbc2: 3.5,
                mapq_histogram: mapq_hist.clone(),
            },
            fragsize: FragsizeSection { total_pairs: 800 },
            tsse: TsseSection { score: 7.5, n_windows: 20 },
            nfr: ScoreSection { n_tss: 10, median_score: 0.5 },
            pt: ScoreSection { n_tss: 10, median_score: 1.2 },
            lib_complexity: LibComplexitySection {
                n_rows: 14,
                distinct_at_1x: 750.0,
            },
        };

        let json = serde_json::to_string(&s).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();

        // Assert top-level keys.
        let obj = v.as_object().unwrap();
        for key in &[
            "schema_version",
            "sample",
            "tool_versions",
            "bamqc",
            "fragsize",
            "tsse",
            "nfr",
            "pt",
            "lib_complexity",
        ] {
            assert!(obj.contains_key(*key), "missing top-level key: {}", key);
        }

        // Assert mapq_histogram round-trips.
        let hist_rt = obj["bamqc"]["mapq_histogram"].as_object().unwrap();
        assert_eq!(hist_rt["30"].as_u64().unwrap(), 500);
        assert_eq!(hist_rt["60"].as_u64().unwrap(), 200);
    }
}
