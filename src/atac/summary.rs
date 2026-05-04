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
    /// Always `"fixed_intervals_v1"` per spec.
    pub split_method: &'static str,
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
    pub atacseqqc_replicates: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BamqcSection {
    pub total_qnames: u64,
    pub duplicate_rate: f64,
    pub mitochondria_rate: f64,
    pub proper_pair_rate: f64,
    pub unmapped_rate: f64,
    pub has_unmapped_mate_rate: f64,
    pub not_passing_qc_rate: f64,
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
    /// Relative path to the per-sample fragsize TSV (e.g. `"fragsize/<sample>.fragsize.tsv"`).
    pub tsv_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TsseSection {
    /// Peak TSSE score (max of the smoothed profile).
    pub score: f64,
    /// Number of windows in the profile (typically 20).
    pub n_windows: u32,
    /// Loess-smoothed normalised signal vector (length == n_windows, typically 20).
    pub values: Vec<f64>,
    /// Relative path to the per-sample TSSE TSV.
    pub tsv_path: String,
}

/// Shared section for NFR and PT scores.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoreSection {
    pub n_tss: u32,
    pub median_score: f64,
    /// Relative path to the per-sample score TSV.
    pub tsv_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibComplexitySection {
    /// Number of rows in the saturation curve (typically 14).
    pub n_rows: u32,
    /// Expected distinct fragments at 1× sequencing depth (relative_size == 1.0).
    /// `None` when the bootstrap produced NaN (e.g. very small fixtures).
    pub extrapolated_total: Option<f64>,
    /// Relative path to the per-sample lib_complexity TSV.
    pub tsv_path: String,
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
                atacseqqc_replicates: "1.36.0".to_string(),
            },
            split_method: "fixed_intervals_v1",
            bamqc: BamqcSection {
                total_qnames: 1000,
                duplicate_rate: 0.05,
                mitochondria_rate: 0.02,
                proper_pair_rate: 0.95,
                unmapped_rate: 0.01,
                has_unmapped_mate_rate: 0.02,
                not_passing_qc_rate: 0.005,
                nrf: 0.8,
                pbc1: 0.9,
                pbc2: 3.5,
                mapq_histogram: mapq_hist.clone(),
            },
            fragsize: FragsizeSection {
                total_pairs: 800,
                tsv_path: "fragsize/test_sample.fragsize.tsv".to_string(),
            },
            tsse: TsseSection {
                score: 7.5,
                n_windows: 20,
                values: vec![1.0; 20],
                tsv_path: "tsse/test_sample.tsse.tsv".to_string(),
            },
            nfr: ScoreSection {
                n_tss: 10,
                median_score: 0.5,
                tsv_path: "nfr/test_sample.nfr.tsv".to_string(),
            },
            pt: ScoreSection {
                n_tss: 10,
                median_score: 1.2,
                tsv_path: "pt/test_sample.pt.tsv".to_string(),
            },
            lib_complexity: LibComplexitySection {
                n_rows: 14,
                extrapolated_total: Some(750.0),
                tsv_path: "lib_complexity/test_sample.libcomplexity.tsv".to_string(),
            },
        };

        let json = serde_json::to_string(&s).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let j = &v;

        // Assert top-level keys per spec.
        for k in ["sample", "tool_versions", "split_method", "bamqc", "fragsize", "tsse", "nfr", "pt", "lib_complexity"] {
            assert!(j.get(k).is_some(), "missing top-level key: {}", k);
        }

        // Assert bamqc fields per spec.
        for k in ["total_qnames", "duplicate_rate", "mitochondria_rate", "proper_pair_rate",
                  "unmapped_rate", "has_unmapped_mate_rate", "not_passing_qc_rate",
                  "nrf", "pbc1", "pbc2", "mapq_histogram"] {
            assert!(j["bamqc"].get(k).is_some(), "bamqc missing {}", k);
        }

        // Assert tool_versions fields per spec.
        assert!(j["tool_versions"].get("rustqc").is_some(), "tool_versions missing rustqc");
        assert!(j["tool_versions"].get("atacseqqc_replicates").is_some(), "tool_versions missing atacseqqc_replicates");

        // Assert split_method value.
        assert_eq!(j["split_method"].as_str().unwrap(), "fixed_intervals_v1");

        // Assert tsse.values is a 20-element array.
        assert_eq!(j["tsse"]["values"].as_array().unwrap().len(), 20);

        // Assert mapq_histogram round-trips.
        let hist_rt = j["bamqc"]["mapq_histogram"].as_object().unwrap();
        assert_eq!(hist_rt["30"].as_u64().unwrap(), 500);
        assert_eq!(hist_rt["60"].as_u64().unwrap(), 200);

        // Assert tsv_path fields exist.
        assert!(j["fragsize"].get("tsv_path").is_some(), "fragsize missing tsv_path");
        assert!(j["tsse"].get("tsv_path").is_some(), "tsse missing tsv_path");
        assert!(j["nfr"].get("tsv_path").is_some(), "nfr missing tsv_path");
        assert!(j["pt"].get("tsv_path").is_some(), "pt missing tsv_path");
        assert!(j["lib_complexity"].get("tsv_path").is_some(), "lib_complexity missing tsv_path");

        // extrapolated_total present (may be null for NaN fixtures).
        assert!(j["lib_complexity"].get("extrapolated_total").is_some(), "lib_complexity missing extrapolated_total");
    }
}
