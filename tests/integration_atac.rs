//! Integration smoke test for `rustqc atac`.
//!
//! Runs the real binary against the GL1 fixture and verifies all expected
//! output files are produced. Numerical correctness is validated in Phase 14
//! against R golden values.

#[test]
fn rustqc_atac_runs_on_gl1_fixture() {
    let outdir = tempfile::tempdir().unwrap();
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .args([
            "atac",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL1.bam"),
            "--gtf",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/atac/GL_stub.gtf"),
            "--outdir",
            outdir.path().to_str().unwrap(),
            "--sample-name",
            "GL1",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "rustqc atac exited non-zero: {:?}", status);

    // Smoke check: every expected metric file exists.
    for sub in [
        "bamqc/GL1.bamqc.tsv",
        "bamqc/GL1.mapq.tsv",
        "fragsize/GL1.fragsize.tsv",
        "fragsize/GL1.fragsize.svg",
        "tsse/GL1.tsse.tsv",
        "tsse/GL1.tsse.svg",
        "nfr/GL1.nfr.tsv",
        "pt/GL1.pt.tsv",
        "lib_complexity/GL1.libcomplexity.tsv",
        "lib_complexity/GL1.libcomplexity.svg",
        "GL1.atac.summary.json",
    ] {
        assert!(
            outdir.path().join(sub).exists(),
            "missing output file: {}",
            sub
        );
    }

    // Verify JSON summary is valid JSON with expected top-level keys.
    let json_path = outdir.path().join("GL1.atac.summary.json");
    let json_str = std::fs::read_to_string(&json_path).unwrap();
    let v: serde_json::Value = serde_json::from_str(&json_str)
        .expect("GL1.atac.summary.json is not valid JSON");
    let obj = v.as_object().expect("JSON root is not an object");
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
        assert!(obj.contains_key(*key), "JSON missing key: {}", key);
    }
    // Sanity check: sample name matches.
    assert_eq!(obj["sample"].as_str().unwrap(), "GL1");
}
