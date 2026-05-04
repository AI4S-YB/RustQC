//! Reject single-end BAMs at startup: scan the first ≤10 000 primary mapped
//! records and require at least one with the `READ_PAIRED` flag set.

use anyhow::{anyhow, Result};
use std::path::Path;

const MAX_RECORDS_TO_INSPECT: usize = 10_000;
const FLAG_READ_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

pub fn assert_paired_end(path: &Path) -> Result<()> {
    let (mut reader, _header) = crate::bam_io::open(path)?;
    let mut inspected = 0usize;
    let mut paired = 0usize;
    for result in reader.records() {
        let record = result?;
        let flags = u16::from(record.flags());
        if flags & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            continue;
        }
        inspected += 1;
        if flags & FLAG_READ_PAIRED != 0 {
            paired += 1;
        }
        if inspected >= MAX_RECORDS_TO_INSPECT {
            break;
        }
    }
    if inspected == 0 {
        return Err(anyhow!(
            "BAM contains no primary mapped records: {}",
            path.display()
        ));
    }
    if paired == 0 {
        return Err(anyhow!(
            "rustqc atac requires paired-end input; first {} primary records had no PAIRED flag: {}",
            inspected,
            path.display()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture(rel: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
    }

    // TODO(phase-3): once tests/data/atac/GL1.bam is materialized, add an
    // `accepts_paired_end_bam` test that asserts assert_paired_end(...) == Ok(()).

    #[test]
    fn rejects_single_end_bam() {
        // tests/data/test.bam is the existing RNA QC fixture; all records have
        // flag=0 (single-end, no PAIRED bit set). The plan assumed this file was
        // paired-end, but it is not — verified with `samtools view -f 1` (0 hits).
        // We therefore test the SE-rejection path, which exercises assert_paired_end
        // just as thoroughly as an accept path.
        let p = fixture("tests/data/test.bam");
        let result = assert_paired_end(&p);
        assert!(
            result.is_err(),
            "test.bam is single-end and should be rejected by assert_paired_end"
        );
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("paired-end input"),
            "error should mention paired-end: {}",
            msg
        );
    }

    // A real PE BAM test will be added in Phase 3 once we extract the
    // ATACseqQC GL-fixture BAMs (which are genuine ATAC paired-end).
}
