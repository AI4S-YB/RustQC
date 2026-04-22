//! BAM/SAM I/O facade over [`noodles`].
//!
//! This module is the single place the rest of the crate imports from when it
//! needs to read alignment records. It exists for two reasons:
//!
//! 1. **Portability.** HTSlib (via `rust-htslib`) depends on a large C/POSIX
//!    stack that does not build cleanly on Windows MSVC. `noodles` is pure
//!    Rust and builds identically on Linux, macOS, and Windows.
//!
//! 2. **Semantic parity.** `noodles` and `rust-htslib` differ subtly on a few
//!    points that the migration has to bridge rather than expose to every
//!    call site:
//!
//!    - **MAPQ `255`.** The SAM spec reserves `255` to mean "mapping quality
//!      unavailable". `rust-htslib` returns the raw `u8`; `noodles` maps
//!      `255` to `None`. [`mapq`] preserves the `rust-htslib` convention so
//!      existing MAPQ-threshold logic keeps behaving the same.
//!    - **Alignment position.** `rust-htslib::Record::pos()` is `0`-based and
//!      returns `-1` for unmapped reads. `noodles` uses a `1`-based
//!      [`Position`] wrapped in `Option<io::Result<…>>`. [`pos_0based`]
//!      yields the `rust-htslib` shape.
//!    - **Aux integer tags.** `rust-htslib` surfaces six distinct `Aux`
//!      variants (`I8/U8/I16/U16/I32/U32`); `noodles` uses a unified `Value`
//!      enum with an [`as_int`]-style convenience. [`get_aux_int`] keeps the
//!      caller ergonomics identical to the old `get_aux_int` helper.
//!
//! Most call sites want to treat records exactly like the previous
//! `rust_htslib::bam::Record`, so the module re-exports the noodles types
//! under short aliases and adds the bridging helpers above.
//!
//! [`Position`]: noodles_core::Position
//! [`as_int`]: noodles_sam::alignment::record::data::field::Value::as_int

use std::fs::File;
use std::io::{self, Cursor, Read};
use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

// --- Public re-exports ---------------------------------------------------
//
// Callers import these instead of reaching into `noodles_*` directly. This
// keeps the migration diff small and lets us change the backing crate later
// without touching every file.

pub use noodles_bam::Record;
pub use noodles_sam::alignment::record::cigar::op::Kind as CigarKind;
pub use noodles_sam::Header;

/// Boxed inner reader so a single `Reader` type can wrap either a `File`
/// (for BAM input) or an in-memory `Cursor` (for SAM input, which is
/// transcoded to BAM bytes once at open time).
type BoxedRead = Box<dyn Read + Send>;

/// A sequential BAM reader. Backed either by a BGZF-decoded file (for BAM
/// input) or by a BGZF-encoded in-memory cursor (for SAM input that has been
/// transcoded at open time).
pub type Reader = bam::io::Reader<bgzf::io::Reader<BoxedRead>>;

/// An indexed BAM reader backed by a BGZF-decoded file with a BAI/CSI index
/// resolved automatically from the standard sibling path.
pub type IndexedReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

// --- Reader constructors --------------------------------------------------

/// Open an alignment file sequentially and read its header.
///
/// Accepts both BAM (binary, BGZF-compressed) and SAM (plain text). SAM input
/// is transcoded to BAM bytes once at open time so downstream code can treat
/// every record as a [`noodles_bam::Record`].
///
/// Returns both the reader and the owned [`Header`] so downstream code can
/// call `reader.query(&header, …)` on an indexed variant without juggling
/// borrows. The header is required for most `noodles` operations.
pub fn open(path: &Path) -> io::Result<(Reader, Header)> {
    let inner: BoxedRead = if is_sam_path(path)? {
        Box::new(transcode_sam_to_bam(path)?)
    } else {
        Box::new(File::open(path)?)
    };
    let mut reader = bam::io::Reader::new(inner);
    let header = reader.read_header()?;
    Ok((reader, header))
}

/// Detect whether `path` looks like a SAM (plain-text) file rather than BAM.
///
/// Uses the file extension first (`.sam`/`.SAM`) and falls back to sniffing
/// the first two bytes for the BGZF magic (`0x1f 0x8b`).
fn is_sam_path(path: &Path) -> io::Result<bool> {
    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
        let lower = ext.to_ascii_lowercase();
        if lower == "sam" {
            return Ok(true);
        }
        if lower == "bam" || lower == "bgz" || lower == "gz" {
            return Ok(false);
        }
    }
    // Fallback: sniff the magic bytes.
    let mut f = File::open(path)?;
    let mut magic = [0u8; 2];
    use std::io::Read as _;
    let n = f.read(&mut magic)?;
    Ok(n < 2 || magic != [0x1f, 0x8b])
}

/// Read a SAM file and re-emit it as a BGZF-compressed BAM byte stream that
/// can then be fed to `bam::io::Reader`. Returns a `Cursor` over the encoded
/// bytes so downstream `read_header()` and `records()` calls behave exactly
/// like a plain-BAM input.
fn transcode_sam_to_bam(path: &Path) -> io::Result<Cursor<Vec<u8>>> {
    use sam::alignment::io::Write as _;

    let file = File::open(path)?;
    let mut sam_reader = sam::io::Reader::new(io::BufReader::new(file));
    let header = sam_reader.read_header()?;

    let mut out = Vec::new();
    {
        let mut bam_writer = bam::io::Writer::new(&mut out);
        bam_writer.write_alignment_header(&header)?;
        for result in sam_reader.record_bufs(&header) {
            let record = result?;
            bam_writer.write_alignment_record(&header, &record)?;
        }
        bam_writer.finish(&header)?;
    }

    Ok(Cursor::new(out))
}

/// Open a BAM file with its BAI/CSI index for region queries.
///
/// The builder automatically resolves `<path>.bai` and falls back to
/// `<path>.csi`, matching the behaviour of `rust_htslib::bam::IndexedReader`.
pub fn open_indexed(path: &Path) -> io::Result<(IndexedReader, Header)> {
    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(path)?;
    let header = reader.read_header()?;
    Ok((reader, header))
}

// --- Semantic bridges over rust-htslib conventions ------------------------

/// Return MAPQ as a raw `u8`, matching the `rust-htslib` convention.
///
/// `noodles` treats `255` as a "missing" sentinel and yields `None`; this
/// helper flattens that back to `255` so filters such as `mapq > cut` behave
/// the way they did before the migration.
#[inline]
pub fn mapq(record: &Record) -> u8 {
    record.mapping_quality().map(u8::from).unwrap_or(255)
}

/// Return the 0-based alignment start position.
///
/// Mirrors `rust_htslib::bam::Record::pos()`: `-1` for unmapped reads,
/// otherwise the POS field minus one. Panics on the (in practice unreachable)
/// case where `noodles` fails to decode an already-parsed alignment start —
/// this matches the infallible shape of the old htslib helper and keeps
/// every call site free of `?` or result plumbing.
#[inline]
pub fn pos_0based(record: &Record) -> i64 {
    match record.alignment_start() {
        Some(result) => usize::from(result.expect("noodles: decode alignment_start")) as i64 - 1,
        None => -1,
    }
}

/// Return the 0-based mate alignment start, matching `Record::mpos()`.
#[inline]
pub fn mpos_0based(record: &Record) -> i64 {
    match record.mate_alignment_start() {
        Some(result) => {
            usize::from(result.expect("noodles: decode mate_alignment_start")) as i64 - 1
        }
        None => -1,
    }
}

/// Return the read name (QNAME) as a byte slice, with `*` for unset names.
///
/// `rust_htslib` returns the raw stored name (which is `b"*"` in BAM when
/// unset); `noodles` returns `None` in that case. This helper keeps the
/// always-a-slice shape so checksum computations that hash QNAME stay
/// bit-identical to the pre-migration behaviour.
#[inline]
pub fn qname(record: &Record) -> &[u8] {
    record.name().map(|n| n.as_ref()).unwrap_or(b"*")
}

/// Return the 0-based, exclusive alignment end position.
///
/// Matches `rust_htslib::bam::record::CigarStringView::end_pos()`: for a read
/// starting at 0-based `pos` with reference-consuming CIGAR length `n`, this
/// returns `pos + n` (i.e. one past the last aligned reference base). Yields
/// `-1` for unmapped reads, mirroring what the old code observed when the
/// CIGAR was empty and `pos()` was `-1`.
pub fn end_pos_0based(record: &Record) -> i64 {
    use noodles_sam::alignment::Record as _;
    match record.alignment_end() {
        Some(result) => usize::from(result.expect("noodles: decode alignment_end")) as i64,
        None => -1,
    }
}

/// Return the reference sequence id (`tid`) as a signed integer.
///
/// Mirrors `rust_htslib::bam::Record::tid()`: unmapped reads yield `-1`.
#[inline]
pub fn tid(record: &Record) -> i32 {
    match record.reference_sequence_id() {
        Some(result) => result.expect("noodles: decode reference_sequence_id") as i32,
        None => -1,
    }
}

/// Return the mate's reference sequence id.
#[inline]
pub fn mtid(record: &Record) -> i32 {
    match record.mate_reference_sequence_id() {
        Some(result) => result.expect("noodles: decode mate_reference_sequence_id") as i32,
        None => -1,
    }
}

/// Extract an integer auxiliary tag value as `i64`.
///
/// Replacement for the previous `rust_htslib`-specific `get_aux_int`: accepts
/// any of the six integer widths (signed/unsigned `8/16/32`) and yields
/// `None` for absent tags or non-integer types. Silently ignores malformed
/// aux fields (consistent with the old behaviour, which also used `.ok()`).
pub fn get_aux_int(record: &Record, tag: &[u8; 2]) -> Option<i64> {
    let value = record.data().get(tag)?.ok()?;
    value.as_int()
}

/// Return the 4-bit encoded IUPAC base code at position `i` in the record's
/// sequence, matching `rust_htslib::bam::record::Seq::encoded_base`.
///
/// Byte layout in BAM: two bases per byte, high nibble first, low nibble
/// second. This preserves the exact hash stream that `hash_sequence_encoded`
/// fed into `DefaultHasher` before the migration.
///
/// Panics on out-of-range `i`; callers are expected to bound by
/// `record.sequence().len()`.
#[inline]
pub fn encoded_base(record: &Record, i: usize) -> u8 {
    let sequence = record.sequence();
    let bytes = sequence.as_ref();
    let b = bytes[i / 2];
    if i.is_multiple_of(2) {
        b >> 4
    } else {
        b & 0x0f
    }
}

/// Decode CIGAR operations into owned `(kind, len)` pairs.
///
/// `noodles` exposes CIGAR as a lazy iterator of `io::Result<Op>`; downstream
/// code historically expected a slice of variants it could match on directly.
/// This helper collects once per record, letting callers keep their old
/// per-op pattern matching with only the variant names changing.
pub fn cigar_ops(record: &Record) -> io::Result<Vec<(CigarKind, u32)>> {
    record
        .cigar()
        .iter()
        .map(|op| op.map(|op| (op.kind(), op.len() as u32)))
        .collect()
}

// --- Header helpers -------------------------------------------------------

/// Flatten `header.reference_sequences()` into a `(name, length)` vector,
/// matching the `(target_name, target_len)` shape that the old
/// `rust-htslib` code used for idxstats and similar reports.
pub fn reference_sequences(header: &Header) -> Vec<(String, u64)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            let name = String::from_utf8_lossy(name.as_ref()).into_owned();
            let length = usize::from(map.length()) as u64;
            (name, length)
        })
        .collect()
}

// --- Tests ----------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn open_test_bam() -> (Reader, Header) {
        open(Path::new("tests/data/test.bam")).expect("open test.bam")
    }

    #[test]
    fn mapq_missing_maps_to_255() {
        let (mut reader, _) = open_test_bam();
        let mut any_missing = false;
        for result in reader.records() {
            let record = result.unwrap();
            if record.mapping_quality().is_none() {
                assert_eq!(mapq(&record), 255);
                any_missing = true;
                break;
            }
        }
        // The test fixture has MAPQ=255 throughout; guard against future
        // fixture changes silently making the assertion vacuous.
        assert!(
            any_missing,
            "expected at least one MAPQ=255 read in fixture"
        );
    }

    #[test]
    fn pos_is_zero_based() {
        let (mut reader, _) = open_test_bam();
        let record = reader.records().next().unwrap().unwrap();
        let pos = pos_0based(&record);
        // First record in the fixture is mapped; exact coordinate is not
        // what we're asserting, only that it's non-negative and one less
        // than the 1-based noodles position.
        let raw = usize::from(record.alignment_start().unwrap().unwrap()) as i64;
        assert_eq!(pos, raw - 1);
        assert!(pos >= 0);
    }

    #[test]
    fn get_aux_int_reads_nh() {
        let (mut reader, _) = open_test_bam();
        let record = reader.records().next().unwrap().unwrap();
        // NH is present on every read in the fixture.
        let nh = get_aux_int(&record, b"NH").expect("NH tag");
        assert!(nh >= 1);
    }

    #[test]
    fn reference_sequences_shape() {
        let (_, header) = open_test_bam();
        let refs = reference_sequences(&header);
        assert!(!refs.is_empty(), "fixture has at least one reference");
        for (name, length) in &refs {
            assert!(!name.is_empty());
            assert!(*length > 0);
        }
    }
}
