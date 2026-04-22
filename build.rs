//! Build script — compiles the C++ RNG shim for preseq compatibility and
//! embeds git/build metadata as compile-time environment variables.
//!
//! The shim wraps `std::mt19937` and `std::binomial_distribution` so that
//! RustQC's preseq bootstrap uses the exact same random number generation
//! as upstream preseq compiled on the same platform.
//!
//! Embedded variables:
//! - `GIT_SHORT_HASH` — short commit hash (e.g. `84ec57f`), or `unknown`
//! - `BUILD_TIMESTAMP` — UTC timestamp of the build (e.g. `2026-03-07T12:34:56Z`)

use std::process::Command;

fn main() {
    // --- C++ RNG shim ---
    cc::Build::new()
        .cpp(true)
        .file("cpp/rng_shim.cpp")
        .std("c++17")
        .warnings(true)
        .compile("rng_shim");

    println!("cargo:rerun-if-changed=cpp/rng_shim.cpp");

    // --- Git short hash ---
    // First check the GIT_SHORT_HASH env var (set via Docker build arg),
    // then fall back to running `git rev-parse` (works in local/CI builds).
    let git_hash = std::env::var("GIT_SHORT_HASH")
        .ok()
        .filter(|s| !s.is_empty() && s != "unknown")
        .unwrap_or_else(|| {
            Command::new("git")
                .args(["rev-parse", "--short", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        });
    println!("cargo:rustc-env=GIT_SHORT_HASH={}", git_hash);

    // --- Build timestamp (UTC, ISO-8601) ---
    // Computed from SystemTime directly so the build works on Windows where
    // POSIX `date` is unavailable.
    let build_ts = match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(d) => format_iso8601_utc(d.as_secs()),
        Err(_) => "unknown".to_string(),
    };
    println!("cargo:rustc-env=BUILD_TIMESTAMP={}", build_ts);

    // Rebuild when HEAD changes (new commits)
    println!("cargo:rerun-if-changed=.git/HEAD");
}

/// Format Unix epoch seconds as `YYYY-MM-DDTHH:MM:SSZ` without pulling in a
/// date crate at build time. Uses Howard Hinnant's civil-calendar algorithm.
fn format_iso8601_utc(secs: u64) -> String {
    let days = secs / 86400;
    let tod = secs % 86400;
    let (y, mo, d) = days_to_ymd(days);
    let h = tod / 3600;
    let m = (tod % 3600) / 60;
    let s = tod % 60;
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

fn days_to_ymd(days: u64) -> (u64, u64, u64) {
    let z = days + 719468;
    let era = z / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };
    (y, mo, d)
}
