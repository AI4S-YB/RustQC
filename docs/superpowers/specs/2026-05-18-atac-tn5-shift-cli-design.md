# RustQC ATAC Tn5 Shift CLI Semantics

- **Date**: 2026-05-18
- **Status**: approved for implementation planning
- **Scope**: clarify how `rustqc atac` decides whether to apply Tn5 +4/-5 shift before
  TSS-dependent QC metrics.

## Problem

ATAC-seq QC mixes two metric classes:

- Basic BAM/fragment metrics such as bamQC, fragment size, NRF/PBC, and library complexity.
  These do not require insertion-site coordinates.
- TSS-dependent metrics such as TSSEscore, NFRscore, and PTscore. These use Tn5 insertion
  sites and should be computed from +4/-5 shifted 5' coordinates when the input BAM is not
  already shifted.

The current implementation always applies an in-memory +4/-5 shift for TSS-dependent metrics.
That matches ATACseqQC for unshifted input, but it can double-shift input BAMs that have already
been shifted upstream. There is also no supported way to request only the basic QC metrics when
the user does not want RustQC to shift unshifted input.

## User-Facing Interface

Add two ATAC-specific CLI/config controls:

```text
--tn5-shift <yes|no>
--input-is-shifted
```

`--tn5-shift` means "should RustQC apply Tn5 shift during this run?". It defaults to `yes`.

`--input-is-shifted` means "the input BAM coordinates are already Tn5 shifted". It defaults to
`false`.

The equivalent YAML fields live under `atac:`:

```yaml
atac:
  tn5_shift: yes
  input_is_shifted: false
```

CLI flags take precedence over YAML values.

## Behavior Matrix

| `tn5_shift` | `input_is_shifted` | Behavior |
|---|---:|---|
| `yes` | `false` | Apply in-memory Tn5 shift for TSS-dependent metrics; emit all QC outputs. |
| `yes` | `true` | Error out because the request would double-shift already shifted input. |
| `no` | `true` | Do not shift; treat input coordinates as insertion-site coordinates; emit all QC outputs. |
| `no` | `false` | Do not shift; skip TSS-dependent metrics and emit only basic QC outputs. |

TSS-dependent metrics are:

- TSSEscore
- NFRscore
- PTscore

Basic QC outputs remain available when TSS-dependent metrics are disabled:

- bamQC TSV and MAPQ TSV
- fragment-size TSV and SVG
- library-complexity TSV and SVG
- JSON summary with clear metadata that TSS-dependent metrics were skipped

## Output Metadata

The ATAC JSON summary should include the effective shift state so downstream tools can interpret
missing or present TSS-dependent metrics:

```json
"tn5_shift": {
  "requested": false,
  "input_is_shifted": false,
  "applied": false,
  "tss_dependent_metrics_enabled": false
}
```

When TSS-dependent metrics are disabled, the `tsse`, `nfr`, and `pt` summary sections should be
present with JSON `null` values. Keeping the top-level keys stable is less disruptive for
downstream consumers than omitting the sections entirely, while the `tn5_shift` metadata explains
why the metric payloads are unavailable.

## Error Handling

`--tn5-shift yes --input-is-shifted` is a hard error. This combination expresses both "please
shift" and "the input has already been shifted", and silently choosing one interpretation risks
quietly producing double-shifted metrics.

The error should mention both flags and explain how to proceed:

- Use `--tn5-shift no --input-is-shifted` for already shifted input.
- Omit `--input-is-shifted` for ordinary unshifted ATAC BAMs.

## Testing

Add focused tests for:

- CLI parsing of `--tn5-shift yes|no`.
- YAML parsing of `atac.tn5_shift` and `atac.input_is_shifted`.
- Resolution precedence from CLI over YAML.
- The four behavior-matrix cases.
- The hard error for `--tn5-shift yes --input-is-shifted`.
- A run with `--tn5-shift no` and unshifted input writes basic QC outputs and skips TSSE/NFR/PT.
- Default behavior remains equivalent to the current ATACseqQC-compatible path.

## Non-Goals

- Do not try to auto-detect whether an input BAM has already been Tn5 shifted. BAM headers and
  coordinate distributions are not reliable enough for that decision.
- Do not make `--emit-shifted-bam` control metric semantics. It remains an output-emission flag.
- Do not shift basic BAM/fragment metrics, because their definitions do not depend on insertion
  site coordinates.
