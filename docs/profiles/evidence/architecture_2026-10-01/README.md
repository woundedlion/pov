# Architecture device-capture evidence (2026-10-01)

Baseline source: `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`.
Checkpoint 1 (MindSplatter and HyperLattice): `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b`.
Checkpoint 2 (LatticeMelt, ChromaticLichen and KaleidoscopeSmooth):
`d25dd85dee17c28d2728f813e1ae00ae64690e1b`. This contains the ranked-instance
change and a WASM-only descriptor change; source identity is also recorded
in every capture provenance and summary. These initial observations remain
preserved below even where the current ranking uses a later capture.

The shipping KaleidoscopeSmooth ranking uses the repeated paired
`baseline-board4-kaleidoscopesmooth-profile` and
`aligned-board4-kaleidoscopesmooth-profile` captures. The latter is base
`0156d0d7490355ea99cab406f4704a128a55ccbf` plus its retained sole alignment
[source patch](aligned-board4-kaleidoscopesmooth-profile_source.diff); the
base also includes WASM-only diagnostics. The final source contains that
alignment at `76ea80495cef34aa347431bb0aaabb5aa9d3ee55`. The O3 report
records its own later capture source and sidecar.

The full-cycle `checkpoint1-kaleidoscopesmooth-profile` capture is additional
attribution evidence for source `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b`;
it does not replace the accepted shipping comparison.

Raw serial captures are copied byte-for-byte to `.txt` with original mtimes.
Each capture has a provenance sidecar, build and environment dumps for both
the measured profile image and its shipping Phantasm attestation image,
source status/diff when available, validation output, and a derived summary.
[Manifest](architecture_capture_manifest.json) records file sizes and SHA-256 hashes. No ELF binaries
are stored here; their hashes and original artifact paths remain in provenance.

Untouched captures passed the standard parser validation before filtering.
Exact runtime cadence excludes only setup frame 1; all later raw frame rows,
including the final partial counter window, transitions and spills, remain. Scope, ISR and preset-window summaries use complete windows
after the startup-containing window. Summaries record these distinct ranges.
Single-owner windows can contain a parameter morph, so they are not claimed
to measure isolated pure holds. Baseline/candidate finite-sample comparisons
are descriptive and do not establish a generalized speedup.
[Paired observations](architecture_comparisons.json) collect exact means,
peaks and spill counts; complete per-capture summaries carry preset attribution.

Rejected clock-validation attempts are excluded from reported captures and
rankings. The successful repeated capture retains the same 5 ppm gate.
