# Convex-face AA candidate profiling evidence, 2026-09-28

Baseline `0c02f3912a98677184d3cf43b31a5d468d4b8ba3` is current source; candidate
`d4b8d7b75f3e2a87ccd2801a946d0b8ad57e6cc6` remains unlanded. The standard profiles
and roster rows use the baseline. The [comparison](../../face_aa_2026-09-28.md)
reports both revisions, leading with exact matched runtime peaks and spills.

Eight captures cover both effects and both build configs before/after.
The supported locked `profile_one.sh` wrapper ran shipping on COM3 and global
O3 on COM4, with each source pair using the same board. HankinSolids uses
260 seconds and its ordered tour; IslamicStars uses 210 seconds and TS4.
All use 16-frame windows and a 4000-revolution epoch.

- `before_*.txt` / `after_*.txt` are untouched raw captures, including setup
  and trailing per-frame rows. The eight `.provenance` records attest compiler,
  source, ELF and environment hashes.
- `candidate.patch` is the complete binary-capable Git diff from baseline to
  candidate; `candidate.json` records both full source SHAs and the patch's
  SHA-256. Apply it to the recorded baseline to reproduce the source change
  without requiring access to the unlanded local branch.
- `*_summary.json` contains exact frame/preset peaks, complete-window scope
  totals, markers, final geometry and raw-log SHA-256. `comparison.json`
  intersects frame numbers and verifies shape ownership before comparing.
- `*_envdump.json` stores exact UTF-8 dump content and its SHA-256; encoding
  `content` as UTF-8 reconstructs the original bytes. `*_build.txt` trims only
  trailing whitespace from readable build output.
- Complete local binary/map/attestation bundles are retained outside the
  worktrees in `C:/work/Holosphere/.git/implement-sessions/face21-captures/`,
  under `before-shipping`, `before-o3`, `after-shipping`, `after-o3`.
  Provenance preserves the original `build/prof/` relative path; its artifact
  basename identifies the same bundle under each archive's `artifacts/`.

Every untouched capture passes the stock parser validation. Additional cycle
checks require all 18 Hankin nodes and 30 landings, and all 23 Islamic recipes
plus the first recipe's return. Setup frame 1 is excluded from runtime
statistics; the whole startup window is excluded from scope/ISR statistics.
Initial Hankin frames belong to tetrahedron, not a fictitious extra preset.
Islamic V/E/F/I comes from `Built Shape`, not its recipe seed at spawn.

Scope trees retain raw duplicate-name rows. Shared or duplicate counters are
not treated as exclusive costs; unavailable pixel-counter rows are not zero.
The candidate measurements do not mark the code-review finding resolved.

Selected native visual comparisons, complete sampled metrics and their
[method](visuals/METHOD.md) are retained in `visuals/`; `manifest.json` hashes
each portable visual artifact. Raw RGB16 captures and diagnostic harness
source remain in the local archive named by the method. The candidate cost
was rejected; the rendering change and review finding remain unlanded/open.
