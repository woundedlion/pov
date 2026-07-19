# Phantasm ITCM spend ledger: `aeba37b5` → `d4816de0` (2026-07-16 … 07-19)

Complete per-commit accounting of `.text.itcm` on the phantasm image across the
window that opened when the arena cut (`aeba37b5`) freed a 32 KiB FlexRAM bank
for ITCM `-O3` promotions.

## Method

Every commit in `aeba37b5^..d4816de0` (140 builds) was checked out in an isolated
worktree and built with `pio run -e phantasm`; `.text.itcm` was read straight
from the linked ELF via `arm-none-eabi-size -A firmware.elf`, **independent of
the teensy gate** (the ELF links before the gate runs, so gate failures do not
perturb the measurement). Perf benefits are quoted verbatim from each commit
body and cross-checked against `docs/profiles/`. Device cadence figures for
MeshFeedback are from the on-device profile capture (`mf_ship_final.log`).

## Headline

| | ITCM | headroom to 196,608 B ceiling |
|---|---|---|
| `aeba37b5` (baseline, post arena cut) | 149,168 B | 47,432 B |
| `d4816de0` (HEAD) | 194,528 B | **2,072 B** |
| **net spend** | **+45,360 B** | |

Gross spend +53,952, gross reclaim −8,592. The "47k headroom before selective-O3"
framing is a conflation: selective-O3 landed *before* this window (07-15); the
~47k was **created** by the arena cut on 07-16 and then deliberately spent down
to fund ITCM promotions — exactly its purpose. There is no next bank available:
DTCM needs all 10 banks for 312,736 B of variables plus the 12,288 B stack floor.

## Spends that bought a cadence tier (74% of the total)

| Δ ITCM | commit | measured benefit | effect |
|---|---|---|---|
| **+19,568** | `d9bd43da` always_inline hot leaf helpers | peak 71.3→**56.8 ms**, spills 25→**0**/2,592, 16 fps locked, below the −O3 ceiling (58.8). **Opened the 6th ITCM bank** — the enabler for `405197d9` | HopfFibration + every HS_O3-region effect |
| **+4,592** | `3c61a33c` BZ raster under selective-O3 | 103→**58.4 ms**, 8→**16 fps**, 0 spills, within 3% of −O3 ceiling; physics 9.2→5.1 ms | BZReactionDiffusion |
| **+4,272** | `dedd3dde` interleave two pixels through the composite | **−6.3 ms peak** (largest single win of the effort), bit-identical; MeshFeedback 1🔴3🟡4🟢 → **0🔴1🟡7🟢**, spill 10.1%→0.21% | MeshFeedback |
| **+3,008** | `405197d9` MeshFeedback flush (R5) | flush 88.6→**45.3 ms** (=−O3 ceiling 45.2); 16 fps coverage ~0%→**54%** | MeshFeedback |
| **+2,080 / +1,648** | `0a5843a3` / `5e1a9fa4` Hopf trail gate + staging | enabling pair for `d9bd43da`'s lock; 16.8k orient+stage calls/frame out of −Os; pixel-identical | HopfFibration |
| **+1,248** | `3b43b353` concave sector walk (K1) | IslamicStars 120.2→**92.2 ms** (−23%), spills 755→363, two solids 8→16 fps; beats −O3 twin on worst shape | IslamicStars |
| **+1,216** | `93a70142` canvas pointer hoist | **−2.9 ms** on MeshFeedback Smoke (quantified in `7a16deaf`) | MeshFeedback |

The `dedd3dde` row corrects an earlier read of this ledger: its commit body quotes
only "output is bit-identical" and no timing, so a message-only audit scored it as
an unquantified spend. The device profile shows it is the **highest-yield 4 KB in
the window**. Lesson: a MeshFeedback restructuring can hide a 6 ms device win
behind a body that mentions only correctness — a commit-message ledger
under-credits algorithmic commits and needs device profiles to be complete.

## Reclaims (algorithmic, all correctness-neutral)

| Δ ITCM | commit | note |
|---|---|---|
| −1,184 | `3ece4f8a` fixed icosahedron | cadence regression accepted for looks (54%→37%); ITCM fell because the carousel machinery went |
| −816 | `6a3557d1` bail trail column cull on unboundable edges | correctness fix, net −816 |
| −464 | `9fd2a4e8` branchless bound in linear_rgb_in_gamut | 48 vmrs→16, 395→327 insns |
| −464 | `d4816de0` portable Fisher-Yates | determinism fix; −464 ride-along |
| −96 | `c93a1dbe` drop unread trail arc-length work | — |
| −48 | `2475eee6` Voronoi block-union shade | while delivering 4.2× host frame |

## Spends that bought nothing recorded

| Δ ITCM | commit | status |
|---|---|---|
| +1,408 | `c2b9ecf8` hoist clip test out of composite loop | no standalone timing, but structurally the enabler for `dedd3dde`'s paired path (hoists the per-cell clip test the interleave depends on) — likely bundled into that 6.3 ms, not separable |

The originally-suspected `dedd3dde` (+4,272) is **not** in this bucket; see above.

## Reverted experiments (measured harmful/dead)

| commit(s) | outcome | ITCM |
|---|---|---|
| `af49fdb6`+`21744baa`+`823fcb0a`, reverted `f7524d19` | per-row spans: probes −10-30% on host, but device span construction cost 31.7 ms/frame vs 7.2 ms saved → worst shape 89.4→108.4 ms, turned a green shape red. **Sized on host probe counts that could not see device construction cost.** | added +3,344, revert −2,928 |
| `54f613ef`, reverted `9761a432` | cubic-form gamut bisection: **measured dead on device** 1,953→1,962 cyc/px (latency-bound on a 16-link serial chain, not flop-bound; 1.13-1.32× on x86, 1.00× on device) | +320, revert −1,312 |

Both experiments died on device after being sized on host. This is the same trap
recorded for host `-O2` CSE in the raymarch-perf memory; it has now cost two
landings in four days. **Standing rule: no device-path lever lands on host timing
alone.**

### Is the per-row-span residual recoverable? No.

The four span-episode deltas sum to +416 B (added 3,344 − removed 2,928), which
reads like leftover residual. Direct inspection at HEAD refutes that:
- The experiment's named machinery (`emit_row_spans`, `raster_rowspan`) is **fully
  gone** from `scan.h`/`sdf.h`. (The `*_row_span` symbols in `plot.h` are the
  unrelated trail-clip edge-span feature; the `span` tokens in `scan.h` are the
  ordinary per-row longitude-interval rasterizer that always existed.)
- `rasterize_face` at HEAD is the bounding-rectangle scan with no row-narrowing.
- **~10 commits rewrote the touched functions after the revert** (`6cab193d`
  inlining, `d59384d6`/`ab8bad4e` branchless walks, `a9a8d500` modulo drop), so
  whatever layout/scheduling drift the episode left has been subsumed and
  re-optimized.

The +416 B is a net-accounting figure — the revert removed 416 B fewer than the
experiment added, measured against neighbors that shifted (`__LINE__`-sensitive
always-on checks, register allocation in shared templates). It is **not a block
of deletable code.** Recovering it would mean re-optimizing already-correct,
already-reverted code for sub-noise bytes; the real headroom lever is `variables`,
not this.

## Correction: `38b76187` (R7) reports its ITCM cost with the wrong sign

The commit body of `38b76187` and the doc commit `9b299107` both claim the
Raymarch march region "returns 208 B" / `RAM1 code 185,784 → 185,576 (−208)`.
**Measured: +752 B** (184,464 → 185,216) — a 960 B discrepancy, wrong direction.
The messages are landed history and cannot be rewritten; the accurate figure is
recorded here. (The only `−208`/`+208` token in tracked docs,
`probe_path_open_items.md:256`, is an unrelated rejected-inline measurement.)

## Everything else

The remaining ~90 commits are docs/tests/tools/CI (0 B delta) or small
algorithmic/correctness/quality commits in the ±32…+512 B noise band (Conway
morph machinery, hankin strap crossfades, RNG-seed plumbing, gamut/OKLab
micro-opts). None individually moves the headroom needle. Full per-commit ITCM
values: see the sweep TSV referenced in the audit session.

## State at HEAD

- ITCM 194,528 B / 196,608 B ceiling — **2,072 B free**, no next bank reachable.
- `tools/teensy_budgets.json` header comment is **stale**: it claims "476,544 B
  of RAM1 with 47,744 B free for locals," true only between `aeba37b5` and
  `d9bd43da` on 07-16. HEAD build reports `variables:312,736 code:194,536
  padding:2,072  free for local variables:14,944`. The bank `aeba37b5` freed was
  handed to ITCM by `d9bd43da` and never returned.
- Any further promotion must be paid for by an equal trim or by shrinking
  `variables` (DTCM), not by finding a spare bank.
