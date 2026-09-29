# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `build/prof/hyperlattice_ship.log`, captured 2026-09-28 18:43 on COM4.
Replaces `profile_hyperlattice_teensy_2026-09-27.md`.

`97eb0bf78` has the same source tree as master `7678630fb`, which reverts an intervening HyperLattice flight change.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 170 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh HyperLattice profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:60864, data:149996, headers:8276` / `RAM1: variables:314944, code:13368, padding:19400, free:176576` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2449–2464 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.3 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `hl_shader_draw` averages 43.61 ms/f; its worst window is 49.96 ms/f (frames 2449–2464). Peak frame render is **56.12 ms** (frame 2482), and **0/2687** frames spilled. Setup frame 1 is excluded from both; it rendered 75.79 ms.

The previous shipping report (2026-09-27 01:02) recorded peak 🟢 55.509 (2) and spilled 🟢 0/1576 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 3 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2481–2496)

```
frame                     62.38 ms   37.43 Mcyc   100%
  pov_preserve_half       137.5 us    82.5 kcyc     0%
  hl_shader_draw          47.68 ms   28.61 Mcyc    76%
  hl_timeline_step          5.4 us     3.3 kcyc     0%
  canvas_clear             85.3 us    51.2 kcyc     0%
  canvas_buffer_wait      12.69 ms    7.62 Mcyc    20%
```

Wall min/avg/max = 52.08/62.38/70.70 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `hl_shader_draw` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `hl_shader_draw` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 2 | — | 🟢 56.12 | 0/1118 | 49.96 | 70/70 |
| 3 | — | 🟢 54.97 | 0/692 | 49.90 | 43/43 |
| 1 | — | 🟢 47.35 | 0/877 | 42.89 | 55/55 |

### Per-pixel figures

This capture carries no `filter_blend` counter, so no per-pixel blend figure is available.

## Column-ISR / DMA marshaling cost

```
isr_wake         1151/frame  min/avg/max 0.6/1.6/11.4 us  cpu 3.03%
isr_pack          144/frame  min/avg/max 6.2/6.7/9.5 us  cpu 1.55%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `hl_shader_draw` — 76% of the peak window, 47.68 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.09 ms/f.
4. `hl_timeline_step` — 0% of the peak window, 0.01 ms/f.

README cells: peak 🟢 56.12 (3), spilled 🟢 0/2687 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: the per-pixel pullback kernels carry `hot` or bare `HS_O3_FN` placement since `97eb0bf78`, in flash via `.text.hot`.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh HyperLattice profile 170 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.

## Supplemental fixed Triangular experiment

The opt-in `experimental-triangular-flight` preset was measured separately on
clean source `167fb6e02`, with a 70-second fixed-preset capture on COM3.
Runtime frames 2–442 render in **126.099 ms mean**,
**134.620 ms peak**, with **441/441 spills (100%)**.
Startup frame 1 (229.849 ms) is excluded. These results are outside
the normal two-preset roster and do not replace its measurements above.

The [paired experimental report](../hyperlattice_triangular_2026-09-27.md)
contains both image sizes, scope trees, ISR costs, exact frame ranges, matched
comparison, and portable evidence.

## Supplemental Octet 3D single-owner correction

Point-in-time snapshot of the corrected single-owner strut renderer.
This fixed experimental preset is separate from the normal HyperLattice cycle.
Its [shipping](#supplemental-octet-3d-single-owner-correction) and
[global-O3](../O3/profile_hyperlattice_teensy_2026-09-27.md#supplemental-octet-3d-single-owner-correction) captures use the same source.
[Raw capture](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship.txt),
[provenance](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship.provenance),
[summary](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship_summary.json),
[validation](../evidence/hyperlattice_octet_single_owner_2026-09-27/ship_validate.txt).
Captured 2026-09-27 22:50 local time.

### Setup

| | |
| --- | --- |
| Hardware | Teensy 4.0, 600 MHz, COM3, flywheel and DMA ISRs live |
| Image | `profile`, selective -O3; shipping retains HS_O3 shader helpers and HS_HOT_FLASH_MEMBER placement |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | HyperLattice, 288×144, fixed `experimental-octet-flight`, clean source `3c0ad6dcaad361db1dea671933e4e931e91faba0` |
| Method | 70 seconds, 16-frame windows; runtime frames 2–549; setup frame 1 excluded; scope windows 17–544 |
| Reproduce | `HS_PROFILE_TREE=<checkout> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh HyperLattice profile 70 16 '-DHS_ENABLE_HYPERLATTICE_EXPERIMENTS=1 -DHS_PROFILE_PRESET=2'` |

```text
teensy_size:   FLASH: code:76184, data:148692, headers:8596   free for files:1798144
teensy_size:    RAM1: variables:315008, code:20376, padding:12392   free for local variables:176512
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 433–448, root 1,199,191,741
cycles / 600 MHz versus 1,998,653 μs wall sum: **0.049 ppm**.
Build logs and environment dumps are retained beside the raw capture.

### Frame cadence

Runtime render minimum/mean/peak: **84.503/88.700/96.869 ms**.
Spills: **548/548 (100.0%)**.
Startup frame 1 took 167.999 ms and is excluded.
Mean wall time: 124.847 ms, approximately 8.01 fps.
The display budget is 62.5 ms. Peak render exceeds it by 34.369 ms.
`canvas_buffer_wait` is alignment idle before the next display flip.

### Phase-by-phase readout

The preset is held while its camera and palette continue moving. No preset
cycle or transition coverage is claimed. Worst complete shader window:

#### Fixed Octet 3D (frames 529–544)

```text
frame                       124.868 ms 74.921 Mcyc 100.0%
  pov_preserve_half           0.138 ms  0.083 Mcyc   0.1%
  hl_shader_draw             94.689 ms 56.814 Mcyc  75.8%
  hl_timeline_step            0.010 ms  0.006 Mcyc   0.0%
  canvas_clear                0.089 ms  0.053 Mcyc   0.1%
  canvas_buffer_wait         28.196 ms 16.917 Mcyc  22.6%
```

Wall minimum/mean/maximum: 124.740/124.868/124.943 ms.
Mean render in this window is 96.672 ms. All listed leaf scopes
run once per frame; their frame costs also give milliseconds per call.
The shader includes traversal, coverage and color without a finer breakdown.

#### Per-pixel figures

The nominal quadrant is 144×72 = 10,368 pixels; the shader margin evaluates
146×73 = 10,658 samples. Across the complete runtime windows, the shader
averages 86.633 ms/frame or 4877.1 cycles/sample.
Direct writes have no `filter_blend` calls. Candidate/layer counts were not captured.

### Column-ISR / DMA marshaling cost

```text
isr_wake        2304.4/f 0.575/1.676/24.428 us 3.09%
isr_pack         288.0/f 6.235/6.997/10.885 us 1.61%
isr_dma_submit   288.0/f 0.616/0.944/1.383 us 0.22%
```

Times are per-call minimum/mean/maximum, followed by CPU share.
Pack averages 6.997 μs versus
0.944 μs for submit. The 600-byte LED
image and black strobe take approximately 400 μs asynchronously at 12 MHz.
ISR share totals 4.92%, leaving approximately 59.425 ms
of foreground time per interval. Render already includes interrupts; its
mean/peak need 1.42×/1.55× reduction to fit.

### Summary ranking

1. `hl_shader_draw`: 86.633 ms/frame, 69.4% of root cycles.
2. `canvas_buffer_wait`: 36.277 ms/frame of display synchronization.
3. Preserve, clear and unscoped preparation account for the remainder.

The previous shipping capture averaged 90.406 ms and peaked at
103.733 ms. Matching frame indices 2–549 gives mean render
90.406 ms before and 88.700 ms after: **1.9% less render time**,
or **1.02× speedup**. Camera motion advances per frame;
the equally long captures can reach different frame indices. The previous
capture is source `39225587612b`; the intervening preset/UI and group-storage
changes mean this is a historical comparison, not an isolated rebuild A/B.
Both use the same fixed Octet 3D settings, board, driver, and compiler.


### Caveats

- All scopes include ISR time. No per-pixel profiling overhead is added.
- Direct writes have no `filter_blend` parenting artifact.
- Shipping uses selective-O3 shader traversal and cached-flash placement;
  global-O3 changes compiler flags, not the placement annotations.
- These measurements cover one authored preset and camera-frame range.
  WASM/native correctness tests are not comparable device timings.
- Both images were built from clean source; no engine or instrumentation
  changes were made for these captures. Setup is excluded, with no warmup cut.

### Harness

`targets/Profile/Profile.ino` supplies the existing HS_PROFILE scopes.
Use the locked reproduce command above; `just profile HyperLattice` without
the experimental/fixed-preset flags measures the normal roster.
