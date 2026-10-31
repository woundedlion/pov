# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot. Replaces the October 11 selective-O3 capture.

**The requested no-spill check fails.** All fully reported frames in all four segment-role captures exceed the 62.5 ms render budget. Startup is reported separately from ongoing animation below.

## Setup

| | |
|---|---|
| Hardware | Two bench Teensy 4.0 boards at 600 MHz: COM3 serial 6362410 for roles 0/2; COM4 serial 6362390 for roles 1/3 |
| Image | `profile`, shipping `-Os` base, newlib-nano and selective `HS_O3` pullback shade, placed pipeline and color helpers; no deep or PC probes |
| Driver | `POVSegmented<288, 4, 480>`; each segment role runs standalone with its actual arm, latitude band, reversed-strip mapping, flywheel ISR and DMA output |
| Effect | KaleidoscopeMelt at 288×144, clean source `1ee8ff397bf45f8e8ec442ad0104aaafa97e03bf`, parameter commit `c4c1e04f7` |
| Method | Four independent 130 s captures, 32-frame windows, natural choreography and temporal sampling. No fixed preset, pause, replay or speed overrides. The harness default one-hour epoch prevents a reset during capture |
| Reproduce | Set `HS_PROFILE_SEGMENT_ID` to 0, 1, 2 and 3 through the locked wrapper; commands below |

The profile wrapper came from the existing uploader-repair checkout at `58b5d7193`, with `HS_PROFILE_TREE` explicitly naming the clean effect checkout. Its headless uploader pins COM to serial and USB location, verifies tycmd 0.9.9 against its SHA-256, and holds the per-board lock through upload and capture. The initial vendor-loader attempts failed before capture and are excluded.

The supplied legacy parameter tuple maps to the current named preset fields. All values already matched except `opacity_low = 0.339f` and `opacity_high = 0.865f`; the authored shader document and generated preset header agree.

| Role | Raw capture | Profile ELF SHA-256 | Captured (Los Angeles) |
|---:|---|---|---|
| 0 | `build/prof/kaleidoscopemelt_role0_20261030_ship.log` | `090e2d709191ca90ee72a64560e01b3ed363c5fb62e29b7f551b69c6618e5943` | 2026-10-30 23:35 |
| 1 | `build/prof/kaleidoscopemelt_role1_20261030_ship.log` | `1ce1900b51ce8b89a573eee89e9fa7cb6d12d0153c6e0ea6dc9254d2a78d28ab` | 2026-10-30 23:33 |
| 2 | `build/prof/kaleidoscopemelt_role2_20261030_ship.log` | `8a337408ec6e379e9ea298d512ed38353d563eee7abff0edf18a26bbc02747af` | 2026-10-30 23:37 |
| 3 | `build/prof/kaleidoscopemelt_role3_20261030_ship.log` | `04adef25ab17d8e5c51cb4cb068954946c3a9467e75a89ceb19aa4c1f7da63a3` | 2026-10-30 23:37 |

Every raw capture retains build/compiler/ABI attestation, exact source revision, build flags and immutable profile/Phantasm ELF copies under its provenance artifact directory.

| Role | FLASH code | FLASH data | RAM1 code | RAM1 variables | RAM1 free | RAM2 variables |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 101,880 | 151,672 | 35,432 | 313,056 | 145,696 | 520,064 |
| 1 | 101,880 | 151,672 | 35,432 | 313,056 | 145,696 | 520,064 |
| 2 | 101,880 | 151,672 | 35,432 | 313,056 | 145,696 | 520,064 |
| 3 | 101,880 | 151,672 | 35,432 | 313,056 | 145,696 | 520,064 |

All four profile images also have FLASH headers 8,592, RAM1 padding 30,104 and RAM2 free 4,224 bytes.

The accompanying shipping Phantasm image passes every region budget and layout gate: FLASH code 1,037,568, data 749,244, headers 8,260; RAM1 variables 312,736, code 192,600, padding 4,008, free 14,944; RAM2 variables 520,064, free 4,224 bytes. Against the pre-change `d7a7a86e8` build, RAM1 code, RAM1 variables and FLASH data deltas are all **0 bytes**. KaleidoscopeMelt remains outside the curated Phantasm playlist; the profiling image measures it explicitly.

Exactness cross-checks use the highest-cost clean window for each role:

| Role | Window frames | Root cycles / 600 (us) | Wall sum (us) | Difference (ppm) |
|---:|---|---:|---:|---:|
| 0 | 225–256 | 5999801.425 | 5999804 | 0.429 |
| 1 | 257–288 | 5998354.280 | 5998355 | 0.120 |
| 2 | 225–256 | 6000323.452 | 6000320 | 0.575 |
| 3 | 257–288 | 5998896.818 | 5998898 | 0.197 |

Each capture passes the parser: correct effect and resolution, continuous frames starting at 1, complete per-frame telemetry, no epoch reset, and root-cycle/wall agreement within 5 ppm. The logged segment banner and `profile role` marker separately match each requested role.

## Frame cadence

A half-revolution at 480 RPM lasts **62.5 ms**, giving a 16 fps render deadline. Each role shades 144×72 = **10,368 pixels**. The horizontal half alternates naturally with the driver; roles 0/2 cover northern rows, roles 1/3 cover southern rows, and roles 2/3 use arm B. Southern strips retain their reversed output order.

| Role | Mapping | Startup frame 1 (ms) | Ongoing peak (ms) | Peak frame | Ongoing spills | Full-pass spills |
|---:|---|---:|---:|---:|---:|---:|
| 0 | A, north +y | 93.505 | 84.662 | 253 | 671/671 | 672/672 |
| 1 | A, south −y | 100.416 | 87.506 | 293 | 671/671 | 672/672 |
| 2 | B, north +y | 93.480 | 86.154 | 243 | 671/671 | 672/672 |
| 3 | B, south −y | 100.442 | 88.522 | 253 | 671/671 | 672/672 |

README cells: peak 🔴 88.522, spilled 🔴 2684/2684 (100.00%).

All ongoing frames spill; startup is not the cause of the ongoing failure. Steady wall cadence is approximately 187.5 ms (three display windows, about 5.33 fps). `canvas_buffer_wait` is display-handoff idle; render telemetry subtracts it. Wall maxima are not substituted for render peaks. Unreported trailing frames are excluded.

## Phase-by-phase readout

The single preset stays active while camera, displacement, generated palette and noise hue modulation evolve naturally. Frame 1 includes initial preparation; subsequent windows share one ongoing rendering regime. Trees below show each role’s highest mean render window, with per-frame averages and original nesting. Percentages use the frame root; `fx_prepare_column` is included within `fx_shader_draw`.

### Role 0 — ongoing frames 225–256

```
frame                       187.494 ms 112.496 Mcyc 100.0%
  pov_preserve_half          139.82 us   83.89 kcyc   0.1% x1
  fx_shader_draw             77.022 ms  46.213 Mcyc  41.1%
    fx_prepare_column        519.64 us  311.79 kcyc   0.3% x144
  fx_prepare_frame            4.309 ms   2.586 Mcyc   2.3%
  fx_advance                 939.09 us  563.45 kcyc   0.5%
    palette_cycle_step       916.96 us  550.18 kcyc   0.5%
  fx_timeline_step            96.36 us   57.82 kcyc   0.1%
  canvas_clear                85.90 us   51.54 kcyc   0.0% x1
  canvas_buffer_wait        104.691 ms  62.815 Mcyc  55.8%
```

Wall min/avg/max = 184.939/187.493/190.231 ms. Across all ongoing frames, render min/avg/max = 75.629/80.477/84.662 ms. The shader dominates foreground work; even the least expensive reported ongoing frame exceeds the deadline.

### Role 1 — ongoing frames 257–288

```
frame                       187.449 ms 112.469 Mcyc 100.0%
  pov_preserve_half          139.26 us   83.55 kcyc   0.1% x1
  fx_shader_draw             79.431 ms  47.659 Mcyc  42.4%
    fx_prepare_column        518.57 us  311.14 kcyc   0.3% x144
  fx_prepare_frame            4.308 ms   2.585 Mcyc   2.3%
  fx_advance                 945.35 us  567.21 kcyc   0.5%
    palette_cycle_step       922.40 us  553.44 kcyc   0.5%
  fx_timeline_step            96.91 us   58.14 kcyc   0.1%
  canvas_clear                86.30 us   51.78 kcyc   0.0% x1
  canvas_buffer_wait        102.230 ms  61.338 Mcyc  54.5%
```

Wall min/avg/max = 184.495/187.448/190.378 ms. Across all ongoing frames, render min/avg/max = 75.222/80.863/87.506 ms. The shader dominates foreground work; even the least expensive reported ongoing frame exceeds the deadline.

### Role 2 — ongoing frames 225–256

```
frame                       187.510 ms 112.506 Mcyc 100.0%
  pov_preserve_half          139.98 us   83.99 kcyc   0.1% x1
  fx_shader_draw             76.991 ms  46.195 Mcyc  41.1%
    fx_prepare_column        521.27 us  312.76 kcyc   0.3% x144
  fx_prepare_frame            4.308 ms   2.585 Mcyc   2.3%
  fx_advance                 941.22 us  564.73 kcyc   0.5%
    palette_cycle_step       918.45 us  551.07 kcyc   0.5%
  fx_timeline_step            98.05 us   58.83 kcyc   0.1%
  canvas_clear                86.79 us   52.07 kcyc   0.0% x1
  canvas_buffer_wait        104.735 ms  62.841 Mcyc  55.9%
```

Wall min/avg/max = 181.180/187.510/194.495 ms. Across all ongoing frames, render min/avg/max = 76.275/80.471/86.154 ms. The shader dominates foreground work; even the least expensive reported ongoing frame exceeds the deadline.

### Role 3 — ongoing frames 257–288

```
frame                       187.466 ms 112.479 Mcyc 100.0%
  pov_preserve_half          139.69 us   83.82 kcyc   0.1% x1
  fx_shader_draw             79.414 ms  47.648 Mcyc  42.4%
    fx_prepare_column        512.49 us  307.49 kcyc   0.3% x144
  fx_prepare_frame            4.309 ms   2.586 Mcyc   2.3%
  fx_advance                 955.09 us  573.06 kcyc   0.5%
    palette_cycle_step       932.79 us  559.67 kcyc   0.5%
  fx_timeline_step            97.18 us   58.31 kcyc   0.1%
  canvas_clear                85.70 us   51.42 kcyc   0.0% x1
  canvas_buffer_wait        102.254 ms  61.352 Mcyc  54.5%
```

Wall min/avg/max = 182.041/187.465/192.466 ms. Across all ongoing frames, render min/avg/max = 75.242/80.869/88.522 ms. The shader dominates foreground work; even the least expensive reported ongoing frame exceeds the deadline.

### Per-pixel figures

The shader assigns 10,368 pixels per frame directly, with 1× quadrant coverage and temporal state prepared for 144 columns. Column preparation averages approximately 3.6 us/call in each selected window. This path has no `filter_blend` counter, so no blended-pixel count or cycles/blend claim is inferred. Shader cycles/pixel in the selected windows are:

| Role | Shader cycles/pixel | Column preparation calls/frame |
|---:|---:|---:|
| 0 | 4457.31 | 144 |
| 1 | 4596.69 | 144 |
| 2 | 4455.50 | 144 |
| 3 | 4595.71 | 144 |

## Column-ISR / DMA marshaling cost

ISR readings use the same selected windows. Pack and submit run inside the inclusive wake scope; they are not added to wake’s CPU share.

Role 0, ISR snapshot interval 6,001,802 us:

```
isr_wake          3457.25/f  0.50/ 1.96/18.02 us 3.61%
  isr_pack         432.12/f  8.72/ 9.45/12.34 us 2.17%
  isr_dma_submit   432.12/f  0.72/ 0.93/ 1.13 us 0.21%
```

Role 1, ISR snapshot interval 6,000,282 us:

```
isr_wake          3456.38/f  0.52/ 1.97/20.74 us 3.62%
  isr_pack         432.03/f  8.81/ 9.47/13.52 us 2.18%
  isr_dma_submit   432.03/f  0.72/ 0.93/ 2.79 us 0.21%
```

Role 2, ISR snapshot interval 6,002,359 us:

```
isr_wake          3457.59/f  0.50/ 1.96/22.42 us 3.61%
  isr_pack         432.19/f  8.75/ 9.46/13.51 us 2.17%
  isr_dma_submit   432.19/f  0.68/ 0.93/ 5.20 us 0.21%
```

Role 3, ISR snapshot interval 6,000,952 us:

```
isr_wake          3456.78/f  0.50/ 1.97/16.54 us 3.62%
  isr_pack         432.06/f  8.81/ 9.47/13.25 us 2.18%
  isr_dma_submit   432.06/f  0.64/ 0.93/ 2.03 us 0.21%
```

Columns are calls/frame, per-call min/avg/max and CPU share. DMA submit costs about one microsecond per column; packing costs several microseconds. The 600-byte framed column transfer runs asynchronously (about 230 us at the requested 24 MHz SPI clock) within the approximately 434 us column period. DMA-completion and other interrupts are not measured.

| Role | Inclusive wake CPU | Foreground budget/window (ms) | Peak / foreground budget |
|---:|---:|---:|---:|
| 0 | 3.61% | 60.244 | 1.405× |
| 1 | 3.62% | 60.237 | 1.453× |
| 2 | 3.61% | 60.244 | 1.430× |
| 3 | 3.62% | 60.237 | 1.470× |

Render already includes ISR interference; the foreground-budget ratios are conservative sizing figures. The unadjusted 62.5 ms deadline alone is sufficient to reject every role.

## Summary ranking

1. `fx_shader_draw` — 42.36% of wall frame, 79.414 ms/frame in role 3’s selected window.
2. `fx_prepare_frame` — 2.30% of wall frame, 4.309 ms/frame in role 3’s selected window.
3. `fx_advance` — 0.51% of wall frame, 0.955 ms/frame in role 3’s selected window.

The October 11 capture used an older pipeline and predates per-column temporal sampling; its 46.26 ms steady peak is historical context, not a matched control for the opacity change. No matched WASM/native performance comparison was run.

## Caveats

- Standalone timing simulates each physical segment’s rendering and output work on two bench boards. It does not test a four-board sync bus, downstream acquisition, cable timing or receiver-interrupt overhead.
- Cycle scopes include live ISR time because CYCCNT free-runs. ISR shares use their own snapshot intervals, including serial-report time.
- No per-pixel profiling scopes were added. Existing frame/column scopes have instrumentation overhead. `filter_blend` parenting can hide that subtree on other paths; this direct-assignment path does not call it.
- Selective `HS_O3` annotations in pullback shading, placed-stage evaluation and color helpers remain active; this is not global `-O3`.
- No dwell compression, parameter forcing or paused animation was used. The single-entry harness retains one effect instance throughout each 130 s capture.
- Captured firmware was built from clean isolated checkouts. Unrelated shared-checkout spec work was preserved and did not enter these firmware builds. The local simulator was installed from the landed code, with that shared spec work still present.

## Harness

`targets/Profile/Profile.ino` uses `HS_PROFILE_TARGET=KaleidoscopeMelt`, `HS_PROFILE_WINDOW=32` and `HS_PROFILE_SEGMENT_ID=<0..3>`. Explicit segment IDs select the existing standalone driver mode and disable the sync transmitter. Omitting the knob preserves strapped segment-0 profiling.

```bash
HS_PROFILE_TREE=/c/work/temp/holosphere-melt-20261030 \
HS_TEENSY_PORT=COM3 HS_UPLOAD_BACKEND=headless \
HS_PROFILE_OUT=/c/work/Holosphere/build/prof/kaleidoscopemelt_role0_20261030_ship.log \
bash /c/work/temp/Holosphere-uploader-repair-20261030/tools/profile_one.sh \
  KaleidoscopeMelt profile 130 32 -D HS_PROFILE_SEGMENT_ID=0
```

Repeat for IDs 1–3 with unique output paths; roles 1/3 used the second isolated checkout on COM4. The uploader wrapper revision is `58b5d7193`; the effect source revision is recorded independently in every capture. `just profile KaleidoscopeMelt` remains the default segment-0 entry point.

Validation: 125/125 native CTest cases passed with `HS_SMOKE_FRAMES=120`; excluded native targets built successfully; 15/15 composed-preset generator tests passed; profiling tools ran 335 tests successfully (one skipped); shipping Phantasm size/layout checks and all four role-image builds passed. The release WASM module was built and installed into the local daydream simulator; a targeted smoke check verified all 23 requested parameter values exactly as binary32 and rendered 120 frames at 288×144.
