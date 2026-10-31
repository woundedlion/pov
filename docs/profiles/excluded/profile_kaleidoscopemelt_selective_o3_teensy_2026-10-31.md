# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-31, **selective -O3**)

**All four ongoing segment roles meet the 58 ms target.** The worst observed render is **55.668 ms**, with **0/8188 ongoing frame spills** against the 62.5 ms display deadline. Temporal sampling remains enabled. This is a 37.11% reduction from the October 30 worst of 88.522 ms.

Startup frame 1 still exceeds the display deadline: the driver primes the complete horizontal width of its latitude band, shading 20,736 pixels instead of the ongoing half-width 10,368. Startup and full-pass counts are explicit below; “no spills” here means the ongoing natural animation.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Two Teensy 4.0 boards at 600 MHz: COM3 serial 6362410 for roles 0/2; COM4 serial 6362390 for roles 1/3 |
| Driver | `POVSegmented<288, 4, 480>`; actual arm, latitude clip, reversed southern strips, flywheel ISR and DMA output; standalone timing for each physical role |
| Image | Shipping `profile`: `-Os` base, newlib-nano, selective `HS_O3`/ITCM helpers; no deep scopes or PC sampler in these four captures |
| Source | Clean `0586276c44359c972d462dd87e4fb3d2e031f4b1`, integrated over `87a5fb2e1` |
| Run | Four independent 130 s captures, 32-frame windows, 288×144, natural choreography, temporal sampling ON, no fixed preset, pause, speed or clock overrides; default one-hour epoch |

The requested 23-field legacy tuple is installed in the authored shader document and generated preset. Only opacity endpoints differed from the prior preset: 0.339 and 0.865. CUP mapping, camera motion, curl displacement and noise hue modulation retain the supplied values.

The locked wrapper is the uploader-repair revision `58b5d7193`, with `HS_PROFILE_TREE` explicitly selecting the measured checkout. It pins tycmd 0.9.9 and the COM/serial/USB identity through upload and capture. The equivalent headless uploader has since landed in the source tree. All raw captures retain compiler/ABI/flag attestation and immutable ELF/map/build/source artifacts.

| Role | Raw capture | Profile ELF SHA-256 |
|---:|---|---|
| 0 | `build/prof/kaleidoscopemelt_integrated_final_role0_20261031_ship.log` | `a0f5f213675a522327690d30c6fc7a67c45daad2f44f15eadfff47b1766f630e` |
| 1 | `build/prof/kaleidoscopemelt_integrated_final_role1_20261031_ship.log` | `601b8ce57365f501d7f9521d5209c9f5e88253fbead795c9e17c793f57e1f63f` |
| 2 | `build/prof/kaleidoscopemelt_integrated_final_role2_20261031_ship.log` | `d228243d2f42d0c9e1ab4294257bf4af542b7c441322c6ea8b39a7a4203e9749` |
| 3 | `build/prof/kaleidoscopemelt_integrated_final_role3_20261031_ship.log` | `137fc34f460c2ca8e592ac5300b88421ac284faaa4cd9de925078365261a2d90` |

## Frame cadence

| Role | Mapping | October 30 peak (ms) | Final ongoing peak (ms) | Peak frame | Startup (ms) | Ongoing spills | Full-pass spills |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0 | A north | 84.662 | 54.917 | 1613 | 76.862 | 0/2047 | 1/2048 |
| 1 | A south | 87.506 | 54.635 | 1806 | 82.212 | 0/2047 | 1/2048 |
| 2 | B north | 86.154 | 55.668 | 1952 | 76.734 | 0/2047 | 1/2048 |
| 3 | B south | 88.522 | 55.493 | 1785 | 82.254 | 0/2047 | 1/2048 |

README cells: peak 🟢 55.668, spilled 🟢 0/8188 (0.00%).

Every capture has continuous, complete frame telemetry starting at 1, the correct effect/resolution and matching segment banner/standalone marker, with no epoch reset. Unreported trailing frames are excluded. Render excludes `canvas_buffer_wait` display-handoff idle and includes live ISR interference; wall maxima are not used as render peaks.

| Role | Highest mean ongoing window | Root cycles / 600 (us) | Wall sum (us) | Difference (ppm) |
|---:|---|---:|---:|---:|
| 0 | 1921–1952 | 1997309.888 | 1997309 | 0.445 |
| 1 | 1761–1792 | 1996585.410 | 1996587 | 0.796 |
| 2 | 1921–1952 | 1997801.053 | 1997798 | 1.528 |
| 3 | 1761–1792 | 1996746.527 | 1996743 | 1.766 |

## Selected implementation

- Fused temporal palette/modulation endpoints and in-place carrier blending remove repeated copies, palette construction and color finalization. Both temporal endpoints and the per-column sampling surface remain active.
- Paired modulation cubes share direction admission and tap coordinates. Endpoint-one shading samples only the current endpoint. Unused hue-rotation tables are omitted for Melt.
- A 256-value × 3-time-row OKLab base cache has two retained owners. It reuses endpoint samples/equal colors during baking and restores the canonical RGB-first base near the discontinuous gamut admission gate.
- The exact full-resolution gamut table is copied into the existing persistent arena. Late palette phases otherwise stall on FLASH table loads. This costs 131,072 arena bytes and was the decisive late-phase improvement.
- Exact range fast paths remove unnecessary periodic/longitude libm calls. Fused unit hues use bounded wrapping. Square lattice endpoints and interior cells skip unused distance work; canonical dodecahedral folds skip a redundant loop/axis test.
- Scalar cube baking inlines callees locally inside its FLASH loop. Gradient sampling keeps the integrated value-return APIs, replaces qualified floating-point control work with integer control, and reads an identical private 1,024-byte gradient table from DTCM.

The Lab cache is admitted only when **both** endpoint color and modulation configurations match the authored recipe. Changed controls, and the first transition restoring them, take the canonical path. Melt does not enable an approximate final RGB/hue cache or a geometry/noise interpolation cache. The independent AshCloud RGB/chamber caches remain intact.

## Accuracy qualification

The finite native cache corpus covers two seeds (0 and 123456789), five chromas, 108 adjacent natural palette pairs per case, 512 field values, 256 hues and nine column fractions: **3,822,059,520 channel samples**. The selected 0.001 linear-RGB gamut guard gives worst raw RGB error 5,139 codes, mean at most 0.3945 code, and worst premultiplied framebuffer error 3,228 codes for the supplied opacity recipe. Existing limits remain 7,000 maximum RGB, 256 mean RGB and 4,000 framebuffer codes; none was relaxed. Maximum observed fallback rate is 4.337%.

Qualification data: `build/prof/kaleidoscopemelt_lab_cache_qualification_20261031.json`, its CSV and guard CSV. This is observed finite-domain coverage, not a proof over continuous inputs or other control recipes. The both-endpoint admission guard confines the approximation to the tested authored controls.

Longitude, bounded hue interpolation, square-distance shortcuts and local cube baking have independent old/new parity checks. The cube-bake check compares 8,294,400 node bytes and the same number of direct floats. DTCM placement preserves all 1,024 lookup bytes and emitted gradient instructions. The integer-control port is exact in normal Debug/Release comparisons; fast-math can reassociate gradient arithmetic (observed value-component difference up to 4.30e-6). Existing renderer qualification tolerances are unchanged.

At the integration base, an optional native fast-math run failed 10,023 strict equality checks in the existing value-return gradient test. That test and its tolerances were preserved. Normal Debug/Release suites and independent fast-math output comparisons were used for this port.

## Code generation and memory

The retained final audit is `build/prof/kaleidoscopemelt_final_codegen_20261031.txt`, tied to the exact role-1 profile and shipping Phantasm ELFs. The transformed gradient leaf is 944 bytes, uses integer control fast paths, and has no explicit local stack load/store sites. Fused color still contains branch-dependent stack traffic and fallback copies. This is not a claim of zero register spills or a whole-program stack bound.

An early global placement experiment was removed after integration: it emitted an unnecessary 1,188-byte MermaidSkin endpoint and a 1,096-byte KaleidoscopeHexOil endpoint in ITCM. Melt already uses `INLINE_ONLY`; inherited dispatch still forwards temporal endpoints. The removal recovered 2,288 linked bytes, preserving the other effects’ previous placement. Temporal color behavior is tested for all three placement modes.

| Shipping Phantasm region | Integrated baseline `87a5fb2e1` | Final | Delta |
|---|---:|---:|---:|
| RAM1 code | 192,568 | 195,784 | +3,216 |
| RAM1 variables | 312,736 | 313,760 | +1,024 |
| FLASH data | 779,156 | 779,156 | 0 |

Final Phantasm FLASH code/data/headers are 1,061,416 / 779,156 / 8,772 bytes. RAM1 padding is 824 bytes, with 13,920 bytes available for locals; RAM2 variables/free are 520,064 / 4,224. Code remains below the unchanged 196,608-byte limit; the 12,288-byte stack floor and 272,384-byte persistent-footprint limit are unchanged. The full gamut copy consumes existing arena capacity, not additional global arena allocation. Its device-only allocation is charged by the compile-time persistent-footprint check; the native arena probe does not allocate that copy. The two Lab node arrays total 18,432 bytes; the 1,024-byte gradient table is the static RAM1-variable increase.

Against the original October 30 source, final RAM1 code is +3,184 bytes, variables +1,024 and FLASH data +29,912. The FLASH-data increase belongs to the intervening shared AshCloud assets; against the actual integration base its delta is zero.

| Role image | FLASH code | FLASH data | RAM1 code | RAM1 variables | RAM1 free |
|---:|---:|---:|---:|---:|---:|
| 0 | 111,280 | 152,016 | 37,864 | 314,112 | 144,640 |
| 1 | 111,280 | 152,016 | 37,864 | 314,112 | 144,640 |
| 2 | 111,280 | 152,016 | 37,864 | 314,112 | 144,640 |
| 3 | 111,280 | 152,016 | 37,864 | 314,112 | 144,640 |

## Diagnostic trials

Trial timing uses clean natural captures unless marked diagnostic. Windows differ with achieved frame count, so individual deltas are observations, not additive predictions. All final decisions are rechecked in the four captures above.

| Trial | Ongoing peak (ms) | Ongoing spills | Decision |
|---|---:|---:|---|
| Original requested parameters, worst role | 88.522 | 671/671 | Baseline failure |
| Guarded Lab, direct finisher and geometry/gradient work, all-four checkpoint | 62.849 | 3/8188 | Continued optimization |
| Exact float modulation-node copies (+28,912 arena bytes) | 62.832 | 3/2047 | Rejected; negligible gain |
| Global O3 control | 62.825 | 3/2047 | Rejected; selective flags retained |
| Exact longitude range shortcut | 61.668 | 0/2047 | Retained |
| Exact full gamut copy, longitude source | 56.468 | 0/2047 | Retained, then integrated |
| Integer controls + unit hue + square interior + local cube bake, no full gamut copy | 59.429 | 0/2047 | Combined with exact gamut copy |
| Final integrated, all four roles | 55.668 | 0/8188 | Below 58 ms |

Other rejected experiments include a final RGB/hue cache (7.7–9.6 ms rebakes outweighed shader savings), raw homogeneous-cube coordinates, blanket gradient unrolling, a curl chamber approximation, and a global noise-sample inline annotation that failed an existing equality check. Temporal sampling was never removed from a shipping candidate. An early non-temporal diagnostic was retired.

The long PC sample at source `7fec51672` covers 137,059 ms and 1,117,824 render samples: generated color 25.1%, geometry pipeline 16.5%, analytic gradient 16.5%, surface integration 5.3%, projection 5.2%, lattice 4.8%, gamut mapping 3.8%, scalar simplex 3.6%, and libm remainder calls 2.6%. It located the late-phase table stall and control overhead. Its instrumentation is excluded from final timing. Data: `build/prof/kaleidoscopemelt_selected_worst_pc_20261031_ship_pc_report.txt`.

The prior October 11 green result predates the October 17 column-time work and later modulation changes. The later sweep covered 36 curated Phantasm effects; Melt was excluded by curation. An old-opacity control (84.474 ms versus 84.662 ms for the requested opacity) ruled out the parameter edit as the cause. No single regression commit was established.

## Phase and ISR breakdown

Trees use each role’s highest mean ongoing window. Times are per-frame averages; percentages use the inclusive frame root, and column preparation is nested in shader draw. Each ongoing frame shades 10,368 pixels directly; no filter/blended-pixel count is inferred.

### Role 0 — frames 1921–1952

```
frame                        62.416 ms  37.450 Mcyc 100.0%
  pov_preserve_half          137.07 us   82.24 kcyc   0.2%
  fx_shader_draw             48.524 ms  29.114 Mcyc  77.7%
    fx_prepare_column        568.39 us  341.04 kcyc   0.9%
  fx_palette_base_bake       349.44 us  209.67 kcyc   0.6%
  fx_prepare_frame            3.001 ms   1.800 Mcyc   4.8%
  fx_advance                 961.59 us  576.96 kcyc   1.5%
    palette_cycle_step       934.83 us  560.90 kcyc   1.5%
  fx_timeline_step            83.88 us   50.33 kcyc   0.1%
  canvas_clear                86.95 us   52.17 kcyc   0.1%
  canvas_buffer_wait          9.084 ms   5.450 Mcyc  14.6%
```

Ongoing render min/mean/max = 44.829/49.463/54.917 ms. Shader cycles/pixel = 2808.08. Column preparations/frame = 144.

ISR snapshot interval: 2,001,039 us. Calls/frame and per-call min/mean/max:

```
isr_wake         1152.66/f   0.54/ 2.02/23.10 us  3.72% CPU
isr_pack          144.06/f   8.75/ 9.48/12.13 us  2.18% CPU
isr_dma_submit    144.06/f   0.65/ 0.94/ 1.15 us  0.21% CPU
```

### Role 1 — frames 1761–1792

```
frame                        62.393 ms  37.436 Mcyc 100.0%
  pov_preserve_half          139.10 us   83.46 kcyc   0.2%
  fx_shader_draw             48.672 ms  29.203 Mcyc  78.0%
    fx_prepare_column        561.07 us  336.64 kcyc   0.9%
  fx_palette_base_bake       347.17 us  208.30 kcyc   0.6%
  fx_prepare_frame            3.002 ms   1.801 Mcyc   4.8%
  fx_advance                 937.82 us  562.69 kcyc   1.5%
    palette_cycle_step       911.25 us  546.75 kcyc   1.5%
  fx_timeline_step            84.56 us   50.74 kcyc   0.1%
  canvas_clear                87.54 us   52.53 kcyc   0.1%
  canvas_buffer_wait          8.926 ms   5.356 Mcyc  14.3%
```

Ongoing render min/mean/max = 44.905/49.630/54.635 ms. Shader cycles/pixel = 2816.68. Column preparations/frame = 144.

ISR snapshot interval: 2,000,286 us. Calls/frame and per-call min/mean/max:

```
isr_wake         1152.22/f   0.51/ 2.02/21.23 us  3.72% CPU
isr_pack          144.00/f   8.79/ 9.52/12.29 us  2.19% CPU
isr_dma_submit    144.00/f   0.61/ 0.94/ 1.12 us  0.21% CPU
```

### Role 2 — frames 1921–1952

```
frame                        62.431 ms  37.459 Mcyc 100.0%
  pov_preserve_half          135.89 us   81.53 kcyc   0.2%
  fx_shader_draw             48.529 ms  29.117 Mcyc  77.7%
    fx_prepare_column        568.29 us  340.98 kcyc   0.9%
  fx_palette_base_bake       349.65 us  209.79 kcyc   0.6%
  fx_prepare_frame            3.001 ms   1.801 Mcyc   4.8%
  fx_advance                 964.96 us  578.98 kcyc   1.5%
    palette_cycle_step       937.65 us  562.59 kcyc   1.5%
  fx_timeline_step            85.12 us   51.07 kcyc   0.1%
  canvas_clear                87.25 us   52.35 kcyc   0.1%
  canvas_buffer_wait          9.084 ms   5.450 Mcyc  14.6%
```

Ongoing render min/mean/max = 45.230/49.459/55.668 ms. Shader cycles/pixel = 2808.39. Column preparations/frame = 144.

ISR snapshot interval: 2,001,577 us. Calls/frame and per-call min/mean/max:

```
isr_wake         1152.97/f   0.51/ 2.02/20.22 us  3.72% CPU
isr_pack          144.12/f   8.73/ 9.48/12.10 us  2.18% CPU
isr_dma_submit    144.12/f   0.61/ 0.94/ 1.08 us  0.21% CPU
```

### Role 3 — frames 1761–1792

```
frame                        62.398 ms  37.439 Mcyc 100.0%
  pov_preserve_half          137.01 us   82.21 kcyc   0.2%
  fx_shader_draw             48.684 ms  29.210 Mcyc  78.0%
    fx_prepare_column        561.93 us  337.16 kcyc   0.9%
  fx_palette_base_bake       350.08 us  210.05 kcyc   0.6%
  fx_prepare_frame            3.001 ms   1.800 Mcyc   4.8%
  fx_advance                 943.74 us  566.24 kcyc   1.5%
    palette_cycle_step       916.20 us  549.72 kcyc   1.5%
  fx_timeline_step            82.77 us   49.66 kcyc   0.1%
  canvas_clear                88.00 us   52.80 kcyc   0.1%
  canvas_buffer_wait          8.920 ms   5.352 Mcyc  14.3%
```

Ongoing render min/mean/max = 44.794/49.634/55.493 ms. Shader cycles/pixel = 2817.37. Column preparations/frame = 144.

ISR snapshot interval: 2,000,463 us. Calls/frame and per-call min/mean/max:

```
isr_wake         1152.34/f   0.51/ 2.03/18.39 us  3.73% CPU
isr_pack          144.03/f   8.80/ 9.52/12.40 us  2.19% CPU
isr_dma_submit    144.03/f   0.65/ 0.94/ 1.08 us  0.21% CPU
```

Pack and DMA submit are included in wake, not added to it. The framed 600-byte transfer is asynchronous; DMA-completion and other interrupts are not separately measured. Cycle scopes already include ISR interference.

## Summary ranking

Role 3's highest mean ongoing window spends 48.684 ms/frame in `fx_shader_draw`, 3.001 ms in `fx_prepare_frame`, and 0.944 ms in `fx_advance`. Shader work remains dominant; the separate Lab-base bake costs 0.350 ms/frame.

## Harness

Final-source Debug default and excluded targets built; **125/125 CTest cases passed** with `HS_SMOKE_FRAMES=120`. Selected Release tests passed for Melt, pullback, composed effects, noise fields, color, projections and interpolation. The affected gamut generator suite passed eight tests. Firmware region/layout gates and all four role images passed. Rebased source formatting/lint/license/include/documentation checks passed; sibling links were checked against the actual daydream checkout.

The landed release WASM engine was built and installed into the local daydream simulator. A smoke check verified all 23 supplied values exactly as binary32 and rendered 120 natural temporal frames at 288×144. Shared spec work and unrelated daydream video files were preserved.

The local bundle records `0586276c44359c972d462dd87e4fb3d2e031f4b1-dirty` because the shared checkout retains unrelated spec work. The device captures use the clean source revision above; the installed engine is a local development snapshot.

```bash
HS_PROFILE_TREE=/c/work/temp/holosphere-melt-final-20261031 \
HS_TEENSY_PORT=COM4 HS_UPLOAD_BACKEND=headless HS_DEVICE_WAIT=600 \
HS_PROFILE_OUT=/c/work/Holosphere/build/prof/kaleidoscopemelt_integrated_final_role3_20261031_ship.log \
bash /c/work/temp/Holosphere-uploader-repair-20261030/tools/profile_one.sh \
  KaleidoscopeMelt profile 130 32 "-D HS_PROFILE_SEGMENT_ID=3"
```

Repeat with IDs 0–3 and unique output paths. Role IDs select the actual standalone driver mapping and disable sync transmission. These runs simulate four physical roles on two boards; they do not validate a four-board sync bus, cables or receiver-interrupt overhead. The captures are finite natural-run observations, not a guarantee for every future palette or arbitrary edited parameters.
