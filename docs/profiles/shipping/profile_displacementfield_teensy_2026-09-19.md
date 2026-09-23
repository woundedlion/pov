# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-09-19, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw captures: after: `build/prof/astra-regression-after-2026-09-19/build/prof/displacementfield_ship.log`, before: `build/prof/astra-regression-before-2026-09-19/build/prof/displacementfield_ship.log`. Replaces the historical 2026-08-26 capture (report no longer retained); the before/after comparison below uses fresh matched captures, not that older report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3 for all four passes; segmented POV, flywheel and DMA ISRs live |
| Image | `profile`; -Os base, shipping HS_O3 regions in DisplacementField field/ring/hue/LUT work and Scan::DistortedRingStack |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, arm A north, 72 physical pixels |
| Effect | DisplacementField 288×144, single-entry playlist, clean tip `8834627d96d747e5646c15f0c72f04681b2dcaac` |
| Method | HS_PROFILE cycle scopes; 70 s, window 32, default 120 s epoch, deep instrumentation off; no dwell, transition or seed override |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh DisplacementField profile 70 32` |

Image size:

- `FLASH: code:78944, data:150112, headers:8512   free for files:1794048`
- `RAM1: variables:315264, code:43848, padding:21688   free for local variables:143488`
- `RAM2: variables:520064  free for malloc/new:4224`

Exactness cross-check: frames 609–640, root 1,199,473,801 cycles ÷ 600 MHz versus measured wall sum 1,999,125 µs agree within **1.0 ppm** (`parse_profile.py ... validate`). Both before/after headers name DisplacementField and the expected configuration, begin at frame 1, and contain 34 complete monotonic windows with no epoch reset. Raw capture mtimes: before `2026-09-19T22:12:45.259485`, after `2026-09-19T22:17:30.294837` (local). Build identity, source status, ELF/map attestations and parser outputs are retained beside the archived logs.

## Frame cadence

**Pass aggregate:** `df_timeline_step` 47.231 ms/frame; worst window 55.581 ms/frame (frames 609–640); peak frame render **58.183 ms**, spilled **0/1088 (0%)**.

A half-revolution display window is 62.5 ms (16 fps). This non-persisting effect does not require a full frame: steady rendering covers one 144×72 quadrant, 10,368 pixels. All observed regimes held 16 fps. `canvas_buffer_wait` is deliberate idle until buffer ownership permits drawing. The first frame now seeds both halves of the segment (288×72, 20,736 pixels); subsequent frames preserve the opposite half.

| Measurement | Before | After | Change |
|---|--:|--:|--:|
| All-frame mean render, ms | 46.913726 | 47.461390 | +0.547664 |
| Matched frames 33–1088 mean render, ms | 47.495682 | 48.047037 | +0.551355 |
| Peak render, ms | 57.382000 | 58.183000 | +0.801000 |
| Startup frame 1 render, ms | 11.084000 | 20.810000 | +9.726000 |
| Spilled frames | 0/1088 | 0/1088 | 0 |
| Direct preservation scope, ms/call | absent | 0.139995 | added |

The matched comparison excludes the entire startup window and pairs identical frame IDs, preserving the same observed phase schedule. Whole-capture mean changes by 1.167%. The direct copy averages **0.139995 ms** across 1087 calls, **0.224%** of the display budget; frame 1 has no copy call. Each steady copy preserves 62,208 bytes in existing OCRAM pixel buffers; no extra pixel buffer is allocated. Startup adds 9.726 ms once for this effect entry. These are measured draw-frame startup times, not boot or effect-constructor timing.

## Phase-by-phase readout

Counter columns are time/frame, cycles/frame, percent of root frame, calls/frame, and ?s/call for leaves.

The initial noise field fades in for 150 frames, holds for 600, then fades out for 150. The following ball-spawning phase lasts 900 frames before draining. This 1088-frame capture covers the initial noise phase and early balls only; it does not establish the maximum over the complete ball phase.

### Noise fade-in (window frames 65–96)

```
frame                       62.72ms 37.63Mcyc 100.0%
  pov_preserve_half        141.20us 84.72kcyc   0.2% x   1.0   141.20us
  df_timeline_step          41.35ms 24.81Mcyc  65.9%
    df_draw_rings           41.28ms 24.77Mcyc  65.8%
      df_hue_table_prep    812.32us 487.39kcyc   1.3% x  19.8    40.94us
      df_lut_bake            8.92ms  5.35Mcyc  14.2% x  33.2   268.28us
      df_chunk_cull        958.69us 575.21kcyc   1.5% x  44.4    21.59us
      df_fused_scan         29.21ms 17.52Mcyc  46.6%
        filter_blend         1.13ms 676.05kcyc   1.8% x9455.8     0.12us
  canvas_clear              84.66us 50.80kcyc   0.1% x   1.0    84.66us
  canvas_buffer_wait        21.14ms 12.68Mcyc  33.7% x   1.0 21136.24us
  df_prepare_fields          0.24us  0.15kcyc   0.0% x   1.0     0.24us
```

Wall min/avg/max = 61.518/62.716/64.149 ms. Coverage and field strength are still rising. The wait remains larger than during the held noise field.

### Held noise field (window frames 609–640)

```
frame                       62.47ms 37.48Mcyc 100.0%
  pov_preserve_half        138.74us 83.25kcyc   0.2% x   1.0   138.74us
  df_timeline_step          55.58ms 33.35Mcyc  89.0%
    df_draw_rings           55.50ms 33.30Mcyc  88.8%
      df_hue_table_prep      1.76ms  1.06Mcyc   2.8% x  21.5    81.97us
      df_lut_bake           10.25ms  6.15Mcyc  16.4% x  43.2   237.11us
      df_chunk_cull          1.11ms 666.22kcyc   1.8% x  45.4    24.47us
      df_fused_scan         40.78ms 24.47Mcyc  65.3%
        filter_blend         1.21ms 728.22kcyc   1.9% x10109.4     0.12us
  canvas_clear              85.31us 51.18kcyc   0.1% x   1.0    85.31us
  canvas_buffer_wait         6.67ms  4.00Mcyc  10.7% x   1.0  6666.03us
  df_prepare_fields          0.24us  0.15kcyc   0.0% x   1.0     0.24us
```

Wall min/avg/max = 60.454/62.472/64.123 ms. The fused scan dominates the render work. The single worst render frame occurs in the preceding 577–608 window, so this representative hold window is not substituted for the exact pass peak.

### Noise fade-out (window frames 801–832)

```
frame                       62.18ms 37.31Mcyc 100.0%
  pov_preserve_half        140.59us 84.36kcyc   0.2% x   1.0   140.59us
  df_timeline_step          44.17ms 26.50Mcyc  71.0%
    df_draw_rings           44.07ms 26.44Mcyc  70.9%
      df_hue_table_prep      1.51ms 903.65kcyc   2.4% x  20.7    72.91us
      df_lut_bake            8.84ms  5.30Mcyc  14.2% x  32.7   270.46us
      df_chunk_cull        989.17us 593.50kcyc   1.6% x  39.5    25.04us
      df_fused_scan         31.45ms 18.87Mcyc  50.6%
        filter_blend         1.15ms 690.93kcyc   1.9% x9619.4     0.12us
  canvas_clear              84.91us 50.95kcyc   0.1% x   1.0    84.91us
  canvas_buffer_wait        17.78ms 10.67Mcyc  28.6% x   1.0 17782.97us
  df_prepare_fields          0.24us  0.15kcyc   0.0% x   1.0     0.24us
```

Wall min/avg/max = 59.287/62.175/65.263 ms. Noise amplitude is falling while rings continue scanning. Buffer wait expands as the render work drops.

### Early ball phase (window frames 1057–1088)

```
frame                       62.44ms 37.46Mcyc 100.0%
  pov_preserve_half        141.32us 84.79kcyc   0.2% x   1.0   141.32us
  df_timeline_step          42.52ms 25.51Mcyc  68.1%
    df_draw_rings           42.32ms 25.39Mcyc  67.8%
      df_hue_table_prep      1.70ms  1.02Mcyc   2.7% x  37.8    44.94us
      df_lut_bake           10.91ms  6.55Mcyc  17.5% x  42.9   254.15us
      df_chunk_cull          1.11ms 663.33kcyc   1.8% x  43.9    25.20us
      df_fused_scan         26.43ms 15.86Mcyc  42.3%
        filter_blend         1.10ms 661.06kcyc   1.8% x9237.6     0.12us
  canvas_clear              84.72us 50.83kcyc   0.1% x   1.0    84.72us
  canvas_buffer_wait        19.68ms 11.81Mcyc  31.5% x   1.0 19677.07us
  df_prepare_fields         11.85us  7.11kcyc   0.0% x   1.0    11.85us
```

Wall min/avg/max = 61.427/62.439/63.870 ms. Ball spawning has begun and the noise fade has completed. Later ball accumulation and draining are outside this capture.

### Per-pixel figures

| Regime window | Blends/frame | Quadrant coverage | Cycles/blend | Fused-scan cycles/blend |
|---|--:|--:|--:|--:|
| 65–96 | 9455.75 | 0.912× | 71.50 | 1853.35 |
| 609–640 | 10109.44 | 0.975× | 72.03 | 2420.07 |
| 801–832 | 9619.44 | 0.928× | 71.83 | 1961.53 |
| 1057–1088 | 9237.59 | 0.891× | 71.56 | 1716.37 |

Blends count actual filter calls, not all scanned pixels; fused scan also includes field evaluation and rejected coverage.

## Column-ISR / DMA marshaling cost

Representative held-noise window, frames 609–640:

```
isr_wake 1152.28/frame 0.563/1.718/15.196 us 3.16%
  isr_pack 144.03/frame 6.401/7.035/9.683 us 1.62%
  isr_dma_submit 144.03/frame 0.716/0.941/1.021 us 0.21%
```

Times are per-call min/average/max; percentages are CPU share over the dump interval.

- Wake includes pack and DMA submission; the nested shares must not be added again. Packing is the larger nested CPU cost.
- A 600-byte composite strobe transfer at 24 MHz occupies the SPI wire for 200 µs asynchronously, within the approximately 434 µs column interval. This is distinct from CPU submission time.
- Wake consumes 3.16% of elapsed time, leaving roughly 60.525 ms of foreground CPU per 62.5 ms window. Render counters already include ISR preemption, so do not subtract it again from measured render time. No speedup is required to meet the deadline in any observed phase.

## Summary ranking

1. `df_fused_scan` — 54.41% of frame wall time, 33.953 ms/frame averaged across the pass.
2. `df_lut_bake` — 15.08% of frame wall time, 9.412 ms/frame averaged across the pass.
3. `df_hue_table_prep` — 1.91% of frame wall time, 1.191 ms/frame averaged across the pass.
4. `df_chunk_cull` — 1.63% of frame wall time, 1.020 ms/frame averaged across the pass.
5. `pov_preserve_half` — 0.22% of frame wall time, 0.140 ms/frame averaged across the pass.

The fused scan remains the main rendering cost. These are live-device measurements; desktop copy benchmarks are not target timing estimates.

## Caveats

- CYCCNT free-runs through interrupts, so every foreground scope includes ISR preemption. Nested counters overlap and are not additive.
- `filter_blend` is parented under `df_fused_scan`; registry parenting can hide a subtree when its parent has zero calls. Per-pixel instrumentation itself has overhead. Deep scopes were disabled.
- Shipping uses selective HS_O3 field/ring/hue/LUT and DistortedRingStack regions; global O3 changes surrounding code too.
- Both trees were clean at capture. Before tip was `b2cefb6139e148c7e0c546b4d9d5d88c81645606`; after includes all 35 fixes and the committed copy scope. No dwell compression, epoch stretching, transition-speed or RNG seed overrides were used. The default harness seed path and frame IDs match; these are single captures, not a repeated-run confidence interval.
- The overall mean delta cannot be assigned entirely to the copy: the images also differ in other fixes, code layout and instrumentation. The direct copy scope isolates its inclusive cost more closely. Startup is a separate full-band render, not steady copying.
- The complete later ball phase was not sampled. This representative non-persisting clipped effect does not exercise the optional sine-distance path (#10) or extreme authored pullback chains (#20). Zero spills here is not proof that every effect or phase remains below budget.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`; `just profile DisplacementField` builds, flashes and captures the shipping image. Exact matched runs used:

```sh
HS_PROFILE_TREE=<before-or-after-tree> HS_TEENSY_PORT=COM3 \
  bash tools/profile_one.sh DisplacementField profile 70 32
```
