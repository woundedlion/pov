# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-09-19, **-O3**)

Shipping sibling: [selective -O3 report](../shipping/profile_displacementfield_teensy_2026-09-19.md).

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw captures: after: `build/prof/astra-regression-after-2026-09-19/build/prof/displacementfield_o3.log`, before: `build/prof/astra-regression-before-2026-09-19/build/prof/displacementfield_o3.log`. Replaces the historical 2026-08-26 capture (report no longer retained); the before/after comparison below uses fresh matched captures, not that older report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM3 for all four passes; segmented POV, flywheel and DMA ISRs live |
| Image | `profile_o3`; global -O3 -ffast-math reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, arm A north, 72 physical pixels |
| Effect | DisplacementField 288×144, single-entry playlist, clean tip `8834627d96d747e5646c15f0c72f04681b2dcaac` |
| Method | HS_PROFILE cycle scopes; 70 s, window 32, default 120 s epoch, deep instrumentation off; no dwell, transition or seed override |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM3 bash tools/profile_one.sh DisplacementField profile_o3 70 32` |

Image size:

- `FLASH: code:104800, data:149872, headers:8496   free for files:1768448`
- `RAM1: variables:315296, code:65992, padding:32312   free for local variables:110688`
- `RAM2: variables:520064  free for malloc/new:4224`

Exactness cross-check: frames 737–768, root 1,199,063,475 cycles ÷ 600 MHz versus measured wall sum 1,998,445 µs agree within **2.9 ppm** (`parse_profile.py ... validate`). Both before/after headers name DisplacementField and the expected configuration, begin at frame 1, and contain 34 complete monotonic windows with no epoch reset. Raw capture mtimes: before `2026-09-19T22:14:55.050687`, after `2026-09-19T22:19:46.075714` (local). Build identity, source status, ELF/map attestations and parser outputs are retained beside the archived logs.

## Frame cadence

**Pass aggregate:** `df_timeline_step` 46.586 ms/frame; worst window 54.605 ms/frame (frames 737–768); peak frame render **57.451 ms**, spilled **0/1088 (0%)**.

A half-revolution display window is 62.5 ms (16 fps). This non-persisting effect does not require a full frame: steady rendering covers one 144×72 quadrant, 10,368 pixels. All observed regimes held 16 fps. `canvas_buffer_wait` is deliberate idle until buffer ownership permits drawing. The first frame now seeds both halves of the segment (288×72, 20,736 pixels); subsequent frames preserve the opposite half.

| Measurement | Before | After | Change |
|---|--:|--:|--:|
| All-frame mean render, ms | 46.794940 | 46.815064 | +0.020124 |
| Matched frames 33–1088 mean render, ms | 47.378811 | 47.392081 | +0.013271 |
| Peak render, ms | 57.472000 | 57.451000 | -0.021000 |
| Startup frame 1 render, ms | 11.671000 | 21.567000 | +9.896000 |
| Spilled frames | 0/1088 | 0/1088 | 0 |
| Direct preservation scope, ms/call | absent | 0.139662 | added |

The matched comparison excludes the entire startup window and pairs identical frame IDs, preserving the same observed phase schedule. Whole-capture mean changes by 0.043%. The direct copy averages **0.139662 ms** across 1087 calls, **0.223%** of the display budget; frame 1 has no copy call. Each steady copy preserves 62,208 bytes in existing OCRAM pixel buffers; no extra pixel buffer is allocated. Startup adds 9.896 ms once for this effect entry. These are measured draw-frame startup times, not boot or effect-constructor timing.

## Phase-by-phase readout

Counter columns are time/frame, cycles/frame, percent of root frame, calls/frame, and ?s/call for leaves.

The initial noise field fades in for 150 frames, holds for 600, then fades out for 150. The following ball-spawning phase lasts 900 frames before draining. This 1088-frame capture covers the initial noise phase and early balls only; it does not establish the maximum over the complete ball phase.

### Noise fade-in (window frames 65–96)

```
frame                       62.71ms 37.62Mcyc 100.0%
  pov_preserve_half        140.27us 84.16kcyc   0.2% x   1.0   140.27us
  df_timeline_step          40.88ms 24.53Mcyc  65.2%
    df_draw_rings           40.82ms 24.49Mcyc  65.1%
      df_hue_table_prep    658.59us 395.15kcyc   1.1% x  19.8    33.19us
      df_lut_bake            8.73ms  5.24Mcyc  13.9% x  33.2   262.62us
      df_chunk_cull        881.90us 529.14kcyc   1.4% x  44.4    19.86us
      df_fused_scan         29.33ms 17.60Mcyc  46.8%
        filter_blend         1.00ms 600.34kcyc   1.6% x9455.7     0.11us
  canvas_clear              84.56us 50.74kcyc   0.1% x   1.0    84.56us
  canvas_buffer_wait        21.60ms 12.96Mcyc  34.4% x   1.0 21598.58us
  df_prepare_fields          0.25us  0.15kcyc   0.0% x   1.0     0.25us
```

Wall min/avg/max = 61.750/62.706/63.977 ms. Coverage and field strength are still rising. The wait remains larger than during the held noise field.

### Held noise field (window frames 609–640)

```
frame                       62.51ms 37.51Mcyc 100.0%
  pov_preserve_half        138.78us 83.27kcyc   0.2% x   1.0   138.78us
  df_timeline_step          54.20ms 32.52Mcyc  86.7%
    df_draw_rings           54.12ms 32.47Mcyc  86.6%
      df_hue_table_prep      1.60ms 962.65kcyc   2.6% x  21.8    73.56us
      df_lut_bake            9.41ms  5.64Mcyc  15.0% x  33.8   278.49us
      df_chunk_cull        911.91us 547.14kcyc   1.5% x  40.3    22.62us
      df_fused_scan         41.02ms 24.61Mcyc  65.6%
        filter_blend         1.08ms 648.06kcyc   1.7% x10127.6     0.11us
  canvas_clear              84.84us 50.90kcyc   0.1% x   1.0    84.84us
  canvas_buffer_wait         8.08ms  4.85Mcyc  12.9% x   1.0  8084.52us
  df_prepare_fields          0.25us  0.15kcyc   0.0% x   1.0     0.25us
```

Wall min/avg/max = 59.161/62.511/65.687 ms. The fused scan dominates the render work. The single worst render frame occurs in the preceding 577–608 window, so this representative hold window is not substituted for the exact pass peak.

### Noise fade-out (window frames 801–832)

```
frame                       62.13ms 37.28Mcyc 100.0%
  pov_preserve_half        139.44us 83.66kcyc   0.2% x   1.0   139.44us
  df_timeline_step          43.61ms 26.17Mcyc  70.2%
    df_draw_rings           43.54ms 26.12Mcyc  70.1%
      df_hue_table_prep      1.34ms 803.26kcyc   2.2% x  20.2    66.32us
      df_lut_bake            9.15ms  5.49Mcyc  14.7% x  41.1   222.58us
      df_chunk_cull        959.70us 575.82kcyc   1.5% x  42.7    22.47us
      df_fused_scan         30.76ms 18.46Mcyc  49.5%
        filter_blend         1.00ms 602.91kcyc   1.6% x9511.3     0.11us
  canvas_clear              84.53us 50.72kcyc   0.1% x   1.0    84.53us
  canvas_buffer_wait        18.30ms 10.98Mcyc  29.4% x   1.0 18295.50us
  df_prepare_fields          0.25us  0.15kcyc   0.0% x   1.0     0.25us
```

Wall min/avg/max = 60.584/62.133/63.747 ms. Noise amplitude is falling while rings continue scanning. Buffer wait expands as the render work drops.

### Early ball phase (window frames 1057–1088)

```
frame                       62.47ms 37.48Mcyc 100.0%
  pov_preserve_half        140.71us 84.42kcyc   0.2% x   1.0   140.71us
  df_timeline_step          42.08ms 25.25Mcyc  67.4%
    df_draw_rings           41.91ms 25.14Mcyc  67.1%
      df_hue_table_prep      1.24ms 744.52kcyc   2.0% x  33.4    37.14us
      df_lut_bake           10.75ms  6.45Mcyc  17.2% x  39.9   269.28us
      df_chunk_cull          1.03ms 618.89kcyc   1.7% x  46.9    22.00us
      df_fused_scan         27.01ms 16.21Mcyc  43.2%
        filter_blend       972.26us 583.35kcyc   1.6% x9212.6     0.11us
  canvas_clear              84.78us 50.87kcyc   0.1% x   1.0    84.78us
  canvas_buffer_wait        20.14ms 12.08Mcyc  32.2% x   1.0 20140.38us
  df_prepare_fields         11.61us  6.96kcyc   0.0% x   1.0    11.61us
```

Wall min/avg/max = 60.331/62.466/64.559 ms. Ball spawning has begun and the noise fade has completed. Later ball accumulation and draining are outside this capture.

### Per-pixel figures

| Regime window | Blends/frame | Quadrant coverage | Cycles/blend | Fused-scan cycles/blend |
|---|--:|--:|--:|--:|
| 65–96 | 9455.72 | 0.912× | 63.49 | 1861.15 |
| 609–640 | 10127.56 | 0.977× | 63.99 | 2430.16 |
| 801–832 | 9511.31 | 0.917× | 63.39 | 1940.40 |
| 1057–1088 | 9212.56 | 0.889× | 63.32 | 1759.31 |

Blends count actual filter calls, not all scanned pixels; fused scan also includes field evaluation and rejected coverage.

## Column-ISR / DMA marshaling cost

Representative held-noise window, frames 609–640:

```
isr_wake 1153.06/frame 0.396/1.591/11.880 us 2.93%
  isr_pack 144.12/frame 6.143/6.936/9.475 us 1.59%
  isr_dma_submit 144.12/frame 0.681/0.930/1.013 us 0.21%
```

Times are per-call min/average/max; percentages are CPU share over the dump interval.

- Wake includes pack and DMA submission; the nested shares must not be added again. Packing is the larger nested CPU cost.
- A 600-byte composite strobe transfer at 24 MHz occupies the SPI wire for 200 µs asynchronously, within the approximately 434 µs column interval. This is distinct from CPU submission time.
- Wake consumes 2.93% of elapsed time, leaving roughly 60.669 ms of foreground CPU per 62.5 ms window. Render counters already include ISR preemption, so do not subtract it again from measured render time. No speedup is required to meet the deadline in any observed phase.

## Summary ranking

1. `df_fused_scan` — 54.49% of frame wall time, 34.005 ms/frame averaged across the pass.
2. `df_lut_bake` — 14.68% of frame wall time, 9.160 ms/frame averaged across the pass.
3. `df_hue_table_prep` — 1.63% of frame wall time, 1.019 ms/frame averaged across the pass.
4. `df_chunk_cull` — 1.51% of frame wall time, 0.944 ms/frame averaged across the pass.
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
  bash tools/profile_one.sh DisplacementField profile_o3 70 32
```

**Global -O3 vs selective -O3:** mean render 47.461390 → 46.815064 ms (1.0138×); exact peak 58.183 → 57.451 ms. Copy cost 139.995 → 139.662 µs/call. FLASH code grows by +25,856 B and ITCM code by +22,144 B. Global O3 crosses a 32 KiB ITCM allocation boundary here: RAM1 local-variable headroom is 110,688 B versus shipping 143,488 B.
