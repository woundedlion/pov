# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/displacementfield_ship.log`. Captured 2026-11-04 22:56 America/Los_Angeles on COM4. Replaces profile_displacementfield_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | DisplacementField, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 121 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh DisplacementField profile 121 32` |

Image size:

```text
teensy_size:   FLASH: code:118640, data:151136, headers:8752   free for files:1753088
teensy_size:    RAM1: variables:314496, code:25512, padding:7256   free for local variables:177024
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `17b8672c3448c4f4097bdac100b19f877c1e69b28b925ec12b915a7b614948e6`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 737–768, root cyc/600 versus measured wall sum differs by **0.118 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   DisplacementField        worst 2 lines/set (0 sets over)  gated 20,480 B  overflow 32,770 B  full-set worst 4
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **44.911 ms** (frame 603), **0/1887 spilled** (0.00%). Setup frame 1: 18.817 ms. Complete windows cover frames 1–1888; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 37.706 ms/frame. `df_timeline_step` averages 37.406 ms/frame over complete windows; its worst window is 43.049 ms/frame (frames 737–768).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 11.551–44.911 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     60.112 ms 36.067 Mcyc 100.0%
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x0.97/frame; 144.556 us/call
  df_timeline_step                        29.276 ms 17.566 Mcyc  48.7%
    df_draw_rings                         29.146 ms 17.488 Mcyc  48.5%
      df_fused_scan                       14.327 ms  8.596 Mcyc  23.8%
        x0.78/frame; 18339.103 us/call
      df_lut_bake                          9.224 ms  5.534 Mcyc  15.3%
        df_hue_table_prep                  0.322 ms  0.193 Mcyc   0.5%
          x12.53/frame; 25.705 us/call
        df_octave_noise                    6.359 ms  3.815 Mcyc  10.6%
          x36.00/frame; 176.634 us/call
      df_chunk_cull                        0.568 ms  0.341 Mcyc   0.9%
        x36.75/frame; 15.465 us/call
      df_lut_bake                          0.616 ms  0.370 Mcyc   1.0%
        df_hue_table_prep                  0.006 ms  0.003 Mcyc   0.0%
          x2.62/frame; 2.133 us/call
      df_chunk_cull                        0.138 ms  0.083 Mcyc   0.2%
        x9.00/frame; 15.311 us/call
      df_fused_scan                        2.487 ms  1.492 Mcyc   4.1%
        ring_stack_table                   1.657 ms  0.994 Mcyc   2.8%
          x1.75/frame; 946.850 us/call
        ring_stack_table                   0.045 ms  0.027 Mcyc   0.1%
          x0.03/frame; 1435.315 us/call
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.180 us/call
  canvas_clear                             0.103 ms  0.062 Mcyc   0.2%
    x1.00/frame; 102.608 us/call
  canvas_buffer_wait                      30.565 ms 18.339 Mcyc  50.8%
    x1.00/frame; 30565.149 us/call
  df_prepare_fields                        0.001 ms  0.000 Mcyc   0.0%
    x1.00/frame; 0.613 us/call
```

Wall min/avg/max: 16.465/60.111/79.053 ms. Window mean render is 29.547 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 577–608)

```text
ring_stack_table                           1.768 ms  1.061 Mcyc   2.8%
  x2.00/frame; 883.759 us/call
frame                                     62.434 ms 37.460 Mcyc 100.0%
  pov_preserve_half                        0.144 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.783 us/call
  df_timeline_step                        42.164 ms 25.298 Mcyc  67.5%
    df_draw_rings                         41.991 ms 25.195 Mcyc  67.3%
      df_fused_scan                       25.003 ms 15.002 Mcyc  40.0%
        x1.00/frame; 25002.842 us/call
      df_lut_bake                         14.420 ms  8.652 Mcyc  23.1%
        df_hue_table_prep                  1.108 ms  0.665 Mcyc   1.8%
          x21.16/frame; 52.358 us/call
        df_octave_noise                    9.982 ms  5.989 Mcyc  16.0%
          x42.59/frame; 234.345 us/call
      df_chunk_cull                        0.724 ms  0.435 Mcyc   1.2%
        x43.91/frame; 16.496 us/call
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.178 us/call
  canvas_clear                             0.101 ms  0.061 Mcyc   0.2%
    x1.00/frame; 101.013 us/call
  canvas_buffer_wait                      19.994 ms 11.997 Mcyc  32.0%
    x1.00/frame; 19994.340 us/call
  df_prepare_fields                        0.000 ms  0.000 Mcyc   0.0%
    x1.00/frame; 0.351 us/call
```

Wall min/avg/max: 60.661/62.433/64.380 ms. Window mean render is 42.439 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 737–768)

```text
ring_stack_table                           1.837 ms  1.102 Mcyc   2.9%
  x2.00/frame; 918.490 us/call
frame                                     62.346 ms 37.407 Mcyc 100.0%
  pov_preserve_half                        0.145 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.536 us/call
  df_timeline_step                        43.049 ms 25.829 Mcyc  69.0%
    df_draw_rings                         42.868 ms 25.721 Mcyc  68.8%
      df_fused_scan                       25.734 ms 15.440 Mcyc  41.3%
        x1.00/frame; 25734.117 us/call
      df_lut_bake                         14.523 ms  8.714 Mcyc  23.3%
        df_hue_table_prep                  1.116 ms  0.670 Mcyc   1.8%
          x21.97/frame; 50.795 us/call
        df_octave_noise                   10.015 ms  6.009 Mcyc  16.1%
          x43.12/frame; 232.224 us/call
      df_chunk_cull                        0.741 ms  0.444 Mcyc   1.2%
        x45.28/frame; 16.356 us/call
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.177 us/call
  canvas_clear                             0.102 ms  0.061 Mcyc   0.2%
    x1.00/frame; 101.959 us/call
  canvas_buffer_wait                      19.019 ms 11.411 Mcyc  30.5%
    x1.00/frame; 19018.763 us/call
  df_prepare_fields                        0.000 ms  0.000 Mcyc   0.0%
    x1.00/frame; 0.352 us/call
```

Wall min/avg/max: 60.984/62.345/64.013 ms. Window mean render is 43.327 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 897–928)

```text
frame                                     62.788 ms 37.673 Mcyc 100.0%
  pov_preserve_half                        0.147 ms  0.088 Mcyc   0.2%
    x1.00/frame; 146.655 us/call
  df_timeline_step                        19.162 ms 11.497 Mcyc  30.5%
    df_draw_rings                         18.953 ms 11.372 Mcyc  30.2%
      df_fused_scan                       13.126 ms  7.876 Mcyc  20.9%
        x0.88/frame; 15001.428 us/call
      df_lut_bake                          2.704 ms  1.622 Mcyc   4.3%
        df_hue_table_prep                  0.072 ms  0.043 Mcyc   0.1%
          x9.38/frame; 7.691 us/call
      df_chunk_cull                        0.178 ms  0.107 Mcyc   0.3%
        x9.41/frame; 18.888 us/call
      df_lut_bake                          0.278 ms  0.167 Mcyc   0.4%
        df_hue_table_prep                  0.003 ms  0.002 Mcyc   0.0%
          x1.50/frame; 2.138 us/call
      df_chunk_cull                        0.062 ms  0.037 Mcyc   0.1%
        x4.38/frame; 14.281 us/call
      df_fused_scan                        1.254 ms  0.753 Mcyc   2.0%
        ring_stack_table                   1.604 ms  0.962 Mcyc   2.6%
          x1.88/frame; 855.546 us/call
    animation_noise_clock                  0.000 ms  0.000 Mcyc   0.0%
      x1.00/frame; 0.179 us/call
  canvas_clear                             0.098 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.082 us/call
  canvas_buffer_wait                      43.333 ms 26.000 Mcyc  69.0%
    x1.00/frame; 43333.449 us/call
  df_prepare_fields                        0.006 ms  0.003 Mcyc   0.0%
    x1.00/frame; 5.783 us/call
```

Wall min/avg/max: 58.395/62.787/63.827 ms. Window mean render is 19.455 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 1  | 29/28 | unavailable | 43.049 | 43.327 | 16.04 | modal-call-count hold |
| 2  | 30/30 | unavailable | 38.797 | 39.100 | 16.02 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 1 | 44.911 | 0/898 | 16.00 |
| 2 | 40.865 | 0/989 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1149.09/frame
  min/avg/max 0.533/1.918/29.568 us; cpu 3.53%
  isr_pack         143.62/frame
    min/avg/max 7.971/8.835/11.460 us; cpu 2.03%
  isr_dma_submit   143.62/frame
    min/avg/max 0.680/0.936/1.040 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1994837 us
isr_wake       n=36771 cyc min/avg/max=320/1151/17741 ns min/avg/max=533/1918/29568 total=70562 us cpu=3.53%
isr_pack       n=4596 cyc min/avg/max=4783/5301/6876 ns min/avg/max=7971/8835/11460 total=40611 us cpu=2.03%
isr_dma_submit n=4596 cyc min/avg/max=408/562/624 ns min/avg/max=680/936/1040 total=4309 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.53% leaves approximately 60.294 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `df_timeline_step` — 42.164 ms/frame, 67.53% of this root window.
2. `df_draw_rings` — 41.991 ms/frame, 67.26% of this root window.
3. `df_fused_scan` — 25.003 ms/frame, 40.05% of this root window.
4. `df_lut_bake` — 14.420 ms/frame, 23.10% of this root window.
5. `df_octave_noise` — 9.982 ms/frame, 15.99% of this root window.

README cells: peak 🟢 44.911 (2), spilled 🟢 0/1887 (0.00%).

README setup frame 1: 18.817 ms. Captured 2026-11-04 22:56 America/Los_Angeles.

This capture measures the generated hot-code layout against the 11-04 matched sweep; see the sweep scorecard for the delta.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=DisplacementField  windows=59
  [PASS] captured >=3 windows (59)
  [PASS] every window names one effect and resolution (DisplacementField 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [PASS] cycle returns to its first index (peak idx 2, 2 distinct)
  [PASS] all 2 presets visited (2)
  [PASS] per-frame render is render, not wall (0 of 59 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 59 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): df_chunk_cull, df_fused_scan, df_hue_table_prep, df_lut_bake, ring_stack_table
  [PASS] scope 'df_timeline_step' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.1 ppm, frames 737-768)

VALID: DisplacementField
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`. `just profile DisplacementField` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
