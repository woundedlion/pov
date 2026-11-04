# Involution on-device profile — Teensy 4.0, segmented mode (2026-11-04, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Involution`).
Raw capture: `build/prof/involution_ship.log` in the capture worktree,
captured **2026-11-04 00:15** America/Los_Angeles on COM3. Replaces the
2026-11-03 report, which measured the S2 projection.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 at 600 MHz, live flywheel and DMA ISRs |
| Image | `profile`, shipping `-Os` with selective `HS_O3`/ITCM shader, noise and color kernels |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Involution 288×144, authored `jewel-melt` preset, cubic atlas at Area Correction 0.7, source `15d4750f8` |
| Method | `HS_PROFILE`, 32-frame windows, 70 seconds, default epoch; no epoch reset |
| Reproduce | `HS_PROFILE_TREE=<tree> bash tools/profile_one.sh Involution profile 70 32` |

Full Phantasm image with Involution in the playlist: `FLASH: code:1081552,
data:788768`; `RAM1: variables:313760, code:194488`; `RAM2: variables:520064`.
The size and layout gate passes.

| Region | Master `9db2be662` (S2) | Final | Delta |
|---|--:|--:|--:|
| RAM1 code | 194,472 | 194,488 | +16 |
| RAM1 variables | 313,760 | 313,760 | 0 |
| FLASH code | 1,079,920 | 1,081,552 | +1,632 |
| FLASH data | 788,124 | 788,768 | +644 |

The RAM1 code delta is the `FieldGate::AREA_CORRECTION` case: 14
`field_gate_open<Spec>` instantiations resident in ITCM each grow 2 B, 16 B
net after alignment.

Exactness cross-check: root cycles/600 MHz differ from the measured wall sum
by **1.5 ppm** in frames 1-32, inside the validator's 5 ppm bound.

## Frame cadence

**Steady peak: 53.399 ms**, frame 616; mean **50.756 ms**,
**0/1087 spilled (0.00%)**, and **0/1087 at or above 58 ms**.
Frames 2-1088 in complete windows are included. Setup frame 1 is
**90.534 ms** and spills separately.

README cells: peak 🟢 53.399, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. Involution shades one 10,368-pixel
quadrant per frame and holds the 16 fps cadence. Setup frame 1 renders the
full canvas without a preserved half. `canvas_buffer_wait` is the
display-sync idle and is subtracted from render.

The S2 build measured a 48.380 ms steady peak and a 46.368 ms mean
(`18601fa2d`, COM4): +5.02 ms peak, +4.39 ms mean.

## Phase-by-phase readout

The single authored preset runs continuously: the camera and projection walk,
the periodic affine lattice turns, and palettes and noise clocks advance. There
are no preset markers.

### Worst steady window (frames 609-640)

```
frame                     62.403 ms 37.442 Mcyc 100.0%
  pov_preserve_half        0.137 ms  0.082 Mcyc   0.2%
  fx_shader_draw          49.722 ms 29.833 Mcyc  79.7%
    fx_prepare_column      5.195 ms  3.117 Mcyc   8.3% x144 36.08 us/col
  fx_palette_base_bake     0.092 ms  0.055 Mcyc   0.1%
  fx_prepare_frame         0.017 ms  0.010 Mcyc   0.0%
  fx_advance               0.956 ms  0.573 Mcyc   1.5%
    palette_cycle_step     0.865 ms  0.519 Mcyc   1.4%
  fx_timeline_step         0.045 ms  0.027 Mcyc   0.1%
  canvas_clear             0.087 ms  0.052 Mcyc   0.1%
  canvas_buffer_wait      11.058 ms  6.635 Mcyc  17.7%
```

Wall min/mean/max: **58.717/62.402/65.886 ms**; render mean/max
51.344/53.399 ms. Against the S2 report's worst steady window, shading outside
column preparation grows 41.80 → 44.53 ms (+2.73 ms), consistent with the
atlas kernel's six face dot products and the correction's two square roots
replacing S2's axis comparisons. `fx_prepare_column` grows 3.10 → 5.20 ms
(21.5 → 36.1 µs per column) although this change does not touch column
preparation; that part is unexplained and may be flash placement.

### Per-pixel figures

At 10,368 shaded pixels per frame, the 49.7 ms shader draw costs about
2,880 cycles per pixel including column preparation.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.58/2.01/23.84 us  cpu 3.70%
isr_pack         144/frame  min/avg/max 8.33/9.42/12.16 us  cpu 2.17%
isr_dma_submit   144/frame  min/avg/max 0.78/0.94/1.04 us  cpu 0.21%
```

- Packing costs about ten times the DMA submission per column.
- The wire transfer itself runs asynchronously on DMA.
- The inclusive `isr_wake` share is 3.70%, leaving about 60.2 ms of each
  62.5 ms window for render; the 53.4 ms steady peak needs no further speedup.
  Pack and submit run nested inside wake and are not added to it.
  DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — 79.7% of the frame, 49.7 ms: per-pixel shading.
2. `fx_prepare_column` — 8.3%, 5.2 ms: 144 interval samples and stage
   preparations per frame.
3. `fx_advance` / `palette_cycle_step` — 1.5%, 1.0 ms: palette rebakes on
   alternate frames.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` scope is active for this shader.
- The shader, simplex noise, lattice and color finishing run in the selective
  `HS_O3`/ITCM regions; the shader body and atlas projection run from flash.
- No dwell-compression knobs were used.
- Captured from the landing worktree at `15d4750f8` with no uncommitted
  changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Involution`,
`HS_PROFILE_WINDOW=32`; `just profile Involution` = build + flash + capture.
