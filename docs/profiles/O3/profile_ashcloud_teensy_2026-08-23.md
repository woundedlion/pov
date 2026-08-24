# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-08-23, **-O3**)

Global-O3 reference twin of the [shipping selective-O3 capture](../shipping/profile_ashcloud_teensy_2026-08-23.md).
Raw capture: `build/prof/ashcloud_o3.log`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference image |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, tip `b203edd365fc` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh AshCloud profile_o3 70 32` |

Image size: `FLASH: code:84,240, data:151,604, headers:8,892` / `RAM1:
variables:315,104, code:38,008, padding:27,528, free:143,648` / `RAM2:
variables:520,064, free:4,224`.

Exactness cross-check: window frames 865–896 root counter cycles ÷ 600 MHz
match the measured wall sum within **3.0 ppm**.

## Frame cadence

**Pass aggregate**: `fx_shader_draw` reaches 41.95 ms/f at the slowest
window; peak frame render is **48.35 ms** at frames 513–544, with
**0/1,088** frames spilled (0%).

The 62.5 ms display window leaves 14.15 ms of worst-frame margin. The effect
renders one quadrant ≈ 10,368 samples and holds 16 fps throughout.

## Phase-by-phase readout

The single preset evolves continuously; frames 513–544 are the worst render
window.

### Worst window (frames 513–544)

```text
frame                  62.46 ms  37.48 Mcyc  100%
  fx_shader_draw       41.67 ms  25.00 Mcyc   66%
  fx_prepare_frame      3.36 ms   2.02 Mcyc    5%
  fx_advance            2.22 ms   1.33 Mcyc    3%
  fx_timeline_step     59.2 us   35.52 kcyc    0%
  canvas_clear         90.6 us   54.36 kcyc    0%
  canvas_buffer_wait   15.05 ms   9.03 Mcyc   24%
```

Wall min/avg/max = 60.52/62.46/64.43 ms. Global O3 reduces the composed draw
by 3.53 ms/frame and preparation by 0.51 ms/frame; that time returns as
display-sync idle.

### Per-pixel figures

The draw covers one 10,368-sample quadrant and spends about 2,411 cycles per
sample in `fx_shader_draw`, including the complete composed pipeline.

## Column-ISR / DMA marshaling cost

```text
isr_wake        1152/frame  min/avg/max 0.56/1.86/16.00 us  cpu 3.43%
isr_pack         144/frame  min/avg/max 6.10/7.04/9.40 us   cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.77/0.94/6.08 us   cpu 0.21%
```

- Packing dominates submit CPU cost; transfer remains asynchronous.
- ISR time is included in the measured render scopes.
- The 48.35 ms peak remains below one display window.

## Summary ranking

1. `fx_shader_draw` — 66% of the worst window, 41.67 ms/frame.
2. `fx_prepare_frame` — 5%, 3.36 ms/frame.
3. `fx_advance` — 3%, 2.22 ms/frame.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- No per-pixel profiling scope was enabled.
- Global O3 is a measurement reference, not a shippable roster image.
- The source tree had no tracked WIP; review/profile documents are ignored
  working-tree artifacts.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=AshCloud`,
`HS_PROFILE_WINDOW=32`.

Global O3 lowers peak render from 52.39 to 48.35 ms (**7.7%**) while adding
**14,344 B FLASH code** and **11,888 B ITCM code** over the shipping image.
