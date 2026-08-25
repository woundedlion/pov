# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-09, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DreamBalls`). Raw capture:
`build/prof/dreamballs_ship.log` in the implementation worktree. This replaces
the 2026-08-05 report after the automatic weave gained parallel-transported
edge framing.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: shipping `-Os` plus selective-O3 raster/filter regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, 72 px/arm |
| Effect | DreamBalls 288×144, single-entry playlist, capture tip unrecoverable — an implementation-worktree commit that was never landed, so no sha in this repository resolves to it; see the `99d03ae6` note below |
| Method | `HS_PROFILE`, window 16, 230 s, `HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile 230 16 -D HS_PROFILE_EPOCH_REVS=2000` |

The landed `99d03ae6` differs from the capture only by placing setup-only solid
construction in flash; the automatic render loop measured here is unchanged.

Image size: `FLASH: code:110240, data:182268, headers:8548` / `RAM1:
variables:315296, code:48856, padding:16680, free:143456` / `RAM2:
variables:520064, free:4224`.

The final shipping roster image passes its FlexRAM gate at `RAM1 code:196504`,
104 bytes below the six-bank ceiling. Exactness cross-check: frames 3329–3344
match root cycles to wall time within **3.3 ppm**.

## Frame cadence

Pass aggregate: peak frame render **44.65 ms**, spilled **0/3,648 frames
(0%)**. Every preset remains below the 62.5 ms display window and holds 16 fps;
the worst frame retains 17.85 ms of render margin. `canvas_buffer_wait` is the
intentional wait for the next display flip.

## Phase-by-phase readout

Five 320-frame presets run in order. The capture contains 11 advance markers,
visits all five presets, and proves two complete wraps with no epoch reset.

### Peak direct weave — Rhombicuboctahedron (frames 3329–3344)

```
frame                    62.57 ms  37.54 Mcyc  100%
  db_timeline_step       41.98 ms  25.19 Mcyc   67%
    db_draw              41.90 ms  25.14 Mcyc   67%
      db_draw_scene      41.90 ms  25.14 Mcyc   67%
        db_mesh_plot     41.09 ms  24.65 Mcyc   66%  x18.0
          filter_blend    5.51 ms   3.31 Mcyc    9%  x55638
        db_warp_orient   366.0 us  219.6 kcyc    1%  x18.0
        db_displace      417.8 us  250.7 kcyc    1%  x18.0
  canvas_clear            88.1 us   52.9 kcyc    0%
  canvas_buffer_wait     20.50 ms  12.30 Mcyc   33%
```

Wall min/avg/max = 57.69/62.57/67.27 ms. Rasterization remains the dominant
cost. The transported framing stages together consume under 0.8 ms/frame in
the cadence-setting 18-copy preset.

### Per-preset table

Rows use each preset's worst clean hold; peak render is the worst single frame
owned by that preset. All five presets hold 16 fps.

| # | Solid | V/E/F | copies | blended px/f | `db_mesh_plot` ms | render ms | peak ms | fps |
|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 | Rhombicuboctahedron | 24/48/26 | 18 | 55,638 | 41.09 | 42.07 | 44.65 | 16 |
| 4 | Icosidodecahedron | 30/60/32 | 10 | 32,211 | 23.99 | 24.75 | 32.61 | 16 |
| 3 | Truncated Cuboctahedron | 48/72/26 | 6 | 26,396 | 20.72 | 21.75 | 26.27 | 16 |
| 2 | Rhombicosidodecahedron | 60/120/62 | 6 | 26,799 | 20.44 | 21.33 | 29.76 | 16 |
| 5 | Snub Cube | 24/60/38 | 4 | 15,002 | 11.69 | 12.36 | 13.00 | 16 |

### Per-pixel figures

The peak hold blends 55,638 samples/frame, 5.37× the 10,368-pixel quadrant,
at about 59 cycles/blend. The remaining `db_mesh_plot` time is geodesic setup,
sampling, clipping, and shading.

## Column-ISR / DMA marshaling cost

```
isr_wake        1154/frame  min/avg/max 0.74/2.00/11.48 us  cpu 3.68%
isr_pack         144/frame  min/avg/max 6.43/7.17/9.24 us   cpu 1.65%
isr_dma_submit   144/frame  min/avg/max 0.66/0.94/1.06 us   cpu 0.21%
```

ISR time is included in the foreground scopes. DMA submission remains much
smaller than packing, and the render retains its full 62.5 ms cadence budget.

## Summary ranking

1. `db_mesh_plot` — 41.09 ms in the peak window; wireframe sampling dominates.
2. `filter_blend` — 5.51 ms, already included beneath mesh plotting.
3. Transported displacement plus warp/orientation — 0.78 ms combined.

The previous shipping report peaked at 44.88 ms. This capture peaks at 44.65
ms, so the framing correction preserves measured cadence and render cost.

## Caveats

- CYCCNT scopes include ISR time.
- `filter_blend` is parented beneath `db_mesh_plot`; calls approximate blended samples.
- The capture used a 2,000-revolution epoch; this changes dwell capacity, not frame cost.
- The capture tip preceded only setup-code flash placement in the landed commit.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, window 16; the
reproduction command above performs the locked build, flash, and capture.
