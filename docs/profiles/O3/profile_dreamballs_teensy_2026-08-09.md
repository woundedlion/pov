# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-09, **-O3**)

See the [shipping report](../shipping/profile_dreamballs_teensy_2026-08-25.md)
for the current selective-O3 capture; its 2026-08-09 twin, taken alongside this
report against the 5-preset roster, has been replaced. Raw capture:
`build/prof/dreamballs_o3.log` in the implementation worktree.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile_o3`: global `-O3 -ffast-math`, single-effect image |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master, 72 px/arm |
| Effect | DreamBalls 288×144, single-entry playlist, capture tip unrecoverable — an implementation-worktree commit that was never landed, so no sha in this repository resolves to it; see the `99d03ae6` note below |
| Method | `HS_PROFILE`, window 16, 230 s, `HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile_o3 230 16 -D HS_PROFILE_EPOCH_REVS=2000` |

The landed `99d03ae6` differs only in setup-time flash placement.

Image size: `FLASH: code:136216, data:182420, headers:9044` / `RAM1:
variables:315328, code:65128, padding:408, free:143424` / `RAM2:
variables:520064, free:4224`. Exactness cross-check: frames 113–128 match root
cycles to wall time within **0.1 ppm**.

## Frame cadence

Peak frame render **42.94 ms**, spilled **0/3,648 frames (0%)**. Every preset
holds 16 fps with at least 19.56 ms of render margin.

## Phase-by-phase readout

The 230-second capture visits all five presets and proves two complete wraps.

### Peak direct weave — Rhombicuboctahedron (frames 3313–3328)

```
frame                    62.51 ms  37.51 Mcyc  100%
  db_timeline_step       36.32 ms  21.79 Mcyc   58%
    db_draw              36.26 ms  21.76 Mcyc   58%
      db_draw_scene      36.26 ms  21.76 Mcyc   58%
        db_mesh_plot     35.38 ms  21.23 Mcyc   57%  x18.0
          filter_blend    5.17 ms   3.10 Mcyc    8%  x49528
        db_warp_orient   411.0 us  246.6 kcyc    1%  x18.0
        db_displace      425.4 us  255.3 kcyc    1%  x18.0
  canvas_clear            87.8 us   52.7 kcyc    0%
  canvas_buffer_wait     26.10 ms  15.66 Mcyc   42%
```

Wall min/avg/max = 48.85/62.51/75.97 ms. Mesh plotting remains dominant;
transported endpoint work stays below 0.9 ms/frame.

### Per-preset table

| # | Solid | V/E/F | copies | blended px/f | `db_mesh_plot` ms | render ms | peak ms | fps |
|---|---|---|--:|--:|--:|--:|--:|--:|
| 1 | Rhombicuboctahedron | 24/48/26 | 18 | 51,082 | 36.51 | 37.56 | 42.94 | 16 |
| 4 | Icosidodecahedron | 30/60/32 | 10 | 35,648 | 25.30 | 26.09 | 31.41 | 16 |
| 2 | Rhombicosidodecahedron | 60/120/62 | 6 | 31,868 | 23.27 | 24.18 | 30.61 | 16 |
| 3 | Truncated Cuboctahedron | 48/72/26 | 6 | 26,486 | 19.82 | 20.88 | 27.12 | 16 |
| 5 | Snub Cube | 24/60/38 | 4 | 14,876 | 11.13 | 11.81 | 12.32 | 16 |

### Per-pixel figures

The peak clean hold blends 51,082 samples/frame, 4.93× quadrant coverage, at
about 63 cycles/blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.62/1.87/11.38 us  cpu 3.44%
isr_pack         144/frame  min/avg/max 6.44/7.15/9.50 us   cpu 1.64%
isr_dma_submit   144/frame  min/avg/max 0.67/0.91/1.13 us   cpu 0.21%
```

DMA submission remains subordinate to packing; ISR time is included in the
render scopes.

## Summary ranking

1. `db_mesh_plot` — 35.38 ms in the peak window.
2. `filter_blend` — 5.17 ms beneath mesh plotting.
3. Transported displacement plus warp/orientation — 0.84 ms combined.

## Caveats

- CYCCNT scopes include ISR time.
- `filter_blend` calls approximate blended samples.
- The epoch extension changes dwell capacity, not per-frame work.
- The final landed change after capture affects setup placement only.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=DreamBalls`, window 16.

## Global -O3 vs selective -O3

Worst-frame render improves from 44.65 ms to 42.94 ms (**1.04×**). Relative to
the shipping profile image, global O3 adds **25,976 B FLASH** and **16,272 B
ITCM**. Both configurations remain green with zero spilled frames.
