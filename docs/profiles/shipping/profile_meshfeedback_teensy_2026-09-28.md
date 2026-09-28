# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MeshFeedback`).
Raw capture: `build/prof/meshfeedback_ship_itcm_levers.log` (A side:
`build/prof/meshfeedback_ship_itcm_base.log`). Replaces
`profile_meshfeedback_teensy_2026-08-26.md`. This capture is the B side of a
same-board A/B against the preceding master tip `54d0eb7fa`, run to check the
five ITCM-reduction commits `181b67168..a6e9d6426`.

## Update: tip `0e45a1fd9` (COM4, 14:46)

Raw capture: `build/prof/meshfeedback_ship_axes_flush.log`. Peak render
**61.85 ms** (frames 4097–4112), spilled **0/6688**. The later ITCM commits
change the random-walk trajectory, so the peak window moved from frames
1201–1216. Same-board (COM4) peaks at each step:

| Tip | Peak ms | Mesh draw ms/f | Blends/f |
|---|--:|--:|--:|
| `574d58674` | 62.24 | 6.68 | 9,049 |
| `574d58674` with finding 1 reverted | 60.02 | 4.88 | 8,391 |
| `83d706add` per-stroke screen-step axes | 61.13 | 5.76 | 9,049 |
| `0e45a1fd9` + flush setup in flash | 61.85 | 5.78 | 9,049 |

The flush commit costs 0.6–0.7 ms/frame of composite and frees a FlexRAM bank
(DTCM free for locals 12,896 → 45,664 B). The sections below describe the
earlier `a6e9d6426` capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 (both A and B) |
| Image | `profile` env: `-Os` base; selective-O3 path crosses pixel-feedback, filter-pipeline, screen-AA, and scan `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `a6e9d6426` (A: `54d0eb7fa`) |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture, `-D HS_PROFILE_EPOCH_REVS=3400` |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"` |

Image size (`profile` env, this effect only): `FLASH: code:114384,
data:190792, headers:9192` / `RAM1: variables:315168, code:50120,
padding:15416, free:143584`. Shipping `phantasm` image at the same tip:
RAM1 code 173,976 B, down from 196,568 B at `54d0eb7fa`.

Exactness cross-check: window frames 3921–3936 root counter cyc ÷ 600 MHz
matches the measured wall sum within **1.6 ppm**
(`tools/parse_profile.py ... validate`, VALID; A side 1.9 ppm).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): peak frame render
61.86 ms (frames 1201–1216), spilled **0/6688** frames. The A side peaked at
62.75 ms in the same window and spilled 1/6688.

A display window is 62.5 ms. Feedback reads outside the segment band, so the
pipeline renders the FULL canvas (41,472 px), not one quadrant. Every preset
holds 16 fps; the tightest, preset 6 (WavyTrails), peaks 0.64 ms under the
window. The `canvas_buffer_wait` scope is the round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: 12 style holds of 241 frames each, each followed by a
parameter transition into the next style; the solid is rebuilt per preset.

### Peak window: WavyTrails transition (window frames 1201–1216)

```
frame                    62.19 ms  37.31 Mcyc  100%
  mf_mesh_draw            7.45 ms   4.47 Mcyc   12%
    filter_blend          1.11 ms   0.66 Mcyc    2%  x10315  64 cyc/blend
  mf_feedback_flush      46.64 ms  27.98 Mcyc   75%
    feedback_composite   38.59 ms  23.16 Mcyc   62%
    feedback_populate     8.02 ms   4.81 Mcyc   13%
    feedback_litscan      0.2 us     139 cyc     0%
  mf_timeline_step       59.3 us   35.6 kcyc     0%
  mf_apply_params        11.7 us    7.0 kcyc     0%
  canvas_clear            0.41 ms   0.25 Mcyc    1%
  canvas_buffer_wait      7.62 ms   4.57 Mcyc   12%
```

Wall min/avg/max = 54.88/62.19/65.60 ms; render avg/max = 54.57/61.86 ms.
Against the A side's same window (render avg/max 55.50/62.75 ms), the mesh
draw fell from 8.68 to 7.45 ms per frame and `filter_blend` from 1.49 to
1.11 ms at an identical 165,041 blends. The flush rose from 46.34 to
46.64 ms: composite +0.13 ms, populate +0.17 ms. The net is 0.9 ms of render
per frame back, which clears the one spilled frame.

### Per-preset table

Each preset's flush cost is the mean over its clean holds; peak render and
spill come from `parse_profile.py ... buckets`, whose preset frames include
the transition that follows. The capture wraps from `Preset: 12/12` to
`1/12` at log line 6146 and every preset has 30–46 clean hold windows.

| rank | preset | style | flush ms A | flush ms B | peak render ms A | peak render ms B | spilled A | spilled B |
|---:|--:|---|--:|--:|--:|--:|--:|--:|
| 1 | 5 | SlowDust | 49.10 | 49.41 | 60.72 | 59.91 | 0/482 | 0/482 |
| 2 | 3 | EnergeticFire | 47.63 | 47.92 | 59.85 | 58.96 | 0/723 | 0/723 |
| 3 | 4 | Smoke | 47.62 | 47.90 | 55.97 | 55.86 | 0/664 | 0/664 |
| 4 | 9 | Miasma | 46.52 | 46.84 | 58.94 | 58.02 | 0/482 | 0/482 |
| 5 | 2 | SlowFire | 46.21 | 46.47 | 54.46 | 54.07 | 0/723 | 0/723 |
| 6 | 8 | MeltingLo | 45.46 | 45.73 | 55.95 | 54.78 | 0/482 | 0/482 |
| 7 | 6 | WavyTrails | 45.07 | 45.34 | 62.75 | 61.86 | 1/482 | 0/482 |
| 8 | 7 | MeltingHi | 44.01 | 44.18 | 51.73 | 51.36 | 0/482 | 0/482 |
| 9 | 10 | LooseWormhole | 43.08 | 43.27 | 52.78 | 51.92 | 0/482 | 0/482 |
| 10 | 1 | ArcingLightning | 42.26 | 42.49 | 55.71 | 54.94 | 0/722 | 0/722 |
| 11 | 11 | TightWormhole | 40.99 | 41.15 | 59.47 | 58.37 | 0/482 | 0/482 |
| 12 | 12 | WigglingWormhole | 39.26 | 39.41 | 61.47 | 59.97 | 0/482 | 0/482 |

The flush is 0.15–0.31 ms per frame slower on every preset (0.4–0.7%), and
peak render is 0.1–1.5 ms lower on every preset. All twelve are 16 fps.

### Per-pixel figures

The peak window blends 10,315 px/frame against the 41,472 px full canvas
(0.25× coverage) at 64 cyc/blend. The composite covers every canvas pixel:
23.16 Mcyc / 41,472 px ≈ 558 cyc/px.

## Column-ISR / DMA marshaling cost

```
isr_wake        1147/frame  min/avg/max 0.56/1.72/11.38 us  cpu 3.16%
isr_pack         143/frame  min/avg/max 6.36/7.22/9.55 us   cpu 1.66%
isr_dma_submit   143/frame  min/avg/max 0.83/0.94/1.01 us   cpu 0.21%
```

- Packing costs 7.7× the DMA submit per column.
- The HD107S wire transfer runs asynchronously on the DMA engine.
- ISR share totals 5.0% (≈ 3.1 ms per 62.5 ms window). The render
  figures already absorb it, so the budget is the whole window, and the
  worst preset now fits with no speedup required.
- ISR figures match the A side within 1%.

## Summary ranking

1. `feedback_composite` — 62% of the frame, 38.6 ms: the per-pixel
   hue-rotating two-pixel composite.
2. `feedback_populate` — 13%, 8.0 ms: coarse warp-field repopulation.
3. `mf_mesh_draw` — 12%, 7.45 ms: wireframe rasterization into the
   Orient + AntiAlias + Feedback pipeline.
4. `canvas_buffer_wait` — 12%, 7.6 ms of idle to the flip.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `mf_mesh_draw`; calls ≈ blended pixels.
- Selective-O3: the pixel-feedback region (`flush`, composite), the
  filter-pipeline sink, screen AA, and scan regions. The identity-fade and
  custom-color composites now run from flash inside flattened
  `HS_FLASH_MEMBER` lambdas; no autoplay preset reaches them.
- The flush slowdown is consistent with the hot composite now inlining its
  per-cell helper and the flush body shifting by 152 B; the mesh-draw gain is
  consistent with the smaller translation unit letting the -Os inliner keep
  more of the anti-aliased plot path inline. Neither attribution was isolated
  by a separate A/B.
- The epoch stretch keeps one effect instance across the full cycle; it does
  not change per-frame cost.
- The A side ran from the detached `hs-itcm-base` worktree at `54d0eb7fa`
  with no local changes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`,
`HS_PROFILE_WINDOW=16`; `just profile MeshFeedback` = build + flash + capture.
