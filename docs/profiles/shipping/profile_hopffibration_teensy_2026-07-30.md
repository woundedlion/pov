# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-07-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HopfFibration`).
Raw capture: `build/prof/hopffibration_ship.log` (23:47 re-capture at tip
`aea5d301`, COM3). Replaces `profile_hopffibration_teensy_2026-07-25.md`,
captured before the trail-gate fusion landed.

Re-capture note: an A/B against a worktree at `fc96824b` (parent of the
`Plot::rasterize` RasterOptions refactor, `ef753ee0`) shows the refactor is
perf-neutral on device — deterministic windows align frame-for-frame, every
window mean within 0.2% (peak frame render 57.62 → 57.74 ms, +0.12 ms),
0/1088 frames spilled on both sides. The tables below are from the 16:57
capture of the same day, and every peak on this page is that capture's: its
peak frame render is 58.06 ms against the 23:47 capture's 57.74 ms, a 0.6%
capture-to-capture spread.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board COM3 |
| Image | `profile` env: `-Os` base + newlib-nano + `HS_PROFILE_ENABLE`; the trail path crosses the `HS_O3` regions on `Plot::rasterize` and `Plot::gate_trail_edges` (`HS_O3_FN` on `render_trails`) |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, tip `d6147b06` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture, default 120 s epoch (no straddle windows) |
| Reproduce | `bash tools/profile_one.sh HopfFibration profile 70 32` |

Image size (`arm-none-eabi-size -A`): `.text.headers 5,120` /
`.text.code 15,420` / `.text.progmem 535,356` / `.text.csf 3,224` /
`.text.itcm 41,248` / `.data 2,752` / `.bss 311,808` / `.bss.dma 519,552`.

Exactness cross-check: window frames 193–224 root counter cyc ÷ 600 MHz
matches the measured wall sum within **0.6 ppm**
(`tools/parse_profile.py ... validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hf_render_trails`
worst window 49.64 ms/f (frames 193–224), peak frame render **58.06 ms**,
spilled **0/1088 frames (0%)**.

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
HopfFibration overrides neither `needs_full_frame()` nor `persists_pixels()`,
so it draws the quadrant clip only. Trail density swings widely as the fibers
tumble through the view: the heavy phase peaks at 58.06 ms render, still inside
one window, and the light phase runs at 13.6 ms. Every frame of the pass holds
16 fps with 4.4 ms of margin at the worst frame. The `canvas_buffer_wait` scope
is the round-up idle to the next display flip, by design — it absorbs the
difference, which is why it swells to 78% of the frame in the light phase.

## Phase-by-phase readout

Phase schedule: one continuous tumble; density varies with how much of the
fibration faces the viewer, so the pass alternates between heavy and light
regimes on a multi-second period rather than a fixed schedule.

### Heavy tumble (window frames 193–224, worst of the capture)

```
frame                   62.85 ms  37.71 Mcyc  100%
  hf_render_trails      49.13 ms  29.48 Mcyc   78%
    hf_trail_raster     38.90 ms  23.34 Mcyc   62%  x148.1  262.6 us/trail
      filter_blend       5.21 ms   3.13 Mcyc    8%  x49310  63.4 cyc/blend
    hf_trail_gate         9.31 ms   5.59 Mcyc   15%  x210     44.4 us/trail
  canvas_buffer_wait    13.42 ms   8.05 Mcyc   21%
  hf_project_record       190 us    114 kcyc    0%
  canvas_clear             88 us   52.8 kcyc    0%
  hf_timeline_step        14.9 us    8.9 kcyc    0%
  hf_advance_tumble        0.5 us     327 cyc    0%
```

Wall min/avg/max = 39.5/62.9/84.5 ms — the raster dominates at 62% of the
frame, and the gate costs 15% on top of it. The 210 gate calls per frame are
one per fiber, against 148 raster calls: 62 fibers per frame are gated away
whole before any rasterization, which is exactly the work the gate exists to
avoid. Everything outside the trail path is noise — projection, tumble
advance, timeline step and canvas clear together account for under 0.5% of the
frame.

### Light tumble (window frames 737–768)

```
frame                   62.62 ms  37.57 Mcyc  100%
  canvas_buffer_wait    49.02 ms  29.41 Mcyc   78%
  hf_render_trails      13.32 ms   7.99 Mcyc   21%
    hf_trail_raster      8.15 ms   4.89 Mcyc   13%  x55.1   147.9 us/trail
      filter_blend       1.07 ms    643 kcyc    2%  x10194  63.0 cyc/blend
    hf_trail_gate         4.26 ms   2.55 Mcyc    7%  x210     20.3 us/trail
  hf_project_record       190 us    114 kcyc    0%
```

Wall min/avg/max = 48.5/62.6/75.9 ms — render collapses to 13.6 ms and the
buffer wait takes the rest. The gate's call count is unchanged at 210/frame
(it is per-fiber, not per-visible-fiber) but its per-call cost halves to
20.3 µs, because a fiber that fails the clip early exits before walking its
polyline. Raster calls fall to 55/frame, so 155 of 210 fibers are gated away
here — the ratio that makes the gate pay for itself.

### A/B against the pre-fusion baseline

Both sides captured on the same board (COM3), same window, `CCACHE_DISABLE=1`
and a wiped `.pio/build/profile` on each side, so no stale object can survive
between them.

| Window | Scope | Baseline | Fused | Δ |
|---|---|--:|--:|--:|
| 193–224 | peak frame render | 58.61 ms | 58.06 ms | −0.55 ms (−0.9%) |
| 193–224 | `hf_render_trails` | 953.0 Mcyc | 943.4 Mcyc | −9.7 Mcyc (−1.0%) |
| 193–224 | `hf_trail_gate` | 194.3 Mcyc | 178.8 Mcyc | −15.5 Mcyc (−8.0%) |
| 193–224 | `hf_trail_raster` | 746.88 Mcyc | 746.91 Mcyc | +0.03 Mcyc (+0.004%) |
| 737–768 | `hf_trail_gate` | 98.5 Mcyc | 81.7 Mcyc | −16.8 Mcyc (−17.1%) |
| 737–768 | peak frame render | 23.89 ms | 23.60 ms | −0.29 ms (−1.2%) |

`filter_blend` calls are **identical on both sides** — 1,577,911 in the heavy
window and 326,201 in the light one — as are the `hf_trail_raster` call counts
(4,740 and 1,764). The fused gate therefore produces bit-identical visibility
bits and the rasterizer does bit-identical work; the saving is purely the
second staging pass that no longer happens.

The win does not move cadence, because the frame was already inside one
display window: the saved render time reappears in `canvas_buffer_wait`
(248.1 → 257.7 Mcyc in the heavy window), leaving total frame cycles flat at
1.2068 Gcyc. It buys headroom, not fps.

### Per-pixel figures

Heavy phase: 49,310 blended px/frame against the 10,368 px quadrant = **4.76×
coverage** (anti-aliased polylines blend a pixel once per crossing edge), at
63.4 cyc/blend and 473 raster cycles per blended pixel. Light phase: 10,194
blended px/frame = 0.98× coverage, 63.0 cyc/blend, 480 raster cycles per
blended pixel. Cycles per blended pixel are flat across a 4.8× density swing,
so the raster cost is linear in blended pixels with no density-dependent term.

## Column-ISR / DMA marshaling cost

```
isr_wake        x1159/frame  min/avg/max 0.586/1.908/11.305 us  cpu 3.51%
isr_pack         x144.9/frame  min/avg/max 6.178/6.968/9.178 us  cpu 1.60%
isr_dma_submit   x144.9/frame  min/avg/max 0.653/0.931/1.206 us  cpu 0.21%
```

- Submit is 7.5× cheaper than pack (0.93 µs against 6.97 µs): the eDMA
  descriptor write is trivial next to marshaling a column's pixels into the
  wire buffer.
- Pack and submit run at the same 144.9/frame rate — one pair per column — while
  wake fires 8× more often; wake is the flywheel's per-edge tick, not per column.
- ISR CPU share is **5.32%**, so a 62.5 ms display window leaves ≈59.2 ms of
  render budget. The worst frame uses 58.06 ms of it, a 1.9% margin. No
  speedup is required to hold 16 fps, but there is no room for another heavy
  phase either.

## Summary ranking

1. `hf_trail_raster` — 62% of the frame, 38.90 ms: anti-aliased polyline
   rasterization of 148 visible fibers, 4.76× quadrant overdraw. This is the
   effect.
2. `hf_trail_gate` — 15% of the frame, 9.31 ms: per-fiber clip test over all
   210 fibers. Pays for itself by removing 62 fibers from the raster in the
   heavy phase and 155 in the light one.
3. `filter_blend` — 8% of the frame, 5.21 ms (nested under raster): 49,310
   blended pixels at 63.4 cyc each.
4. Everything else — under 0.5% combined.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs); the 5.32% ISR share above is
  included in every figure on this page.
- `filter_blend` parents under `hf_trail_raster`; its subtree is hidden in
  windows where that parent had 0 calls. Its calls ≈ blended pixels.
- Selective -O3: `render_trails` carries `HS_O3_FN`, so both the raster and the
  gate run under the `HS_O3` regions on the `-Os` base image. The regions
  activate only on the device build, so this capture measures the shipping
  selective-O3 configuration by construction.
- No dwell-compression knobs were used — HopfFibration is not a cycling effect,
  so the capture is a single continuous regime and the 70 s window never
  crosses the 120 s epoch.
- Both captures ran on committed trees with clean working copies; the fused
  side was `606e0d85` in a worktree, landed as `d6147b06`.
- The parallel first attempt at this A/B failed: two `profile_one.sh` runs on
  two different boards collided inside the Teensy loader, which is a single
  shared instance. Per-board device locks do not make flashing parallel-safe —
  serialize the flashes.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HopfFibration`,
`HS_PROFILE_WINDOW=32`; `just profile HopfFibration` = build + flash + capture.
