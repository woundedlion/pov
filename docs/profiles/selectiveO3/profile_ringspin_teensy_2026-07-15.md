# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-07-15, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`). Raw
capture: `build/profile_ringspin_selo3.log`; same-tip baselines:
`build/profile_ringspin_base_os.log` (pre-region `-Os`),
`build/profile_ringspin_base_o3.log` (`profile_o3` global-O3 twin). This is
the first report in `selectiveO3/`; the `Os/` and `O3/` 2026-07-15 reports
remain the per-level references for the pre-region tree.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: shipping Phantasm flags (`-Os`, newlib-nano, DMA LEDs, N=4) + `HS_PROFILE_ENABLE` + the landed `HS_O3` regions (blend sink, `Scan::RingGroup`, sdf.h `HS_O3_FN`s) |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist |
| Method | `HS_PROFILE` cycle scopes, window = 64 frames, 120 s capture |
| Reproduce | `just profile RingSpin` (the `-Os` `profile` env measures the shipping selective-O3 config) |

Image size: `FLASH: code:44052, data:77104, headers:8884` / `RAM1:
variables:345696, code:30168, padding:2600, free:145824` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 1729–1792 root counter
2,394,631,381 cyc ÷ 600 MHz = 3,991,052 µs vs measured wall sum
3,991,068 µs — **3.9 ppm**.

## Frame cadence

A display window (half-revolution at 480 RPM) is 62.5 ms; `draw_frame` wall
time snaps up to whole windows. The effect renders one quadrant ≈ 10,368 px.
**Every one of the 29 captured windows averages 62.1–62.8 ms — 16 fps
hard-locked across the whole trail-density cycle**, with `rs_buffer_wait`
(the round-up idle, by design) at 33–43% of the frame. At `-Os` on the same
tree the heavy windows averaged 71–98 ms (8↔16 jitter); at the global-O3
ceiling the cadence is the same 16 fps lock seen here.

## Phase-by-phase readout

RingSpin sweeps trail density continuously; the two regimes below bracket it.

### Steady regime (window frames 1729–1792, median scan)

```
frame                   62.36 ms  37.42 Mcyc  100%
  rs_draw_rings         35.63 ms  21.38 Mcyc   57%
    rs_ring_scan        34.89 ms  20.94 Mcyc   56%  x76    459 us/ring
      filter_blend       5.43 ms   3.26 Mcyc    9%  x38642  84 cyc/blend
  rs_timeline_step      66.1 us     40 kcyc     0%
  rs_buffer_wait        26.67 ms  16.00 Mcyc   43%
```

Wall min/avg/max = 47.6/62.4/76.9 ms. The fused group scan is 98% of the
draw; the region's -O3 codegen puts the median scan within 4% of the
global-O3 ceiling (33.3 ms) from an `-Os` baseline of 46.3 ms.

### Heavy regime (window frames 1537–1600, max scan)

```
frame                   62.59 ms  37.55 Mcyc  100%
  rs_draw_rings         41.78 ms  25.07 Mcyc   67%
    rs_ring_scan        41.04 ms  24.62 Mcyc   66%  x76    540 us/ring
      filter_blend       6.52 ms   3.91 Mcyc   10%  x46478  84 cyc/blend
  rs_timeline_step      69.0 us     41 kcyc     0%
  rs_buffer_wait        20.73 ms  12.44 Mcyc   33%
```

Wall min/avg/max = 41.1/62.6/80.8 ms. Even the densest trail window keeps
20.7 ms of idle — the max-window figure at `-Os` was 57.5 ms of scan
(zero-idle, 2-window frames); the ceiling's max is 39.2 ms. Individual
frames still range to 80.8 ms wall (a frame that misses one window boundary
picks up the next), but the window *average* stays locked.

### Per-pixel figures

Steady: 38,642 blended px/frame vs the 10,368 px quadrant = 3.7× overdraw
(the per-blend-chain bound the fused/slim RingGroup work left in place);
84 cyc/blend at the -O3 sink (vs 149–156 at `-Os`, 103 measured at global
-O3); total scan cost ≈ 542 cyc per blended pixel.

## Column-ISR / DMA marshaling cost

```
isr_wake        1149/frame  min/avg/max 0.58/1.66/147 us   cpu 3.05%
isr_pack         144/frame  min/avg/max 4.3/5.7/146 us     cpu 1.30%
isr_dma_submit   144/frame  min/avg/max 0.64/0.96/1.2 us   cpu 0.22%
```

- Submit is negligible next to pack; the wire transfer itself is async DMA.
- ISR CPU share ≈ 4.6% ⇒ ~59.6 ms of render budget per 62.5 ms window;
  the steady render (35.7 ms) and heavy render (41.9 ms) both fit with
  margin — which is exactly the 16 fps lock observed.

## Summary ranking

1. `rs_ring_scan` — 56–66% of frame, 34.9–41.0 ms: the fused ring-group
   scan, now -O3; within 4–5% of the global-O3 ceiling.
2. `rs_buffer_wait` — 33–43%: display-sync idle (headroom, not cost).
3. `filter_blend` — 9–10% (inside the scan): 84 cyc/blend at the -O3 sink.

Selective -O3 closes ~90% of the -Os → -O3 gap on this effect's hot scope
and delivers the ceiling's cadence outright.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs); that is the point of
  profiling under the live driver.
- `filter_blend` parents under whichever scope first enters it; its
  calls/cycles ≈ blended pixels.
- Per-call figures are window totals ÷ calls; epoch-boundary windows were
  skipped when picking representatives.
- Shipping config is `-Os` + selective `HS_O3` regions — this report *is*
  the shipping image's profile, unlike the `O3/` reports (unshippable
  single-effect global-O3 images).
- Captured on a clean landed tip (`0f7b0616`); no uncommitted working-tree
  state.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingSpin`,
`HS_PROFILE_WINDOW=64`; `just profile RingSpin` = build + flash + capture.
