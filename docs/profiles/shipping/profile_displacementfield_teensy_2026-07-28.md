# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-07-28, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`. Replaces
`profile_displacementfield_teensy_2026-07-25.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board COM4 |
| Image | `profile` env: `-Os` + newlib-nano + DMA LEDs + `HS_PROFILE_ENABLE`; `Scan::DistortedRingStack::draw` (R2 fused driver), the ring/polyline `distance()` chain in `sdf.h`, and the effect's own `draw_rings` run inside `HS_O3` regions (docs/selective_o3_spec.md) |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `542a5b49` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `just profile DisplacementField` |

Image size: `FLASH: code:67740, data:539204, headers:8472` / `RAM1:
variables:314592, code:44024, padding:21512, free:144160` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 897–928 root counter cyc ÷ 600 MHz
(1,206,802,482 ÷ 600 = 2,011,337 µs) matches the measured wall sum
(2,011,344 µs) within **3.2 ppm** (`tools/parse_profile.py ... validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `df_timeline_step`
19.2–55.5 ms/f over the dwell cycle, peak frame render **58.71 ms**
(frames 673–704), spilled **0/1088 frames (0%)**.

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px. The
whole pass holds 16 fps, but this is the tightest green on the roster: the peak
frame leaves **3.8 ms** against the 62.5 ms window, and only 0.7 ms once the
5.00% ISR share is subtracted from the budget. The `canvas_buffer_wait` scope is
the round-up idle to the next display flip, by design — it collapses from
43.6 ms at the dwell trough to 7.0 ms at the peak, which is the margin being
consumed.

## Phase-by-phase readout

Phase schedule: one continuous regime whose cost breathes with the displacement
field's dwell cycle — the ring population entering the stack scan rises and
falls, taking the per-ring LUT bakes with it. The capture covers a full breath:
cost ramps 27 → 55 ms/f over windows 1–416, holds near the peak through window
800, then falls to a 19.2 ms/f trough at frames 897–928 before rebuilding.

### Dwell peak — worst of the capture (window frames 673–704)

```
frame                    62.58 ms  37.55 Mcyc  100%
  df_timeline_step       55.48 ms  33.29 Mcyc   89%
    df_draw_rings        55.36 ms  33.22 Mcyc   88%
      df_hue_table_prep   0.88 ms   0.53 Mcyc    1%  x22.8  38.7 us/call
      df_lut_bake        10.01 ms   6.01 Mcyc   16%  x45.4  220 us/bake
      df_chunk_cull       1.29 ms   0.78 Mcyc    2%  x47.2  27.4 us/call
      df_fused_scan      42.07 ms  25.24 Mcyc   67%  x1
        filter_blend      0.89 ms   0.53 Mcyc    1%  x10397 51.2 cyc/blend
  canvas_clear           89.47 us  53685 cyc     0%
  canvas_buffer_wait      7.01 ms   4.21 Mcyc   11%
  df_prepare_fields       0.09 us     67 cyc     0%
```

Wall min/avg/max = 59.22/62.58/64.94 ms — a narrow spread, because this effect's
render evolves smoothly frame to frame and wall tracks the render *derivative*
(see the wall-metric caveat below), not the display, which is rigid.
`df_fused_scan` at 42.07 ms is 67% of the frame and 76% of the ring draw: the
fused ring-stack rasterizer is the effect, and everything else is preparation
for it. `df_lut_bake` is the clear second at 10.01 ms over 45.4 bakes/frame
(220 µs each) — the bake count, not the bake cost, is what the dwell cycle
moves. `filter_blend` is a thin 1%: 10,397 blends is 1.00× the quadrant, so
there is essentially no overdraw and the scan's probe work dominates.

### Dwell trough (window frames 897–928)

```
frame                    62.85 ms  37.71 Mcyc  100%
  df_timeline_step       19.20 ms  11.52 Mcyc   31%
    df_draw_rings        19.07 ms  11.44 Mcyc   30%
      df_hue_table_prep   0.65 ms   0.39 Mcyc    1%  x11.2  58.5 us/call
      df_lut_bake         1.47 ms   0.88 Mcyc    2%  x13.1  112 us/bake
      df_chunk_cull       0.36 ms   0.21 Mcyc    1%  x13.6  26.1 us/call
      df_fused_scan      15.87 ms   9.52 Mcyc   25%  x1
        filter_blend      0.78 ms   0.47 Mcyc    1%  x9153  51.4 cyc/blend
  canvas_clear           85.47 us  51296 cyc     0%
  canvas_buffer_wait     43.57 ms  26.14 Mcyc   69%
  df_prepare_fields       1.88 us   1130 cyc     0%
```

Wall min/avg/max = 58.83/62.85/67.27 ms — indistinguishable from the peak
window's wall, as it must be: every frame occupies one display window either
way, and only the idle split moves. Render is 2.9× cheaper than at the peak, and the cut is
almost entirely structural: `df_lut_bake` drops 45.4 → 13.1 calls/frame
(10.01 → 1.47 ms) and `df_chunk_cull` 47.2 → 13.6, while `filter_blend`'s
per-blend cost is flat at 51.4 cyc and its count barely moves (10,397 → 9,153).
Fewer rings reach the stack scan; the ones that do cost the same.

### Per-pixel figures

Dwell peak: 10,397 blended px/frame against the 10,368 px quadrant = **1.00×
coverage**, at **51.2 cyc/blend**. `df_fused_scan` spends 25.24 Mcyc over those
blends = **2,428 cyc per blended pixel**, so ~97% of the scan is probe
evaluation that never reaches a blend. The trough holds the same ratio
(9,153 blends at 51.4 cyc/blend), which is why the dwell cycle moves render
cost without moving coverage: the effect is probe-bound, and the field's ring
population sets the probe count.

## Column-ISR / DMA marshaling cost

```
isr_wake        x1154/frame  min/avg/max 0.69/1.78/10.78 us  cpu 3.27%
isr_pack        x144/frame   min/avg/max 6.11/6.60/8.77 us   cpu 1.52%
isr_dma_submit  x144/frame   min/avg/max 0.66/0.93/1.17 us   cpu 0.21%
```

(Window frames 673–704, the peak.)

- Submit is 14% of pack's CPU cost — the column marshal dominates, the DMA kick
  is nearly free.
- The wire transfer is asynchronous; submit returns in 0.93 µs and the engine
  runs against the render.
- Total ISR share is 5.00%, leaving ≈ 59.4 ms of render budget per 62.5 ms
  window. The peak frame uses 58.71 ms, so it clears by **0.7 ms** — no speedup
  is required today, but this is the roster's thinnest margin and a workload
  increase of ~1% would put it over.

## Summary ranking

1. `df_fused_scan` — 67% of the peak frame, 42.07 ms: the fused ring-stack
   rasterizer, 2,428 cyc per blended pixel. Already inside an `HS_O3` region.
2. `df_lut_bake` — 16%, 10.01 ms: 45.4 per-ring displacement-LUT bakes at
   220 µs each; the count scales with the dwell cycle, so this is the scope
   that decides whether the peak is 58.7 ms or something worse.
3. `df_chunk_cull` — 2%, 1.29 ms: 47.2 chunk tests/frame at 27.4 µs, the
   `field_bound` cull that keeps the scan off empty space.
4. `df_hue_table_prep` — 1%, 0.88 ms.
5. `df_prepare_fields` — 0.09 µs/frame; the field stack setup is free.

Against the ledger (docs/selective_o3_spec.md): the fused-tip `-O3` ceiling for
`df_fused_scan` was 46.0 ms against 59.1 ms at `-Os`, and the selective build
now measures 42.07 ms at the dwell peak — the region is carrying its full
crossing.

## Caveats

- **Wall min/avg/max is not a cadence signal.** `buffer_wait` sits at the head of
  `draw_frame`, so a measured frame runs from the end of the previous render to
  the end of this one and wall picks up the render *derivative*:
  `wall(N) = P + render(N) - render(N-1)`, which reproduces every frame of this
  capture to a 7 us median / 0.22 ms max residual (P = 62.48 ms). The display
  flip itself is rigid and `buffer_wait` absorbs the remainder exactly; total
  pass wall is 0.999x `frames x 62.5 ms`, i.e. one display window per frame with
  no doubled frames. A volatile per-frame cost therefore shows a wide wall
  spread while sitting perfectly on cadence — read `spilled` and peak render,
  never wall max.
- All scopes absorb ISR time (CYCCNT free-runs) — the 5.00% ISR share is inside
  every number above.
- `filter_blend` parents under `df_fused_scan`; its subtree is hidden in windows
  where that parent had 0 calls. Its calls ≈ blended pixels regardless.
- `Scan::DistortedRingStack::draw` (R2), the `sdf.h` ring/polyline `distance()`
  chain, and `draw_rings` sit inside `HS_O3` regions, which activate only on the
  `-Os` device build — the `profile` env therefore measures the shipping
  selective-O3 configuration.
- No dwell-compression knobs: the capture runs the effect's own cycle at
  shipping speed, and 70 s covers a full breath inside one 120 s epoch.
- Captured at tip `542a5b49`. The working tree carried uncommitted KiCad/doc
  edits under `hardware/phantasm/` and `docs/` only; nothing on the firmware
  build path was dirty.
- **Blend cost fell sharply since 2026-07-25** (`32478115`): `filter_blend` was
  169 cyc/blend over 12,120 blends (3.42 ms/f) and now reads 51.2 cyc/blend over
  10,397 (0.89 ms/f), taking peak render 60.59 → 58.71 ms. The core/color batch
  in that range (findings 363–374, notably the `Pixel16` → `CRGB` split-decode)
  is the obvious candidate; this capture establishes the new number but does not
  isolate the cause — an A/B would be needed to attribute it.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`,
`HS_PROFILE_WINDOW=32`; `just profile DisplacementField` = build + flash +
capture.
