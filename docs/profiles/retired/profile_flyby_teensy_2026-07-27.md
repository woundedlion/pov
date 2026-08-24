# Flyby on-device profile — Teensy 4.0, segmented mode (2026-07-27, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Flyby`).
Raw capture: `build/prof/flyby_ship.log`. Replaces
`profile_flyby_teensy_2026-07-25.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, board COM4 |
| Image | `profile` env: `-Os` + newlib-nano + DMA LEDs + `HS_PROFILE_ENABLE`; the stereographic shader runs inside an `HS_O3` region |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Flyby 288×144, single-entry playlist, tip `ece0955b` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 310 s capture, `HS_PROFILE_EPOCH_REVS=2560` |
| Reproduce | `just profile Flyby` |

Image size: `FLASH: code:46148, data:1062396, headers:8632` / `RAM1:
variables:314368, code:27880, padding:4888, free:177152` / `RAM2:
variables:519552, free:4736`.

Exactness cross-check: window frames 1457–1472 root counter cyc ÷ 600 MHz
matches the measured wall sum within **3.4 ppm**
(`tools/parse_profile.py ... validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fly_shader_draw`
37.3–44.0 ms/f across the 5 presets, worst preset 2 at 44.00 ms/f, peak frame
render 45.05 ms (frames 2961–2976), spilled 0/4928 frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px. Every
preset holds 16 fps with 17 ms of margin at the worst frame — the widest
headroom of any effect in the upper half of the roster. The
`canvas_buffer_wait` scope is the round-up idle to the next display flip, by
design; it tracks the shader cost inversely, 18.4 ms at preset 2 against
27.3 ms at preset 0.

## Phase-by-phase readout

Phase schedule: five presets chained by a continuous LERP — Flyby never holds a
preset still, it interpolates between them, so there is no steady-state hold
window. Per the Method note below, each preset's figure is the peak window
nearest its anchor. 11 `Preset:` markers over the pass, all 5 presets visited,
cycle wrapped.

### Preset 2 anchor — worst of the capture (window frames 2961–2976)

```
frame                   62.50 ms  37.50 Mcyc  100%
  fly_shader_draw       44.00 ms  26.40 Mcyc   70%
  fly_timeline_step     15.81 us   9491 cyc     0%
  canvas_clear          91.88 us  55145 cyc     0%
  canvas_buffer_wait    18.38 ms  11.03 Mcyc   29%
```

Wall min/avg/max = 60.67/62.50/63.92 ms — a 3.25 ms spread across 16 frames,
the flattest in the capture. Flyby is a single flat shader scope with no
sub-tree: there is no cull, no per-primitive setup, no blend filter. Every
pixel of the quadrant is evaluated every frame, so the cost is set by shader
complexity at the current LERP point and nothing else. `fly_timeline_step` is
15.8 µs — the preset interpolation that drives the whole effect is 0.03% of
the frame.

### Preset 0 anchor — lightest of the capture (window frames 3825–3840)

```
frame                   62.45 ms  37.47 Mcyc  100%
  fly_shader_draw       34.99 ms  21.00 Mcyc   56%
  fly_timeline_step     20.13 us  12106 cyc     0%
  canvas_clear          91.44 us  54870 cyc     0%
  canvas_buffer_wait    27.34 ms  16.40 Mcyc   44%
```

Wall min/avg/max = 61.87/62.45/62.97 ms — 9.0 ms/frame cheaper than preset 2
at an identical pixel count, so the delta is purely per-pixel shader work along
the LERP path. Render swings 35.0→44.0 ms across the cycle while wall stays
pinned at one display window throughout: the cadence never changes, only how
much of the window is spent.

### Per-preset table

**Method note:** Flyby has no clean-hold regime — the presets are chained by a
continuous LERP, so the figure for each is the peak window nearest its anchor
rather than a steady hold. Cycle wrap to preset 0 confirmed (5 distinct
presets, peak index 4). Windows per preset: 60–68.

| # | holds | fly_shader_draw ms | peak render ms | spilled | fps |
|---|--:|--:|--:|--:|--:|
| 2 | 60 | 44.00 | 45.05 | 0/960 | 16 |
| 3 | 60 | 43.51 | 45.03 | 0/960 | 16 |
| 1 | 68 | 43.15 | 44.15 | 0/1088 | 16 |
| 4 | 60 | 42.12 | 43.93 | 0/960 | 16 |
| 0 | 60 | 37.29 | 38.20 | 0/960 | 16 |

The span is only 1.18× (37.3–44.0 ms) — the tightest per-preset spread of any
cycler on the roster. Four of the five presets sit within 2 ms of each other;
preset 0 is the single outlier.

### Per-pixel figures

`filter_blend` does not appear: the shader writes every pixel of the quadrant
directly, with no blended-pixel path. Coverage is therefore exactly 1.0× —
10,368 px/frame — and the per-pixel cost is 44.00 ms ÷ 10,368 = **4.24 µs/px
(2,546 cyc/px)** at preset 2, against 3.37 µs/px (2,025 cyc/px) at preset 0.
That per-pixel figure is what the `HS_O3` region on the shader path buys down;
it is the whole cost of the effect.

## Column-ISR / DMA marshaling cost

```
isr_wake        x1153/frame  min/avg/max 0.62/1.79/11.51 us  cpu 3.30%
isr_pack        x144/frame   min/avg/max 6.35/7.01/9.50 us   cpu 1.61%
isr_dma_submit  x144/frame   min/avg/max 0.78/0.93/1.06 us   cpu 0.21%
```

- Submit is 13% of pack's CPU cost — marshaling dominates, the DMA kick is
  nearly free.
- The wire transfer is asynchronous; submit returns in 0.93 µs and the engine
  runs against the render.
- Total ISR share is 5.12%, leaving ≈ 59.3 ms of render budget per 62.5 ms
  window. The worst frame uses 45.05 ms, so Flyby holds 16 fps with 14.2 ms of
  budget unspent and needs no speedup.

## Summary ranking

1. `fly_shader_draw` — 70% of the frame, 44.00 ms: the entire effect. One flat
   scope, 2,546 cyc per quadrant pixel.
2. `canvas_clear` — 0.15%, 91.9 µs: the clip-scoped buffer clear.
3. `fly_timeline_step` — 0.03%, 15.8 µs: preset LERP.

There is no second hot spot to rank — Flyby is a pure per-pixel shader, so the
only lever on its cost is cycles per pixel.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs) — the 5.12% ISR share is inside
  every number above.
- `filter_blend` does not parent under any Flyby scope in this capture (the
  shader path does not enter it), so no blend subtree appears.
- The stereographic shader sits inside an `HS_O3` region, which activates only
  on the `-Os` device build — the `profile` env therefore measures the shipping
  selective-O3 configuration.
- `HS_PROFILE_EPOCH_REVS=2560` stretches the epoch so one 2,400-frame preset
  cycle fits a single effect instance. It changes epoch length, never per-frame
  render cost.
- Captured at tip `ece0955b` with a clean working tree.
- **Metric boundary changed at `a1b06170`**: `*_buffer_wait` now covers the
  `buffer_free()` spin only, where it previously covered the whole `Canvas`
  constructor. Render therefore now includes `canvas_clear` (91.9 µs/frame) and
  the buffer advance, which the 2026-07-25 report charged to the wait. Peak
  render is not directly comparable across that boundary.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Flyby`,
`HS_PROFILE_WINDOW=16`; `just profile Flyby` = build + flash + capture.
