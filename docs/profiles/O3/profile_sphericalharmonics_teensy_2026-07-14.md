# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode, **-O3** (2026-07-14)

Point-in-time snapshot (regenerate with `just profile`). Raw capture:
`build/profile_sphericalharmonics_o3.log` (64-frame windows, ~75 s single pass).

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 (IMXRT1062 @ 600 MHz), bench-attached over USB, no motor/LEDs |
| Image | `profile_o3` PlatformIO env = Phantasm shipping flags (newlib-nano, `USE_DMA_LEDS`, `tools/phantasm.ld`) + `-D HS_PROFILE_ENABLE`, but with base **`-O3`** (`-ffast-math -fno-finite-math-only`) instead of `-Os` — the -O3 twin of the shipping `-Os` profile |
| Driver | `POVSegmented<288, 4, 480>` — shipping 4-segment config, this board strapped as **segment 0 (master)**, flywheel + DMA LED ISRs live |
| Effect | `SphericalHarmonics<288, 144>` only (single-entry playlist), current working-tree state |
| Method | DWT `CYCCNT` cycle counters (`HS_PROFILE` RAII scopes, 600 cyc/µs) + `micros()` wall clock per `draw_frame`, dumped every 64 frames then reset |
| Reproduce | `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (no `just` recipe for the -O3 env) |

Image size: FLASH code **47,516 B**; ITCM (RAM1 code) **34,504 B**; RAM2 free 4,736 B.
-O3 inflates both code sections vs the shipping -Os image (33,508 / 21,400 B):
FLASH +14,008 B (+42 %), ITCM +13,104 B (+61 %) — see the `## -O3 vs -Os` section.

**Exactness cross-check** — the cycle counter and wall clock agree to sub-ppm:
window frames 513–576 root counter = 2,400,398,084 cyc = 4,000,663.47 µs vs
measured `micros()` window sum 4,000,670 µs (Δ ≈ 1.6 ppm).

## Frame cadence

At 480 RPM a display window (half-revolution) is **62.5 ms**; the driver
releases one `draw_frame` per window, each rendering this segment's quadrant
(72-row band × the window's 144-column half ≈ **10,368 px**). SphericalHarmonics
is a full-quadrant mesh/field rasterizer: a single `Scan::rasterize` sweeps the
whole quadrant and evaluates the harmonic field per pixel in the shader.
`filter_blend` confirms it shades essentially the **entire** quadrant every
frame (a constant 10,658 blends/frame, identical to the -Os run). Wall time snaps
up to a whole number of 62.5 ms windows.

At -O3 the render (rasterize + timeline + the un-scoped per-frame field build)
costs only **~14.5 ms/frame**, deeper inside one window than the -Os build's
~20 ms, so the cadence is a **rock-steady 16 fps (62.5 ms/frame)** across the
entire pass — the render cost breathes (12.2–17.1 ms as the field animates) but
comes nowhere near the window boundary, so it never flips to 8 fps.

`sh_buffer_wait` (the `buffer_free()` spin inside the `Canvas` ctor waiting for
the display flip) is timed separately and is exactly the round-up idle —
~48 ms/frame, the bulk of the window the (now cheaper) render leaves empty.

## Phase-by-phase readout (64-frame windows, per-frame averages)

SphericalHarmonics has **no discrete phases**: every steady window looks alike
(`sh_rasterize` 19–27 % of frame, `sh_buffer_wait` 73–80 %). One representative
window suffices; several others (frames 65–128, 129–192, 193–256, 641–704,
705–768) were checked and agree. The first window (`frames 1-64`, min ≈ 13.6 ms
epoch-start) is skipped. Columns: time/frame, cycles/frame, % of frame; leaves
add xN calls/frame and per-call cost.

### Steady render — 16 fps (representative: frames 513–576)

```
frame                 62.51 ms  37.51 Mcyc  100%
  sh_rasterize        14.46 ms   8.68 Mcyc   23%
    filter_blend       1.55 ms   0.93 Mcyc        x10658   87 cyc/blend
  sh_buffer_wait      48.04 ms  28.82 Mcyc   77%
  sh_timeline_step      12.2 us    7.3 kcyc    0%
```

Wall: min 60.4 / avg 62.5 / max 65.6 ms. The single `sh_rasterize` leaf (the
full-sphere `Scan::rasterize`, per-pixel harmonic-field eval in the shader) is
the whole render cost; `sh_timeline_step` only advances animations and is
negligible. The per-frame field build (`decode_lm` ×2 + `HarmonicField` ctor) is
cheap un-scoped setup — the un-instrumented remainder is **~0.7 µs/frame**.
This regime holds across the whole pass: `sh_rasterize` swings 12.2–17.1 ms
(19–27 % of frame) as the harmonic field animates, but coverage is fixed, so
the cadence never leaves 16 fps.

### Per-pixel figures

`filter_blend` (the pre-existing per-pixel counter in filter.h): a **constant**
10,658 blended pixels/frame in every window — the harmonic field covers
essentially the **entire 10,368-px quadrant** every frame (the swing in
`sh_rasterize` is per-pixel cost, not coverage) — at **87 cyc (0.15 µs) per
blend** (includes its own scope overhead). `sh_rasterize` ≈ 8.68 Mcyc over
10,658 blended px ⇒ **~814 cyc per blended pixel**, of which the blend is only
~87 cyc; the remaining **~727 cyc/px is the spherical-harmonic field
evaluation** in the shader (the dominant cost). Across windows the field eval
breathes from ~600 (frames 705–768) to ~872 cyc/blended px (frames 193–256).

## Column-ISR / DMA marshaling cost (`build/profile_sphericalharmonics_o3.log`)

`HS_ISR_PROFILE` raw-CYCCNT accumulators, representative steady window
(frames 513–576, 4.0 s). The pack and submit run inside the flywheel wake, on
the 1-in-8 wakes that render a column. Columns: rate, per-call min/avg/max,
CPU share:

```
isr_wake         18419/s  0.50 / 1.59 / ~228 us  cpu 2.92%
  isr_pack        2304/s  4.65 / 6.54 / ~227 us  cpu 1.50%
  isr_dma_submit  2304/s  0.60 / 0.91 / 1.23 us  cpu 0.20%
```

`isr_wake` = whole flywheel ISR (18,432/s is the COLUMN_US/8 wake grid);
`isr_pack` = the 72× `packPixel` marshal, once per column; `isr_dma_submit` =
`submitFrame` (overrun check + dcache flush + eDMA kick). The ISR machinery is
plain C with no fast-math surface, so it is unchanged from the -Os run (~2.9 %).

- **The DMA kick itself is ~1 µs of CPU** (submit floor 0.60 µs); the pixel-pack
  marshal around it is the real per-column cost (floor 4.65 µs). Whole
  per-column ISR CPU ≈ **13 µs of the 434 µs column period (~2.9 %)**. The wire
  transfer runs asynchronously on eDMA and costs no CPU.
- The avg/max figures inflate (pack avg 6.54 µs, max ~227 µs) — this is
  preemption by the serial dumps themselves, not real ISR growth; the `min`
  floor (pack 4.65 µs / 2789 cyc, submit ~0.6 µs / 360 cyc) is stable across
  every window (this window prints a 64-cyc submit min, a lone dump-boundary
  fluke).
- Net: the ISR machinery steals **~2.9 % of the chip**. Against the 62.5 ms
  window that's ~1.8 ms, leaving **~60 ms of render budget per window** — the
  ~14.5 ms render now fits ~4× over, hence the unwavering 16 fps.

## Summary ranking (steady window, share of the 62.5 ms frame)

1. display-window sync (`sh_buffer_wait`) — **77%** (48.0 ms), **idle by
   design**: the cheap render leaves most of the window empty, and wall time
   rounds up to the whole 62.5 ms window. -O3 makes render cheaper, so this idle
   grows to absorb the saved time.
2. `Scan::rasterize` (full-quadrant harmonic-field shade + blend) — **23%**
   (14.5 ms) — the **entire algorithmic render cost**; `filter_blend` (~1.55 ms)
   is just the composite, ~727 cyc/px is the field eval.
3. everything else (`sh_timeline_step` + un-scoped field build) — <1% (~13 µs).

SphericalHarmonics is a cheap, fully render-bound effect that never stresses the
window: the only real lever is the per-pixel harmonic-field math (`sh_rasterize`
at ~814 cyc/blended px), and even halving it would not change the 16 fps
cadence — it is already idle-bound, not compute-bound.

## -O3 vs -Os

Same effect, same driver, same 10,658-blend coverage — the **only** variable is
the optimization level (`profile_o3` vs the shipping `profile`/-Os env). All
figures are the representative steady window (frames 513–576) in each run.

| metric | -Os | -O3 | change |
|---|---|---|---|
| `sh_rasterize` (render) | 20.14 ms / 12.08 Mcyc | 14.46 ms / 8.68 Mcyc | **1.39× faster, −5.68 ms/frame** |
| `sh_rasterize` % of frame | 32 % | 23 % | −9 pts |
| cyc / blended px (total) | ~1,134 | ~814 | 1.39× |
| field eval cyc / px | ~964 | ~727 | 1.33× |
| `filter_blend` per blend | 170 cyc | 87 cyc | **1.95×** |
| frame wall / cadence | 62.5 ms / 16 fps | 62.5 ms / 16 fps | unchanged |
| `sh_buffer_wait` (idle) | 42.35 ms / 68 % | 48.04 ms / 77 % | +5.7 ms idle |
| FLASH code | 33,508 B | 47,516 B | +14,008 B (+42 %) |
| ITCM (RAM1 code) | 21,400 B | 34,504 B | +13,104 B (+61 %) |
| RAM2 free | 4,736 B | 4,736 B | unchanged |

-O3 buys a real **1.39× on the render** (−5.68 ms/frame), with the blend itself
gaining a disproportionate **1.95×** (`-ffast-math` + inlining collapses the
per-pixel composite) and the harmonic field eval **1.33×**. But the effect is
already idle-bound at -Os, so the win is **invisible in the frame cadence**: it
stays 16 fps and the saved 5.7 ms simply moves into `sh_buffer_wait`. **-O3 is
NOT the shipping config** — the full 26-effect Phantasm image overflows FlexRAM
(ITCM) at -O3 (a +61 % ITCM blow-up on this one effect illustrates why); that is
exactly why the ship build is -Os. -O3 is only viable here as a single-effect
profile image.

## Caveats

- **ISR time is included**: `CYCCNT` free-runs, so every scope also absorbs the
  flywheel/DMA/USB ISRs firing inside it — the real shipping condition, not a
  pure-CPU algorithm cost.
- **`filter_blend` parenting**: it nests under `sh_rasterize` (the only scope
  that enters it) and its printed % is parent-relative (~11 % of rasterize,
  ~2.5 % of frame); its calls/cycles are correct and constant every window.
- Per-pixel `filter_blend` scope overhead is folded into the ~87 cyc/blend and
  slightly inflates `sh_rasterize`; the coarse `sh_*` scopes are negligible.
- The first window (`frames 1-64`, `min ≈ 13.6 ms` epoch-start) mixes the
  instance's cold first frames and is skipped when picking the representative.
- **-O3 build** (the `profile_o3` env): this does **not** ship — the shipping
  Phantasm image is -Os because the full roster overflows ITCM at -O3. Use this
  report only as a single-effect optimization-headroom reference.
- Run with uncommitted `HS_PROFILE` instrumentation in
  `effects/SphericalHarmonics.h` (three scopes: `sh_buffer_wait`,
  `sh_timeline_step`, `sh_rasterize`); the per-frame field build
  (`decode_lm` ×2 + `HarmonicField` ctor) is left un-scoped.

## Harness

- `targets/Profile/Profile.ino` — generic single-effect wrapper
  (`-D HS_PROFILE_TARGET=SphericalHarmonics`, `-D HS_PROFILE_WINDOW=<frames>`);
  also dumps the `HS_ISR_PROFILE` column-ISR accumulators (flywheel wake, pixel
  pack, DMA submit) each window.
- `pio run -e profile_o3 -t upload` + `python tools/profile_capture.py` (the -O3
  twin of the shipping `-Os` `just profile SphericalHarmonics`).
