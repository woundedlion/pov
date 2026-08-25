# KaleidoscopeHexOil on-device profile — Teensy 4.0, segmented mode (2026-08-17, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexOil`).
Raw capture: `build/prof/kaleidoscopehexoil_ship.log`. First profile of this effect;
no prior report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live |
| Image | `profile` env: Phantasm shipping flags (`-Os` + newlib-nano) + `HS_PROFILE_ENABLE`; the `Scan::Shader` draw path crosses the landed `HS_O3` regions in core/render/scan.h |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master (COM3) |
| Effect | KaleidoscopeHexOil 288×144, single-entry playlist, tip `fd982e2e` (branch `kaleidoscopehexoil/promote`) |
| Method | Shared composed-effect `HS_PROFILE` scopes, window = 16 frames, 140 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` (200 s epoch so the 2-preset cycle never crosses it) |
| Reproduce | `bash tools/profile_one.sh KaleidoscopeHexOil profile 140 16 '-D HS_PROFILE_EPOCH_REVS=1600'` |

Image size: `FLASH: code:75648, data:149492, headers:8332` / `RAM1:
variables:315232, code:40360, padding:25176, free:143520` / `RAM2:
variables:520064, free:4224`.

Exactness cross-check: window frames 1553–1568 root counter 599,328,742 cyc ÷
600 MHz = 998,881 µs against the measured wall sum 998,886 µs — within
**4.8 ppm** (`tools/parse_profile.py build/prof/kaleidoscopehexoil_ship.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw`
35.0–36.9 ms/f across clean holds, worst window 40.0 ms/f mean (frames
1553–1568), peak frame render **40.52 ms**, spilled **0/2208** frames (0%).

A display window is 62.5 ms; the effect shades one quadrant ≈ 10,368 px
through the fully inlined pullback pipeline. Every frame of the capture lands
in exactly one display window — steady 16 fps through both holds and both
crossfades, with ≥ 22 ms of window margin at the peak. The
`canvas_buffer_wait` scope is the round-up idle to the next display flip, by
design.

## Phase-by-phase readout

Phase schedule: 600-frame hold on `kaleidoscope-hex-oil` → 480-frame parameter
crossfade → 600-frame hold on `kaleidoscope-hex-oil-2` → crossfade back (wrap
confirmed by the `Preset: 1/2` marker at ~frame 1680). Render cost also
drifts ±3 ms within a hold as the projection/camera wander sweeps the view
across the folded noise field, so the regimes below bracket the capture.

### kaleidoscope-hex-oil-2 hold (window frames 1553–1568, worst of the capture)

```
frame                  62.43 ms  37.46 Mcyc  100%
  fx_shader_draw       36.88 ms  22.13 Mcyc   59%  x1  36.9 ms/call
  canvas_buffer_wait   22.44 ms  13.47 Mcyc   35%
  fx_advance            2.20 ms   1.32 Mcyc    3%
  fx_prepare_frame      0.75 ms   0.45 Mcyc    1%
  canvas_clear          0.09 ms   0.05 Mcyc    0%
  fx_timeline_step      0.06 ms   0.03 Mcyc    0%
```

Wall min/avg/max = 61.29/62.43/63.61 ms — the densest orientation of the
tighter (scale 3.66) displacement field. Render peaks at 40.52 ms here; the
rest of the frame is display-sync idle.

### Crossfade toward kaleidoscope-hex-oil (window frames 1841–1856)

```
frame                  62.45 ms  37.47 Mcyc  100%
  fx_shader_draw       32.11 ms  19.27 Mcyc   51%  x1  32.1 ms/call
  canvas_buffer_wait   26.98 ms  16.19 Mcyc   43%
  fx_advance            2.19 ms   1.31 Mcyc    3%
  fx_prepare_frame      0.89 ms   0.53 Mcyc    1%
  canvas_clear          0.09 ms   0.05 Mcyc    0%
  fx_timeline_step      0.11 ms   0.07 Mcyc    0%
```

Wall min/avg/max = 61.88/62.45/63.01 ms — mid-crossfade the surface scale and
mapping frequency interpolate between the endpoints and the cost sits between
the two holds. The per-frame hue-rotation LUT rebuild while
`hue_shift_amount` interpolates is inside `fx_prepare_frame` and stays under
1 ms/f, so the transition adds no regime of its own.

### Per-preset table

Wrap-to-0 confirmed (`Preset: 2/2` at ~frame 600, `Preset: 1/2` at ~frame
1680); 38/67/33 windows attribute to the initial hold, preset 2, and the
post-wrap return respectively (a preset's windows include the transition that
follows it).

| # | preset | blended px/f | fx_shader_draw ms | peak render ms | fps |
|---|---|--:|--:|--:|--:|
| 2 | kaleidoscope-hex-oil-2 | full quadrant | 36.88 | 40.52 | 16 |
| 1 | kaleidoscope-hex-oil | full quadrant | 35.01 | 39.63 | 16 |

### Per-pixel figures

The shader shades the full 10,368-px quadrant every frame (no sparse
raster, no `filter_blend` pass). At the 36.88 ms clean-hold cost that is
3.56 µs ≈ **2,134 cycles per pixel** through camera → hex-prism fold →
direct-noise displacement → stereographic projection → spiral → material →
palette.

## Column-ISR / DMA marshaling cost

```
isr_wake        1152/frame  min/avg/max 0.7/2.0/25.2 us  cpu 3.6%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.2 us   cpu 1.6%
isr_dma_submit   144/frame  min/avg/max 0.8/0.9/5.5 us   cpu 0.2%
```

- Submit is ~7× cheaper than pack per call; pack dominates the marshaling
  cost at ~15.5 µs-equivalent per column pair.
- Wire transfer runs async under DMA; only the submit setup is CPU.
- ISR share ≈ 5.4% ⇒ ~59 ms of render budget per 62.5 ms window; the
  40.5 ms peak leaves ~31% headroom, so no speedup is required.

## Summary ranking

1. `fx_shader_draw` — 51–59% of the frame, 35.0–36.9 ms clean-hold (40.5 ms
   peak): the per-pixel pullback pipeline; the tighter noise field of
   `kaleidoscope-hex-oil-2` costs ~1.9 ms/f more than `kaleidoscope-hex-oil`.
2. `canvas_buffer_wait` — 35–43%, 22–27 ms: display-sync idle, not work.
3. `fx_advance` — 3%, 2.2 ms: clocks, camera walks, palette cycler step.
4. `fx_prepare_frame` — 1%, 0.75–0.89 ms: frame snapshot + hue-LUT upkeep
   (rebuilt per frame only while a crossfade interpolates the hue shift).

No WASM/native ledger entry exists for this effect yet.

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- No `filter_blend` counter appears: the effect is a full-quadrant shader
  with no blended-raster pass.
- Selective-O3 shipping config: the `Scan::Shader` draw path crosses the
  landed `HS_O3` regions in core/render/scan.h; no effect-specific region
  exists.
- `HS_PROFILE_EPOCH_REVS=1600` stretched the epoch so the 1,680-frame cycle
  fits one instance; the stretch changes hold lengths, never per-frame cost.
- Captured from worktree branch `kaleidoscopehexoil/promote` (`fd982e2e`), not yet
  landed on master at capture time; the tree was clean at the profiled tip.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeHexOil`,
`HS_PROFILE_WINDOW=16`; `just profile KaleidoscopeHexOil` = build + flash + capture.
