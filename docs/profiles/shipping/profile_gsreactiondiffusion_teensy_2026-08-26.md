# GSReactionDiffusion on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile GSReactionDiffusion`).
Raw capture: `build/prof/gsreactiondiffusion_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_gsreactiondiffusion_teensy_2026-08-09.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses the GS physics/shade `HS_O3_FN` loops and `ReactionDiffusionBase` accumulation helpers |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | GSReactionDiffusion 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 130 s capture, `-D HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh GSReactionDiffusion profile 130 32 "-D HS_PROFILE_EPOCH_REVS=1200"` |

Image size:

```text
FLASH: code:61,128, data:239,948, headers:9,196
       free for files:1,721,344
RAM1:  variables:315,232, code:29,528, padding:3,240
       free for local variables:176,288
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 449–480 root counter
cycles ÷ 600 MHz match the measured wall sum within **0.6 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `grd_render` averages
48.00 ms/frame; its worst window is 53.77 ms/frame
(frames 449–480). Peak frame render is
**55.26 ms** (frames 481–512);
spilled **0/2048 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 7.24 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 449–480)

```text
frame                         62.46 ms  37.47 Mcyc  100%
  grd_render                  53.77 ms  32.26 Mcyc   86%
    grd_rasterize             41.27 ms  24.76 Mcyc   66%
      grd_shader_draw         38.33 ms  23.00 Mcyc   61%  x1.0 38329.7us/c
      grd_cull_flags           2.57 ms   1.54 Mcyc    4%  x1.0 2573.7us/c
      grd_orient              366.7 us 220.05 kcyc    1%  x1.0 366.7us/c
    grd_simulate              11.96 ms   7.18 Mcyc   19%  x1.0 11958.5us/c
  rd_timeline_step             27.9 us  16.76 kcyc    0%  x1.0 27.9us/c
  canvas_clear                 89.8 us  53.90 kcyc    0%  x1.0 89.8us/c
  canvas_buffer_wait           8.57 ms   5.14 Mcyc   14%  x1.0 8568.2us/c
```

Wall min/avg/max = 61.22/62.45/63.59 ms. `grd_render`
accounts for 86.1% of this measured frame at 53.77 ms/frame.
The complete render is 53.89 ms; `canvas_buffer_wait` contributes
8.57 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 993–1024)

```text
frame                         62.36 ms  37.42 Mcyc  100%
  grd_render                  27.65 ms  16.59 Mcyc   44%
    grd_palette_rebake         66.4 us  39.86 kcyc    0%  x0.0 2125.0us/c
    grd_rasterize             15.57 ms   9.34 Mcyc   25%
      grd_shader_draw         11.48 ms   6.89 Mcyc   18%  x1.0 11483.1us/c
      grd_cull_flags           3.71 ms   2.23 Mcyc    6%  x1.0 3710.8us/c
      grd_orient              372.8 us 223.72 kcyc    1%  x1.0 372.8us/c
    grd_simulate              11.95 ms   7.17 Mcyc   19%  x1.0 11952.3us/c
  rd_timeline_step             28.2 us  16.93 kcyc    0%  x1.0 28.2us/c
  canvas_clear                 86.0 us  51.61 kcyc    0%  x1.0 86.0us/c
  canvas_buffer_wait          34.60 ms  20.76 Mcyc   55%  x1.0 34600.3us/c
```

Wall min/avg/max = 54.73/62.36/68.25 ms. `grd_render`
accounts for 44.3% of this measured frame at 27.65 ms/frame.
The complete render is 27.76 ms; `canvas_buffer_wait` contributes
34.60 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `grd_render` uses 3111.6 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1160.3/f  0.43/1.72/20.41 us  cpu 3.19%
isr_pack          143.9/f  6.24/7.41/16.50 us  cpu 1.71%
isr_dma_submit    143.9/f  0.58/0.94/8.29 us  cpu 0.22%
```

- Pack plus submit costs 1.20 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 5.11% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `grd_render` — 76.9% of aggregated root time, 48.00 ms/frame: inclusive measured scope in the live driver.
2. `grd_rasterize` — 57.6% of aggregated root time, 35.93 ms/frame: inclusive measured scope in the live driver.
3. `grd_shader_draw` — 51.5% of aggregated root time, 32.17 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches the GS physics/shade `HS_O3_FN` loops and `ReactionDiffusionBase` accumulation helpers; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `-D HS_PROFILE_EPOCH_REVS=1200`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=GSReactionDiffusion` and
`HS_PROFILE_WINDOW=32`; `just profile GSReactionDiffusion` performs the
locked build, flash, capture, and artifact attestation.
