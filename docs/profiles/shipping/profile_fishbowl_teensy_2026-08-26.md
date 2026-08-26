# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Fishbowl`).
Raw capture: `build/prof/fishbowl_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_fishbowl_teensy_2026-08-02.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses mesh scan, transformer, and raster `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Fishbowl profile 70 32` |

Image size:

```text
FLASH: code:70,520, data:149,696, headers:9,160
       free for files:1,802,240
RAM1:  variables:315,136, code:45,720, padding:19,816
       free for local variables:143,616
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 481–512 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fish_build_vertices` averages
12.51 ms/frame; its worst window is 13.23 ms/frame
(frames 481–512). Peak frame render is
**23.22 ms** (frames 193–224);
spilled **0/1088 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 39.28 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 481–512)

```text
frame                         62.48 ms  37.49 Mcyc  100%
  fish_multiline_draw          8.28 ms   4.97 Mcyc   13%
    filter_blend              852.1 us 511.29 kcyc    1%  x6757.9 0.1us/c
  fish_build_vertices         13.23 ms   7.94 Mcyc   21%  x1.0 13229.5us/c
  fish_noise_prepare            0.3 us     171 cyc    0%  x1.0 0.3us/c
  fish_timeline_step          127.8 us  76.71 kcyc    0%  x1.0 127.8us/c
  canvas_clear                 85.2 us  51.16 kcyc    0%  x1.0 85.2us/c
  canvas_buffer_wait          40.76 ms  24.46 Mcyc   65%  x1.0 40762.0us/c
```

Wall min/avg/max = 61.05/62.48/63.90 ms. `fish_build_vertices`
accounts for 21.2% of this measured frame at 13.23 ms/frame.
The complete render is 21.72 ms; `canvas_buffer_wait` contributes
40.76 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 1–32)

```text
frame                         58.78 ms  35.27 Mcyc  100%
  fish_multiline_draw          1.32 ms 793.15 kcyc    2%
    filter_blend              146.2 us  87.72 kcyc    0%  x1177.4 0.1us/c
  fish_build_vertices          1.93 ms   1.16 Mcyc    3%  x1.0 1926.9us/c
  fish_noise_prepare            0.2 us     164 cyc    0%  x1.0 0.2us/c
  fish_timeline_step          127.8 us  76.66 kcyc    0%  x1.0 127.8us/c
  canvas_clear                 85.7 us  51.44 kcyc    0%  x1.0 85.7us/c
  canvas_buffer_wait          55.31 ms  33.19 Mcyc   94%  x1.0 55311.2us/c
```

Wall min/avg/max = 0.55/58.77/63.22 ms. `fish_build_vertices`
accounts for 3.3% of this measured frame at 1.93 ms/frame.
The complete render is 3.46 ms; `canvas_buffer_wait` contributes
55.31 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-pixel figures

The peak dominant-scope window blends 6,757.9 pixels/frame, 0.65× the 10,368-pixel quadrant. `filter_blend` costs 75.7 cycles/blend; `fish_build_vertices` uses 1174.6 cycles per blended pixel.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.8/f  0.47/1.61/28.14 us  cpu 2.99%
isr_pack          143.8/f  6.23/6.69/9.65 us  cpu 1.54%
isr_dma_submit    143.8/f  0.60/0.94/2.50 us  cpu 0.22%
```

- Pack plus submit costs 1.10 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.74% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `fish_build_vertices` — 20.1% of aggregated root time, 12.51 ms/frame: inclusive measured scope in the live driver.
2. `fish_multiline_draw` — 12.5% of aggregated root time, 7.78 ms/frame: inclusive measured scope in the live driver.
3. `filter_blend` — 1.3% of aggregated root time, 0.80 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches mesh scan, transformer, and raster `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- No dwell-compression knob applies to this steady capture. Capture knobs were `none`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=Fishbowl` and
`HS_PROFILE_WINDOW=32`; `just profile Fishbowl` performs the
locked build, flash, capture, and artifact attestation.
