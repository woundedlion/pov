# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-08-26, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile SphericalHarmonics`).
Raw capture: `build/prof/sphericalharmonics_ship.log`, sourced from the isolated sweep tree.
Replaces `profile_sphericalharmonics_teensy_2026-08-21.md` with the full 2026-08-26 sweep capture.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base; selective-O3 path crosses mesh scan, face-SDF, shading, and transformer `HS_O3` regions |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048` |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size:

```text
FLASH: code:39,656, data:146,432, headers:8,472
       free for files:1,837,056
RAM1:  variables:315,008, code:20,664, padding:12,104
       free for local variables:176,512
RAM2:  variables:520,064
       free for malloc/new:4,224
```

Exactness cross-check: window frames 2449–2464 root counter
cycles ÷ 600 MHz match the measured wall sum within **2.1 ppm**
(`tools/parse_profile.py ... validate`). The parser also confirms monotonic
frames, the expected header, and no epoch reset.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `sh_rasterize` averages
9.39 ms/frame; its worst window is 12.46 ms/frame
(frames 2449–2464). Peak frame render is
**12.64 ms** (frames 913–928);
spilled **0/3488 frames**
(0.0%).

A display window is 62.5 ms. The segmented master renders one quadrant, approximately 10,368 pixels. The peak retains 49.86 ms of the display window. The effect's
`*_buffer_wait` scope is round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: every parsed window contributes to the aggregate. The sections
below show the dominant-scope extremes that bound the measured regimes; preset
ownership and clean holds are tabulated separately when markers are present.

### Peak measured regime (window frames 2449–2464)

```text
frame                         62.44 ms  37.46 Mcyc  100%
  sh_rasterize                12.46 ms   7.47 Mcyc   20%  x1.0 12455.8us/c
  sh_timeline_step             22.5 us  13.51 kcyc    0%  x1.0 22.5us/c
  canvas_clear                 84.9 us  50.96 kcyc    0%  x1.0 84.9us/c
  canvas_buffer_wait          49.87 ms  29.92 Mcyc   80%  x1.0 49869.1us/c
```

Wall min/avg/max = 62.23/62.44/62.55 ms. `sh_rasterize`
accounts for 19.9% of this measured frame at 12.46 ms/frame.
The complete render is 12.57 ms; `canvas_buffer_wait` contributes
49.87 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Lowest-cost measured regime (window frames 2801–2816)

```text
frame                         62.38 ms  37.43 Mcyc  100%
  sh_rasterize                 7.52 ms   4.51 Mcyc   12%  x1.0 7520.8us/c
  sh_timeline_step             24.6 us  14.78 kcyc    0%  x1.0 24.6us/c
  canvas_clear                 84.9 us  50.99 kcyc    0%  x1.0 84.9us/c
  canvas_buffer_wait          54.75 ms  32.85 Mcyc   88%  x1.0 54750.1us/c
```

Wall min/avg/max = 61.66/62.38/62.66 ms. `sh_rasterize`
accounts for 12.1% of this measured frame at 7.52 ms/frame.
The complete render is 7.63 ms; `canvas_buffer_wait` contributes
54.75 ms of flip-alignment idle. The remaining spread is live-driver
and flywheel timing already absorbed by the counters.

### Per-preset table

The parser observed a closed cycle: the marker sequence returned to its first
entry. Rows are ranked by the modal-call-count clean-hold `sh_rasterize` cost;
cadence peaks use per-frame ownership, including the following transition.
The `windows` column records modal-call-count clean windows over owned windows.

| rank | preset/shape/mode | geometry | windows | blended px/f | `sh_rasterize` ms | render ms | fps |
|---:|---|---|--:|--:|--:|--:|--:|
| 1 | `21` | — | 8/8 | — | 12.46 | 12.64 | 16.0 |
| 2 | `20` | — | 8/8 | — | 12.45 | 12.61 | 16.0 |
| 3 | `19` | — | 8/8 | — | 11.29 | 11.50 | 16.0 |
| 4 | `22` | — | 8/8 | — | 11.28 | 11.48 | 16.0 |
| 5 | `13` | — | 10/10 | — | 10.79 | 10.98 | 16.0 |
| 6 | `12` | — | 12/12 | — | 10.78 | 10.98 | 16.0 |
| 7 | `23` | — | 8/8 | — | 9.94 | 10.09 | 16.0 |
| 8 | `18` | — | 8/8 | — | 9.92 | 10.09 | 16.0 |
| 9 | `14` | — | 8/8 | — | 9.54 | 9.74 | 16.0 |
| 10 | `11` | — | 12/12 | — | 9.51 | 9.72 | 16.0 |
| 11 | `24` | — | 8/8 | — | 9.19 | 9.36 | 16.0 |
| 12 | `17` | — | 8/8 | — | 9.18 | 9.32 | 16.0 |
| 13 | `7` | — | 12/12 | — | 9.02 | 9.23 | 16.0 |
| 14 | `6` | — | 8/8 | — | 8.99 | 9.21 | 16.0 |
| 15 | `16` | — | 8/8 | — | 8.85 | 9.03 | 16.0 |
| 16 | `10` | — | 12/12 | — | 8.78 | 8.93 | 16.0 |
| 17 | `15` | — | 8/8 | — | 8.77 | 8.93 | 16.0 |
| 18 | `9` | — | 12/12 | — | 8.47 | 8.65 | 16.0 |
| 19 | `1` | — | 8/8 | — | 8.44 | 8.59 | 16.0 |
| 20 | `8` | — | 12/12 | — | 8.38 | 8.55 | 16.0 |
| 21 | `5` | — | 8/8 | — | 8.37 | 8.52 | 16.0 |
| 22 | `4` | — | 8/8 | — | 8.06 | 8.25 | 16.0 |
| 23 | `3` | — | 8/8 | — | 7.88 | 8.08 | 16.0 |
| 24 | `2` | — | 8/8 | — | 7.86 | 8.05 | 16.0 |

### Per-pixel figures

No `filter_blend` population appears in the peak window; the path writes directly. `sh_rasterize` uses 720.8 cycles per configured sample over 10,368 samples/frame.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1154.1/f  0.43/1.63/15.88 us  cpu 3.01%
isr_pack          143.9/f  6.24/6.82/9.65 us  cpu 1.57%
isr_dma_submit    143.9/f  0.58/0.93/8.39 us  cpu 0.21%
```

- Pack plus submit costs 1.11 ms of CPU per frame.
- LED wire transfer remains asynchronous; the table measures CPU marshaling only.
- Aggregate ISR CPU share is 4.79% of one 62.5 ms window. The peak
  render requires 1.00× the one-window budget, so the cadence target
  needs no speedup.

## Summary ranking

1. `sh_rasterize` — 15.0% of aggregated root time, 9.39 ms/frame: inclusive measured scope in the live driver.
2. `canvas_clear` — 0.1% of aggregated root time, 0.08 ms/frame: inclusive measured scope in the live driver.
3. `sh_timeline_step` — 0.0% of aggregated root time, 0.02 ms/frame: inclusive measured scope in the live driver.

No matched current WASM/native capture was used; this ranking is the on-device
segmented-driver result.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` parents under whichever measured scope reaches it first; its
  subtree is hidden in windows where that parent has zero calls. Calls are an
  approximation of blended pixels, and the scope itself adds per-pixel
  instrumentation overhead only in profile images.
- Shipping selective O3 reaches mesh scan, face-SDF, shading, and transformer `HS_O3` regions; the rest of the image retains the
  `-Os` base policy.
- Marker-driven dwell/transition compression changes time spent in each preset, not its per-frame cost. Capture knobs were `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048`.
- Provenance attests a clean source tree at `63268c3768ae91c8e2b47acd728db641570179eb`; no
  uncommitted source state was profiled.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=SphericalHarmonics` and
`HS_PROFILE_WINDOW=16`; `just profile SphericalHarmonics` performs the
locked build, flash, capture, and artifact attestation.
