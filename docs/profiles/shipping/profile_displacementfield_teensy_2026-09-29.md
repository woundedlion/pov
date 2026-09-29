# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-09-29 09:09 on COM3.
Replaces `profile_displacementfield_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; the fused ring-stack scan (`Scan::DistortedRingStack::draw`) and `draw_rings` compile at -O3 in ITCM; the stack's table build and both bake helpers run from hot flash |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `e1fd2e1ec` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` so one instance covers a full NOISE → BALLS → NOISE cycle |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 150 32 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size (shipping `phantasm` image): `RAM1: variables:314784, code:172136, padding:24472, free for local variables:12896`.

Exactness cross-check: window frames 2177–2208 root counter cyc ÷ 600 MHz matches the measured wall sum within **1.7 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `df_timeline_step` averages 28.47 ms/f; its worst window is 33.63 ms/f (frames 2177–2208). Peak frame render is **35.27 ms** (frames 2177–2208), and **0/2368** frames spilled.

A same-board capture of the previous tip (`29faba01d`, taken immediately before this one on COM3) recorded peak 60.41 ms and 0/2368 spilled.

A display window is 62.5 ms, so render at or under it holds 16 fps; every phase now leaves ≥ 27 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: NOISE fade-in (150 frames) → NOISE dwell (600) → fade-out (150) → BALLS spawning (900) → drain, then the next NOISE phase from frame ~1950.

### NOISE dwell (window frames 2177–2208, worst of the capture)

```
frame                     62.46 ms   37.47 Mcyc   100%
  pov_preserve_half       142.3 us    85.4 kcyc     0%
  df_timeline_step        33.63 ms   20.18 Mcyc    54%
    df_draw_rings         33.49 ms   20.10 Mcyc    54%
      df_lut_bake          8.99 ms    5.40 Mcyc    14%  x46.1  195 us/c
        df_octave_noise    5.07 ms    3.04 Mcyc     8%  x46.1  110 us/c
        df_hue_table_prep  1.15 ms    0.69 Mcyc     2%  x23.6   49 us/c
      df_chunk_cull       600.1 us    0.36 Mcyc     1%  x47.7   13 us/c
      df_fused_scan       23.10 ms   13.86 Mcyc    37%
        ring_stack_table   1.82 ms    1.09 Mcyc     3%
        filter_blend      813.7 us    0.49 Mcyc     1%  x10303  47 cyc/c
  canvas_clear             84.9 us    51.0 kcyc     0%
  canvas_buffer_wait      28.60 ms   17.16 Mcyc    46%
  df_prepare_fields         0.6 us     0.4 kcyc     0%
```

Wall min/avg/max = 60.77/62.46/64.13 ms. The fused scan is two thirds of the render; the octave-grid noise evaluation is the largest bake cost.

### BALLS steady (window frames 1537–1568)

```
frame                     62.39 ms   37.44 Mcyc   100%
  pov_preserve_half       143.3 us    86.0 kcyc     0%
  df_timeline_step        27.44 ms   16.47 Mcyc    44%
    df_draw_rings         27.21 ms   16.32 Mcyc    44%
      df_lut_bake          6.27 ms    3.76 Mcyc    10%  x40.3  155 us/c
        df_hue_table_prep 326.1 us    0.20 Mcyc     1%  x36.8    9 us/c
      df_chunk_cull       516.5 us    0.31 Mcyc     1%  x40.3   13 us/c
      df_fused_scan       19.46 ms   11.67 Mcyc    31%
        ring_stack_table   1.51 ms    0.90 Mcyc     2%
        filter_blend      741.1 us    0.44 Mcyc     1%  x9336   48 cyc/c
  canvas_clear             84.4 us    50.6 kcyc     0%
  canvas_buffer_wait      34.70 ms   20.82 Mcyc    56%
  df_prepare_fields        11.9 us     7.1 kcyc     0%
```

Wall min/avg/max = 60.85/62.39/63.81 ms. Balls bake only across the azimuth arc each cap covers, and the hue table fills only to the largest reached shift, so the ball phase now costs less than the noise dwell.

### Per-pixel figures

`filter_blend` ran 10,303 times per frame in the NOISE peak window (0.99× the 10,368 px quadrant) at 47 cycles per blend; the fused scan spends 1,345 cycles per blended pixel there.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.5/1.7/11.5 us  cpu 3.11%
isr_pack         144/frame  min/avg/max 6.4/7.1/9.8 us   cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/1.2 us   cpu 0.21%
```

- Submit costs a seventh of pack per column.
- The ISR share (~4.9%) leaves ~59.4 ms of each 62.5 ms window for render; every phase fits with ≥ 26 ms to spare.

## Summary ranking

1. `df_fused_scan` — 37% of the NOISE peak window, 23.10 ms (table build 1.82 ms).
2. `df_lut_bake` — 14%, 8.99 ms, of which octave noise 5.07 ms.
3. `df_chunk_cull` — 1%, 0.60 ms.

Same board, previous tip → this tip: NOISE dwell render 54.99 → 31.97 ms/f, BALLS 40.92 → 26.66 ms/f, peak 60.41 → 35.27 ms.

README cells: peak 🟢 35.27, spilled 🟢 0/2368 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `df_fused_scan`; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures.
- The epoch stretch keeps one effect instance across the whole cycle; it does not change per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh DisplacementField profile 150 32 -D HS_PROFILE_EPOCH_REVS=1600` builds, flashes and captures under the device lock.
