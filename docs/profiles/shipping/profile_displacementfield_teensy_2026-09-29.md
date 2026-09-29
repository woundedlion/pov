# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-09-29, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile DisplacementField`).
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-09-29 00:06 on COM3.
Replaces `profile_displacementfield_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; the fused ring-stack scan (`Scan::DistortedRingStack::draw`), `draw_rings` and both bake helpers compile at -O3 |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `8031d4d6f` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 150 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` so one instance covers a full NOISE → BALLS → NOISE cycle |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 150 32 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size (shipping `phantasm` image): `RAM1: variables:314784, code:176920, padding:19688, free for local variables:12896`.

Exactness cross-check: window frames 2177–2208 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.2 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `df_timeline_step` averages 28.06 ms/f; its worst window is 33.10 ms/f (frames 2177–2208). Peak frame render is **34.70 ms** (frames 2177–2208), and **0/2368** frames spilled.

A same-board capture of the previous tip (`29faba01d`, taken immediately before this one on COM3) recorded peak 60.41 ms and 0/2368 spilled.

A display window is 62.5 ms, so render at or under it holds 16 fps; every phase now leaves ≥ 27 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: NOISE fade-in (150 frames) → NOISE dwell (600) → fade-out (150) → BALLS spawning (900) → drain, then the next NOISE phase from frame ~1950.

### NOISE dwell (window frames 2177–2208, worst of the capture)

```
frame                     62.48 ms   37.49 Mcyc   100%
  pov_preserve_half       142.3 us    85.4 kcyc     0%
  df_timeline_step        33.10 ms   19.86 Mcyc    53%
    df_draw_rings         32.99 ms   19.79 Mcyc    53%
      df_lut_bake          8.91 ms    5.34 Mcyc    14%  x46.1  193 us/c
        df_octave_noise    4.97 ms    2.98 Mcyc     8%  x46.1  108 us/c
        df_hue_table_prep  1.15 ms    0.69 Mcyc     2%  x23.6   49 us/c
      df_chunk_cull       595.6 us    0.36 Mcyc     1%  x47.7   12 us/c
      df_fused_scan       22.70 ms   13.62 Mcyc    36%
        ring_stack_table   1.83 ms    1.10 Mcyc     3%
        filter_blend      777.9 us    0.47 Mcyc     1%  x10303  45 cyc/c
  canvas_clear             84.3 us    50.6 kcyc     0%
  canvas_buffer_wait      29.14 ms   17.49 Mcyc    47%
  df_prepare_fields         0.3 us     0.2 kcyc     0%
```

Wall min/avg/max = 60.86/62.48/64.03 ms. The fused scan is two thirds of the render; the octave-grid noise evaluation is the largest bake cost.

### BALLS steady (window frames 1537–1568)

```
frame                     62.41 ms   37.45 Mcyc   100%
  pov_preserve_half       145.4 us    87.2 kcyc     0%
  df_timeline_step        27.11 ms   16.27 Mcyc    43%
    df_draw_rings         26.90 ms   16.14 Mcyc    43%
      df_lut_bake          6.29 ms    3.77 Mcyc    10%  x40.3  156 us/c
        df_hue_table_prep 326.6 us    0.20 Mcyc     1%  x36.8    9 us/c
      df_chunk_cull       521.8 us    0.31 Mcyc     1%  x40.3   13 us/c
      df_fused_scan       19.13 ms   11.48 Mcyc    31%
        ring_stack_table   1.52 ms    0.91 Mcyc     2%
        filter_blend      715.1 us    0.43 Mcyc     1%  x9336   46 cyc/c
  canvas_clear             84.6 us    50.8 kcyc     0%
  canvas_buffer_wait      35.05 ms   21.03 Mcyc    56%
  df_prepare_fields        10.6 us     6.3 kcyc     0%
```

Wall min/avg/max = 60.87/62.41/63.79 ms. Balls bake only across the azimuth arc each cap covers, and the hue table fills only to the largest reached shift, so the ball phase now costs less than the noise dwell.

### Per-pixel figures

`filter_blend` ran 10,303 times per frame in the NOISE peak window (0.99× the 10,368 px quadrant) at 45 cycles per blend; the fused scan spends 1,322 cycles per blended pixel there.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.6/1.7/12.0 us  cpu 3.11%
isr_pack         144/frame  min/avg/max 6.4/7.1/9.8 us   cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/4.4 us   cpu 0.21%
```

- Submit costs a seventh of pack per column.
- The ISR share (~4.9%) leaves ~59.4 ms of each 62.5 ms window for render; every phase fits with ≥ 26 ms to spare.

## Summary ranking

1. `df_fused_scan` — 36% of the NOISE peak window, 22.70 ms (table build 1.83 ms).
2. `df_lut_bake` — 14%, 8.91 ms, of which octave noise 4.97 ms.
3. `df_chunk_cull` — 1%, 0.60 ms.

Same board, previous tip → this tip: NOISE dwell render 54.99 → 31.43 ms/f, BALLS 40.92 → 26.34 ms/f, peak 60.41 → 34.70 ms.

README cells: peak 🟢 34.70, spilled 🟢 0/2368 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `df_fused_scan`; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures.
- The epoch stretch keeps one effect instance across the whole cycle; it does not change per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh DisplacementField profile 150 32 -D HS_PROFILE_EPOCH_REVS=1600` builds, flashes and captures under the device lock.
