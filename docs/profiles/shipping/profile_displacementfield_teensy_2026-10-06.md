# DisplacementField on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/displacementfield_ship.log`, captured 2026-10-06 19:24 on COM3.
Replaces `profile_displacementfield_teensy_2026-09-29.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean main checkout at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DisplacementField 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 160 s capture, `-D HS_PROFILE_EPOCH_REVS=1600` so one instance covers a full NOISE → BALLS → NOISE cycle |
| Reproduce | `bash tools/profile_one.sh DisplacementField profile 160 32 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size (`profile` env, this effect only): `FLASH: code:80756, data:155052, headers:8924` / `RAM1: variables:315264, code:33752, padding:31784, free:143488` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2177–2208 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.5 ppm** (`tools/parse_profile.py ... validate`, VALID).

Method note: `tools/profile_sweep.sh` budgets DisplacementField at 70 s with no epoch override (group `g1_ship`). That covers only the first NOISE phase and the start of BALLS and misses the second NOISE dwell, where this pass peaks. This report uses a 160 s capture under a stretched epoch instead.

## Frame cadence

**Pass aggregate**: `df_timeline_step` averages 27.47 ms/f; its worst window is 32.16 ms/f (frames 2177–2208). Peak frame render is **33.66 ms** (frame 2199), and **0/2527** frames spilled. Setup frame 1 is excluded from both; it rendered 19.11 ms.

The previous shipping report (2026-09-29 09:09) recorded peak 🟢 35.27 and spilled 🟢 0/2368 (0.00%). Both passes peak in the same window of the second NOISE dwell (frames 2177–2208): its `df_timeline_step` reads 33.63 → 32.16 ms/f and the peak 35.27 → 33.66 ms (−4.6%).

A display window is 62.5 ms, so render at or under it holds 16 fps; every phase leaves at least 28.84 ms of margin. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: NOISE fade-in (frames 1–~160, 21.97 → 31.12 ms/f) → NOISE dwell (to frame 800, 29.40–31.70 ms/f) → fade-out (801–896) → BALLS spawning (897–928, 15.04 ms/f) → BALLS steady (to ~1824, 23.8–26.8 ms/f) → BALLS drain into NOISE (1857–1920, bottoming at 12.47 ms/f) → second NOISE fade-in and dwell from frame ~1921 (30.2–32.2 ms/f).

### Second NOISE dwell (frames 2177–2208, worst of the capture)

```
frame                         62.48 ms    37.49 Mcyc   100%
  pov_preserve_half           144.7 us     86.8 kcyc     0%
  df_timeline_step            32.16 ms    19.30 Mcyc    51%
    df_draw_rings             32.04 ms    19.22 Mcyc    51%
      df_lut_bake              8.96 ms     5.38 Mcyc    14%  x46.1  194 us/c
        df_octave_noise        5.12 ms     3.07 Mcyc     8%  x46.1  111 us/c
        df_hue_table_prep      1.15 ms    689.7 kcyc     2%  x23.6  49 us/c
      df_chunk_cull           702.2 us    421.3 kcyc     1%  x47.7  15 us/c
      df_fused_scan           21.51 ms    12.90 Mcyc    34%
        filter_blend          806.5 us    483.9 kcyc     1%  x10303  47 cyc/c
        ring_stack_table       1.29 ms    774.6 kcyc     2%
  canvas_clear                 86.2 us     51.7 kcyc     0%
  canvas_buffer_wait          30.09 ms    18.05 Mcyc    48%
  df_prepare_fields             0.3 us      192 cyc      0%
```

Wall min/avg/max = 60.96/62.48/63.95 ms. Per-frame values are window averages; `xN` is calls per frame. The fused scan is two thirds of the render and octave noise the largest bake cost. Against the previous report's same window, the fused scan reads 23.10 → 21.51 ms/f and its table build 1.82 → 1.29 ms/f, while the bake is flat (8.99 → 8.96).

### First NOISE dwell (frames 577–608)

```
frame                         62.51 ms    37.51 Mcyc   100%
  pov_preserve_half           144.0 us     86.4 kcyc     0%
  df_timeline_step            31.11 ms    18.67 Mcyc    50%
    df_draw_rings             30.99 ms    18.60 Mcyc    50%
      df_lut_bake              8.67 ms     5.20 Mcyc    14%  x41.7  208 us/c
        df_octave_noise        5.12 ms     3.07 Mcyc     8%  x41.7  123 us/c
        df_hue_table_prep      1.04 ms    622.3 kcyc     2%  x20.9  50 us/c
      df_chunk_cull           657.9 us    394.7 kcyc     1%  x43.4  15 us/c
      df_fused_scan           20.83 ms    12.50 Mcyc    33%
        filter_blend          794.0 us    476.4 kcyc     1%  x10160  47 cyc/c
        ring_stack_table       1.12 ms    669.8 kcyc     2%
  canvas_clear                 86.5 us     51.9 kcyc     0%
  canvas_buffer_wait          31.17 ms    18.70 Mcyc    50%
  df_prepare_fields             0.4 us      231 cyc      0%
```

Wall min/avg/max = 60.91/62.51/64.18 ms. This window holds the first dwell's peak frame (33.14 ms). It bakes 41.7 rings per frame against the second dwell's 46.1, so it runs about 1 ms/f cheaper.

### BALLS spawning (frames 897–928)

```
frame                         62.64 ms    37.58 Mcyc   100%
  pov_preserve_half           146.9 us     88.1 kcyc     0%
  df_timeline_step            15.04 ms     9.03 Mcyc    24%
    df_draw_rings             14.89 ms     8.94 Mcyc    24%
      df_lut_bake              1.66 ms    995.9 kcyc     3%  x14.2  117 us/c
        df_hue_table_prep      83.3 us     50.0 kcyc     0%  x12.4  7 us/c
      df_chunk_cull           212.2 us    127.3 kcyc     0%  x15.2  14 us/c
      df_fused_scan           12.36 ms     7.41 Mcyc    20%
        filter_blend          726.8 us    436.1 kcyc     1%  x9273  47 cyc/c
        ring_stack_table      995.8 us    597.5 kcyc     2%
  canvas_clear                 86.0 us     51.6 kcyc     0%
  canvas_buffer_wait          47.35 ms    28.41 Mcyc    76%
  df_prepare_fields             3.8 us      2.3 kcyc     0%
```

Wall min/avg/max = 59.41/62.64/65.71 ms. Octave noise is gone and only 14 rings bake per frame while the first balls appear.

### BALLS steady (frames 1025–1056)

```
frame                         62.39 ms    37.43 Mcyc   100%
  pov_preserve_half           143.3 us     86.0 kcyc     0%
  df_timeline_step            26.78 ms    16.07 Mcyc    43%
    df_draw_rings             26.55 ms    15.93 Mcyc    43%
      df_lut_bake              6.68 ms     4.01 Mcyc    11%  x44.0  152 us/c
        df_hue_table_prep     324.6 us    194.7 kcyc     1%  x38.4  8 us/c
      df_chunk_cull           660.8 us    396.5 kcyc     1%  x44.1  15 us/c
      df_fused_scan           18.12 ms    10.87 Mcyc    29%
        filter_blend          732.8 us    439.7 kcyc     1%  x9317  47 cyc/c
        ring_stack_table       1.20 ms    722.9 kcyc     2%
  canvas_clear                 86.3 us     51.8 kcyc     0%
  canvas_buffer_wait          35.36 ms    21.21 Mcyc    57%
  df_prepare_fields            11.2 us      6.7 kcyc     0%
```

Wall min/avg/max = 60.63/62.39/63.64 ms. The costliest BALLS window: bake 2.28 ms/f and scan 3.39 ms/f below the second NOISE dwell. BALLS windows span 23.80–26.78 ms/f through frame 1824.

### BALLS drain into NOISE (frames 1889–1920)

```
frame                         62.77 ms    37.66 Mcyc   100%
  pov_preserve_half           148.9 us     89.3 kcyc     0%
  df_timeline_step            12.47 ms     7.48 Mcyc    20%
    df_draw_rings             12.32 ms     7.39 Mcyc    20%
      df_lut_bake              1.20 ms    719.0 kcyc     2%  x13.0  92 us/c
        df_octave_noise       264.0 us    158.4 kcyc     0%  x2.3  114 us/c
        df_hue_table_prep      46.2 us     27.7 kcyc     0%  x6.8  7 us/c
      df_chunk_cull           190.3 us    114.2 kcyc     0%  x13.2  14 us/c
      df_fused_scan           10.33 ms     6.20 Mcyc    16%
        filter_blend          722.1 us    433.3 kcyc     1%  x9256  47 cyc/c
        ring_stack_table      940.8 us    564.5 kcyc     1%
  canvas_clear                 86.3 us     51.8 kcyc     0%
  canvas_buffer_wait          50.05 ms    30.03 Mcyc    80%
  df_prepare_fields             9.6 us      5.8 kcyc     0%
```

Wall min/avg/max = 62.05/62.77/70.63 ms. The cheapest window of the capture: 13 rings bake per frame as the balls drain and octave noise restarts (2.3 calls per frame). The 70.63 ms wall maximum is display-sync wait, not render; its peak frame render is 22.72 ms.

### Per-pixel figures

`filter_blend` ran 10,303 times per frame in the second NOISE dwell window (0.99× the 10,368 px quadrant) at 47 cycles per blend; the fused scan spends 1,252 cycles per blended pixel there.

## Column-ISR / DMA marshaling cost

```
isr_wake        1153/frame  min/avg/max 0.5/1.7/14.7 us  cpu 3.05%
isr_pack         144/frame  min/avg/max 6.4/7.1/10.0 us  cpu 1.62%
isr_dma_submit   144/frame  min/avg/max 0.6/0.9/7.1 us  cpu 0.21%
```

Second NOISE dwell window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `df_timeline_step` — 51% of the second NOISE dwell window, 32.16 ms/f: `df_fused_scan` 21.51 ms (table build 1.29), `df_lut_bake` 8.96 ms (octave noise 5.12), `df_chunk_cull` 0.70 ms.
2. `pov_preserve_half` — 0% of the second NOISE dwell window, 0.14 ms/f.
3. `canvas_clear` — 0% of the second NOISE dwell window, 0.09 ms/f.
4. `df_prepare_fields` — 0% of the second NOISE dwell window, under 0.01 ms/f.

Previous report → this one: second NOISE dwell window 33.63 → 32.16 ms/f, BALLS window 1537–1568 27.44 → 26.16 ms/f, peak 35.27 → 33.66 ms.

README cells: peak 🟢 33.66, spilled 🟢 0/2527 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under `df_fused_scan`; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- The epoch stretch keeps one effect instance across the whole cycle; it does not change per-frame cost. `validate` reports 0 epoch resets.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `e1fd2e1ec` (also on COM3), and the delta between them is not attributed here.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DisplacementField`, `HS_PROFILE_WINDOW=32`, `HS_PROFILE_EPOCH_REVS=1600`; `bash tools/profile_one.sh DisplacementField profile 160 32 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes and captures under the device lock.
