# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-09-28, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/dreamballs_ship.log`, captured 2026-09-28 18:45 on COM3.
Replaces `profile_dreamballs_teensy_2026-08-26.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, tip `97eb0bf78` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture, `-D HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` |

Image size (`profile` env, this effect only): `FLASH: code:108904, data:191220, headers:9124` / `RAM1: variables:315296, code:44248, padding:21288, free:143456` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2673–2688 root counter cyc ÷ 600 MHz matches the measured wall sum within **2.8 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `db_timeline_step` averages 25.43 ms/f; its worst window is 40.84 ms/f (frames 2673–2688). Peak frame render is **42.66 ms** (frame 2765), and **0/3647** frames spilled. Setup frame 1 is excluded from both; it rendered 59.68 ms.

The previous shipping report (2026-08-26 02:21) recorded peak 🟢 38.73 (11) and spilled 🟢 0/3648 (0.0%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 10 preset entries; each owns its hold and the transition that follows it. The block below is the window holding the pass's peak frame.

### Peak window (frames 2753–2768)

```
frame                     62.42 ms   37.45 Mcyc   100%
  pov_preserve_half       142.7 us    85.6 kcyc     0%
  db_timeline_step        40.18 ms   24.11 Mcyc    64%
    db_draw               40.10 ms   24.06 Mcyc    64%
      db_draw_scene       40.10 ms   24.06 Mcyc    64%
        db_mesh_plot      38.64 ms   23.18 Mcyc    62%  x6.0  6439 us/c
          filter_blend     5.64 ms    3.38 Mcyc     9%  x43118  78 cyc/c
        db_orient         259.4 us   155.6 kcyc     0%  x6.0  43 us/c
        db_displace        1.02 ms   614.1 kcyc     2%  x6.0  171 us/c
  canvas_clear             84.6 us    50.8 kcyc     0%
  canvas_buffer_wait      22.02 ms   13.21 Mcyc    35%
```

Wall min/avg/max = 59.56/62.42/65.86 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `db_timeline_step` is the costliest modal-call-count window of each entry.

| Entry | Meta | Peak render ms | Spilled/frames | Clean `db_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 9 | — | 🟢 42.66 | 0/320 | 40.84 | 20/20 |
| 1 | — | 🟢 39.88 | 0/638 | 38.40 | 40/40 |
| 10 | — | 🟢 32.17 | 0/320 | 29.79 | 20/20 |
| 8 | — | 🟢 31.91 | 0/320 | 30.35 | 20/20 |
| 4 | — | 🟢 26.09 | 0/320 | 23.45 | 20/20 |
| 3 | — | 🟢 23.82 | 0/320 | 21.77 | 20/20 |
| 2 | — | 🟢 22.55 | 0/449 | 21.11 | 28/28 |
| 7 | — | 🟢 21.01 | 0/320 | 19.75 | 20/20 |
| 6 | — | 🟢 15.82 | 0/320 | 13.92 | 20/20 |
| 5 | — | 🟢 14.41 | 0/320 | 13.66 | 20/20 |

### Per-pixel figures

`filter_blend` ran 43,118 times per frame in the peak window at 78 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake         1152/frame  min/avg/max 0.6/1.7/11.0 us  cpu 3.06%
isr_pack          144/frame  min/avg/max 6.2/6.8/9.2 us  cpu 1.56%
isr_dma_submit    144/frame  min/avg/max 0.8/0.9/1.0 us  cpu 0.21%
```

The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `db_timeline_step` — 64% of the peak window, 40.18 ms/f.
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 42.66 (10), spilled 🟢 0/3647 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately.
- Selective-O3: `97eb0bf78` changed placement only for the per-pixel pullback, noise and projection kernels; this effect's own `HS_O3` regions are unchanged.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` builds, flashes and captures under the device lock.
