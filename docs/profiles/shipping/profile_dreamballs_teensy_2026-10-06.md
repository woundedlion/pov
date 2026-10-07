# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-10-06, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw capture: `build/prof/dreamballs_ship.log`, captured 2026-10-06 18:28 on COM3.
Replaces `profile_dreamballs_teensy_2026-09-28.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base with the landed `HS_O3` regions; clean tree at the tip |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, tip `e2f5b0a3d` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture, `-D HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` |

Image size (`profile` env, this effect only): `FLASH: code:114516, data:194156, headers:8764` / `RAM1: variables:315296, code:44936, padding:20600, free:143456` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 2673–2688 root counter cyc ÷ 600 MHz matches the measured wall sum within **0.9 ppm** (`tools/parse_profile.py ... validate`, VALID).

## Frame cadence

**Pass aggregate**: `db_timeline_step` averages 24.96 ms/f; its worst window is 40.00 ms/f (frames 2673–2688). Peak frame render is **41.44 ms** (frame 2734), and **0/3647** frames spilled. Setup frame 1 is excluded from both; it rendered 60.16 ms.

The previous shipping report (2026-09-28 18:45) recorded peak 🟢 42.66 (10) and spilled 🟢 0/3647 (0.00%).

A display window is 62.5 ms, so render at or under it holds 16 fps. The `canvas_buffer_wait` scope is the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: 10 preset entries; each owns its hold and the transition that follows it. The capture runs entry 1 through 10, wraps back to entry 1 and ends on entry 2. The block below is the window holding the pass's peak frame.

### Peak window (frames 2721–2736)

```
frame                    62.40 ms  37.44 Mcyc   100%
  pov_preserve_half      144.9 us   87.0 kcyc     0%
  db_timeline_step       39.85 ms  23.91 Mcyc    64%
    db_draw              39.78 ms  23.87 Mcyc    64%
      db_draw_scene      39.78 ms  23.87 Mcyc    64%
        db_mesh_plot     38.29 ms  22.97 Mcyc    61%  x6.0  6382 us/c
          filter_blend    5.71 ms   3.43 Mcyc     9%  x43459  79 cyc/c
        db_orient        265.6 us  159.4 kcyc     0%  x6.0  44 us/c
        db_displace       1.03 ms  617.3 kcyc     2%  x6.0  171 us/c
  canvas_clear            84.3 us   50.6 kcyc     0%
  canvas_buffer_wait     22.32 ms  13.39 Mcyc    36%
```

Wall min/avg/max = 61.13/62.40/64.04 ms. Per-frame values are window averages; `xN` is calls per frame.

### Per-preset table

Buckets from the per-frame owner stamps, setup frame 1 excluded; clean-hold `db_timeline_step` is from `parse_profile.py ... presets` (the costliest modal-call-count window of each entry, windows straddling an advance excluded). The cycle wraps back to entry 1 within the capture (validate: cycle returns to its first index).

| Entry | Meta | Peak render ms | Spilled/frames | Clean `db_timeline_step` ms/f | Clean windows |
|---|---|--:|--:|--:|--:|
| 9 | — | 🟢 41.44 | 0/320 | 40.00 | 19/20 |
| 1 | — | 🟢 38.63 | 0/638 | 37.59 | 38/40 |
| 10 | — | 🟢 31.22 | 0/320 | 29.16 | 19/20 |
| 8 | — | 🟢 30.96 | 0/320 | 29.74 | 19/20 |
| 4 | — | 🟢 25.85 | 0/320 | 22.93 | 19/20 |
| 3 | — | 🟢 23.13 | 0/320 | 21.32 | 19/20 |
| 2 | — | 🟢 21.87 | 0/449 | 20.68 | 27/28 |
| 7 | — | 🟢 20.40 | 0/320 | 19.35 | 19/20 |
| 6 | — | 🟢 15.57 | 0/320 | 13.64 | 19/20 |
| 5 | — | 🟢 14.27 | 0/320 | 13.38 | 19/20 |

### Per-pixel figures

`filter_blend` ran 43,459 times per frame in the peak window at 79 cycles per blend.

## Column-ISR / DMA marshaling cost

```
isr_wake        1151/frame  min/avg/max 0.5/1.7/12.0 us  cpu 3.07%
isr_pack         144/frame  min/avg/max 6.2/6.8/9.5 us  cpu 1.56%
isr_dma_submit   144/frame  min/avg/max 0.7/0.9/1.1 us  cpu 0.21%
```

Read from the peak window. The ISR share is included in every scope above, since CYCCNT free-runs.

## Summary ranking

1. `db_timeline_step` — 64% of the peak window, 39.85 ms/f (`db_mesh_plot` 38.29).
2. `pov_preserve_half` — 0% of the peak window, 0.14 ms/f.
3. `canvas_clear` — 0% of the peak window, 0.08 ms/f.

README cells: peak 🟢 41.44 (10), spilled 🟢 0/3647 (0.00%).

## Caveats

- All scopes absorb ISR time (CYCCNT free-runs).
- `filter_blend` parents under whichever scope first enters it; calls ≈ blended pixels.
- Setup frame 1 is excluded from the peak and spill figures and reported separately; the worst-window figure skips the window that holds it.
- Selective-O3: shipping `profile` image at `e2f5b0a3d`; the previous report ran at `97eb0bf78`, and the delta between them is not attributed here.
- Dwell-compression knobs change how long an entry holds, not its per-frame cost.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DreamBalls`, `HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh DreamBalls profile 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` builds, flashes and captures under the device lock.
