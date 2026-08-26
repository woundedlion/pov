# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_shapeshifter_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh ShapeShifter profile_o3 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"`). Raw capture:
`build/prof/shapeshifter_o3.log`, captured 2026-08-26 03:11 local. This replaces `profile_shapeshifter_teensy_2026-08-08.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | ShapeShifter 288×144, single-entry playlist, tip `20ca3cb48892795a4575dd9d16d31a699f79df75` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 155 s capture, flags `-D HS_PROFILE_EPOCH_REVS=1600` |
| Reproduce | `bash tools/profile_one.sh ShapeShifter profile_o3 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` |

Image size: `FLASH: code:104,248, data:150,500, headers:8,420` / `RAM1: variables:315,104, code:72,184, padding:26,120, free:110,880` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2289–2304
root counter cycles ÷ 600 MHz match the measured wall sum within **2.3 ppm**
(`tools/parse_profile.py build/prof/shapeshifter_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `ss_draw_all` avg
23.47 ms/f, worst window 51.13 ms/f
(frames 2289–2304),
peak frame render 59.80 ms, spilled 0/2448
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `ss_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 2289–2304)

```text
frame                          62.27ms  37.36Mcyc 100.0%
  ss_draw_all                  51.13ms  30.68Mcyc  82.1%
    ss_plot_dispatch           50.93ms  30.56Mcyc  81.8% x137.9 369us/c
  ss_timeline_step              42.1us   25.3kcyc   0.1% x1.0 42us/c
  ss_buffer_wait               11.09ms   6.66Mcyc  17.8%
    canvas_clear                88.9us   53.3kcyc   0.1% x1.0 89us/c
    canvas_buffer_wait         11.00ms   6.60Mcyc  17.7% x1.0 11003us/c
```

Wall min/avg/max = 50.49/62.27/76.32 ms. The `ss_draw_all` scope averages 51.13 ms/f in this window, while the exact frame-render peak is 59.80 ms. That is 2.70 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 545–560)

```text
frame                          62.50ms  37.50Mcyc 100.0%
  ss_draw_all                   4.90ms   2.94Mcyc   7.8%
    ss_plot_dispatch            4.80ms   2.88Mcyc   7.7% x24.1 199us/c
  ss_timeline_step              35.1us   21.1kcyc   0.1% x1.0 35us/c
  ss_buffer_wait               57.56ms  34.53Mcyc  92.1%
    canvas_clear                84.5us   50.7kcyc   0.1% x1.0 84us/c
    canvas_buffer_wait         57.47ms  34.48Mcyc  92.0% x1.0 57472us/c
```

Wall min/avg/max = 61.24/62.50/63.73 ms. The `ss_draw_all` scope averages 4.90 ms/f in this window, while the exact frame-render peak is 6.10 ms. That is 56.40 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `ss_draw_all` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `1` | — | 15/15 | — | 51.13 | 59.80 | 16.0 |
| `0` | — | 15/15 | — | 44.65 | 51.41 | 16.0 |
| `9` | — | 15/15 | — | 35.52 | 36.05 | 16.0 |
| `4` | — | 15/15 | — | 34.77 | 35.69 | 16.0 |
| `8` | — | 15/15 | — | 24.27 | 25.64 | 16.0 |
| `7` | — | 15/15 | — | 20.12 | 22.63 | 16.0 |
| `6` | — | 15/15 | — | 17.15 | 18.23 | 16.0 |
| `2` | — | 18/18 | — | 12.90 | 15.26 | 16.0 |
| `5` | — | 15/15 | — | 9.86 | 11.24 | 16.0 |
| `3` | — | 15/15 | — | 7.63 | 8.71 | 16.0 |

138 advance markers visit 9/9 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1148.8/frame  min/avg/max 0.34/1.58/10.85 us  cpu 2.91%
isr_pack          143.6/frame  min/avg/max 6.10/6.89/9.19 us  cpu 1.58%
isr_dma_submit    143.6/frame  min/avg/max 0.63/0.93/1.00 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.70% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `ss_draw_all` — 37.6% of measured root time, 23.47 ms/f average.
2. `ss_plot_dispatch` — 37.4% of measured root time, 23.36 ms/f average.
3. `canvas_clear` — 0.1% of measured root time, 0.09 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- No effect-local selective-O3 boundary is exposed by this counter tree; the paired selective-O3 capture is the shipping reference, while global O3 compiles every translation unit.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=1600` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `20ca3cb48892795a4575dd9d16d31a699f79df75` and paired shipping source
  `20ca3cb48892795a4575dd9d16d31a699f79df75`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=ShapeShifter`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh ShapeShifter profile_o3 155 16 "-D HS_PROFILE_EPOCH_REVS=1600"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 58.13 ms versus
59.80 ms here: global O3 raises the peak by 1.66 ms (2.9%). O3-image minus shipping-image
size deltas are **FLASH code +27,272 B** and **ITCM code
+23,440 B**. Spill fractions are compared rather than raw counts:
shipping 0/2448 (0%)
versus O3 0/2448 (0%).
