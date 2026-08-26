# Comets on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_comets_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh Comets profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"`). Raw capture:
`build/prof/comets_o3.log`, captured 2026-08-26 02:02 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Comets 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 260 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2400` |
| Reproduce | `bash tools/profile_one.sh Comets profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` |

Image size: `FLASH: code:80,096, data:148,888, headers:8,584` / `RAM1: variables:315,104, code:48,136, padding:17,400, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 1921–1936
root counter cycles ÷ 600 MHz match the measured wall sum within **0.4 ppm**
(`tools/parse_profile.py build/prof/comets_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `cm_draw_trail` avg
12.51 ms/f, worst window 22.61 ms/f
(frames 1921–1936),
peak frame render 28.71 ms, spilled 0/4128
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 3521–3536)

```text
frame                          62.72ms  37.63Mcyc 100.0%
  cm_draw_trail                19.62ms  11.77Mcyc  31.3%
    filter_blend                1.08ms  647.0kcyc   1.7% x11981.4 0us/c
  cm_wipe_rebake                2.17ms   1.30Mcyc   3.5% x1.0 2170us/c
  cm_timeline_step              97.2us   58.3kcyc   0.2% x1.0 97us/c
  canvas_clear                  84.8us   50.9kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           40.75ms  24.45Mcyc  65.0% x1.0 40748us/c
```

Wall min/avg/max = 50.91/62.72/74.70 ms. The `cm_draw_trail` scope averages 19.62 ms/f in this window, while the exact frame-render peak is 28.71 ms. That is 33.79 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 3985–4000)

```text
frame                          62.47ms  37.48Mcyc 100.0%
  cm_draw_trail                 1.75ms   1.05Mcyc   2.8%
    filter_blend                36.2us   21.7kcyc   0.1% x393.8 0us/c
  cm_wipe_rebake                 0.0us      29cyc   0.0% x1.0 0us/c
  cm_timeline_step              68.6us   41.2kcyc   0.1% x1.0 69us/c
  canvas_clear                  84.4us   50.7kcyc   0.1% x1.0 84us/c
  canvas_buffer_wait           60.56ms  36.34Mcyc  96.9% x1.0 60559us/c
```

Wall min/avg/max = 61.52/62.47/63.46 ms. The `cm_draw_trail` scope averages 1.75 ms/f in this window, while the exact frame-render peak is 2.44 ms. That is 60.06 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `cm_draw_trail` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `1` | — | 20/20 | 14,629.8 | 22.61 | 28.67 | 16.0 |
| `11` | — | 20/20 | 11,981.4 | 19.62 | 28.71 | 16.0 |
| `10` | — | 20/20 | 11,026.6 | 19.57 | 23.83 | 16.0 |
| `2` | — | 28/28 | 12,295.8 | 19.54 | 23.25 | 16.0 |
| `12` | — | 20/20 | 11,537.2 | 19.42 | 27.38 | 16.0 |
| `8` | — | 20/20 | 11,040.8 | 18.97 | 25.02 | 16.0 |
| `4` | — | 20/20 | 12,154.9 | 18.75 | 25.08 | 16.0 |
| `5` | — | 20/20 | 10,652.2 | 18.36 | 26.88 | 16.0 |
| `9` | — | 20/20 | 8,604.9 | 17.05 | 22.95 | 16.0 |
| `3` | — | 20/20 | 9,209.1 | 14.82 | 17.76 | 16.0 |
| `6` | — | 20/20 | 5,886.4 | 12.60 | 17.41 | 16.0 |
| `0` | — | 10/10 | 8,353.8 | 11.68 | 13.83 | 16.0 |
| `7` | — | 20/20 | 6,760.9 | 11.11 | 14.46 | 16.0 |

248 advance markers visit 12/12 presets; wrap-to-0 is confirmed.

### Per-pixel figures

`filter_blend` averages 11,981.4 blended px/frame versus the 10,368-px quadrant (1.16× coverage), at 54.0 cyc/blend. `cm_draw_trail` contributes 982.5 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1156.9/frame  min/avg/max 0.54/1.55/14.72 us  cpu 2.86%
isr_pack          144.6/frame  min/avg/max 5.99/6.60/8.95 us  cpu 1.52%
isr_dma_submit    144.6/frame  min/avg/max 0.69/0.94/7.38 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.59% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `cm_draw_trail` — 20.0% of measured root time, 12.51 ms/f average.
2. `cm_wipe_rebake` — 1.0% of measured root time, 0.63 ms/f average.
3. `filter_blend` — 0.9% of measured root time, 0.59 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `cm_draw_trail`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- No effect-local selective-O3 boundary is exposed by this counter tree; the paired selective-O3 capture is the shipping reference, while global O3 compiles every translation unit.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=2400` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Comets`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh Comets profile_o3 260 16 "-D HS_PROFILE_EPOCH_REVS=2400"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 33.91 ms versus
28.71 ms here: global O3 lowers the peak by 5.20 ms (15.3%). O3-image minus shipping-image
size deltas are **FLASH code +15,432 B** and **ITCM code
+12,944 B**. Spill fractions are compared rather than raw counts:
shipping 0/4128 (0%)
versus O3 0/4128 (0%).
