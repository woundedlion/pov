# SphericalHarmonics on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_sphericalharmonics_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh SphericalHarmonics profile_o3 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"`). Raw capture:
`build/prof/sphericalharmonics_o3.log`, captured 2026-08-26 01:57 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | SphericalHarmonics 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 220 s capture, flags `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048` |
| Reproduce | `bash tools/profile_one.sh SphericalHarmonics profile_o3 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` |

Image size: `FLASH: code:47,672, data:146,580, headers:8,500` / `RAM1: variables:315,040, code:27,064, padding:5,704, free:176,480` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2449–2464
root counter cycles ÷ 600 MHz match the measured wall sum within **1.3 ppm**
(`tools/parse_profile.py build/prof/sphericalharmonics_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `sh_rasterize` avg
8.69 ms/f, worst window 11.53 ms/f
(frames 2449–2464),
peak frame render 11.70 ms, spilled 0/3488
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 913–928)

```text
frame                          62.45ms  37.47Mcyc 100.0%
  sh_rasterize                 11.52ms   6.91Mcyc  18.4% x1.0 11515us/c
  sh_timeline_step              19.5us   11.7kcyc   0.0% x1.0 20us/c
  canvas_clear                  84.6us   50.8kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           50.83ms  30.50Mcyc  81.4% x1.0 50826us/c
```

Wall min/avg/max = 62.20/62.45/62.58 ms. The `sh_rasterize` scope averages 11.52 ms/f in this window, while the exact frame-render peak is 11.70 ms. That is 50.80 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 2801–2816)

```text
frame                          62.44ms  37.46Mcyc 100.0%
  sh_rasterize                  6.76ms   4.05Mcyc  10.8% x1.0 6756us/c
  sh_timeline_step              21.4us   12.8kcyc   0.0% x1.0 21us/c
  canvas_clear                  84.3us   50.6kcyc   0.1% x1.0 84us/c
  canvas_buffer_wait           55.58ms  33.35Mcyc  89.0% x1.0 55576us/c
```

Wall min/avg/max = 62.26/62.44/62.62 ms. The `sh_rasterize` scope averages 6.76 ms/f in this window, while the exact frame-render peak is 6.98 ms. That is 55.52 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `sh_rasterize` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `21` | — | 8/8 | — | 11.53 | 11.70 | 16.0 |
| `20` | — | 8/8 | — | 11.32 | 11.47 | 16.0 |
| `22` | — | 8/8 | — | 10.46 | 10.64 | 16.0 |
| `19` | — | 8/8 | — | 10.46 | 10.65 | 16.0 |
| `13` | — | 10/10 | — | 10.05 | 10.22 | 16.0 |
| `12` | — | 12/12 | — | 9.84 | 10.02 | 16.0 |
| `18` | — | 8/8 | — | 9.29 | 9.45 | 16.0 |
| `23` | — | 8/8 | — | 9.23 | 9.38 | 16.0 |
| `14` | — | 8/8 | — | 8.87 | 9.04 | 16.0 |
| `11` | — | 12/12 | — | 8.86 | 9.05 | 16.0 |
| `17` | — | 8/8 | — | 8.47 | 8.60 | 16.0 |
| `24` | — | 8/8 | — | 8.46 | 8.62 | 16.0 |
| `7` | — | 12/12 | — | 8.43 | 8.61 | 16.0 |
| `6` | — | 8/8 | — | 8.25 | 8.43 | 16.0 |
| `16` | — | 8/8 | — | 8.21 | 8.37 | 16.0 |
| `10` | — | 12/12 | — | 8.17 | 8.30 | 16.0 |
| `15` | — | 8/8 | — | 8.10 | 8.25 | 16.0 |
| `9` | — | 12/12 | — | 7.88 | 8.04 | 16.0 |
| `1` | — | 8/8 | — | 7.78 | 7.92 | 16.0 |
| `8` | — | 12/12 | — | 7.75 | 7.89 | 16.0 |
| `5` | — | 8/8 | — | 7.74 | 7.88 | 16.0 |
| `4` | — | 8/8 | — | 7.47 | 7.63 | 16.0 |
| `3` | — | 8/8 | — | 7.31 | 7.49 | 16.0 |
| `2` | — | 8/8 | — | 7.00 | 7.17 | 16.0 |

218 advance markers visit 24/24 presets; wrap-to-0 is confirmed.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.1/frame  min/avg/max 0.40/1.49/10.61 us  cpu 2.74%
isr_pack          144.0/frame  min/avg/max 5.99/6.69/8.94 us  cpu 1.54%
isr_dma_submit    144.0/frame  min/avg/max 0.59/0.93/1.06 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.49% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `sh_rasterize` — 13.9% of measured root time, 8.69 ms/f average.
2. `canvas_clear` — 0.1% of measured root time, 0.08 ms/f average.
3. `sh_timeline_step` — 0.0% of measured root time, 0.02 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=SphericalHarmonics`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh SphericalHarmonics profile_o3 220 16 "-D HS_PROFILE_ORDERED_CYCLE -D HS_PROFILE_EPOCH_REVS=2048"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 12.64 ms versus
11.70 ms here: global O3 lowers the peak by 0.95 ms (7.5%). O3-image minus shipping-image
size deltas are **FLASH code +8,016 B** and **ITCM code
+6,400 B**. Spill fractions are compared rather than raw counts:
shipping 0/3488 (0%)
versus O3 0/3488 (0%).
