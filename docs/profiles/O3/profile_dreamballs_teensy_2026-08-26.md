# DreamBalls on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_dreamballs_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh DreamBalls profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"`). Raw capture:
`build/prof/dreamballs_o3.log`, captured 2026-08-26 02:19 local. This replaces `profile_dreamballs_teensy_2026-08-09.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | DreamBalls 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 230 s capture, flags `-D HS_PROFILE_EPOCH_REVS=2000` |
| Reproduce | `bash tools/profile_one.sh DreamBalls profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` |

Image size: `FLASH: code:137,200, data:187,004, headers:8,596` / `RAM1: variables:315,360, code:62,360, padding:3,176, free:143,392` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 2849–2864
root counter cycles ÷ 600 MHz match the measured wall sum within **0.1 ppm**
(`tools/parse_profile.py build/prof/dreamballs_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `db_timeline_step` avg
20.53 ms/f, worst window 33.02 ms/f
(frames 2849–2864),
peak frame render 34.32 ms, spilled 0/3648
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 2593–2608)

```text
frame                          62.46ms  37.48Mcyc 100.0%
  db_timeline_step             32.59ms  19.55Mcyc  52.2%
    db_draw                    32.54ms  19.52Mcyc  52.1%
      db_draw_scene            32.54ms  19.52Mcyc  52.1%
        db_mesh_plot           31.05ms  18.63Mcyc  49.7%
          filter_blend          4.14ms   2.48Mcyc   6.6% x42711.4 0us/c
        db_orient              271.0us  162.6kcyc   0.4% x6.0 45us/c
        db_displace             1.02ms  611.5kcyc   1.6% x6.0 170us/c
  canvas_clear                  89.3us   53.6kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait           29.78ms  17.87Mcyc  47.7% x1.0 29780us/c
```

Wall min/avg/max = 60.79/62.46/64.79 ms. The `db_timeline_step` scope averages 32.59 ms/f in this window, while the exact frame-render peak is 34.32 ms. That is 28.18 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 1585–1600)

```text
frame                          61.82ms  37.09Mcyc 100.0%
  db_timeline_step             10.21ms   6.13Mcyc  16.5%
    db_draw                    10.16ms   6.10Mcyc  16.4%
      db_draw_scene            10.16ms   6.10Mcyc  16.4%
        db_mesh_plot            9.82ms   5.89Mcyc  15.9%
          filter_blend          1.35ms  809.5kcyc   2.2% x13782.0 0us/c
        db_orient               57.9us   34.8kcyc   0.1% x3.8 15us/c
        db_displace            211.6us  127.0kcyc   0.3% x3.8 56us/c
  canvas_clear                  88.9us   53.4kcyc   0.1% x1.0 89us/c
  canvas_buffer_wait           51.52ms  30.91Mcyc  83.3% x1.0 51519us/c
```

Wall min/avg/max = 51.66/61.82/63.35 ms. The `db_timeline_step` scope averages 10.21 ms/f in this window, while the exact frame-render peak is 11.41 ms. That is 51.09 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `db_timeline_step` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `9` | — | 20/20 | 42,791.1 | 33.02 | 34.32 | 16.0 |
| `1` | — | 20/20 | 42,575.4 | 31.02 | 32.19 | 16.0 |
| `0` | — | 20/20 | 42,443.4 | 30.87 | 31.75 | 16.0 |
| `10` | — | 20/20 | 32,851.2 | 24.62 | 26.50 | 16.0 |
| `8` | — | 20/20 | 32,679.8 | 23.93 | 25.66 | 16.0 |
| `4` | — | 20/20 | 25,855.2 | 18.86 | 21.12 | 16.0 |
| `3` | — | 20/20 | 23,123.8 | 17.30 | 18.91 | 16.0 |
| `2` | — | 28/28 | 23,110.1 | 17.15 | 18.07 | 16.0 |
| `7` | — | 20/20 | 21,496.6 | 16.25 | 16.79 | 16.0 |
| `6` | — | 20/20 | 14,432.3 | 11.16 | 12.60 | 16.0 |
| `5` | — | 20/20 | 14,755.8 | 11.00 | 11.67 | 16.0 |

208 advance markers visit 10/10 presets; wrap-to-0 is confirmed.

### Per-pixel figures

`filter_blend` averages 42,711.4 blended px/frame versus the 10,368-px quadrant (4.12× coverage), at 58.1 cyc/blend. `db_timeline_step` contributes 457.8 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.3/frame  min/avg/max 0.42/1.55/11.32 us  cpu 2.86%
isr_pack          144.1/frame  min/avg/max 6.02/6.76/9.52 us  cpu 1.55%
isr_dma_submit    144.1/frame  min/avg/max 0.84/0.95/1.03 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.62% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `db_timeline_step` — 32.9% of measured root time, 20.53 ms/f average.
2. `db_draw` — 32.8% of measured root time, 20.49 ms/f average.
3. `db_draw_scene` — 32.8% of measured root time, 20.48 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `db_mesh_plot`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=2000` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=DreamBalls`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh DreamBalls profile_o3 230 16 "-D HS_PROFILE_EPOCH_REVS=2000"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 38.73 ms versus
34.32 ms here: global O3 lowers the peak by 4.42 ms (11.4%). O3-image minus shipping-image
size deltas are **FLASH code +26,896 B** and **ITCM code
+12,352 B**. Spill fractions are compared rather than raw counts:
shipping 0/3648 (0%)
versus O3 0/3648 (0%).
