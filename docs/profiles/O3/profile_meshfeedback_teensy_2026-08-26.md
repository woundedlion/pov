# MeshFeedback on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_meshfeedback_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh MeshFeedback profile_o3 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"`). Raw capture:
`build/prof/meshfeedback_o3.log`, captured 2026-08-26 02:10 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | MeshFeedback 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 16 frames, 420 s capture, flags `-D HS_PROFILE_EPOCH_REVS=3400` |
| Reproduce | `bash tools/profile_one.sh MeshFeedback profile_o3 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"` |

Image size: `FLASH: code:150,944, data:186,388, headers:8,780` / `RAM1: variables:315,264, code:85,528, padding:12,776, free:110,720` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 977–992
root counter cycles ÷ 600 MHz match the measured wall sum within **0.7 ppm**
(`tools/parse_profile.py build/prof/meshfeedback_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `mf_feedback_flush` avg
40.35 ms/f, worst window 47.02 ms/f
(frames 977–992),
peak frame render 58.32 ms, spilled 0/6688
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: Marker-defined presets own their following transitions; the first section is the worst exact-render window and the second is the lowest-cost marked hold.

### Worst exact-render regime (window frames 1201–1216)

```text
frame                          62.04ms  37.22Mcyc 100.0%
  mf_mesh_draw                  5.44ms   3.27Mcyc   8.8%
    filter_blend                1.33ms  797.1kcyc   2.1% x8985.1 0us/c
  mf_feedback_flush            43.52ms  26.11Mcyc  70.2%
    feedback_composite         37.10ms  22.26Mcyc  59.8% x1.0 37103us/c
    feedback_populate           6.40ms   3.84Mcyc  10.3% x1.0 6404us/c
    feedback_litscan             0.2us     120cyc   0.0% x1.0 0us/c
  mf_timeline_step              38.1us   22.9kcyc   0.1% x1.0 38us/c
  mf_apply_params               10.5us    6.3kcyc   0.0% x1.0 10us/c
  canvas_clear                 398.4us  239.1kcyc   0.6% x1.0 398us/c
  canvas_buffer_wait           12.62ms   7.57Mcyc  20.3% x1.0 12622us/c
```

Wall min/avg/max = 54.06/62.04/67.00 ms. The `mf_feedback_flush` scope averages 43.52 ms/f in this window, while the exact frame-render peak is 58.32 ms. That is 4.18 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 5617–5632)

```text
frame                          62.57ms  37.54Mcyc 100.0%
  mf_mesh_draw                  8.80ms   5.28Mcyc  14.1%
    filter_blend                2.13ms   1.28Mcyc   3.4% x14448.4 0us/c
  mf_feedback_flush            32.60ms  19.56Mcyc  52.1%
    feedback_composite         32.58ms  19.55Mcyc  52.1% x1.0 32580us/c
    feedback_populate            0.0us      21cyc   0.0% x1.0 0us/c
    feedback_litscan             0.2us     115cyc   0.0% x1.0 0us/c
  mf_timeline_step              30.1us   18.1kcyc   0.0% x1.0 30us/c
  mf_apply_params                0.1us      71cyc   0.0% x1.0 0us/c
  canvas_clear                 392.3us  235.4kcyc   0.6% x1.0 392us/c
  canvas_buffer_wait           20.74ms  12.44Mcyc  33.1% x1.0 20736us/c
```

Wall min/avg/max = 61.93/62.57/63.15 ms. The `mf_feedback_flush` scope averages 32.60 ms/f in this window, while the exact frame-render peak is 43.14 ms. That is 19.36 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Per-preset table

| Preset/shape | Geometry | Hold windows | Blended px/f | `mf_feedback_flush` ms | Peak render ms | fps |
|---|---|---:|---:|---:|---:|---:|
| `5` | — | 30/30 | 9,235.1 | 47.02 | 58.29 | 16.0 |
| `3` | — | 45/45 | 9,612.9 | 45.12 | 54.42 | 16.0 |
| `9` | — | 30/30 | 10,537.0 | 44.88 | 57.40 | 16.0 |
| `10` | — | 30/30 | 8,211.3 | 44.56 | 54.56 | 16.0 |
| `2` | — | 46/46 | 5,294.0 | 44.11 | 49.99 | 16.0 |
| `4` | — | 41/41 | 3,711.1 | 43.60 | 50.79 | 16.0 |
| `8` | — | 30/30 | 10,228.1 | 42.58 | 48.54 | 16.0 |
| `6` | — | 31/31 | 8,918.3 | 41.86 | 58.32 | 16.0 |
| `7` | — | 30/30 | 8,258.3 | 40.75 | 50.69 | 16.0 |
| `1` | — | 30/30 | 6,700.6 | 39.42 | 47.72 | 16.0 |
| `0` | — | 15/15 | 8,113.9 | 39.12 | 45.72 | 16.0 |
| `12` | — | 30/30 | 12,466.1 | 38.67 | 53.69 | 16.0 |
| `11` | — | 30/30 | 12,855.0 | 38.13 | 54.15 | 16.0 |

403 advance markers visit 12/12 presets; wrap-to-0 is confirmed.

### Per-pixel figures

`filter_blend` averages 8,985.1 blended px/frame versus the 10,368-px quadrant (0.87× coverage), at 88.7 cyc/blend. `mf_feedback_flush` contributes 2906.4 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1144.7/frame  min/avg/max 0.41/1.61/11.12 us  cpu 2.96%
isr_pack          143.1/frame  min/avg/max 6.10/7.00/9.21 us  cpu 1.61%
isr_dma_submit    143.1/frame  min/avg/max 0.78/0.94/1.02 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.78% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `mf_feedback_flush` — 64.6% of measured root time, 40.35 ms/f average.
2. `feedback_composite` — 56.3% of measured root time, 35.16 ms/f average.
3. `feedback_populate` — 8.3% of measured root time, 5.17 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `mf_mesh_draw`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- The capture uses `-D HS_PROFILE_EPOCH_REVS=3400` to keep the full cycle inside one effect epoch; it changes dwell/epoch length, not per-frame render cost.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `20ca3cb48892795a4575dd9d16d31a699f79df75`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=MeshFeedback`,
`HS_PROFILE_WINDOW=16`; `bash tools/profile_one.sh MeshFeedback profile_o3 420 16 "-D HS_PROFILE_EPOCH_REVS=3400"` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 58.30 ms versus
58.32 ms here: global O3 raises the peak by 0.02 ms (0.0%). O3-image minus shipping-image
size deltas are **FLASH code +34,144 B** and **ITCM code
+21,776 B**. Spill fractions are compared rather than raw counts:
shipping 0/6688 (0%)
versus O3 0/6688 (0%).
