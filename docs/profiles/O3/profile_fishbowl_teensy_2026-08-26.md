# Fishbowl on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_fishbowl_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh Fishbowl profile_o3 70 32`). Raw capture:
`build/prof/fishbowl_o3.log`, captured 2026-08-26 01:18 local. This replaces `profile_fishbowl_teensy_2026-08-02.md`.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Fishbowl 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Fishbowl profile_o3 70 32` |

Image size: `FLASH: code:98,376, data:149,460, headers:9,188` / `RAM1: variables:315,136, code:72,120, padding:26,184, free:110,848` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 481–512
root counter cycles ÷ 600 MHz match the measured wall sum within **2.2 ppm**
(`tools/parse_profile.py build/prof/fishbowl_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fish_build_vertices` avg
11.50 ms/f, worst window 12.16 ms/f
(frames 481–512),
peak frame render 21.29 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 961–992)

```text
frame                          62.44ms  37.47Mcyc 100.0%
  fish_multiline_draw           7.70ms   4.62Mcyc  12.3%
    filter_blend               718.8us  431.3kcyc   1.2% x7291.4 0us/c
  fish_build_vertices          12.12ms   7.27Mcyc  19.4% x1.0 12124us/c
  fish_noise_prepare             0.2us     160cyc   0.0% x1.0 0us/c
  fish_timeline_step            96.2us   57.7kcyc   0.2% x1.0 96us/c
  canvas_clear                  87.8us   52.7kcyc   0.1% x1.0 88us/c
  canvas_buffer_wait           42.43ms  25.46Mcyc  68.0% x1.0 42434us/c
```

Wall min/avg/max = 60.44/62.44/64.46 ms. The `fish_build_vertices` scope averages 12.12 ms/f in this window, while the exact frame-render peak is 21.29 ms. That is 41.21 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 481–512)

```text
frame                          62.50ms  37.50Mcyc 100.0%
  fish_multiline_draw           7.46ms   4.48Mcyc  11.9%
    filter_blend               671.2us  402.7kcyc   1.1% x6760.5 0us/c
  fish_build_vertices          12.16ms   7.30Mcyc  19.5% x1.0 12160us/c
  fish_noise_prepare             0.2us     161cyc   0.0% x1.0 0us/c
  fish_timeline_step            96.6us   58.0kcyc   0.2% x1.0 97us/c
  canvas_clear                  87.5us   52.5kcyc   0.1% x1.0 88us/c
  canvas_buffer_wait           42.69ms  25.62Mcyc  68.3% x1.0 42692us/c
```

Wall min/avg/max = 61.11/62.50/63.89 ms. The `fish_build_vertices` scope averages 12.16 ms/f in this window, while the exact frame-render peak is 20.43 ms. That is 42.07 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 7,291.4 blended px/frame versus the 10,368-px quadrant (0.70× coverage), at 59.1 cyc/blend. `fish_build_vertices` contributes 997.6 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.7/frame  min/avg/max 0.40/1.48/11.50 us  cpu 2.73%
isr_pack          143.9/frame  min/avg/max 5.99/6.69/8.79 us  cpu 1.54%
isr_dma_submit    143.9/frame  min/avg/max 0.67/0.93/1.12 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.48% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fish_build_vertices` — 18.4% of measured root time, 11.50 ms/f average.
2. `fish_multiline_draw` — 11.2% of measured root time, 7.01 ms/f average.
3. `filter_blend` — 1.0% of measured root time, 0.64 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `fish_multiline_draw`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- No effect-local selective-O3 boundary is exposed by this counter tree; the paired selective-O3 capture is the shipping reference, while global O3 compiles every translation unit.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Fishbowl`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Fishbowl profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 23.22 ms versus
21.29 ms here: global O3 lowers the peak by 1.93 ms (8.3%). O3-image minus shipping-image
size deltas are **FLASH code +27,856 B** and **ITCM code
+26,400 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
