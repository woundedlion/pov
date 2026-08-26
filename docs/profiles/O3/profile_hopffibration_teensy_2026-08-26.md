# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_hopffibration_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh HopfFibration profile_o3 70 32`). Raw capture:
`build/prof/hopffibration_o3.log`, captured 2026-08-26 01:27 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | HopfFibration 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh HopfFibration profile_o3 70 32` |

Image size: `FLASH: code:80,440, data:146,584, headers:8,496` / `RAM1: variables:315,200, code:59,176, padding:6,360, free:143,552` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 193–224
root counter cycles ÷ 600 MHz match the measured wall sum within **2.5 ppm**
(`tools/parse_profile.py build/prof/hopffibration_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `hf_render_trails` avg
23.38 ms/f, worst window 39.62 ms/f
(frames 193–224),
peak frame render 46.52 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 193–224)

```text
frame                          62.78ms  37.67Mcyc 100.0%
  hf_render_trails             39.62ms  23.77Mcyc  63.1%
    hf_trail_raster            30.85ms  18.51Mcyc  49.1%
      filter_blend              5.38ms   3.23Mcyc   8.6% x49280.1 0us/c
    hf_trail_gate               7.95ms   4.77Mcyc  12.7% x210.0 38us/c
  hf_project_record            157.6us   94.6kcyc   0.3% x1.0 158us/c
  hf_advance_tumble              0.4us     265cyc   0.0% x1.0 0us/c
  hf_timeline_step              16.3us    9.8kcyc   0.0% x1.0 16us/c
  canvas_clear                  86.7us   52.0kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           22.89ms  13.74Mcyc  36.5% x1.0 22895us/c
```

Wall min/avg/max = 44.54/62.78/79.57 ms. The `hf_render_trails` scope averages 39.62 ms/f in this window, while the exact frame-render peak is 46.52 ms. That is 15.98 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 49,280.1 blended px/frame versus the 10,368-px quadrant (4.75× coverage), at 65.5 cyc/blend. `hf_render_trails` contributes 482.4 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1158.1/frame  min/avg/max 0.40/1.57/11.35 us  cpu 2.90%
isr_pack          144.8/frame  min/avg/max 6.08/7.00/9.35 us  cpu 1.61%
isr_dma_submit    144.8/frame  min/avg/max 0.65/0.93/1.10 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.72% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `hf_render_trails` — 37.5% of measured root time, 23.38 ms/f average.
2. `hf_trail_raster` — 26.5% of measured root time, 16.55 ms/f average.
3. `hf_trail_gate` — 9.6% of measured root time, 6.02 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `hf_trail_raster`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=HopfFibration`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh HopfFibration profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 48.50 ms versus
46.52 ms here: global O3 lowers the peak by 1.97 ms (4.1%). O3-image minus shipping-image
size deltas are **FLASH code +19,616 B** and **ITCM code
+18,272 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
