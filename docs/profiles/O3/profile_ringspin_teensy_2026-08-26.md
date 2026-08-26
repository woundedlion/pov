# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_ringspin_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh RingSpin profile_o3 70 32`). Raw capture:
`build/prof/ringspin_o3.log`, captured 2026-08-26 01:38 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingSpin profile_o3 70 32` |

Image size: `FLASH: code:67,624, data:148,000, headers:8,632` / `RAM1: variables:315,104, code:45,432, padding:20,104, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 737–768
root counter cycles ÷ 600 MHz match the measured wall sum within **1.5 ppm**
(`tools/parse_profile.py build/prof/ringspin_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rs_draw_rings` avg
29.26 ms/f, worst window 40.53 ms/f
(frames 737–768),
peak frame render 50.81 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 769–800)

```text
frame                          62.20ms  37.32Mcyc 100.0%
  rs_draw_rings                29.62ms  17.77Mcyc  47.6%
    rs_ring_scan               29.10ms  17.46Mcyc  46.8%
      filter_blend              4.16ms   2.50Mcyc   6.7% x41275.7 0us/c
  rs_timeline_step              67.4us   40.5kcyc   0.1% x1.0 67us/c
  canvas_clear                  86.3us   51.8kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           32.43ms  19.46Mcyc  52.1% x1.0 32430us/c
```

Wall min/avg/max = 33.19/62.20/89.70 ms. The `rs_draw_rings` scope averages 29.62 ms/f in this window, while the exact frame-render peak is 50.81 ms. That is 11.69 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 737–768)

```text
frame                          62.03ms  37.22Mcyc 100.0%
  rs_draw_rings                40.53ms  24.32Mcyc  65.3%
    rs_ring_scan               40.01ms  24.01Mcyc  64.5%
      filter_blend              5.89ms   3.54Mcyc   9.5% x58593.5 0us/c
  rs_timeline_step              69.4us   41.6kcyc   0.1% x1.0 69us/c
  canvas_clear                  87.0us   52.2kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           21.34ms  12.80Mcyc  34.4% x1.0 21338us/c
```

Wall min/avg/max = 42.57/62.03/79.46 ms. The `rs_draw_rings` scope averages 40.53 ms/f in this window, while the exact frame-render peak is 49.91 ms. That is 12.59 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 41,275.7 blended px/frame versus the 10,368-px quadrant (3.98× coverage), at 60.5 cyc/blend. `rs_draw_rings` contributes 430.5 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1147.4/frame  min/avg/max 0.40/1.54/11.68 us  cpu 2.84%
isr_pack          143.4/frame  min/avg/max 5.99/6.68/9.23 us  cpu 1.53%
isr_dma_submit    143.4/frame  min/avg/max 0.77/0.93/4.24 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.58% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `rs_draw_rings` — 46.9% of measured root time, 29.26 ms/f average.
2. `rs_ring_scan` — 46.1% of measured root time, 28.75 ms/f average.
3. `filter_blend` — 6.7% of measured root time, 4.18 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `rs_ring_scan`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingSpin`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh RingSpin profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 49.42 ms versus
50.81 ms here: global O3 raises the peak by 1.39 ms (2.8%). O3-image minus shipping-image
size deltas are **FLASH code +15,080 B** and **ITCM code
+13,088 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
