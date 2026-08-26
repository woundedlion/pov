# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_petalflow_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh PetalFlow profile_o3 70 32`). Raw capture:
`build/prof/petalflow_o3.log`, captured 2026-08-26 01:32 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | PetalFlow 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh PetalFlow profile_o3 70 32` |

Image size: `FLASH: code:78,096, data:147,032, headers:8,344` / `RAM1: variables:315,136, code:56,120, padding:9,416, free:143,616` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 961–992
root counter cycles ÷ 600 MHz match the measured wall sum within **0.8 ppm**
(`tools/parse_profile.py build/prof/petalflow_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `pf_draw_rings` avg
9.70 ms/f, worst window 10.07 ms/f
(frames 961–992),
peak frame render 10.77 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 289–320)

```text
frame                          62.48ms  37.49Mcyc 100.0%
  pf_draw_rings                 9.84ms   5.91Mcyc  15.8%
    pf_ring_scan                8.70ms   5.22Mcyc  13.9%
      filter_blend             728.4us  437.0kcyc   1.2% x6428.8 0us/c
    pf_ring_build               1.05ms  632.3kcyc   1.7% x22.7 47us/c
  pf_timeline_step              19.1us   11.5kcyc   0.0% x1.0 19us/c
  canvas_clear                  84.8us   50.9kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           52.53ms  31.52Mcyc  84.1% x1.0 52533us/c
```

Wall min/avg/max = 61.68/62.48/63.53 ms. The `pf_draw_rings` scope averages 9.84 ms/f in this window, while the exact frame-render peak is 10.77 ms. That is 51.73 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 6,428.8 blended px/frame versus the 10,368-px quadrant (0.62× coverage), at 68.0 cyc/blend. `pf_draw_rings` contributes 918.7 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1152.6/frame  min/avg/max 0.40/1.42/21.48 us  cpu 2.61%
isr_pack          144.1/frame  min/avg/max 5.99/6.45/8.79 us  cpu 1.48%
isr_dma_submit    144.1/frame  min/avg/max 0.62/0.93/11.92 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.30% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `pf_draw_rings` — 15.6% of measured root time, 9.70 ms/f average.
2. `pf_ring_scan` — 13.7% of measured root time, 8.53 ms/f average.
3. `pf_ring_build` — 1.7% of measured root time, 1.07 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `pf_ring_scan`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses selective-O3 mesh/scan/filter hot paths identified by the counter tree; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=PetalFlow`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh PetalFlow profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 11.85 ms versus
10.77 ms here: global O3 lowers the peak by 1.08 ms (9.1%). O3-image minus shipping-image
size deltas are **FLASH code +22,272 B** and **ITCM code
+21,008 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
