# RingShower on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_ringshower_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh RingShower profile_o3 70 32`). Raw capture:
`build/prof/ringshower_o3.log`, captured 2026-08-26 01:36 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | RingShower 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh RingShower profile_o3 70 32` |

Image size: `FLASH: code:78,024, data:145,660, headers:8,764` / `RAM1: variables:315,104, code:50,840, padding:14,696, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 385–416
root counter cycles ÷ 600 MHz match the measured wall sum within **0.3 ppm**
(`tools/parse_profile.py build/prof/ringshower_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rsh_draw_rings` avg
0.79 ms/f, worst window 1.85 ms/f
(frames 385–416),
peak frame render 3.86 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 193–224)

```text
frame                          62.42ms  37.45Mcyc 100.0%
  rsh_draw_rings                1.28ms  769.9kcyc   2.1%
    rsh_ring_plot               1.28ms  769.4kcyc   2.1%
      filter_blend             112.7us   67.6kcyc   0.2% x1026.8 0us/c
  rsh_timeline_step            105.2us   63.1kcyc   0.2% x1.0 105us/c
  canvas_clear                  84.6us   50.8kcyc   0.1% x1.0 85us/c
  canvas_buffer_wait           60.94ms  36.57Mcyc  97.6% x1.0 60943us/c
```

Wall min/avg/max = 61.21/62.42/64.40 ms. The `rsh_draw_rings` scope averages 1.28 ms/f in this window, while the exact frame-render peak is 3.86 ms. That is 58.64 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 385–416)

```text
frame                          62.38ms  37.43Mcyc 100.0%
  rsh_draw_rings                1.85ms   1.11Mcyc   3.0%
    rsh_ring_plot               1.85ms   1.11Mcyc   3.0%
      filter_blend             159.8us   95.9kcyc   0.3% x1449.3 0us/c
  rsh_timeline_step             56.1us   33.7kcyc   0.1% x1.0 56us/c
  canvas_clear                  84.3us   50.6kcyc   0.1% x1.0 84us/c
  canvas_buffer_wait           60.40ms  36.24Mcyc  96.8% x1.0 60395us/c
```

Wall min/avg/max = 60.34/62.38/63.85 ms. The `rsh_draw_rings` scope averages 1.85 ms/f in this window, while the exact frame-render peak is 3.42 ms. That is 59.08 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 1,026.8 blended px/frame versus the 10,368-px quadrant (0.10× coverage), at 65.9 cyc/blend. `rsh_draw_rings` contributes 749.8 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.5/frame  min/avg/max 0.40/1.42/11.06 us  cpu 2.61%
isr_pack          143.9/frame  min/avg/max 5.99/6.45/8.61 us  cpu 1.48%
isr_dma_submit    143.9/frame  min/avg/max 0.74/0.93/1.01 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.30% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `rsh_draw_rings` — 1.3% of measured root time, 0.79 ms/f average.
2. `rsh_ring_plot` — 1.3% of measured root time, 0.79 ms/f average.
3. `canvas_clear` — 0.1% of measured root time, 0.08 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `rsh_ring_plot`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- No effect-local selective-O3 boundary is exposed by this counter tree; the paired selective-O3 capture is the shipping reference, while global O3 compiles every translation unit.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=RingShower`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh RingShower profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 3.98 ms versus
3.86 ms here: global O3 lowers the peak by 0.13 ms (3.2%). O3-image minus shipping-image
size deltas are **FLASH code +16,336 B** and **ITCM code
+15,136 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
