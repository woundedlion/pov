# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_raymarch_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh Raymarch profile_o3 70 32`). Raw capture:
`build/prof/raymarch_o3.log`, captured 2026-08-26 01:34 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | Raymarch 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh Raymarch profile_o3 70 32` |

Image size: `FLASH: code:118,120, data:183,224, headers:8,928` / `RAM1: variables:315,072, code:44,104, padding:21,432, free:143,680` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 513–544
root counter cycles ÷ 600 MHz match the measured wall sum within **0.8 ppm**
(`tools/parse_profile.py build/prof/raymarch_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `rm_shader_draw` avg
46.93 ms/f, worst window 49.90 ms/f
(frames 513–544),
peak frame render 58.80 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 321–352)

```text
frame                          62.38ms  37.43Mcyc 100.0%
  rm_shader_draw               48.64ms  29.18Mcyc  78.0%
    filter_blend               418.0us  250.8kcyc   0.7% x5707.3 0us/c
  rm_timeline_step             284.5us  170.7kcyc   0.5% x1.0 285us/c
  canvas_clear                  86.6us   52.0kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait           10.51ms   6.31Mcyc  16.8% x1.0 10510us/c
```

Wall min/avg/max = 54.07/62.38/74.14 ms. The `rm_shader_draw` scope averages 48.64 ms/f in this window, while the exact frame-render peak is 58.80 ms. That is 3.70 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 513–544)

```text
frame                          62.43ms  37.46Mcyc 100.0%
  rm_shader_draw               49.90ms  29.94Mcyc  79.9%
    filter_blend               420.2us  252.2kcyc   0.7% x5818.2 0us/c
  rm_timeline_step             281.9us  169.2kcyc   0.5% x1.0 282us/c
  canvas_clear                  87.1us   52.3kcyc   0.1% x1.0 87us/c
  canvas_buffer_wait            9.29ms   5.58Mcyc  14.9% x1.0 9292us/c
```

Wall min/avg/max = 59.79/62.43/65.19 ms. The `rm_shader_draw` scope averages 49.90 ms/f in this window, while the exact frame-render peak is 55.99 ms. That is 6.51 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

`filter_blend` averages 5,707.3 blended px/frame versus the 10,368-px quadrant (0.55× coverage), at 43.9 cyc/blend. `rm_shader_draw` contributes 5113.4 cycles per blended pixel in the selected window; this ratio is attribution, not an isolated microbenchmark.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1150.8/frame  min/avg/max 0.48/1.56/12.44 us  cpu 2.87%
isr_pack          143.8/frame  min/avg/max 6.00/6.68/9.01 us  cpu 1.53%
isr_dma_submit    143.8/frame  min/avg/max 0.64/0.93/1.57 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.61% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `rm_shader_draw` — 75.2% of measured root time, 46.93 ms/f average.
2. `filter_blend` — 0.6% of measured root time, 0.39 ms/f average.
3. `rm_timeline_step` — 0.5% of measured root time, 0.28 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `rm_shader_draw`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=Raymarch`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh Raymarch profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 62.22 ms versus
58.80 ms here: global O3 lowers the peak by 3.41 ms (5.5%). O3-image minus shipping-image
size deltas are **FLASH code +10,216 B** and **ITCM code
+6,688 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
