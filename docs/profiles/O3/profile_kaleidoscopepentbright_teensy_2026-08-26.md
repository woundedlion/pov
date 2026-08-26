# KaleidoscopePentBright on-device profile — Teensy 4.0, segmented mode (2026-08-26, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_kaleidoscopepentbright_teensy_2026-08-26.md).
Point-in-time snapshot (regenerate with `bash tools/profile_one.sh KaleidoscopePentBright profile_o3 70 32`). Raw capture:
`build/prof/kaleidoscopepentbright_o3.log`, captured 2026-08-26 02:46 local. This is the first archived global-O3 report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM4 |
| Image | `profile_o3` env: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopePentBright 288×144, single-entry playlist, tip `63268c3768ae91c8e2b47acd728db641570179eb` |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture |
| Reproduce | `bash tools/profile_one.sh KaleidoscopePentBright profile_o3 70 32` |

Image size: `FLASH: code:83,776, data:152,136, headers:8,824` / `RAM1: variables:315,104, code:37,848, padding:27,688, free:143,648` / `RAM2: variables:520,064, free:4,224`.

Exactness cross-check: window frames 417–448
root counter cycles ÷ 600 MHz match the measured wall sum within **3.1 ppm**
(`tools/parse_profile.py build/prof/kaleidoscopepentbright_o3.log validate`).

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw` avg
24.76 ms/f, worst window 26.40 ms/f
(frames 417–448),
peak frame render 29.53 ms, spilled 0/1088
frames (0%).

A display window is 62.5 ms; the effect renders one quadrant ≈ 10,368 px.
Every measured frame holds 16 fps. The `canvas_buffer_wait` scope is round-up idle to the next
display flip, by design.

## Phase-by-phase readout

Phase schedule: No discrete preset markers are present; the selected windows bound the observed continuous regime.

### Worst exact-render regime (window frames 865–896)

```text
frame                          62.42ms  37.45Mcyc 100.0%
  fx_shader_draw               25.73ms  15.44Mcyc  41.2% x1.0 25732us/c
  fx_prepare_frame             695.1us  417.1kcyc   1.1% x1.0 695us/c
  fx_advance                    2.00ms   1.20Mcyc   3.2% x1.0 1995us/c
  fx_timeline_step              44.4us   26.6kcyc   0.1% x1.0 44us/c
  canvas_clear                  86.1us   51.6kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           33.86ms  20.31Mcyc  54.2% x1.0 33857us/c
```

Wall min/avg/max = 60.47/62.42/64.48 ms. The `fx_shader_draw` scope averages 25.73 ms/f in this window, while the exact frame-render peak is 29.53 ms. That is 32.97 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.
### Comparison regime (window frames 417–448)

```text
frame                          62.41ms  37.45Mcyc 100.0%
  fx_shader_draw               26.40ms  15.84Mcyc  42.3% x1.0 26401us/c
  fx_prepare_frame             664.3us  398.6kcyc   1.1% x1.0 664us/c
  fx_advance                    1.93ms   1.16Mcyc   3.1% x1.0 1932us/c
  fx_timeline_step              31.4us   18.9kcyc   0.1% x1.0 31us/c
  canvas_clear                  86.0us   51.6kcyc   0.1% x1.0 86us/c
  canvas_buffer_wait           33.29ms  19.98Mcyc  53.3% x1.0 33294us/c
```

Wall min/avg/max = 61.83/62.41/63.03 ms. The `fx_shader_draw` scope averages 26.40 ms/f in this window, while the exact frame-render peak is 29.40 ms. That is 33.10 ms of margin against one 62.5 ms display interval; display-sync wait is idle rather than render work.

### Per-pixel figures

This capture has no `filter_blend` calls in the selected worst window, so blended-pixel coverage, cycles/blend, and scan cycles per blended pixel are not defined. The dominant-scope timing above remains the device-level per-frame measurement.

## Column-ISR / DMA marshaling cost

```text
isr_wake         1151.8/frame  min/avg/max 0.39/1.53/14.77 us  cpu 2.82%
isr_pack          144.0/frame  min/avg/max 6.06/6.96/9.76 us  cpu 1.60%
isr_dma_submit    144.0/frame  min/avg/max 0.67/0.93/4.76 us  cpu 0.21%
```

- Packing and submit are CPU marshaling costs; their relative totals are reported above without treating asynchronous LED transfer as render work.
- LED wire transfer runs asynchronously; display synchronization remains isolated in the effect's `*_buffer_wait` scope.
- The selected window's ISR counters total 4.63% CPU. The worst render requires 0.0% speedup to fit one 62.5 ms interval (zero means it already fits).

## Summary ranking

1. `fx_shader_draw` — 39.7% of measured root time, 24.76 ms/f average.
2. `fx_advance` — 3.1% of measured root time, 1.96 ms/f average.
3. `fx_prepare_frame` — 1.2% of measured root time, 0.75 ms/f average.

No matched WASM/native datum is asserted here; the paired on-device shipping
capture is the comparison baseline.

## Caveats

- All scopes absorb ISR time because `CYCCNT` free-runs.
- `filter_blend` parents under `the first entering scope`; its subtree is hidden in
  windows where that parent has zero calls, and calls approximate blended pixels.
- The paired shipping image traverses the selective-O3 `Scan::Shader` and shared math/color hot paths used by this effect; global O3 compiles the entire single-effect image.
- No dwell-compression or ordered-cycle override was used.
- Provenance records O3 source `63268c3768ae91c8e2b47acd728db641570179eb` and paired shipping source
  `63268c3768ae91c8e2b47acd728db641570179eb`; no uncommitted working-tree state is recorded.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopePentBright`,
`HS_PROFILE_WINDOW=32`; `bash tools/profile_one.sh KaleidoscopePentBright profile_o3 70 32` builds, flashes, and captures
the global-O3 twin through the device lock.

## Global -O3 vs selective -O3

The paired shipping capture peaks at 33.37 ms versus
29.53 ms here: global O3 lowers the peak by 3.84 ms (11.5%). O3-image minus shipping-image
size deltas are **FLASH code +14,464 B** and **ITCM code
+11,712 B**. Spill fractions are compared rather than raw counts:
shipping 0/1088 (0%)
versus O3 0/1088 (0%).
