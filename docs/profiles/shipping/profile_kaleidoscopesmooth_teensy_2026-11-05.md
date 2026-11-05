# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeSmooth`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/kaleidoscopesmooth_ship.log`. Captured 2026-11-05 00:32 America/Los_Angeles on COM4. Replaces profile_kaleidoscopesmooth_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeSmooth, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 240 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh KaleidoscopeSmooth profile 240 16` |

Image size:

```text
teensy_size:   FLASH: code:103608, data:153336, headers:8272   free for files:1766400
teensy_size:    RAM1: variables:314144, code:28664, padding:4104   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `9560317e50f818fa75beaa1926d297ba781a77940403b48e0c6cf8dc48672040`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1057–1072, root cyc/600 versus measured wall sum differs by **9.031 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   KaleidoscopeSmooth       worst 2 lines/set (0 sets over)  gated 20,478 B  overflow  3,036 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **36.376 ms** (frame 1054), **0/3807 spilled** (0.00%). Setup frame 1: 56.382 ms. Complete windows cover frames 1–3808; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 28.461 ms/frame. `fx_shader_draw` averages 25.725 ms/frame over complete windows; its worst window is 31.266 ms/frame (frames 1057–1072).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 23.264–36.376 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     59.802 ms 35.881 Mcyc 100.0%
  pov_preserve_half                        0.134 ms  0.081 Mcyc   0.2%
    x0.94/frame; 143.321 us/call
  fx_shader_draw                          27.314 ms 16.388 Mcyc  45.7%
    fx_prepare_column                      0.779 ms  0.467 Mcyc   1.3%
      x144.00/frame; 5.410 us/call
  fx_palette_base_bake                     0.190 ms  0.114 Mcyc   0.3%
    fx_palette_rgb_bake                    0.072 ms  0.043 Mcyc   0.1%
      x1.00/frame; 71.587 us/call
  fx_prepare_frame                         0.585 ms  0.351 Mcyc   1.0%
    x1.00/frame; 585.337 us/call
  fx_advance                               1.244 ms  0.746 Mcyc   2.1%
    palette_cycle_step                     1.083 ms  0.650 Mcyc   1.8%
      x1.00/frame; 1083.314 us/call
  fx_timeline_step                         0.091 ms  0.055 Mcyc   0.2%
    x1.00/frame; 91.187 us/call
  canvas_clear                             0.108 ms  0.065 Mcyc   0.2%
    x1.00/frame; 107.535 us/call
  canvas_buffer_wait                      29.904 ms 17.943 Mcyc  50.0%
    x1.00/frame; 29904.181 us/call
```

Wall min/avg/max: 26.830/59.801/66.009 ms. Window mean render is 29.898 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 1041–1056)

```text
frame                                     62.458 ms 37.475 Mcyc 100.0%
  pov_preserve_half                        0.144 ms  0.087 Mcyc   0.2%
    x1.00/frame; 144.362 us/call
  fx_shader_draw                          30.998 ms 18.599 Mcyc  49.6%
    fx_prepare_column                      0.959 ms  0.575 Mcyc   1.5%
      x144.00/frame; 6.659 us/call
  fx_palette_base_bake                     0.762 ms  0.457 Mcyc   1.2%
    fx_palette_rgb_bake                    0.638 ms  0.383 Mcyc   1.0%
      x1.00/frame; 638.353 us/call
  fx_prepare_frame                         0.369 ms  0.221 Mcyc   0.6%
    x1.00/frame; 368.677 us/call
  fx_advance                               1.142 ms  0.685 Mcyc   1.8%
    palette_cycle_step                     1.121 ms  0.672 Mcyc   1.8%
      x1.00/frame; 1120.661 us/call
  fx_timeline_step                         0.163 ms  0.098 Mcyc   0.3%
    x1.00/frame; 162.776 us/call
  canvas_clear                             0.101 ms  0.061 Mcyc   0.2%
    x1.00/frame; 100.912 us/call
  canvas_buffer_wait                      28.603 ms 17.162 Mcyc  45.8%
    x1.00/frame; 28603.368 us/call
```

Wall min/avg/max: 57.955/62.457/66.741 ms. Window mean render is 33.855 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 1057–1072)

```text
frame                                     62.369 ms 37.421 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.281 us/call
  fx_shader_draw                          31.266 ms 18.760 Mcyc  50.1%
    fx_prepare_column                      0.948 ms  0.569 Mcyc   1.5%
      x144.00/frame; 6.586 us/call
  fx_palette_base_bake                     0.755 ms  0.453 Mcyc   1.2%
    fx_palette_rgb_bake                    0.631 ms  0.378 Mcyc   1.0%
      x1.00/frame; 630.819 us/call
  fx_prepare_frame                         0.365 ms  0.219 Mcyc   0.6%
    x1.00/frame; 364.535 us/call
  fx_advance                               1.153 ms  0.692 Mcyc   1.8%
    palette_cycle_step                     1.130 ms  0.678 Mcyc   1.8%
      x1.00/frame; 1130.078 us/call
  fx_timeline_step                         0.157 ms  0.094 Mcyc   0.3%
    x1.00/frame; 157.475 us/call
  canvas_clear                             0.101 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.801 us/call
  canvas_buffer_wait                      28.258 ms 16.955 Mcyc  45.3%
    x1.00/frame; 28258.337 us/call
```

Wall min/avg/max: 58.188/62.368/66.378 ms. Window mean render is 34.110 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 593–608)

```text
frame                                     62.693 ms 37.616 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 143.319 us/call
  fx_shader_draw                          29.796 ms 17.878 Mcyc  47.5%
    fx_prepare_column                      0.875 ms  0.525 Mcyc   1.4%
      x144.00/frame; 6.076 us/call
  fx_palette_base_bake                     0.188 ms  0.113 Mcyc   0.3%
    fx_palette_rgb_bake                    0.065 ms  0.039 Mcyc   0.1%
      x1.00/frame; 64.686 us/call
  fx_prepare_frame                         0.308 ms  0.185 Mcyc   0.5%
    x1.00/frame; 308.136 us/call
  fx_advance                               1.078 ms  0.647 Mcyc   1.7%
    palette_cycle_step                     1.052 ms  0.631 Mcyc   1.7%
      x1.00/frame; 1051.510 us/call
  fx_timeline_step                         0.126 ms  0.076 Mcyc   0.2%
    x1.00/frame; 125.981 us/call
  canvas_clear                             0.100 ms  0.060 Mcyc   0.2%
    x1.00/frame; 100.486 us/call
  canvas_buffer_wait                      30.754 ms 18.452 Mcyc  49.1%
    x1.00/frame; 30754.148 us/call
```

Wall min/avg/max: 58.385/62.692/66.445 ms. Window mean render is 31.939 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 2  | 67/66 | unavailable | 31.266 | 34.110 | 16.03 | modal-call-count hold |
| 1  | 38/37 | unavailable | 29.008 | 31.121 | 16.02 | modal-call-count hold |
| 3  | 68/67 | unavailable | 26.666 | 30.146 | 15.97 | modal-call-count hold |
| 4  | 65/65 | unavailable | 26.465 | 29.414 | 16.03 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 2 | 36.376 | 0/1079 | 16.00 |
| 1 | 36.079 | 0/598 | 16.00 |
| 3 | 33.267 | 0/1079 | 16.00 |
| 4 | 32.871 | 0/1051 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.62/frame
  min/avg/max 0.533/1.963/26.101 us; cpu 3.62%
  isr_pack         144.06/frame
    min/avg/max 8.511/9.233/11.636 us; cpu 2.12%
  isr_dma_submit   144.06/frame
    min/avg/max 0.781/0.936/1.095 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1000504 us
isr_wake       n=18442 cyc min/avg/max=320/1178/15661 ns min/avg/max=533/1963/26101 total=36221 us cpu=3.62%
isr_pack       n=2305 cyc min/avg/max=5107/5540/6982 ns min/avg/max=8511/9233/11636 total=21283 us cpu=2.12%
isr_dma_submit n=2305 cyc min/avg/max=469/562/657 ns min/avg/max=781/936/1095 total=2162 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.62% leaves approximately 60.237 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 30.998 ms/frame, 49.63% of this root window.
2. `fx_advance` — 1.142 ms/frame, 1.83% of this root window.
3. `palette_cycle_step` — 1.121 ms/frame, 1.79% of this root window.
4. `fx_prepare_column` — 0.959 ms/frame, 1.54% of this root window.
5. `fx_palette_base_bake` — 0.762 ms/frame, 1.22% of this root window.

README cells: peak 🟢 36.376 (4), spilled 🟢 0/3807 (0.00%).

README setup frame 1: 56.382 ms. Captured 2026-11-05 00:32 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=KaleidoscopeSmooth  windows=238
  [PASS] captured >=3 windows (238)
  [PASS] every window names one effect and resolution (KaleidoscopeSmooth 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 4, 3 distinct)
  [FAIL] all 4 presets visited (3)
  [PASS] per-frame render is render, not wall (0 of 238 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 238 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (9.0 ppm, frames 1057-1072)

INVALID: KaleidoscopeSmooth
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 4, 3 distinct); [FAIL] all 4 presets visited (3); [FAIL] root cyc/600 vs wall sum within 5 ppm (9.0 ppm, frames 1057-1072). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`. `just profile KaleidoscopeSmooth` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
