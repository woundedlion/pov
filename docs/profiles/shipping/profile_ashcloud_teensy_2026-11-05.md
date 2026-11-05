# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/ashcloud_ship.log`. Captured 2026-11-05 00:10 America/Los_Angeles on COM4. Replaces profile_ashcloud_teensy_2026-10-31.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh AshCloud profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:122216, data:182684, headers:8444   free for files:1718272
teensy_size:    RAM1: variables:314144, code:39976, padding:25560   free for local variables:144608
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `a2ae4872e8d8ca6e487fb3a2d9a95b125c5d7af6354c3632028a51bbb31ee457`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 897–928, root cyc/600 versus measured wall sum differs by **5.042 ppm**. This fails the 5 ppm check; see validator diagnostics below.

Hot-layout gate in the measured profile build:

```text
  ok   AshCloud                 worst 2 lines/set (0 sets over)  gated 20,480 B  overflow 11,014 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **56.969 ms** (frame 920), **0/927 spilled** (0.00%). Setup frame 1: 81.521 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 49.044 ms/frame. `fx_shader_draw` averages 45.832 ms/frame over complete windows; its worst window is 49.213 ms/frame (frames 897–928).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 38.404–56.969 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.373 ms 37.424 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.467 us/call
  fx_shader_draw                          41.366 ms 24.819 Mcyc  66.3%
    fx_prepare_column                      1.716 ms  1.030 Mcyc   2.8%
      x144.00/frame; 11.916 us/call
  fx_palette_base_bake                     0.202 ms  0.121 Mcyc   0.3%
    fx_palette_rgb_bake                    0.049 ms  0.029 Mcyc   0.1%
      x1.00/frame; 48.826 us/call
  fx_prepare_frame                         0.895 ms  0.537 Mcyc   1.4%
    x1.00/frame; 895.238 us/call
  fx_advance                               1.201 ms  0.721 Mcyc   1.9%
    palette_cycle_step                     1.075 ms  0.645 Mcyc   1.7%
      x1.00/frame; 1075.026 us/call
  fx_timeline_step                         0.111 ms  0.067 Mcyc   0.2%
    x1.00/frame; 111.278 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.910 us/call
  canvas_buffer_wait                      18.103 ms 10.862 Mcyc  29.0%
    x1.00/frame; 18102.702 us/call
```

Wall min/avg/max: 38.479/62.372/81.525 ms. Window mean render is 44.270 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 897–928)

```text
frame                                     62.439 ms 37.464 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 136.597 us/call
  fx_shader_draw                          49.213 ms 29.528 Mcyc  78.8%
    fx_prepare_column                      1.853 ms  1.112 Mcyc   3.0%
      x144.00/frame; 12.866 us/call
  fx_palette_base_bake                     0.636 ms  0.382 Mcyc   1.0%
    fx_palette_rgb_bake                    0.488 ms  0.293 Mcyc   0.8%
      x1.00/frame; 487.875 us/call
  fx_prepare_frame                         1.081 ms  0.649 Mcyc   1.7%
    x1.00/frame; 1081.375 us/call
  fx_advance                               1.190 ms  0.714 Mcyc   1.9%
    palette_cycle_step                     1.111 ms  0.667 Mcyc   1.8%
      x1.00/frame; 1111.275 us/call
  fx_timeline_step                         0.112 ms  0.067 Mcyc   0.2%
    x1.00/frame; 112.345 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.352 us/call
  canvas_buffer_wait                       9.728 ms  5.837 Mcyc  15.6%
    x1.00/frame; 9727.562 us/call
```

Wall min/avg/max: 56.603/62.439/68.348 ms. Window mean render is 52.712 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 65–96)

```text
frame                                     62.329 ms 37.398 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 139.114 us/call
  fx_shader_draw                          42.781 ms 25.668 Mcyc  68.6%
    fx_prepare_column                      1.740 ms  1.044 Mcyc   2.8%
      x144.00/frame; 12.085 us/call
  fx_palette_base_bake                     0.267 ms  0.160 Mcyc   0.4%
    fx_palette_rgb_bake                    0.112 ms  0.067 Mcyc   0.2%
      x1.00/frame; 111.666 us/call
  fx_prepare_frame                         1.070 ms  0.642 Mcyc   1.7%
    x1.00/frame; 1069.586 us/call
  fx_advance                               1.173 ms  0.704 Mcyc   1.9%
    palette_cycle_step                     1.115 ms  0.669 Mcyc   1.8%
      x1.00/frame; 1114.992 us/call
  fx_timeline_step                         0.111 ms  0.067 Mcyc   0.2%
    x1.00/frame; 111.496 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.439 us/call
  canvas_buffer_wait                      16.435 ms  9.861 Mcyc  26.4%
    x1.00/frame; 16434.770 us/call
```

Wall min/avg/max: 58.970/62.329/66.161 ms. Window mean render is 45.894 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1153.69/frame
  min/avg/max 0.495/2.171/28.761 us; cpu 4.00%
  isr_pack         144.22/frame
    min/avg/max 8.831/9.543/12.321 us; cpu 2.19%
  isr_dma_submit   144.22/frame
    min/avg/max 0.616/0.930/4.508 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=2002836 us
isr_wake       n=36918 cyc min/avg/max=297/1303/17257 ns min/avg/max=495/2171/28761 total=80219 us cpu=4.00%
isr_pack       n=4615 cyc min/avg/max=5299/5726/7393 ns min/avg/max=8831/9543/12321 total=44045 us cpu=2.19%
isr_dma_submit n=4615 cyc min/avg/max=370/558/2705 ns min/avg/max=616/930/4508 total=4295 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 4.00% leaves approximately 60.000 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 49.213 ms/frame, 78.82% of this root window.
2. `fx_prepare_column` — 1.853 ms/frame, 2.97% of this root window.
3. `fx_advance` — 1.190 ms/frame, 1.91% of this root window.
4. `palette_cycle_step` — 1.111 ms/frame, 1.78% of this root window.
5. `fx_prepare_frame` — 1.081 ms/frame, 1.73% of this root window.

README cells: peak 🟢 56.969, spilled 🟢 0/927 (0.00%).

README setup frame 1: 81.521 ms. Captured 2026-11-05 00:10 America/Los_Angeles.

Pullback shader: also affected by the 1a3244da8..dc2954193 cached-phase/color fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=AshCloud  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (AshCloud 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [FAIL] root cyc/600 vs wall sum within 5 ppm (5.0 ppm, frames 897-928)

INVALID: AshCloud
```

Validation flags: [FAIL] root cyc/600 vs wall sum within 5 ppm (5.0 ppm, frames 897-928). The reported cadence is supported by direct per-frame render telemetry; unobserved wrap and cycle-counter exactness are not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=AshCloud`, `HS_PROFILE_WINDOW=32`. `just profile AshCloud` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
