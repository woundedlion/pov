# MobiusGrid on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MobiusGrid`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/mobiusgrid_ship.log`. Captured 2026-11-04 23:03 America/Los_Angeles on COM4. Replaces profile_mobiusgrid_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MobiusGrid, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 105 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh MobiusGrid profile 105 16` |

Image size:

```text
teensy_size:   FLASH: code:96184, data:152816, headers:9048   free for files:1773568
teensy_size:    RAM1: variables:314112, code:25480, padding:7288   free for local variables:177408
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `cfdd7a6cf4140319ec4b2318026478752f79e2066492915f9f5903a9a5a08c16`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 689–704, root cyc/600 versus measured wall sum differs by **0.068 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   MobiusGrid               worst 2 lines/set (0 sets over)  gated 18,776 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **32.874 ms** (frame 700), **0/1647 spilled** (0.00%). Setup frame 1: 53.328 ms. Complete windows cover frames 1–1648; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 28.993 ms/frame. `fx_shader_draw` averages 27.000 ms/frame over complete windows; its worst window is 29.196 ms/frame (frames 689–704).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 26.874–32.874 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     59.610 ms 35.766 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x0.94/frame; 148.348 us/call
  fx_shader_draw                          27.744 ms 16.646 Mcyc  46.5%
    fx_prepare_column                      0.648 ms  0.389 Mcyc   1.1%
      x144.00/frame; 4.502 us/call
  fx_palette_base_bake                     0.109 ms  0.065 Mcyc   0.2%
    x1.00/frame; 108.939 us/call
  fx_prepare_frame                         0.280 ms  0.168 Mcyc   0.5%
    x1.00/frame; 279.642 us/call
  fx_advance                               1.003 ms  0.602 Mcyc   1.7%
    palette_cycle_step                     0.856 ms  0.513 Mcyc   1.4%
      x1.00/frame; 855.720 us/call
  fx_timeline_step                         0.123 ms  0.074 Mcyc   0.2%
    x1.00/frame; 122.564 us/call
  canvas_clear                             0.104 ms  0.062 Mcyc   0.2%
    x1.00/frame; 103.675 us/call
  canvas_buffer_wait                      29.839 ms 17.904 Mcyc  50.1%
    x1.00/frame; 29839.476 us/call
```

Wall min/avg/max: 26.931/59.609/65.183 ms. Window mean render is 29.770 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 689–704)

```text
frame                                     62.383 ms 37.430 Mcyc 100.0%
  pov_preserve_half                        0.144 ms  0.086 Mcyc   0.2%
    x1.00/frame; 144.076 us/call
  fx_shader_draw                          29.196 ms 17.517 Mcyc  46.8%
    fx_prepare_column                      1.936 ms  1.162 Mcyc   3.1%
      x144.00/frame; 13.445 us/call
  fx_palette_base_bake                     0.128 ms  0.077 Mcyc   0.2%
    x1.00/frame; 128.253 us/call
  fx_prepare_frame                         0.292 ms  0.175 Mcyc   0.5%
    x1.00/frame; 292.222 us/call
  fx_advance                               0.980 ms  0.588 Mcyc   1.6%
    palette_cycle_step                     0.948 ms  0.569 Mcyc   1.5%
      x1.00/frame; 947.534 us/call
  fx_timeline_step                         0.218 ms  0.131 Mcyc   0.3%
    x1.00/frame; 217.903 us/call
  canvas_clear                             0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 94.011 us/call
  canvas_buffer_wait                      31.117 ms 18.670 Mcyc  49.9%
    x1.00/frame; 31116.822 us/call
```

Wall min/avg/max: 59.614/62.383/65.150 ms. Window mean render is 31.267 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 593–608)

```text
frame                                     62.652 ms 37.591 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.514 us/call
  fx_shader_draw                          27.908 ms 16.745 Mcyc  44.5%
    fx_prepare_column                      1.334 ms  0.800 Mcyc   2.1%
      x144.00/frame; 9.261 us/call
  fx_palette_base_bake                     0.117 ms  0.070 Mcyc   0.2%
    x1.00/frame; 117.051 us/call
  fx_prepare_frame                         0.296 ms  0.177 Mcyc   0.5%
    x1.00/frame; 295.721 us/call
  fx_advance                               0.979 ms  0.587 Mcyc   1.6%
    palette_cycle_step                     0.935 ms  0.561 Mcyc   1.5%
      x1.00/frame; 935.276 us/call
  fx_timeline_step                         0.177 ms  0.106 Mcyc   0.3%
    x1.00/frame; 177.323 us/call
  canvas_clear                             0.096 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.215 us/call
  canvas_buffer_wait                      32.697 ms 19.618 Mcyc  52.2%
    x1.00/frame; 32696.894 us/call
```

Wall min/avg/max: 59.522/62.651/66.965 ms. Window mean render is 29.955 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 2  | 65/65 | unavailable | 29.196 | 31.267 | 16.03 | modal-call-count hold |
| 1  | 38/37 | unavailable | 27.744 | 29.770 | n/a | startup window (includes setup frame 1) |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 2 | 32.874 | 0/1049 | 16.00 |
| 1 | 30.604 | 0/598 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.19/frame
  min/avg/max 0.630/1.980/19.708 us; cpu 3.65%
  isr_pack         144.00/frame
    min/avg/max 8.348/9.313/11.875 us; cpu 2.14%
  isr_dma_submit   144.00/frame
    min/avg/max 0.653/0.940/9.888 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1000105 us
isr_wake       n=18435 cyc min/avg/max=378/1188/11825 ns min/avg/max=630/1980/19708 total=36506 us cpu=3.65%
isr_pack       n=2304 cyc min/avg/max=5009/5588/7125 ns min/avg/max=8348/9313/11875 total=21459 us cpu=2.14%
isr_dma_submit n=2304 cyc min/avg/max=392/564/5933 ns min/avg/max=653/940/9888 total=2166 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.65% leaves approximately 60.219 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 29.196 ms/frame, 46.80% of this root window.
2. `fx_prepare_column` — 1.936 ms/frame, 3.10% of this root window.
3. `fx_advance` — 0.980 ms/frame, 1.57% of this root window.
4. `palette_cycle_step` — 0.948 ms/frame, 1.52% of this root window.
5. `fx_prepare_frame` — 0.292 ms/frame, 0.47% of this root window.

README cells: peak 🟢 32.874 (2), spilled 🟢 0/1647 (0.00%).

README setup frame 1: 53.328 ms. Captured 2026-11-04 23:03 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=MobiusGrid  windows=103
  [PASS] captured >=3 windows (103)
  [PASS] every window names one effect and resolution (MobiusGrid 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 103 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 103 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.1 ppm, frames 689-704)

INVALID: MobiusGrid
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 2, 1 distinct); [FAIL] all 2 presets visited (1). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MobiusGrid`, `HS_PROFILE_WINDOW=16`. `just profile MobiusGrid` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
