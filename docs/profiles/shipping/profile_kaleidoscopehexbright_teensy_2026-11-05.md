# KaleidoscopeHexBright on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeHexBright`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/kaleidoscopehexbright_ship.log`. Captured 2026-11-05 00:35 America/Los_Angeles on COM4. Replaces profile_kaleidoscopehexbright_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeHexBright, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 105 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh KaleidoscopeHexBright profile 105 32` |

Image size:

```text
teensy_size:   FLASH: code:101984, data:152596, headers:8588   free for files:1768448
teensy_size:    RAM1: variables:314144, code:27752, padding:5016   free for local variables:177376
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `877b24b841bafa52e0cce0e026c52caf3e4612f9d676744df9d11c8f93ad0f9d`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 385–416, root cyc/600 versus measured wall sum differs by **0.366 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   KaleidoscopeHexBright    worst 2 lines/set (0 sets over)  gated 20,474 B  overflow  3,140 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **36.811 ms** (frame 407), **0/1631 spilled** (0.00%). Setup frame 1: 53.101 ms. Complete windows cover frames 1–1632; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 30.876 ms/frame. `fx_shader_draw` averages 25.339 ms/frame over complete windows; its worst window is 27.543 ms/frame (frames 385–416).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 25.999–36.811 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     61.031 ms 36.618 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x0.97/frame; 143.175 us/call
  fx_shader_draw                          24.615 ms 14.769 Mcyc  40.3%
    fx_prepare_column                      0.620 ms  0.372 Mcyc   1.0%
      x144.00/frame; 4.307 us/call
  fx_palette_base_bake                     0.224 ms  0.134 Mcyc   0.4%
    fx_palette_rgb_bake                    0.104 ms  0.062 Mcyc   0.2%
      x1.00/frame; 104.000 us/call
  fx_prepare_frame                         3.361 ms  2.017 Mcyc   5.5%
    x1.00/frame; 3361.133 us/call
  fx_advance                               1.072 ms  0.643 Mcyc   1.8%
    palette_cycle_step                     0.965 ms  0.579 Mcyc   1.6%
      x1.00/frame; 965.349 us/call
  fx_timeline_step                         0.097 ms  0.058 Mcyc   0.2%
    x1.00/frame; 96.824 us/call
  canvas_clear                             0.102 ms  0.061 Mcyc   0.2%
    x1.00/frame; 101.952 us/call
  canvas_buffer_wait                      31.209 ms 18.725 Mcyc  51.1%
    x1.00/frame; 31208.796 us/call
```

Wall min/avg/max: 28.673/61.030/65.214 ms. Window mean render is 29.822 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 385–416)

```text
frame                                     62.396 ms 37.437 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.084 Mcyc   0.2%
    x1.00/frame; 140.634 us/call
  fx_shader_draw                          27.543 ms 16.526 Mcyc  44.1%
    fx_prepare_column                      0.664 ms  0.398 Mcyc   1.1%
      x144.00/frame; 4.609 us/call
  fx_palette_base_bake                     1.225 ms  0.735 Mcyc   2.0%
    fx_palette_rgb_bake                    1.100 ms  0.660 Mcyc   1.8%
      x1.00/frame; 1099.759 us/call
  fx_prepare_frame                         3.329 ms  1.998 Mcyc   5.3%
    x1.00/frame; 3329.242 us/call
  fx_advance                               1.069 ms  0.641 Mcyc   1.7%
    palette_cycle_step                     1.020 ms  0.612 Mcyc   1.6%
      x1.00/frame; 1020.398 us/call
  fx_timeline_step                         0.099 ms  0.059 Mcyc   0.2%
    x1.00/frame; 99.046 us/call
  canvas_clear                             0.099 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.554 us/call
  canvas_buffer_wait                      28.686 ms 17.212 Mcyc  46.0%
    x1.00/frame; 28686.372 us/call
```

Wall min/avg/max: 56.747/62.395/67.935 ms. Window mean render is 33.709 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 577–608)

```text
frame                                     62.431 ms 37.459 Mcyc 100.0%
  pov_preserve_half                        0.142 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.706 us/call
  fx_shader_draw                          23.954 ms 14.373 Mcyc  38.4%
    fx_prepare_column                      0.757 ms  0.454 Mcyc   1.2%
      x144.00/frame; 5.258 us/call
  fx_palette_base_bake                     0.169 ms  0.102 Mcyc   0.3%
    fx_palette_rgb_bake                    0.053 ms  0.032 Mcyc   0.1%
      x1.00/frame; 53.116 us/call
  fx_prepare_frame                         3.282 ms  1.969 Mcyc   5.3%
    x1.00/frame; 3282.070 us/call
  fx_advance                               1.037 ms  0.622 Mcyc   1.7%
    palette_cycle_step                     0.985 ms  0.591 Mcyc   1.6%
      x1.00/frame; 985.316 us/call
  fx_timeline_step                         0.111 ms  0.066 Mcyc   0.2%
    x1.00/frame; 110.750 us/call
  canvas_clear                             0.098 ms  0.059 Mcyc   0.2%
    x1.00/frame; 98.428 us/call
  canvas_buffer_wait                      33.439 ms 20.064 Mcyc  53.6%
    x1.00/frame; 33439.173 us/call
```

Wall min/avg/max: 57.667/62.430/67.076 ms. Window mean render is 28.992 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 1  | 19/18 | unavailable | 27.543 | 33.709 | 16.03 | modal-call-count hold |
| 2  | 32/32 | unavailable | 26.969 | 32.424 | 16.03 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 1 | 36.811 | 0/598 | 16.00 |
| 2 | 34.698 | 0/1033 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.03/frame
  min/avg/max 0.496/1.991/29.876 us; cpu 3.67%
  isr_pack         144.00/frame
    min/avg/max 8.845/9.540/12.206 us; cpu 2.19%
  isr_dma_submit   144.00/frame
    min/avg/max 0.668/0.941/5.576 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999889 us
isr_wake       n=36865 cyc min/avg/max=298/1195/17926 ns min/avg/max=496/1991/29876 total=73442 us cpu=3.67%
isr_pack       n=4608 cyc min/avg/max=5307/5724/7324 ns min/avg/max=8845/9540/12206 total=43966 us cpu=2.19%
isr_dma_submit n=4608 cyc min/avg/max=401/565/3346 ns min/avg/max=668/941/5576 total=4342 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.67% leaves approximately 60.206 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `fx_shader_draw` — 27.543 ms/frame, 44.14% of this root window.
2. `fx_prepare_frame` — 3.329 ms/frame, 5.34% of this root window.
3. `fx_palette_base_bake` — 1.225 ms/frame, 1.96% of this root window.
4. `fx_palette_rgb_bake` — 1.100 ms/frame, 1.76% of this root window.
5. `fx_advance` — 1.069 ms/frame, 1.71% of this root window.

README cells: peak 🟢 36.811 (2), spilled 🟢 0/1631 (0.00%).

README setup frame 1: 53.101 ms. Captured 2026-11-05 00:35 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=KaleidoscopeHexBright  windows=51
  [PASS] captured >=3 windows (51)
  [PASS] every window names one effect and resolution (KaleidoscopeHexBright 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 2, 1 distinct)
  [FAIL] all 2 presets visited (1)
  [PASS] per-frame render is render, not wall (0 of 51 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 51 windows are short)
  [PASS] scope 'fx_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.4 ppm, frames 385-416)

INVALID: KaleidoscopeHexBright
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 2, 1 distinct); [FAIL] all 2 presets visited (1). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeHexBright`, `HS_PROFILE_WINDOW=32`. `just profile KaleidoscopeHexBright` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
