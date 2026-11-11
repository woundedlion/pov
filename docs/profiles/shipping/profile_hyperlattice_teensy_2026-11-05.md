# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HyperLattice`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/hyperlattice_ship.log`. Captured 2026-11-05 00:03 America/Los_Angeles on COM4. Replaces profile_hyperlattice_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 300 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh HyperLattice profile 300 16` |

Image size:

```text
teensy_size:   FLASH: code:145984, data:148696, headers:8424   free for files:1728512
teensy_size:    RAM1: variables:313024, code:13400, padding:19368   free for local variables:178496
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `1b3c6fc3fca1e9b2eb5bc5ce426339dece314f3e5b0ba96f0457999217859730`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 3089–3104, root cyc/600 versus measured wall sum differs by **0.164 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   HyperLattice             worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 65,162 B  full-set worst 6
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **53.649 ms** (frame 3209), **0/4767 spilled** (0.00%). Setup frame 1: 43.112 ms. Complete windows cover frames 1–4768; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 30.081 ms/frame. `hl_shader_draw` averages 27.733 ms/frame over complete windows; its worst window is 43.340 ms/frame (frames 3089–3104).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 21.039–53.649 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     58.857 ms 35.314 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.082 Mcyc   0.2%
    x0.94/frame; 145.300 us/call
  hl_shader_draw                          21.631 ms 12.978 Mcyc  36.8%
    x0.94/frame; 23072.686 us/call
  hl_shader_draw                           2.654 ms  1.592 Mcyc   4.5%
    x0.06/frame; 42462.030 us/call
  palette_cycle_step                       1.461 ms  0.877 Mcyc   2.5%
    x1.00/frame; 1461.022 us/call
  hl_timeline_step                         0.018 ms  0.011 Mcyc   0.0%
    x1.00/frame; 17.503 us/call
  canvas_clear                             0.095 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.824 us/call
  canvas_buffer_wait                      32.519 ms 19.511 Mcyc  55.3%
    x1.00/frame; 32518.930 us/call
```

Wall min/avg/max: 24.671/58.856/64.111 ms. Window mean render is 26.338 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 3201–3216)

```text
frame                                     62.015 ms 37.209 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.519 us/call
  hl_shader_draw                          42.542 ms 25.525 Mcyc  68.6%
    x1.00/frame; 42541.922 us/call
  palette_cycle_step                       1.614 ms  0.968 Mcyc   2.6%
    x1.00/frame; 1613.907 us/call
  hl_timeline_step                         0.032 ms  0.019 Mcyc   0.1%
    x1.00/frame; 31.778 us/call
  canvas_clear                             0.094 ms  0.056 Mcyc   0.2%
    x1.00/frame; 93.840 us/call
  canvas_buffer_wait                      17.243 ms 10.346 Mcyc  27.8%
    x1.00/frame; 17242.756 us/call
```

Wall min/avg/max: 51.574/62.015/74.075 ms. Window mean render is 44.773 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 3089–3104)

```text
frame                                     62.429 ms 37.457 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.319 us/call
  hl_shader_draw                          43.340 ms 26.004 Mcyc  69.4%
    x1.00/frame; 43340.438 us/call
  palette_cycle_step                       1.578 ms  0.947 Mcyc   2.5%
    x1.00/frame; 1577.769 us/call
  hl_timeline_step                         0.017 ms  0.010 Mcyc   0.0%
    x1.00/frame; 17.126 us/call
  canvas_clear                             0.093 ms  0.056 Mcyc   0.1%
    x1.00/frame; 92.660 us/call
  canvas_buffer_wait                      16.910 ms 10.146 Mcyc  27.1%
    x1.00/frame; 16909.661 us/call
```

Wall min/avg/max: 59.875/62.428/64.900 ms. Window mean render is 45.519 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 2673–2688)

```text
frame                                     63.161 ms 37.896 Mcyc 100.0%
  hl_shader_draw                           2.639 ms  1.583 Mcyc   4.2%
    x0.06/frame; 42219.273 us/call
  pov_preserve_half                        0.140 ms  0.084 Mcyc   0.2%
    x1.00/frame; 139.904 us/call
  hl_shader_draw                          37.094 ms 22.256 Mcyc  58.7%
    x0.94/frame; 39567.002 us/call
  palette_cycle_step                       1.589 ms  0.953 Mcyc   2.5%
    x1.00/frame; 1589.024 us/call
  hl_timeline_step                         0.038 ms  0.023 Mcyc   0.1%
    x1.00/frame; 38.453 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.838 us/call
  canvas_buffer_wait                      21.213 ms 12.728 Mcyc  33.6%
    x1.00/frame; 21213.226 us/call
```

Wall min/avg/max: 59.505/63.160/75.966 ms. Window mean render is 41.947 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 6  | 35/34 | unavailable | 43.340 | 45.519 | 16.02 | modal-call-count hold |
| 9  | 26/25 | unavailable | 32.819 | 35.001 | 16.01 | modal-call-count hold |
| 5  | 43/42 | unavailable | 29.550 | 31.765 | 16.01 | modal-call-count hold |
| 8  | 42/42 | unavailable | 29.444 | 31.613 | 15.98 | modal-call-count hold |
| 3  | 35/34 | unavailable | 29.173 | 31.345 | 15.85 | modal-call-count hold |
| 4  | 27/26 | unavailable | 28.168 | 30.347 | 16.00 | modal-call-count hold |
| 2  | 43/42 | unavailable | 25.663 | 27.793 | 16.03 | modal-call-count hold |
| 7  | 27/26 | unavailable | 25.308 | 27.454 | 15.96 | modal-call-count hold |
| 1  | 20/18 | unavailable | 23.790 | 25.941 | 16.05 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 6 | 53.649 | 0/559 | 16.00 |
| 3 | 39.088 | 0/559 | 16.00 |
| 9 | 36.225 | 0/416 | 16.00 |
| 8 | 34.626 | 0/679 | 16.00 |
| 4 | 33.252 | 0/439 | 16.00 |
| 5 | 32.354 | 0/679 | 16.00 |
| 7 | 30.523 | 0/439 | 16.00 |
| 2 | 29.997 | 0/679 | 16.00 |
| 1 | 26.763 | 0/318 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1153.00/frame
  min/avg/max 0.613/1.791/19.801 us; cpu 3.30%
  isr_pack         144.12/frame
    min/avg/max 6.025/7.815/10.331 us; cpu 1.80%
  isr_dma_submit   144.12/frame
    min/avg/max 0.726/0.931/3.536 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1000802 us
isr_wake       n=18448 cyc min/avg/max=368/1075/11881 ns min/avg/max=613/1791/19801 total=33063 us cpu=3.30%
isr_pack       n=2306 cyc min/avg/max=3615/4689/6199 ns min/avg/max=6025/7815/10331 total=18023 us cpu=1.80%
isr_dma_submit n=2306 cyc min/avg/max=436/559/2122 ns min/avg/max=726/931/3536 total=2149 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.30% leaves approximately 60.438 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `hl_shader_draw` — 42.542 ms/frame, 68.60% of this root window.
2. `palette_cycle_step` — 1.614 ms/frame, 2.60% of this root window.
3. `pov_preserve_half` — 0.139 ms/frame, 0.22% of this root window.
4. `canvas_clear` — 0.094 ms/frame, 0.15% of this root window.
5. `hl_timeline_step` — 0.032 ms/frame, 0.05% of this root window.

README cells: peak 🟢 53.649 (9), spilled 🟢 0/4767 (0.00%).

README setup frame 1: 43.112 ms. Captured 2026-11-05 00:03 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=HyperLattice  windows=298
  [PASS] captured >=3 windows (298)
  [PASS] every window names one effect and resolution (HyperLattice 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 298 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 298 windows are short)
  [INFO] not exclusive costs (MIXED-PARENT/DUPLICATE-NAME): hl_shader_draw
  [PASS] scope 'hl_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.2 ppm, frames 3089-3104)

INVALID: HyperLattice
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 9, 8 distinct); [FAIL] all 9 presets visited (8). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`. `just profile HyperLattice` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
