# ShapeShifter on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile ShapeShifter`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/shapeshifter_ship.log`. Captured 2026-11-04 23:13 America/Los_Angeles on COM4. Replaces profile_shapeshifter_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ShapeShifter, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 269 s, 16-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh ShapeShifter profile 269 16` |

Image size:

```text
teensy_size:   FLASH: code:199976, data:152420, headers:9076   free for files:1670144
teensy_size:    RAM1: variables:314048, code:39128, padding:26408   free for local variables:144704
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `1c639d9d21ac05f3f8e10f85c21a825a7dda6be9660523aabae0ba88b5f3d5d4`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 3553–3568, root cyc/600 versus measured wall sum differs by **0.487 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   ShapeShifter             worst 2 lines/set (0 sets over)  gated 20,476 B  overflow 103,074 B  full-set worst 8
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **57.031 ms** (frame 3554), **0/4271 spilled** (0.00%). Setup frame 1: 78.200 ms. Complete windows cover frames 1–4272; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 38.126 ms/frame. `ss_draw_all` averages 37.775 ms/frame over complete windows; its worst window is 54.015 ms/frame (frames 3553–3568).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 17.581–57.031 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–16)

```text
frame                                     62.084 ms 37.251 Mcyc 100.0%
  pov_preserve_half                        0.128 ms  0.077 Mcyc   0.2%
    x0.94/frame; 136.208 us/call
  ss_draw_all                             43.882 ms 26.329 Mcyc  70.7%
    ss_surface_knots                       5.351 ms  3.210 Mcyc   8.6%
      ss_surface_mask                      3.641 ms  2.185 Mcyc   5.9%
        x20.25/frame; 179.814 us/call
  ss_timeline_step                         0.111 ms  0.066 Mcyc   0.2%
    x1.00/frame; 110.810 us/call
  ss_buffer_wait                          17.911 ms 10.747 Mcyc  28.8%
    canvas_clear                           0.107 ms  0.064 Mcyc   0.2%
      x1.00/frame; 106.647 us/call
    canvas_buffer_wait                    17.766 ms 10.660 Mcyc  28.6%
      x1.00/frame; 17765.954 us/call
```

Wall min/avg/max: 44.044/62.084/78.200 ms. Window mean render is 44.173 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 3553–3568)

```text
frame                                     62.206 ms 37.324 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.341 us/call
  ss_draw_all                             54.015 ms 32.409 Mcyc  86.8%
    x1.00/frame; 54014.849 us/call
  ss_timeline_step                         0.067 ms  0.040 Mcyc   0.1%
    x1.00/frame; 66.626 us/call
  ss_buffer_wait                           7.974 ms  4.784 Mcyc  12.8%
    canvas_clear                           0.087 ms  0.052 Mcyc   0.1%
      x1.00/frame; 87.196 us/call
    canvas_buffer_wait                     7.847 ms  4.708 Mcyc  12.6%
      x1.00/frame; 7846.528 us/call
```

Wall min/avg/max: 56.826/62.206/67.455 ms. Window mean render is 54.232 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Preset-transition regime (frames 3345–3360)

```text
frame                                     62.585 ms 37.551 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 137.967 us/call
  ss_draw_all                             46.080 ms 27.648 Mcyc  73.6%
    x1.00/frame; 46080.297 us/call
  ss_timeline_step                         0.074 ms  0.044 Mcyc   0.1%
    x1.00/frame; 73.833 us/call
  ss_buffer_wait                          16.279 ms  9.768 Mcyc  26.0%
    canvas_clear                           0.087 ms  0.052 Mcyc   0.1%
      x1.00/frame; 86.721 us/call
    canvas_buffer_wait                    16.153 ms  9.692 Mcyc  25.8%
      x1.00/frame; 16153.061 us/call
```

Wall min/avg/max: 53.341/62.585/71.879 ms. Window mean render is 46.306 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-preset table

Clean-hold windows use modal draw-call counts and exclude straddling ownership where possible. When no hold exists, the parser uses the nearest available anchor window. The cadence buckets include transitions and exclude setup frame 1 in README cells. Blend counts are unavailable in this timing image. Wrap and visited-preset diagnostics are retained below. The parser counts explicit advance markers for validation but infers the initial preset for ownership tables; an inferred initial row does not count as an explicit visit or a confirmed wrap.

Rows are ranked by dominant-scope cost; each row states its observation method.

| Preset / geometry | Windows / clean | Blended px/f | Dominant ms/f | Render ms/f | fps | Method |
|---|--:|--:|--:|--:|--:|---|
| 8  | 29/29 | unavailable | 54.015 | 54.232 | 16.08 | modal-call-count hold |
| 7  | 30/29 | unavailable | 51.234 | 51.452 | 15.95 | modal-call-count hold |
| 1  | 30/29 | unavailable | 48.577 | 48.876 | 16.03 | modal-call-count hold |
| 4  | 30/29 | unavailable | 44.819 | 45.051 | 16.07 | modal-call-count hold |
| 6  | 30/29 | unavailable | 43.482 | 43.699 | 16.02 | modal-call-count hold |
| 5  | 30/29 | unavailable | 40.764 | 40.995 | 16.15 | modal-call-count hold |
| 9  | 28/28 | unavailable | 35.639 | 35.890 | 16.18 | modal-call-count hold |
| 3  | 30/29 | unavailable | 31.528 | 31.753 | 16.08 | modal-call-count hold |
| 2  | 30/29 | unavailable | 25.471 | 25.694 | 15.97 | modal-call-count hold |

| Preset | Peak render ms | Spilled/frames | Cadence fps |
|---|--:|--:|--:|
| 8 | 57.031 | 0/479 | 16.00 |
| 1 | 54.816 | 0/470 | 16.00 |
| 7 | 53.838 | 0/479 | 16.00 |
| 5 | 48.673 | 0/479 | 16.00 |
| 4 | 47.398 | 0/479 | 16.00 |
| 6 | 45.040 | 0/479 | 16.00 |
| 3 | 38.853 | 0/479 | 16.00 |
| 9 | 38.113 | 0/448 | 16.00 |
| 2 | 29.893 | 0/479 | 16.00 |

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.75/frame
  min/avg/max 0.630/2.423/29.481 us; cpu 4.46%
  isr_pack         144.06/frame
    min/avg/max 7.495/8.453/10.776 us; cpu 1.94%
  isr_dma_submit   144.06/frame
    min/avg/max 0.626/0.935/2.086 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1000594 us
isr_wake       n=18444 cyc min/avg/max=378/1454/17689 ns min/avg/max=630/2423/29481 total=44725 us cpu=4.46%
isr_pack       n=2305 cyc min/avg/max=4497/5072/6466 ns min/avg/max=7495/8453/10776 total=19487 us cpu=1.94%
isr_dma_submit n=2305 cyc min/avg/max=376/561/1252 ns min/avg/max=626/935/2086 total=2157 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 4.46% leaves approximately 59.712 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `ss_draw_all` — 54.015 ms/frame, 86.83% of this root window.
2. `pov_preserve_half` — 0.137 ms/frame, 0.22% of this root window.
3. `canvas_clear` — 0.087 ms/frame, 0.14% of this root window.
4. `ss_timeline_step` — 0.067 ms/frame, 0.11% of this root window.

README cells: peak 🟢 57.031 (9), spilled 🟢 0/4271 (0.00%).

README setup frame 1: 78.200 ms. Captured 2026-11-04 23:13 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=ShapeShifter  windows=267
  [PASS] captured >=3 windows (267)
  [PASS] every window names one effect and resolution (ShapeShifter 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [FAIL] cycle returns to its first index (peak idx 9, 8 distinct)
  [FAIL] all 9 presets visited (8)
  [PASS] per-frame render is render, not wall (0 of 267 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 267 windows are short)
  [PASS] scope 'ss_draw_all' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 3553-3568)

INVALID: ShapeShifter
```

Validation flags: [FAIL] cycle returns to its first index (peak idx 9, 8 distinct); [FAIL] all 9 presets visited (8). The reported cadence is supported by direct per-frame render telemetry; an unobserved wrap is not claimed.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ShapeShifter`, `HS_PROFILE_WINDOW=16`. `just profile ShapeShifter` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
