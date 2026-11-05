# HopfFibration on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile HopfFibration`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/hopffibration_ship.log`. Captured 2026-11-04 23:01 America/Los_Angeles on COM4. Replaces profile_hopffibration_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HopfFibration, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh HopfFibration profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:128272, data:145832, headers:8520   free for files:1748992
teensy_size:    RAM1: variables:312992, code:30152, padding:2616   free for local variables:178528
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `6d7343fd94c09127b18a251107b92d9940297604f96695b6168111da1e0039e8`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 193–224, root cyc/600 versus measured wall sum differs by **0.515 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   HopfFibration            worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 51,496 B  full-set worst 5
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **57.199 ms** (frame 824), **0/927 spilled** (0.00%). Setup frame 1: 0.646 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 27.459 ms/frame. `hf_render_trails` averages 26.894 ms/frame over complete windows; its worst window is 47.534 ms/frame (frames 193–224).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 1.870–57.199 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     59.336 ms 35.602 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.084 Mcyc   0.2%
    x0.97/frame; 143.792 us/call
  hf_render_trails                        11.091 ms  6.655 Mcyc  18.7%
    x1.00/frame; 11091.023 us/call
  hf_project_record                        0.166 ms  0.099 Mcyc   0.3%
    x1.00/frame; 165.775 us/call
  hf_advance_tumble                        0.006 ms  0.004 Mcyc   0.0%
    x1.00/frame; 6.330 us/call
  hf_timeline_step                         0.077 ms  0.046 Mcyc   0.1%
    x1.00/frame; 77.155 us/call
  canvas_clear                             0.093 ms  0.056 Mcyc   0.2%
    x1.00/frame; 92.851 us/call
  canvas_buffer_wait                      47.742 ms 28.645 Mcyc  80.5%
    x1.00/frame; 47741.842 us/call
```

Wall min/avg/max: 0.646/59.335/69.904 ms. Window mean render is 11.594 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 801–832)

```text
frame                                     62.680 ms 37.608 Mcyc 100.0%
  pov_preserve_half                        0.139 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.572 us/call
  hf_render_trails                        40.792 ms 24.475 Mcyc  65.1%
    x1.00/frame; 40792.367 us/call
  hf_project_record                        0.193 ms  0.116 Mcyc   0.3%
    x1.00/frame; 192.990 us/call
  hf_advance_tumble                        0.005 ms  0.003 Mcyc   0.0%
    x1.00/frame; 5.161 us/call
  hf_timeline_step                         0.094 ms  0.056 Mcyc   0.1%
    x1.00/frame; 93.621 us/call
  canvas_clear                             0.091 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.475 us/call
  canvas_buffer_wait                      21.339 ms 12.804 Mcyc  34.0%
    x1.00/frame; 21339.268 us/call
```

Wall min/avg/max: 31.413/62.679/95.292 ms. Window mean render is 41.340 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 193–224)

```text
frame                                     62.834 ms 37.700 Mcyc 100.0%
  pov_preserve_half                        0.137 ms  0.082 Mcyc   0.2%
    x1.00/frame; 137.129 us/call
  hf_render_trails                        47.534 ms 28.521 Mcyc  75.7%
    x1.00/frame; 47534.329 us/call
  hf_project_record                        0.196 ms  0.118 Mcyc   0.3%
    x1.00/frame; 196.041 us/call
  hf_advance_tumble                        0.005 ms  0.003 Mcyc   0.0%
    x1.00/frame; 4.923 us/call
  hf_timeline_step                         0.091 ms  0.055 Mcyc   0.1%
    x1.00/frame; 90.948 us/call
  canvas_clear                             0.092 ms  0.055 Mcyc   0.1%
    x1.00/frame; 91.790 us/call
  canvas_buffer_wait                      14.754 ms  8.852 Mcyc  23.5%
    x1.00/frame; 14753.578 us/call
```

Wall min/avg/max: 41.433/62.833/82.505 ms. Window mean render is 48.080 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 129–160)

```text
frame                                     62.331 ms 37.399 Mcyc 100.0%
  pov_preserve_half                        0.143 ms  0.086 Mcyc   0.2%
    x1.00/frame; 142.773 us/call
  hf_render_trails                        13.341 ms  8.005 Mcyc  21.4%
    x1.00/frame; 13341.221 us/call
  hf_project_record                        0.170 ms  0.102 Mcyc   0.3%
    x1.00/frame; 169.745 us/call
  hf_advance_tumble                        0.008 ms  0.005 Mcyc   0.0%
    x1.00/frame; 8.072 us/call
  hf_timeline_step                         0.076 ms  0.046 Mcyc   0.1%
    x1.00/frame; 75.855 us/call
  canvas_clear                             0.088 ms  0.053 Mcyc   0.1%
    x1.00/frame; 87.877 us/call
  canvas_buffer_wait                      48.487 ms 29.092 Mcyc  77.8%
    x1.00/frame; 48487.222 us/call
```

Wall min/avg/max: 41.354/62.331/84.019 ms. Window mean render is 13.844 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.97/frame
  min/avg/max 0.568/1.733/19.730 us; cpu 3.19%
  isr_pack         143.97/frame
    min/avg/max 5.991/7.310/10.371 us; cpu 1.68%
  isr_dma_submit   143.97/frame
    min/avg/max 0.681/0.925/1.211 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999848 us
isr_wake       n=36863 cyc min/avg/max=341/1040/11838 ns min/avg/max=568/1733/19730 total=63912 us cpu=3.19%
isr_pack       n=4607 cyc min/avg/max=3595/4386/6223 ns min/avg/max=5991/7310/10371 total=33681 us cpu=1.68%
isr_dma_submit n=4607 cyc min/avg/max=409/555/727 ns min/avg/max=681/925/1211 total=4268 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.19% leaves approximately 60.506 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `hf_render_trails` — 40.792 ms/frame, 65.08% of this root window.
2. `hf_project_record` — 0.193 ms/frame, 0.31% of this root window.
3. `pov_preserve_half` — 0.139 ms/frame, 0.22% of this root window.
4. `hf_timeline_step` — 0.094 ms/frame, 0.15% of this root window.
5. `canvas_clear` — 0.091 ms/frame, 0.15% of this root window.

README cells: peak 🟢 57.199, spilled 🟢 0/927 (0.00%).

README setup frame 1: 0.646 ms. Captured 2026-11-04 23:01 America/Los_Angeles.

Also affected by the 1a3244da8..dc2954193 placement/render fixes; this delta does not isolate layout.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=HopfFibration  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (HopfFibration 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'hf_render_trails' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 193-224)

VALID: HopfFibration
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HopfFibration`, `HS_PROFILE_WINDOW=32`. `just profile HopfFibration` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
