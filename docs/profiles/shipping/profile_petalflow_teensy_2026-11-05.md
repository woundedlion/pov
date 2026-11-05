# PetalFlow on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile PetalFlow`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/petalflow_ship.log`. Captured 2026-11-04 23:04 America/Los_Angeles on COM4. Replaces profile_petalflow_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | PetalFlow, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-20261105-push-sweep HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh PetalFlow profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:112400, data:145780, headers:9084   free for files:1764352
teensy_size:    RAM1: variables:312992, code:30552, padding:2216   free for local variables:178528
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `f2542b447ca3094c7cb19889b005aa8bb3b6ee7adb8145ac3ce33c1eda351e0b`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 321–352, root cyc/600 versus measured wall sum differs by **0.676 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   PetalFlow                worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 35,152 B  full-set worst 4
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **25.961 ms** (frame 309), **0/447 spilled** (0.00%). Setup frame 1: 20.617 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 21.393 ms/frame. `pf_draw_rings` averages 21.078 ms/frame over complete windows; its worst window is 23.606 ms/frame (frames 321–352).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 17.451–25.961 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     59.713 ms 35.828 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.100 us/call
  pf_draw_rings                           19.273 ms 11.564 Mcyc  32.3%
    pf_ring_scan                          18.280 ms 10.968 Mcyc  30.6%
      x23.34/frame; 783.067 us/call
    pf_ring_build                          0.852 ms  0.511 Mcyc   1.4%
      x23.34/frame; 36.496 us/call
  pf_timeline_step                         0.063 ms  0.038 Mcyc   0.1%
    x1.00/frame; 62.853 us/call
  canvas_clear                             0.094 ms  0.057 Mcyc   0.2%
    x1.00/frame; 94.353 us/call
  canvas_buffer_wait                      40.124 ms 24.074 Mcyc  67.2%
    x1.00/frame; 40123.534 us/call
```

Wall min/avg/max: 18.065/59.713/67.219 ms. Window mean render is 19.590 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 289–320)

```text
frame                                     62.519 ms 37.511 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x1.00/frame; 135.026 us/call
  pf_draw_rings                           23.548 ms 14.129 Mcyc  37.7%
    pf_ring_scan                          22.563 ms 13.538 Mcyc  36.1%
      x23.34/frame; 966.536 us/call
    pf_ring_build                          0.851 ms  0.510 Mcyc   1.4%
      x23.34/frame; 36.444 us/call
  pf_timeline_step                         0.060 ms  0.036 Mcyc   0.1%
    x1.00/frame; 60.289 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.535 us/call
  canvas_buffer_wait                      38.662 ms 23.197 Mcyc  61.8%
    x1.00/frame; 38662.322 us/call
```

Wall min/avg/max: 58.484/62.518/66.016 ms. Window mean render is 23.856 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Worst dominant-scope regime (frames 321–352)

```text
frame                                     62.397 ms 37.438 Mcyc 100.0%
  pov_preserve_half                        0.136 ms  0.081 Mcyc   0.2%
    x1.00/frame; 135.551 us/call
  pf_draw_rings                           23.606 ms 14.164 Mcyc  37.8%
    pf_ring_scan                          22.615 ms 13.569 Mcyc  36.2%
      x23.34/frame; 968.763 us/call
    pf_ring_build                          0.849 ms  0.509 Mcyc   1.4%
      x23.34/frame; 36.359 us/call
  pf_timeline_step                         0.060 ms  0.036 Mcyc   0.1%
    x1.00/frame; 60.294 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.853 us/call
  canvas_buffer_wait                      38.482 ms 23.089 Mcyc  61.7%
    x1.00/frame; 38481.929 us/call
```

Wall min/avg/max: 59.091/62.397/65.036 ms. Window mean render is 23.915 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 33–64)

```text
frame                                     62.442 ms 37.465 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.281 us/call
  pf_draw_rings                           19.454 ms 11.672 Mcyc  31.2%
    pf_ring_scan                          18.464 ms 11.078 Mcyc  29.6%
      x23.34/frame; 790.956 us/call
    pf_ring_build                          0.851 ms  0.510 Mcyc   1.4%
      x23.34/frame; 36.438 us/call
  pf_timeline_step                         0.058 ms  0.035 Mcyc   0.1%
    x1.00/frame; 57.628 us/call
  canvas_clear                             0.089 ms  0.054 Mcyc   0.1%
    x1.00/frame; 89.337 us/call
  canvas_buffer_wait                      42.679 ms 25.607 Mcyc  68.3%
    x1.00/frame; 42678.738 us/call
```

Wall min/avg/max: 57.934/62.442/66.824 ms. Window mean render is 19.764 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1152.03/frame
  min/avg/max 0.618/1.678/18.328 us; cpu 3.09%
  isr_pack         144.00/frame
    min/avg/max 6.113/6.960/9.530 us; cpu 1.60%
  isr_dma_submit   144.00/frame
    min/avg/max 0.663/0.928/1.028 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999926 us
isr_wake       n=36865 cyc min/avg/max=371/1007/10997 ns min/avg/max=618/1678/18328 total=61897 us cpu=3.09%
isr_pack       n=4608 cyc min/avg/max=3668/4176/5718 ns min/avg/max=6113/6960/9530 total=32075 us cpu=1.60%
isr_dma_submit n=4608 cyc min/avg/max=398/557/617 ns min/avg/max=663/928/1028 total=4283 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.09% leaves approximately 60.569 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `pf_draw_rings` — 23.548 ms/frame, 37.67% of this root window.
2. `pf_ring_scan` — 22.563 ms/frame, 36.09% of this root window.
3. `pf_ring_build` — 0.851 ms/frame, 1.36% of this root window.
4. `pov_preserve_half` — 0.135 ms/frame, 0.22% of this root window.
5. `canvas_clear` — 0.090 ms/frame, 0.14% of this root window.

README cells: peak 🟢 25.961, spilled 🟢 0/447 (0.00%).

README setup frame 1: 20.617 ms. Captured 2026-11-04 23:04 America/Los_Angeles.

This capture measures the generated hot-code layout against the 11-04 matched sweep; see the sweep scorecard for the delta.

Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=PetalFlow  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (PetalFlow 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'pf_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.7 ppm, frames 321-352)

VALID: PetalFlow
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=PetalFlow`, `HS_PROFILE_WINDOW=32`. `just profile PetalFlow` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
