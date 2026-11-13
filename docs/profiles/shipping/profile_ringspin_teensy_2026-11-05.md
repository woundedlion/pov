# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/ringspin_ship.log`. Captured 2026-11-04 23:07 America/Los_Angeles on COM4. Replaces profile_ringspin_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | RingSpin, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 60 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh RingSpin profile 60 32` |

Image size:

```text
teensy_size:   FLASH: code:91896, data:148164, headers:8772   free for files:1782784
teensy_size:    RAM1: variables:313984, code:27576, padding:5192   free for local variables:177536
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `7025c0416ab48027e32407da55c455383f08cc1108bf233f6e50d7d53145b886`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 737–768, root cyc/600 versus measured wall sum differs by **0.509 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   RingSpin                 worst 2 lines/set (0 sets over)  gated 20,478 B  overflow 18,930 B  full-set worst 3
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **53.875 ms** (frame 761), **0/927 spilled** (0.00%). Setup frame 1: 3.399 ms. Complete windows cover frames 1–928; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 38.715 ms/frame. `rs_draw_rings` averages 38.371 ms/frame over complete windows; its worst window is 41.485 ms/frame (frames 737–768).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 3.415–53.875 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     59.506 ms 35.704 Mcyc 100.0%
  pov_preserve_half                        0.135 ms  0.081 Mcyc   0.2%
    x0.97/frame; 139.108 us/call
  rs_draw_rings                           25.079 ms 15.047 Mcyc  42.1%
    rs_ring_scan                          24.358 ms 14.615 Mcyc  40.9%
      x54.62/frame; 445.919 us/call
  rs_timeline_step                         0.073 ms  0.044 Mcyc   0.1%
    x1.00/frame; 72.611 us/call
  canvas_clear                             0.089 ms  0.054 Mcyc   0.2%
    x1.00/frame; 89.324 us/call
  canvas_buffer_wait                      34.122 ms 20.473 Mcyc  57.3%
    x1.00/frame; 34121.667 us/call
```

Wall min/avg/max: 3.400/59.505/89.089 ms. Window mean render is 25.384 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Peak-render regime (frames 737–768)

```text
frame                                     62.616 ms 37.569 Mcyc 100.0%
  pov_preserve_half                        0.141 ms  0.085 Mcyc   0.2%
    x1.00/frame; 141.271 us/call
  rs_draw_rings                           41.485 ms 24.891 Mcyc  66.3%
    rs_ring_scan                          40.460 ms 24.276 Mcyc  64.6%
      x76.00/frame; 532.374 us/call
  rs_timeline_step                         0.071 ms  0.043 Mcyc   0.1%
    x1.00/frame; 71.320 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.363 us/call
  canvas_buffer_wait                      20.824 ms 12.495 Mcyc  33.3%
    x1.00/frame; 20824.433 us/call
```

Wall min/avg/max: 43.428/62.615/81.381 ms. Window mean render is 41.791 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 673–704)

```text
frame                                     62.389 ms 37.434 Mcyc 100.0%
  pov_preserve_half                        0.138 ms  0.083 Mcyc   0.2%
    x1.00/frame; 138.236 us/call
  rs_draw_rings                           36.012 ms 21.607 Mcyc  57.7%
    rs_ring_scan                          35.005 ms 21.003 Mcyc  56.1%
      x76.00/frame; 460.595 us/call
  rs_timeline_step                         0.071 ms  0.043 Mcyc   0.1%
    x1.00/frame; 70.842 us/call
  canvas_clear                             0.085 ms  0.051 Mcyc   0.1%
    x1.00/frame; 85.004 us/call
  canvas_buffer_wait                      26.075 ms 15.645 Mcyc  41.8%
    x1.00/frame; 26075.157 us/call
```

Wall min/avg/max: 53.080/62.389/72.396 ms. Window mean render is 36.314 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1149.16/frame
  min/avg/max 0.560/1.863/16.988 us; cpu 3.43%
  isr_pack         143.66/frame
    min/avg/max 6.365/8.343/11.883 us; cpu 1.92%
  isr_dma_submit   143.66/frame
    min/avg/max 0.655/0.936/2.336 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1994963 us
isr_wake       n=36773 cyc min/avg/max=336/1118/10193 ns min/avg/max=560/1863/16988 total=68548 us cpu=3.43%
isr_pack       n=4597 cyc min/avg/max=3819/5006/7130 ns min/avg/max=6365/8343/11883 total=38361 us cpu=1.92%
isr_dma_submit n=4597 cyc min/avg/max=393/562/1402 ns min/avg/max=655/936/2336 total=4309 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.43% leaves approximately 60.356 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `rs_draw_rings` — 41.485 ms/frame, 66.25% of this root window.
2. `rs_ring_scan` — 40.460 ms/frame, 64.62% of this root window.
3. `pov_preserve_half` — 0.141 ms/frame, 0.23% of this root window.
4. `canvas_clear` — 0.085 ms/frame, 0.14% of this root window.
5. `rs_timeline_step` — 0.071 ms/frame, 0.11% of this root window.

README cells: peak 🟢 53.875, spilled 🟢 0/927 (0.00%).

README setup frame 1: 3.399 ms. Captured 2026-11-04 23:07 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=RingSpin  windows=29
  [PASS] captured >=3 windows (29)
  [PASS] every window names one effect and resolution (RingSpin 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 29 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 29 windows are short)
  [PASS] scope 'rs_draw_rings' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.5 ppm, frames 737-768)

VALID: RingSpin
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`. `just profile RingSpin` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
