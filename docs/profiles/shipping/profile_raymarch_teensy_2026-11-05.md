# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-11-05, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `C:/work/temp/holosphere-full-sweep-20261105/build/prof/raymarch_ship.log`. Captured 2026-11-04 23:06 America/Los_Angeles on COM4. Replaces profile_raymarch_teensy_2026-10-23.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: shipping `-Os` base, selective `HS_O3` hot paths |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | Raymarch, 288×144, single-entry playlist, source `e00f48da5cb4dec7b67bd32e5f6179186904f3bb` |
| Method | 30 s, 32-frame counter windows; stock sweep flags |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=COM4 HS_DEVICE_WAIT=3600 bash tools/profile_one.sh Raymarch profile 30 32` |

Image size:

```text
teensy_size:   FLASH: code:111832, data:186644, headers:8724   free for files:1724416
teensy_size:    RAM1: variables:313920, code:25880, padding:6888   free for local variables:177600
teensy_size:    RAM2: variables:520064  free for malloc/new:4224
```

Profile ELF SHA-256: `ef9e3a3bbc5df44c3128cce5332a3c4157864d1e2ba7211aaf9647d4ab199856`. Phantasm ELF SHA-256: `e8b226ce326cecb75a2ea9f2eaf8cd9a75d2704b9e28ac61b3854f34eb9c42c4`. Both archived ELF files are rehashed while generating this report; the adjacent provenance and envdump files retain the compiler/ABI and exact flags.

Exactness cross-check: frames 1–32, root cyc/600 versus measured wall sum differs by **0.355 ppm**. This passes the 5 ppm check.

Hot-layout gate in the measured profile build:

```text
  ok   Raymarch                 worst 2 lines/set (0 sets over)  gated 18,210 B  overflow      0 B  full-set worst 2
```

## Frame cadence

**Steady render, frame 1 excluded:** peak **52.693 ms** (frame 10), **0/447 spilled** (0.00%). Setup frame 1: 93.551 ms. Complete windows cover frames 1–448; trailing rows without a complete counter window are excluded. Transitions remain included.

Pass aggregate: steady mean render 47.617 ms/frame. `rm_shader_draw` averages 44.176 ms/frame over complete windows; its worst window is 45.814 ms/frame (frames 1–32).

A display window is 62.5 ms (480 RPM); one quadrant is 10,368 pixels. Render excludes `canvas_buffer_wait`, the round-up idle to the next display flip. The observed steady render range is 42.269–52.693 ms. All measured steady frames fit 16 fps.

## Phase-by-phase readout

Phase schedule follows the shipping timeline. The startup, peak-render and dominant-scope windows below bound the observed regimes. Scopes are inclusive; nested and tagged scopes overlap. Counter blocks retain every logged instance and its original indentation. Repeated names are separate counter instances; the validator lists MIXED-PARENT/DUPLICATE-NAME warnings, and these costs are not additive.

### Startup (frames 1–32)

```text
frame                                     62.820 ms 37.692 Mcyc 100.0%
  pov_preserve_half                        0.130 ms  0.078 Mcyc   0.2%
    x0.97/frame; 133.836 us/call
  rm_shader_draw                          45.814 ms 27.488 Mcyc  72.9%
    x1.00/frame; 45813.561 us/call
  rm_timeline_step                         0.388 ms  0.233 Mcyc   0.6%
    x1.00/frame; 388.322 us/call
  canvas_clear                             0.090 ms  0.054 Mcyc   0.1%
    x1.00/frame; 90.199 us/call
  canvas_buffer_wait                      13.476 ms  8.085 Mcyc  21.5%
    x1.00/frame; 13475.805 us/call
```

Wall min/avg/max: 51.271/62.819/93.552 ms. Window mean render is 49.344 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Lowest-cost steady regime (frames 129–160)

```text
frame                                     62.232 ms 37.339 Mcyc 100.0%
  pov_preserve_half                        0.133 ms  0.080 Mcyc   0.2%
    x1.00/frame; 132.843 us/call
  rm_shader_draw                          42.639 ms 25.583 Mcyc  68.5%
    x1.00/frame; 42638.698 us/call
  rm_timeline_step                         0.388 ms  0.233 Mcyc   0.6%
    x1.00/frame; 388.082 us/call
  canvas_clear                             0.086 ms  0.052 Mcyc   0.1%
    x1.00/frame; 85.840 us/call
  canvas_buffer_wait                      16.044 ms  9.627 Mcyc  25.8%
    x1.00/frame; 16044.259 us/call
```

Wall min/avg/max: 54.453/62.232/70.245 ms. Window mean render is 46.188 ms. The remainder of wall time includes display synchronization; scope means describe this regime rather than the single-frame peak.

### Per-pixel figures

This shipping timing image does not enable per-pixel blend counters. Blended pixels/frame, cycles/blend and scan cycles/blended pixel are unavailable; geometric coverage is not substituted.

## Column-ISR / DMA marshaling cost

Last complete counter window:

```text
isr_wake        1151.78/frame
  min/avg/max 0.543/1.800/16.461 us; cpu 3.32%
  isr_pack         143.97/frame
    min/avg/max 6.011/7.871/11.055 us; cpu 1.81%
  isr_dma_submit   143.97/frame
    min/avg/max 0.616/0.936/1.128 us; cpu 0.21%
```

Raw ISR readout and separately measured ISR interval:

```text
ISR window=1999465 us
isr_wake       n=36857 cyc min/avg/max=326/1080/9877 ns min/avg/max=543/1800/16461 total=66392 us cpu=3.32%
isr_pack       n=4607 cyc min/avg/max=3607/4723/6633 ns min/avg/max=6011/7871/11055 total=36267 us cpu=1.81%
isr_dma_submit n=4607 cyc min/avg/max=370/562/677 ns min/avg/max=616/936/1128 total=4316 us cpu=0.21%
```

`isr_wake` is inclusive: packing and submission run inside it. Use wake alone for measured interrupt CPU share; never add pack and submit. Submit records CPU cost, not asynchronous wire-transfer duration, which is not measured. DMA-completion and other interrupts are not separately measured. Render scopes already include interrupt time, so the measured 62.5 ms cadence budget needs no second ISR subtraction. The peak regime requires a 1.000× render speedup to fit one window (1× means it already fits).

The inclusive wake share of 3.32% leaves approximately 60.425 ms of foreground CPU time per 62.5 ms display window. This describes CPU availability; the observed render telemetry already includes that ISR work.

## Summary ranking

1. `rm_shader_draw` — 45.814 ms/frame, 72.93% of this root window.
2. `rm_timeline_step` — 0.388 ms/frame, 0.62% of this root window.
3. `pov_preserve_half` — 0.130 ms/frame, 0.21% of this root window.
4. `canvas_clear` — 0.090 ms/frame, 0.14% of this root window.

README cells: peak 🟢 52.693, spilled 🟢 0/447 (0.00%).

README setup frame 1: 93.551 ms. Captured 2026-11-04 23:06 America/Los_Angeles.


Native/WASM results are not a matched hardware timing baseline; the unchanged comparison set remains in the aggregate index.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs.
- `filter_blend` first-caller parenting and per-pixel scope overhead apply when enabled.
- Shipping selective `HS_O3` draw paths are active; no additional dwell-compression knobs were used.
- Clean detached source worktree; the default one-hour epoch is used.

Parser validation:

```text
effect=Raymarch  windows=14
  [PASS] captured >=3 windows (14)
  [PASS] every window names one effect and resolution (Raymarch 288x144)
  [PASS] frame numbers monotonic (no epoch reset) (0 resets)
  [INFO] no preset markers - non-cycling capture
  [PASS] per-frame render is render, not wall (0 of 14 windows have render == wall)
  [PASS] per-frame rows cover every frame of their window (0 of 14 windows are short)
  [PASS] scope 'rm_shader_draw' appears in at least one window
  [PASS] richest window carries a root 'frame' counter and a wall sum
  [PASS] root cyc/600 vs wall sum within 5 ppm (0.4 ppm, frames 1-32)

VALID: Raymarch
```

All parser checks pass.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Raymarch`, `HS_PROFILE_WINDOW=32`. `just profile Raymarch` builds, flashes and captures through the shared-device lock; use the exact command above for this sweep duration.
