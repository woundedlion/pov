# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-09-24, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: [raw capture](../capture20260924ringspin/data/ringspin_ship.log.txt), captured 2026-09-24 20:15 local on COM4.
Replaces `profile_ringspin_teensy_2026-08-26.md`. That older capture used a different source revision and included setup in its aggregates; this is not a controlled before/after comparison.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile`: `-Os` base, with shape-scan and ring-SDF selective-O3 regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `9442dbc204d556612f50980a5e41b082ca1d7961` |
| Method | 70 s capture, 32-frame windows, default 120 s epoch; exact runtime frames 2–1088, setup frame 1 excluded; scope/ISR summaries use complete windows 33–1088 |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE="$PWD" bash tools/profile_one.sh RingSpin profile 70 32` |

Image size: `FLASH: code:52968, data:148496, headers:8456, free:1821696` / `RAM1: variables:315072, code:33208, padding:32328, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 737–768, root 1,205,533,233 cycles ÷ 600 MHz versus measured wall sum 2,009,228 us differs by **2.96 ppm**. The parser validates the effect, monotonic frames, complete per-frame telemetry and absence of epoch resets. Firmware SHA-256 and source SHA match the archived provenance.

## Frame cadence

**Runtime aggregate**, setup excluded: mean render **30.179 ms/frame**, peak render **49.920 ms**, spilled **0/1087 live frames (0.0%)**. Mean wall time is 62.429 ms/frame, or 16.018 observed fps.
Startup setup render is **1.028 ms** for frame 1; it is excluded from all runtime denominators and peaks.

Every live frame fits the 62.5 ms display interval (16 fps cadence), with 12.580 ms peak margin. RingSpin renders one quadrant, approximately 10,368 pixels. `canvas_buffer_wait` is the round-up idle until the next display flip; wall extrema include that synchronization and do not represent render peaks.

Complete post-setup windows give `rs_draw_rings` a mean of 30.296 ms/frame and a worst window mean of 35.600 ms/frame at 737–768.

## Phase-by-phase readout

Phase schedule: four rings wander continuously; their 19-sample trails fill at startup and their projected coverage changes throughout the pass. There are no preset transitions. These two complete windows bound the observed ring-scan regimes.

### Highest-cost measured regime (frames 737–768)

```text
frame                       62.79 ms 37.67 Mcyc 100.0%
  pov_preserve_half         143.6 us  86.2 kcyc  0.2% x1.0 143.6us/call
  rs_draw_rings             35.60 ms 21.36 Mcyc 56.7%
    rs_ring_scan            35.01 ms 21.01 Mcyc 55.8%
      filter_blend           4.76 ms  2.86 Mcyc  7.6% x46520.1 61.4cyc/blend
  rs_timeline_step           72.5 us  43.5 kcyc  0.1% x1.0 72.5us/call
  canvas_clear               84.4 us  50.6 kcyc  0.1% x1.0 84.4us/call
  canvas_buffer_wait        26.89 ms 16.13 Mcyc 42.8% x1.0 26886.3us/call
```

Wall min/avg/max = 44.608/62.788/80.940 ms. `rs_draw_rings` averages 35.600 ms/frame, mostly `rs_ring_scan` at 35.009 ms/frame. Complete render averages 35.902 ms/frame; the changing ring coverage shifts work between render and display-sync idle.

### Lowest-cost measured regime (frames 65–96)

```text
frame                       62.42 ms 37.45 Mcyc 100.0%
  pov_preserve_half         141.2 us  84.7 kcyc  0.2% x1.0 141.2us/call
  rs_draw_rings             24.06 ms 14.43 Mcyc 38.5%
    rs_ring_scan            23.48 ms 14.09 Mcyc 37.6%
      filter_blend           3.32 ms  1.99 Mcyc  5.3% x32577.6 61.1cyc/blend
  rs_timeline_step           75.9 us  45.6 kcyc  0.1% x1.0 75.9us/call
  canvas_clear               84.2 us  50.5 kcyc  0.1% x1.0 84.2us/call
  canvas_buffer_wait        38.06 ms 22.84 Mcyc 61.0% x1.0 38061.4us/call
```

Wall min/avg/max = 52.159/62.421/71.557 ms. `rs_draw_rings` averages 24.057 ms/frame, mostly `rs_ring_scan` at 23.478 ms/frame. Complete render averages 24.360 ms/frame; the changing ring coverage shifts work between render and display-sync idle.

### Per-pixel figures

The high-cost window blends 46,520.1 pixels/frame (4.49× the quadrant), with 61.43 cycles/blend. `rs_ring_scan` uses 451.54 cycles per blended pixel. Overlap from the ring trails explains coverage above one quadrant.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.4/f  0.49/1.66/20.65 us  CPU 3.06%
  isr_pack        144.0/f  6.23/6.80/9.82 us  CPU 1.57%
  isr_dma_submit  144.0/f  0.62/0.94/9.22 us  CPU 0.22%
```

Each row gives calls/frame and min/avg/max time per call, across the complete post-setup windows. `isr_wake` is inclusive: `flywheel_isr()` calls the measured pack and submit blocks through `run_wake_sequence()`. The indented child costs are already inside the wake cost and must not be added to it.

- Pack plus submit consumes 1.114 ms of CPU per rendered frame; DMA submission starts an asynchronous transfer.
- At 24 MHz, the 600-byte image-plus-black strobe packet takes 200 us on the wire per column. This transfer overlaps foreground work and is not CPU marshaling time.
- Measured inclusive flywheel ISR share is 3.06%, leaving approximately 60.587 ms per 62.5 ms interval after that measured ISR cost. DMA-completion ISR cost is not separately quantified, so this is not a complete foreground CPU budget. Render counters already include ISR time: the 49.920 ms live peak fits the actual interval, so no speedup is needed.

## Summary ranking

1. `rs_draw_rings` — 48.50% of root time, 30.296 ms/frame; inclusive measured scope.
2. `rs_ring_scan` — 47.56% of root time, 29.708 ms/frame; inclusive measured scope.
3. `filter_blend` — 6.50% of root time, 4.058 ms/frame; inclusive measured scope.

No matched current WASM/native timing capture was used. Nested scope costs overlap; they must not be added.

## Caveats

- CYCCNT free-runs, so all foreground scopes include ISR time.
- `filter_blend` parents under `rs_ring_scan`; its subtree can be hidden when its first parent has no calls. It adds per-pixel scope overhead in both profiling configurations.
- Shipping crosses shape-scan and ring-SDF `HS_O3` regions; global O3 changes the rest of the image too. No deep instrumentation or dwell-compression flags were used.
- No uncommitted source changes were profiled. Archived source status records documentation-only profile README work in progress.
- Runtime statistics use every captured live frame. Scope trees and ISR summaries exclude the entire first window because window totals cannot isolate setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`; `just profile RingSpin` routes through the locked wrapper. The capture and provenance are published in [the evidence archive](../capture20260924ringspin/evidence.md); local build logs and ELF files are not included.
