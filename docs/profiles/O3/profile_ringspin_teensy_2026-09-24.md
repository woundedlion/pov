# RingSpin on-device profile — Teensy 4.0, segmented mode (2026-09-24, **-O3**)

Global-O3 twin of the [shipping report](../shipping/profile_ringspin_teensy_2026-09-24.md).

Point-in-time snapshot (regenerate with `just profile RingSpin`).
Raw capture: [raw capture](../capture20260924ringspin/data/ringspin_o3.log.txt), captured 2026-09-24 20:18 local on COM4.
Replaces `profile_ringspin_teensy_2026-08-26.md`. That older capture used a different source revision and included setup in its aggregates; this is not a controlled before/after comparison.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | `profile_o3`: global `-O3 -ffast-math` reference |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | RingSpin 288×144, single-entry playlist, tip `9442dbc204d556612f50980a5e41b082ca1d7961` |
| Method | 70 s capture, 32-frame windows, default 120 s epoch; exact runtime frames 2–1088, setup frame 1 excluded; scope/ISR summaries use complete windows 33–1088 |
| Reproduce | `HS_TEENSY_PORT=COM4 HS_PROFILE_TREE="$PWD" bash tools/profile_one.sh RingSpin profile_o3 70 32` |

Image size: `FLASH: code:68368, data:148632, headers:8280, free:1806336` / `RAM1: variables:315072, code:46088, padding:19448, free:143680` / `RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 737–768, root 1,190,980,513 cycles ÷ 600 MHz versus measured wall sum 1,984,969 us differs by **0.74 ppm**. The parser validates the effect, monotonic frames, complete per-frame telemetry and absence of epoch resets. Firmware SHA-256 and source SHA match the archived provenance.

## Frame cadence

**Runtime aggregate**, setup excluded: mean render **29.428 ms/frame**, peak render **50.750 ms**, spilled **0/1087 live frames (0.0%)**. Mean wall time is 62.416 ms/frame, or 16.021 observed fps.
Startup setup render is **1.022 ms** for frame 1; it is excluded from all runtime denominators and peaks.

Every live frame fits the 62.5 ms display interval (16 fps cadence), with 11.750 ms peak margin. RingSpin renders one quadrant, approximately 10,368 pixels. `canvas_buffer_wait` is the round-up idle until the next display flip; wall extrema include that synchronization and do not represent render peaks.

Complete post-setup windows give `rs_draw_rings` a mean of 29.563 ms/frame and a worst window mean of 40.324 ms/frame at 737–768.

## Phase-by-phase readout

Phase schedule: four rings wander continuously; their 19-sample trails fill at startup and their projected coverage changes throughout the pass. There are no preset transitions. These two complete windows bound the observed ring-scan regimes.

### Highest-cost measured regime (frames 737–768)

```text
frame                       62.03 ms 37.22 Mcyc 100.0%
  pov_preserve_half         141.4 us  84.9 kcyc  0.2% x1.0 141.4us/call
  rs_draw_rings             40.32 ms 24.19 Mcyc 65.0%
    rs_ring_scan            39.81 ms 23.89 Mcyc 64.2%
      filter_blend           5.89 ms  3.53 Mcyc  9.5% x58593.7 60.3cyc/blend
  rs_timeline_step           69.8 us  41.9 kcyc  0.1% x1.0 69.8us/call
  canvas_clear               84.4 us  50.6 kcyc  0.1% x1.0 84.4us/call
  canvas_buffer_wait        21.41 ms 12.85 Mcyc 34.5% x1.0 21410.5us/call
```

Wall min/avg/max = 42.646/62.030/79.428 ms. `rs_draw_rings` averages 40.324 ms/frame, mostly `rs_ring_scan` at 39.814 ms/frame. Complete render averages 40.620 ms/frame; the changing ring coverage shifts work between render and display-sync idle.

### Lowest-cost measured regime (frames 993–1024)

```text
frame                       62.74 ms 37.64 Mcyc 100.0%
  pov_preserve_half         142.5 us  85.5 kcyc  0.2% x1.0 142.5us/call
  rs_draw_rings             22.23 ms 13.34 Mcyc 35.4%
    rs_ring_scan            21.73 ms 13.04 Mcyc 34.6%
      filter_blend           3.32 ms  1.99 Mcyc  5.3% x32988.7 60.3cyc/blend
  rs_timeline_step           69.5 us  41.7 kcyc  0.1% x1.0 69.5us/call
  canvas_clear               84.1 us  50.5 kcyc  0.1% x1.0 84.1us/call
  canvas_buffer_wait        40.21 ms 24.13 Mcyc 64.1% x1.0 40212.1us/call
```

Wall min/avg/max = 53.771/62.737/70.347 ms. `rs_draw_rings` averages 22.229 ms/frame, mostly `rs_ring_scan` at 21.731 ms/frame. Complete render averages 22.525 ms/frame; the changing ring coverage shifts work between render and display-sync idle.

### Per-pixel figures

The high-cost window blends 58,593.7 pixels/frame (5.65× the quadrant), with 60.28 cycles/blend. `rs_ring_scan` uses 407.69 cycles per blended pixel. Overlap from the ring trails explains coverage above one quadrant.

## Column-ISR / DMA marshaling cost

```text
isr_wake       1152.3/f  0.32/1.53/26.31 us  CPU 2.82%
  isr_pack        144.0/f  5.99/6.66/9.84 us  CPU 1.53%
  isr_dma_submit  144.0/f  0.62/0.93/6.22 us  CPU 0.21%
```

Each row gives calls/frame and min/avg/max time per call, across the complete post-setup windows. `isr_wake` is inclusive: `flywheel_isr()` calls the measured pack and submit blocks through `run_wake_sequence()`. The indented child costs are already inside the wake cost and must not be added to it.

- Pack plus submit consumes 1.093 ms of CPU per rendered frame; DMA submission starts an asynchronous transfer.
- At 24 MHz, the 600-byte image-plus-black strobe packet takes 200 us on the wire per column. This transfer overlaps foreground work and is not CPU marshaling time.
- Measured inclusive flywheel ISR share is 2.82%, leaving approximately 60.736 ms per 62.5 ms interval after that measured ISR cost. DMA-completion ISR cost is not separately quantified, so this is not a complete foreground CPU budget. Render counters already include ISR time: the 50.750 ms live peak fits the actual interval, so no speedup is needed.

## Summary ranking

1. `rs_draw_rings` — 47.33% of root time, 29.563 ms/frame; inclusive measured scope.
2. `rs_ring_scan` — 46.52% of root time, 29.058 ms/frame; inclusive measured scope.
3. `filter_blend` — 6.79% of root time, 4.243 ms/frame; inclusive measured scope.

No matched current WASM/native timing capture was used. Nested scope costs overlap; they must not be added.

## Caveats

- CYCCNT free-runs, so all foreground scopes include ISR time.
- `filter_blend` parents under `rs_ring_scan`; its subtree can be hidden when its first parent has no calls. It adds per-pixel scope overhead in both profiling configurations.
- Shipping crosses shape-scan and ring-SDF `HS_O3` regions; global O3 changes the rest of the image too. No deep instrumentation or dwell-compression flags were used.
- No uncommitted source changes were profiled. Archived source status records documentation-only profile README work in progress.
- Runtime statistics use every captured live frame. Scope trees and ISR summaries exclude the entire first window because window totals cannot isolate setup.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=RingSpin`, `HS_PROFILE_WINDOW=32`; `just profile RingSpin` routes through the locked wrapper. The capture and provenance are published in [the evidence archive](../capture20260924ringspin/evidence.md); local build logs and ELF files are not included.

## Global -O3 vs selective -O3

Mean live render is 30.179 ms shipping versus 29.428 ms global O3 (1.026× ratio). Peaks are 49.920/50.750 ms; both have zero live spills and hold 16 fps. Global O3 adds +15,400 B FLASH code and +12,880 B ITCM. This pair measures the same source and board, with separate 70-second captures.
