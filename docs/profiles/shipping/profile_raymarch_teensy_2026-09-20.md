# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-09-20, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `build/prof/raymarch_implementation_20260920/raymarch_stage2_final_ship.log`.
Replaces the earlier September 20 baseline report. The separate implementation
measurement archive is no longer retained.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 COM3, 600 MHz, live flywheel and DMA interrupts |
| Image | `profile`, -Os with selective O3 volume scan, SDF and Raymarch surface/shader |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Raymarch 288×144, default controls, `fb68b57885a76c604f82ef2304586c7aadfb4a9c` |
| Method | 110 s, 32-frame windows, deterministic seed; runtime frames 2–1737; frame 1 excluded; scope summaries use complete windows from frame 33 |
| Reproduce | `bash tools/profile_one.sh Raymarch profile 110 32` at the recorded source |

```text
FLASH: code:107264, data:183564, headers:9204   free for files:1731584
RAM1: variables:315072, code:37832, padding:27704   free for local variables:143680
RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 321–352, root cycles / 600 versus
measured wall sum agree within **3.21 ppm**. Raw-log parser validation passed.

## Frame cadence

**Peak live render: 56.071 ms at frame 349;
spilled 0/1736 live frames.** Headroom against 62.5 ms:
6.429 ms. Runtime mean render
48.577 ms is descriptive only;
acceptance and rankings use the peak and spill count.

Startup setup: 96.727 ms at frame 1, drawn before publication;
excluded from all runtime statistics, denominators and README rankings.
Startup renders 20,736 pixels; each live quadrant contains 10,368 pixels.
`canvas_buffer_wait` is synchronization idle to the next display flip.

## Phase-by-phase readout

Phase schedule: continuous rotation, default 26-volume placement, twist 2 and
18 primary steps. No preset dwell compression or epoch crossing.

### Window containing the live peak (frames 321–352)

Window averages below attribute cost; they do not substitute for peak timing.

```text
frame                       62.403 ms 37.442 Mcyc 100.0%
  pov_preserve_half          0.132 ms  0.079 Mcyc   0.2%
  rm_shader_draw            47.246 ms 28.348 Mcyc  75.7%
    filter_blend             0.520 ms  0.312 Mcyc   0.8%
  rm_timeline_step           0.375 ms  0.225 Mcyc   0.6%
  canvas_clear               0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait        11.439 ms  6.864 Mcyc  18.3%
```

Wall min/avg/max: 54.779/62.403/71.345 ms.
Geometry and shading remain the dominant work. Render includes interrupts;
the buffer wait records the spare interval before display publication.

### Per-pixel figures

5810.1 blends/frame versus 10,368 quadrant pixels;
53.7 cycles/blend and
4879.1 inclusive shader cycles per blended pixel.
Rejected rays and occlusion probes are included in the shader ratio.

## Column-ISR / DMA marshaling cost

Rates and min/mean/max service time for the same complete window:

```text
isr_wake        1151.0/frame  0.60/1.68/13.25 us  3.09%
isr_pack         143.9/frame  6.24/6.84/9.49 us  1.57%
isr_dma_submit   143.9/frame  0.61/0.93/1.10 us  0.21%
```

Packing dominates DMA submit CPU time. LED wire transfer is asynchronous.
Interrupt time is already included in the measured render budget; do not
subtract it a second time. The observed peak requires no further speedup to
meet the 62.5 ms live deadline.

## Summary ranking

1. `rm_shader_draw`: 47.246 ms/window frame,
   75.7% of wall time, inclusive march and shading.
2. `rm_timeline_step`: 0.375 ms/window frame.
3. Preserve/clear work: 0.217 ms/window frame.

Native framebuffer comparisons validate images; native timings are not used
to predict device deadlines. Earlier deep-stage numbers describe the baseline,
not this optimized image.

## Caveats

- All cycle scopes absorb interrupt time. No per-pixel deep instrumentation
  was enabled in acceptance captures.
- `filter_blend` counts blends, not unique pixels; parent selection can vary.
- Shipping selective O3 includes the compact Raymarch surface/shader functions.
- Captures used committed source with attested profile and full-roster ELFs.
- These are observed peaks for the default seeded trajectory, not a proof
  over every possible control setting or duration.
- Full-roster and single-effect firmware can make different inlining choices.

## Harness

`targets/Profile/Profile.ino`: `HS_PROFILE_TARGET=Raymarch`,
`HS_PROFILE_WINDOW=32`; `just profile Raymarch` builds, flashes and captures.
