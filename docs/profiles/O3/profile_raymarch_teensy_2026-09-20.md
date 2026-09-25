# Raymarch on-device profile — Teensy 4.0, segmented mode (2026-09-20, **-O3**)

Point-in-time snapshot (regenerate with `just profile Raymarch`).
Raw capture: `build/prof/raymarch_implementation_20260920/raymarch_final_o3.log`.
Replaces the earlier September 20 baseline report. The separate implementation
measurement archive is no longer retained.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 COM4, 600 MHz, live flywheel and DMA interrupts |
| Image | `profile_o3`, global -O3 -ffast-math reference |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Raymarch 288×144, default controls, `fb68b57885a76c604f82ef2304586c7aadfb4a9c` |
| Method | 110 s, 32-frame windows, deterministic seed; runtime frames 2–1737; frame 1 excluded; scope summaries use complete windows from frame 33 |
| Reproduce | `bash tools/profile_one.sh Raymarch profile_o3 110 32` at the recorded source |

```text
FLASH: code:119832, data:183744, headers:8744   free for files:1719296
RAM1: variables:315104, code:45736, padding:19800   free for local variables:143648
RAM2: variables:520064  free for malloc/new:4224
```

Exactness cross-check: frames 321–352, root cycles / 600 versus
measured wall sum agree within **1.60 ppm**. Raw-log parser validation passed.

## Frame cadence

**Peak live render: 56.044 ms at frame 351;
spilled 0/1736 live frames.** Headroom against 62.5 ms:
6.456 ms. Runtime mean render
47.936 ms is descriptive only;
acceptance and rankings use the peak and spill count.

Startup setup: 96.849 ms at frame 1, drawn before publication;
excluded from all runtime statistics, denominators and README rankings.
Startup renders 20,736 pixels; each live quadrant contains 10,368 pixels.
`canvas_buffer_wait` is synchronization idle to the next display flip.

## Phase-by-phase readout

Phase schedule: continuous rotation, default 26-volume placement, twist 2 and
18 primary steps. No preset dwell compression or epoch crossing.

### Window containing the live peak (frames 321–352)

Window averages below attribute cost; they do not substitute for peak timing.

```text
frame                       62.381 ms 37.429 Mcyc 100.0%
  pov_preserve_half          0.132 ms  0.079 Mcyc   0.2%
  rm_shader_draw            46.182 ms 27.709 Mcyc  74.0%
    filter_blend             0.466 ms  0.280 Mcyc   0.7%
  rm_timeline_step           0.288 ms  0.173 Mcyc   0.5%
  canvas_clear               0.085 ms  0.051 Mcyc   0.1%
  canvas_buffer_wait        12.841 ms  7.704 Mcyc  20.6%
```

Wall min/avg/max: 54.536/62.381/73.462 ms.
Geometry and shading remain the dominant work. Render includes interrupts;
the buffer wait records the spare interval before display publication.

### Per-pixel figures

5707.3 blends/frame versus 10,368 quadrant pixels;
49.0 cycles/blend and
4855.0 inclusive shader cycles per blended pixel.
Rejected rays and occlusion probes are included in the shader ratio.

## Column-ISR / DMA marshaling cost

Rates and min/mean/max service time for the same complete window:

```text
isr_wake        1150.8/frame  0.51/1.56/18.71 us  2.87%
isr_pack         143.8/frame  5.99/6.64/9.29 us  1.53%
isr_dma_submit   143.8/frame  0.63/0.93/1.01 us  0.21%
```

Packing dominates DMA submit CPU time. LED wire transfer is asynchronous.
Interrupt time is already included in the measured render budget; do not
subtract it a second time. The observed peak requires no further speedup to
meet the 62.5 ms live deadline.

## Summary ranking

1. `rm_shader_draw`: 46.182 ms/window frame,
   74.0% of wall time, inclusive march and shading.
2. `rm_timeline_step`: 0.288 ms/window frame.
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

## Global -O3 vs selective -O3

Live peaks: 56.044 versus 56.071 ms. Global O3 changes profile FLASH code by +12,568 B and ITCM by +7,904 B; it remains a single-effect reference.
