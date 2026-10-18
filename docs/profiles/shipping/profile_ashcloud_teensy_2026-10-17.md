# AshCloud on-device profile — Teensy 4.0, segmented mode (2026-10-17, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile AshCloud`).
Raw capture: `C:/work/temp/rollout-ashcloud.log`, captured on COM4 at
2026-10-17 17:12 America/Los_Angeles. Replaces the October 6 shipping report.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | Shipping `profile` environment: `-Os`, newlib-nano, selective hot pullback shading/color regions and hot flash column traversal |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | AshCloud 288×144, single-entry playlist, temporal sampling enabled |
| Source | Clean `282d28354391cb71e2024942a5904aafc4e18ec5` |
| Method | 70 s capture, 32-frame windows; no epoch stretch or dwell compression |
| Reproduce | Pin `HS_PROFILE_TREE` and `HS_TEENSY_PORT=COM4`; run `bash tools/profile_one.sh AshCloud profile 70 32` |

Profile image size:
`FLASH: code:86772, data:154736, headers:8344` /
`RAM1: variables:313024, code:22184, padding:10584, free:178496` /
`RAM2: variables:520064, free:4224`.

The full shipping `phantasm` size/layout gate also passes. This is the memory
change for the entire 17-effect rollout, not AshCloud alone:

| Memory | Optimized GridSpace-only bytes | Full rollout bytes | Delta bytes |
|---|---:|---:|---:|
| RAM1 code | 181,704 | 186,104 | +4,400 |
| RAM1 variables | 312,736 | 312,736 | 0 |
| FLASH data | 851,388 | 851,388 | 0 |
| RAM2 variables | 520,064 | 520,064 | 0 |
| Stack allowance | 14,944 | 14,944 | 0 |

ITCM remains within its existing 196,608-byte allocation. Retained temporal
resources use the existing reserved arena; the fixed RAM1 allocation did not grow.
All ten firmware environments pass their build gates. The cold verbose warning
check covers 40/40 first-party translation units with zero first-party warnings.

Parser validation passes: 34 complete windows, 1,088 frames, correct effect and
resolution, monotonic numbering, and no epoch reset. For frames 545–576, root
cycles 1,197,616,413 / 600 give 1,996,027.355 us versus a 1,996,030 us wall sum:
**1.3 ppm** agreement.

The wrapper attested GCC 15.2.1 and matching profile/shipping compiler and ABI.
Profile ELF SHA-256:
`3acc34f933a1204b37eb6077d78ce66a26c31c0ba82beee6bd70252795999ff6`.
Raw logs, build logs, environment dumps, ELF files, and provenance are retained
locally under `C:/work/temp`, with immutable image artifacts in
`artifacts/rollout-ashcloud_282d28354391_3acc34f933a1`.

## Frame cadence

**Initialization/priming draw: 71.809 ms.** The driver primes a full-width
288×72 band before publishing the effect. This construction cost is excluded
from live-frame statistics; it is not a missed live display frame.

**Live frames 2–1,088:** render averages **43.581 ms**, peaks at **46.270 ms**
on frame 534, and spills **0/1,087 frames (0.00%)**. A half-turn is 62.5 ms;
each live frame emits one 144×72 quadrant (10,368 pixels). The observed cadence
remains **16 fps**, with **16.230 ms** between the worst live render and the
deadline. `canvas_buffer_wait` is the idle time until the display flip.

The unfiltered parser footer includes setup frame 1: peak 71.809 ms and
1/1,088 spills. The ranked indexes exclude that setup draw. Complete steady
windows (frames 33–1,088) average **37.034 ms/frame** in `fx_shader_draw`;
the worst shading window is **39.143 ms/frame**, frames 545–576.

The October 6 capture peaked at 42.67 ms with zero live spills, but used
`e2f5b0a3d` on COM3. This is not a matched temporal-off baseline, so the
difference does not isolate temporal sampling overhead.

## Phase-by-phase readout

Phase schedule: one continuously animated preset; no preset markers or cycle
compression. Initialization is reported separately above.

### Steady peak window (frames 513–544)

```
frame                    62.391 ms  37.435 Mcyc 100.0%
  pov_preserve_half     140.47 us   84.30 kcyc   0.2% x1 140.5 us/call
  fx_shader_draw         38.803 ms  23.282 Mcyc  62.2% x1 38803.4 us/call
    fx_prepare_column   604.50 us  362.71 kcyc   1.0% x146 4.14 us/call
  fx_prepare_frame        3.698 ms   2.219 Mcyc   5.9% x1 3698.3 us/call
  fx_advance              2.308 ms   1.385 Mcyc   3.7% x1 2307.6 us/call
  fx_timeline_step       67.81 us   40.69 kcyc   0.1% x1 67.8 us/call
  canvas_clear           85.72 us   51.43 kcyc   0.1% x1 85.7 us/call
  canvas_buffer_wait     16.985 ms  10.191 Mcyc  27.2% x1 16984.6 us/call
```

Wall min/avg/max = 60.448/62.391/64.294 ms. Wall variation includes sync wait
and is not render spill. Column preparation is included in shading; its 146
calls include border columns. Retaining endpoint resources lies outside the
named shading/preparation scopes and is included in total render. Across all
complete steady windows, column preparation averages 0.621 ms/frame.

AshCloud uses the shared prepared temporal clock/rotation and color sampling
paths. Its surface displacement and kaleidoscope lens precede projection, so
this pipeline does not use GridSpace's factored Rotate-to-Project column basis.

### Per-pixel figures

There is no `filter_blend` counter: this shader writes pixels directly.
Steady shading averages **2,143 cycles per nominal output pixel**, dividing
the inclusive scope by 10,368 emitted pixels. This includes column preparation,
border traversal, and ISR interference; it is not an isolated instruction count.

## Column-ISR / DMA marshaling cost

Peak-window rates are calls per 32 rendered frames. CPU shares use the separately
logged 2,000,255 us ISR interval:

```
counter          calls/frame  min/avg/max per call     CPU
isr_wake        1152.22/frame  0.61/2.03/19.95 us  3.74%
isr_pack         144.00/frame  8.86/9.57/12.40 us  2.20%
isr_dma_submit   144.00/frame  0.63/0.94/6.82 us   0.21%
```

- DMA submission averages 0.940 us/call; packing averages 9.566 us/call.
- The 72-LED image plus black strobe frame occupies 600 bytes. At the configured
  24 MHz SPI rate, the driver budgets 230 us including byte framing, within the
  approximately 434 us column period. The transfer is asynchronous.
- Inclusive wake consumes 3.74%, leaving approximately 60.16 ms of foreground
  CPU time per 62.5 ms window. Pack and submit run inside wake and are not added
  to it. Render measurements already include ISR interference; no speedup is
  required to sustain the observed cadence. DMA-completion and other interrupts
  are not measured by these counters.

## Summary ranking

README cells: peak 🟢 46.27, spilled 🟢 0/1087 (0.00%).

1. `fx_shader_draw` — 62.2% of the peak window, 38.803 ms/frame, including
   0.605 ms of column preparation.
2. `fx_prepare_frame` — 5.9%, 3.698 ms/frame.
3. `fx_advance` — 3.7%, 2.308 ms/frame.
4. `pov_preserve_half` — 0.2%, 0.140 ms/frame.

The newly enabled temporal path remains below the live deadline throughout
this capture. Native/WASM checks validate functionality; their timings are not
substituted for this on-device measurement.

## Caveats

- CYCCNT runs through interrupts, so foreground scopes include ISR time.
- Frame 1 primes the effect before publication and is initialization cost.
- There are no per-pixel scopes. Column-preparation instrumentation remains
  included in the measured image.
- The capture used a clean source commit, without dwell compression or an epoch
  override. One 70-second trajectory does not bound all GUI parameter settings.
- Historical shipping and global-O3 captures use different source and board
  contexts. The O3 report and its paired code-size deltas were not refreshed.
- This device capture validates AshCloud; it does not establish timing for all
  17 newly enabled effects or measure perceptual seam visibility.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=AshCloud` and
`HS_PROFILE_WINDOW=32`. `just profile AshCloud 70` invokes the supported locked
build, flash, and capture wrapper. Pin `HS_PROFILE_TREE` and `HS_TEENSY_PORT`
for the selected source and board. Validate with
`python tools/parse_profile.py <log> validate`.
