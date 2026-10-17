# GridSpace on-device profile — Teensy 4.0, segmented mode (2026-10-17, **selective -O3**)

Point-in-time snapshot (regenerate with the Reproduce command below).
Raw captures: `build/prof/gridspace_temporal_before_20261017.log` and
`build/prof/gridspace_temporal_after_20261017.log`, captured 2026-10-17 14:43
and 2026-10-17 14:47 America/Los_Angeles on COM4. Replaces profile_gridspace_teensy_2026-10-06.md.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, COM4, flywheel + DMA ISRs live |
| Image | Shipping `profile` environment: `-Os`, newlib-nano, selective `HS_O3` shading, `HS_HOT_FLASH_MEMBER` column traversal |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | GridSpace 288×144, single-entry playlist, clean candidate `5cbdd7d91c493c5af7d84d5ecadf61d1bcd4cc22` |
| Baseline | Same board and flags, clean `0f8a086aa6f76b6f415f90d13ffe1cc89ca29a7b`, captured first |
| Method | One 70 s capture per image, 32-frame windows, no epoch stretch or dwell compression |
| Reproduce | Set `HS_PROFILE_TREE` to the selected checkout and `HS_TEENSY_PORT=COM4`; run `bash tools/profile_one.sh GridSpace profile 70 32` |

Candidate image size (`profile`, this effect only): `FLASH: code:72604, data:155260, headers:8676` / `RAM1: variables:312992, code:19976, padding:12792, free:178528` / `RAM2: variables:520064, free:4224`.
Baseline image size: `FLASH: code:69324, data:154376, headers:8744` / `RAM1: variables:312992, code:19624, padding:13144, free:178528` / `RAM2: variables:520064, free:4224`.

The full shipping `phantasm` size/layout gate passes:

| Memory | Baseline bytes | Temporal bytes | Delta bytes |
|---|---:|---:|---:|
| RAM1 code | 179,640 | 179,816 | +176 |
| RAM1 variables | 312,736 | 312,736 | 0 |
| FLASH data | 851,396 | 851,388 | −8 |
| FLASH total | 1,496,060 | 1,498,108 | +2,048 |
| GridSpace persistent arena, host probe | 13,584 | 25,536 | +11,952 |

The arena comparison uses `arena_measure` at 288×144 for 120 frames. These
are host high-water marks, not device allocations; the host ABI has larger
pointers. The additional retained state lives inside the existing global arena,
so RAM1's statically reserved variable footprint is unchanged.

Both logs pass `tools/parse_profile.py` validation: 34 complete windows,
1,088 frames, monotonic numbering, correct effect/resolution, and no epoch
reset. For frames 1–32, candidate root cycles 1,182,914,060 / 600 give
1,971,523.433 us versus wall sum 1,971,530 us: **3.3 ppm**. Baseline agrees
within **2.4 ppm**.

The wrapper attested GCC 15.2.1 and matching profile/shipping compiler and ABI.
Candidate profile ELF SHA-256 is
`6964f7dc68cb825591e80f75361fa6693e7b0c444c932abb4dffb0fe52866311`;
baseline is `9414a42b158e89818795c08f3f6d3d2aa31d851340361ddcbc99e95421200927`.
Raw logs, build logs, environment dumps, maps, ELFs, and provenance are retained
under C:/work/temp/grid-temporal-profile-20261017, in separate baseline and
candidate artifact directories.

## Frame cadence

**Initialization/priming draw:** **62.977 ms**, versus baseline **40.398 ms**.
`hardware/pov_segment_frame.h` deliberately draws a full-width 288×72 band
(20,736 pixels) in `prime_segment_frame()`. This runs inside `rebuild_effect()`
before `handoff.publish()`, while the replacement effect is not yet live.
It is construction cost, not a missed live display frame.

The unfiltered parser footer reports one of 1,088 draw calls above 62.5 ms.
That arithmetic threshold includes the priming draw; it must not be described
as a dropped display frame. No live render in this capture exceeded the budget.

**Steady display, frames 2–1,088:** render averages **36.130 ms**,
peak **36.615 ms** at frame 399, and **0/1,087 spilled (0.00%)**.
The ranked indexes retain their existing convention of excluding setup frame 1.
Shader scope averages **29.657 ms/frame** over complete steady
windows (frames 33–1,088); its worst window is **29.770 ms/frame**,
frames 65–96.

A half-turn display window is 62.5 ms. Each steady frame shades one 144×72
quadrant (10,368 pixels), keeping **16 fps** with 25.885 ms between the worst
observed steady render and the deadline. `canvas_buffer_wait` is the round-up
idle until the display flip.

### Same-board baseline comparison

| Work | Baseline ms/frame | Temporal ms/frame | Delta ms | Delta |
|---|---:|---:|---:|---:|
| Render | 24.222 | 36.130 | +11.907 | +49.2% |
| Pixel shading + column preparation | 17.838 | 29.657 | +11.819 | +66.3% |
| Frame preparation | 3.762 | 3.811 | +0.049 | +1.3% |
| Animation advancement | 2.267 | 2.291 | +0.024 | +1.1% |
| Steady peak render | 24.642 | 36.615 | +11.973 | +48.6% |

Render means use per-frame rows 2–1,088. Scope means exclude the entire startup
window because its counters cannot separate individual frames. The added
**11.907 ms/frame (49.2%)** remains within the display window.
Almost all the increase lies inside `fx_shader_draw`, which includes column
preparation, traversal, and the two endpoint color evaluations. This capture
does not isolate those three contributions. Baseline also had zero steady spills.

## Phase-by-phase readout

Phase schedule: one steady preset with continuous camera, warp, and color
animation. Initialization is accounted for separately above.

### Steady peak window (frames 385–416)

```
frame                     62.41 ms  37.45 Mcyc 100.0%
  pov_preserve_half       140.2 us   84.1 kcyc   0.2% x1 140.2 us/call
  fx_shader_draw          29.60 ms  17.76 Mcyc  47.4% x1 29604.5 us/call
  fx_prepare_frame         3.96 ms   2.37 Mcyc   6.3% x1 3955.9 us/call
  fx_advance               2.27 ms   1.36 Mcyc   3.6% x1 2267.8 us/call
  fx_timeline_step         71.7 us   43.0 kcyc   0.1% x1 71.7 us/call
  canvas_clear             85.6 us   51.4 kcyc   0.1% x1 85.6 us/call
  canvas_buffer_wait      26.22 ms  15.73 Mcyc  42.0% x1 26215.8 us/call
```

Wall min/avg/max = 61.611/62.414/63.092 ms. The wall variation includes display
wait and does not indicate render spills. Shading occupies about 47% of the
window; the remaining sync wait averages 26.216 ms in this window. Retaining
the current endpoint happens outside the named shading and preparation scopes
and remains included in total render.

### Per-pixel figures

There is no `filter_blend` counter: this shader writes pixels directly.
The steady shading average over 10,368 output samples is
**1716 cycles/output pixel**, versus
1032 before the change.
These inclusive figures include column preparation and ISR interference;
they are not isolated pixel-shader instruction counts.

## Column-ISR / DMA marshaling cost

Peak-window ISR snapshot; rates are calls per 32 rendered frames, and CPU
shares use the separately logged 2,000,234 us ISR window:

```
counter          calls/frame  min/avg/max per call     CPU
isr_wake        1152.22/frame  0.62/2.02/27.92 us  3.72%
isr_pack         144.00/frame  8.39/9.19/12.17 us  2.11%
isr_dma_submit   144.00/frame  0.64/0.94/5.25 us  0.21%
```

- DMA submission averages 0.941 us/call; packing averages 9.188 us/call.
- The 72-LED image plus black strobe frame occupies 600 bytes. At the configured
  24 MHz SPI rate, the driver budgets 230 us including byte framing, within
  the approximately 434 us column period. The transfer is asynchronous.
- Inclusive `isr_wake` consumes 3.72% of the snapshot interval, corresponding
  to about 60.18 ms of foreground CPU time per 62.5 ms window. Pack and submit
  run inside wake and are not added to it. The measured render already includes
  ISR interference; no speedup is required for the observed steady regime.
  DMA-completion and other interrupts are not measured by these counters.

## Summary ranking

README cells: peak 🟢 36.62, spilled 🟢 0/1087 (0.00%).

1. `fx_shader_draw` — 47.4% of the steady peak window, 29.605 ms/frame.
2. `fx_prepare_frame` — 6.3%, 3.956 ms/frame.
3. `fx_advance` — 3.6%, 2.268 ms/frame.
4. `pov_preserve_half` — 0.2%, 0.140 ms/frame.

Temporal sampling adds roughly 11.9 ms/frame while preserving the observed
steady 16 fps cadence. This is a same-board firmware comparison; no native or
WASM timing comparison was collected. The simulator was installed before
device capture, and both images were measured through the supported locked
wrapper on the same COM4 board.

## Caveats

- CYCCNT runs through interrupts, so all foreground scopes include ISR time.
- Frame 1 is the full-width priming draw before effect publication. It is
  measured and reported as initialization cost, then excluded from live-frame
  peak/spill statistics under the archive's existing convention.
- There are no per-pixel profiling scopes or `filter_blend` parenting effects
  in this capture. Adding per-pixel instrumentation would alter these costs.
- Selective-O3 shipping behavior includes the hot flash column traversal and
  the existing hot pullback shading/color operations; the global-O3 image was
  not remeasured.
- No dwell compression, epoch override, or uncommitted source state was used.
- Endpoint interpolation samples retained frame states and color resources;
  this timing run does not measure visual seam quality or model the simulator.
- One sequential baseline/candidate capture measures this 70-second trajectory;
  it does not establish a worst-case bound for every GUI parameter setting.

## Harness

`targets/Profile/Profile.ino` with `HS_PROFILE_TARGET=GridSpace` and
`HS_PROFILE_WINDOW=32`. `just profile GridSpace 70` runs the supported locked
build, flash, and capture wrapper. Pin `HS_PROFILE_TREE` and `HS_TEENSY_PORT`
when reproducing a particular checkout and board; validate the result with
`python tools/parse_profile.py <log> validate`.
