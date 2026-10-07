# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-01, **selective -O3**)

Point-in-time architecture snapshot. Raw capture: checkpoint1-mindsplatter-profile.txt (supporting artifact removed); capture provenance (supporting artifact removed). Superseded in the ranking by [profile_mindsplatter_teensy_2026-10-06.md](profile_mindsplatter_teensy_2026-10-06.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile`: `-Os`, `-ffast-math -fno-finite-math-only`; Plot point raster and culling HS_O3 regions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, source `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b` |
| Method | `HS_PROFILE`, 16-frame windows, 110 s capture; runtime 2–1737, setup frame 1 excluded; counter summaries 17–1728, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh MindSplatter profile 110 16` |

Image size: `FLASH: code:68864, data:551076, headers:8796` /
`RAM1: variables:315200, code:37464, padding:28072, free:143552` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 993–1008 root cycles / 600 MHz
match the wall sum within **0.34 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are no longer retained in the repository.

## Frame cadence

README cells: peak 🟢 54.457 (8), spilled 🟢 0/1736 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**23.733 ms/f**, peak render **54.457 ms**, spilled
**0/1736** live frames (**0.00%**).
Window summaries (17–1728): `msp_draw_particles` mean 21.277 ms/f;
worst window 47.302 ms/f (993–1008).

Startup setup render: **0.286 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
8.043 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Eight authored presets, 112-frame dwell and 48-frame linear parameter departure; the initial preset is held for 160 frames.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 993–1008, richest in this regime)

```text
frame                   62.162 ms 37.297 M 100.0%
  pov_preserve_half        0.137 ms  0.082 M   0.2% x1.0 137.22 us
  msp_draw_particles      47.302 ms 28.381 M  76.1%
    msp_particle_scan       47.295 ms 28.377 M  76.1%
      plot_ps_raster          34.239 ms 20.543 M  55.1% x576.0 59.44 us
      plot_ps_deferred         0.568 ms  0.341 M   0.9% x576.0 0.99 us
      plot_ps_gate             6.768 ms  4.061 M  10.9%
        plot_ps_cartesian_gate   0.728 ms  0.437 M   1.2% x1537.0 0.47 us
      plot_ps_tween            4.924 ms  2.954 M   7.9% x1537.0 3.20 us
  msp_particle_step        2.942 ms  1.765 M   4.7% x1.0 2942.35 us
  msp_timeline_step        0.070 ms  0.042 M   0.1% x1.0 69.92 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 84.76 us
  canvas_buffer_wait      11.625 ms  6.975 M  18.7% x1.0 11624.63 us
```

Wall min/avg/max = 59.141/62.161/64.809 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 945–960, richest in this regime)

```text
frame                   62.854 ms 37.712 M 100.0%
  pov_preserve_half        0.135 ms  0.081 M   0.2% x1.0 135.30 us
  msp_draw_particles      40.518 ms 24.311 M  64.5%
    msp_particle_scan       40.511 ms 24.307 M  64.5%
      plot_ps_raster          25.667 ms 15.400 M  40.8% x553.1 46.41 us
      plot_ps_deferred         2.275 ms  1.365 M   3.6% x553.1 4.11 us
      plot_ps_gate             6.457 ms  3.874 M  10.3%
        plot_ps_cartesian_gate   0.763 ms  0.458 M   1.2% x1657.1 0.46 us
      plot_ps_tween            5.269 ms  3.161 M   8.4% x1657.1 3.18 us
  msp_particle_step        7.558 ms  4.535 M  12.0% x1.0 7558.45 us
  msp_timeline_step        0.056 ms  0.033 M   0.1% x1.0 55.67 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 84.79 us
  canvas_buffer_wait      14.498 ms  8.699 M  23.1% x1.0 14498.40 us
```

Wall min/avg/max = 59.042/62.854/69.802 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # (one-based marker order) | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 7 | Cube, friction 0.7465 | 10 / 9 | not recorded | 47.302 | 50.537 | 16.09 |
| 6 | Dodecahedron | 10 / 9 | not recorded | 39.943 | 47.719 | 15.98 |
| 2 | Cube, strong well | 20 / 18 | not recorded | 29.946 | 32.448 | 16.14 |
| 8 | Cube, no Mobius | 10 / 9 | not recorded | 29.892 | 32.230 | 15.88 |
| 1 | Cube, friction 0.85 | 19 / 17 | not recorded | 27.726 | 30.052 | 15.93 |
| 3 | Cube, friction 0.9645 | 18 / 17 | not recorded | 27.359 | 29.689 | 16.03 |
| 4 | Octahedron | 10 / 9 | not recorded | 27.309 | 30.535 | 16.05 |
| 5 | Tetrahedron, no Mobius | 10 / 9 | not recorded | 16.689 | 18.469 | 16.07 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner (one-based marker order): preset 7 (Cube, friction 0.7465). Per-preset peaks span 22.224–54.457 ms.

| # (one-based marker order) | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 7 | Cube, friction 0.7465 | 54.457 | 0/159 (0.00%) |
| 6 | Dodecahedron | 51.354 | 0/159 (0.00%) |
| 8 | Cube, no Mobius | 44.393 | 0/159 (0.00%) |
| 3 | Cube, friction 0.9645 | 36.706 | 0/306 (0.00%) |
| 2 | Cube, strong well | 35.114 | 0/318 (0.00%) |
| 1 | Cube, friction 0.85 | 34.123 | 0/317 (0.00%) |
| 4 | Octahedron | 33.511 | 0/159 (0.00%) |
| 5 | Tetrahedron, no Mobius | 22.224 | 0/159 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. Particle, gate and raster calls count plot work, not blended pixels; dividing them into a per-pixel cost would misstate the measurement.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–1728; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.33/f 0.52/1.66/20.72 us 3.06% CPU
isr_pack        144.03/f 6.25/6.87/9.86 us 1.58% CPU
isr_dma_submit  144.03/f 0.61/0.94/1.24 us 0.22% CPU
```

- DMA submit averages 0.94 us/event; packing averages 6.87 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Measured inclusive flywheel ISR (`isr_wake`) share is 3.06%, leaving approximately 60.588 ms of CPU-only work per 62.5 ms window; `isr_pack` and `isr_dma_submit` are nested inside it, not additional. DMA-completion and other interrupts are not measured. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `msp_draw_particles`: 47.302 ms/f, 76.1% of the richest runtime window.
2. `msp_particle_step`: 2.942 ms/f, 4.7% of the richest runtime window.
3. `pov_preserve_half`: 0.137 ms/f, 0.2% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 24.342 ms/f,
peak 55.278 ms, spills 0/1736.
Candidate `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b`: mean render 23.733 ms/f
(-2.501% in these captures), peak 54.457 ms,
spills 0/1736. Baseline and candidate include 1736 and
1736 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses Plot point raster and culling HS_O3 regions; the global-O3 twin optimizes all compiled code.
- No ordered-cycle or transition-speed override was used; dwell and transition settings remain authored.
- Captured source identity is recorded above; see the archive source-reachability note. Build flags, source status and hashes are not retained in the repository; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:15**.
