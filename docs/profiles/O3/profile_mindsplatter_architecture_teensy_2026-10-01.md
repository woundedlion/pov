# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-01, **-O3**)

Shipping twin: [MindSplatter selective -O3](../shipping/profile_mindsplatter_architecture_teensy_2026-10-01.md).

Point-in-time architecture snapshot. Raw capture: [checkpoint1-mindsplatter-profile_o3.txt](../evidence/architecture_2026-10-01/checkpoint1-mindsplatter-profile_o3.txt); [capture provenance](../evidence/architecture_2026-10-01/checkpoint1-mindsplatter-profile_o3.provenance). Updates the current ranking alongside the preserved earlier capture [profile_mindsplatter_teensy_2026-08-26.md](profile_mindsplatter_teensy_2026-08-26.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile_o3`: `-O3`, `-ffast-math -fno-finite-math-only`; global optimization |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, source `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b` |
| Method | `HS_PROFILE`, 16-frame windows, 110 s capture; runtime 2–1737, setup frame 1 excluded; counter summaries 17–1728, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh MindSplatter profile_o3 110 16` |

Image size: `FLASH: code:90840, data:550952, headers:8448` /
`RAM1: variables:315232, code:54664, padding:10872, free:143520` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 993–1008 root cycles / 600 MHz
match the wall sum within **1.41 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are in the [evidence directory](../evidence/architecture_2026-10-01/README.md).

## Frame cadence

README cells: peak 🟢 54.539 (8), spilled 🟢 0/1736 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**22.132 ms/f**, peak render **54.539 ms**, spilled
**0/1736** live frames (**0.00%**).
Window summaries (17–1728): `msp_draw_particles` mean 19.898 ms/f;
worst window 43.821 ms/f (993–1008).

Startup setup render: **0.298 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
7.961 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Eight authored presets, 112-frame dwell and 48-frame linear parameter departure; the initial preset is held for 160 frames.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 993–1008, richest in this regime)

```text
frame                   62.266 ms 37.359 M 100.0%
  pov_preserve_half        0.137 ms  0.082 M   0.2% x1.0 136.58 us
  msp_draw_particles      43.821 ms 26.293 M  70.4%
    msp_particle_scan       43.816 ms 26.290 M  70.4%
      plot_ps_raster          34.750 ms 20.850 M  55.8% x563.0 61.72 us
      plot_ps_deferred         0.452 ms  0.271 M   0.7% x563.0 0.80 us
      plot_ps_gate             5.027 ms  3.016 M   8.1%
        plot_ps_cartesian_gate   0.637 ms  0.382 M   1.0% x1536.9 0.41 us
      plot_ps_tween            3.204 ms  1.922 M   5.1% x1536.9 2.08 us
  msp_particle_step        2.874 ms  1.725 M   4.6% x1.0 2874.37 us
  msp_timeline_step        0.055 ms  0.033 M   0.1% x1.0 54.68 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 84.52 us
  canvas_buffer_wait      15.294 ms  9.176 M  24.6% x1.0 15293.67 us
```

Wall min/avg/max = 50.920/62.265/73.045 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 945–960, richest in this regime)

```text
frame                   62.387 ms 37.432 M 100.0%
  pov_preserve_half        0.135 ms  0.081 M   0.2% x1.0 134.92 us
  msp_draw_particles      33.693 ms 20.216 M  54.0%
    msp_particle_scan       33.687 ms 20.212 M  54.0%
      plot_ps_raster          23.617 ms 14.170 M  37.9% x484.2 48.77 us
      plot_ps_deferred         1.742 ms  1.045 M   2.8% x484.2 3.60 us
      plot_ps_gate             4.493 ms  2.696 M   7.2%
        plot_ps_cartesian_gate   0.667 ms  0.400 M   1.1% x1657.1 0.40 us
      plot_ps_tween            3.436 ms  2.062 M   5.5% x1657.1 2.07 us
  msp_particle_step        5.931 ms  3.559 M   9.5% x1.0 5930.95 us
  msp_timeline_step        0.052 ms  0.031 M   0.1% x1.0 52.21 us
  canvas_clear             0.084 ms  0.051 M   0.1% x1.0 84.46 us
  canvas_buffer_wait      22.489 ms 13.493 M  36.0% x1.0 22489.14 us
```

Wall min/avg/max = 58.652/62.387/67.300 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 7 | Cube, friction 0.7465 | 10 / 9 | not recorded | 43.821 | 46.972 | 16.06 |
| 2 | Cube, strong well | 20 / 18 | not recorded | 35.715 | 38.081 | 15.72 |
| 6 | Dodecahedron | 10 / 9 | not recorded | 35.663 | 41.808 | 15.88 |
| 4 | Octahedron | 10 / 9 | not recorded | 25.558 | 28.161 | 16.05 |
| 8 | Cube, no Mobius | 10 / 9 | not recorded | 24.822 | 27.115 | 16.03 |
| 1 | Cube, friction 0.85 | 19 / 17 | not recorded | 24.796 | 27.041 | 15.98 |
| 3 | Cube, friction 0.9645 | 18 / 17 | not recorded | 23.551 | 25.811 | 15.99 |
| 5 | Tetrahedron, no Mobius | 10 / 9 | not recorded | 15.471 | 16.779 | 16.00 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner: preset 7 (Cube, friction 0.7465). Per-preset peaks span 20.411–54.539 ms.

| # | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 7 | Cube, friction 0.7465 | 54.539 | 0/159 (0.00%) |
| 2 | Cube, strong well | 48.413 | 0/318 (0.00%) |
| 6 | Dodecahedron | 46.617 | 0/159 (0.00%) |
| 4 | Octahedron | 30.677 | 0/159 (0.00%) |
| 8 | Cube, no Mobius | 30.438 | 0/159 (0.00%) |
| 3 | Cube, friction 0.9645 | 29.629 | 0/306 (0.00%) |
| 1 | Cube, friction 0.85 | 28.777 | 0/317 (0.00%) |
| 5 | Tetrahedron, no Mobius | 20.411 | 0/159 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. Particle, gate and raster calls count plot work, not blended pixels; dividing them into a per-pixel cost would misstate the measurement.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–1728; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.31/f 0.41/1.49/18.02 us 2.75% CPU
isr_pack        144.03/f 5.99/6.65/9.47 us 1.53% CPU
isr_dma_submit  144.03/f 0.60/0.93/1.22 us 0.21% CPU
```

- DMA submit averages 0.93 us/event; packing averages 6.65 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.50%, leaving approximately 59.689 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `msp_draw_particles`: 43.821 ms/f, 70.4% of the richest runtime window.
2. `msp_particle_step`: 2.874 ms/f, 4.6% of the richest runtime window.
3. `pov_preserve_half`: 0.137 ms/f, 0.2% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 22.129 ms/f,
peak 54.516 ms, spills 0/1736.
Candidate `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b`: mean render 22.132 ms/f
(+0.014% in these captures), peak 54.539 ms,
spills 0/1736. Baseline and candidate include 1736 and
1736 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses Plot point raster and culling HS_O3 regions; the global-O3 twin optimizes all compiled code.
- No ordered-cycle or transition-speed override was used; dwell and transition settings remain authored.
- Captured source is the exact committed SHA above. Build flags, source status and hashes are retained as evidence; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:22**.

Global -O3 versus selective -O3: mean render 22.132 versus 23.733 ms/f (ratio 1.072× for these samples). The O3 image adds +21,976 B FLASH code and +17,200 B ITCM. Live frame counts differ; compare spill fractions and peaks alongside these means.
