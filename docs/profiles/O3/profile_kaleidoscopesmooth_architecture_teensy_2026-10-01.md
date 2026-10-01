# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-10-01, **-O3**)

Shipping twin: [KaleidoscopeSmooth selective -O3](../shipping/profile_kaleidoscopesmooth_architecture_teensy_2026-10-01.md).

Point-in-time architecture snapshot. Raw capture: [final-kaleidoscopesmooth-profile_o3.txt](../evidence/architecture_2026-10-01/final-kaleidoscopesmooth-profile_o3.txt); [capture provenance](../evidence/architecture_2026-10-01/final-kaleidoscopesmooth-profile_o3.provenance). Updates the current ranking alongside the preserved earlier capture [profile_kaleidoscopesmooth_teensy_2026-08-26.md](profile_kaleidoscopesmooth_teensy_2026-08-26.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile_o3`: `-O3`, `-ffast-math -fno-finite-math-only`; global optimization |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, source `76ea80495cef34aa347431bb0aaabb5aa9d3ee55` |
| Method | `HS_PROFILE`, 16-frame windows, 260 s capture; `HS_PROFILE_EPOCH_REVS=2400` (300 s epoch); runtime 2–4138, setup frame 1 excluded; counter summaries 17–4128, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh KaleidoscopeSmooth profile_o3 260 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image size: `FLASH: code:85736, data:156616, headers:8528` /
`RAM1: variables:315040, code:28600, padding:4168, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 865–880 root cycles / 600 MHz
match the wall sum within **1.62 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are in the [evidence directory](../evidence/architecture_2026-10-01/README.md).

## Frame cadence

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**25.481 ms/f**, peak render **30.314 ms**, spilled
**0/4137** live frames (**0.00%**).
Window summaries (17–4128): `fx_shader_draw` mean 22.196 ms/f;
worst window 26.438 ms/f (865–880).

Startup setup render: **47.534 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 10 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
32.186 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Four authored presets, 600-frame dwell and 480-frame eased parameter departure.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 865–880, richest in this regime)

```text
frame                   62.374 ms 37.424 M 100.0%
  pov_preserve_half        0.145 ms  0.087 M   0.2% x1.0 145.28 us
  fx_shader_draw          26.438 ms 15.863 M  42.4% x1.0 26437.87 us
  fx_prepare_frame         0.727 ms  0.436 M   1.2% x1.0 726.60 us
  fx_advance               2.136 ms  1.281 M   3.4% x1.0 2135.63 us
  fx_timeline_step         0.152 ms  0.091 M   0.2% x1.0 151.58 us
  canvas_clear             0.084 ms  0.051 M   0.1% x1.0 84.37 us
  canvas_buffer_wait      32.654 ms 19.592 M  52.4% x1.0 32653.99 us
```

Wall min/avg/max = 61.601/62.373/63.121 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 593–608, richest in this regime)

```text
frame                   62.526 ms 37.516 M 100.0%
  pov_preserve_half        0.142 ms  0.085 M   0.2% x1.0 141.92 us
  fx_shader_draw          25.555 ms 15.333 M  40.9% x1.0 25554.73 us
  fx_prepare_frame         0.440 ms  0.264 M   0.7% x1.0 440.12 us
  fx_advance               1.961 ms  1.176 M   3.1% x1.0 1960.74 us
  fx_timeline_step         0.134 ms  0.081 M   0.2% x1.0 134.21 us
  canvas_clear             0.084 ms  0.051 M   0.1% x1.0 84.34 us
  canvas_buffer_wait      34.178 ms 20.507 M  54.7% x1.0 34178.29 us
```

Wall min/avg/max = 60.564/62.526/64.372 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 2 | direct-grid | 67 / 66 | not recorded | 26.438 | 29.720 | 16.03 |
| 1 | coupled-grid | 55 / 54 | not recorded | 25.575 | 28.554 | 16.03 |
| 3 | double-map | 68 / 67 | not recorded | 22.968 | 26.049 | 16.02 |
| 4 | stretched-grid | 67 / 66 | not recorded | 22.958 | 26.458 | 16.03 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner: preset 2 (direct-grid). Per-preset peaks span 27.012–30.314 ms.

| # | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 2 | direct-grid | 30.314 | 0/1079 (0.00%) |
| 1 | coupled-grid | 29.498 | 0/900 (0.00%) |
| 4 | stretched-grid | 27.171 | 0/1079 (0.00%) |
| 3 | double-map | 27.012 | 0/1079 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–4128; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.09/f 0.39/1.57/30.44 us 2.90% CPU
isr_pack        144.00/f 5.98/6.88/12.34 us 1.58% CPU
isr_dma_submit  144.00/f 0.58/0.93/10.87 us 0.22% CPU
```

- DMA submit averages 0.93 us/event; packing averages 6.88 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.70%, leaving approximately 59.565 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 26.438 ms/f, 42.4% of the richest runtime window.
2. `fx_advance`: 2.136 ms/f, 3.4% of the richest runtime window.
3. `fx_prepare_frame`: 0.727 ms/f, 1.2% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 25.542 ms/f,
peak 30.299 ms, spills 0/4136.
Candidate `76ea80495cef34aa347431bb0aaabb5aa9d3ee55`: mean render 25.481 ms/f
(-0.238% in these captures), peak 30.314 ms,
spills 0/4137. Baseline and candidate include 4136 and
4137 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 853.200 us/f versus baseline 847.679 us/f. Shading averages 22.196 versus 22.266 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; the global-O3 twin optimizes all compiled code.
- Only the epoch is stretched to 300 s; dwell and transition settings remain authored.
- Captured source is the exact committed SHA above. Build flags, source status and hashes are retained as evidence; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 11:19**.

Global -O3 versus selective -O3: mean render 25.481 versus 26.830 ms/f (ratio 1.053× for these samples). The O3 image adds +15,048 B FLASH code and +8,624 B ITCM. Live frame counts differ; compare spill fractions and peaks alongside these means.
