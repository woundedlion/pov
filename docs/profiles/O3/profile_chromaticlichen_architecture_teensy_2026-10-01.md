# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-10-01, **-O3**)

Shipping twin: [ChromaticLichen selective -O3](../shipping/profile_chromaticlichen_architecture_teensy_2026-10-01.md).

Point-in-time architecture snapshot. Raw capture: checkpoint2-chromaticlichen-profile_o3.txt (supporting artifact removed); capture provenance (supporting artifact removed). Updates the current ranking alongside the preserved earlier capture [profile_chromaticlichen_teensy_2026-08-26.md](profile_chromaticlichen_teensy_2026-08-26.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile_o3`: `-O3`, `-ffast-math -fno-finite-math-only`; global optimization |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, source `d25dd85dee17c28d2728f813e1ae00ae64690e1b` |
| Method | `HS_PROFILE`, 32-frame windows, 70 s capture; runtime 2–1097, setup frame 1 excluded; counter summaries 33–1088, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh ChromaticLichen profile_o3 70 32` |

Image size: `FLASH: code:84328, data:154936, headers:8544` /
`RAM1: variables:315040, code:28456, padding:4312, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 1–32 root cycles / 600 MHz
match the wall sum within **0.32 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are no longer retained in the repository.

## Frame cadence

README cells: peak 🟢 35.218, spilled 🟢 0/1096 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**34.600 ms/f**, peak render **35.218 ms**, spilled
**0/1096** live frames (**0.00%**).
Window summaries (33–1088): `fx_shader_draw` mean 28.597 ms/f;
worst window 28.796 ms/f (961–992).

Startup setup render: **61.387 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
27.282 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: One authored preset; no ownership-changing cycle.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 961–992, richest in this regime)

```text
frame                   62.430 ms 37.458 M 100.0%
  pov_preserve_half        0.142 ms  0.085 M   0.2% x1.0 141.61 us
  fx_shader_draw          28.796 ms 17.278 M  46.1% x1.0 28796.00 us
  fx_prepare_frame         3.860 ms  2.316 M   6.2% x1.0 3860.26 us
  fx_advance               1.946 ms  1.167 M   3.1% x1.0 1945.61 us
  fx_timeline_step         0.072 ms  0.043 M   0.1% x1.0 72.04 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 84.59 us
  canvas_buffer_wait      27.513 ms 16.508 M  44.1% x1.0 27513.08 us
```

Wall min/avg/max = 62.183/62.429/62.765 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 33–1088; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.07/f 0.36/1.56/26.43 us 2.87% CPU
isr_pack        144.00/f 5.98/6.66/9.36 us 1.54% CPU
isr_dma_submit  144.00/f 0.58/0.93/11.59 us 0.22% CPU
```

- DMA submit averages 0.93 us/event; packing averages 6.66 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.62%, leaving approximately 59.613 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 28.796 ms/f, 46.1% of the richest runtime window.
2. `fx_prepare_frame`: 3.860 ms/f, 6.2% of the richest runtime window.
3. `fx_advance`: 1.946 ms/f, 3.1% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 34.860 ms/f,
peak 35.415 ms, spills 0/1096.
Candidate `d25dd85dee17c28d2728f813e1ae00ae64690e1b`: mean render 34.600 ms/f
(-0.746% in these captures), peak 35.218 ms,
spills 0/1096. Baseline and candidate include 1096 and
1096 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 3671.752 us/f versus baseline 3664.534 us/f. Shading averages 28.597 versus 28.852 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; noise-contour source; the global-O3 twin optimizes all compiled code.
- No ordered-cycle or transition-speed override was used; dwell and transition settings remain authored.
- Captured source identity is recorded above; see the archive source-reachability note. Build flags, source status and hashes are not retained in the repository; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:51**.

Global -O3 versus selective -O3: mean render 34.600 versus 35.274 ms/f (ratio 1.019× for these samples). The O3 image adds +15,264 B FLASH code and +8,672 B ITCM. Live frame counts differ; compare spill fractions and peaks alongside these means.
