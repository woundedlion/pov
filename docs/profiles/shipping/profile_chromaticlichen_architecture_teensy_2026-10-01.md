# ChromaticLichen on-device profile — Teensy 4.0, segmented mode (2026-10-01, **selective -O3**)

Point-in-time architecture snapshot. Raw capture: [checkpoint2-chromaticlichen-profile.txt](../evidence/architecture_2026-10-01/checkpoint2-chromaticlichen-profile.txt); [capture provenance](../evidence/architecture_2026-10-01/checkpoint2-chromaticlichen-profile.provenance). Updates the current ranking alongside the preserved earlier capture [profile_chromaticlichen_teensy_2026-09-28.md](profile_chromaticlichen_teensy_2026-09-28.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile`: `-Os`, `-ffast-math -fno-finite-math-only`; the composed shade and generated-palette HS_O3 functions; noise-contour source |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | ChromaticLichen 288×144, single-entry playlist, source `d25dd85dee17c28d2728f813e1ae00ae64690e1b` |
| Method | `HS_PROFILE`, 32-frame windows, 70 s capture; runtime 2–1097, setup frame 1 excluded; counter summaries 33–1088, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh ChromaticLichen profile 70 32` |

Image size: `FLASH: code:69064, data:155048, headers:8336` /
`RAM1: variables:315040, code:19784, padding:12984, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 1–32 root cycles / 600 MHz
match the wall sum within **3.00 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are in the [evidence directory](../evidence/architecture_2026-10-01/README.md).

## Frame cadence

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**35.274 ms/f**, peak render **35.773 ms**, spilled
**0/1096** live frames (**0.00%**).
Window summaries (33–1088): `fx_shader_draw` mean 28.930 ms/f;
worst window 29.256 ms/f (801–832).

Startup setup render: **62.118 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
26.727 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: One authored preset; no ownership-changing cycle.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 801–832, richest in this regime)

```text
frame                   62.458 ms 37.475 M 100.0%
  pov_preserve_half        0.138 ms  0.083 M   0.2% x1.0 138.01 us
  fx_shader_draw          29.256 ms 17.554 M  46.8% x1.0 29256.28 us
  fx_prepare_frame         3.874 ms  2.324 M   6.2% x1.0 3874.11 us
  fx_advance               1.958 ms  1.175 M   3.1% x1.0 1958.05 us
  fx_timeline_step         0.070 ms  0.042 M   0.1% x1.0 69.55 us
  canvas_clear             0.086 ms  0.052 M   0.1% x1.0 86.31 us
  canvas_buffer_wait      27.046 ms 16.227 M  43.3% x1.0 27045.55 us
```

Wall min/avg/max = 62.189/62.457/62.680 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 33–1088; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.07/f 0.59/1.68/27.74 us 3.10% CPU
isr_pack        144.00/f 6.25/6.87/11.76 us 1.58% CPU
isr_dma_submit  144.00/f 0.60/0.93/11.03 us 0.22% CPU
```

- DMA submit averages 0.93 us/event; packing averages 6.87 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.90%, leaving approximately 59.438 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 29.256 ms/f, 46.8% of the richest runtime window.
2. `fx_prepare_frame`: 3.874 ms/f, 6.2% of the richest runtime window.
3. `fx_advance`: 1.958 ms/f, 3.1% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 35.331 ms/f,
peak 35.875 ms, spills 0/1096.
Candidate `d25dd85dee17c28d2728f813e1ae00ae64690e1b`: mean render 35.274 ms/f
(-0.160% in these captures), peak 35.773 ms,
spills 0/1096. Baseline and candidate include 1096 and
1096 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 4038.536 us/f versus baseline 4031.874 us/f. Shading averages 28.930 versus 28.961 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; noise-contour source; the global-O3 twin optimizes all compiled code.
- No ordered-cycle or transition-speed override was used; dwell and transition settings remain authored.
- Captured source is the exact committed SHA above. Build flags, source status and hashes are retained as evidence; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=ChromaticLichen`, `HS_PROFILE_WINDOW=32`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:40**.
