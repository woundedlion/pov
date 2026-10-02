# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-01, **selective -O3**)

Point-in-time architecture snapshot. Raw capture: checkpoint2-latticemelt-profile.txt (supporting artifact removed); capture provenance (supporting artifact removed). Updates the current ranking alongside the preserved earlier capture [profile_latticemelt_teensy_2026-09-28.md](profile_latticemelt_teensy_2026-09-28.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile`: `-Os`, `-ffast-math -fno-finite-math-only`; the composed shade and generated-palette HS_O3 functions; sphere curl-noise displacement |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, source `d25dd85dee17c28d2728f813e1ae00ae64690e1b` |
| Method | `HS_PROFILE`, 16-frame windows, 110 s capture; `HS_PROFILE_EPOCH_REVS=1200` (150 s epoch); runtime 2–1737, setup frame 1 excluded; counter summaries 17–1728, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh LatticeMelt profile 110 16 -D HS_PROFILE_EPOCH_REVS=1200` |

Image size: `FLASH: code:71600, data:156284, headers:8660` /
`RAM1: variables:315040, code:19960, padding:12808, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 1–16 root cycles / 600 MHz
match the wall sum within **1.06 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are no longer retained in the repository.

## Frame cadence

README cells: peak 🟢 36.583 (2), spilled 🟢 0/1736 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**36.120 ms/f**, peak render **36.583 ms**, spilled
**0/1736** live frames (**0.00%**).
Window summaries (17–1728): `fx_shader_draw` mean 32.448 ms/f;
worst window 32.602 ms/f (385–400).

Startup setup render: **68.510 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
25.917 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Two authored presets, 600-frame dwell and 480-frame eased parameter departure.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 385–400, richest in this regime)

```text
frame                   62.419 ms 37.451 M 100.0%
  pov_preserve_half        0.139 ms  0.083 M   0.2% x1.0 138.77 us
  fx_shader_draw          32.602 ms 19.561 M  52.2% x1.0 32602.18 us
  fx_prepare_frame         1.180 ms  0.708 M   1.9% x1.0 1179.78 us
  fx_advance               2.275 ms  1.365 M   3.6% x1.0 2275.28 us
  fx_timeline_step         0.088 ms  0.053 M   0.1% x1.0 88.08 us
  canvas_clear             0.090 ms  0.054 M   0.1% x1.0 90.18 us
  canvas_buffer_wait      25.998 ms 15.599 M  41.7% x1.0 25998.04 us
```

Wall min/avg/max = 62.231/62.419/62.578 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 1665–1680, richest in this regime)

```text
frame                   62.432 ms 37.459 M 100.0%
  pov_preserve_half        0.142 ms  0.085 M   0.2% x1.0 141.91 us
  fx_shader_draw          32.455 ms 19.473 M  52.0% x1.0 32455.24 us
  fx_prepare_frame         1.000 ms  0.600 M   1.6% x1.0 1000.48 us
  fx_advance               2.311 ms  1.386 M   3.7% x1.0 2310.63 us
  fx_timeline_step         0.089 ms  0.053 M   0.1% x1.0 88.80 us
  canvas_clear             0.090 ms  0.054 M   0.1% x1.0 89.97 us
  canvas_buffer_wait      26.300 ms 15.780 M  42.1% x1.0 26299.90 us
```

Wall min/avg/max = 62.158/62.431/62.674 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # (one-based marker order) | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 1 | open-curl | 40 / 39 | not recorded | 32.602 | 36.421 | 16.02 |
| 2 | dense-curl | 67 / 66 | not recorded | 32.510 | 36.066 | 16.02 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner (one-based marker order): preset 1 (open-curl). Per-preset peaks span 36.535–36.583 ms.

| # (one-based marker order) | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 1 | open-curl | 36.583 | 0/657 (0.00%) |
| 2 | dense-curl | 36.535 | 0/1079 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–1728; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.07/f 0.62/1.70/19.16 us 3.13% CPU
isr_pack        144.00/f 6.26/7.01/9.85 us 1.62% CPU
isr_dma_submit  144.00/f 0.58/0.94/11.54 us 0.22% CPU
```

- DMA submit averages 0.94 us/event; packing averages 7.01 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.96%, leaving approximately 59.402 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 32.602 ms/f, 52.2% of the richest runtime window.
2. `fx_advance`: 2.275 ms/f, 3.6% of the richest runtime window.
3. `fx_prepare_frame`: 1.180 ms/f, 1.9% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 36.374 ms/f,
peak 36.908 ms, spills 0/1736.
Candidate `d25dd85dee17c28d2728f813e1ae00ae64690e1b`: mean render 36.120 ms/f
(-0.698% in these captures), peak 36.583 ms,
spills 0/1736. Baseline and candidate include 1736 and
1736 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 981.633 us/f versus baseline 974.352 us/f. Shading averages 32.448 versus 32.722 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; sphere curl-noise displacement; the global-O3 twin optimizes all compiled code.
- Only the epoch is stretched to 150 s; dwell and transition settings remain authored.
- Captured source identity is recorded above; see the archive source-reachability note. Build flags, source status and hashes are not retained in the repository; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=LatticeMelt`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:38**.
