# KaleidoscopeSmooth on-device profile — Teensy 4.0, segmented mode (2026-10-01, **selective -O3**)

Point-in-time architecture snapshot. Raw capture: [aligned-board4-kaleidoscopesmooth-profile.txt](../evidence/architecture_2026-10-01/aligned-board4-kaleidoscopesmooth-profile.txt); [capture provenance](../evidence/architecture_2026-10-01/aligned-board4-kaleidoscopesmooth-profile.provenance). Updates the current ranking alongside the preserved earlier capture [profile_kaleidoscopesmooth_teensy_2026-09-28.md](profile_kaleidoscopesmooth_teensy_2026-09-28.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile`: `-Os`, `-ffast-math -fno-finite-math-only`; the composed shade and generated-palette HS_O3 functions |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | KaleidoscopeSmooth 288×144, single-entry playlist, source `0156d0d7490355ea99cab406f4704a128a55ccbf` plus [retained source patch](../evidence/architecture_2026-10-01/aligned-board4-kaleidoscopesmooth-profile_source.diff) |
| Method | `HS_PROFILE`, 16-frame windows, 260 s capture; `HS_PROFILE_EPOCH_REVS=2400` (300 s epoch); runtime 2–4137, setup frame 1 excluded; counter summaries 17–4128, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh KaleidoscopeSmooth profile 260 16 -D HS_PROFILE_EPOCH_REVS=2400` |

Image size: `FLASH: code:70688, data:156636, headers:8196` /
`RAM1: variables:315040, code:19976, padding:12792, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 3953–3968 root cycles / 600 MHz
match the wall sum within **4.46 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are in the [evidence directory](../evidence/architecture_2026-10-01/README.md).

## Frame cadence

README cells: peak 🟢 32.230 (4), spilled 🟢 0/4136 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**26.830 ms/f**, peak render **32.230 ms**, spilled
**0/4136** live frames (**0.00%**).
Window summaries (17–4128): `fx_shader_draw` mean 22.884 ms/f;
worst window 27.799 ms/f (3953–3968).

Startup setup render: **50.393 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
30.270 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Four authored presets, 600-frame dwell and 480-frame eased parameter departure.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 3953–3968, richest in this regime)

```text
frame                   62.440 ms 37.464 M 100.0%
  pov_preserve_half        0.141 ms  0.084 M   0.2% x1.0 140.82 us
  fx_shader_draw          27.799 ms 16.680 M  44.5% x1.0 27799.46 us
  fx_prepare_frame         1.012 ms  0.607 M   1.6% x1.0 1012.49 us
  fx_advance               2.093 ms  1.256 M   3.4% x1.0 2093.41 us
  fx_timeline_step         0.154 ms  0.092 M   0.2% x1.0 153.92 us
  canvas_clear             0.088 ms  0.053 M   0.1% x1.0 88.27 us
  canvas_buffer_wait      31.115 ms 18.669 M  49.8% x1.0 31115.25 us
```

Wall min/avg/max = 61.207/62.440/63.746 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 593–608, richest in this regime)

```text
frame                   62.533 ms 37.520 M 100.0%
  pov_preserve_half        0.142 ms  0.085 M   0.2% x1.0 142.47 us
  fx_shader_draw          26.439 ms 15.864 M  42.3% x1.0 26439.40 us
  fx_prepare_frame         0.600 ms  0.360 M   1.0% x1.0 599.51 us
  fx_advance               1.928 ms  1.157 M   3.1% x1.0 1928.27 us
  fx_timeline_step         0.116 ms  0.070 M   0.2% x1.0 115.87 us
  canvas_clear             0.087 ms  0.052 M   0.1% x1.0 87.19 us
  canvas_buffer_wait      33.187 ms 19.912 M  53.1% x1.0 33186.69 us
```

Wall min/avg/max = 60.491/62.533/64.327 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # (one-based marker order) | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 1 | coupled-grid | 55 / 54 | not recorded | 27.799 | 31.325 | 16.02 |
| 2 | direct-grid | 67 / 66 | not recorded | 27.177 | 31.023 | 16.02 |
| 3 | double-map | 68 / 67 | not recorded | 23.848 | 26.893 | 16.02 |
| 4 | stretched-grid | 67 / 66 | not recorded | 23.799 | 28.317 | 16.02 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner (one-based marker order): preset 1 (coupled-grid). Per-preset peaks span 28.432–32.230 ms.

| # (one-based marker order) | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 1 | coupled-grid | 32.230 | 0/899 (0.00%) |
| 2 | direct-grid | 31.439 | 0/1079 (0.00%) |
| 4 | stretched-grid | 28.874 | 0/1079 (0.00%) |
| 3 | double-map | 28.432 | 0/1079 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–4128; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.09/f 0.62/1.70/24.15 us 3.13% CPU
isr_pack        144.00/f 6.24/7.04/22.30 us 1.62% CPU
isr_dma_submit  144.00/f 0.61/0.93/11.77 us 0.22% CPU
```

- DMA submit averages 0.93 us/event; packing averages 7.04 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.97%, leaving approximately 59.392 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 27.799 ms/f, 44.5% of the richest runtime window.
2. `fx_advance`: 2.093 ms/f, 3.4% of the richest runtime window.
3. `fx_prepare_frame`: 1.012 ms/f, 1.6% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 26.336 ms/f,
peak 31.997 ms, spills 0/4136.
Candidate `0156d0d7490355ea99cab406f4704a128a55ccbf`: mean render 26.830 ms/f
(+1.878% in these captures), peak 32.230 ms,
spills 0/4136. Baseline and candidate include 4136 and
4136 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 1536.384 us/f versus baseline 1021.429 us/f. Shading averages 22.884 versus 22.932 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; the global-O3 twin optimizes all compiled code.
- Only the epoch is stretched to 300 s; dwell and transition settings remain authored.
- Captured source is the base SHA above plus the exact retained source patch; the source diff hash is 783e073d58c6daf2ebaac200edc4f0fbee222ef9a9e5d462a84ec36d4d08b7d3. This capture does not represent an unmodified committed tree. Build flags, source status and hashes are retained as evidence; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=KaleidoscopeSmooth`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 11:05**.
