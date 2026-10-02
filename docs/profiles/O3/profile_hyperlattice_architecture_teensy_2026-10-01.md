# HyperLattice on-device profile — Teensy 4.0, segmented mode (2026-10-01, **-O3**)

Shipping twin: [HyperLattice selective -O3](../shipping/profile_hyperlattice_architecture_teensy_2026-10-01.md).

Point-in-time architecture snapshot. Raw capture: checkpoint1-hyperlattice-profile_o3.txt (supporting artifact removed); capture provenance (supporting artifact removed). Updates the current ranking alongside the preserved earlier capture [profile_hyperlattice_teensy_2026-09-27.md](profile_hyperlattice_teensy_2026-09-27.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile_o3`: `-O3`, `-ffast-math -fno-finite-math-only`; global optimization |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | HyperLattice 288×144, single-entry playlist, source `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b` |
| Method | `HS_PROFILE`, 16-frame windows, 170 s capture; `HS_PROFILE_EPOCH_REVS=1600` (200 s epoch); runtime 2–2697, setup frame 1 excluded; counter summaries 17–2688, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh HyperLattice profile_o3 170 16 -D HS_PROFILE_EPOCH_REVS=1600` |

Image size: `FLASH: code:83872, data:152828, headers:9060` /
`RAM1: variables:314976, code:20504, padding:12264, free:176544` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 1185–1200 root cycles / 600 MHz
match the wall sum within **2.03 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are no longer retained in the repository.

## Frame cadence

README cells: peak 🟢 36.975 (3), spilled 🟢 0/2696 (0.00%).

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**25.296 ms/f**, peak render **36.975 ms**, spilled
**0/2696** live frames (**0.00%**).
Window summaries (17–2688): `hl_shader_draw` mean 23.282 ms/f;
worst window 28.318 ms/f (1185–1200).

Startup setup render: **42.528 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
25.525 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Three normal-roster presets, 320-frame dwell, with 240-frame eased parameter morphs or fades selected by the authored departure.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 1185–1200, richest in this regime)

```text
frame                   62.851 ms 37.710 M 100.0%
  pov_preserve_half        0.143 ms  0.086 M   0.2% x1.0 142.54 us
  hl_shader_draw          28.318 ms 16.991 M  45.1% x1.0 28317.91 us
  hl_timeline_step         0.004 ms  0.002 M   0.0% x1.0 3.50 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 85.35 us
  canvas_buffer_wait      32.552 ms 19.531 M  51.8% x1.0 32551.65 us
```

Wall min/avg/max = 57.405/62.850/68.400 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 2673–2688, richest in this regime)

```text
frame                   62.514 ms 37.508 M 100.0%
  pov_preserve_half        0.145 ms  0.087 M   0.2% x1.0 144.53 us
  hl_shader_draw          25.453 ms 15.272 M  40.7% x1.0 25453.33 us
  hl_timeline_step         0.011 ms  0.007 M   0.0% x1.0 11.24 us
  canvas_clear             0.085 ms  0.051 M   0.1% x1.0 85.40 us
  canvas_buffer_wait      35.033 ms 21.020 M  56.0% x1.0 35033.16 us
```

Wall min/avg/max = 61.996/62.513/63.353 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # (one-based marker order) | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 3 | hypercube-flight | 35 / 34 | not recorded | 28.318 | 30.299 | 15.91 |
| 2 | cubic-wide-flight | 86 / 84 | not recorded | 24.920 | 26.927 | 16.09 |
| 1 | cubic-flight | 46 / 44 | not recorded | 23.265 | 25.277 | 15.99 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner (one-based marker order): preset 3 (hypercube-flight). Per-preset peaks span 26.078–36.975 ms.

| # (one-based marker order) | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 3 | hypercube-flight | 36.975 | 0/581 (0.00%) |
| 2 | cubic-wide-flight | 30.313 | 0/1358 (0.00%) |
| 1 | cubic-flight | 26.078 | 0/757 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–2688; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.10/f 0.41/1.52/12.48 us 2.81% CPU
isr_pack        144.00/f 5.98/6.55/9.68 us 1.51% CPU
isr_dma_submit  144.00/f 0.60/0.94/4.21 us 0.22% CPU
```

- DMA submit averages 0.94 us/event; packing averages 6.55 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.54%, leaving approximately 59.665 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `hl_shader_draw`: 28.318 ms/f, 45.1% of the richest runtime window.
2. `pov_preserve_half`: 0.143 ms/f, 0.2% of the richest runtime window.
3. `canvas_clear`: 0.085 ms/f, 0.1% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 25.316 ms/f,
peak 37.004 ms, spills 0/2696.
Candidate `ac4870adb95eecb7130f86dbfd36a86d0b1fb96b`: mean render 25.296 ms/f
(-0.077% in these captures), peak 36.975 ms,
spills 0/2696. Baseline and candidate include 2696 and
2696 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses `HS_HOT_FLASH_MEMBER` (`HS_O3_FN` plus `hot`) for the cached shader, with no ITCM placement; the global-O3 twin optimizes all compiled code.
- Only the epoch is stretched to 200 s; dwell and transition settings remain authored.
- Captured source identity is recorded above; see the archive source-reachability note. Build flags, source status and hashes are not retained in the repository; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=HyperLattice`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:26**.

Global -O3 versus selective -O3: mean render 25.296 versus 25.652 ms/f (ratio 1.014× for these samples). The O3 image adds +9,616 B FLASH code and +4,912 B ITCM. Live frame counts differ; compare spill fractions and peaks alongside these means.
