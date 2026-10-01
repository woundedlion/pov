# LatticeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-01, **-O3**)

Shipping twin: [LatticeMelt selective -O3](../shipping/profile_latticemelt_architecture_teensy_2026-10-01.md).

Point-in-time architecture snapshot. Raw capture: [checkpoint2-latticemelt-profile_o3.txt](../evidence/architecture_2026-10-01/checkpoint2-latticemelt-profile_o3.txt); [capture provenance](../evidence/architecture_2026-10-01/checkpoint2-latticemelt-profile_o3.provenance). Updates the current ranking alongside the preserved earlier capture [profile_latticemelt_teensy_2026-08-26.md](profile_latticemelt_teensy_2026-08-26.md). Earlier reports used different source tips and are not the architecture baseline.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, live flywheel and DMA ISRs |
| Image | `profile_o3`: `-O3`, `-ffast-math -fno-finite-math-only`; global optimization |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master |
| Effect | LatticeMelt 288×144, single-entry playlist, source `d25dd85dee17c28d2728f813e1ae00ae64690e1b` |
| Method | `HS_PROFILE`, 16-frame windows, 110 s capture; `HS_PROFILE_EPOCH_REVS=1200` (150 s epoch); runtime 2–1737, setup frame 1 excluded; counter summaries 17–1728, excluding the entire startup-containing window |
| Reproduce | `tools/profile_one.sh LatticeMelt profile_o3 110 16 -D HS_PROFILE_EPOCH_REVS=1200` |

Image size: `FLASH: code:97352, data:156460, headers:8332` /
`RAM1: variables:315040, code:39432, padding:26104, free:143712` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window 1–16 root cycles / 600 MHz
match the wall sum within **3.80 ppm**. The untouched baseline and candidate captures
both passed `tools/parse_profile.py <capture> validate` before filtering; their validation
outputs and build/envdump records are in the [evidence directory](../evidence/architecture_2026-10-01/README.md).

## Frame cadence

**Runtime aggregate**, exact per-frame telemetry, setup excluded: mean render
**38.156 ms/f**, peak render **38.705 ms**, spilled
**0/1736** live frames (**0.00%**).
Window summaries (17–1728): `fx_shader_draw` mean 34.667 ms/f;
worst window 34.832 ms/f (385–400).

Startup setup render: **72.880 ms** for frame 1, drawn before publication.
Every subsequent raw per-frame record is included, including warm-up, departures, fades and spills. The final 9 records beyond the last dumped counter window remain in exact render statistics and preset attribution.
One 62.5 ms display window covers a quadrant of approximately 10,368 pixel positions.
The measured worst live render fits that window by
23.795 ms. `canvas_buffer_wait` is display-sync idle to the next flip.

## Phase-by-phase readout

Phase schedule: Two authored presets, 600-frame dwell and 480-frame eased parameter departure.

Trees show per-frame averages. `M` means million cycles, `x` means calls/frame,
and the final microsecond figure on each leaf is its average cost/call.
Percentages divide by the root frame counter, rather than the parent counter.

### Single-owner window (frames 385–400, richest in this regime)

```text
frame                   62.438 ms 37.463 M 100.0%
  pov_preserve_half        0.138 ms  0.083 M   0.2% x1.0 138.06 us
  fx_shader_draw          34.832 ms 20.899 M  55.8% x1.0 34832.27 us
  fx_prepare_frame         1.112 ms  0.667 M   1.8% x1.0 1111.68 us
  fx_advance               2.234 ms  1.340 M   3.6% x1.0 2233.95 us
  fx_timeline_step         0.074 ms  0.044 M   0.1% x1.0 73.66 us
  canvas_clear             0.086 ms  0.051 M   0.1% x1.0 85.64 us
  canvas_buffer_wait      23.937 ms 14.362 M  38.3% x1.0 23936.97 us
```

Wall min/avg/max = 62.156/62.438/62.702 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Ownership-changing window (frames 593–608, richest in this regime)

```text
frame                   62.452 ms 37.471 M 100.0%
  pov_preserve_half        0.141 ms  0.085 M   0.2% x1.0 141.36 us
  fx_shader_draw          34.719 ms 20.831 M  55.6% x1.0 34718.95 us
  fx_prepare_frame         0.999 ms  0.599 M   1.6% x1.0 999.13 us
  fx_advance               2.143 ms  1.286 M   3.4% x1.0 2142.51 us
  fx_timeline_step         0.123 ms  0.074 M   0.2% x1.0 122.99 us
  canvas_clear             0.086 ms  0.051 M   0.1% x1.0 85.58 us
  canvas_buffer_wait      24.217 ms 14.530 M  38.8% x1.0 24216.96 us
```

Wall min/avg/max = 58.968/62.451/65.723 ms. The scope tree includes interrupt preemption; display-sync idle accounts for the remainder of the frame. Ownership markers distinguish these windows, but a single-owner window can still contain a parameter morph. These trees are measured regime examples rather than isolated pure-hold estimates.

### Per-preset table

All authored presets were visited and the marker sequence wrapped to the first index. Rows use the richest complete subsequent window with one owner and the modal draw-call count. They are ranked by scope cost; parameter morphs can remain within a single-owner window.

| # | Preset | Windows / clean | blended px/f | Scope ms | Render ms | fps |
|---|---|--:|--:|--:|--:|--:|
| 1 | open-curl | 40 / 39 | not recorded | 34.832 | 38.501 | 16.02 |
| 2 | dense-curl | 67 / 66 | not recorded | 34.745 | 38.524 | 16.02 |

Cadence buckets in the READMEs use all owner-attributed live frames, including the departure following that preset, rather than only these selected windows.

Worst live owner: preset 1 (open-curl). Per-preset peaks span 38.682–38.705 ms.

| # | Live owner | Peak render ms | Spilled / live frames |
|---|---|--:|--:|
| 1 | open-curl | 38.705 | 0/657 (0.00%) |
| 2 | dense-curl | 38.682 | 0/1079 (0.00%) |

### Per-pixel figures

This capture has no `filter_blend` counter and records no exact blended-pixel count. A quadrant has 10,368 nominal shader sample positions. Scope time divided by that count describes cost per nominal position, not measured cycles per successful blend.

## Column-ISR / DMA marshaling cost

Complete counter windows 17–1728; rates divide total events by their frame count, CPU shares divide ISR time by total measured window duration.

```text
isr_wake        1152.08/f 0.39/1.59/27.54 us 2.93% CPU
isr_pack        144.00/f 5.98/6.86/19.09 us 1.58% CPU
isr_dma_submit  144.00/f 0.59/0.94/11.84 us 0.22% CPU
```

- DMA submit averages 0.94 us/event; packing averages 6.86 us/event. Wire transmission continues asynchronously.
- The 24 MHz requested SPI clock, 600-byte composite frame and configured LPSPI divider/delays model **230 us** per transfer; this is a framing-model duration, not an independently measured wire trace.
- Total measured ISR CPU share is 4.73%, leaving approximately 59.544 ms of CPU-only work per 62.5 ms window. Render scopes already include interrupts, so this budget must not be subtracted from measured render a second time. No measured regime requires a speedup to hold 16 fps.

## Summary ranking

1. `fx_shader_draw`: 34.832 ms/f, 55.8% of the richest runtime window.
2. `fx_advance`: 2.234 ms/f, 3.6% of the richest runtime window.
3. `fx_prepare_frame`: 1.112 ms/f, 1.8% of the richest runtime window.

Architecture baseline `b6d2c0100e10d57b18f9bd749d5ab28b34574dd8`: mean render 38.163 ms/f,
peak 38.672 ms, spills 0/1736.
Candidate `d25dd85dee17c28d2728f813e1ae00ae64690e1b`: mean render 38.156 ms/f
(-0.017% in these captures), peak 38.705 ms,
spills 0/1736. Baseline and candidate include 1736 and
1736 live frames respectively. These finite samples can differ in camera/noise
phase and transition timing; a lower observed average is not a general speedup claim.
Host or WASM timings are not substituted for device measurements.

Preparation over the complete post-startup windows averages 872.996 us/f versus baseline 872.491 us/f. Shading averages 34.667 versus 34.680 ms/f. These separate counters distinguish preparation cost from shader cost.

## Caveats

- CYCCNT free-runs, so every scope absorbs ISR preemption and nested scopes overlap their parents.
- No `filter_blend` subtree or exact blend count was recorded in these captures.
- Shipping uses the composed shade and generated-palette HS_O3 functions; sphere curl-noise displacement; the global-O3 twin optimizes all compiled code.
- Only the epoch is stretched to 150 s; dwell and transition settings remain authored.
- Captured source is the exact committed SHA above. Build flags, source status and hashes are retained as evidence; later documentation edits do not change that source identity.
- Counter summaries omit the startup-containing window; exact cadence excludes only actual setup frame 1 and retains all following live frames.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=LatticeMelt`, `HS_PROFILE_WINDOW=16`;
the reproduce command builds, flashes and captures using the real segmented driver.
Raw log mtime (America/Los_Angeles): **2026-10-01 09:49**.

Global -O3 versus selective -O3: mean render 38.156 versus 36.120 ms/f (ratio 0.947× for these samples). The O3 image adds +25,752 B FLASH code and +19,472 B ITCM. Live frame counts differ; compare spill fractions and peaks alongside these means.
