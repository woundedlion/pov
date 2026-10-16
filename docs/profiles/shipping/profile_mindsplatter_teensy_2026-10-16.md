# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-16, **selective -O3**)

Point-in-time snapshot. Replaces profile_mindsplatter_teensy_2026-10-06.md. Captured 2026-10-16 08:55 PDT on COM3.
Raw captures: `build/prof/mindsplatter_hd107_before_20261016.log` and
`build/prof/mindsplatter_hd107_after_20261016.log`; adjacent provenance files
record source, compiler and ELF hashes. Their artifact directories retain both
profile and shipping ELFs, maps, build logs and clean source status.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0, 600 MHz, same COM3 board for both passes |
| Image | Shipping selective-O3 `profile` environment, `-Os` base, live DMA/flywheel ISRs |
| Driver | `POVSegmented<288,4,480>`, segment 0 master, 72 LEDs |
| Effect | MindSplatter, 288×144, deterministic effect seed, full preset cycle |
| Method | Two 110-second captures, 16-frame windows; no dwell compression |
| Source | Before `2dc0422e8658`; HD packing `c883a255a0b9`; clean isolated worktrees |
| Reproduce | Set `HS_PROFILE_TREE` to each source checkout and `HS_TEENSY_PORT=COM3`; run `bash tools/profile_one.sh MindSplatter profile 110 16` with distinct `HS_PROFILE_OUT` paths |

HD profile image: FLASH code 70,044 B, data 549,708 B, headers 8,980 B;
RAM1 code 36,600 B, variables 313,152 B, padding 28,936 B, free 145,600 B;
RAM2 variables 520,064 B, free 4,224 B.
Versus the baseline profile image: RAM1 code +304 B, RAM1 variables +0 B,
FLASH data +0 B. The full Phantasm image's earlier +432 B code delta is a
different image and does not substitute for these single-effect measurements.

Exactness cross-check, frames 993–1008: HD root
596542709 cycles / 600 = 994237.848 µs,
wall sum 994240 µs, difference **2.16 ppm**. HD capture passes
all parser checks. The baseline passes identity, completeness, monotonic frames
and full-cycle checks, but its corresponding root/wall difference is
**9.77 ppm**, failing the 5 ppm exactness gate. Baseline comparison
figures are retained with this limitation; the ISR comparison uses independent
ISR snapshot intervals, not the root/wall subtraction.

## Frame cadence

MindSplatter was the highest peak-render entry in the existing shipping ranking.
Both captures visit all eight presets and wrap back to entry 1. A display window
is 62.5 ms; each frame renders a 10,368-pixel quadrant. Both passes sustain
16 fps with zero spilled frames. The `canvas_buffer_wait` scope is display-sync
idle, not render work.

| Metric | Before | HD packing | Difference |
|---|---:|---:|---:|
| Peak frame render | 55.382 ms | 55.714 ms | +0.332 ms |
| Mean render, frames 17–1728 | 24.345 ms | 24.480 ms | +0.135 ms |
| Spilled, frame 1 excluded | 0/1727 | 0/1727 | 0 |

HD peak margin is 6.786 ms. These are single-pass observations,
not a confidence interval or a claim that every render delta is packing alone.

## Phase-by-phase readout

Eight preset entries include their transitions. The heaviest clean window is
entry 7; entry 6 is the next largest. Entry 5 has the lowest draw cost.

### Peak regime, frames 993–1008

Per-frame averages; M denotes million cycles. Tree percentages use frame root.

```
frame                             62.14 ms 37.28 M 100.0%
  pov_preserve_half                0.13 ms  0.08 M   0.2%
  msp_draw_particles              48.46 ms 29.08 M  78.0%
    msp_particle_scan             48.45 ms 29.07 M  78.0%
      plot_ps_raster              35.62 ms 21.37 M  57.3% x576
      plot_ps_deferred             0.56 ms  0.34 M   0.9% x576
      plot_ps_gate                 6.53 ms  3.92 M  10.5% x1537
      plot_ps_tween                4.94 ms  2.96 M   8.0% x1537
  msp_particle_step                3.01 ms  1.80 M   4.8%
  msp_timeline_step                0.06 ms  0.04 M   0.1%
  canvas_clear                     0.08 ms  0.05 M   0.1%
  canvas_buffer_wait              10.39 ms  6.23 M  16.7%
```

Wall min/avg/max = 59.035/62.140/64.774 ms.
Particle drawing dominates; its worst-window mean rises from 48.130 to 48.459 ms.
The raster subtree contains most of the work. Leaf per-call costs are window
cycles divided by calls; `plot_ps_raster` averages
61.84 µs/call.

### Per-preset table

Clean-hold draw costs use the parser's modal call-count windows; transition
straddles are excluded. Peak/spill columns use per-frame owner buckets including
transitions. The cycle wraps to entry 1. Rows rank by clean draw cost.

| Entry | Clean draw ms/frame | Peak render ms | Spilled/frames | Clean/hold windows | fps |
|---|---:|---:|---:|---:|---:|
| 7 | 48.459 | 55.714 | 0/159 | 9/10 | 16 |
| 6 | 40.956 | 52.707 | 0/159 | 9/10 | 16 |
| 2 | 30.652 | 35.964 | 0/318 | 18/20 | 16 |
| 8 | 30.539 | 45.425 | 0/159 | 9/10 | 16 |
| 1 | 28.437 | 35.035 | 0/317 | 17/19 | 16 |
| 3 | 27.991 | 37.422 | 0/297 | 17/18 | 16 |
| 4 | 27.942 | 34.396 | 0/159 | 9/10 | 16 |
| 5 | 17.100 | 22.829 | 0/159 | 9/10 | 16 |

### Per-pixel figures

No `filter_blend` counter is present, so blended-pixel counts and cycles per blend
cannot be derived. Amortized packing overhead is 17.63 cycles per LED
(29.38 ns); this is a column measurement divided by 72, not a per-pixel timer.

## Column-ISR / DMA marshaling cost

Weighted over complete windows 17–1728; initial setup window omitted from both
passes. ISR counters use their own snapshot interval, including report output.
Both runs contain 246,585 column-pack calls in approximately 107.025 seconds.

| Metric | Before | HD packing | Difference |
|---|---:|---:|---:|
| Pack per 72-LED column | 5.714 µs | 7.830 µs | +2.115 µs (+37.0%) |
| Pack CPU share | 1.317% | 1.804% | +0.487 percentage points |
| Inclusive wake CPU share | 2.810% | 3.271% | +0.461 percentage points |
| DMA submit per column | 0.952 µs | 0.944 µs | -0.008 µs |

HD ISR tree: rate, per-call min/average/max, CPU share.

```
isr_wake          18433.2/s  0.51/ 1.77/24.22 us  3.27%
  isr_pack         2304.0/s  6.00/ 7.83/11.67 us  1.80%
  isr_dma_submit   2304.0/s  0.59/ 0.94/ 3.36 us  0.22%
```

`isr_wake` includes pack and submit; do not add its children to it.
DMA-completion and other interrupts are not measured by these counters.
The measured wake share leaves approximately
60.456 ms per 62.5 ms window for foreground
work before other interrupts. Packing adds 0.305 ms per display window
at 144 columns/window. No speedup is required to retain 16 fps in this capture.
Wire format and clock are unchanged: 300-byte image frame at 24 MHz = 100 µs;
the optional trailing black frame makes 600 bytes = 200 µs. Transfers are async.

## Summary ranking

1. Particle drawing remains the bottleneck: 48.459 ms/frame
   in the heaviest window, 78.0% of frame wall time.
2. HD packing costs +2.115 µs/column, about +0.49 percentage points of CPU.
3. Peak render rises by 0.332 ms in this paired observation; no cadence loss.

README cells: peak 🟢 55.71 (8), spilled 🟢 0/1727 (0.00%).
No fresh WASM/native performance comparison was run; host timing cannot measure
the physical LED ISR. Historical O3 results are different-source comparisons.

## Caveats

- Baseline exactness gate failure is reported above; HD validation passes.
- One full-cycle pass per version on one board; no statistical repeat study.
- Main-loop scopes include interrupt time because CYCCNT free-runs.
- No `filter_blend` counter is present; its usual first-parent artifact prevents
  inferring blend counts from a missing subtree. Per-pixel timers were not added.
- Shipping selective-O3 hot regions remain enabled, including particle scan and
  column packing. No global-O3 or dwell-compression override was used.
- Both builds had clean source state; the only functional delta is HD packing.
- The default strip and Candle gains remain enabled in both images.
- The isolated Python 3.12 environment used pinned PlatformIO requirements because
  system Python 3.14 could not load SSL. Compiler/ABI attestation passed.
- The board is left running the HD MindSplatter profiling image.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`,
`HS_PROFILE_WINDOW=16`. The supported `tools/profile_one.sh` wrapper held the
COM3 device lock through build, flash and capture. Default epoch 120 seconds
exceeds the 110-second capture; no epoch reset occurred.
Use `just profile MindSplatter 110` for a fresh default-window capture, or the
Setup command for this exact 16-frame protocol.
