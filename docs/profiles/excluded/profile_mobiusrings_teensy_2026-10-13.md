# MobiusRings on-device profile — Teensy 4.0, segmented mode (2026-10-13, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MobiusRings`).
Raw captures and ELF/source provenance: `C:/work/temp/review1013-ab-*.log`.
This is the first archived MobiusRings segmented profile. MobiusRings remains
excluded from the Phantasm playlist as a low-resolution-only effect.

## Setup

| | |
|---|---|
| Hardware | Two Teensy 4.0 boards at 600 MHz; COM3 serial 6362410, COM4 serial 6362390 |
| Image | `profile`, selective O3; filter/pipeline hot regions active |
| Driver | `POVSegmented<288, 4, 480>`, segment 0 master; flywheel and DMA ISRs live |
| Effect | MobiusRings 288×144; baseline `1930ffccf89a`, candidate `c2c8e0fd6812` |
| Method | 70 s captures, 32-frame windows, epoch 1200 revolutions; fixed effect seed |
| Reproduce | `HS_PROFILE_TREE=<tree> HS_TEENSY_PORT=<COMn> bash tools/profile_one.sh MobiusRings profile 70 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

Candidate image: FLASH code 80284, data 153716, headers 8684 bytes;
RAM1 variables 315072, code 39656, padding 25880, free 143680 bytes;
RAM2 variables 520064, free 4224 bytes. Relative to the profile baseline:
RAM1 code +608, RAM1 variables 0, FLASH data +96 bytes. The shipping Phantasm
image has zero change in RAM1 code, RAM1 variables and FLASH data because this
effect is excluded from that playlist.

Exactness cross-check: frames 129–160, root cycles / 600 MHz
and measured wall sum differ by 3.8 ppm. All four logs pass
`tools/parse_profile.py ... validate`.

## Frame cadence

A display window is 62.5 ms. The effect renders a quadrant of 10,368 pixels.
Both versions hold the 16 fps display cadence without spilled frames.
`canvas_buffer_wait` is display synchronization idle, removed from render cost.
Candidate pass peak render is 12.30 ms; spilled 0/1088 frames.
README cells: peak 🟢 12.30, spilled 🟢 0/1088 (0.00%).

### Finding 6 matched A/B

COM3 ran baseline then candidate; COM4 ran candidate then baseline. Each comparison
uses identical frame indices starting at frame 65; startup frames are discarded.
The source and ELF hashes are retained beside every capture. Each board runs both
versions, controlling board variation. The cached cosine test retains vector-length
correction, avoiding the visual error caused by assuming perfectly unit sampler vectors.

| Board | Matched frames | Baseline mean ms | Candidate mean ms | Saving | Peak old / new ms | Spilled old / new |
|---|--:|--:|--:|--:|--:|--:|
| COM4 | 1034 | 8.2036 | 8.0830 | 1.470% | 12.479 / 12.299 | 0 / 0 |
| COM3 | 1034 | 8.2037 | 8.0828 | 1.474% | 12.461 / 12.297 | 0 / 0 |

Both same-board comparisons show a positive mean render saving. This is a modest
CPU saving; display frame rate stays 16 fps. Native parity checks cover radius
updates, near-unit vectors, disabled holes and radii beyond a half turn.

## Phase-by-phase readout

The animation is continuous: warp repeats every 160 frames, ring/line density
every 320, orientation every 400, and palette wiping every 120. The captures
cover multiple repeats of every individual period, but not their combined cycle.

### Continuous animation, worst render window

Frames 129–160; per-frame averages:

```
frame                   62.441 ms  37.465 Mcyc 100.0%
  pov_preserve_half        0.144 ms   0.086 Mcyc   0.2%
  mr_draw_grid             8.229 ms   4.937 Mcyc  13.2%
    mr_lines_draw            7.579 ms   4.547 Mcyc  12.1%
    mr_rings_draw            0.649 ms   0.390 Mcyc   1.0%
  mr_wipe_rebake           0.887 ms   0.532 Mcyc   1.4%
  mr_timeline_step         0.066 ms   0.040 Mcyc   0.1%
  canvas_clear             0.084 ms   0.051 Mcyc   0.1%
  canvas_buffer_wait      53.023 ms  31.814 Mcyc  84.9%
```

Wall min/avg/max: 59.075 / 62.441 /
65.752 ms. Line rendering dominates the draw grid cost.
Most elapsed frame time is synchronization idle. No separate preset table applies.

### Per-pixel figures

Blended-pixel counts were not enabled in these timing images. The nominal quadrant
contains 10,368 pixels; actual blended-pixel and cycles-per-blend figures are unavailable.

## Column-ISR / DMA marshaling cost

Representative ISR readings from the same frame window; the ISR interval is
2.000414 s. Rates divide call counts by that interval; per-call times use the
logged nanosecond readings.

```
scope               calls/s      min/avg/max us          CPU
isr_wake            18433.18     0.626/1.633/17.310       3.01%
  isr_pack           2304.02     6.231/6.605/9.096        1.52%
  isr_dma_submit     2304.02     0.775/0.935/1.011        0.21%
```

`isr_wake` includes pack and submit; its CPU share is the measured flywheel ISR
share. Pack costs more than submit. LED wire transfer is asynchronous; its exact
duration was not measured. DMA completion and other interrupts are not measured.
The roughly 3% inclusive wake share leaves about 60.6 ms of the 62.5 ms window
for rendering; neither version requires a speedup to meet that budget.

## Summary ranking

1. `mr_draw_grid`, especially `mr_lines_draw`, dominates active rendering.
2. `canvas_buffer_wait` dominates wall time and contributes no render work.
3. Half-buffer preservation and clearing have smaller costs.

No matched WASM/native timing comparison was captured. The two-board device A/B
supports landing the cached cap classification without a performance regression.

## Caveats

ISR time is included in the cycle scopes. Filter-blend parenting can attribute
work beneath the active caller; per-pixel scope instrumentation adds overhead.
This measures selective-O3 profile flags rather than global O3. The epoch stretch
changes effect lifetime, not per-frame work. Only temporary host upload/capture
scripts were uncommitted; compiled effect/filter code is pinned to each recorded
commit. Fresh loader processes and explicit ports prevent stale GUI selection.
The comparison covers about 1,000 matched frames per board, not every possible
parameter combination or the combined animation period.

## Harness

`targets/Profile/Profile.ino` supplies fixed seeds, per-frame render telemetry,
32-frame cycle windows and ISR readings. `just profile MobiusRings 70` reproduces
the default run; the exact A/B commands additionally set epoch 1200 and pin the
source worktree and board as shown above.
