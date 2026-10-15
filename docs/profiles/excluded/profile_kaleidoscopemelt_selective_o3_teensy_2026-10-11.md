# KaleidoscopeMelt on-device profile — Teensy 4.0, segmented mode (2026-10-11, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile KaleidoscopeMelt`).
Raw capture: `build/prof/kaleidoscopemelt_ship.log`, captured on COM3 at
2026-10-11 11:38 America/Los_Angeles. First report for this effect.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, POV segmented mode, flywheel + DMA ISRs live, COM3 |
| Image | `profile` env: `-Os` base + newlib-nano, shipping selective `HS_O3` shader traversal, composed shading, curl surface displacement, lattice distance and palette colorization |
| Driver | `POVSegmented<288, 4, 480>`, board = segment 0 master |
| Effect | KaleidoscopeMelt 288×144, single-entry playlist, base tip `2ae741709553ce016e00e0751f1b1013555ce3de` plus the uncommitted roster replacement |
| Method | `HS_PROFILE` cycle scopes, window = 32 frames, 70 s capture; single preset, no dwell-compression or epoch override |
| Reproduce | `HS_PROFILE_TREE=/c/work/temp/Holosphere-melt-20261011 HS_TEENSY_PORT=COM3 bash tools/profile_one.sh KaleidoscopeMelt profile 70 32` |

Image size (`profile` env): `FLASH: code:70180, data:155252, headers:9060` /
`RAM1: variables:315040, code:19704, padding:13064, free:176480` /
`RAM2: variables:520064, free:4224`.

Exactness cross-check: window frames 545–576 root counter
1,199,175,292 cycles ÷ 600 MHz = 1,998,625.487 us, versus the measured wall
sum 1,998,632 us, within **3.3 ppm**
(`tools/parse_profile.py ... validate`, VALID). All 34 complete windows
name KaleidoscopeMelt at 288×144; frame numbers are monotonic, per-frame rows
are complete and no epoch reset occurs. The wrapper verified the configuration,
fresh first frame 1 and matching saved ELF provenance.

## Frame cadence

**Pass aggregate** (`parse_profile.py ... windows` footer): `fx_shader_draw`
averages **37.80 ms/f**, with worst window **39.94 ms/f** at frames 545–576.
The full parser aggregate is peak render **76.53 ms**, spilled **1/1088**
frames (0.09%): setup frame 1 accounts for the only spill.

For consistency with the archive's display-cadence tables, setup frame 1 is
excluded from the README figures. Frames 2–1088 peak at **46.26 ms**, with
**0/1087** spilled (0.00%). README cells: peak 🟢 46.26, spilled 🟢 0/1087 (0.00%).

A display window is 62.5 ms at 480 RPM. The effect shades one 144×72 quadrant,
10,368 pixels, and holds **16 fps** after setup, with at least 16.24 ms of
render margin in this capture. Setup renders for 76.534 ms once; it does not
represent a recurring preset transition. The `canvas_buffer_wait` scope is
the round-up idle to the next display flip, by design.

## Phase-by-phase readout

Phase schedule: one initial setup frame followed by a single continuously
animated preset. Camera and curl-field evolution vary shader cost; there is
no discrete preset cycle. The blocks use per-frame averages across 32 frames,
with nesting matching the capture. Each non-root scope runs once per frame,
except the startup window's 31 preservation calls; its per-call cost is its
listed total divided by that call rate.

### Startup window (frames 1–32)

```
frame                    62.25 ms  37.35 Mcyc 100.0%
  pov_preserve_half      135.6 us   81.4 kcyc   0.2% x0.97 140 us/call
  fx_shader_draw         37.54 ms  22.52 Mcyc  60.3% x1 37539 us/call
  fx_prepare_frame        3.92 ms   2.35 Mcyc   6.3% x1 3925 us/call
  fx_advance              1.95 ms   1.17 Mcyc   3.1% x1 1954 us/call
  fx_timeline_step        63.6 us   38.2 kcyc   0.1% x1 64 us/call
  canvas_clear            90.0 us   54.0 kcyc   0.1% x1 90 us/call
  canvas_buffer_wait     18.51 ms  11.11 Mcyc  29.7% x1 18509 us/call
```

Wall min/avg/max = **43.004/62.252/76.534 ms**. Frame 1 initializes the
first-frame rendering path and causes the only spill. The following frames
already fit one display window; averaging this window distributes that
one-time cost over 32 frames.

### Highest shader-cost steady window (frames 545–576)

```
frame                    62.46 ms  37.47 Mcyc 100.0%
  pov_preserve_half      138.8 us   83.3 kcyc   0.2% x1 139 us/call
  fx_shader_draw         39.94 ms  23.97 Mcyc  64.0% x1 39944 us/call
  fx_prepare_frame        3.57 ms   2.14 Mcyc   5.7% x1 3567 us/call
  fx_advance              1.99 ms   1.19 Mcyc   3.2% x1 1991 us/call
  fx_timeline_step        72.7 us   43.6 kcyc   0.1% x1 73 us/call
  canvas_clear            85.9 us   51.5 kcyc   0.1% x1 86 us/call
  canvas_buffer_wait     16.63 ms   9.98 Mcyc  26.6% x1 16628 us/call
```

Wall min/avg/max = **61.653/62.457/63.185 ms**. Render averages 45.829 ms
and peaks at 46.262 ms, leaving every frame within the 62.5 ms render budget.
The shader grows to 39.944 ms/frame versus 37.539 ms in the startup window,
while buffer wait falls to 16.628 ms/frame. Small wall variations around
62.5 ms reflect synchronization and readout timing, without a render spill.

### Per-pixel figures

The shader directly assigns all **10,368 pixels/frame**, one sample per
pixel, with 1× quadrant coverage. At the highest-cost steady window,
`fx_shader_draw` costs **2,311.58 cycles/pixel** (3.853 us/pixel), including
traversal and color assignment. This path does not call `filter_blend`, so
there is no measured blended-pixel count or cycles/blend figure.

## Column-ISR / DMA marshaling cost

The following readings are from frames 545–576. Pack and submit execute inside
the inclusive wake scope; rates are calls per rendered frame.

```
isr_wake       1152.56/f  min/avg/max 0.595/1.721/18.230 us  cpu 3.17%
  isr_pack      144.06/f  min/avg/max 6.436/7.268/10.306 us  cpu 1.67%
  isr_dma_submit 144.06/f  min/avg/max 0.590/0.940/1.045 us   cpu 0.21%
```

- Submit consumes 4,334 us per ISR window versus pack's 33,511 us: roughly
  13% of the pack CPU cost. DMA transfers asynchronously after submission.
- The segment sends a 600-byte image-plus-black composite at a requested
  24 MHz SPI clock. Driver byte framing gives a **230 us** transfer per
  column, below the approximately 434 us column interval.
- Inclusive `isr_wake` CPU share is **3.17%**, leaving approximately
  **60.52 ms** per 62.5 ms window for foreground work. Steady render's
  46.26 ms peak needs no speedup to hold 16 fps; setup's 76.53 ms would need
  about **1.26×** reduction to fit this remaining budget. Measured render
  already includes ISR interference, so this is a conservative comparison.
  DMA-completion and other interrupts are not measured.

## Summary ranking

1. `fx_shader_draw` — **63.96%**, **39.94 ms/frame** in the highest-cost
   steady window: dodecahedral folding, curl displacement, equirectangular
   projection, lattice sampling and palette colorization dominate.
2. `fx_prepare_frame` — **5.71%**, **3.57 ms/frame**: prepared frame state
   and lookup work before traversal.
3. `fx_advance` — **3.19%**, **1.99 ms/frame**: clocks, camera, palette
   and runtime advancement.

`canvas_buffer_wait` is 26.62%, 16.63 ms/frame, and represents cadence idle.
No measured WASM/native timing ledger entry exists for the newly added effect;
this capture supplies the device numbers without inferring host performance.

## Caveats

- All scopes absorb ISR time because CYCCNT free-runs under the live driver.
- `filter_blend` parenting can hide a subtree under an inactive first parent;
  this direct-assignment shader path emits no `filter_blend` counter.
- No per-pixel profiling scope was added. Its instrumentation overhead would
  distort the pixel cost; the quoted figure divides one frame-level scope.
- The shipping selective-O3 path includes the typed shader traversal,
  composed shade entry, curl surface displacement, lattice distance and
  generated-palette colorization. The lens/displacement/projection stage
  group retains the effect's out-of-line flash placement.
- No dwell-compression knobs were used. The single preset runs continuously.
- Captured from the isolated worktree's uncommitted KaleidoscopeMelt roster
  replacement. The saved artifact directory retains ELF files, build logs,
  build flags, source status and the source diff. Compiler and ARM ABI match
  the independently built, budget-passing shipping Phantasm image.
- README cadence figures exclude setup frame 1 consistently with the archive;
  the full parser aggregate including setup is reported above. The final nine
  individual frames after the last complete window are not in that aggregate.

## Harness

`targets/Profile/Profile.ino` + `HS_PROFILE_TARGET=KaleidoscopeMelt`,
`HS_PROFILE_WINDOW=32`; `just profile KaleidoscopeMelt` builds, flashes and
captures through the shared-device lock. Set `HS_PROFILE_TREE` explicitly when
profiling a worktree. The report's 70-second capture stays within one epoch.
