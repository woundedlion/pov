# MindSplatter on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with `just profile MindSplatter`).
Raw capture: `build/prof/mindsplatter_ship.log`, copied with its original mtime from `C:/work/temp/flash-20261030-mindsplatter-after-com3.log`. This replaces `profile_mindsplatter_teensy_2026-10-29.md`.

## Setup

| Item | Configuration |
|---|---|
| Hardware | Teensy 4.0, COM3, 600 MHz; flywheel and DMA ISRs live |
| Image | Shipping selective-O3 `profile`; existing MindSplatter trail/point/lattice and Plot hot regions; no deep, PC, replay or count flags |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | MindSplatter 288×144, single-entry playlist, source `c87e96523a1b7cc97e8e7bb5529cc52e78ad1503` |
| Method | 90 s, window 16; normal preset dwells; one complete eight-preset cycle and wrap, within the default 120 s profiling epoch |
| Reproduce | `HS_PROFILE_TREE=C:/work/temp/Holosphere-flash-20261030 HS_TEENSY_PORT=COM3 HS_TEENSY_TOOLS=C:/work/temp/teensy-tools-loopback-20261030 bash tools/profile_one.sh MindSplatter profile 90 16` |

Profile image size: `FLASH: code:164216, data:546336, headers:8296   free for files:1312768 / RAM1: variables:313152, code:60568, padding:4968   free for local variables:145600 / RAM2: variables:520064  free for malloc/new:4224`.

Shipping image: FLASH code 1,039,544 B, data 840,984 B, headers 8,752 B, total 1,889,280 B. Against its parent, total Flash falls 3,072 B; RAM1 variables 312,736 B and code 183,496 B are unchanged, as are RAM2 variables 520,064 B and both RAM allocations.

Exactness cross-check: frames 993–1008 root cycles / 600 MHz = 996784.837 us against measured wall sum 996782 us, a 2.85 ppm delta. `parse_profile.py ... validate` passes markers, wrap, frame coverage and epoch checks. Captured 2026-10-30 09:29 America/Los_Angeles.

The change shares setup-only integer registration and option-label checks. Same-board COM3 baseline (`73e418e2c`, `flash-20261030-mindsplatter-before-full.log`) peaks at 57.761 ms; this capture peaks at 57.717 ms (-0.076%), with zero spilled frames in both 1,408-frame counter-window sets. The later parent also includes an independent pullback color change; MindSplatter uses its existing renderer. A supplementary pre-rebase COM4 candidate capture peaks at 57.848 ms with 0/1,408 spills.

## Frame cadence

Pass aggregate (`parse_profile.py ... windows`): peak frame render **57.717 ms**, spilled **0/1408 (0%)**. Setup frame 1 is 0.410 ms; the README separates setup from the remaining 1,407 frames.

At 480 RPM, a half revolution is 62.5 ms, giving 16 fps. The clipped quadrant is 10,368 pixels. Every captured preset and transition stays within one display window, with 4.783 ms at the peak. `canvas_buffer_wait` is the round-up idle until the next display flip.

README cells: peak 🟢 57.717 (8), spilled 🟢 0/1407 (0.00%).

## Phase-by-phase readout

Preset dwell is 112 frames, departure LERP 48 frames, with an extended initial dwell. Markers stamp adoption at departure start. The table below selects only windows fully beyond the 48-frame blend and before the next marker; cadence buckets also include transitions.

### Initial population (frames 1–16)

```
frame                     54.744 ms  32.846 Mcyc 100.0%
  pov_preserve_half       140.25 us   84.15 kcyc   0.3% x0.94 150 us/call
  msp_draw_particles      823.43 us  494.06 kcyc   1.5%
    msp_particle_scan       0.07 us    0.04 kcyc   0.0% x0.06 1 us/call
  msp_particle_step       126.45 us   75.87 kcyc   0.2% x1 126 us/call
  msp_timeline_step        82.85 us   49.71 kcyc   0.2% x1 83 us/call
  canvas_clear             96.97 us   58.18 kcyc   0.2% x1 97 us/call
  canvas_buffer_wait      53.463 ms  32.078 Mcyc  97.7% x1 53463 us/call
```

Wall min/avg/max = 0.412/54.743/62.733 ms. The particle population grows from startup. This window is dominated by display synchronization rather than drawing; registration finishes before these frame scopes.

### Held preset 7 (frames 1009–1024)

```
frame                     61.449 ms  36.869 Mcyc 100.0%
  pov_preserve_half       134.62 us   80.77 kcyc   0.2% x1 135 us/call
  msp_draw_particles      39.508 ms  23.705 Mcyc  64.3% x1 39508 us/call
  msp_particle_step        2.646 ms   1.588 Mcyc   4.3% x1 2646 us/call
  msp_timeline_step       114.99 us   69.00 kcyc   0.2% x1 115 us/call
  canvas_clear             92.40 us   55.44 kcyc   0.2% x1 92 us/call
  canvas_buffer_wait      18.933 ms  11.360 Mcyc  30.8% x1 18933 us/call
```

Wall min/avg/max = 49.653/61.449/72.559 ms. This is the costliest fully held draw window. Particle drawing dominates the foreground cost; all its frames retain 16 fps.

### Parameter transition (frames 961–976)

```
frame                     62.162 ms  37.297 Mcyc 100.0%
  pov_preserve_half       134.43 us   80.66 kcyc   0.2% x1 134 us/call
  msp_draw_particles      33.926 ms  20.356 Mcyc  54.6% x1 33926 us/call
  msp_particle_step        7.773 ms   4.664 Mcyc  12.5% x1 7773 us/call
  msp_timeline_step       135.66 us   81.39 kcyc   0.2% x1 136 us/call
  canvas_clear             91.59 us   54.95 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait      20.071 ms  12.043 Mcyc  32.3% x1 20071 us/call
```

Wall min/avg/max = 56.297/62.161/68.000 ms. This window lies wholly within a logged LERP interval. Its particle-update and draw costs remain inside the display cadence budget.

### Peak transition-to-hold anchor (frames 993–1008)

```
frame                     62.299 ms  37.379 Mcyc 100.0%
  pov_preserve_half       136.34 us   81.81 kcyc   0.2% x1 136 us/call
  msp_draw_particles      49.677 ms  29.806 Mcyc  79.7% x1 49677 us/call
  msp_particle_step        2.959 ms   1.775 Mcyc   4.7% x1 2959 us/call
  msp_timeline_step       134.09 us   80.46 kcyc   0.2% x1 134 us/call
  canvas_clear             92.32 us   55.39 kcyc   0.1% x1 92 us/call
  canvas_buffer_wait       9.275 ms   5.565 Mcyc  14.9% x1 9275 us/call
```

Wall min/avg/max = 58.157/62.298/66.068 ms. The maximum individual frame occurs at the preset-7 anchor. It sets the pass cadence result; a window average is not used as the peak.

### Per-preset table

All eight presets are visited and the serial log wraps from 8 back to 1. Rows use fully held windows, ranked by particle-draw cost. Blended pixels are not counted by this timing image. Per-preset cadence buckets, including transitions, are all green.

| Preset | Held windows | Representative frames | Blended px/f | Draw ms/f | Render ms/f | fps |
|---|--:|---|--:|--:|--:|--:|
| 7 | 6 | 1009–1024 | — | 39.508 | 42.516 | 16 |
| 6 | 6 | 929–944 | — | 38.718 | 46.847 | 16 |
| 2 | 6 | 209–224 | — | 32.953 | 35.464 | 16 |
| 4 | 6 | 529–544 | — | 32.841 | 36.245 | 16 |
| 8 | 6 | 1233–1248 | — | 27.929 | 30.380 | 16 |
| 1 | 14 | 1329–1344 | — | 27.643 | 30.076 | 16 |
| 3 | 6 | 369–384 | — | 26.728 | 29.092 | 16 |
| 5 | 6 | 705–720 | — | 17.454 | 19.141 | 16 |

### Per-pixel figures

No complete blended-pixel scope is emitted by the current direct point/column route, so cycles per blended pixel cannot be derived from this timing capture. The display quadrant contains 10,368 pixels. Counter scope costs include their instrumentation overhead.

## Column-ISR / DMA marshaling cost

Peak-frame window, using the separately snapshotted ISR interval:

```
isr_wake          1151.1/f 0.615/1.918/21.656 us 3.53%
  isr_pack         143.9/f 6.490/8.536/11.540 us 1.96%
  isr_dma_submit   143.9/f 0.703/0.943/6.208 us 0.21%
```

- Rows show call rate, min/avg/max CPU time and CPU share. Pack is costlier than DMA submission; submission measures CPU setup.
- The 72-pixel strobed composite is 600 bytes. At the configured 24 MHz SPI clock, `dma::transfer_us` includes byte-framing gaps and budgets 230 us asynchronously; that wire time is not submit CPU time.
- Inclusive wake share is 3.53%, leaving about 60.294 ms of foreground capacity per 62.5 ms display window. The measured inclusive render peak already fits; no cadence speedup is required. Pack and submit nest inside wake and are not added to it. DMA-completion and other interrupts are not measured by these rows.

## Summary ranking

1. `msp_draw_particles` — 79.7% of the frame, 49.677 ms/f.
2. `msp_particle_step` — 4.7% of the frame, 2.959 ms/f.
3. `pov_preserve_half` — 0.2% of the frame, 0.136 ms/f.

The same-board before/after peak has no measured regression and no new spill. Native checks passed 124/125 CTest tests (including 120-frame smoke); the assertion-reachability check for `verify_export` also fails on the untouched baseline. Release canvas/parameter checks pass 1,044 assertions. Shipping Holosphere, Phantasm and rev-1.2 builds pass size/layout gates.

## Caveats

- CYCCNT free-runs, so all foreground scopes include ISR time.
- `filter_blend` parenting depends on its first caller; it is absent from this direct-route timing capture, so no complete blend count is inferred.
- No per-pixel scopes were added. Existing instrumentation cost remains in the timings.
- The shipping selective-O3 regions, ITCM trail/point/lattice entries and cached-flash Plot leaves remain active. The refactor changes cold registration only.
- No dwell-compression knob is used. This 90-second run extends beyond the normal 80-second show entry to confirm wrap, but remains within one 120-second profiling instance. It is not a natural-handoff capture.
- The measured source worktree was clean. A rebased parent contains the independent pullback color update; the matched COM3 baseline predates it.
- The host uploader uses isolated copies under `C:/work/temp/teensy-tools-loopback-20261030`. Their TCP-table listener-state preflight was corrected because it falsely rejected live localhost listeners. Loopback endpoint/greeting checks, pinned board selection and the supported wrapper's locks remain active. Installed vendor tools and measured firmware are unmodified. Byte offsets and original/corrected helper hashes are retained in `loopback-preflight.json` there.

## Harness

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=MindSplatter`, `HS_PROFILE_WINDOW=16`; `just profile MindSplatter` runs the supported build/flash/capture wrapper. Use the 90-second/window-16 reproduce command above for this complete preset comparison. Raw logs retain source, compiler, ELF hashes and immutable artifact paths.
