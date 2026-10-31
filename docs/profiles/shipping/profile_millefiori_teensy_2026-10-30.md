# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with the command below) of clean source `d3bfdbc02c0ec88e14658e6fdae2d3d0d7fb7872`. Raw capture:
`C:/work/temp/ripple-unify-study-20261030/millefiori-arena-lazy-d3bfdbc02.log`. This replaces the earlier 2026-10-30 snapshot.

**Ongoing ≤58 ms target: MET.** Fully reported frames
2–2048: peak **56.720 ms**
at frame 1627; 0 frames above 58 ms;
**0/2047 ongoing display spills**.
Setup frame 1 is 81.618 ms. Full pass peak is
81.618 ms at frame 1, with **1/2048 spills**.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV, flywheel and DMA ISRs live; COM3, serial 6362410 |
| Image | Shipping `profile`: `-Os` base with selective `HS_O3` shader/scan/color helpers; PC=False, deep=False |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Millefiori 288×144, single-entry playlist, source `d3bfdbc02c0e` |
| Method | 32-frame windows; clock sampling **paired**; `HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

`profile`: FLASH code 105,920, data 153,376, headers 8,992; RAM1 variables 313,056, code 32,216, padding 552, free local 178,464; RAM2 variables 520,064, free 4,224 bytes.

`phantasm`: FLASH code 1,039,200, data 749,240, headers 8,680; RAM1 variables 312,736, code 196,600, padding 8, free local 14,944; RAM2 variables 520,064, free 4,224 bytes.

Canonical parser validation: **VALID**, root checksum 0.9 ppm, window [1953, 1984]. Wrapper validation=True; clean artifact=True. 64 complete windows; 8 unreported tail rows excluded.

The root counter’s cycles divided by 600 MHz match measured wall sum within the reported checksum error. Paired wall/cycle reads use one interrupt mask and restore its preceding state.

### Frozen provenance

|Artifact|SHA-256/value|
|---|---|
|Source|`d3bfdbc02c0ec88e14658e6fdae2d3d0d7fb7872`|
|Profile ELF|`de81c4a6e1203762568415f32d3e514200a6ee7cccefd6783dc0288dc7f17cb5`|
|Shipping ELF|`98c06ee6023ba5cc461028d6511dc14339be794dbbc497e1d8d8c483481ff1f9`|
|Capture|`891dcbcc41a13fcc01130fc21214c07370ac7dcf3709a56180b62ad4a9834535`|
|Profile envdump|`b65b3058b70410bc8667000cfa1ff2d180edf58e6bcc9620264c29b2019c967a`|
|Compiler|`GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`|

Artifacts: `C:/work/temp/ripple-unify-study-20261030/artifacts/millefiori-arena-lazy-d3bfdbc02_d3bfdbc02c0e_de81c4a6e120`. ELF and envdump hashes are verified against the archive; source status and diff, PC/deep flags, excluded windows and exact counter trees are retained in the companion JSON.

## Frame cadence

| Scope | Frames | Peak render | >58 ms | Display spills (>62.5 ms) |
|---|---:|---:|---:|---:|
| Setup | 1 | 81.618 ms | 1 | 1/1 |
| Ongoing acceptance | 2–2048 | 56.720 ms, frame 1627 | 0 | 0/2047 |
| Full pass | 1–2048 | 81.618 ms, frame 1 | 1 | 1/2048 |

README cells: peak 🟢 56.720, spilled 🟢 0/2047 (0.00%).
These cells use ongoing frames; setup is reported separately.

Acceptance uses exact per-frame **render**, never wall. Ongoing mean render is
50.432 ms. The ongoing peak has
1.280 ms margin against 58 ms and
5.780 ms against the 62.5 ms display window.
Setup remains included in full-pass figures. No shorter matched prefix replaces
the ALL-ongoing peak.

Nominal display quadrant: 144×72=10,368 pixels. Clip margins shade 146×73=
**10,658 pixels/ongoing frame**; prime frame shades 288×73=**21,024 pixels**.
Both prepare 144 shared sweep-time contexts/frame. `canvas_buffer_wait` is
display synchronization idle and is removed from render.

## Phase-by-phase readout

This is continuous animation; windows describe workload, not distinct presets.
Trees are per-frame averages: milliseconds, million cycles, root-cycle share,
calls/frame and microseconds/call. Parents include children.

### Startup-containing window (frames 1–32)

```text
frame                     62.540ms 37.52M 100.0%
  pov_preserve_half        0.135ms  0.08M  0.2% x0.969 139.26us
  fx_shader_draw          46.246ms 27.75M 73.9% x1 46246.50us
    fx_prepare_column      2.225ms  1.33M  3.6% x144 15.45us
  fx_prepare_frame         0.640ms  0.38M  1.0% x1 639.53us
  fx_advance               0.987ms  0.59M  1.6% x1 986.97us
    palette_cycle_step     0.902ms  0.54M  1.4% x1 902.47us
  fx_timeline_step         0.080ms  0.05M  0.1% x1 80.31us
  canvas_clear             0.098ms  0.06M  0.2% x1 97.53us
  canvas_buffer_wait      13.980ms  8.39M 22.4% x1 13980.47us
```

Wall min/avg/max=43.357/62.540/81.618 ms; mean render=48.560 ms; exact peak render=81.618 ms.

### Ongoing peak window (frames 1601–1632)

```text
frame                     62.457ms 37.47M 100.0%
  pov_preserve_half        0.137ms  0.08M  0.2% x1 136.81us
  fx_shader_draw          52.424ms 31.45M 83.9% x1 52424.28us
    fx_prepare_column      2.239ms  1.34M  3.6% x144 15.55us
  fx_prepare_frame         0.863ms  0.52M  1.4% x1 863.44us
  fx_advance               0.970ms  0.58M  1.6% x1 969.66us
    palette_cycle_step     0.941ms  0.56M  1.5% x1 941.47us
  fx_timeline_step         0.083ms  0.05M  0.1% x1 83.34us
  canvas_clear             0.095ms  0.06M  0.2% x1 94.66us
  canvas_buffer_wait       7.524ms  4.51M 12.0% x1 7523.91us
```

Wall min/avg/max=58.493/62.457/66.422 ms; mean render=54.934 ms; exact peak render=56.720 ms.

### Highest mean shader cost (frames 1953–1984)

```text
frame                     62.386ms 37.43M 100.0%
  pov_preserve_half        0.137ms  0.08M  0.2% x1 137.09us
  fx_shader_draw          52.914ms 31.75M 84.8% x1 52914.00us
    fx_prepare_column      2.253ms  1.35M  3.6% x144 15.65us
  fx_prepare_frame         0.846ms  0.51M  1.4% x1 846.09us
  fx_advance               0.953ms  0.57M  1.5% x1 952.69us
    palette_cycle_step     0.923ms  0.55M  1.5% x1 922.53us
  fx_timeline_step         0.084ms  0.05M  0.1% x1 84.03us
  canvas_clear             0.094ms  0.06M  0.2% x1 93.72us
  canvas_buffer_wait       6.997ms  4.20M 11.2% x1 6996.97us
```

Wall min/avg/max=61.154/62.386/63.739 ms; mean render=55.389 ms; exact peak render=56.165 ms.

### Lowest mean shader cost (frames 737–768)

```text
frame                     62.355ms 37.41M 100.0%
  pov_preserve_half        0.138ms  0.08M  0.2% x1 138.41us
  fx_shader_draw          42.030ms 25.22M 67.4% x1 42030.34us
    fx_prepare_column      2.198ms  1.32M  3.5% x144 15.26us
  fx_prepare_frame         0.596ms  0.36M  1.0% x1 596.25us
  fx_advance               0.953ms  0.57M  1.5% x1 952.84us
    palette_cycle_step     0.924ms  0.55M  1.5% x1 923.97us
  fx_timeline_step         0.088ms  0.05M  0.1% x1 87.56us
  canvas_clear             0.094ms  0.06M  0.2% x1 93.78us
  canvas_buffer_wait      18.101ms 10.86M 29.0% x1 18101.41us
```

Wall min/avg/max=59.143/62.354/65.332 ms; mean render=44.253 ms; exact peak render=46.768 ms.

### Per-pixel figures

Direct pixel shading has no `filter_blend` subtree. These inclusive shader figures include preparation and ISR interruptions.

|Window|Actual pixels/window|Cycles/pixel|µs/pixel|
|---|---:|---:|---:|
|1–32|351,422|2526.7|4.211|
|1601–1632|341,056|2951.3|4.919|
|1953–1984|341,056|2978.8|4.965|
|737–768|341,056|2366.1|3.944|

## Column-ISR / DMA marshaling cost

Frames 1601–1632, ISR interval 2,001,777 µs.

|Scope|Calls/frame|Calls/s|Min/avg/max µs|CPU share|
|---|---:|---:|---:|---:|
|`isr_wake`|1153.09|18433.1|0.473/2.013/28.981|3.71%|
|`isr_pack`|144.12|2304.0|8.755/9.503/12.821|2.18%|
|`isr_dma_submit`|144.12|2304.0|0.715/0.930/1.035|0.21%|

Inclusive wake corresponds to 2.321 ms per 62.5 ms display window.

Frames 1953–1984, ISR interval 1,999,729 µs.

|Scope|Calls/frame|Calls/s|Min/avg/max µs|CPU share|
|---|---:|---:|---:|---:|
|`isr_wake`|1151.94|18433.5|0.473/2.008/24.801|3.70%|
|`isr_pack`|144.00|2304.3|8.670/9.455/12.638|2.17%|
|`isr_dma_submit`|144.00|2304.3|0.578/0.930/1.060|0.21%|

Inclusive wake corresponds to 2.315 ms per 62.5 ms display window.

Wake includes nested pack and submit: never add their shares. Rates use the ISR interval, not the frame wall sum. Submit is CPU setup; async wire transfer, DMA-completion and other interrupts are not measured here. Render already includes interruptions.

## Summary ranking

- `fx_shader_draw`: 48.084 ms/frame, 77.05% of root.
- `fx_prepare_column`: 2.224 ms/frame, 3.56% of root.
- `fx_prepare_frame`: 0.740 ms/frame, 1.19% of root.
- `fx_advance`: 0.949 ms/frame, 1.52% of root.
- `canvas_buffer_wait`: 11.964 ms/frame, 19.17% of root.

### Shared kernel and arena preparation

`core/math/ripple.h` shares the wave equations, support tests, rotation and
radial-transfer math between pullback and IslamicStars. Scheduling and sampling
remain in their adapters. IslamicStars keeps its sequential wave semantics;
Millefiori keeps its approved grouped radial-table path. The shared kernel alone
reduced shipping instruction code by 96 B and FLASH data by 4 B.

Millefiori's prepared frame lives in the existing persistent arena on every
target. One outlined construction helper serves both temporal columns and static
frames. The frame owns its copied context; column shading completes before the
arena slot is replaced. The portable storage option remains opt-in for other
effects.

The exact burst preparer resolves a center once per nonempty burst. In the
2,048-draw native replay, fixed-center evaluations fell from 853,235 to 11,753
and random-center evaluations from 843,695 to 321,935. The unchanged table baker
accounts for the remaining fixed-center evaluations. These counts are not
device-cycle measurements.

For eligible radial-table columns, the surface now leaves exact wave kernels
unprepared. Every table rejection prepares the exact kernels locally and
restarts the grouped evaluator from the original input. Literal endpoints,
missing/disabled tables and custom or derived providers retain eager behavior.
No wave equation, table sample, support threshold or temporal interpolation
changed. The prepared state remains immutable and contains no borrowed context
or mutable cache.

The native replay measured 9,773 transfer rejections among 21,528,576 eligible
table evaluations (0.0454%); these are fixture counts, not device telemetry.
Debug and Release tests pin all rejection reasons, original-input restart and
provider behavior. Full 600-frame RGB16 comparisons are bit-exact with both
stack and arena storage.

Actual final ARM code skips wave construction, center/hash evaluation and sync
on the successful table preparation path. One 584 B clear and two 584 B tuple
copies remain; their removal is an unlanded experiment and is not part of this
report's performance claim. The complete accepted arena configuration measures
56.720 ms ongoing peak. Individual source changes are not assigned isolated
timing gains from the separate-binary captures.

### Controlled captures and storage experiments

Arena-backed frame storage and the stack both reside in DTCM on Teensy 4.0. A code-generation difference is not evidence that either memory region is faster. A controlled same-binary experiment over 4,096 constructor pairs measured 39.327009 µs to stack and 39.328889 µs to arena: 0.004779%, about 1.13 cycles. This found no meaningful destination-address penalty. Earlier separate-build comparisons are not causal evidence of an arena penalty. An arena constructor proxy preserved the factory call but eliminated no additional context copy; tuple-assignment probes retained the large preparation temporary and added initialization. Only measured captures establish performance.

|Capture|Source|Clock sampling|ALL-ongoing peak ms|Peak frame|>58 ms|Spills|
|---|---|---|---:|---:|---:|---:|
|millefiori-combined-1bb25bcb2|`1bb25bcb2c`|independent / unverified|59.614|1872|63|0/2047|
|millefiori-after-e14c265f9|`e14c265f92`|independent / unverified|58.053|1874|2|0/2047|
|millefiori-final-9d2fcbc64|`9d2fcbc642`|paired|58.592|1874|5|0/2047|
|millefiori-controlled-baseline-fff5d338d|`fff5d338d8`|paired|58.032|1874|1|0/2047|
|millefiori-arena-3f1272108|`3f12721084`|paired|60.123|1621|113|0/2047|
|millefiori-center-129ee42da|`129ee42da8`|paired|58.135|1623|4|0/2047|
|millefiori-arena-lazy-d3bfdbc02|`d3bfdbc02c`|paired|56.720|1627|0|0/2047|

Older independent-clock captures remain historical evidence. Current acceptance and controlled comparisons use interrupt-matched wall/cycle samples. Timing from different clock methods is not presented as a clean algorithmic delta.

### Preserved original optimization evidence

The following records the original 50e659aa5 study and its clock method; it does not establish current target acceptance.

The largest measured region remains the shader. The clean equal-coverage
precursor `1921ee80b` peaks at 61.908 ms, frame 1621, over frames 2–2048, with
412 frames above 58 ms and 0 display spills. `50e659aa5` reduces that peak by
**4.420 ms (7.14%)** to 57.488 ms and removes all >58 ms frames. Ongoing mean
falls from 54.746 to
50.325 ms. This is a same-coverage comparison;
older shorter captures are retained in the external comparison JSON and do not
replace any candidate’s ALL-ongoing peak.

### Diagnostic PC sampling and direct carrier writes

The separate precursor PC run `palette_lab_pc_1921ee80b.log` recorded 620,600
ticks over 62,060 ms, including 553,232 render samples. Its report attributes
22,376 samples (4.0% of render samples) to `memcpy`, 8.6% to the generated
palette path, and 27.2% to the inlined Millefiori shader. Statistical PC bins
are not exclusive phase timing or proof that all memcpy belongs to one caller.
The PC capture **fails strict checksum validation** (144.4 ppm, frames 193–224)
and is diagnostic only; it is excluded from target acceptance.

Source change `50e659aa5` writes cached palette carrier fields directly in
`GeneratedPalette::apply_shared_base` instead of assigning a temporary aggregate,
then stores cached Lab and its validity flag. It retains the field, opacity,
palette coordinate, hue/lift defaults, and brightness gain. The clean no-PC
before/after captures establish the net speed improvement; sampling explains
the investigation without claiming an isolated memcpy speedup.

ARM disassembly confirms that the cached-base helper shrinks from 784 to
752 bytes and its local stack allocation from 92 to 12 bytes. Its 36-byte
`memset` and 69-byte `memcpy` calls disappear. Every carrier member is set
explicitly; a dirty-output test catches an omitted lift reset. The final
rendering change is bit-exact over 600 native frames. The palette cache itself
uses two 64-by-3 OKLab buffers (4,608 bytes of node payload) in the existing
arena. Compared with the preceding phase-only renderer, its 600-frame RGB16
difference has maximum 3,229 and mean 3.12 channel codes; ripple geometry is
unchanged.


### Validation results

User approved the radial-table appearance. The renderer retains deferred temporal color composition and one finish per final pixel.

- Native `composed_effect`: 358,777 passed, 0 failed.
- Native `pullback`: 1,363,013 passed, 0 failed.
- Native `shader_chain`: 5,173,555 passed, 0 failed.
- Native `profiling`: 483 passed, 0 failed.
- Native `transformers`: 2,329 passed, 0 failed.
- Native `effects`: 24,584,729 passed, 0 failed.
- Render goldens: 175 passed, 0 failed.
- Cold Phantasm and Holosphere builds: all ten first-party translation units compiled; zero first-party warnings; size and layout gates passed.
- Local validation covers the six relevant suites and rendering goldens; the full matrix runs in remote CI.
- WASM smoke: 120 frames; 0 failures.
- WASM Millefiori 96×20 flat measured stack 1236/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 96×20 temporal measured stack 1252/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 288×144 flat measured stack 1236/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 288×144 temporal measured stack 1252/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.

## Caveats

- This finite capture covers one default continuously animated effect, not every control setting or future frame.
- All fully reported ongoing frames count; setup and full-pass spills remain explicit. Unreported tails are excluded.
- CYCCNT includes ISR and instrumentation; counter parents and children overlap.
- This shipping selective-O3 configuration is not an all-O3 image. No per-pixel scopes were added.
- `HS_PROFILE_EPOCH_REVS=1200` extends the playlist lifetime; it does not reduce per-frame work.
- Color fidelity to the old framebuffer was relaxed; ripple geometry/motion was qualified and the user approved the look.
- Assertion reachability and promoted document bindings pass their checks. The typed document API is exercised by headless smoke; the Daydream compiled-apply consumer is unchanged.

## Harness

```sh
bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'
```

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`, `HS_PROFILE_WINDOW=32`.
The supported wrapper owns the device lock. This report generator performs no
build, upload, install or capture action.
