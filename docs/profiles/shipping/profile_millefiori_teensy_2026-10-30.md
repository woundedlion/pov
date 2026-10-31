# Millefiori on-device profile — Teensy 4.0, segmented mode (2026-10-30, **selective -O3**)

Point-in-time snapshot (regenerate with the command below) of clean source `86497f3b02e0e78909a36907b94845a35e34b948`. Raw capture:
`C:/work/temp/ripple-unify-study-20261030/millefiori-size-86497f3b0.log`. This replaces the earlier 2026-10-30 timing snapshot and records the subsequent ITCM audit.

**Ongoing ≤58 ms target: MET.** Fully reported frames
2–2048: peak **56.200 ms**
at frame 1623; 0 frames above 58 ms;
**0/2047 ongoing display spills**.
Setup frame 1 is 82.100 ms. Full pass peak is
82.100 ms at frame 1, with **1/2048 spills**.

## Setup

| | |
|---|---|
| Hardware | Teensy 4.0 @ 600 MHz, segmented POV, flywheel and DMA ISRs live; COM3, serial 6362410 |
| Image | Shipping `profile`: `-Os` base with selective `HS_O3` shader/scan/color helpers; PC=False, deep=False |
| Driver | `POVSegmented<288,4,480>`, segment 0 master |
| Effect | Millefiori 288×144, single-entry playlist, source `86497f3b02e0` |
| Method | 32-frame windows; clock sampling **paired**; `HS_PROFILE_EPOCH_REVS=1200` |
| Reproduce | `bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'` |

`profile`: FLASH code 105,912, data 153,376, headers 9,000; RAM1 variables 313,056, code 29,848, padding 2,920, free local 178,464; RAM2 variables 520,064, free 4,224 bytes.

`phantasm`: FLASH code 1,037,664, data 749,244, headers 9,188; RAM1 variables 312,736, code 192,632, padding 3,976, free local 14,944; RAM2 variables 520,064, free 4,224 bytes.

Canonical parser validation: **VALID**, root checksum 0.6 ppm, window [1953, 1984]. Wrapper validation=True; clean artifact=True. 64 complete windows; 8 unreported tail rows excluded.

The root counter’s cycles divided by 600 MHz match measured wall sum within the reported checksum error. Paired wall/cycle reads use one interrupt mask and restore its preceding state.

### Frozen provenance

|Artifact|SHA-256/value|
|---|---|
|Source|`86497f3b02e0e78909a36907b94845a35e34b948`|
|Profile ELF|`9f3ff19ed883b65e2f0843b1dd5007667c09c28a5936fdf55cdc59c09fb2a55d`|
|Shipping ELF|`2c709a9716a9c605656b1f2b8f81f3987c42f864d8dc83590776415b4891fa09`|
|Capture|`bc55d7c61c1ca82f9b96b93a5399a678ff3453d7a853f6fbac9ac303f500679d`|
|Profile envdump|`82220a805ab73106cfaae38c360b378b35aa4f8367e8d9658046e2c3ee8482eb`|
|Compiler|`GCC: (Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)) 15.2.1 20251203`|

Artifacts: `C:/work/temp/ripple-unify-study-20261030/artifacts/millefiori-size-86497f3b0_86497f3b02e0_9f3ff19ed883`. ELF and envdump hashes are verified against the archive; source status and diff, PC/deep flags, excluded windows and exact counter trees are retained in the companion JSON.

## Frame cadence

| Scope | Frames | Peak render | >58 ms | Display spills (>62.5 ms) |
|---|---:|---:|---:|---:|
| Setup | 1 | 82.100 ms | 1 | 1/1 |
| Ongoing acceptance | 2–2048 | 56.200 ms, frame 1623 | 0 | 0/2047 |
| Full pass | 1–2048 | 82.100 ms, frame 1 | 1 | 1/2048 |

README cells: peak 🟢 56.200, spilled 🟢 0/2047 (0.00%).
These cells use ongoing frames; setup remains separate.

Acceptance uses exact per-frame **render**, never wall. Ongoing mean render is
49.820 ms. The ongoing peak has
1.800 ms margin against 58 ms and
6.300 ms against the 62.5 ms display window.
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
frame                     62.526ms 37.52M 100.0%
  pov_preserve_half        0.136ms  0.08M  0.2% x0.969 139.90us
  fx_shader_draw          45.516ms 27.31M 72.8% x1 45516.41us
    fx_prepare_column      1.771ms  1.06M  2.8% x144 12.30us
  fx_prepare_frame         0.635ms  0.38M  1.0% x1 635.09us
  fx_advance               1.135ms  0.68M  1.8% x1 1134.66us
    palette_cycle_step     1.034ms  0.62M  1.7% x1 1034.50us
  fx_timeline_step         0.101ms  0.06M  0.2% x1 101.19us
  canvas_clear             0.093ms  0.06M  0.1% x1 93.22us
  canvas_buffer_wait      14.517ms  8.71M 23.2% x1 14517.34us
```

Wall min/avg/max=42.431/62.525/82.101 ms; mean render=48.008 ms; exact peak render=82.100 ms.

### Ongoing peak window (frames 1601–1632)

```text
frame                     62.454ms 37.47M 100.0%
  pov_preserve_half        0.136ms  0.08M  0.2% x1 136.34us
  fx_shader_draw          51.764ms 31.06M 82.9% x1 51763.91us
    fx_prepare_column      1.797ms  1.08M  2.9% x144 12.48us
  fx_prepare_frame         0.861ms  0.52M  1.4% x1 861.31us
  fx_advance               1.066ms  0.64M  1.7% x1 1065.66us
    palette_cycle_step     1.031ms  0.62M  1.7% x1 1031.34us
  fx_timeline_step         0.099ms  0.06M  0.2% x1 99.19us
  canvas_clear             0.089ms  0.05M  0.1% x1 88.69us
  canvas_buffer_wait       8.055ms  4.83M 12.9% x1 8055.38us
```

Wall min/avg/max=58.611/62.453/66.349 ms; mean render=54.399 ms; exact peak render=56.200 ms.

### Highest mean shader cost (frames 1953–1984)

```text
frame                     62.390ms 37.43M 100.0%
  pov_preserve_half        0.136ms  0.08M  0.2% x1 136.19us
  fx_shader_draw          52.295ms 31.38M 83.8% x1 52294.56us
    fx_prepare_column      1.803ms  1.08M  2.9% x144 12.52us
  fx_prepare_frame         0.845ms  0.51M  1.4% x1 845.44us
  fx_advance               0.933ms  0.56M  1.5% x1 932.69us
    palette_cycle_step     0.899ms  0.54M  1.4% x1 899.06us
  fx_timeline_step         0.097ms  0.06M  0.2% x1 97.31us
  canvas_clear             0.090ms  0.05M  0.1% x1 89.59us
  canvas_buffer_wait       7.601ms  4.56M 12.2% x1 7601.16us
```

Wall min/avg/max=61.085/62.390/63.707 ms; mean render=54.789 ms; exact peak render=55.524 ms.

### Lowest mean shader cost (frames 737–768)

```text
frame                     62.330ms 37.40M 100.0%
  pov_preserve_half        0.139ms  0.08M  0.2% x1 139.44us
  fx_shader_draw          41.352ms 24.81M 66.3% x1 41352.19us
    fx_prepare_column      1.779ms  1.07M  2.9% x144 12.36us
  fx_prepare_frame         0.604ms  0.36M  1.0% x1 603.84us
  fx_advance               0.975ms  0.59M  1.6% x1 975.06us
    palette_cycle_step     0.945ms  0.57M  1.5% x1 944.56us
  fx_timeline_step         0.111ms  0.07M  0.2% x1 110.69us
  canvas_clear             0.088ms  0.05M  0.1% x1 88.22us
  canvas_buffer_wait      18.687ms 11.21M 30.0% x1 18687.03us
```

Wall min/avg/max=59.096/62.329/65.426 ms; mean render=43.643 ms; exact peak render=46.238 ms.

### Per-pixel figures

Direct pixel shading has no `filter_blend` subtree. These inclusive shader figures include preparation and ISR interruptions.

|Window|Actual pixels/window|Cycles/pixel|µs/pixel|
|---|---:|---:|---:|
|1–32|351,422|2486.8|4.145|
|1601–1632|341,056|2914.1|4.857|
|1953–1984|341,056|2944.0|4.907|
|737–768|341,056|2328.0|3.880|

## Column-ISR / DMA marshaling cost

Frames 1601–1632, ISR interval 2,001,717 µs.

|Scope|Calls/frame|Calls/s|Min/avg/max µs|CPU share|
|---|---:|---:|---:|---:|
|`isr_wake`|1153.06|18433.2|0.488/2.041/30.855|3.76%|
|`isr_pack`|144.12|2304.0|8.755/9.501/12.301|2.18%|
|`isr_dma_submit`|144.12|2304.0|0.770/0.933/9.050|0.21%|

Inclusive wake corresponds to 2.353 ms per 62.5 ms display window.

Frames 1953–1984, ISR interval 1,999,831 µs.

|Scope|Calls/frame|Calls/s|Min/avg/max µs|CPU share|
|---|---:|---:|---:|---:|
|`isr_wake`|1151.97|18433.1|0.488/2.033/28.568|3.74%|
|`isr_pack`|143.97|2303.7|8.665/9.453/12.800|2.17%|
|`isr_dma_submit`|143.97|2303.7|0.648/0.928/1.075|0.21%|

Inclusive wake corresponds to 2.344 ms per 62.5 ms display window.

Wake includes nested pack and submit: never add their shares. Rates use the ISR interval, not the frame wall sum. Submit is CPU setup; async wire transfer, DMA-completion and other interrupts are not measured here. Render already includes interruptions.

## Summary ranking

- `fx_shader_draw`: 47.411 ms/frame, 75.98% of root.
- `fx_prepare_column`: 1.788 ms/frame, 2.87% of root.
- `fx_prepare_frame`: 0.741 ms/frame, 1.19% of root.
- `fx_advance`: 0.973 ms/frame, 1.56% of root.
- `canvas_buffer_wait`: 12.566 ms/frame, 20.14% of root.

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
56.200 ms ongoing peak after the size cleanup. Individual source changes are not assigned isolated
timing gains from the separate-binary captures.

### ITCM audit and recovery

The historical audit rebuilds clean pre-promotion source `9183a1438` with the
same GCC 15.2.1 shipping flags, linker script and build hooks, and compares 16
retained image milestones. Actual ITCM means the ELF `.text.itcm` section.
The stock size tool adds an 8 B FLASH unwind table to its reported RAM1 code;
that reported number is not the actual instruction-section size.

| Configuration | Actual ITCM B | Actual headroom B | Observed ongoing peak |
|---|---:|---:|---|
| Before Millefiori promotion, `9183a1438` | 183,808 | 12,800 | Not in playlist |
| First promoted configuration, `5e5f76362` | 184,288 | 12,320 | 116.378 ms, frames 2-352 |
| Arena target before size cleanup, `d3bfdbc02` | 196,592 | 16 | 56.720 ms, frames 2-2048 |
| Smaller configuration, `86497f3b0` | **192,624** | **3,984** | **56.200 ms**, frames 2-2048 |

The earliest timing used a promotion patch later committed as `5e5f76362`.
Capture coverage and clock sampling differ across the historical study; these
are milestone observations, not isolated timing gains attributable to each
commit. Current acceptance uses the clean, interrupt-matched final capture.

Before cleanup, the full pre-Millefiori-to-target regression was 12,784 B of
actual ITCM. The cleanup recovers 3,968 B; the remaining increase is 8,816 B
(4.8%). All configurations reserve six 32 KiB ITCM banks. RAM1 globals remain
312,736 B and linker free local space remains 14,944 B; this audit concerns
instruction headroom inside the existing bank allocation. Total FLASH falls
from 1,894,400 to 1,796,096 B, a separate 96 KiB reduction.

The largest retained code-growth buckets before cleanup were ripple helpers
(4,000 B), the Millefiori pixel shader (2,776 B), paired palette preparation
(2,290 B) and deferred finishing (1,354 B). The palette bucket contained eleven
identical 126 B endpoint lambdas. One typed non-template endpoint helper now
shares that work. Its compact GCC body remains in ITCM; Clang can inline it.
Nineteen effects remain bit-exact over 600-frame comparisons, and the cached
Millefiori path is unchanged by this deduplication.

Only the exact sequential and grouped ripple wrappers move to hot FLASH. Their
math, guards and O3 settings are unchanged; the radial-table success kernel
stays in ITCM. The native replay counted 308,781 grouped evaluations out of
21,827,584 pixel evaluations (1.4146%), including 299,008 canonical endpoint
evaluations. It counted 6,285 nested sequential fallbacks (0.0288% of all
evaluations). These counts do not describe every control setting.

### Table-disabled control comparison

Two separate clean diagnostic branches change only the initial ripple period
to 30, which disables radial tables. They preserve the same controls and
rendering math, and compare grouped code in ITCM versus FLASH. Both captures
use 70 s, 32-frame windows, the same board and the supported locked uploader.
They are control experiments, not the default target capture.

| Period 30, frames 2-1088 | Grouped in ITCM | Grouped in FLASH |
|---|---:|---:|
| Ongoing peak at frame 914 | 55.989 ms | 56.733 ms |
| Ongoing mean | 47.896 ms | 48.333 ms |
| Ongoing spills | 0/1087 | 0/1087 |

The observed FLASH cost is 0.744 ms at this control's peak and 0.437 ms on its
mean. Both remain below 58 ms over the sampled prefix. Source revisions are
`80edc16cf6b8ecafeb27be124c3beb0eb403bb39` and
`fd3a37d89b2b299ac1208db45dea8a66bd9ab8b9`; the source, retained ELF and capture
hashes pass integrity checks. Other table-disabled parameter combinations and
future effects are not claimed to have the same cost.

The full historical sizes, symbol buckets and counter details are retained in
`C:/work/temp/ripple-ci-firmware-study-20261030/millefiori-full-itcm-history.json`
and the period-30 comparison in
`C:/work/temp/ripple-unify-study-20261030/period30-comparison.json`.

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
|millefiori-size-86497f3b0|`86497f3b02`|paired|56.200|1623|0|0/2047|

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

- Native `composed_effect`: 358,803 passed, 0 failed.
- Native `pullback`: 1,363,013 passed, 0 failed.
- Native `shader_chain`: 5,173,555 passed, 0 failed.
- Native `profiling`: 483 passed, 0 failed.
- Native `transformers`: 2,329 passed, 0 failed.
- Native `effects`: 24,584,729 passed, 0 failed.
- Native `kaleidoscope_melt`: 96 passed, 0 failed.
- Render goldens: 175 passed, 0 failed.
- Cold Holosphere, Phantasm and eight-board builds: all fifteen first-party translation units compiled; zero first-party warnings; all applicable size/layout checks passed.
- Static CTest: 10 passed, 0 failed.
- Local validation covers seven relevant native suites, ten death/static checks and rendering goldens; the full matrix runs in remote CI.
- WASM smoke: 120 frames; 0 failures.
- WASM Millefiori 96×20 flat measured stack 1236/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 96×20 temporal measured stack 1252/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 288×144 flat measured stack 1236/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.
- WASM Millefiori 288×144 temporal measured stack 1252/8192 B, init 1088 B; verified against the rebuilt module. Existing 2,048 B creep gate remains unchanged.

Final source `decb445af3b73f78ddf047b7244b1e2be2e59f6e` includes a GCC-equivalent placement-macro follow-up that restores Clang inlining. The validation record explicitly covers captured source `86497f3b0`; fresh final-module stack and smoke results are reported above.

GridSpace's separate CI repair selects canonical coordinate evaluation to preserve its existing 16-code color bound at a gamut-mapping discontinuity. It saves 192 B of ITCM and measured 13.88% slower in a host replay; this is not a device timing measurement. Millefiori's factored path remains enabled.

## Caveats

- This finite capture covers one default continuously animated effect, not every control setting or future frame.
- All fully reported ongoing frames count; setup and full-pass spills remain explicit. Unreported tails are excluded.
- CYCCNT includes ISR and instrumentation; counter parents and children overlap.
- This shipping selective-O3 configuration is not an all-O3 image. No per-pixel scopes were added.
- `HS_PROFILE_EPOCH_REVS=1200` extends the playlist lifetime; it does not reduce per-frame work.
- Color fidelity to the old framebuffer was relaxed; ripple geometry/motion was qualified and the user approved the look.
- Assertion reachability, document bindings, enum-guard death coverage and stateless fixture layout pass their checks. The typed document API is exercised by headless smoke; the Daydream compiled-apply consumer is unchanged.

## Harness

```sh
bash tools/profile_one.sh Millefiori profile 130 32 '-D HS_PROFILE_EPOCH_REVS=1200'
```

`targets/Profile/Profile.ino`, `HS_PROFILE_TARGET=Millefiori`, `HS_PROFILE_WINDOW=32`.
The supported wrapper owns the device lock. This report generator performs no
build, upload, install or capture action.
