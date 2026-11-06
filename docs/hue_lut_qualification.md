# Hue lookup qualification

Hue modulation retains the 64-value by 16-hue grid and RGB16 channel-error
budgets of 7000 maximum and 256 mean. Each table owns a cell mask. Cells whose
gamut classification changes across nodes, or whose corner RGB16 span exceeds
14000 codes, use procedural hue conversion. Temporal resources, retained
frames and chain samples carry the matching mask with its table.

Cached hue rows store RGB16 nodes and use the existing Q16 pixel interpolation
along value, hue and time. Marked cached cells resolve the interpolated Lab
base procedurally. Near the gamut-admission boundary, the sampler restores
the base from its owning palette endpoints before conversion; the linear-RGB admission limits remain -0.0001 and 1.0001. The separate
canonical-base restoration guard has width 0.001. Neutral hue and alpha retain their separate
contracts. Equivalent retained/current hue tables and masks share the current
view for interior columns; frame endpoints retain their own views.

The optimized bake reuses its pre-mapping linear RGB for cell classification.
Column views retain one mask when both RGB rows share it; distinct masks
remain independently guarded. These changes remove duplicate work while
preserving cell classification, canonical-base restoration and allocation
sizes. Accuracy remains the selected contract. Final integrated checks and stock
peak/spill measurements are recorded below; canonical CI remains required.

Palette-cache views borrow their canonical endpoint palettes as well as Lab
and RGB rows and gamut masks. Sampling must finish before those resources
change. Post-draw replay retains the original interval endpoints and rebinds
a reference cache to them; frame-history capture can overwrite the engine's
previous-palette storage after the draw has consumed its views.

## Retained hue views

Interval preparation compares palette entries, hue pixels and gamut masks.
Matching resources share a hue view through modulation-recipe changes;
differing resources retain procedural fallback. The comparison occurs once
per interval, not per pixel or column.

A same-board GridWorms modulation-edit stress pair over 1088 matched frames
measures mean 54.536 to 52.727 ms and live peak 60.470 to 60.714 ms, with no
spills in either image. Both images include the accurate-color implementation.
These are custom-recipe results, not stock-preset measurements. The qualified
source pair is `5fa946f9` and `3f1c4143`; observed maxima from one pair do not
establish a universal worst-case bound.

## Accuracy checks and their scope

The normal `unit_pullback` module registers two independent procedural-color
comparisons. Their corpus sizes must not be confused with the broader
diagnostic sweeps:

| Check | Corpus | RGB samples | Compared channels | Accurate-control maximum |
| --- | --- | ---: | ---: | ---: |
| `test_pullback_modulation_lut_qualification` | Initial bank, chroma 0.788, TRIADIC, 8 frame targets, 512 values, 128 hues | 524288 | 1572864 | 2184 |
| Broader ordinary diagnostic | 3 sequence offsets, 8 chromas, 3 harmonies, 17 targets, 512 values, 128 hues | 80216064 | 240648192 | 5975 |
| `test_pullback_cached_hue_qualification` | 3 offsets, 4 chromas, 3 harmonies, 11 targets, 3 interior times, 256 values, 128 hues | 38928384 | 116785152 | 5449 |

Offsets are 0, 17 and 83. The cached chromas are 0.314, 0.718, 0.788 and 1;
the ordinary diagnostic also includes 0, 0.05, 0.15 and 0.5. Harmonies are
TRIADIC, COMPLEMENTARY and ANALOGOUS. Restored generated-palette snapshots
are validated before sampling. This establishes valid, restorable states,
not that every sampled coordinate occurs in a stock preset or follows the
same accumulated floating-point trajectory as a forward run.

Clang 23 strict, fast-math and NDEBUG accuracy-control builds satisfy the
unchanged budgets. The cached maximum case mean is below 84.872 channel codes;
a denser ordinary diagnostic reaches maximum 6280. Removing cell admission
produces maximum 17641; disabling canonical-base restoration produces maximum
7914. Both negative controls fail the 7000 limit. A matching float-interpolation
control differs from Q16 by at most one channel code on the expanded cached
corpus. Existing cache tests cover both node axes, hue wrap, temporal endpoints,
Q16 compatibility, owner copies and bake reuse.

The optimized implementation passes 78383939 assertions in each strict and
fast-math qualification run. Strict output matches the accurate control across
38928384 cached RGB samples. Fast-math changes 12372 channels by one code, with
no alpha-bit changes or larger differences in that corpus. Distinct-mask and
raw-classification negative controls fail as intended. These finite corpora
do not prove an error bound over every floating-point input.

## Earlier renderer diagnostics

The earlier maximum of 11040 covered only the initial-bank, chroma-0.788,
TRIADIC, eight-target ordinary corpus. It was not a bound on the full supported
palette domain. On the broader diagnostic, the pre-fix renderer reaches 18152
ordinary and 17641 cached channel error. Strict, fast-math and NDEBUG controls
reproduce those maxima. Maximum per-case mean remains below 184.831 ordinary
and 182.886 cached; these are distinct from aggregate means.

In the uniform diagnostic grids, 15110 of 80216064 ordinary RGB samples
(0.01884%) and 13762 of 38928384 cached samples (0.03535%) have a channel error
above 7000. The accurate control has none. These frequencies describe the
constructed grids, not observed shipping pixels or perceptual visibility.
Errors are absolute unpremultiplied linear RGB16 channel-code differences,
not a perceptual color-distance metric. High-chroma custom states are included.

## Final integrated stock comparisons

The qualified snapshot `492885884` includes the retained-view optimization.
Its control is the current pre-color renderer `8dca144c6`. Each stock preset-0
pair uses the same board, flags, 70-second duration and 32-frame windows. All
six canonical validators pass. Frame 1 is setup; every subsequent captured
frame is retained in the live statistics.

| Effect | Control/accurate live peak ms | Live spills | Live frames at or above 59 ms | Control/accurate setup ms |
| --- | ---: | ---: | ---: | ---: |
| AlienBrain | 31.639 / 35.855 | 0 / 0 | 0 / 0 | 52.600 / 59.585 |
| JewelMelt | 47.712 / 51.121 | 0 / 0 | 0 / 0 | 73.707 / 77.258 |
| AshCloud | 55.686 / 56.982 | 0 / 0 | 0 / 0 | 82.588 / 89.998 |

Each side contains 1087 live frames. JewelMelt and AshCloud each have one
setup spill in both images. Mean live times, as diagnostics, are
26.916/30.274, 41.860/43.974 and 48.633/49.844 ms respectively. These results
measure the final accuracy cost against the current renderer. They do not
establish a bound for every preset or custom recipe.

The shipping gate passes with RAM1 code 189016 to 192328 bytes and FLASH code
1172240 to 1180440 bytes. RAM1 variables remain 313760 bytes, FLASH data
790784 bytes, and free DTCM 13920 bytes. Both images use six ITCM banks; the
accurate image retains 4280 bytes of code headroom. RAM2 variables/free remain
520064/4224 bytes. No budgets were raised.

The snapshot passes 77 native checks and 162586557 assertions, plus full
Release WASM smoke and 166 render goldens. ABI catalogs travel with the
qualified module. Canonical CI covers the remaining platform and math-mode
matrix.

Two further targeted JewelMelt trials retain the 51.121 ms accurate control.
Moving the guarded sampler to ITCM raises live peak to 53.205 ms. A two-tier
guard passes numerical qualification but raises peak to 51.626 ms. Both remain
unlanded. The selected implementation is the best measured candidate in these
trials; the search does not establish that its cost is unavoidable.

## Measured optimization

Fresh pairs compare accurate control
`0afb0a264056787674e0cb43ab332fdce8872975` with work-elimination candidate
`5fa946f9539a4cb6916b8a5a12aee80c56d8bbb9`. All six captures use the same
Teensy 4.0, serial 6362410, at 600 MHz, GCC 15.2.1, the shipping selective-O3
Profile configuration, 288x144, four synchronized segments, 480 RPM, preset 0,
70 seconds and 32-frame windows. Each image contributes 1088 matched frames.
Every canonical validator exits zero and ends in VALID.

Observed peak frame time and spills are the primary performance criteria.
Mean differences remain useful diagnostics of repeated work:

| Effect | Accurate/optimized peak ms | Spills | Accurate/optimized mean ms | Mean change |
| --- | ---: | ---: | ---: | ---: |
| AlienBrain | 60.395 / 59.726 | 0 / 0 | 30.775534 / 30.327124 | -0.448410 ms (-1.4570%) |
| JewelMelt | 78.431 / 77.689 | 1 / 1 | 44.616056 / 44.226133 | -0.389923 ms (-0.8740%) |
| AshCloud | 86.774 / 86.140 | 1 / 1 | 50.455200 / 49.984269 | -0.470931 ms (-0.9334%) |

JewelMelt and AshCloud spills occur at startup; all steady captures have zero
spills. The 64-frame block intervals for paired mean differences are
[-486.987, -409.294], [-426.125, -352.759] and [-579.485, -361.672] microseconds,
respectively. These single-pair intervals describe observed paired variation,
not independent-run drift. The comparisons isolate the optimization against
the accurate control. They do not measure the final merged renderer against
the earlier inaccurate renderer.

The restored Deep `gamut_clip` scope preserves all nine allocated sections'
bytes, addresses and sizes in the shipping image and in each of the three
normal Profile images. The normal images use identical capture flags; their
only full-ELF difference is a nonallocated `__rtc_localtime` symbol value.
Their recorded normal timing remains applicable. A separate Deep build proves
that the restored counter is registered and referenced by the out-of-gamut
bake branch. This is static target-code reachability, not an observed device
counter count. Other `gamut_clip` sites retain the same label.

## Historical original-to-accurate study

The earlier six-capture study compared the original `5a5cd76f` runtime with
the accurate control, using the same validated ITCM clock sampler in each
image. Every individual pair used the same board and settings. AlienBrain ran
on serial 6362410; AshCloud and JewelMelt ran on serial 6362390.

| Effect | Original/accurate peak ms | Spills | Original/accurate mean ms |
| --- | ---: | ---: | ---: |
| AlienBrain | 52.495 / 60.374 | 0 / 0 | 26.867842 / 30.777451 |
| AshCloud | 83.198 / 86.729 | 1 / 1 | 48.995972 / 50.458306 |
| JewelMelt | 73.935 / 78.517 | 1 / 1 | 41.718107 / 44.618477 |

The AshCloud and JewelMelt spills are startup-only. These historical pairs and
the fresh optimization pairs remain separate. Subtracting results across
their different runs or boards would not produce a measured final cost.

## Memory, ABI and remaining integration checks

| Region | Original 5a bytes | Accurate control bytes | Optimized bytes |
| --- | ---: | ---: | ---: |
| RAM1 code | 192216 | 193656 | 195496 |
| RAM1 variables | 313760 | 313760 | 313760 |
| FLASH code | 1173688 | 1178960 | 1180800 |
| FLASH data | 791624 | 791624 | 791624 |
| DTCM stack free | 13920 | 13920 | 13920 |

All isolated shipping gates pass without raising caps. The optimized image
has 1112 bytes below the 196608-byte, six-bank ITCM ceiling. All 46 native
footprint rows pass with unchanged optimized allocation sizes. The accurate
control's smallest native margin is 45986 bytes. Its focused Release WASM
draw peaks for AlienBrain/AshCloud/JewelMelt are 464/964/464 bytes; initialization
peaks at 1120 bytes, below 2048. These are source-specific results. The final integrated study above records current memory, smoke, stack and
device results; canonical CI supplies the remaining math-mode matrix. Isolated
deltas cannot establish merged headroom or performance.

Catalogs regenerated from the merged source change only generated-palette
prepared-state size: native 88 to 96 bytes and WASM32 64 to 68 bytes. They must
travel with the module generated from that source.

Windows Clang 23 exact-node probes retain a separate known limitation:
`test_palette_rgb_cache_matches_hue_rotation_table` has three one-code equality
failures at `-O2` and seven at `-Os`, with fast-math and finite-math-only disabled.
The same failures reproduce in the original pre-review `0ba616d4fe71` and
current pre-fix source. The accurate and optimized `-Os` controls add none.
Assertions and tolerances remain unchanged. These controls do not establish
the result of canonical Linux Clang 22 fast-math CI, which remains required.

HankinSolids and IslamicStars do not call the affected hue-table or cached-RGB
palette paths in the tested full sequences. Native counters observe zero
bakes and samples over 4000 frames of each effect; the two-frame AlienBrain
positive control observes one hue bake, 3840 hue samples, two RGB-cache bakes
and 3760 RGB-cache samples. Actual target full-cycle captures independently
show no composed prepare/bake/draw scope calls, while their own timeline scopes
execute for all 3648 Hankin and 3328 Islamic frames. This does not exclude
indirect costs from global code placement, caches or memory. The final merged
peak/spill comparisons above cover AlienBrain, JewelMelt and AshCloud; no additional Hankin/Islamic capture is scheduled for this change.

[Machine-readable qualification and artifact hashes](hue_lut_qualification.json)
preserve each study's source and capture provenance. Profile clock annotation
and reach instrumentation are excluded from the production patch.
