# ShaderBall 59 ms quest plan

Status: **complete.** The 12-preset roster landed at `43bd0b41`; source
`c06c5fbc` passes the two-board shipping gate.

The shipping gate is a peak render time below **59.00 ms** for every authored
ShaderBall preset and every transition frame in the selective-O3 Phantasm
image. Development targets **56.00 ms** so normal orientation, interrupt, and
measurement variation cannot turn a nominal pass into a field failure. Final
acceptance is a complete 12-preset on-device cycle on each board with zero
spilled frames, not a host benchmark or a single favorable window.

This plan supersedes the completion claim for the former 23-preset bank in
`shaderball_optimization_plan.md`. That document remains the historical
optimization ledger and exact-oracle record.

## Roster reduction and numbering

The current bank removes former presets 5, 6, 7, 9, and 10. Historical
profiles and experiment rows retain their captured indices; they are never
relabelled as measurements of the renumbered programs.

| Former preset | Current preset |
|---:|---:|
| 0 | 0 |
| 1 | 1 |
| 2 | 2 |
| 3 | 3 |
| 4 | 4 |
| 5 | Retired |
| 6 | Retired |
| 7 | Retired |
| 8 | 5 |
| 9 | Retired |
| 10 | Retired |
| 11 | 6 |
| 12 | 7 |
| 13 | 8 |
| 14 | 9 |
| 15 | 10 |
| 16 | 11 |

## Historical evidence and baseline integrity

The 2026-08-13 shipping report is the initial symptom map. It covers all 17
presets, wraps 16 to 0, agrees with the DWT root counter within 0.2 ppm, and
reports a 147.35 ms peak with 664 of 1,472 frames spilled. Five buckets are
green and twelve are red.

The initial report pair was not a matched compiler comparison: its shipping
log named `f0821319`, its O3 log named `d0848a94`, and both reports named
`6f2e7e77`. The final reports replace that ambiguity with four provenance-
attested captures from source `c06c5fbc`, two boards per configuration.

| Preset | Authored topology | Peak, ms | Clean shader, ms | Gap to 59 ms |
|---:|---|---:|---:|---:|
| 16 | Gnomonic dodecahedral vector-noise + mirror grid | 147.35 | 133.79 | 88.36 |
| 0 | Glitch wave-shear grid | 146.53 | 35.50 | 87.54* |
| 13 | Sinusoidal curl lattice, fine field | 107.89 | 85.34 | 48.90 |
| 14 | Sinusoidal curl lattice, coarse field | 97.91 | 84.98 | 38.92 |
| 15 | Stereographic prism polar-wave lattice | 97.44 | 68.77 | 38.45 |
| 7 (retired) | Dodecahedral grid + mirror | 82.27 | 69.73 | 23.28 |
| 8 | Peirce dodecahedral grid | 81.58 | 67.64 | 22.59 |
| 9 (retired) | Dodecahedral noise grid | 79.59 | 51.93 | 20.60* |
| 11 | Gnomonic dodecahedral wave/mirror grid | 79.08 | 66.34 | 20.09 |
| 12 | Gnomonic affine lattice contour | 78.86 | 43.91 | 19.87* |
| 10 (retired) | Dodecahedral noise lattice + mirror | 78.25 | 56.80 | 19.26 |
| 6 (retired) | Kaleidoscope noise grid + edge fade | 75.36 | 62.68 | 16.37 |
| 5 (retired) | Peirce kaleidoscope lattice | 49.35 | 30.21 | pass |
| 4 | Bonne kaleidoscope lattice + mirror | 49.30 | 37.94 | pass |
| 1 | Kaleidoscope twin-wave + inner mirror | 47.81 | 19.63 | pass |
| 3 | Gnomonic glitch grid + mirror | 33.36 | 21.33 | pass |
| 2 | Gnomonic kaleidoscope grid + mirror | 32.89 | 21.13 | pass |

`*` The clean hold is already below 59 ms. The cycle bucket inherits an
expensive adjacent endpoint or contains a phase spike, so a fixed-preset run
and per-frame transition trace must distinguish the two before code changes.

The attested full-roster image currently reports 195,640 of 196,608 ITCM
bytes, only **968 bytes** of headroom. RAM2 has 4,224 bytes free, only 128 bytes
above its 4 KiB allocator floor. RAM1 leaves 12,800 bytes for locals, 512 bytes
above the stack floor. Flash has ample room. Large lookup tables, template
clones, and persistent fields therefore require an equal-or-larger reclaim or
reuse of existing bounded storage.

## Historical measurement campaign

### Establish a trustworthy baseline

1. Build and capture shipping selective-O3 and global-O3 from one clean
   current-master worktree through `tools/profile_one.sh`, pinned to one board.
2. Use 16-frame windows, fast-cycle choreography, and a 1,400-revolution epoch.
   Validate all 17 markers, the 16-to-0 wrap, no reset, fresh frame numbers,
   artifact hashes, log mtimes, and DWT/wall exactness.
3. Record the full Phantasm FLASH, ITCM, RAM1, RAM2, and stack gates from each
   build. The single-effect profile image is never a memory acceptance proxy.
4. Capture fixed 70-second, 32-frame shipping runs for presets 0, 6-16. Fixed
   runs expose orientation and animation phases hidden by compressed holds.
5. Repeat any apparent pass within 2 ms of the gate. A result is accepted only
   when the worst valid frame from both runs is below 59 ms.

### Deep stage attribution

Use the existing low-distortion ShaderBall DWT accumulator, extended only when
the current counters cannot answer a decision. Capture per-preset totals for:

- outer camera and sphere lookup;
- lens, including reflection count distribution and maximum;
- surface noise by basis, vector/scalar evaluations, and integrator steps;
- projection and projection metadata;
- planar warp 1 and 2 separately, including noise evaluations;
- source, material, value transfer, coverage, and color;
- hue/value field build versus lookup;
- shader dispatch and scan-loop overhead; and
- endpoint preparation, clear, buffer wait, and ISR/DMA share.

Counters are plain CYCCNT deltas accumulated at stage boundaries. No registry
scope, serial output, or branch is added per pixel to a timing image. Diagnostic
subtraction builds may temporarily replace one stage with an identity or a
recorded input, but those images only price a stage; they cannot prove a ship
result.

For every red bucket, retain a per-frame trace around entry, hold, exit, and
wrap. This resolves transition ownership and catches rare orientation peaks
that window means conceal.

### Generated-code examination

Archive the exact shipping and O3 ELFs and use the same GCC 15.2.1
`arm-none-eabi-nm`, `objdump`, `readelf`, and `size` binaries recorded in their
provenance. For every hot program:

1. Map the `Scan::Shader::draw -> FrameShader -> InversePipeline` call chain,
   demangled symbol size, section, address, and call target.
2. Disassemble all 16 program descriptors and the shared noise, lens,
   projection, warp, source, material, and color leaves.
3. Count indirect calls, flash calls per sample, stack spills, table loads,
   conditional branches, transcendental library calls, divides, square roots,
   conversions, and duplicated template bodies.
4. Compare shipping and O3 instruction topology rather than assuming O3 is
   faster. The historical global-O3 Peirce/dodecahedral regression makes code
   layout a measured variable.
5. Produce targeted `-Os`, `-O2`, `-O3`, and `-O3 -ffast-math` object variants
   for isolated leaves. Inspect assembly before spending device time.
6. Correlate instruction deltas with DWT stage deltas. A smaller symbol or
   lower synthetic cycle count is not accepted unless the live driver improves.

The analysis ledger records ELF hashes, symbol bytes, section movement, full
roster bank delta, device frame delta, and the decision for every experiment.

## Optimization levers and expected wins

The ranges below are hypotheses for the presets that exercise each path. They
overlap and are not additive. Each lever is retained only after a fixed-preset
A/B and the full-roster memory gate.

| Lever | Primary presets | Expected frame win | Cost and decision rule |
|---|---|---:|---|
| Bake surface Curl into a coarse sphere field and interpolate tangent vectors | 13, 14 | 30-55 ms | Reuse bounded ShaderBall storage or arena memory; preserve poles, seam, unit length, and angular-error oracle |
| Bake or fuse projected vector-noise warp samples with scan-coherent interpolation | 16 | 25-60 ms | Highest-value hypothesis; field must follow the animated coordinate frame without visible swimming |
| Share multi-output noise traversal, gradients, octave weights, and prepared phase | 6, 7, 9, 10, 13, 14, 16 | 5-25 ms | Prefer arithmetic removal over ITCM promotion; retain scalar exact oracle |
| Row/ring-coherent sphere, frame, and projection recurrence | 7-16 | 4-18 ms | Replace repeated matrix/trig setup with incremental coordinates; periodically renormalize and test drift |
| Branchless bounded polyhedral chamber fold specialized by solid | 7-11, 16 | 7-20 ms | Use measured reflection maxima; exact fallback and dense boundary oracle remain |
| Specialized planar-warp chains with identity/zero terms erased | 11, 12, 15, 16 | 5-25 ms | Fuse the authored two-stage chain without changing stage order or GUI admission |
| Prism polar-chart and wave-shear algebra reduction | 15 | 10-25 ms | Hoist frame constants, replace generic remainder/trig where exact periodic folds allow it |
| Projection-specific kernels and metadata pruning | 8, 12-16 | 3-15 ms | Peirce fast path already exists; remove unused edge/fade/layout work per topology |
| Move program dispatch outside the per-pixel loop | all, especially 13-16 | 1-6 ms | Clone only the scan shell or group compatible programs; reject if ITCM growth exceeds the reclaim |
| Rebalance inlining and flash/ITCM placement per hot leaf | all red presets | 2-15 ms | Favor small hot loops in ITCM and cold/exact oracles in flash; global O3 is not a shortcut |
| Eliminate redundant normalize, clamp, wrap, and coordinate conversions across fused stages | 6-16 | 3-12 ms | Prove range invariants and finite behavior with dense host scans |
| Reuse and extend the persistent value/hue field for compatible source/color work | 6, 9, 10, 13, 14, 16 | 5-20 ms | No second large field; quantify framebuffer error against the exact path |
| Pack small hot tables or constants into reclaimed DTCM and arrange cache-line locality | noise and color families | 1-8 ms | RAM2 is effectively closed; every byte is charged against roster/stack floors |
| Cortex-M7 DSP/fixed-point kernels for mirror, lattice, interpolation, and color | remaining holdouts | 5-20 ms | Higher-risk fallback; exact float oracle, saturation, and seam tests required |
| Approximate transcendental or reciprocal kernels with bounded domains | remaining projection/lens holdouts | 2-15 ms | Exhaust algebraic removal first; record max/mean error and boundary behavior |
| Re-author a preset with a cheaper equivalent topology | irreducible holdout only | unbounded | Last resort; requires explicit rendered-look review and does not redefine the 59 ms gate |

### Historical preset-family attack order

1. **Preset 16:** split dodecahedral lens, vector-noise warp, projection,
   mirror, and shared shading. It needs about a 2.5x peak-render improvement,
   so a local compiler tweak cannot close it. Attack sample-count reduction and
   scan-coherent field construction first.
2. **Presets 13/14:** identical Curl cost at two scales strongly suggests a
   fixed per-pixel noise/integrator bill. Price direct samples, field build,
   interpolation, normalization, projection, and color independently.
3. **Preset 15:** isolate prism fold, polar chart, wave shear, projection, and
   lattice/color work. It has enough intrinsic cost to need a structural win.
4. **Dodecahedral family 7-11:** one fold/code-layout improvement can close
   several presets. Use preset 5 as the cheap Peirce control and compare 7
   versus 8 to price projection separately from the shared solid lens.
5. **Preset 6 and near-gate preset 10:** shared noise and color cleanup should
   close them after the structural work.
6. **Presets 0, 9, and 12:** re-evaluate after their expensive neighbors pass.
   Optimize only an intrinsic fixed-run spike, not a parser ownership artifact.

## Experiment protocol

Every logical experiment gets its own commit in the isolated worktree.

1. Record the source SHA, compiler flags, exact ELF hash, symbol/section delta,
   Phantasm bank usage, and expected affected presets.
2. Run native mathematical, ShaderBall oracle, framebuffer, serialization, and
   topology tests before flashing.
3. Capture the target fixed preset on the same board and configuration as its
   baseline. Keep only a repeatable improvement outside run noise.
4. Immediately reject correctness drift, a memory-floor failure, or a speed
   regression. Do not stack speculative changes to hide a failed A/B.
5. After each accepted shared-kernel change, recheck a cheap control preset and
   one transition so optimization does not merely move the red bucket.
6. Rebase and land accepted changes fast-forward-only. Rebuild the simulator
   after landing; profiling reports remain working-tree artifacts.

The quest ledger will use this shape:

| SHA | Experiment | Preset/run | Peak before -> after | Ship/O3 symbol delta | Phantasm ITCM/RAM delta | Tests | Decision |
|---|---|---|---:|---:|---:|---|---|

## Quest outcome and experiment ledger

The final current-source shipping captures on COM3 and COM4 each contain 648
16-frame windows, all 12 presets, an 11→0 wrap, no reset, and zero spills.
Both peak at 55.69 ms. The matched global-O3 captures each contain 558 windows;
five presets are red, each board spills 1,440/8,928 frames, and the worst peak
is 84.64 ms. The full per-preset evidence is in the
[shipping](profiles/shipping/profile_shaderball_teensy_2026-08-14.md) and
[global-O3](profiles/O3/profile_shaderball_teensy_2026-08-14.md) reports.

`N/R` means the per-experiment artifact was not retained. Cumulative timing is
shown only where it is the narrowest honest attribution; it is not presented
as an isolated A/B. All symbol sizes are bytes. `RAM Δ0` means the measured
RAM1/RAM2 floors were unchanged.

| SHA | Experiment | Preset/run | Peak before -> after | Ship/O3 symbol delta | Phantasm ITCM/RAM delta | Tests | Decision |
|---|---|---|---:|---|---|---|---|
| `43bd0b41` | Retire former presets 5, 6, 7, 9, and 10; renumber retained bank | Current 0–11 full cycle, dual boards | Historical pre-reduction 58.23 -> final-source 55.69 ms; isolated A/B N/R | Final-source ship named symbols 41,076 B; per-commit delta N/R | Cumulative current vs pre-reduction: -1,296 B ITCM, RAM Δ0, flash -5,416 B | Phantasm gate; two validated 648-window shipping cycles | Keep; current bank passes with 3.31 ms headroom |
| `c06c5fbc` | Compare source-matched global O3 with shipping selective O3 | Current 0–11 full cycle, dual boards | Ship 55.69 -> O3 84.64 ms | O3 - ship named symbols: +1,172 B flash / +976 B ITCM (+2,148 B total); whole profile +13,000/+12,352 B | Phantasm ITCM/RAM Δ0 | Four validated cycles; O3 has 5 red presets and 1,440/8,928 spills per board | Reject global O3; keep selective O3 |
| `387c046f` | Derive Simplex curl analytically | 13/14 fixed | Chain: 107.89/97.91 -> 48.41/49.04 ms; isolated A/B N/R | Ship shared curl leaf: N/R; O3: N/R | Per-commit N/R | Native curl/oracle suite | Keep |
| `d7e9303a` | Specialize the Simplex curl surface path | 13/14 fixed | Chain result above; isolated A/B N/R | Ship specialized leaf: N/R; O3: N/R | Per-commit N/R | Native ShaderBall oracle | Keep |
| `e367a443` | Fuse the Simplex Euler surface step | 13/14 fixed | Chain result above; isolated A/B N/R | Ship fused leaf: N/R; O3: N/R | Per-commit N/R | Native ShaderBall oracle | Keep |
| `ff2ea771` | Skip hidden palette rebakes | Transition preparation | N/R | Ship palette cycler: N/R; O3: N/R | Per-commit N/R | Palette/native suite | Keep |
| `b92ac517` | Defer hidden palette boundary rebakes | Transition preparation | N/R | Ship palette cycler: N/R; O3: N/R | Per-commit N/R | Palette/native suite | Keep |
| `67577f49` | Cache spherical hue noise | Noise-colored presets | N/R | Ship hue-field builder/sampler: N/R; O3: N/R | RAM2 Δ0; reused bounded field storage | Noise-field and framebuffer oracles | Keep |
| `8cb5bf5a` | Couple Simplex vector channels | 16 fixed | Chain: 147.35 -> about 55.8 ms after `a2e39bc0`; isolated A/B N/R | Ship vector-noise leaf: N/R; O3: N/R | Per-commit N/R | FastNoise and ShaderBall oracles | Keep |
| `36265e3e` | Bypass redundant palette RGB lookup | Shared color path | N/R | Ship color leaf: N/R; O3: N/R | Per-commit N/R | Color and ShaderBall oracles | Keep |
| `357be008` | Split the hue-field sampler | Hue-noise presets | N/R | Ship sampler: N/R; O3: N/R | RAM2 Δ0 | Hue-field and framebuffer oracles | Keep |
| `1a9644ee` | Pin preset 0 pipeline in ITCM; move cold ShaderBall frame/setup paths to flash | 0 fixed | Flash layout -> 25.8–27.0 ms | Ship p0: 0x298 flash -> ITCM; O3: N/R | Net per-commit N/R | Native suite, Phantasm gate, device fixed run | Keep |
| `a079e9e3` | Temporarily pin preset 16 pipeline in ITCM | 16 fixed | No repeatable useful win | Ship p16: flash -> ITCM; O3: N/R | Net per-commit N/R | Native suite, Phantasm gate, device A/B | Superseded by arithmetic work and `304caf17` |
| `a2e39bc0` | Prepare vector-warp trig and loop offsets | 16 fixed | 147.35 -> about 55.8 ms cumulative | Ship p16: specialized body delta N/R; O3: N/R | Per-commit N/R | Native ShaderBall oracle, device fixed run | Keep |
| `304caf17` | Pin preset 10; return preset 16 to flash; move cold handlers to flash | 10/16 fixed | p10 about 82 -> 51.40 ms; p16 stays about 55.8 ms | Ship p10: 0x42c flash -> ITCM; p16: ITCM -> flash; O3: N/R | Net per-commit N/R | Native suite, Phantasm gate, device A/B | Keep |
| `056ad3f3` | Prepare exact wave-shear direction | 11 fixed | 66.34 -> about 44.6 ms | Ship p11 body delta: N/R; O3: N/R | Per-commit N/R | Native ShaderBall oracle, device fixed run | Keep |
| `ddc13795` | Trade preset 0 ITCM residence for preset 8 | 0/8 fixed | p8 about 76 -> 50.9 ms; p0 regressed in flash | Ship p0: ITCM -> flash; p8: 0x2b4 flash -> ITCM; O3: N/R | Placement swap; net N/R | Native suite, Phantasm gate, device A/B | Superseded by `7fe8dc4c` |
| `7fe8dc4c` | Keep presets 0, 8, and 10 in ITCM; move cold setup code to flash | 0/8/10 fixed | 25.83/50.94/51.40 ms | Ship p0/p8/p10: 0x298/0x2b4/0x42c in ITCM; O3: N/R | Measured final-class footprint: 196,440 B ITCM; RAM floors unchanged | Full hook, Phantasm gate, fixed device runs | Keep |
| `70e31a4c` | Stop discarded transition-source runtime after midpoint | Full transition cycle | p10 transition 64.51 -> about 64.5 ms; no isolated timing win | Ship `draw_frame`/runtime control: N/R; O3: N/R | Per-commit N/R | Transition/runtime tests, full device cycle | Keep; exact work removal |
| `cd79895e` | Update only the visible transition palette | 9->10->11 | p10 transition about 64.6 -> 50.2 ms | Ship palette transition path: N/R; O3: N/R | RAM Δ0; ITCM N/R | Palette/transition tests, dual targeted runs | Keep |
| `b31cda1c` | Resume normal choreography after profile preset selection | Harness, start at 9 | Timing N/A | Profile-only symbol delta: N/R; O3: N/A | Phantasm Δ0 | Profile harness tests, targeted device capture | Keep; measurement only |
| `58335255` | Pin preset 11; move cold boot/control code to flash | 9->10->11 | p11 about 77 -> 50.99 ms; p9 regressed about 58.7 -> 65–67 ms | Ship p11: 0x2d4 flash -> ITCM; cold flash moves: 0x3c + 0x44 + 0x40 + 0x214; O3: N/R | +0x2d4 hot / -0x2d4 cold; measured 196,440 B ITCM, 168 B free; RAM floors unchanged | Effect/pov-sync tests, 71-test hook, Phantasm gate, dual device run | Keep; p9 required follow-up |
| `3af89af1` | Isolate preset-9 shader on a 1 KiB flash boundary | 9->10->11, dual boards | 65–67 -> 59.32/60.12 ms | Ship p9: 0x37c flash body aligned 1024; O3: N/R | ITCM/RAM Δ0; flash padding only | Full hook, dual targeted runs | Superseded; missed gate on both-board maximum |
| `0125443d` | Page-align preset-9 shader | 9->10->11 targeted; final full cycle | 60.12 -> about 42.1 ms targeted p9 phase; final full-cycle p9 58.23 ms | Ship p9: 0x37c at `0x60035000`, aligned 4096; O3: N/R | ITCM/RAM Δ0; flash code before p4 wrapper 427,496 B | Full hook, dual targeted runs, dual full cycles | Keep; 0.77 ms final headroom |
| `c5cb0bb4` | Page-align preset-4 shader | 4->5->6->7 targeted; final full cycle | 65.98 -> 45.69 ms targeted sequence; final p4/p5 buckets 41.64/42.31 ms | Ship p4: 0x300 at `0x60036000`, aligned 4096; O3: N/R | ITCM/RAM Δ0; final flash 427,976 B | Full hook, dual targeted runs, dual 648-window full cycles | Keep; closes gate |
| `94fbacdf` | Replace curl normalization with reciprocal square root | 13/14 fixed | No accepted improvement; exact A/B N/R | Ship curl leaf: N/R; O3: N/R | N/R | Native oracle and device A/B | Reject; numerical/benefit tradeoff |
| `c752085c` | Add frame-preparation attribution probes | Diagnostic image | Acceptance timing N/A | Diagnostic-only instrumentation: N/R | N/R | Diagnostic capture | Reject from shipping; attribution only |
| `e5abfd18` | Prepare the older vector phase/layout | 16 fixed | Superseded by exact offset preparation; A/B N/R | Ship p16 body: N/R; O3: N/R | N/R | Native oracle and device A/B | Reject |
| `c716070b` | Quantize hue lookup to nearest samples | Hue-noise presets | A/B N/R | Ship hue sampler: N/R; O3: N/R | RAM Δ0 | Hue/framebuffer oracle | Reject; quality/benefit tradeoff |
| `894b89d3` | Combine split hue sampler with nearest lookup | Hue-noise presets | A/B N/R | Ship hue sampler: N/R; O3: N/R | RAM Δ0 | Hue/framebuffer oracle | Reject |
| Uncommitted | Compile preset 10 shader with global O3 | 10 fixed | Regression; exact peak N/R | Ship/O3 body delta N/R | N/R | Build, disassembly, device A/B | Reject |
| Uncommitted | Move `SolidBuilder` cold wrappers to flash | Full roster size | Timing N/A | Ship COMDAT: about +4.6 KiB | ITCM benefit defeated by duplicate code; RAM Δ0 | Full-roster build gate | Reject |

The final Phantasm footprint is 195,144/196,608 B ITCM with 1,464 B padding,
314,880 B RAM1 variables with 12,800 B local free, and 520,064 B RAM2
variables with 4,224 B allocator free. The roster reduction reclaimed 1,296 B
of ITCM and 5,416 B of
flash from the attested pre-reduction image; RAM1/RAM2 floors are unchanged.
Historical `O3: N/R` cells still mean the required per-experiment matched
artifact was not retained. The final source-matched compiler comparison is the
explicit `c06c5fbc` row above.

## Acceptance gates

- Two validated current-source shipping cycles visit all 12 presets, wrap
  11-to-0, have zero epoch resets, and report **every peak render <59.00 ms**.
- Every current preset bucket on both boards is below 59.00 ms; anything above
  56.00 ms requires a repeat and a low-margin callout.
- No transition frame, including the 11-to-0 wrap, exceeds 59.00 ms.
- The full Phantasm image passes FLASH, ITCM, DTCM/stack, RAM2, and allocator
  floors. The profile image alone cannot satisfy this gate.
- Exact reference paths remain for every approximate geometry/noise/color
  kernel. Dense pole, seam, chamber-boundary, period-wrap, finite-value,
  framebuffer, and error-budget tests pass.
- Native CTest, ShaderBall configuration/serialization tests, release WASM
  smoke at 96x20 and 288x144, and the installed simulator pass.
- Final shipping and matched O3 reports are regenerated from attested logs, and
  all three ranked profile READMEs are updated without source-provenance drift.

## Stop conditions

Do not buy a timing pass by weakening the render telemetry, shortening a cycle
before its wrap, excluding transitions, changing the 59 ms threshold, using a
single-effect memory result, or publishing an unmatched shipping/O3 pair.
The five-preset retirement is an explicit roster decision. Further content
changes require a recorded re-authoring decision with visual review; they are
never smuggled in as an optimization.
