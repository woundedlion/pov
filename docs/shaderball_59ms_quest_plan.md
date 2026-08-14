# ShaderBall 59 ms quest plan

Status: **planned; no performance implementation has begun.**

The shipping gate is a peak render time below **59.00 ms** for every authored
ShaderBall preset and every transition frame in the selective-O3 Phantasm
image. Development targets **56.00 ms** so normal orientation, interrupt, and
measurement variation cannot turn a nominal pass into a field failure. The
final evidence is a complete 17-preset on-device cycle with zero spilled
frames, not a host benchmark or a single favorable window.

This plan supersedes the completion claim for the former 23-preset bank in
`shaderball_optimization_plan.md`. That document remains the historical
optimization ledger and exact-oracle record.

## Existing evidence and baseline integrity

The 2026-08-13 shipping report is the initial symptom map. It covers all 17
presets, wraps 16 to 0, agrees with the DWT root counter within 0.2 ppm, and
reports a 147.35 ms peak with 664 of 1,472 frames spilled. Five buckets are
green and twelve are red.

The current report pair is not a matched compiler comparison. The attested
shipping log names source `f0821319`, the O3 log names `d0848a94`, and the
reports name `6f2e7e77`. All three revisions precede current master
`74120cb4`. The shipping numbers are useful for prioritization, but neither
the shipping/O3 delta nor either report's source label is acceptance evidence.
A clean current-master matched pair is the first measurement task.

| Preset | Authored topology | Peak, ms | Clean shader, ms | Gap to 59 ms |
|---:|---|---:|---:|---:|
| 16 | Gnomonic dodecahedral vector-noise + mirror grid | 147.35 | 133.79 | 88.36 |
| 0 | Glitch wave-shear grid | 146.53 | 35.50 | 87.54* |
| 13 | Sinusoidal curl lattice, fine field | 107.89 | 85.34 | 48.90 |
| 14 | Sinusoidal curl lattice, coarse field | 97.91 | 84.98 | 38.92 |
| 15 | Stereographic prism polar-wave lattice | 97.44 | 68.77 | 38.45 |
| 7 | Dodecahedral grid + mirror | 82.27 | 69.73 | 23.28 |
| 8 | Peirce dodecahedral grid | 81.58 | 67.64 | 22.59 |
| 9 | Dodecahedral noise grid | 79.59 | 51.93 | 20.60* |
| 11 | Gnomonic dodecahedral wave/mirror grid | 79.08 | 66.34 | 20.09 |
| 12 | Gnomonic affine lattice contour | 78.86 | 43.91 | 19.87* |
| 10 | Dodecahedral noise lattice + mirror | 78.25 | 56.80 | 19.26 |
| 6 | Kaleidoscope noise grid + edge fade | 75.36 | 62.68 | 16.37 |
| 5 | Peirce kaleidoscope lattice | 49.35 | 30.21 | pass |
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

## Measurement campaign

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

### Preset-family attack order

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

## Acceptance gates

- Two validated current-source shipping cycles visit all 17 presets, wrap
  16-to-0, have zero epoch resets, and report **every peak render <59.00 ms**.
- Every fixed run for presets 0 and 6-16 is below 59.00 ms; anything above
  56.00 ms is repeated and called out as low margin.
- No transition frame, including the 16-to-0 wrap, exceeds 59.00 ms.
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
Content changes are allowed only as an explicit last-resort re-authoring
decision with visual review; they are never smuggled in as an optimization.
