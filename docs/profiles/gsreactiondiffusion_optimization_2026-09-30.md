# GSReactionDiffusion optimization campaign — 2026-09-30

The original optimization campaign retained four SSAA samples, seven stencil nodes, six physics substeps, the original ordered graph, and the existing palettes. This report records twelve main candidate captures, a matched packed-coordinate trial, diagnostic instrumentation, and final shipping/global-O3 captures.

Current shading, measured at `7baf3cc430753e1a3693134150b282f212aa4e9e` on COM4, has shipping peak **313.695 ms** and global-O3 peak **279.686 ms**. Shipping spills **443/443 live frames**, and O3 spills **538/538**, after excluding setup frame 1. The earlier 36.528/36.477 ms results predate pigment blending, hue rotation and shimmer and do not characterize the current artwork. [Current-state measurements](#current-shading-measurements) include raw evidence and current image sizes.

The original optimization campaign finished with shipping peak **36.528 ms**, down from **55.413 ms**: **18.885 ms (34.08%) faster**. The independent selected-candidate pass measured 36.549 ms; campaign global-O3 measured **36.477 ms**. Both campaign captures have **0/2047 live-frame spills**. The controlled optimization series adds **288 B of full-shipping ITCM**, **9,760 B of FLASH data**, and **0 B of RAM variables**. That campaign's rebased Phantasm image has 166,680 B RAM1 code, **336 B above the original baseline**; 48 B of that increase appeared when incorporating concurrent master changes. All 15 final shader probes stay within **2/65535 per channel** with **zero coverage changes**. Final native physics hashes match all **2,049 baseline states and 33 parameter cases**.

Each main capture used the supported `tools/profile_one.sh` wrapper and locked Teensy on **COM3**, 600 MHz, four real POV segments, 480 RPM, 288x144 effect resolution, `ship` configuration, **130 seconds**, **32-frame windows**, and `HS_PROFILE_EPOCH_REVS=1200`. All twelve main captures contain 64 complete windows, frames 1-2048, without an epoch reset. Statistics use **2,047 live frames, 2-2048**, excluding frame 1. Peak is the maximum per-frame render `r` column, not the wall time that includes display-buffer waiting. Every main capture validates, and every one has **zero frames exceeding 62.5 ms** (`spilled=0`); this is a deadline statistic, not a claim about CPU register spills.

The captured compiler is Arm GNU Toolchain 15.2.Rel1, GCC 15.2.1. Source revisions, dirty diffs, compiler identity, and SHA-256 hashes of profile/shipping ELFs and environment dumps are recorded with each capture. Baseline source revision is `c1837c0973a5276e3af34a7ec9e30c5d1ba3b568`. These are single captures per candidate; sub-millisecond differences should not be represented as confidence intervals or guaranteed worst-case bounds.

Negative timing deltas mean faster. Parent identifies the actual retained starting point: graph physics starts from symmetry after branchless was rejected. ITCM costs below are full shipping `phantasm` RAM1 code, not the much smaller single-effect profile image.

| Experiment | Parent | Peak ms | Delta ms parent | Delta ms baseline | Shipping RAM1 code B | Delta ITCM B parent / baseline | Decision |
|---|---|---:|---:|---:|---:|---:|---|
| [baseline](evidence/gs_optimization_2026-09-30/baseline/summary.json) | -- | 55.413 | +0.000 | +0.000 | 166,344 | +0 / +0 | Reference |
| [interchange](evidence/gs_optimization_2026-09-30/interchange/summary.json) | baseline | 52.778 | -2.635 | -2.635 | 166,312 | -32 / -32 | Keep, then extend |
| [symmetry](evidence/gs_optimization_2026-09-30/symmetry/summary.json) | interchange | 50.574 | -2.204 | -4.839 | 166,344 | +32 / +0 | Keep |
| [branchless](evidence/gs_optimization_2026-09-30/branchless/summary.json) | symmetry | 52.277 | +1.703 | -3.136 | 166,344 | +0 / +0 | Reject: slower |
| [graphphysics](evidence/gs_optimization_2026-09-30/graphphysics/summary.json) | symmetry | 45.662 | -4.912 | -9.751 | 166,408 | +64 / +64 | Superseded by ordered physics |
| [graphcull](evidence/gs_optimization_2026-09-30/graphcull/summary.json) | graphphysics | 42.954 | -2.708 | -12.459 | 166,568 | +160 / +224 | Keep cull; replace inherited physics |
| [ordered](evidence/gs_optimization_2026-09-30/ordered/summary.json) | graphcull | 42.659 | -0.295 | -12.754 | 166,728 | +160 / +384 | Keep constrained physics |
| [graphshader](evidence/gs_optimization_2026-09-30/graphshader/summary.json) | ordered | 39.425 | -3.234 | -15.988 | 166,744 | +16 / +400 | Keep |
| [scaled](evidence/gs_optimization_2026-09-30/scaled/summary.json) | graphshader | 37.517 | -1.908 | -17.896 | 166,664 | -80 / +320 | Keep |
| [floatpalette](evidence/gs_optimization_2026-09-30/floatpalette/summary.json) | scaled | 36.819 | -0.698 | -18.594 | 167,336 | +672 / +992 | Keep palette arithmetic; add loop control |
| [roundfold](evidence/gs_optimization_2026-09-30/roundfold/summary.json) | floatpalette | 37.263 | +0.444 | -18.150 | 167,448 | +112 / +1,104 | Reject: slower; closer rounding behavior |
| [paletteloop](evidence/gs_optimization_2026-09-30/paletteloop/summary.json) | floatpalette | 36.549 | -0.270 | -18.864 | 166,632 | -704 / +288 | Selected; confirmation pending |

| Experiment | Profile RAM1 code B | Delta profile B parent / baseline | Shipping FLASH data B | Delta FLASH data B parent / baseline | Profile FLASH data B |
|---|---:|---:|---:|---:|---:|
| baseline | 23,784 | +0 / +0 | 736,472 | +0 / +0 | 242,108 |
| interchange | 23,640 | -144 / -144 | 736,472 | +0 / +0 | 242,108 |
| symmetry | 23,672 | +32 / -112 | 736,472 | +0 / +0 | 242,108 |
| branchless | 23,672 | +0 / -112 | 736,472 | +0 / +0 | 242,108 |
| graphphysics | 23,736 | +64 / -48 | 738,552 | +2,080 / +2,080 | 244,188 |
| graphcull | 23,896 | +160 / +112 | 738,552 | +0 / +2,080 | 244,188 |
| ordered | 23,896 | +0 / +112 | 738,552 | +0 / +2,080 | 244,188 |
| graphshader | 23,928 | +32 / +144 | 746,232 | +7,680 / +9,760 | 251,868 |
| scaled | 23,832 | -96 / +48 | 746,232 | +0 / +9,760 | 251,868 |
| floatpalette | 24,504 | +672 / +720 | 746,232 | +0 / +9,760 | 251,868 |
| roundfold | 24,632 | +128 / +848 | 746,232 | +0 / +9,760 | 251,868 |
| paletteloop | 23,800 | -704 / +16 | 746,232 | +0 / +9,760 | 251,868 |

All twelve main captures keep full-shipping RAM1 variables at **314,784 B**, profile RAM1 variables at **315,168 B**, and RAM2 variables at **520,064 B**. Shipping free local-variable space remains **12,896 B** because these code changes stay within the same ITCM bank allocation; code padding changes instead. Profile free local-variable space remains 176,352 B. The graph run representation contains 148 ordered runs x14 B plus a 4 B count =2,076 B logical payload, linked as a 2,080 B FLASH-data increment. Adding the 7,680 B node-to-run map brings logical compact tables to 9,756 B and linked total data growth to 9,760 B. The original 92,160 B adjacency table remains for other consumers. GS shader refinement/stencil reads use the compact decoder, while BZ retains its original default access path.

Interchange makes the seven stencil nodes the outer loop and updates the four SSAA accumulators per node. It avoids repeatedly loading gathered positions/concentrations for each sample. The ARM shader prologue's local stack reservation drops from 136 B to 52 B; saved register sets also change, so those figures are not total stack usage. Peak improves 2.635 ms and full-shipping ITCM falls 32 B.

Symmetry represents each horizontal SSAA pair as midpoint plus/minus an offset. Distances share a squared base and tangent cross term, reducing repeated geometry work. It saves another 2.204 ms for 32 B ITCM relative to interchange, returning total shipping ITCM to baseline. This algebra changes floating-point rounding slightly; image evidence below bounds the tested effect.

Branchless replaces the `u > 0` support test with a clamp and unconditional weight accumulation. It regresses 1.703 ms with no measured ITCM change. Keep the support branch: skipping work outside the kernel is valuable on this workload.

Graph physics sweeps 148 contiguous ranges of identical ordered neighbor offsets instead of streaming 92 KB of absolute indices on each of six substeps. It saves 4.912 ms versus symmetry and costs 64 B shipping ITCM plus the compact table. This first version allowed different floating-point reassociation and is not the retained physics implementation. Graph cull applies the same exact offsets to both directed dilation passes, saving a further 2.708 ms for 160 B shipping ITCM. Its integer threshold/boolean result is exact; the capture still inherited the unconstrained physics and is superseded as a combined candidate.

Ordered keeps the compact run sweep but matches the baseline Cortex-M7 Laplacian evaluation: add neighbors 0 and 1; fuse subtraction of six times the center into that sum; then add neighbors 2, 3, 4, and 5 in order. Empty inline-assembly register constraints prevent reassociation without emitting arithmetic instructions. Actual baseline and ordered ARM kernel disassemblies were compared through diffusion, reaction, timestep, and clamp operations. Ordered improves another 0.295 ms in this capture while costing 160 B more full-shipping ITCM than graphcull; profile RAM1 code happens to be unchanged.

Graphshader adds one uint8 run ID per node so shader refinement and stencil gathering decode `node + neighbor_runs[neighbor_run_index[node]].delta[k]`. This extends the exact sequential representation to random accesses, reducing the active adjacency footprint to about 9.8 KB. It saves 3.234 ms for 16 B shipping ITCM and 7,680 B FLASH data. Ordered tuples and node indices remain unchanged.

Scaled folds the inverse support radius into the shared horizontal-pair calculation: compute `base_u = 1 - base * INV_R2` once per row, scale the tangent cross term once, and obtain the two support coordinates as `base_u - cross` and `base_u + cross`. It saves 1.908 ms while reducing shipping ITCM by 80 B. This retains all four SSAA samples and all seven stencil nodes; it changes evaluation/rounding order rather than reducing samples.

Floatpalette interpolates the existing palette as scaled floating RGB for every contributing SSAA sample, accumulates the quarter-scale colors, then rounds/clamps once at the pixel boundary. It avoids the per-sample Q16 interpolation/quantization and integer averaging path without averaging palette coordinates or dropping samples. It saves 0.698 ms against scaled at a cost of 672 B shipping ITCM. The tradeoff is tiny rounding differences on many lit pixels, including a constant-field bias of up to 2 channel units, not a field evolution change.

Roundfold retains floating palette interpolation but restores an integer contribution for each sample using `uint32_t(rgb + 0.625f)` for the quarter-scale channels. The 0.625 term folds round-to-nearest followed by `(channel + 2) / 4`; exact preservation still depends on the interpolation arithmetic. It reduces the observed average error substantially, but loses 0.444 ms versus floatpalette and costs another 112 B shipping ITCM. These data support choosing based on the explicit rounding tradeoff rather than claiming either is strictly better in every metric.

Paletteloop starts from floatpalette, discards the roundfold alternative, and applies `#pragma GCC unroll 1` to the four-sample palette loop. This keeps the compact loop in the GCC device build instead of duplicating the interpolation/accumulation body. The measured result is another 0.270 ms peak improvement and **704 B less shipping ITCM** relative to floatpalette, with no data or variable allocation change. The selected ARM shader is **1,664 B** and reserves **120 B total stack**: 28 B saved general registers, 56 B saved floating registers, and 36 B local storage. These figures come from `paletteloop-shader-disassembly.txt`; the 36 B `sub sp` alone is not total stack use. The source guard leaves Clang host builds unaffected.

Scope values below are averages over each capture's own peak-containing 32-frame window, not individual-frame scope maxima. The peak moves from early frame 489 to frame 1450 after graphshader, so rows with different windows cannot isolate a phase speedup by direct subtraction.

| Experiment | Window | Simulate ms/frame | Cull ms/frame | Shader draw ms/frame | Peak frame |
|---|---|---:|---:|---:|---:|
| baseline | 481-512 | 11.975 | 3.272 | 37.451 | 489 |
| interchange | 481-512 | 11.969 | 3.290 | 34.491 | 489 |
| symmetry | 481-512 | 11.977 | 3.278 | 32.445 | 489 |
| branchless | 481-512 | 11.971 | 3.283 | 34.303 | 491 |
| graphphysics | 481-512 | 6.996 | 3.282 | 32.524 | 493 |
| graphcull | 481-512 | 6.996 | 0.471 | 32.465 | 489 |
| ordered | 481-512 | 6.713 | 0.473 | 32.485 | 489 |
| graphshader | 1441-1472 | 6.714 | 0.344 | 30.812 | 1450 |
| scaled | 1441-1472 | 6.714 | 0.343 | 28.877 | 1450 |
| floatpalette | 1441-1472 | 6.714 | 0.344 | 28.270 | 1450 |
| roundfold | 1441-1472 | 6.714 | 0.343 | 28.674 | 1450 |
| paletteloop | 449-480 | 6.716 | 0.322 | 28.099 | 464 |

Quality evidence comes from the standalone host harness, not captured device images. `quality.cpp` compares the candidate shader against the original four-sample shader using the same state, orientation, and palette over all 288x144 pixels. Each candidate has 15 probes: 11 evolved frames through frame 2000 and four synthetic fields (constant, gradient, hard hemisphere, and oscillating field). Baseline and interchange are pixel-identical on all probes. Symmetry has maximum channel error **2/65535**, zero lit/unlit coverage changes, maximum per-probe MAE **0.004854681** and RMSE **0.069905874** in 16-bit channel units, and minimum global SSIM **0.999999999962**. Branchless also stays within 2/65535 on these probes but is rejected on speed. Global SSIM here is the harness's whole-image/channel statistic, not a local-window perceptual certification.

Physics evidence is separate from shader evidence. Comparing baseline and unconstrained graphphysics standalone-host CSVs finds changed field/lifecycle hashes in **2,048 of 2,049** recorded states (first change at frame 1) and changed float hashes in **25 of 33** parameter corner cases. The ordered host implementation matches all baseline field/lifecycle hashes for **2,049 states** and all float hashes for **33 corners**. These are host standalone-harness results. The device-side assurance is inspection of the actual ARM kernel's arithmetic sequence; no device full-state hash comparison was captured, so this does not establish an on-device bit-identical trajectory.

The compact graph expands exactly to all **46,080 ordered neighbor indices** and covers every node without gaps. Native regression additionally compares both cull arrays to directed original-table gathers for **36 full-node cases**, including sparse/dense fields, constant/threshold boundaries, and isolated north/south pole activity. Sentinel destinations and guard bytes detect missing output nodes and writes outside the arrays. The cull test passes against both the original implementation and the final candidate.

The packed-coordinate experiment is a separate COM4 branch from graphshader, not from floatpalette. It has its own same-board baseline and standard 130-second captures with deep scopes disabled. Both validate with 2,047 live frames and no deadline spills. The main COM3 graphshader result is 39.425 ms; the COM4 repeat is 39.452 ms. The experimental cost is measured against that COM4 baseline, not against a later main candidate.

| COM4 experiment | Peak ms | Delta peak vs own baseline | Shipping ITCM B | Delta ITCM own / original baseline | Profile RAM1 code B | Shipping FLASH data B |
|---|---:|---:|---:|---:|---:|---:|
| [packed-baseline](evidence/gs_optimization_2026-09-30/packed-baseline/summary.json) | 39.452 | +0.000 | 166,744 | +0 / +400 | 23,928 | 746,232 |
| [packed-candidate](evidence/gs_optimization_2026-09-30/packed-candidate/summary.json) | 41.889 | +2.437 | 167,416 | +672 / +1,072 | 24,504 | 746,232 |

Packed stores oriented XYZ as signed Q15 with B in an 8-byte node, reducing world scratch from 92,160 B to 61,440 B and raster scratch to 76,800 B. The physics scratch maximum remains 122,880 B, so allocated RAM does not decrease. Saturating `QSUB16` differences and `SMUAD` dual multiply-adds are present in actual ARM shader code. Saturation occurs only for distances already beyond kernel support; the largest three squared saturated lanes sum to 3,221,225,472 and fit uint32. The implementation converts each relevant result to unsigned before summing to avoid signed overflow.

Despite those DSP instructions, packing/conversion and shader overhead lose **2.437 ms (6.18%)** peak and **2.387105 ms** mean relative to the COM4 baseline. At the matched frames 449-480, shader time grows 30.94359 to 32.82119 ms/frame and orientation/packing grows 0.37138 to 1.01325 ms/frame. Shader code grows 1,784 to 2,316 B; the archived notes report stack usage falling 136 to 112 B. Orientation method code grows 188 to 584 B. The prototype is rejected on performance and is not committed.

Packed host image evidence has maximum error 167/65535, maximum per-probe MAE 2.124115869 and RMSE 4.451995538, no pixels above 4,096 error, and at most two lit/unlit differences per 41,472-pixel probe. The smallest global RGB SSIM is 0.999999844405. Coordinate quantization bounds are retained in `packed-candidate/quantization-bounds.json`; they bound geometric distance error, not final color error across a hard coverage threshold. ARM uses fused quantization arithmetic while a host may leave multiply/add separate, so host images are not device bit-exact certification. Selected sRGB8 comparisons and fixed-scale error heatmaps are retained in portable evidence; larger RGB16/NPY files remain in the local artifact archive.

The deep experiment is diagnostic only: **COM4**, 600 MHz, **50 seconds**, 32-frame windows, **24 windows / 767 live frames**. Per-pixel RAII scopes perturb both measurement overhead and register allocation, so its 39.676 ms diagnostic maximum is not a comparable optimization row. At frames 481-512, refinement costs 2.457662 ms/frame (139.31 cycles/call), weights 17.025159 ms/frame (965.07 cycles/call), and palette 6.058112 ms/frame (343.40 cycles/call), with about 10,584.84 calls/frame. At its own peak window 449-480 those phases are 2.559856, 17.229380, and 6.914233 ms/frame. This identifies weight evaluation and palette work as useful targets; it does not establish their exact uninstrumented fraction. See `deep/summary.json` and the scoped source diff.

Host shader error statistics for the added experiments use the same 15-probe methodology; rows report maxima across probes and minimum global RGB SSIM. Errors are in original 16-bit linear channel units. These are empirical bounds for sampled cases, not a proof over every possible palette, state, and orientation.

| Variant | Max channel | Max MAE | Max RMSE | Max coverage changes | Min global RGB SSIM |
|---|---:|---:|---:|---:|---:|
| graphshader | 2 | 0.004854681 | 0.069905874 | 0 | 0.999999999962 |
| scaled | 2 | 0.004669817 | 0.068336062 | 0 | 0.999999999960 |
| floatpalette | 2 | 1.333333333 | 1.414213562 | 0 | 0.999999491041 |
| roundfold | 2 | 0.004951132 | 0.070478414 | 0 | 0.999999999960 |
| selected final host quality | 2 | 1.333333333 | 1.414213562 | 0 | 0.999999491041 |
| packed | 167 | 2.124115869 | 4.451995538 | 2 | 0.999999844405 |
| palette-average | 10950 | 252.092488104 | 745.381686507 | 0 | 0.998831476226 |
| palette-guard002 | 325 | 1.333333333 | 7.920785052 | 0 | 0.999999491041 |

The palette-average host-only experiment replaces four color lookups with one lookup at the average contributing palette coordinate. Palette nonlinearity makes the operations unequal: `palette(mean(t))` need not equal `mean(palette(t))`. Its worst channel error is 10,950/65535, maximum per-probe MAE 252.092488104, and up to 1,231 pixels per probe exceed 4,096 error. A guard requiring a palette-coordinate span at most 0.002 reduces errors but still allows 325 channel units. Both fail the quality standard set by the retained <=2-unit candidates. No device timing or ITCM measurement for either average-palette variant is present, so neither receives an invented performance row.

## Original optimization validation and provenance

The final effects and reaction-graph checks pass. The rebased full native suite ran 103 tests: **101 passed, one skipped, one failed**. The failure is the pre-existing guard-coverage bookkeeping for newly added guards in `targets/wasm/engine_bindings.h` (14 versus 13 allowed) and `workbench/shader/shader_host.h` (9 versus 8 allowed), outside this optimization. Its assertions remain enabled. The earlier generated `Profile.ino.cpp` census race was removed by rebuilding the generated census after PlatformIO finished; it is not the remaining failure. A concurrently introduced ambiguous `FamilyRank` lookup prevented compilation and was fixed with explicit namespace qualification in a separate commit.

Both `phantasm` and `holosphere` pass the firmware region-budget and layout gates. Final Phantasm RAM1 variables are 314,784 B, RAM1 code 166,680 B, FLASH data 746,232 B, and RAM2 variables 520,064 B. Free RAM1 local-variable space remains 12,896 B. The final host fast-math path recognizes Windows `_M_FP_FAST`, which remains defined when Clang's `-fno-finite-math-only` removes `__FAST_MATH__`; this restores the baseline host evaluation sequence without changing ARM code. The standalone state comparison was rerun after that correction and matches all baseline hashes.

The spatial-error population test excludes at most two RGB16 units of rounding while preserving raw changed-pixel counts, coverage, hard-error, maximum-error, and mean-error checks. A separate fixed-stencil reference test isolates numerical shader errors from the existing shared-stencil approximation. Shader quality remains an empirical host comparison, not captured LED imagery.

The final shipping capture used source `0fa3bede8`; the O3 capture includes the later namespace and host-only fast-math fixes. Captures contain exact source revision, compiler identity, build flags, working diff, and ELF hashes. Build artifacts used GCC 15.2.1. Full ELF/map files remain in the local `C:/work/hs-gs-support` archive; portable raw captures, build summaries, source patches, focused disassembly, quality/state CSVs, and reproduction harnesses are checked in alongside this report. Harness source is stored with a `.txt` suffix to keep it out of production source scans.

The full Phantasm `.text.itcm` bytes are identical before and after those two follow-up fixes: SHA-256 `2a76a4735763b5fa2769588e20287ef94270025f4ddbc78196d720a9c76f1507`. Thus neither changes the shipping ARM instruction stream measured in the final capture.

Portable text copies normalize line endings and trailing whitespace. Source diffs are stored verbatim in each capture's source JSON file, under the `patch` field; decoding that string restores the patch. [The portable manifest](evidence/gs_optimization_2026-09-30/portable_manifest.json) records input and retained hashes. ELF and environment hashes in capture footers identify the original local artifacts.

The original campaign captures remain in `evidence/gs_optimization_2026-09-30/finalship/` and `evidence/gs_optimization_2026-09-30/finalo3/`. Their global-O3 image added 13,856 B FLASH code and 8,896 B ITCM. The linked profile reports now describe the current shading captures below.

The measured bottleneck remains four-sample shading. Branchless support evaluation and packed Q15/DSP geometry were slower, while averaging palette coordinates produced excessive color error. No tested larger approximation improved both speed and the retained image-quality bounds.

[Final image metrics](evidence/gs_optimization_2026-09-30/quality-final.csv), [physics states](evidence/gs_optimization_2026-09-30/state-quality-final-states.csv), and [parameter corners](evidence/gs_optimization_2026-09-30/state-quality-final-corners.csv).

## Current shading measurements

Both current images use clean source `7baf3cc430753e1a3693134150b282f212aa4e9e`: pigment diffusion and two-palette accumulation, hue/shimmer shading and the per-frame noise LUT are enabled. This candidate also includes the noise-modifier default-speed and pigment scratch-lifetime corrections. The supported device-lock wrapper captured both images on COM4 for 130 seconds with 32-frame windows and `HS_PROFILE_EPOCH_REVS=1200`. No epoch reset occurred. A failed shipping capture attempt on COM3 and the wrapper's stale-image retry are not used as timing evidence.

| Configuration | Captured | Live frames | Mean render ms | Peak render ms | Spilled | Scope/ISR frames |
|---|---|---|--:|--:|--:|---|
| Shipping | 2026-09-30 19:21 | 2–444 | 262.576 | 313.695 | 443/443 (100.00%) | 33–416 |
| Global O3 | 2026-09-30 19:15 | 2–539 | 212.279 | 279.686 | 538/538 (100.00%) | 33–512 |

Every complete individual `f` row is retained in runtime statistics, including the unfinished final window; only setup frame 1 is excluded. Counter and ISR summaries instead use complete windows after the first. The measured shipping peak is 251.195 ms over the 62.5 ms deadline, requiring 5.02× render speedup to fit it. Both compiler configurations miss the deadline throughout their captured live intervals. These captures cover fewer simulation frames than the earlier 16 fps campaign, so they do not establish a full-lifecycle worst case.

Current full-roster Phantasm passes its firmware size/layout gate: RAM1 code **171,928 B**, RAM1 variables **314,784 B**, FLASH data **746,852 B**, RAM2 variables **520,064 B**, and local-variable space **12,896 B**. Relative to the campaign's final rebased image, these are +5,248 B RAM1 code, +0 B RAM1 variables and +620 B FLASH data. They are whole-revision deltas, not isolated costs of one shading feature. Current single-effect global-O3 versus shipping image deltas are **+16,024 B FLASH code** and **+10,832 B ITCM**.

The historical campaign's host-image error bounds and state hashes describe its earlier shaders and physics. No equivalent current-artwork image comparison is claimed here. Both current timing captures use COM4; the old campaign used COM3, a further reason not to interpret the old/new timing difference as a controlled optimization experiment.

Current reports: [shipping](shipping/profile_gsreactiondiffusion_teensy_2026-09-30.md) and [global O3](O3/profile_gsreactiondiffusion_teensy_2026-09-30.md). Portable evidence includes [shipping summary](evidence/gs_shading_2026-09-30/ship/summary.json), [O3 summary](evidence/gs_shading_2026-09-30/o3/summary.json), their complete raw captures, parser validation, compiler/build/environment records and source state. Original ELF/map files remain in the local artifact directories named by each provenance file; manifests record original and retained hashes.
