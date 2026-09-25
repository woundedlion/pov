# Peirce pole-cap classification A/B/A

Experimental Teensy 4.0 diagnostic, COM4, 600 MHz, real segmented POV driver. All arms use identical target, stable effect seed, fixed preset parameters, instrumentation, 70-second capture, and 32-frame windows. Runs are sequential on the same board. A is instrumentation-only commit `3156f636ff9c247370ee2bbfd861efbfa7d49aea`; B is `bf4fa1223439ab092855f10b41259e1d2cd68cdb`, adding only the classification fix and its regression test. These were isolated commits at capture time. The production correction was subsequently landed as e91e74d82 ; the diagnostic instrumentation remains separate.

The workload is **nonshipping PeirceProbe**, because Shader is deliberately simulator-only and no production roster effect instantiates Peirce. It uses the existing composed-effect lifecycle, pullback pipeline, and quadrant raster. A deliberate metadata-dependent region tint retains the classification result in the final pixel. This measures the added check in a representative consumer, not actual Shader firmware. Production feature gates remain unchanged.

## Matched results

The steady columns exclude frames 1–32 from every arm. Render = root minus canvas wait; shader is the existing `fx_shader_draw` scope. Peaks and spills cover every frame, including cold start. All captures contain 1,088 frames (34 windows).

| Config / arm | Steady shader ms | Steady render ms | All-frame render mean ms | Render peak ms | Spilled | Captured local |
|---|--:|--:|--:|--:|--:|---|
| ship / A1 | 42.929568 | 49.389869 | 49.390857 | 87.016 | 1/1088 | 2026-09-20T00:18:33 |
| ship / B | 43.077910 | 49.579500 | 49.582734 | 87.407 | 1/1088 | 2026-09-20T00:23:26 |
| ship / A2 | 42.930425 | 49.390750 | 49.391584 | 87.002 | 1/1088 | 2026-09-20T00:30:07 |
| o3 / A1 | 38.542038 | 44.771594 | 44.781631 | 79.095 | 1/1088 | 2026-09-20T00:20:59 |
| o3 / B | 38.632608 | 44.806141 | 44.815099 | 79.091 | 1/1088 | 2026-09-20T00:25:40 |
| o3 / A2 | 38.544239 | 44.775106 | 44.785234 | 78.936 | 1/1088 | 2026-09-20T00:32:08 |

## Repeat bounds and interpretation

**ship:** shader B minus mean(A1,A2) = +0.147913 ms (+0.345%); baseline repeat span 0.000857 ms. render B minus mean(A1,A2) = +0.189190 ms (+0.383%); baseline repeat span 0.000881 ms.

**o3:** shader B minus mean(A1,A2) = +0.089470 ms (+0.232%); baseline repeat span 0.002201 ms. render B minus mean(A1,A2) = +0.032791 ms (+0.073%); baseline repeat span 0.003512 ms.

Both configurations show a measurable slowdown in this diagnostic: the candidate differences exceed the observed baseline repeat spans. This does not support classifying the patch as having no measurable performance sacrifice. Both preserve the 16 fps steady display tier. These are observed repeat bounds, not confidence intervals. Only two baseline captures and one candidate were collected per config. Counter scope differences include interrupts and code-layout/cache effects. The absence of a cadence-tier change does not imply zero CPU cost.

## Correctness and generated code

A native sweep of 768 near-pole samples (384 with rounded `abs(y)==1`) finds 280 baseline edge-class mismatches versus the general projection and zero candidate mismatches. `unit_projections` passes with the new regression. The raster capture does not establish that any tiny pole-cap sample occurred on-device; the correctness sweep is separate.

Both ARM binaries retain the added `vabs.f32`, `vcmpe.f32`, `vmrs`, `it ge`, and `movge r4,#0` before storing edge_class. Shipping candidate addresses are `0x600078c6`–`0x600078f4`; global-O3 candidate `0x60008490`–`0x600084be`. Caller assembly loads edge_class and converts it for the tint, proving it affects output rather than being optimized away. See `*_actual_kernel_asm.txt`, full `*_disassembly.txt`, and `correctness-probe.txt` under `build/prof/review_20260920/analysis/115/`.

Baseline A1/A2 disassembly is identical within each config (ELF hashes differ in non-instruction metadata); saved assembly diffs contain only the input file path.

FLASH code grows **16 bytes** in each config; ITCM code and RAM allocations are unchanged. Alignment headers shrink 16 bytes, so total occupied flash is unchanged in this particular image. Global-O3 versus shipping is +15,104 B FLASH code and +11,984 B ITCM in both arms. The experimental target is not evidence that a full global-O3 roster fits.

## Validation and artifacts

The logs, JSON summaries, assembly dumps and patches named below are local, gitignored capture artifacts. They are not distributed with this repository; the tables above and linked profiles are the retained results.

All accepted captures pass `parse_profile.py validate` and have exact per-frame render telemetry, clean source provenance, correct effect/config headers, fresh initial frame numbering, and no epoch crossing. Each profile invocation also builds production Phantasm with profiling disabled. Raw logs, build/env logs, source state, ELF, map, and provenance are preserved under `C:/work/Holosphere/build/prof/review_20260920/115/`; immutable copies of full build artifacts are in its `artifacts/<capture>/` subdirectories.

The first unsupported Shader build and invalid-parameter probe boot are excluded. They produced no accepted measurements. Candidate native validation passes: 87 tests passed, one pre-existing replay-form test skipped, zero failures of 88. See candidate-native-build.log and candidate-ctest-final.log. Standard report drafts are [shipping](shipping/profile_peirceprobe_teensy_2026-09-20.md) and [global-O3](O3/profile_peirceprobe_teensy_2026-09-20.md). `metrics.json`, `summary.json`, and `comparison.json` provide machine-readable results.

Instrumentation and fix remain separate commits. `instrumentation.patch` and `finding-115-tested.patch` preserve their exact content externally. Only the isolated experimental tree was edited; canonical reports and production source remain the integrator's responsibility.

Full analysis files: `build/prof/review_20260920/analysis/115/`.
