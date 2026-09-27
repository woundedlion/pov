# Spherical perspective implementation and admission (2026-09-26)

The [spherical perspective specification](../specs/spherical_perspective_spec.md)
is implemented. HyperLattice retains its two authored analytic configurations;
the additional patterns are core demonstrators with explicit experimental
status. Final firmware and device measurements use clean source
`8441a47efd4570cc09914b1b43c03c758606f058`. The subsequent standalone-profiler
header move and documentation commits do not change either shipping image.

## Architecture and compatibility

The shared [ray core](../../core/render/ray.h) owns camera/domain preparation,
bounded analytic event merging, certified first-boundary search, appearance,
and compositing. Pattern definitions and adapters live under
[SDF](../../core/render/sdf.h). A reusable
[pullback ray stage](../../core/render/pullback/ray.h) binds the renderer to
sphere samples. Geometry and numerical tests compile without effects.

The existing volume renderer uses the extracted scalar kernel without changing
its closest-surface policy. Spherical torus and warped-torus demonstrations
reuse the same shape types and kernel. Certified surface search is a separate
policy: proximity does not certify a hit, verification queries consume budget,
and unresolved, invalid, unsupported, and exhausted searches remain explicit.
Verified subray averaging preserves failed-search status and premultiplied
coverage. A near-slice 4D ball does not create a phantom certified intersection.

HyperLattice owns choreography, controls, palette selection, and admission.
Its public effect ID and `cubic-flight` / `hypercube-flight` preset IDs survive.
Dimensional rift is removed. Schema version 10 rejects version-9 snapshots
without mutation. Configuration changes adopt one complete geometry block;
cross-configuration transitions switch at their midpoint, while common near
fade interpolates. The selector is dense, Lattice Planes replaces Shells, and
4D Spin becomes read-only in 3D with schema-generation notification.

## Admission

| Demonstrator | Decision | Evidence or limit |
| --- | --- | --- |
| Cubic / 3D analytic coverage | Shipping, authored preset and transition sweep | Fixed and cycling device captures; compatibility regression images |
| Cubic / 4D slice analytic coverage | Shipping, authored preset and transition sweep | Fixed and cycling device captures; compatibility regression images |
| Triangular-prism framework | Experimental; excluded from selector | Native nonorthogonal-stream, ordering, and rendered-contribution tests; no device admission capture |
| Torus and warped torus with spherical rays | Experimental; excluded from selector | Native shared-kernel and placement comparisons; no device admission capture |
| Marched 3D / 4D lattice field | Experimental; excluded from selector | Native query, slice, and geometry comparisons; no device admission capture |
| Cosine and gyroid nodal surfaces | Experimental; excluded from selector | [Native quality sweep](spherical_periodic_native_2026-09-26.md) and [device cost sweep](spherical_periodic_device_2026-09-26.md) fail quality/deadline admission |

The shipping policy intentionally preserves bounded plane-crossing coverage
and horizon fading. These are declared approximations, not exact ray/cylinder
intersections. At most three planes per axis produce nine candidates in 3D or
twelve in 4D; the authored two-plane presets produce at most six or eight.
The generic trace cap is twelve candidates and 32 layers, above these finite
adapter bounds. Each ray buffers one merge identity; all lattice candidates
share that identity. No surface marcher or hidden query budget runs here.
Reference tests cover normal completion, opacity saturation, exhaustion,
invalid payloads, partial contribution preservation, and strict cursor progress.

Device admission covers the following authored values and both 240-frame
transitions, with 320-frame holds. Motion and palette cycling remain enabled.
Manual control extrema have native validation but are not a measured exhaustive
device performance envelope.

| Parameter | Cubic | 4D slice |
| --- | ---: | ---: |
| Radial start / cell size | 1 / 1 | 0 / 1 |
| Wire radius / softness, in cells | 0.055 / 0.08 | 0.03546 / 0.029612 |
| Far distance | 4.198 | 8 |
| Near fade / AA strength | 0.5 / 1 | 0.5 / 1 |
| Flight speed | 0.05 | 0.03 |
| 3D / 4D spin | 0.015 / 0 | 0.01089 / 0.015 |
| Planes / appearance | 2 / depth | 2 / depth |

The [ARM audit](evidence/spherical_codegen_2026-09-26/arm-audit.txt) measures a
964-byte effect object plus 4,432 bytes of persistent palette payload from an
aligned arena start (4,452-byte declared bound). The parameter block is 48
bytes. FrameState / PreparedTrace are 100 / 156 bytes; per-ray Events is 104
bytes and its single buffered Contribution is 40 bytes. These object sizes
are not additive stack requirements: the optimizer scalarizes or removes them.
There is no per-ray allocation. The specialized shader's local reservation
falls from 340 to 176 bytes; its own frame including saved registers is 248
bytes, excluding callees and interrupts. Global RAM1 variables and RAM2 usage
are unchanged from baseline.

## Code generation and placement

ARM ELF symbols and disassembly, rather than source-level inline guesses,
guided these changes:

- `HS_COLD_MEMBER` keeps runtime preset lookup and schema refresh in flash,
  prevents inline/IPA duplication, and retains compile-time preset evaluation.
- `HS_FLASH_MEMBER` places the measured frame entry and animation update in
  cached flash. Hot per-ray helpers remain integrated into the shader loop.
- Two- and three-plane specialized aliases share one compiled pipeline.
- The lattice adapter declares a single merge-identity slot instead of paying
  for four general-purpose contributions on every ray.
- Unused cursor fields are initialized only when a stream becomes active.
  This removes the remaining per-ray `memset` without reading inactive fields.
- Static depth shading removes runtime feature-palette selection and integer
  feature division from the specialized path.
- One terminal flush path replaces repeated compositing epilogues. Invalid
  query status still takes precedence when flushing saturates the consumer.
- Preparation uses references and initializes/scales the result matrix directly,
  removing a 64-byte intermediate. Measured preparation stack including saved
  registers drops from 352 to 264 bytes. Schema lookup occurs on mode changes.

The capacity change alone reduced the two shader bodies from 6,666 / 7,112
bytes to 4,476 / 4,748 bytes and reduced saved floating-point registers.
Subsequent shared terminal code reduced palette-emission sites from five to
two. Final specialized / generic shader bodies are 3,936 / 4,440 bytes in
flash, each with zero `memset` or `memcpy` calls. Small one-caller preparation
leaves were left inline where extracting
them would add calls without removing duplicate machine code.

| Phantasm checkpoint | RAM1 code | Distance below 196,608-byte ceiling |
| --- | ---: | ---: |
| Original `79a93c6f8` | 195,304 | 1,304 |
| Initial extracted implementation | 201,480 | -4,872 |
| First fitting extraction `b9078314a` | 196,568 | 40 |
| Final measured source `8441a47ef` | 195,544 | 1,064 |

The final pass saves 1,024 bytes of ITCM versus the first fitting extraction,
and 5,936 versus the initial extraction. The gate reserves **zero** ITCM padding
(`min_headroom_bytes: 0`); its ceiling follows the existing FlexRAM bank layout
and stack floor. No budget was raised to fit the implementation.

## Firmware and runtime validation

All seven environments build cold: holosphere, holosphere_dma, phantasm,
profile, profile_o3, bench, and phantasm8. Both configured size/layout gates
pass. The warning audit covers all 28 first-party translation units with zero
warnings. [Metrics and ELF hashes](evidence/spherical_codegen_2026-09-26/firmware-summary.json),
[cold build](evidence/spherical_codegen_2026-09-26/firmware-cold.txt), and
[warning audit](evidence/spherical_codegen_2026-09-26/firmware-warnings.txt)
retain the evidence.

| Full Phantasm region | Baseline | Final | Change |
| --- | ---: | ---: | ---: |
| RAM1 code | 195,304 | 195,544 | +240 |
| RAM1 variables | 314,784 | 314,784 | 0 |
| FLASH code | 501,392 | 502,072 | +680 |
| FLASH data | 725,460 | 725,448 | -12 |
| RAM2 variables | 520,064 | 520,064 | 0 |

The bank-rounded RAM1 allocation still leaves 12,896 bytes for stack/locals.
Experimental pattern definitions do not add emitted geometry to Phantasm.

Native validation finishes with 100 passing tests and one intentional replay
skip. The [full run](evidence/spherical_codegen_2026-09-26/native-full.txt)
initially found the standalone profiling header in the unit-test directory;
moving it to tools resolves the [include-policy rerun](evidence/spherical_codegen_2026-09-26/native-roster.txt).
Ray contracts, event merging, pattern definitions, demonstrators, HyperLattice
compatibility/control tests, volume regressions, memory/stack checks, and the
120-frame effect smoke suite pass. The release WASM build and
[full WASM smoke](evidence/spherical_codegen_2026-09-26/wasm-smoke.txt) also pass.

Final captures use the real 288-by-144 four-segment driver at 600 MHz, live
flywheel/DMA interrupts, and a 62.5 ms display deadline. The following figures
include every captured individual runtime frame, including trailing partial
windows; frame 1 is setup and is separately excluded.

| Capture | Live frames | Mean render ms | Peak render ms | Spills |
| --- | ---: | ---: | ---: | ---: |
| Shipping fixed cubic, COM3 | 1,096 | 42.826 | 47.077 | 0 |
| Shipping fixed 4D slice, COM4 | 1,096 | 49.725 | 57.797 | 0 |
| Shipping full cycle, COM4 | 1,896 | 45.975 | 55.579 | 0 |
| Global-O3 full cycle, COM3 | 1,576 | 44.883 | 54.879 | 0 |

Setup renders are 77.287, 67.786, 77.272, and 76.882 ms respectively. Both
cycling captures visit both presets and wrap, with monotone frame numbering.
The [shipping report](shipping/profile_hyperlattice_teensy_2026-09-26.md) and
[O3 report](O3/profile_hyperlattice_teensy_2026-09-26.md) provide counter trees,
complete-window statistics, ISR accounting, and raw/provenance links.

The pre-optimization extracted 4D image averaged 70.988 ms and spilled on
531 of 567 runtime frames. The final fixed capture averages 49.725 ms with
zero spills, approximately 30% faster. It remains slower than the original
fixed-preset baseline (44.200 ms); cubic likewise rises from 37.669 to
42.826 ms. The shared architecture fits and meets the measured deadline, but
does not claim zero abstraction cost or performance equivalence for marchers.
