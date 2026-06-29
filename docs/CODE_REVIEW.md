# Holosphere — Code Quality Review

**Date:** 2026-06-29
**Scope:** The Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`, `targets/`, `tests/`, build/CI tooling) and the daydream web simulator (`daydream/*.js`, `tools/`). Out of scope per request: `core/effects_legacy.h`, `core/rotate.h`, `targets/Holosphere/Holosphere.ino`, and all vendored third-party code (`three.js/`, `node_modules/`, `FastNoiseLite.h`).
**Method:** Thirteen component sub-agents read every in-scope file line-by-line against the architecture described in the README. Every substantive finding was then handed to an independent verifier that re-read the cited code and returned a CONFIRMED / REFUTED / ALREADY-HANDLED verdict with a corrected severity. Findings that did not survive verification were dropped; the dropped claims are noted at the end for transparency.

---

## Executive Summary

Holosphere is an exceptionally well-engineered codebase. It is a compile-time-resolution-templated C++17/20 rendering engine that runs bit-faithfully across three targets (Teensy 4 firmware, Emscripten/WASM, and native test host), backed by a deterministic arena allocator, a fail-fast invariant discipline, and a native test suite whose rigor (oracle-based property tests, an out-of-process death harness, energy/topology invariants) is well above the norm for graphics or embedded code. The documentation — both the 2,100-line README and the inline doxygen — is among the most accurate and genuinely explanatory this review has encountered; comments consistently justify *why* a non-obvious decision is load-bearing rather than restating the code.

The review surfaced **one High-severity latent defect** (a stale-buffer-reuse path in `MeshState::clone` that is not currently triggered by live callers but is real and corruption-class), a modest set of Medium issues (a handful of genuine edge cases, two avoidable hot-path costs, several real test-coverage gaps, and one CI gate that doesn't verify what it claims to), and a long tail of Low-severity quality, documentation, and robustness nits. No defect was found that produces wrong output on common input in the shipped configuration.

---

## Letter Grades by Quality Dimension

| Dimension | Grade | Rationale |
|---|---|---|
| **Correctness & Robustness** | **A−** | Rigorous edge-case handling throughout (NaN-folds, clamp-before-cast, divisor floors, seam/pole degeneracy). Held back by one confirmed High latent corruption path (`MeshState::clone`) and a few Medium edge cases (Feedback seam re-center, MindSplatter underflow clamp). |
| **Architecture & Design Elegance** | **A** | Compile-time `<W,H>` specialization, the persistent+scratch arena model with explicit `Arena&` threading, the World/Screen/Pixel filter-domain pipeline, and the bounded type-erasure boundaries are a coherent, principled whole. Three targets share one engine with no `#ifdef` sprawl thanks to `platform.h`. |
| **API / Interface Expressiveness** | **A** | Fluent `SolidBuilder`, designated-initializer param structs, concept-checked callables, RAII guards (Canvas, ScratchScope, Persist, correction guards), and single-source-of-truth X-macros (`HS_RESOLUTIONS`, `HS_EFFECT_LIST`, `MESHOP_LIST`) that make drift a build error. |
| **Code Quality & Readability** | **A** | Dense but disciplined; uniform naming; comments are load-bearing and accurate. Minor duplication (IslamicStars/HankinSolids carousel) and a few open-coded helpers that bypass canonical ones. |
| **Memory Safety** | **A−** | The arena model (generation stamps, LIFO scope enforcement, rebind counters, fail-fast on OOM/overflow) is excellent and converts most misuse into a trap. The `MeshState::clone` stale-reuse path is the one blemish. |
| **Concurrency / Real-Time Safety** | **A** | The ISR/main-loop single-writer model, the position-from-time flywheel (structurally cannot drop columns), the 1-wire symbol codec (degrades to "missed, never wrong"), and the daydream worker-pool generation fence are all rigorously reasoned and bounded by `static_assert`/`HS_CHECK`. |
| **Performance & Efficiency** | **A−** | Strong per-pixel discipline (LUT-driven trig, baked palettes, premultiplied SSAA, no hidden heap). Docked for two avoidable hot-path costs (RingShower per-fragment `acos`, RD per-neighbor fixed-point→float reconversion) and several sliders that expose genuinely expensive worst cases with no power-aware ceiling. |
| **Testing & Verification** | **A−** | Oracle-based property/golden tests, a verified death harness, dual-resolution effect coverage, and self-gating CI tooling. Real Medium gaps: complex/large-mesh rasterization, Hankin topology, non-Ring SDF pixel output, OKLab vs external reference, and large-N k>1 KDTree. |
| **Build System & CI** | **A−** | Drift is engineered out, not policed; CI runs runtime smoke, stack-margin, provenance, and self-tested size gates. Docked for a param-order gate that checks length not order, and a fast-math compile/link asymmetry unverified on the shipped artifact. |
| **Documentation** | **A** | README and doxygen are exceptionally accurate and explanatory, including subtle host/device divergences. A few comments drift from code (noted in Low findings). |
| **Error Handling & Fail-Fast Discipline** | **A** | The `HS_CHECK` doctrine is applied at exactly the right cold seams, with reject-and-log (never trap) at the untyped JS/WASM boundary. Consistent and well-reasoned. |
| **Portability / Platform Abstraction** | **A** | `platform.h` reproduces FastLED integer semantics bit-faithfully; every intentional host/device divergence is documented at its site. |
| **Overall** | **A−** | A mature, cohesive, expertly-documented system whose defects are concentrated in latent edge cases, perf ceilings, and test-coverage gaps rather than shipped-path bugs. |

---

## Prioritized Fix List

Each item is numbered sequentially. Items are grouped by priority; fix lower numbers first.

### Priority 1 — High

1. ✅ **`core/spatial.h` `MeshState::clone` retains stale geometry on a reused destination.** When the destination `MeshState` is reused (the documented "bounce between arenas / compaction" use case) and the *source* has an empty-but-bound `face_offsets` or `topology`, the conditional copies (`if (fo_size > 0) copy_vector(...)` and `if (!src.topology.is_empty()) copy_vector(...)`) are skipped, so the destination keeps its previous bound array and size. Because the owned-first accessor `get_face_offsets_size()` discriminates on `is_bound()` (not emptiness), the reused mesh then reports stale offsets/topology that don't match its new faces — corrupt face indexing / potential OOB on the next rasterization. Currently latent (live callers placement-new fresh `Transients`, and `Persist` resets the target first), but it is a general-purpose primitive advertised as safe for exactly this reuse. Note: the generic `MeshOps::clone` in `core/mesh.h` is *not* affected — it calls `copy_vector` unconditionally (rebinding even at size 0). Fix: in the `MeshState` specialization, unconditionally rebind `face_offsets`/`topology` to the source size (including 0), or explicitly clear the owned vectors in the empty branch so the accessors fall through to the cleared state.

### Priority 2 — Medium

2. ✅ **`core/filter.h` `Feedback::flush` re-centers neighbor taps against `d00` only.** The three neighbor displacements `d10/d01/d11` are wrapped to within `HALF_WQ` of the `d00` anchor; if `d00` is itself the seam outlier (sits just past the θ=0 wrap while the other three cluster on the far side), all three are pulled toward the outlier and the bilinear blend warps in a one-pixel band at the seam. Fix: re-center against the median of the four taps, or anchor on the tap nearest the sampled pixel's own column.

3. ✅ **`effects/MindSplatter.h` clamp is a no-op in exactly the case it guards.** `p_idx = std::min<size_t>(p_idx, particle_system.active_count - 1)` underflows to `SIZE_MAX` when `active_count == 0`, so the clamp cannot bound the index — and the preceding `assert` is stripped under `NDEBUG` on device, making this the sole device-side guard. Latent OOB read of `pool[SIZE_MAX]` if the "fragments only exist for live particles" precondition is ever violated. Fix: `if (active_count) p_idx = std::min<size_t>(p_idx, active_count - 1);`.

4. ✅ **`scripts/wasm_smoke.mjs` param-order seam check verifies length, not order.** The gate asserts only `values.length === defs.length`. A reordering bug in `param_marshal.h` — the single source of parameter order whose entire purpose is to prevent index drift — would keep equal lengths while transposing every slider and still pass green. Fix: also assert `values[i]` matches `defs[i]` by name/value (accounting for the bool→`>0.5` collapse) for each `i`, which actually proves the two streams are zippable.

5. ✅ **`targets/wasm/wasm.cpp` exposes an unguarded, untested empty-mesh `MeshOps` constructor.** `.constructor<>()` lets JS do `new Module.MeshOps()`, yielding a wrapper over a default-constructed empty `PolyMesh`; operators run `check_live()` (generation guard) but do not reject an empty mesh, so a zero-vertex mesh flows straight into `MeshOps`. The smoke test only ever constructs via `fromSolidName`, leaving this path untested. Fix: drop the public default constructor (all real construction goes through `fromSolidName`), or reject an empty mesh at the operator boundary.

6. ✅ **`daydream/gui.js` destroying a child folder cancels the shared parent URL-writer.** Child folders share the root's single `urlWriter` (`wrapped.urlWriter = this.urlWriter`); `DeepLinkGUI.destroy()` calls `urlWriter.cancel()`, which `clearTimeout`s the shared 200 ms debounce and clears all pending writes — so tearing down one child drops sibling/parent URL writes queued in that window. Fix: only the root owns the writer lifecycle (guard `cancel()` on `this.parent === null`), or ref-count it.

7. ✅ **`effects/RingShower.h` computes a per-fragment `acos` on a rasterizer-bound path.** `fragment_shader` calls `angle_between(z, v) / PI` per ring fragment with `z` constant, purely to index a palette by radial coordinate. An inverse-trig per fragment is the most expensive way to get a monotonic radial key. Fix: index the palette from `dot(z, v)` directly (or a cos-domain remap); reserve `acos` for when the true angle is needed.

8. ✅ **`effects/GSReactionDiffusion.h` / `BZReactionDiffusion.h` reconvert each lattice value ~6× per substep.** The neighbor loop calls `from_q16`/`from_q8` on every neighbor read, so each node's fixed-point value is converted to float once for itself and again each time it appears as a neighbor — ~6–7× redundant conversions across the lattice, over RD_N nodes × 16 substeps × every frame. Fix: pre-convert the current generation into a float scratch buffer once per substep, then read floats in the neighbor loop.

9. ✅ **Test gap: mesh rasterization is verified only on a single octahedron.** `tests/test_mesh_raster.h` exercises `Plot::Mesh::draw` and `Scan::Mesh`/`SDF::Face` solid fill only on the 8-triangle octahedron; no cube/dodecahedron, no non-triangular faces, and none of the high-face-count Conway/Hankin/IslamicStars meshes (the actual effect payload, and the rasterizer-bound hot path) are ever drawn to asserted pixels. Mixed-face and near-pole face-walk/bounding-cull paths are pixel-unverified. Fix: add a cube/dodecahedron and at least one large Conway/Goldberg mesh to the raster test, asserting full sphere tiling + edge-arc membership.

10. ✅ **Test gap: Hankin output meshes get no topology check.** Unlike Conway and solids (Euler + manifold edge-degree + winding), `tests/test_hankin.h` validates compiled Hankin output only with structural smoke (face-count consistency, index range, loose unit-length). Hankin builds the most likely of the three families to produce a non-manifold or inconsistently-wound mesh, and that exact class of bug is unguarded. Fix: run `check_euler_genus0` and a winding check on compiled Hankin output for the cube and icosahedron seeds.

11. ✅ **Test gap: non-Ring SDF primitives have no pixel-level test, and there is no overdraw/composite test.** Only `Scan::Ring` is rasterized end-to-end into asserted pixels; `Scan::Star`, `Scan::PlanarPolygon`, and `Scan::Flower` appear only in `distance()`/cull-coverage checks, never drawn to a Canvas. No multi-shape overlapping-stroke composited-color pixel test exists anywhere. Fix: add a placement-oracle pixel test for Star/Polygon/Flower and a two-stroke overlap blend check.

12. ✅ **Test gap: OKLab/OKLCH are never pinned against external reference coordinates.** `tests/test_color.h` validates the perceptual space only by sRGB→space→sRGB round-trip and named-palette goldens — a wrong-but-invertible M1/M2 matrix pair would pass every assertion. Fix: pin a few published sRGB→OKLab/OKLCH reference triples (e.g. pure red → its documented L,a,b).

13. ✅ **Test gap: KDTree k>1 is not cross-checked against brute force on generic geometry.** Brute-force agreement is tested at k=1 (N=16) and at k>1 only on coincident/collinear points where the tree can't really err. There is no large-N (~100) random-distinct-point k>1 sorted-neighbor match, so bbox-pruning/recursion correctness for k>1 on generic geometry is untested. Fix: add a k>1 sorted-set brute-force match on ~100 randomly-placed distinct points.

### Priority 3 — Low

14. **`daydream/segment_controller.js` `composite()` self-heals only one of two view aliases.** Unlike the single-engine path (which re-points both `Daydream.pixels` and `dotMesh.instanceColor.array`), `composite()` re-points only `Daydream.pixels`. No reachable trigger was found (composite routes through `host.refresh()`, which re-points both aliases together), so this is a latent defensive-symmetry gap rather than a live bug. Fix: mirror both aliases in the divergence branch for consistency with the single-engine heal.

15. **`effects/ChaoticStrings.h` binds `vertices` to exactly `MAX_FRAGMENTS` with zero slack.** `deep_tween` can push up to `TRAIL_LENGTH * ORIENTATION_SUBSTEPS` fragments into a buffer of exactly that capacity; any off-by-one in the emitted sample count trips the arena fail-fast trap. Fix: bind with a small slack (`+ ORIENTATION_SUBSTEPS`) or assert the post-loop tween count.

16. **`effects/Dynamo.h` `color()` has no local bound; the size invariant lives in `draw_frame`.** `color()` reads `baked_palettes_[i+1]` whose in-bounds-ness is guaranteed only by an `HS_CHECK` in `draw_frame`. The in-frame ordering is currently safe (the check precedes `draw_nodes`), so this is an invariant-locality smell, not a live OOB. Fix: add `HS_CHECK(i + 1 < baked_palettes_.size())` inside `color()` to make the bound local.

17. ✅ **`hardware/dma_led.h` cache-flush comment misattributes the flush, and the raw-pointer API has no enforced precondition.** The `transferComplete_` doxygen claims the `arm_dcache_flush` happens inside `transmitAsync()` before `dma_.enable()`, but the flush is actually in the caller `submitFrame()` via `frames_[back].flush()`. `transmitAsync(ptr, len)` takes a raw pointer a future caller could invoke without flushing. Fix: correct the comment to state the flush is the caller's contract (document it on the signature), or move an idempotent flush into `transmitAsync()`.

18. **`effects/SphericalHarmonics.h` linearly blends raw field values mid-morph, causing transient dimming.** The morph does `val += (val2 - val) * blend` on the scalar harmonic values, then shades from the blended `val`; when the two modes have opposite-sign or differing-magnitude fields, the blend passes through near-zero mid-transition and dims the sphere. Fix (optional): crossfade the two fully-shaded colors instead of the raw field, or document the value-blend as intentional.

19. ✅ **`effects/Raymarch.h` Fresnel rim term uses N·L instead of N·V.** `fresnel = 1 - clamp(ndotl, 0, 1)` reuses the diffuse light·normal term; a physically-correct rim keys off the view angle (N·V). It is correct only under the headlight model, but `light` is explicitly tilted off-axis, so the rim follows the light rather than the silhouette. Fix: compute Fresnel from `dot(normal_w, view_dir)`.

20. **Several effects expose expensive slider worst cases with no power-aware ceiling.** `Moire` draws up to 2×H `DistortedRing` rasterizations/frame at the density cap (≈288 at H=144); `Raymarch` runs up to 20 sphere-traces × 30 steps + per-halo occlusion probes; `DreamBalls` transforms up to copies×vertices per frame. All are correct and bounded but reach genuinely heavy worst cases. Fix: cap totals across layers/copies, or gate the heaviest probes behind a quality factor; at minimum document the worst case at each cap site.

21. **`effects/MindSplatter.h` mixes `assert` and `HS_CHECK` for the same invariant family.** The per-fragment index bound uses a stripped `assert` (line ~278) while a sibling capacity check uses always-on `HS_CHECK` (line ~293). Per the fail-fast doctrine, the index bound that actually guards `pool[p_idx]` should be the always-on trap. Fix: promote the per-fragment bound to `HS_CHECK` (or back it with a real clamp — see item 3).

22. ✅ **`core/sdf.h` `Subtract` empty-B fast path skips seam normalization that its sibling branches apply.** The empty-B passthrough emits A's intervals raw while the other two branches run `normalize_intervals_to_range<W>` first; harmless because `scan_region` re-normalizes downstream, but the asymmetry invites a future one-sided edit. Fix: route the passthrough through the same normalization, or document `scan_region` as the single seam authority.

23. ✅ **`core/scan.h` full-row detection comment overstates its guarantee.** `scan_region` tests `iv.second - iv.first >= W` on raw pre-normalize intervals, which catches a single full-circle span but not full coverage assembled from abutting pieces; the comment implies the latter. Output stays correct (the slow path paints every covered column); it is a missed fast-path / doc overstatement only. Fix: tighten the comment, or detect full coverage post-coalesce.

24. ✅ **`core/geometry.h` `Orientation::upsample` silently produces no motion smear from a single source frame.** When `old_num_frames == 1` (the common post-`set()`/`collapse()` state), every output frame collapses to `old[0]`. This is mathematically the only well-defined result (no motion to interpolate), so it is correct behavior. Fix: document that a single-frame trail upsamples to a flat smear; real motion blur requires ≥2 pushed frames.

25. ✅ **`core/conway.h` `gyro` carries two unreconciled doc strings.** It is implemented as `dual(snub(...))` and documented "dual of snub," while the `SolidBuilder` wrapper calls it "pentagonal chiral subdivision." The descriptions are compatible (recipe vs geometric effect) but read as different operations. Fix: cross-reference the two, or align the wording.

26. **`daydream/vendor-importmap.js` "CDN fallback" naming is inaccurate.** The file header and README §10.8 describe a "local-first / CDN-fallback" mechanism, but vendor resolution is baked at build time to a single source with no runtime fallback — if the chosen source is unreachable, every import hard-fails. The code comment itself notes the old runtime probing was removed. Fix: rename the concept to "build-time-baked vendor resolution" in the README and doc comment.

27. ✅ **`scripts/capture_screenshots.mjs` has a dead "WRONG EFFECT" log branch.** Control reaches the save/log line only when `honored === true` (the `!honored` case `continue`s earlier), so the inline `${honored ? '' : ', WRONG EFFECT (fell back)'}` can never print. Fix: drop the dead ternary, or restructure so a tolerated fallback can reach it.

28. **`tools/teensy_gate.py` does not assert exactly one definition per layout symbol.** The match filter excludes zero-size and UND rows but never checks `len(matches) == 1`, so a benign duplicate definition (vague-linkage across TUs) yields duplicate violations rather than a surfaced "unexpected duplicate." Fix: assert a single match (or dedup by value/size).

29. **`.github/workflows/ci.yml` `pio-objs` cache key guarantees a primary-key miss every run.** The key embeds `${{ github.sha }}`, so it never restores an exact hit (always falls to `restore-keys`) and saves a fresh full-size entry per commit. This is a valid rolling-cache idiom but is undocumented and doubles storage churn. Fix: document the rolling intent like the toolchain cache, or key on a source content hash.

30. **`CMakeLists.txt` passes `-ffast-math -fno-finite-math-only` at compile but not at link, and the wasm LTO NaN-survival path is unverified in CI.** The in-file comment correctly concedes both points; the shipped artifact's finite-math correctness rides on per-function IR attributes surviving the LTO link, which no test exercises. Fix: add the pair to `target_link_options` too (a documented no-op today that removes the failure mode), or add a wasm-side smoke assertion that a NaN coordinate is clipped, not wrapped.

31. **`daydream/driver.js` `strobeColumns` is uninitialized, causing a one-frame wrong-mode paint.** The constructor never sets `this.strobeColumns` (only `setStrobeColumns()` does), so the first render before any effect binds reads `undefined === false` and takes the strobe branch. Fix: initialize `this.strobeColumns` in the constructor.

32. **`daydream/driver.js` device-pixel-ratio is not re-applied on resize.** `setPixelRatio(min(dpr, 1))` is set once in the constructor; `setCanvasSize()` calls `setSize` but not `setPixelRatio`, so moving the window to a different-DPR monitor keeps the stale ratio until reload. Fix: re-call `setPixelRatio` in `setCanvasSize()` (or confirm the cap-at-1 is intentional and document it).

33. **`daydream/recorder.js` leaks a capture track if `MediaRecorder` construction throws.** In the `requestFrame`-absent fallback, a new `captureStream` track is created; if `new MediaRecorder(stream, options)` then throws (unsupported options), `start()` unwinds with that track still live. Fix: wrap the construction in try/catch that stops the stream's tracks on failure.

34. **`daydream/state.js` `roundUrlNumber` has no non-finite guard.** `parseFloat(value.toFixed(4))` on `NaN`/`Infinity` serializes the literal `"NaN"`/`"Infinity"` into the URL. Latent (no current caller passes non-finite). Fix: `return Number.isFinite(value) ? parseFloat(value.toFixed(4)) : 0;`.

35. **`daydream/tools/solids.html` Hankin op seeds a degenerate 0° angle from a null mesh.** `addOp('hankin')` derives its angle from `computeInternalAngle(currentMesh)`, which returns 0 for a null/load-failed mesh, unconditionally overriding the `OP_DEFS.hankin` default of 54°. Fix: fall back to the default when `computeInternalAngle` returns 0 or the mesh is null.

36. **`effects/IslamicStars.h` and `effects/HankinSolids.h` duplicate the carousel/palette/topology-draw machinery.** Both declare the same `palettes_slots[2]` pattern and near-identical `MeshCarousel` + `shuffle_indices` + `classify_faces_by_topology` + `Persist`-compaction + `draw_shape`/`draw_mesh` bodies. Fix: extract a shared `TopologyMeshDrawer` helper.

37. **`effects/DistortedRing.h` documents a `thickness` default that is never observed.** `Params::thickness` is doc'd as default `1.0f`, but `init()` overwrites it to `4.0f * px` before `registerParam` snaps the slider default. Fix: drop the misleading `= 1.0f` or comment it as reseeded in `init()`.

38. **`effects/Voronoi.h` open-codes a 3-channel cross-cell color lerp.** The per-channel `sec + (best - sec) * t` with manual `uint16_t` casts duplicates the canonical `Color4::lerp`/`Pixel` blend and risks divergence from it. Fix: use the existing color lerp helper.

39. **`effects/GnomonicStars.h` sources the warp speed twice.** The pinned warp is spawned with `params.warp_speed` and then `draw_frame` mirrors `warp_->speed = params.warp_speed` every frame, so the spawn-time value is immediately superseded. Fix: pass a don't-care at spawn and rely on the per-frame mirror, or comment the redundancy.

40. **`effects/PetalFlow.h` has a misplaced doxygen block.** The block describing `init_timeline` sits above `build_shift_table`, so Doxygen mis-associates it and `init_timeline` gets no doc of its own. Fix: move the block down to directly precede `init_timeline()`.

41. **`core/memory.cpp` arena-failure log tags are inconsistent.** `configure_arenas` logs over-subscription as `[FATAL]` while `Arena::allocate` logs OOM as `[OOM]` and mesh/Hankin traps use bare messages. Fix: standardize the prefix across arena failure paths.

42. **Test gap: `fib_spiral` distribution property is unverified.** `tests/test_geometry.h` checks unit-length, determinism, and opposite-hemisphere endpoints but not the golden-angle spacing / even-distribution that is the spiral's defining purpose; a wrong-but-unit-length spiral would pass. Fix: add a nearest-neighbor-spacing or golden-angle-increment assertion.

43. **Test gap: POV strip interior ordering and N>4 segment layouts are untested.** The POV map tests assert endpoint LEDs and an exactly-once bijection but not monotonic interior ordering within a half-strip (a scrambled-but-bijective interior would pass), and segmented tiling is only checked for N∈{2,4}. Fix: assert monotonic interior LED indices, and add an N=6/8 tiling sweep if the driver supports it.

44. **Minor documentation drifts (batch).** (a) `core/3dmath.h` / `core/geometry.h`: `Spherical::theta` is doc'd `(-π,π]` but `fast_atan2` attains `-π` → `[-π, π]`. (b) `core/3dmath.h` `make_rotation` precondition is over-strict — any non-zero axis works (the final `normalized()` corrects magnitude); only a zero axis traps. (c) `core/platform.h`: the host/device divide-by-zero parity comments depend on the Cortex-M7 `DIV_0_TRP` bit staying clear — cross-reference or assert it at startup. (d) `core/animation.h`: `kLargestConcreteAnimSize` is a hand-maintained list a new non-templated animation can silently escape; consider driving it from an X-macro like `HS_EFFECT_LIST`. (e) `hardware/pov_sync.h` `SymbolEmitter::tick()` interleaves two "queue empty" encodings; express "no burst in flight" once.

---

## Claims Raised but Dropped After Verification

For transparency, the independent verification pass refuted or reclassified these initially-reported items, and they are **not** counted as defects:

- **`core/presets.h` `int`/`size_t` accessor mismatch** — refuted; indices are always within a small `std::array` extent, so the narrowing is benign/cosmetic.
- **`core/animation.h` `Lerp` "illusory" generic easing param** — refuted; the `template<typename Easing>` is the standard convertible-callable idiom, not a promise of capturing-lambda support.
- **`core/3dmath.h:1222` "stray backslash" breaking a comment** — refuted; no backslash exists (the search artifact was the UTF-8 em-dash byte sequence).
- **`core/styles.h` `downsample` divide-by-zero** — already handled; `filter.h` guards it with both a `static_assert` and a runtime `HS_CHECK` that fail-fast on 0.
- **`daydream/recorder.js` segmented frame-sync drop/dup** — refuted; `advanced` and the composite share one clock (`stepSimulation` runs only when `advanced`), so no real frame is dropped.
- **`daydream/tools/palette_math.js` LCG precision loss** — refuted; the `state * 1664525 + increment` product stays below 2⁵³ for every uint32 input, so it matches a true uint32 LCG (a `Math.imul` rewrite would be idiomatic but changes nothing).
- **`daydream/segment_controller.js:518` GPU-aliasing bug** — reclassified from High to Low; no reachable trigger exists (see item 14).
