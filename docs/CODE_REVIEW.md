# Holosphere / daydream — Code Quality Review

**Date:** 2026-06-30
**Scope:** The Holosphere C++ rendering engine + firmware (`core/`, `effects/`, `hardware/`, `targets/wasm`, `targets/Phantasm`, native test suite, build/CI tooling) and the daydream web simulator (app/bridge, Three.js renderer, state/GUI, recorder, segmented-POV worker pool, geometry tools, JS tests).
**Out of scope:** `core/effects_legacy.h`, `targets/Holosphere/Holosphere.ino`, `core/rotate.h`; vendored/generated artifacts (`FastNoiseLite.h`, `core/color_luts.h`, `holosphere_wasm.{js,wasm}`, `three.js/`, `node_modules/`).

## Methodology

The codebase was partitioned into **24 components**, each audited by a dedicated reviewer that read every in-scope file against the architecture described in the README. Each candidate finding was then handed to a **fresh, independent validator** that re-read the cited code from scratch and ruled the finding **Confirmed**, **Adjusted** (severity corrected), or **Rejected** (non-issue, already handled, out of scope, or a disproportionate fix), applying the eligibility bar from the `code-review-fix` skill.

- **128** candidate findings gathered → **85 rejected** on independent validation → **43 retained**.
- Severity of retained findings: **0 Critical · 0 High · 3 Medium · 40 Low.**

The high rejection rate is itself a signal: the great majority of "issues" a first pass surfaces in this codebase turn out to be already-guarded, already-documented, or intentional design tradeoffs. What remains is almost entirely latent-robustness hardening, documentation precision, and test-coverage gaps — not defects on shipped paths.

---

## Overall Grade: **A−**

A mature, unusually disciplined codebase. The engineering is defensive by construction: a heap-free partitioned arena with RAII scoping, an always-on fail-fast trap discipline (`HS_CHECK`) that is itself verified by a subprocess death-test harness, compile-time resolution specialization, a lock-free double-buffer with a rigorously reasoned index protocol, and a C++↔WASM↔JS boundary held together by single-source-of-truth macros with anti-drift CI checks. Documentation is exceptional — nearly every epsilon, sentinel, and hot-path contract is justified inline. No correctness, memory-safety, or concurrency defect was found on any shipped path. The gap between this and an unqualified **A** is made up of the residue below: a handful of latent tripwires that fire only under future extension, a few doc/comment inaccuracies, and some specific test-coverage holes (device-resolution arena peaks, the volumetric occluder path, the runtime importmap builder).

### Quality Dimensions

| Dimension | Grade | Rationale |
|---|---|---|
| **Architecture & Design** | **A** | Clean generate→transform→rasterize→filter pipeline; variadic compile-time filter chains; shared engine across three targets; two-repo product with a build that keeps both in sync. Elegant and consistent across the board. |
| **Interface Expressiveness** | **A−** | Named factories (`from_spherical`, `normalized_or`), enum-driven dispatch, and explicit `Arena&` parameters make contracts hard to misuse. Docked for a few dead/misleading parameters (`RingSpin` normal, `MobiusGrid` Z-axis seed) and the `needs_full_frame` opt-in footgun. |
| **Correctness & Robustness** | **A−** | Singularity handling (poles, antipodes, degenerate axes/quaternions, θ=0 seams) is meticulous and well-tested. Residual issues are latent threshold-coarseness and one real UX defect (Liquid2D pause). |
| **Performance & Efficiency** | **A** | Split trig LUTs (~145× memory cut), fast-trig/cbrt hot paths, branchless wrap, instanced single-draw GL rendering, zero-copy WASM memory views. Tradeoffs are measured and documented. |
| **Memory Safety & Resource Mgmt** | **A** | Heap-free arena, trivially-copyable value types, placement-new + launder, fixed-capacity containers with generation tracking, `HS_CHECK` bounds at cold seams. A couple of value-init/guard nits only. |
| **Concurrency & Real-Time** | **A−** | ISR double-buffer and Phantasm flywheel sync are model real-time engineering; the JS worker pool uses a generation fence + single fault latch to break every deadlock path. Docked for documented single-thread LUT-init assumptions and one blocking-in-ISR path on the legacy single-board target. |
| **Error Handling** | **A** | Principled split between always-on fail-fast (`HS_CHECK`, survives `NDEBUG`) for invariant violations and soft-degrade for transient conditions; the WASM boundary rejects-and-returns rather than aborting. |
| **Testing & Verification** | **A−** | 45 cross-pinned native modules + death harness + effect smoke/determinism sweep; 276 passing JS tests wired into the deploy gate. Docked for specific holes: device-H arena peaks, `Scan::Volume` occluder logic, the runtime importmap IIFE, and a few streaming-recorder edge cases. |
| **Documentation** | **A** | Among the best seen in a hobby-scale engine: a 2,000-line architectural README plus inline rationale for nearly every non-obvious constant and contract. A few comments have drifted from the code (docked accordingly). |
| **Readability & Naming** | **A** | Clear naming, single-source helpers, terse-but-complete doxygen; heavy comment density is justified by the numerical subtlety. |
| **Portability** | **A−** | Clean Arduino/WASM/Desktop abstraction with bit-exact host mocks. Minor latent traps: a host stub signature diverging from the Teensy prototype, an IEEE-754-layout bit-hack. |
| **Security** | **A** | Limited attack surface (local art tool + static web page). Input validation at the WASM/JS boundary is consistent; no injection, unsafe eval, or credential handling of concern. |
| **Maintainability** | **A−** | Single sources of truth (rosters, param ordering, provenance hashes) with anti-drift CI. Docked for scattered magic-constant duplication, a twin include/roster list, and a couple of duplicated helpers in the JS tools. |

### Per-Component Grades

| Component | Grade | Component | Grade |
|---|---|---|---|
| `core` math / geometry | A− | `effects` (a–d) | A− |
| `core` color / palettes | A− | `hardware` drivers | A |
| `core` memory / containers | A | `targets/wasm` bridge | A |
| `core` sdf / scan | A− | native test suite | A |
| `core` plot / filter | A− | build / CI / codegen | A |
| `core` animation | A− | daydream app / bridge | A− |
| `core` mesh / conway / spatial | A− | daydream renderer | A− |
| `core` hankin / solids / generators | A− | daydream state / GUI | A− |
| `core` reaction-graph / noise | A | daydream recorder | A− |
| `core` engine / registry / platform | A | daydream segmented-POV | A |
| | | daydream tools + JS tests | A− |

---

## Prioritized Fix List

Every retained defect is listed below under its priority section, numbered sequentially. Items are phrased so the `code-review-fix` skill can mark each one ✅ (fixed) or ❌ (rejected) in place. There are **no Critical or High findings**; the list begins at Medium.

### Priority: Critical

_None._

### Priority: High

_None._

### Priority: Medium

1. ✅ **Device-resolution arena high-water is unverified by CI** (`Testing`, `effects/HankinSolids.h:31-38`) — The scratch_a (24KB)/scratch_b (32KB) budgets are hand-tuned and the comment states the device H=144 high-water "isn't in CI; confirm the peak on hardware". An under-sized arena traps (fail-fast) at runtime on device, a class of regression CI cannot catch for this effect. _Fix:_ Add a host test exercising the heaviest hankin mesh at device H to gate the scratch peaks, or a static worst-case bound.

2. ✅ **Liquid2D time accumulators keep advancing while animations are paused** (`Correctness`, `effects/Liquid2D.h:60-65, 108-121`) — The accumulated_time and cycle_phase Drivers are added without the `&anims_paused_` gate, and sin_phase/cos_phase are advanced unconditionally in draw_frame by params.time_speed. Touching the animated Time/Cycle sliders engages "Pause Animation", yet the warp and trig motion never freeze, so the user cannot tweak a static frame. Raymarch's Drivers pass `&anims_paused_`. _Fix:_ Pass `&anims_paused_` to both Drivers and gate the sin_phase/cos_phase advance on `!animationsPaused()`, matching Raymarch.

3. ❌ **PiP camera renders the identical main-view image** (`Correctness`, `daydream/driver.js:604-606`) — Rejected: working as intended. The PiP deliberately mirrors the main camera pose to render a scaled-down thumbnail of the main view, so the operator can see how the piece reads from a distance; it is not meant to show a fixed/opposite orientation. The behavior is correct — only the README's "fixed orientation" / "front and back visible simultaneously" description was inaccurate, and it has been corrected to match.

### Priority: Low

4. ✅ **Complex::operator/ divisor threshold can emit STEREO_INF for merely small (not singular) denominators** (`Correctness`, `core/3dmath.h:669-681`) — denom < EPS_LEN_SQ (1e-6) triggers the infinity-sentinel branch. A denominator at |b| ~ 3e-4 (denom ~ 1e-7) is force-snapped to STEREO_INF magnitude, discarding a finite quotient of magnitude ~O(3e3). Mobius chains near a pole may quantize a valid finite point to infinity. _Fix:_ Gate the sentinel on the resulting magnitude exceeding STEREO_INF rather than on denom alone, or document the pole-radius this threshold implies.

5. ❌ **wrap() template NDEBUG fallback masks a non-positive modulo base instead of trapping** (`Correctness`, `core/util.h:45-60`) — Rejected: hot path. wrap() is the branchless SDF angular-repeat inner loop whose comment explicitly forbids a branch; an always-on `HS_CHECK(m > 0)` (survives NDEBUG) inserts a per-call branch into the hottest path — a measurable perf cost to catch a caller logic error the existing debug assert already traps under test builds. The modulo base is a compile-time-ish shape constant, never caller-derived on this path, so the fmax NaN-floor is the correct release behavior and the assert the correct debug diagnostic.

6. ✅ **hue_rotate(amount) overload omits the chroma-preserving renormalization sync_hue applies** (`Correctness`, `core/color.h:780-783`) — The (ca,sa) hue_rotate is fed non-orthonormal fast trig; Style::sync_hue renormalizes (1/hypot) precisely because fast trig scales chroma and "compounds per frame under feedback". The amount overload calls fast_cosf/fast_sinf raw with no renormalization, so any feedback-loop caller of it silently compounds chroma drift. Currently only Flyby (single-pass) uses it, so latent. _Fix:_ Renormalize (ca,sa) inside the amount overload too, or document it as single-pass-only and forbid feedback use.

7. ✅ **Style test pins hue_fade exactly equal to the un-renormalized hue_rotate overload** (`Testing`, `tests/test_styles.h:289-292`) — After s.sync_hue() (which renormalizes hue_ca/hue_sa), the test asserts hue_fade == hue_rotate(...,s.hue_shift) with HS_EXPECT_EQ. The reference overload uses raw non-renormalized fast trig, so the two rotations use different (ca,sa). Exact equality holds only because fast-trig error plus 16-bit quantization round identically at this angle — a fragile pin, not guaranteed math. _Fix:_ Compare against hue_rotate(c, s.hue_ca, s.hue_sa), or use a tolerance, not exact EQ.

8. ✅ **frac_to_q16 doc contract "[0,1)" contradicts Color4::lerp which passes clamped [0,1]** (`Documentation`, `core/color.h:218-220`) — frac_to_q16 documents frac as "assumed in [0,1) (caller-guaranteed; not clamped)", but Color4::lerp passes frac_to_q16(ct) with ct clamped to [0,1] inclusive, and BakedPalette/Gradient endpoints special-case idx>=max before calling. The value at 1.0 (65535.5->65535) is safe, but the stated half-open contract is inaccurate and could mislead a future caller into assuming 1.0 is excluded. _Fix:_ Document the contract as closed [0,1] (which is already safe), matching actual call sites.

9. ✅ **BrightnessProfile value ranges never reach 255 due to half-open rand_int** (`Correctness`, `core/color.h:1124-1145`) — hs::rand_int is half-open [min,max) (per animation_scalars.h:36 comment). Ranges like rand_int(204,255) and rand_int(178,255) can never emit the intended top value 255, so ASCENDING/DESCENDING/CUP keys top out at 254. Purely aesthetic (authoring-time, generative), but the literals read as inclusive and mildly bias brightness down. _Fix:_ Use 256 as the exclusive upper bound where 255 is meant to be reachable, or accept and note it.

10. ✅ **SmoothUnion::get_horizontal_intervals reads TrigLUT before ensuring initialization** (`Correctness`, `core/sdf.h:1023`) — SmoothUnion computes pad_px from TrigLUT<W,H>::sin_phi[y] at the top of get_horizontal_intervals without the `if(!initialized) init()` guard every leaf shape performs. Safe on the live path (scan_region inits first), but fragile: a standalone call (e.g. a future unit test) reads an uninitialized static LUT (all zeros) and computes pad_px=W silently. _Fix:_ Add the standard `if (!TrigLUT<W,H>::initialized) TrigLUT<W,H>::init();` guard at the top, matching the leaf shapes.

11. ✅ **Scan::Volume::draw has no direct unit test for the occluder/AA probe logic** (`Testing`, `core/scan.h:1199-1334`) — test_scan.h covers rasterize, shader, seam de-dup, and one TransformedVolume roundtrip, but the substantial probe_occluder / self-occlusion-AA / soft-corner-fill logic in Volume::draw (the most intricate branchy code in scan.h) has no behavioral test asserting the occluded-edge blend or the local-minimum soft coverage. Regressions there would surface only visually. _Fix:_ Add a Volume::draw test with a two-surface shape asserting the foreground edge antialiases over the background rather than fading to black.

12. ✅ **steps_cache capacity backstop silently degrades sampling with no diagnostic** (`Maintainability`, `core/plot.h:555-560`) — The 2*W steps_cache cap breaks out of the sim loop on pathological segments and stretches cached steps, but there is no counter/trace when this fires. On-device a persistently-hit backstop (e.g. a mis-sized shape wrapping the sphere) would silently coarsen every such curve with no way to notice it in a live show. _Fix:_ Add a debug-only counter/trace (HS_PROFILE or trap under a debug flag) when the backstop trips.

13. ✅ **MobiusWarpEvolving uint32->float phase freeze is undocumented** (`Correctness`, `core/animation_mesh.h:419`) — step() computes `float time = t * speed_` with t a uint32_t converted to float. Past t==2^24 float can no longer represent consecutive frames, so the phase quantizes/freezes — the same limit Noise and RandomWalk document with an explicit caveat, but here it is silent, and this is a perpetual (duration=-1) animation that runs forever. _Fix:_ Add the same 2^24 precision caveat as Noise/RandomWalk, or accumulate a float phase per step instead of scaling raw t.

14. ✅ **Channel-offset magic number inconsistency between noise warps** (`Maintainability`, `core/transformers.h:467-469`) — noise_transform() names its decorrelation shifts kChannelYOffset/kChannelZOffset (100/200) for clarity, but stereo_noise_warp() in the same file re-hardcodes a bare 100.0f for the same purpose. The named-constant convention introduced two functions up is not carried through, inviting drift if one is tuned. _Fix:_ Reuse a shared named offset constant in stereo_noise_warp instead of the literal 100.0f.

15. ✅ **MeshMorph guards empty source but not empty dest** (`Correctness`, `core/animation_mesh.h:278-320`) — The ctor HS_CHECKs `!source.vertices.is_empty()` but not dest. An empty dest silently produces a degenerate morph: correspondence and interpolation loops run zero times, step() SLERPs nothing, and the incoming half draws an empty mesh — no crash but a silently blank crossfade half that is hard to diagnose. _Fix:_ Add `HS_CHECK(!dest.vertices.is_empty())` alongside the source check.

16. ✅ **Boundary-mesh behavior differs silently between relax and other operators** (`Documentation`, `core/conway.h:1189`) — relax tolerates a boundary mesh (skips unpaired-twin vertices, partial relaxation) while dual/ambo/truncate/expand/chamfer/snub call require_closed_manifold and trap. The divergence is documented on relax but there is no single place stating the per-operator manifold contract, making the asymmetry easy to trip over when composing. _Fix:_ Add a one-line manifold-precondition summary table in the operator-block header comment.

17. ✅ **No static_assert that RD_N fits in int16_t table element type** (`Correctness`, `core/reaction_graph.h:27-33, 79`) — neighbors[] is int16_t and node indices are stored as int16_t, but only RD_N>=2 and the D_AVG sync are asserted. Bumping RD_N past 32767 (INT16_MAX) would silently truncate/overflow generated indices. Sibling subsystems (mesh.h, solids.h, spatial.h) guard exactly this INT16_MAX bound. _Fix:_ Add `static_assert(RD_N <= INT16_MAX, "node index must fit int16_t");` next to the RD_N asserts.

18. ❌ **needs_full_frame override is an easily-forgotten footgun with no compile-time enforcement** (`Correctness`, `core/canvas.h:112`) — A stateful effect that reads pixels outside its band but forgets to override needs_full_frame() renders wrong under segmented clipping, with no diagnostic. The base default silently favors the clipping win. Known API hazard; no cheap forget-proof test exists (the only forget-proof fix is a disproportionate CRTP, and a roster-gate test already catches known effects). _Fix:_ None cheap; retain as a documented hazard (likely ❌ on the eligibility bar). _Rejected:_ no cheap forget-proof enforcement exists — the only compile-time fix is a disproportionate CRTP, and the roster-gate test already covers every shipped effect; kept as a documented hazard.

19. ✅ **BZ scratch-arena size is unguarded, unlike ChaoticStrings** (`MemorySafety`, `effects/BZReactionDiffusion.h:84-86`) — BZ/GS static_assert only the persistent budget. Scratch_a is given the remainder (GLOBAL_ARENA_SIZE - persistent) with no static_assert that its per-frame 115KB (BZ) / 92KB (GS) fit. On device that remainder is ~165KB/156KB; a future persistent-arena bump could silently starve scratch and trap only at runtime allocate(). _Fix:_ Add a static_assert that the 6 render() scratch buffers fit in (GLOBAL_ARENA_SIZE - kPersistentBytes), mirroring ChaoticStrings' SCRATCH_A_BYTES guard.

20. ✅ **advance_substeps is a non-static member using no instance state** (`Style`, `effects/ReactionDiffusionBase.h:160-171`) — advance_substeps() references only its parameters and RD_N; it touches no member. land_back() beside it is already static. Keeping it an instance method obscures that the substep driver is pure and invites a false assumption of per-object state. _Fix:_ Make advance_substeps static for parity with land_back and to document its purity.

21. ✅ **Stale TODO preset comment left in shipped source** (`Maintainability`, `effects/Flyby.h:45`) — A commented "TODO: Good preset for later" float list sits in init(). Per the repo's remove-dead-over-document convention, an unused candidate preset embedded as a TODO comment is dead data that should be dropped or promoted into the presets array. _Fix:_ Delete the TODO comment or fold the candidate into the Presets<Params,4> table if it is actually wanted.

22. ✅ **HoleRef filters seeded on Z_AXIS but poles are the Y axis** (`Maintainability`, `effects/MobiusGrid.h:33, 84-87`) — The ctor initializes holeN(Z_AXIS)/holeS(-Z_AXIS), but draw_frame transforms the Y_AXIS poles (n_in=Y_AXIS) and draws rings around Y_AXIS. The Z_AXIS seed is dead: holeN/holeS are recomputed at the top of every draw_frame before any rasterization, so the initial value never affects output but misleads readers about the pole axis. _Fix:_ Seed the holes with Y_AXIS/-Y_AXIS (or drop the seed and document the recompute) so the axis convention is consistent.

23. ✅ **cos_eh array left partially uninitialized** (`MemorySafety`, `effects/MindSplatter.h:251-255, 279-289`) — std::array<float,AttractSolid::NUM_VERTS> cos_eh is default-initialized (indeterminate) then filled only for [0, attractors.size()). Entries beyond the live attractor count hold garbage. They are never read (the vertex shader loops to attractors.size()), so it is currently safe, but relies on that invariant holding. _Fix:_ Value-initialize (`cos_eh{}`) or size to attractors.size() to avoid indeterminate values.

24. ✅ **Ring normal parameter is effectively dead; every ring spawns with Y_AXIS** (`API`, `effects/RingSpin.h:157`) — Ring's ctor takes a plane normal n, but spawn_ring always passes Y_AXIS, so all four rings share the same great-circle plane and initial orientation seed, diverging only via per-ring noise. The parameter implies per-ring plane variety that never exists, and rings start nearly coincident. _Fix:_ Pass a random_vector() normal (or per-ring axis) at spawn, or drop the unused ctor parameter.

25. ✅ **Single-board DMA (USE_DMA_LEDS) path is never instantiated or tested** (`Testing`, `hardware/pov_single.h:26-33,159-172,215-228`) — The shipped Holosphere target builds the FastLED/WS2801 branch (platformio [env:holosphere] has no `-D USE_DMA_LEDS`; Holosphere.ino never invokes HS_DEFINE_POV_SINGLE_LED_CONTROLLER). The entire DMA half of POVDisplay — packPixel loop, submitFrame, the DMAMEM specialization macro — compiles on no target, so it drifts untested against dma_led.h changes. _Fix:_ Either wire a CI/target that exercises the single-board DMA path, or delete the dead USE_DMA_LEDS branch (both options are heavier than a Low fix warrants; likely ❌ or deferred).

26. **Host arm_dcache_flush stub signature diverges from Teensy prototype** (`Portability`, `hardware/hd107s_frame.h:33-35`) — The non-ARDUINO stub is `void arm_dcache_flush(void*, size_t)`, but Teensy's real declaration is `arm_dcache_flush(void*, uint32_t)`. On LP64 hosts size_t is 64-bit while callers pass an int COMPOSITE_SIZE; it works today but a signature mismatch across the two builds is a latent trap if the stub is ever taken as a function pointer or overload-resolved. _Fix:_ Match the Teensy prototype exactly: second parameter uint32_t.

27. **FastLED single-board path issues two blocking FastLED.show() calls per column in the column ISR** (`Performance`, `hardware/pov_single.h:173-184`) — In the non-DMA (shipped) path show_col() calls FastLED.show() then, when strobe_columns, FastLED.showColor(black) — both synchronous bit-banged strip writes — from inside the IntervalTimer ISR. This blocks the ISR for the full strip transfer (the very stall dma_led.h exists to avoid); acceptable only because the 96x20 rig's ~1.3 ms column budget is generous, but it is a hard real-time coupling worth noting. _Fix:_ Document the blocking-in-ISR tradeoff at show_col, or budget-assert the column period against strip transfer time.

28. **bakeLut return view detachment contract is documented on the struct but not on the method callers read** (`Documentation`, `targets/wasm/wasm.cpp:1188-1224`) — PaletteOps::bakeLut returns a typed_memory_view aliasing `lut`, sharing getPixels' detach-on-growth + read-before-next-call contract. The contract is noted on the `lut` member comment but the method's own @return says only "JS Uint8Array view", without the must-consume-before-next-call warning that getPixels/getParamValues carry inline. Easy to miss when reading the method. _Fix:_ Add the memory-view/read-before-next-call caveat to bakeLut's @return, matching getPixels' doc.

29. **getParameterDefinitions bool value uses v.value>0.5f but getParamValues streams raw float for the same param** (`API`, `targets/wasm/wasm.cpp:594-611`) — For an is_bool param, getParameterDefinitions emits a JS boolean (value>0.5), while getParamValues streams the raw float in the same index slot. The two parallel streams therefore carry different JS types for bool params. Consumers already key off the definition type (the value stream is a Float32Array and physically cannot carry a bool), so this is a documentation gap, not a defect. _Fix:_ Document in param_marshal.h that the value stream carries raw floats even for bool params; JS keys off the definition type.

30. **No validation that BLANK_FLOOR <= MIN_LIT, allowing a saved frame flagged BLANK** (`Correctness`, `scripts/capture_screenshots.mjs:41-44`) — MIN_LIT (stop-retry threshold) and BLANK_FLOOR (blank-warning threshold) are independent env-tunable numbers. If an operator sets BLANK_FLOOR > MIN_LIT, a frame that cleared MIN_LIT and stopped early can still be marked "STILL BLANK" and fail the run, a self-contradiction. _Fix:_ Clamp/assert `BLANK_FLOOR = Math.min(BLANK_FLOOR, MIN_LIT)` or warn on inversion.

31. **Doxygen build job has no timeout-minutes; CI convention broken** (`Maintainability`, `.github/workflows/docs.yml:44-75`) — ci.yml bounds its slow generator job with timeout-minutes: 10, citing hang protection. docs.yml's Doxygen build (external tarball download + graph rendering via graphviz) has no job timeout, so a network stall on the curl download or a doxygen hang can occupy a runner up to the 6h default. _Fix:_ Add timeout-minutes to the docs build (and deploy) jobs.

32. **Generator main() ignores self-test check() on the normal emit path** (`Testing`, `scripts/generate_luts.py:152-168`) — check() (monotonicity + round-trip) runs only when --check is passed. The default `python generate_luts.py > color_luts.h` path emits without validating, so a maintainer regenerating locally gets no self-test unless they separately run --check. CI runs both, but the generator's own happy path could emit a table that fails its own invariants. _Fix:_ Run check() before emitting on the normal path and abort on failure.

33. **Test-All ticker can livelock on a persistently-rejected effect** (`Correctness`, `daydream/daydream.js:617-632, 619-627`) — The interval advances via nextIndex = (indexOf(currentEffect)+1)%len. When the next effect's setEffect() is rejected, the subscriber reverts appState to `old`, so currentEffect stays put and the ticker recomputes the same rejected next index every second — cycling between two effects instead of skipping the bad one. _Fix:_ On a rejected effect during Test-All, advance past the failed index (track/skip it) rather than reverting to the same predecessor.

34. **OrbitControls lacks minDistance; zoom clips sphere at CAMERA_NEAR** (`Correctness`, `daydream/driver.js:183, 101-102`) — The main controls set no minDistance/maxDistance (default minDistance=0). With CAMERA_NEAR=100 and the sphere surface some scene units out, a user can zoom the camera close enough that the near plane clips into the sphere front, cutting away the model. _Fix:_ Set controls.minDistance to keep the camera outside CAMERA_NEAR + SPHERE_RADIUS (and a maxDistance under CAMERA_FAR).

35. **controls.update() comment cites damping/auto-rotate that are never enabled** (`Documentation`, `daydream/driver.js:408`) — The comment "Must run every frame for damping/auto-rotate; emits change" overstates behavior: enableDamping and autoRotate are never set (both default false), so controls.update() only does work during active pointer interaction. The per-frame call is near-inert and the rationale is misleading. _Fix:_ Either enable damping/auto-rotate to match the comment, or trim the comment to reflect that update() only services live interaction.

36. **Fallback makeUrlParamWriter commit path lacks a location.hash guard** (`Correctness`, `daydream/gui.js:50-66`) — The standalone-page fallback commit() rebuilds the URL from pathname+query only, dropping any location.hash — unlike URLSync.flush/reset which explicitly append window.location.hash. A tool page using a fragment loses it on the first GUI change. _Fix:_ Append window.location.hash in commit() to match URLSync's hash-preserving behavior.

37. **Timed-fallback captureStream samples the offscreen before any frame is blitted** (`Correctness`, `daydream/recorder.js:131-143`) — On the requestFrame-less fallback, captureStream(fps) begins sampling the offscreen on wall-clock immediately at recorder.start(), but the offscreen is only filled by captureFrame() blits. The leading interval encodes a blank (cleared) offscreen until the first captureFrame arrives. _Fix:_ Do one blit into the offscreen before recorder.start() on the fallback path, or drop leading blank frames.

38. **Streaming finish() catch can silently lose data when close() (not write) throws** (`Correctness`, `daydream/recorder.js:388-407`) — If streaming succeeded (failed=false) but the final writable.close() rejects, control falls to the catch which only downloads chunks if chunks.length. On the streaming path chunks is empty, so a close failure yields no download and only a console.warn, without the truncation error the write-failure path emits. _Fix:_ Report the close failure as a possible-truncation error like the mid-stream write path does.

39. **No test covers picker-cancelled-after-N-buffered-chunks streaming fallback** (`Testing`, `daydream/recorder.js:342-410`) — The streaming sink relies on a promise chain ordering writes before finish. Tests await a single microtask flush and never exercise the picker-cancelled mid-session fallback (writable null after some writes were already queued and buffered into chunks). _Fix:_ Add a test for picker-cancelled-after-N-chunks buffering fallback and drain the chain deterministically.

40. **Runtime importmap builder (IIFE) has no unit test** (`Testing`, `daydream/vendor-importmap.js:22-72`) — generate_importmap.test.js covers only the build-time generator. The browser-side IIFE — self-path detection via document.currentScript, core-key-wins EXTRA merge, local-vs-CDN URL assembly — has zero direct coverage. Its logic is non-trivial and a regression (e.g. dropping the core-key precedence) would ship silently. _Fix:_ Add a Node test that evals the IIFE with a stub document/currentScript and asserts the imports map for cdn and local VENDOR variants.

41. **generateFuncAndRecipe assumes item.ops is present** (`Correctness`, `daydream/tools/solid_codegen.js:95`) — item.base is validated as a C++ identifier, but item.ops is dereferenced via .forEach without a guard. A spec missing ops (or with ops non-array) throws an opaque "Cannot read properties of undefined (reading forEach)" instead of the module's clear, named errors used everywhere else. _Fix:_ Guard: `if (!Array.isArray(item.ops)) throw new Error('generateFuncAndRecipe: item.ops must be an array')`.

42. **Duplicated `/* Navigation Tabs */` comment** (`Style`, `daydream/styles/index.css:11-12`) — The "Navigation Tabs" section comment is repeated on two consecutive lines. Harmless but a copy-paste artifact. _Fix:_ Delete the duplicate comment line.

43. **showFatalError duplicated between palettes.html and shared.js** (`Maintainability`, `daydream/tools/palettes.html:388-398`) — palettes.html inlines its own showFatalError (to avoid pulling three via shared.js), duplicating the banner logic and inline styles from shared.js. The two can drift. Extracting the banner into the dependency-free clipboard.js (already imported here) would let both share one implementation. _Fix:_ Move showFatalError into a THREE-free module (e.g. clipboard.js or a new banner.js); import it in both shared.js and palettes.html.

---

## Closing Note

This review found **no defect that corrupts output, crashes a shipped path, or races**. The 43 retained items are, almost without exception, hardening a codebase that is already correct: converting latent tripwires into compile-time asserts, tightening documentation that has drifted a line or two from the code, and closing specific test-coverage gaps. That the independent validation pass rejected two-thirds of candidate findings — overwhelmingly because the "problem" was already guarded, documented, or intentional — is the clearest single measure of the codebase's maturity. The engineering discipline on display (fail-fast invariants verified by a death harness, single-source-of-truth macros with anti-drift CI, bit-exact host mocks, a heap-free deterministic arena) is well above what the hobby/art-installation category typically exhibits.
