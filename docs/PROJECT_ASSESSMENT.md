# Holosphere + daydream — Project Assessment

**Date:** 2026-07-26
**Basis:** a 28-agent audit of both repositories, covering ~55k lines of hand-written C++, ~49k lines
of native tests, ~14k lines of JavaScript, ~6k lines of build/profiling tooling, and the hardware
generation pipeline. 4,575 commits since 2015-07-27.

---

## 1. What this actually is

It would be easy to file this under "LED art project." That classification is wrong by roughly two
orders of magnitude of engineering content.

Holosphere is **a purpose-built real-time 3-D rendering engine whose output device is the surface of
a sphere**, plus the firmware to drive that surface from four synchronized microcontrollers, plus a
browser simulator that runs the identical C++ compiled to WebAssembly, plus the toolchain to profile,
gate, and reproduce all of it. The LED hardware is the display; the engine is the work.

The engine contains, as first-class subsystems: a signed-distance-field shape library with analytic
scanline bounds defined on the sphere rather than the plane; an adaptive-step curve rasterizer with
2-D screen-velocity tracking; a compile-time variadic filter pipeline that lifts and projects between
world, screen, and pixel domains automatically; a 16-bit linear-light color system with OKLCH
perceptual interpolation and first-exit gamut mapping; a complete Conway polyhedron operator algebra
with convergent relaxation and baked reproducibility proofs; a Hankin-method Islamic star-pattern
compiler that runs on arbitrary polyhedra and morphs continuously between them; a K-nearest-neighbour
Fibonacci-lattice graph for reaction-diffusion on the sphere; and a partitioned arena allocator with
three-stamp use-after-free detection. All of it runs in 298 KiB on a 600 MHz Cortex-M7 with no heap
allocation on the render path.

---

## 2. Magnitude of the achievement

Four things individually would each be a respectable project. Their combination is what is unusual.

**The rendering architecture.** Sphere-native SDF rasterization is not a common technique — most
spherical displays render to an equirectangular texture with planar tools and accept the pole
distortion. Here the distance functions, the scanline bounds, and the anti-aliasing kernels are
defined in radians on the unit sphere, with each conservative bound carrying its derivation in the
comment above it. The filter pipeline's automatic domain conversion means an effect author writes
`filters.plot(canvas, world_vector, color)` and the compiler decides where the projection happens.
That is a level of abstraction usually found in engines with teams behind them.

**The dual-target discipline.** The same C++ compiles to Teensy firmware and to WebAssembly, and the
project treats *bit-identical output between them* as an enforced invariant — banning `std::shuffle`,
pinning the generator and runtime to double-precision lockstep, recompiling NaN-clamp assertions
under the shipping `-ffast-math` flags in a dedicated translation unit, and maintaining an explicit
device/host divergence ledger. Most projects that ship a simulator let it drift within a month.

**The multi-board synchronization.** The Phantasm 1-wire sync protocol — count-coded, odd-only,
distance-2 symbol alphabet over a shared pin, disciplining a per-board cycle-counter flywheel, with
epoch beacons and a position-weighted checksum — is real protocol design, documented as a signal
datasheet and validated by an event-driven four-board simulator with per-board ppm skew, EMI
injection, and mid-show reboot. This is the part a professional embedded team would be proud of.

**The verification and performance regime.** A 0.90:1 test-to-source ratio; a death-test harness that
verifies the fail-fast traps actually fire; provenance gates that re-derive generated tables rather
than re-hashing them; a device-lock and profile-sweep regime producing per-effect cycle reports; an
ITCM byte ledger; and — most tellingly — written records of *dead levers*, optimizations that were
measured and rejected, so they are not re-attempted. That last habit is rare even in industry.

---

## 3. Technical merit in context

| Compared against | Verdict |
|---|---|
| **Peer hobbyist POV / LED-sphere projects** | Not comparable. The published POV globes and spinning-LED builds are, almost universally, an Arduino sketch driving a fixed bitmap or a small palette animation, with 8-bit sRGB blending and no test suite. Holosphere has a rendering engine, a shader interface, a memory architecture, a simulator, and CI. The nearest peers are FastLED-ecosystem art frameworks (WLED, Pixelblaze), which are excellent at what they do but operate on 1-D/2-D strips with a scripting layer, not a templated 3-D engine with compile-time resolution specialization. |
| **Commercial / professional offerings** | Consumer "holographic fan" displays render pre-authored video onto a POV disc — the hardware is comparable, the software is a decoder. Professional spherical displays (Pufferfish, OmniGlobe, science-museum globes) are projector-based and drive equirectangular video from a PC; they solve a distribution problem, not a real-time-synthesis-on-an-MCU problem. On *engine* sophistication Holosphere sits closer to demoscene-grade real-time graphics or to an indie game engine's renderer than to any LED product. Where it falls short of commercial software is the surrounding product surface: no versioned releases, no plugin API, no user-facing documentation beyond the README. |
| **Academic research** | The individual mathematical components are drawn from the literature rather than novel: Hankin's construction and its computational treatment (Kaplan), Conway operator algebra, the Hopf fibration, Gray–Scott and Belousov–Zhabotinsky systems, Fibonacci sphere lattices, Ottosson's OKLab. The novelty is **integration and constraint** — running a Hankin pattern compiler and continuous Conway morphs at interactive rates inside 298 KiB on an MCU, with sphere-native anti-aliased rasterization, is a systems result that no paper I would expect to find has demonstrated. As a research artifact it would be a credible SIGGRAPH/Bridges *systems or art* contribution, not a theory contribution. |
| **State of the art in real-time rendering** | On a desktop GPU everything here is trivially achievable; that comparison is not the interesting one. Judged as *state of the art in what can be rendered on a 600 MHz Cortex-M7 with 298 KiB and no FPU-heavy budget to spare*, this is at or beyond the frontier I am aware of. The selective-O3 region system, the ITCM budgeting, and the per-effect cycle accounting are the kind of work that only happens when someone has genuinely run out of headroom and refused to lower the ambition. |

**Where it is weakest technically** is not in the engine but at its edges: four CI gates that pass on
empty input, a documentation layer that has drifted from the code it describes, and a published API
reference that silently omits the entire JS-facing surface. These are all cheap to fix and none of
them touch the architecture. They matter because this codebase's own standard is higher than that —
the project built `shard-coverage` precisely to catch the shard gap it currently has.

---

## 4. Artistic merit

The effect catalogue is not a grab-bag. Read together — Islamic star patterns on morphing polyhedra,
the Hopf fibration, spherical harmonics, Möbius transformations, reaction–diffusion, Voronoi cells,
gnomonic star fields — it is a coherent program: **mathematical structures that are natural on a
sphere, rendered on a sphere.** The sphere is not a canvas the work happens to use; it is the domain
the work is about. That is a genuine artistic thesis, and it is unusual for a technical artist to
sustain one over eleven years.

The craft decisions support it. Blending in linear light and interpolating palettes through OKLCH is
not a checkbox — it is the difference between colour that mixes the way a painter expects and colour
that goes muddy, and on a device where every pixel is an emissive point against black, that choice is
visible. The feedback style system (nineteen named presets over composable space/colour transforms)
is essentially a small instrument for a particular species of visual motion. The continuous Conway
morphing — where a truncated icosahedron becomes its dual becomes an ambo form, with the star pattern
tracking the topology change — is the strongest single idea in the catalogue, and it is one that
could not be pre-rendered.

Artistically the honest limitation is that the work lives in the tradition of generative/mathematical
art rather than extending it conceptually: it is superb execution of a well-established aesthetic
lineage (Islamic geometry, flow fields, cellular systems) on an unusual and well-chosen substrate.
The substrate and the real-time constraint are what make it distinctive, not the visual vocabulary.

---

## 5. Value to the community

**High as a reference; limited as a dependency.** The parts other people would benefit most from are:

- Sphere-native SDF rasterization with analytic scanline bounds — I am not aware of another open
  implementation.
- The multi-board POV synchronization protocol and its host-testable core; this is directly reusable
  by anyone building a segmented POV device, and the datasheet-style documentation makes it
  legible.
- The arena model with generation + rebind + rewind staleness stamping, as a pattern for embedded
  code that cannot use a heap.
- The dual-target bit-identity discipline, as a worked example of keeping a simulator honest.
- The performance regime — the ledgers of rejected optimizations are, unusually, more valuable to a
  reader than the accepted ones.

What limits the value is licensing and packaging. The split license (PolyForm Noncommercial for the
engine, all-rights-reserved for `effects/`) means the most transferable parts cannot be depended on
commercially, and the daydream `LICENSE` does not even carry the carve-out for the effects code its
shipped `.wasm` links. There is no release tagging, no `CONTRIBUTING`, and nothing extracted as a
standalone library. The live WebAssembly demo is the project's best community asset by a wide margin:
it makes the work immediately experienceable without hardware, which is exactly the right decision.

---

## 6. Honest limitations

- **Single-author bus factor.** The design coherence is a direct product of one mind holding the
  whole system; so is the risk. The comment discipline and the ledgers mitigate this better than most
  projects do, but the README drift shows what happens to the parts that only the author reads.
- **Documentation is the weakest deliverable.** A newcomer following the README will hit a deleted
  class, two wrong parameter lists, a wrong effect count, and an API reference missing the whole
  browser-facing bridge.
- **The gates need one afternoon of attention.** Four tests running in zero CI jobs and three gates
  that go green on empty input are not architectural problems, but they silently erode the
  verification story that is otherwise this project's strongest claim.
- **Scope has outgrown its scaffolding in two files.** `sdf.h` (4,589 lines) and `animation/mesh.h`
  (2,711 lines) are internally well-organized but are now large enough that a reader must trust the
  section comments rather than hold the file.

---

## 7. Verdict

> **An outstanding piece of independent engineering — top-percentile for a self-directed project, and
> credible against professional work in its specific niche.**

Judged as a *rendering engine under hard real-time and memory constraints*, this is excellent work by
professional standards: the abstractions are well-chosen, the invariants are enforced by the compiler
where possible and by traps where not, and the performance engineering is genuinely rigorous rather
than folkloric. Judged as a *product*, it is an unreleased one — the documentation and gate layers
have not kept pace with the engine, and the licensing forecloses most downstream reuse.

Judged as *art*, it is a sustained eleven-year exploration of mathematical structure on a spherical
substrate, executed with a colour pipeline and a motion vocabulary far more considered than the genre
usually receives — and the decision to publish a browser simulator that runs the identical code is
what turns a private object into something the public can actually encounter.

The gap between what this codebase *is* and what its documentation *claims* is the single highest-
leverage thing to close. The engineering underneath does not need defending.
