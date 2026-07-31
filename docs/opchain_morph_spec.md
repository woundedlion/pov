# OpChainMorph — Op-by-Op Morphing Generation of Solids — Design Spec

Status: **LANDED.** Implemented across merges 810db7c6, 1138ff8e and 661f1a78;
`effects/IslamicStars.h`, `core/animation/mesh.h` and `core/mesh/recipe.h` cite
its sections as the shipped design. §3 and §12 record constraints and dead ends
established against the renderer; read them before §10.

Goal: a solid should not appear, it should be **built**. Instead of fading a
finished polyhedron in with `Segue::TerminatorSweep`, the effect starts from the
recipe's base solid and animates **every Conway operator in its chain, one at a
time**, until the final solid stands complete — then fades out with
TerminatorSweep as it does today. Coloration must stay seamless across every op
boundary, to the standard `HankinSolids` already meets: no palette pop, no
position jump, no crease flash.

`IslamicStars` is the test case: its 24 recipes
(`core/mesh/solids.h:1096-1160`) are 2–6 op chains over 9 distinct seeds, so it
exercises composition depth, mid-chain `hankin`, mid-chain `relax`, and the two
ops that need genuinely new machinery (`kis`, `dual`). The design generalizes to
`catalan_registry` (13 entries, every one a single `.dual()` leg) and to
`simple_registry`, whose chains `HankinSolids` already sweeps.

**The load-bearing result** (§2.3): **13 of the 24 recipes contain no partition
op**, so they sit entirely outside §3's rasterizer constraints. Nine need
`kis`/`dual`, whose mechanism (§3.3) is unproven and measurement-gated; two more
are blocked by a clamp bound (§2.1). §10 phases the work along that split so the
feature ships on the tractable majority first. The 13 are not free — see §2.3 for
what they still need and §2.4 for the assumption eight of them rest on.

This spec is the sequel to `docs/conway_morph_spec.md` (LANDED 2026-07-16,
`06c95e93..190a91d8`). Read that first — it establishes the leg discipline,
the palette-provenance rules, and the bookend contract this design extends. Where
the two disagree, the landed ConwayMorph behaviour wins; this document only adds.

---

## 1. Facts this design builds on

1. **A solid is already a seed plus an op chain.** Only the five Platonic solids
   are explicit vertex data (`core/mesh/solids.h:56-170`). Every Archimedean,
   Catalan and Islamic entry is a `Solids::SolidBuilder` chain
   (`core/mesh/solids.h:203-364`), e.g.
   `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` (`solids.h:885-895`) =
   dodeca → hankin(35°) → ambo → hankin(62°) → ambo → relax(100) → hankin(42°).
   The chain exists; it is just *imperative* today, buried in a `FLASHMEM static
   PolyMesh f(Arena&, Arena&)` function body. §4 makes it data.

2. **Sweeping one op with constant topology is a solved problem.**
   `Animation::ConwayMorph` (`core/animation/mesh.h:177-703`) re-runs one Conway
   operator per frame at an interpolated parameter, compiles, and draws — 173 µs
   per frame for the 18-solid roster, 0.3 % of a frame, 0 spilled frames over two
   full tours. There is no vertex correspondence and no crossfade: topology is
   constant along a leg (`core/animation/mesh.h:435`), so vertex identity is implicit in the
   operator's deterministic emission order.

3. **Topology changes only at bookends, and only under clamp.** `clamp_param`
   (`core/animation/mesh.h:298-305`) keeps the sweep inside the topology-constant open interval:
   `T_EPS = 0.02` at the zero end, `T_EPS_AMBO = 0.005` below truncate's 0.5
   short-circuit (`conway_graph.h:113-124`). Newborn faces at the ε end are
   **zero-area and culled at draw time** by `SDF::Face` (|signed area| <
   1e-5·radius²) — see `docs/conway_morph_spec.md:48-52`. That cull is what makes
   a bookend swap invisible.

4. **Coloration is per-face, gradient-mapped, and provenance-tracked.** The
   overload IslamicStars actually uses on the hoisted path
   (`shading.h:158-167`, called at `IslamicStars.h:213`) is

   ```cpp
   float t = clamp(fragment_edge_dist(f) * gain, 0, 1);
   float cover = segue.fill(t, phase);        // t by NON-CONST ref: may remap or cull
   if (cover <= 0.0f) return Color4();
   Color4 c = segue.grade(palette.get(t), phase);
   c.alpha = cover * segue.opacity(phase);
   ```

   — an inset gradient keyed on the fragment's distance to its **own face's**
   boundary, normalized by face size (`fragment_edge_dist`, `shading.h:63`), then
   routed through the segue policy. The plain `palette.get(clamp(edge_dist·gain))`
   form is the *other* overload (`shading.h:101-111`) and is not what runs here;
   §3.3's note on policy dependence turns on the difference. Palette selection is per topology
   class through `MeshPaletteBank` (`core/color/palettes.h:111-175`, N = 5).
   Seamlessness across a leg comes from `build_palette_mapping`
   (`core/animation/mesh.h:595-698`): geometric provenance by departed-face centroid, with a
   checked bijection inside `PROVENANCE_TOL_SQ = 0.15²`.

5. **`hankin` is already a sweepable op.** `compile_hankin` bakes the
   angle-independent topology once; `update_hankin` repositions only the star
   points (`core/mesh/hankin.h:86-340`). At angle → 0 every star point collapses
   onto its corner via the `p_corner` fallback (`hankin.h:285-302`), so the
   hankin mesh at angle 0 **is** the base solid with midpoint vertices.

6. **`IslamicStars` cuts, it does not morph.** `spawn_shape`
   (`effects/IslamicStars.h:236-318`) walks the registry round-robin, generates
   the next solid whole into the carousel's back slot, flips front, and hides the
   discontinuity behind `Segue::TerminatorSweep`
   (`core/animation/mesh.h:884-926`). Per-shape choreography is
   `fade + still + burst_span + still + fade` (`IslamicStars.h:299`).

7. **The effect is rasterizer- and probe-bound, not geometry-bound.** Per
   `docs/islamicstars_next_levers.md` and the landed sector-walk work, ship peak
   is ~92 ms (8 fps) on the heaviest recipes, dominated by SDF probes on rejected
   pixels. Worst recipe is F = 1082. The per-frame Conway op is noise against
   that — see §8.

---

## 2. The op taxonomy — the central finding

A leg can be swept iff the operator has a **zero-area birth limit**: a parameter
value at which the faces the op *creates* have zero area and the faces it
*preserves* are unchanged. That is the only condition under which the opening
bookend swap changes no pixels, which is the only reason ConwayMorph is seamless.

Auditing `core/mesh/conway.h` against that condition splits the operator set
cleanly in two, and the split is the whole design:

| Op | Signature | Class | Leg treatment |
|---|---|---|---|
| `truncate(t)` | `conway.h:564` | **inflate** | sweep `T_EPS → t★` — but see §2.1, four recipes have unreachable `t★` |
| `expand(t)` | `conway.h:665` | **inflate** | sweep `T_EPS → t★` |
| `snub(t, twist)` | `conway.h:929` | **inflate** | sweep `T_EPS → t★`, twist lerped |
| `chamfer(t)` | `conway.h:736` | **inflate** | sweep `T_EPS → t★` |
| `hankin(θ)` | `hankin.h:388` | **inflate** | `compile_hankin` once, `update_hankin` per frame, sweep `θ_EPS → θ★` |
| `relax(n)` | `conway.h:813` | **continuous** | existing settle slerp (`core/animation/mesh.h:322-338, 422-427`) |
| `ambo` | `conway.h:458` | **inflate-equivalent** | truncate sweep to `0.5 − T_EPS_AMBO`, then clean swap to `ambo` |
| `bevel(t)` = `ta` | `conway.h:1082` | composite | ambo leg, then truncate leg |
| `kis` | `conway.h:405` | **partition** | §3.3 gated swap, no sweep (no renderable parameter, §3.2) |
| `zip` = `dk` | `conway.h:1069` | composite | kis leg, then dual leg |
| `needle` = `kd` | `conway.h:1057` | composite | dual leg, then kis leg |
| `gyro` = `ds` | `conway.h:1020` | composite | snub leg, then dual leg |
| `meta` = `kda` | `conway.h:1045` | composite | ambo leg, dual leg, kis leg |
| `dual` | `conway.h:355` | **reflection** | §3.3 gated swap + provenance mapping |

**Everything reduces to eight primitives: `truncate`, `expand`, `snub`,
`chamfer`, `hankin`, `relax`, `kis`, `dual`.**

Of those eight, **three are device-proven sweeps**: `ConwayGraph::MorphOp` is
exactly `{TRUNCATE, EXPAND, SNUB}` (`conway_graph.h:97`) and `run_op` has exactly
those three cases (`core/animation/mesh.h:495-506`). The rest each need work:

- `chamfer` is **deliberately unused** by the landed graph — "no simple-registry
  endpoint is a chamfered form" (`conway_graph.h`, echoed at
  `conway_morph_spec.md:423`). It has never been swept, never been through
  `clamp_param`, and has no `T_EPS` characterization. One Islamic recipe needs it
  (`truncatedIcosahedron_hk58_chamfer63`, `solids.h:683`) — and it runs on a
  hankin mesh, so see §2.4.
- `hankin` is swept by `HankinSolids` through an angle `Mutation`, not by
  `run_op`; the mechanism exists but has to be brought into the leg scheduler.
- `relax` is not a sweep at all — it is a settle-slerp bolted to a leg's end, and
  §2.2 shows it cannot always be bolted on.
- `kis` and `dual` need §3.3, which is unproven.

So the honest count is **three proven sweeps, one adaptable mechanism, four open
items**. Three of those four open items (`chamfer`, `hankin`, `relax`) fall
*inside* §2.3's pure-inflate set; only `kis`/`dual` are outside it.

The reduction is a **lowering step, not a storage format** — recipes keep the
authored op (`gyro`, not `snub, dual`) and the leg scheduler lowers on
construction. §4.2 explains why that distinction is load-bearing.

Consequences worth stating outright:

- **`ambo` is not a new problem, but not for the reason it looks.** `truncate` at
  exactly `0.5` short-circuits to `ambo` (`conway.h:567`), and face count is
  `F + V` on both sides (`conway.h:578`, `conway.h:466`). What *does* change
  across the boundary is **vertex count (2E → E)** and **every primary face's
  degree (2n → n)** — the `face_counts` array is not preserved. That, not
  emission order, is why a leg cannot cross 0.5 and why `T_EPS_AMBO` exists.
  (`conway_graph.h:111-112`'s comment is loose on this; do not inherit it.)
  ConwayMorph does already land legs on the cuboctahedron and icosidodecahedron
  through this boundary, so the ambo leg is a configuration of a device-proven
  one.

- **Decomposing composites is a feature, not a tax.** `icosahedron_kis_gyro`
  becomes icosa → **kis** → **snub** → **dual**: three legible legs where the
  recipe reads as two opaque ones. This is literally the requested "op-by-op".
  (All five decompositions and their leg orders were re-verified against
  `conway.h:1020, 1046, 1057, 1069, 1085`; composition is right-to-left and the
  table above reflects it.)

### 2.1 Four recipes have targets the clamp cannot reach

`T_EPS = 0.02`, `T_EPS_AMBO = 0.005` (`conway_graph.h:113, 119`). The shipping
registry violates both bounds:

| recipe | call | problem |
|---|---|---|
| `truncatedIcosidodecahedron_truncate50d_ambo_dual` (`solids.h:745-750`) | `.truncate(50·D2R)` = **0.8727** | far side of the ambo point |
| `truncatedIcosahedron_truncate50d_ambo_dual` (`solids.h:854-860`) | same | far side of the ambo point |
| `truncatedIcosahedron_ambo_relax_truncate001_hankin59` (`solids.h:918-924`) | `.truncate(0.01f)` | **below `T_EPS`** |
| `..._truncate001_hankin73` (`solids.h:934-940`) | `.truncate(0.01f)` | below `T_EPS` |

**Neither failure is a trap or a CI failure — both are silent, on-screen.**
`clamp_param` (`core/animation/mesh.h:298-307`) applies to `t_start` *and*
`t_end`, and it lives only in `ConwayMorph`: `build_recipe` and `entry.generate`
never see it, so §9.1's bitwise gate compares two unclamped chains and is
unaffected. The leg is what goes wrong:

- **`truncate001`:** both endpoints clamp to `T_EPS = 0.02`, so
  `t_start == t_end` and **the leg is a still image** for its whole duration.
  Its final frame sits at `t = 0.02` against a held solid built at `t = 0.01`.
- **`truncate50d`:** `t_end` clamps down to `0.5 − T_EPS_AMBO = 0.495`, and
  `tr.topo` is classified from `run_op` at that *clamped* arrival
  (`core/animation/mesh.h:315-343`), so face counts agree every frame and there
  is **no trap**. The leg sweeps to a cuboctahedron-like form and then clean-swaps
  to the self-intersecting `truncate(0.873)` mesh the recipe actually specifies —
  a full-screen geometry pop. (`conway.h:555-560` documents that `t > 0.5`
  deliberately crosses the cut points to produce those self-intersecting faces.)

**Treatment:** all four are non-morphable on their truncate leg and both failures
are §9.5 seam items, not gate failures — which means CI will *not* catch them and
they must be checked by eye. Either give the four a null recipe, or add a
`t > 0.5` leg kind seeded from `ambo(seed)` with its own clamp on the far side of
the ambo point. The first is correct for Phase 2; the second is speculative and
should not be scheduled without a reason.

**Fifth case, benign but must be handled:** `bevel(0.5)` in
`truncatedIcosidodecahedron_bevel5_relax_hk77` (`solids.h:964`) lowers to a
truncate leg whose target is *exactly* `0.5` — the one value a leg may never land
on. `expand_to_primitives` (§4.2) must special-case `t == 0.5f` in a bevel to emit
`AMBO` instead of a truncate step.

### 2.2 `relax` cannot always fold into the preceding leg

The settle mechanism slerps per-vertex between the swept mesh and
`relax(arrival)` and therefore requires vertex-count identity —
`HS_CHECK(swept.vertices.size() == tr.relaxed.size())` (`core/animation/mesh.h:423`, with the
construction check at `:336`).

The Islamic registry has **nine** `relax` calls (`solids.h:697, 872, 892, 921,
937, 965, 978, 993, 1009`), and **six** of them follow an ambo-class leg —
`:697, 892, 921, 937`, plus `:965` (`bevel(0.5)` lowers to ambo, §2.1) and
`:993`. One follows `snub` (`:872`); two follow the truncate leg of a
`bevel(0.2)` / `bevel(0.33)` (`:978, 1009`).

An ambo leg is a truncate sweep whose swept mesh has **2E** vertices
(`conway.h:576`), while `relax(ambo(seed))` has **E** (`conway.h:465`). The slerp
cannot run during that sweep. **Relax after ambo needs its own leg**, starting
after the clean swap to `ambo`, with its own frames. Only relax following a
vertex-count-preserving leg folds for free — three of the nine.

The Phase-2 (pure-inflate) standalone-relax sites are `solids.h:892, 965, 993`.

### 2.3 The inflate-only subset — the most useful fact in this document

Classifying the 24 Islamic recipes by whether they contain any partition op
(`kis`, `dual`, or a composite containing one — `gyro`, `meta`, `needle`, `zip`):

- **Pure-inflate (13):** `dodecahedron_hk62_ambo_hk62`,
  `truncatedIcosahedron_hk58_chamfer63`, `dodecahedron_ambo_bevel33_relax_hk66`,
  `truncatedIcosahedron_ambo_relax_truncate33_hk64`,
  `truncatedIcosidodecahedron_bevel5_relax_hk77`,
  `icosahedron_ambo_truncate033_hankin59`,
  `dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`, `octahedron_hk17_ambo_hk73`,
  `octahedron_hk34_ambo_hk72`, `rhombicuboctahedron_hk63_ambo_hk63`,
  `truncatedIcosahedron_hk54_ambo_hk72`, `dodecahedron_hk54_ambo_hk72`,
  `icosahedron_snub_relax_truncate033_hankin62`.
- **Pure-inflate but blocked by §2.1 (2):** the two `truncate001` recipes.
- **Needs §3.3 (9):** the rest — every `_dual`, `_gyro`, `_kis`, `_needle` entry.

**Within those 9, `kis` and `dual` are not equal partners.** Lowering every
composite and counting primitives:

| primitive | Islamic recipes | also unlocks |
|---|---|---|
| `dual` | **9 of 9** | all 13 Catalans (each is a single `.dual()` leg) |
| `kis` | **3** — `..._hk54_needle`, `icosahedron_kis_gyro`, `truncatedOctahedron_gyro_kis_hk17` | nothing |

All three `kis` recipes also contain `dual`, so **`kis` unlocks no recipe on its
own**. `dual` alone unlocks 6 Islamic entries plus all 13 Catalans — 19 solids;
`kis` adds 3 more, and only on top of `dual`. §10 Phase 3 is split accordingly.

Note `bevel` is *not* a partition op: `bevel = ta` decomposes to an ambo leg plus
a truncate leg, both inflate.

**This reframes the whole project.** Thirteen of twenty-four recipes — the
majority, and including the visually richest multi-hankin chains — need **no
partition machinery**, so the unproven §3.3 gate is off their critical path. It
is required for nine.

Be precise about what the 13 still need, because it is not nothing: all 13 need
the `hankin` leg kind brought into the scheduler, five need §2.2's standalone
relax leg, one needs a `chamfer` characterization, and eight need §2.4 resolved.
Those are Phase-2 scope. What they do *not* need is `kis`, `dual`, or anything
downstream of §3's rasterizer constraints.

### 2.4 Eight of the 13 apply their Conway op to a *hankin* mesh

The largest unexamined assumption in this document. These recipes run `ambo` or
`chamfer` on the output of a `hankin`, not on a simple solid:
`dodecahedron_hk62_ambo_hk62`, `truncatedIcosahedron_hk58_chamfer63`,
`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`, `octahedron_hk17_ambo_hk73`,
`octahedron_hk34_ambo_hk72`, `rhombicuboctahedron_hk63_ambo_hk63`,
`truncatedIcosahedron_hk54_ambo_hk72`, `dodecahedron_hk54_ambo_hk72`.

The whole §2 taxonomy — zero-area birth limits, `T_EPS` sizing, "newborn faces
are zero-area and culled at draw time" — is characterized against **simple
solids with ≤ 92 faces** (`conway_graph.h:114-118`). A hankin mesh has hundreds
of rosette and strap faces at wildly varying scales, and three things follow that
nobody has checked:

1. **`T_EPS = 0.02` births are not sub-pixel relative to a *small* rosette
   face.** The cull threshold is a fraction of each face's own radius
   (`SDF::Face`, |signed area| < 1e-5·radius²), but `T_EPS` is a global constant
   tuned at ≤ 92 faces. On the small faces the opening bookend may be visible.
2. **`MeshOps::compile` drops degenerate faces size-dependently**, so
   `compiled.face_counts.size()` **can** move mid-leg — which *is* the real
   `HS_CHECK` trap at `core/animation/mesh.h:435`. §2.1's cases turned out not to
   trap; this one plausibly does.
3. **`require_closed_manifold` must hold for the hankin output under sweep.** It
   holds today at the chains' fixed parameters, since they ship — but nothing
   pins it across a swept `t`.

Both of §10 Phase 1's recipes are ambo-on-hankin, so this is the assumption the
plan rests on first. **§10 Phase 0 exists to answer it**: a host test stepping an
ambo leg on a hankin seed frame by frame, asserting constant compiled face count
and manifold-ness, before any animation work. If it fails, `T_EPS` needs a
per-mesh derivation — scaled from the smallest face radius rather than a global
constant — and the whole pure-inflate set inherits the fix.

- **`kis` has no zero-area limit and never will.** It *partitions* an n-gon into
  n triangles that always sum to the parent's area. Lowering the apex into the
  face plane makes the fan coplanar with the parent, which is geometrically
  identity — but **not** shading identity, because each child restarts the
  `fragment_edge_dist` gradient at its own boundary. A coplanar fan pops as a
  starburst of n bright seams. Same for `dual`, plus relocation on top. §3 is the
  answer.

---

## 3. Partition ops: rasterizer constraints and what remains

Two constraints bound every option for `kis` and `dual`. Both are properties of
the **renderer**, not of the mesh or animation layers where operator design
naturally happens, and neither is visible from those layers.

### 3.1 Constraint: coplanar partitions always seam

`Scan::Mesh::draw` rasterizes one `SDF::Face` at a time and antialiases every
face against **its own** boundary (`core/render/scan.h:1238-1241`):

```cpp
float alpha = 1.0f;
if (d > -pixel_width) {                    // within a pixel of ANY edge of THIS face
  float t_aa = 0.5f - d / (2.0f * pixel_width);
  alpha = quintic_kernel(clamp(t_aa, 0, 1));
}
```

and `plot` composites **over**, not additively (`core/render/filter.h:151` ->
`lerp16`, `core/color/color.h:415`). On a pixel along an interior spoke of a
coplanar fan both children report `d ~ 0` -> `alpha = 0.5` each
(`quintic_kernel(0.5) = 0.5` exactly, `core/math/3dmath.h:72-75`), and
over-blending onto black gives `1 - (1-0.5)^2 = 0.75*c` where the parent gave
`1.0*c`. At `CANVAS_W = 288` / `CANVAS_H = 144`, `pixel_width = 2*PI/288` equals
the row spacing, so the band is isotropic.

**Every interior edge of a partition renders as a ~2 px band up to 25 % dark, at
any `gain`, including 0.** Near a parent vertex it is worse: the adjacent fan
triangles have different edge geometry from the parent, so their coverages do not
sum to the parent's at all.

**MEASURED-DEAD: no shading-domain trick removes this.** Collapsing the shading
`gain` to zero removes the *gradient* discontinuity and leaves the *coverage*
discontinuity untouched — the two are independent. The only exact fix is a
coverage-accumulating plot path for mesh faces, a global rasterizer change
touching every effect, out of scope here.

Same defect class as `conway_morph_spec.md:47-51`, overrule 5: zero-area straps
still rasterized as ~1 px lines along every edge.

### 3.2 Constraint: radial vertex motion moves shading, not silhouette

`SDF::Face::setup_frame_and_polygon` projects gnomonically
(`core/render/sdf.h:2489-2496`) and derives edge normals as
`cross(v_i, v_{i+1})` normalized (`:2513-2517`). Both are invariant under
positive scaling of a vertex, and great circles map to straight lines under any
gnomonic center, so **the drawn silhouette does not move**: inside/outside
classification is a function of vertex *directions* only.

What is **not** invariant is everything that gives the face its shading. The
frame is built from `center = sum of vertices` **un-normalized**, then normalized
(`sdf.h:2472-2478`), so *non-uniform* radial motion rotates `basis_u/v/w`. Every
magnitude downstream of that basis moves with it:

- `Face::distance` projects the sample point into the same basis and reports
  `res.dist = atan(plane_dist) - thickness` (`sdf.h:3221-3223, 3296-3297`) — a
  2-D Euclidean distance in the gnomonic plane, hence basis-dependent;
- `radius` / `max_dist_sq` (`sdf.h:2497-2499`) and the gnomonic inradius `size`
  (`compute_inradius`) change;
- `fragment_edge_dist = -v1/size` (`shading.h:63-65`) therefore rescales the
  whole shading field — the defect class that killed the pyritohedron leg
  (`conway_morph_spec.md:523-525`);
- IslamicStars' per-face segue offset is also a raw vertex **sum**
  (`IslamicStars.h:194-199`), so face phases shift too.

**Consequence for `kis`:** an overload sweeping the apex radius from the face
plane to the sphere produces **a gradient shimmer with no motion** — the fan's
outline never changes. Useless as an animation. `compute_phi_extent` additionally
takes `min`/`max` of the raw `vertices[idx].y` (`sdf.h:2436-2440`), so a sub-unit
apex shrinks the scanline band without shrinking the drawn region and can drop
rows.

**MEASURED-DEAD: apex-radius sweeps, and radial "grow it outward" motion
generally, for `SDF::Face`-rendered geometry.** Scope the rule precisely:

- *Uniform* radial scaling of a face is fully invisible.
- *Per-vertex* radial motion changes shading but not silhouette.
- The existing inflate sweeps are unaffected: they move vertices **tangentially**
  (and `normalize` puts everything back on the sphere, changing directions),
  which is exactly why they read as motion.
- IslamicStars' ripple is **not** a counterexample: `ripple_transform` is a
  rotation about `cross(center, v)` (`core/engine/transformers.h:459-499`) —
  tangential and magnitude-preserving.
- Outside mesh rasterization (plot, raymarch, `FieldTransformer` effects) radius
  matters entirely. This rule is about `SDF::Face` only.

### 3.3 What remains: the gated swap

Both ops are partitions with no continuous parameter and an unavoidable coverage
seam. They reduce to the same three-step shape, with **no sweep segment**:

```
GATED SWAP (per partition leg)

  1. CLOSE   F_gate frames: gain 1 -> gain_min on affected parent faces.
  2. SWAP    1 frame: replace the mesh with op(mesh); children inherit
             the parent's colour by geometric provenance.
  3. OPEN    F_gate frames: gain_min -> 1 on the children toward their
             leg-target palettes.
```

The gate is **not** exact. What it buys is bounded but real: it reduces the
swap's visible delta from *gradient restart plus coverage seam* to *coverage seam
only*. Whether that is good enough is empirical — §9.6 measures it, and §10 makes
the measurement a gate before any of this is scheduled.

Two costs to price in:

- **`gain_min` must not be 0.** At `gain = 0` every fragment resolves to
  `palette.get(0)`, and the five bank palettes are IQ cosine ramps whose `t = 0`
  values differ wildly in luminance (`core/color/palettes.h:27-90`): EMBERS
  (0.058, 0.035, 0.186) near-black, LAVENDER_LAKE (0.787, 0.112, 0.832) magenta,
  BRIGHT_SUNRISE (1.0, 0.827, 0.341) near-white. The mesh would collapse to five
  flat, high-contrast patches with all Hankin star structure gone, twice per
  partition leg — the dissolve-to-flat the landed design deliberately deleted
  (`conway_morph_spec.md:9-11`). Tune `gain_min` as the shallowest gate that
  hides the seam.
- **The gate is segue-policy-dependent.** The shading overload in play
  (`shading.h:158-167`) passes `t` to `segue.fill(t, phase)` by **non-const
  reference**, which may remap or cull it (`core/animation/mesh.h:798`).
  `TerminatorSweep` inherits the identity `fill` so a low `gain` survives, but
  `IrisBloom` and `Lace` (`:845, :866`) would cull small `t` outright. Any effect
  adopting the gate must re-check its segue.

**The two ops are different visual events, and §9.6 must measure them
separately.** Neither can fade in — coverage AA is geometric, so a new edge
either exists in a frame or does not:

- **`kis` is additive.** It keeps every original vertex and edge and adds an apex
  plus `n` spokes per face (`conway.h:405-449`). The swap frame reveals `n·F` new
  lines over an otherwise unchanged edge network.
- **`dual` is wholesale.** Its vertices are the source's face centroids and none
  of the source's edges survive (`conway.h:355-396`). The swap frame replaces the
  entire edge network at once.

§3.1's coplanar-partition analysis applies directly to `kis`; for `dual` the
question is not seams but whether a whole-network replacement under flat colour
reads as a transformation or as a cut.

**`dual`'s provenance mapping.** Each dual face `d` corresponds to source vertex
`v` and takes the palette of the orbit face whose centroid is nearest `v`,
matching the nearest-centroid discipline in `build_palette_mapping`
(`core/animation/mesh.h:640-666`). The orbit is already walked by
`emit_vertex_orbit_faces<'P'>` (`conway.h:150`), so the mapping is a by-product
of the op rather than a second pass. This gives colour-*locality* across the swap
— a dual face opens in the colour already painted where it lands. Pixel-identity
is not available (§3.1).

### 3.4 Named fallback

If §9.6 rejects the gated swap, the fallback is a short local
`Segue::Crossfade` (`core/animation/mesh.h:815-832`, the only overlapping policy
in the codebase) across the swap frame. An opacity blend of two complete meshes
has no coverage-seam problem, at the cost of drawing two meshes for the window's
duration on an effect already at 8 fps for its heaviest recipes. Six frames of
double-draw at the small end of a chain is affordable; at the large end it is
not, which likely means a window length scaling inversely with face count.

It also cuts against `conway_morph_spec.md:9-11` ("no crossfade or dissolve
anywhere in the design"). That was a deliberate aesthetic call on the landed
work; re-opening it for two ops needs the owner's agreement, not a spec's
assertion.

### 3.5 Gate cost

`F_gate` initial proposal **6 frames** each side (matching `SHAPE_FRAMES`,
`HankinSolids.h`). A partition leg costs `6 + 1 + 6 = 13` frames — there is no
sweep segment (§3.2). A recipe with three partition ops pays 39 frames.

---

## 4. Recipes as data

`Solids::Entry::generate` (`solids.h:1025-1030`) is an opaque function pointer.
The build needs to step the chain, so the chain must be inspectable. So does the
solids tool (§4.4), which today can only *emit* chains, never read one back.

### 4.1 The table

```cpp
// core/mesh/recipe.h  (new)
namespace Solids {

// Authored ops, including the composites. expand_to_primitives() lowers these
// to the eight of section 2 for the leg scheduler; the tool and codegen consume
// them unlowered.
enum class Op : uint8_t { TRUNCATE, EXPAND, SNUB, CHAMFER, HANKIN, RELAX,
                          KIS, DUAL, AMBO,                 // primitives
                          BEVEL, GYRO, META, NEEDLE, ZIP };  // composites

struct OpStep {
  Op    op;
  float param = 0.0f;   // t / angle(rad); iteration count for RELAX
  float twist = 0.0f;   // SNUB only
};

struct Recipe {
  uint8_t       seed;   // index into simple_registry (Platonic | Archimedean)
  const OpStep *steps;
  uint8_t       count;
};

} // namespace Solids
```

Rules:

- **`Recipe` is `constexpr` table data**, one per morphable registry entry,
  living beside the registry it mirrors.
- **The seed is a `simple_registry` index**, not a new enum. That registry is
  `[Platonic | Archimedean]` and its order is already load-bearing and pinned by
  18 `static_assert`s in `ConwayGraph` (`conway_graph.h:46-90`). It also yields a
  **name** through `get_entry(i).name`, which §4.3 needs. Recipes are therefore
  *shallow*: a Catalan is `seed = snubDodecahedron, steps = [DUAL]`, not a
  flattening down to a Platonic. That is exactly the leg decomposition §2 wants
  and exactly the base+ops model the tool already has.
- **`Solids::Entry` gains `const Recipe *recipe`**, nullable. A null recipe means
  "not morphable" and the effect falls back to today's whole-generate + full
  TerminatorSweep. This is the incremental-adoption lever: land the machinery
  with two recipes populated, extend the table over time.
- **`build_recipe(const Recipe&, Arena&, Arena&) -> PolyMesh`** replays the chain
  through `SolidBuilder`. It is not used at runtime by the effect; it exists so
  §9.1 can gate `build_recipe(entry.recipe) == entry.generate()` **bitwise**.
  That equivalence gate is the anchor of the whole feature: it is what stops the
  recipe table and the shipping geometry from silently diverging.

Deliberately **not** doing: deriving recipes from the generator functions, or
replacing the generators with the recipes. The generators stay the source of
truth for shipping geometry; the recipes are a parallel declaration that CI
proves equal. Two independent statements that must agree catch authoring errors
that one clever statement cannot.

Equally deliberately **not** doing: parsing the registry *name*. The names encode
the chain (`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42`) but lossily —
`hk35` is a rounded degree, `truncate50d` is `50·D2R = 0.873` used as a raw `t`,
`relax` omits its iteration count, and `chamfer63` means `0.63`. A name parser
would be a second, weaker source of truth competing with the recipe table. The
name stays a display label.

### 4.2 Authored ops, not lowered ops

`OpStep::op` stores what the author **wrote** — `GYRO`, not `SNUB, DUAL`. The
leg scheduler lowers composites through a pure

```cpp
size_t expand_to_primitives(const Recipe&, OpStep *out, size_t cap);
```

per §2's table. Two reasons the lowering must live at the scheduler and not at
authoring time:

1. **Codegen round-trips.** The tool's `generateFuncAndRecipe`
   (`solid_codegen.js:102-183`) emits `.gyro()` and names the solid `_gyro`. If
   the recipe stored `SNUB, DUAL`, loading `icosahedron_kis_gyro` and
   re-exporting it would produce `icosahedron_kis_snub_dual` — a different name
   for identical geometry, and registry churn on every round trip.
2. **Authored intent is the reviewable artifact.** A recipe that reads like the
   generator it mirrors is checkable by eye against `solids.h`; a pre-lowered one
   is not, which weakens §9.1's gate to a machine-only check.

The morph runtime still never decomposes anything at draw time — lowering is
once per shape, in the cold constructor.

### 4.3 Crossing the WASM boundary

One new embind class function beside the existing `getRegistry`
(`targets/wasm/wasm.cpp:1064-1076`, bound at `:1284-1302`):

```cpp
// MeshOps.getRecipe(name) -> {seed: string, ops: [{op: string, param, twist}]}
//                         -> null when the entry has no recipe
static val getRecipe(const std::string &name);
```

Shape notes:

- **`seed` is a name string**, resolved through `get_entry(recipe->seed).name`,
  because the tool's whole model is "load a base by name, then apply ops"
  (`solids.html:1844`, `applyOp` `:1652-1663`). Handing it an index would force
  the tool to keep its own index→name table.
- **`op` is the lowercase op string** already used by the tool's `data-op`
  attributes and `OP_DEFS` (`solids.html:465-504`, `:720-734`). It does **not**
  drop straight into `state.ops` — §4.4 lists the translation required.
- **Params cross in engine-native units** (radians for `hankin`, raw `t` for
  `truncate`), matching `MeshOps`' own bindings. The tool converts on the way in
  (§4.4); putting the conversion in the binding would make the payload disagree
  with every other engine surface.
- Follows `fromSolidName`'s precedent (`wasm.cpp:792-803`): **return null for an
  unknown or recipe-less name**, never abort. The engine is `-fno-exceptions` and
  traps fatally; the tool's own comment at `solids.html:1826-1838` says aborts
  are not catchable.
- Pure table read — no arena, no mesh, no `clearToolingMemory()` pairing.

`scripts/wasm_smoke.mjs:277-375` already gates this surface against signature
drift and gets one more case.

### 4.4 The solids tool reconstructs the chain on load

Today `tools/solids.html` loads a solid as an **opaque base**: `fromSolidName(name)`
returns a finished mesh and the editable op list starts empty. The 14-button op
palette, the draggable chain, the parameter sliders and the C++ codegen
(`solid_codegen.js:102-183`) all exist — but only for chains you build from
scratch. You cannot load `dodecahedron_hk62_ambo_hk62`, nudge one angle, and
re-export; you have to retype the chain from the name.

**With `getRecipe`, loading a solid populates the editor:**

```
load("dodecahedron_hk62_ambo_hk62")
  -> base = "dodecahedron"
     ops  = [ hankin(1.0821 rad), ambo, hankin(1.0821 rad) ]
```

which the tool's existing render path (`fromSolidName(base)` then `applyOp` per
step) reproduces — **once the payload is translated.** The translation is not
optional and is the fiddliest part of this section:

- **`state.ops` entries are `{op, params:{…}}` with op-specific key names** —
  `truncate:{t}`, `hankin:{angle}`, `relax:{iter}`, `snub:{t,twist}`,
  `bevel:{t}` (`OP_DEFS`, `solids.html:720-734`; `applyOp`, `:1652-1663`). The
  flat `{op, param, twist}` payload needs a per-op key map.
- **`hankin.angle` is in DEGREES in the tool** — `applyOp` does
  `mesh.hankin(o.params.angle * (Math.PI / 180))` (`solids.html:1658`). A payload
  carrying radians must be converted on load, or `hk62` loads as 1.08°.
- **Slider domains reject shipping values.** `truncate.t` and `bevel.t` are both
  `[0.01, 0.5]` (`solids.html:725, 733`), so `truncate(0.8727)` — §2.1's
  `truncate50d` pair — is unrepresentable in the editor. Those recipes load with
  a value the UI cannot display or restore. Either widen the domains to match the
  engine's `[0, 1]` or mark the op read-only when the value is out of range;
  silently clamping would corrupt the chain on any subsequent edit.

Behavioural changes:

- `state.base` becomes the **recipe seed**, not the entry name; `state.ops` is
  seeded from `recipe.ops`. Entries with a null recipe keep today's behaviour
  (base = entry name, ops = []). **This is the default, not a toggle** (owner
  decision, §11) — clearing the ops recovers the opaque view, so nothing is lost
  and no second load mode has to be maintained.
- **Mark recipe-less entries in the UI** — a badge on the thumbnail strip
  (`generateThumbnails`, `solids.html:933-1087`). This turns the tool into the
  live **coverage dashboard** for recipe adoption across §10's phases, which is
  worth more than the badge costs.
- **`state.base` is also the thumbnail identity key**, not just a build input:
  it drives `.active` / `aria-pressed` on the thumbnail strip
  (`solids.html:1027-1028`, `:1059`), and titles saved solids via
  `formatSolidName(state.base)` (`:1120`, `:1301-1302`). Seeding it with the
  *seed* means clicking `dodecahedron_hk62_ambo_hk62` highlights the
  **dodecahedron** thumbnail and saves as "Dodecahedron". Fix by tracking the
  loaded entry name separately from the build base and keying the strip and the
  save title on the former. This is the one behavioural change that is not small.
- **Keep the validation gate on reconstructed chains.** It is tempting to skip
  the `withValidator` / `chainIsValid` round trip (`solids.html:1665-1676`) on the
  grounds that §9.1 already proved the chain valid. It did not prove the same
  thing: §9.1 proves the chain builds under the **native/device** arenas, while
  the gate exists because a trap **permanently wedges the WASM module**, whose
  `tooling_arena` is a different budget entirely.
- **Thumbnails keep calling `fromSolidName` directly.** They want the finished
  mesh, not the chain; replaying 55 chains offscreen would be strictly slower for
  an identical result.
- Optionally, an **"expand seed"** affordance: recipes are shallow (§4.1), so
  loading a Catalan gives `seed = snubDodecahedron, ops = [dual]`. One click
  substitutes the seed's own recipe to walk further down. Nice-to-have, not
  required.

Two things this buys beyond convenience:

1. **It is the cheapest possible validator for the recipe table.** Every solid in
   the roster becomes a visual A/B: reconstructed chain versus
   `fromSolidName`. Authoring errors that §9.1's bitwise gate catches as a red CI
   line are caught here as a visibly wrong shape, with the offending op sitting
   in an editable list. This lands in Phase 1 and de-risks every later phase.
2. **It closes the authoring loop.** Load an existing solid → tweak → save →
   codegen → paste into `solids.h`. Today that loop only runs forward.

---

## 5. The leg scheduler

### 5.1 Shape

Do not build a monolithic multi-op animation. Chain single-leg animations on the
`Timeline`, exactly as `HankinSolids` chains morph legs through
`start_morph_cycle` / `finish_morph_cycle` (`HankinSolids.h:690-842`).

```cpp
// core/animation/mesh.h
class OpLeg : public AnimationBase<OpLeg> {  // generalizes ConwayMorph
  // Three leg kinds, dispatched once at construction:
  //   CONWAY_SWEEP  - run_op(op, seed, t(frame))          [existing path]
  //   HANKIN_SWEEP  - update_hankin(compiled, out, θ(frame))
  //   GATED_SWAP    - gain gate around kis / dual, no sweep [§3.3]
};
```

Leg construction runs `expand_to_primitives` (§4.2) once per shape, so the
scheduler only ever sees the eight primitives and the draw path never
decomposes anything.

**`HankinSolids` does not adopt `OpLeg` in this work** (owner decision, §11).
The refactor below must be behaviour-preserving for its existing consumer: after
it lands, `HankinSolids` still drives graph-edge legs exactly as it does today,
and its device profile is unchanged. Converging the two effects onto the chain
scheduler is a separate, later change. Concretely this means the graph-edge
constructor stays a supported entry point alongside the recipe-step one, rather
than being reworked into it.

`ConwayMorph` (`core/animation/mesh.h:177`) is the CONWAY_SWEEP kind with the graph-edge
constructor. **Refactor it into `OpLeg`'s first kind rather than forking it** —
the palette machinery (`PaletteHandoff`, `BookendClasses`, `Landing`,
`build_palette_mapping`, `bake_palette_blend` pre-blend, `MAX_BLEND_PAIRS`) is
kind-agnostic and must not be duplicated. `run_op` (`core/animation/mesh.h:495-506`) grows the
missing primitives.

### 5.2 State carried between legs

Leg `k`'s seed is leg `k−1`'s **completed output**, held as a persistent
`PolyMesh`. Per-frame cost is therefore **one op**, never a chain replay —
identical to ConwayMorph today. What crosses the boundary:

| Carried | Type | Source |
|---|---|---|
| seed mesh | `PolyMesh` (persistent, cloned) | previous leg's final `run_op` |
| departed provenance | `ConwayMorph::PaletteHandoff` | previous leg's `Landing` + centroids |
| arrival classification | `ConwayMorph::BookendClasses` | `classify_faces_by_topology` on the next leg's clean endpoint |
| slot → palette | `std::array<uint8_t, 5>` | previous `Landing::to_palette` |
| gate state | `float gain` per face | gated swap, §3.3 |

`Landing` (`core/animation/mesh.h:244-251`) already carries `topology`, `faces`,
`primary_faces`, `to_palette` and is the existing leg→effect handoff. It needs
one addition: the **final swept `PolyMesh`**, or a cheap way to reconstruct it
(op + final parameter + seed), so leg `k+1` can seed from it without a
whole-chain replay.

### 5.3 Palette targets across a chain

`HankinSolids` re-shuffles palette slots once per solid. For a chain build the
shuffle happens **once per shape, before the first leg** — the whole build is one
shape. Newborn classes inherit their first face's nearest-departed palette
(`core/animation/mesh.h:657-666`), which keeps `T_EPS`-wide births from popping
in target-coloured.

**Convergence is spread across the whole build** (owner decision, §11): the shape
reaches its final slot assignment on the *last* leg, not the first, so colour and
geometry arrive together. Each leg carries a share of the total colour distance
weighted by its frame count.

**Spreading is specified but not implemented.** `OpLeg` carries no per-leg blend
range; every leg applies its own `blend_fn` weight directly against the bank
endpoints, so each leg converges over its own span. The rest of this section is
the design a chained implementation would follow.

**Mechanics.** `ramp_from` / `ramp_to` are `uint8_t` **bank indices**, consumed
as `tr.bank->entries[tr.ramp_from[r]]` (`core/animation/mesh.h:445-451,
486-487`). `bake_palette_blend` (`core/color/composition.h:1042-1051`) resolves
to `from.get(t).lerp(to.get(t), w)` -> `Color4::lerp` -> `Pixel::lerp16`
(`core/color/color.h:153-162`): **a linear lerp on 16-bit linear-light
channels**, no perceptual space anywhere. So

```
lerp(lerp(A, B, c0), lerp(A, B, c1), w)  ==  lerp(A, B, c0 + w*(c1 - c0))
```

and leg `k` spanning cumulative colour fractions `c(k-1) -> c(k)` is expressed by
feeding a remapped weight against the **unchanged bank endpoints**. No
intermediate palettes are ever materialized and no new storage is needed.

Two costs that are *not* free, both easy to miss:

- **The aliasing fast path stops firing.** `bake_palette_blend` aliases the
  endpoint LUT with no allocation and no bake when `w <= 0` or `w >= 1`
  (`composition.h:1037-1040`). Under spreading the remapped weight is `c(k-1)` /
  `c(k)`, which are never 0 or 1 on any intermediate leg — so every frame of the
  20 % plateaus at both ends of every leg now bakes `num_ramps x 256` `Color4`
  entries into `scratch_arena_b` instead of aliasing. With `MAX_BLEND_PAIRS = 8`
  and ~3.5 KB of headroom in `b` (§8.2), this lands on exactly the frames the
  arena analysis assumed were cheapest. §9.3 must exercise it.
  (The `ramp_from[r] == ramp_to[r]` fast path at `core/animation/mesh.h:448` is a separate,
  index-based test and is unaffected by re-basing.)
- **Bookend swaps stay exact, but only to a few LSB.** The same `w'` on both
  sides of a swap yields the same LUT, so the swap itself is exact. The nested
  lerp, however, quantizes: the weight goes through `frac_to_q16` (uint16),
  channels are uint16, and `lerp16`'s shift-divide is documented as within 1 LSB
  of `round(x/65535)`. Rebasing drifts a few LSB. That sits inside
  `conway_morph_spec.md:50-51`'s +/-2 LSB bookend tolerance, but §9.5's
  "inflate-op bookends <= 0.5 %" threshold must be stated as tolerating it rather
  than as exactness.

**The visible cost is a stutter.** `blend_weight` is exactly 0 for `p <= 0.2` and
exactly 1 for `p >= 0.8` (`core/animation/mesh.h:560-568`) — plateaus that exist
to make each leg's bookends exact swaps. Under spreading, colour is therefore
**frozen for the first and last 20 % of every leg**: it advances in a
move/freeze/move/freeze rhythm, 40 % dead per leg, with the period varying by leg
length rather than staying constant.

Mitigation if that reads badly: give colour its own clock over the whole build,
sampled per frame independent of leg boundaries, and keep `blend_weight`'s
plateaus only for the *geometric* bookends. That decouples the two and is a small
change, but it gives up exact colour swaps at leg boundaries — measure before
adopting. Eyeball-gate at §9.10.

---

## 6. IslamicStars integration

### 6.1 Choreography

Today (`IslamicStars.h:299`):

```
fade_in(TerminatorSweep) → still → ripple burst → still → fade_out(TerminatorSweep)
```

Proposed:

```
fade_in(seed solid, TerminatorSweep)      fade frames
build: leg 1 … leg n                      Σ leg frames      ← new
still                                     still frames
ripple burst                              burst_span
still                                     still frames
fade_out(final solid, TerminatorSweep)    fade frames
```

**The seed sweeps in on TerminatorSweep** (owner decision, §11). It has to arrive
from somewhere, and the framing is explicit: *start with the base solid*, build,
*then* fade out. The seed is small (a Platonic or low-F Archimedean solid), so its
sweep-in is the cheapest phase in the whole cycle — and the visual reading is
right: the sweep delivers a bare seed, the build does the work, the sweep takes
away a finished solid.

Consequence for `retarget`: `TerminatorSweep::retarget(v)` rolls a fresh sweep
axis and `fade_seed` per transition (`IslamicStars.h:279-280`, `core/animation/mesh.h:~905`).
It must keep firing **once per shape, before the seed sweeps in** — not per leg.
The build phase sits on the phase-1 plateau where the axis is unused, so a
mid-build retarget would be invisible on the build and then wrong on the way out.

`MeshCarousel::schedule_segue` (`core/animation/mesh.h:1126-1129`) assumes a symmetric
`(duration, window)` sprite. The build phase sits on the phase-1 plateau where
`TerminatorSweep::opacity(1) = 1` and `fill` is identity, so **the simplest
integration is to leave the carousel contract alone and lengthen `duration` by
the build span**, drawing build frames through the same `draw_fn` at phase 1.
Preferred over adding an asymmetric-window API for one caller.

`draw_shape` (`IslamicStars.h:151-228`) needs one change: the per-face hoist loop
(`:194-206`) must **build** a `face_gains` array alongside the palette pointers it
already resolves, and the fragment shader must read it in place of the constant
`1.0f` currently passed to `shade_mesh_topology` (`:213`, `:218`). One more
parallel array on a loop that already exists; no per-fragment cost beyond the
lookup.

### 6.2 Ripple interaction

The ripple burst fires at `fade + still` (`IslamicStars.h:303`) and must move to
`fade + build_span + still`. `RippleTransformer` deforms per frame in
`MeshOps::transform` (`IslamicStars.h:163`) and is orthogonal to the build — but
**do not overlap them**. A rippling mesh whose topology is simultaneously
changing under a gain gate muddles the visual reading and makes the seam metric
in §9.5 unmeasurable — and §9.5 is now the gate that decides whether partition
ops ship at all (§10 Phase 3).
Keep the stages disjoint, as the existing comment at `IslamicStars.h:282-289`
already insists.

### 6.3 Round-robin ordering

Unchanged: `solid_idx = (solid_idx + 1) % solids.size()` (`IslamicStars.h:238`).
See §10 Phase 4 for the seed-affinity ordering that would let consecutive shapes
share a seed and drop the fades entirely.

---

## 7. Timing budget

Leg length must scale with the op's visual weight, not be constant. Proposal,
in frames at the effect's nominal cadence, all divided by `params.trans_speed`
like every other stage (`IslamicStars.h:290-298`):

| Leg kind | Frames | Rationale |
|---|---|---|
| `truncate` / `expand` / `snub` / `chamfer` | 24 | matches ConwayMorph's landed legs (13.8–21.8 ms/frame there) |
| `ambo` (truncate to 0.5) | 24 | same sweep |
| `hankin` | 32 | the star growth is the recipe's signature moment |
| `kis` (gated) | 6 + 1 + 6 = 13 | no renderable parameter (§3.2); gates are the whole leg |
| `dual` (gated) | 6 + 1 + 6 = 13 | no sweep |
| `relax` after a vertex-count-preserving leg | folded into that leg's settle | existing mechanism |
| `relax` after `ambo` | 16, standalone | vertex counts differ, §2.2 |

Worst recipe by leg count after composite expansion:
`dodecahedron_hk35_ambo_hk62_ambo_relax_hk42` → hankin, ambo, hankin, ambo,
relax (standalone, §2.2), hankin = 32+24+32+24+16+32 = **160 frames ≈ 10 s at
16 fps**. `icosahedron_kis_gyro` → kis, snub, dual = 13+24+13 = **50 frames**.

Add that to the existing ~`72 + 16 + 128 + 16 + 72` frame cycle and shapes get
roughly 1.15–1.55× longer. **Accepted** (owner decision, §11): the stretch is fine,
and leg lengths are tuned for legibility rather than to hit a duration target.
The `Trans Speed` slider already divides every stage
(`IslamicStars.h:290-298`), so profiling and impatience both keep an escape
hatch — and leg frames must go through that same divisor, with the same `>= 1`
floor every other stage uses.

One knock-on: the full 24-recipe round-robin now takes roughly twice as long to
cycle. Nothing depends on the period, but it does mean §9.4's soak and any
capture that wants roster coverage need their frame counts doubled.

---

## 8. Memory and performance

### 8.1 Per-frame cost

ConwayMorph measures **173 µs** for op + `compile` on the 18-solid roster
(F ≤ 92). Op cost scales roughly with `I` (Σ face degrees). The Islamic worst
case is F = 1082, so **an order-of-magnitude estimate is ~2 ms/frame** on the
final legs — against a ship peak of ~92 ms/frame that IslamicStars already pays
on its heaviest recipes (rasterizer- and probe-bound, per
`docs/islamicstars_next_levers.md`).

**Early legs are cheap** because the mesh is small; the expensive legs are the
last ones, which are also the ones closest to today's already-measured hold
frame. ConwayMorph found the same shape — "the morph is the effect's lightest
phase".

**Do not predict that the build phase stays under the hold frame** — it cannot,
by construction. A build frame does everything a hold frame does
(`MeshOps::transform` plus `Scan::Mesh::draw` on a mesh that, by the last leg, is
essentially the final solid) *plus* one Conway op, one `MeshOps::compile`, and
the palette bake. The last leg is necessarily slower than the hold.

**The prediction to gate (§9.10): build frame ≤ hold frame + op cost, with the op
cost measured, not extrapolated.** The `173 µs → ~2 ms` figure is
complexity-defensible (`build_half_edge_mesh` is linear-plus-pairing,
`core/mesh/mesh.h:243`) but its constant is untested, and it silently assumes
`MeshOps::compile` also scales linearly. Do not predict 16 fps; the heavy Islamic
recipes run at 8 fps today and this feature does not change that.

### 8.2 Arena

`IslamicStars::init` configures `GLOBAL − 194 KB` persistent / 114 KB scratch_a /
80 KB scratch_b (`IslamicStars.h:36-37`), sized to a measured generation
high-water of 112.3 / 76.5 KB. A morph frame needs, concurrently:

- the leg seed `PolyMesh` — **persistent**, new (previously nothing was held
  mid-shape beyond the two carousel slots),
- `swept` (scratch_a) + `compile` scratch (scratch_b) — as ConwayMorph,
- `MAX_BLEND_PAIRS × 3 KB` blended LUTs (scratch_b) — as ConwayMorph,
- the carousel's two `MeshState` slots — persistent, unchanged.

The seed hold is the new persistent pressure and it peaks on the *last* leg's
seed, which is nearly the final solid. **This is the single largest technical
risk in the feature** on a 298 KB device arena
(`core/engine/memory.h:31`). Mitigations, in order of preference:

1. Drop to a **single** carousel slot during the build — the outgoing shape is
   gone by then (`TerminatorSweep` is sequential: `schedule_sequential` returns
   `duration`, `core/animation/mesh.h:757`; `IslamicStars.h:277-278` already relies on this), so
   slot reuse is legal and buys back a full `MeshState`. **Not free, though:**
   `compact_keep_front` (`core/animation/mesh.h:1164-1170`) is built to keep the *front* slot
   alive across the back slot's regeneration, and `spawn_shape` compacts → 
   generates → flips (`IslamicStars.h:253-274`). Going to one slot restructures
   that ordering, and the front slot is what `Persist` evacuates through
   `scratch_arena_b` — 80 KB against a measured front-slot evacuation of up to
   63.7 KB
   (`IslamicStars.h:33-37`). Evacuating a near-final F = 1082 mesh through that
   path is itself near budget.
2. Restrict the morphable roster (null `recipe` per §4) to recipes whose peak
   fits, measured not guessed.
3. Reconstruct the leg seed from `(seed₀, ops[0..k])` on leg transition only,
   trading a one-off replay for the hold. **Not cheap:** replaying a chain
   containing `.relax(100)` on a 1082-face mesh is a substantial cost —
   `HankinSolids`' own `load_shape` spike (`conway_morph_spec.md:380-383`) is the
   existing evidence for that class — and it would land as a visible hitch at each
   leg boundary unless amortized.

§9.3 makes this a hard CI gate, not a hope. Note also that `Arena` OOM **traps**
(`memory.h`), so a budget miss is a crash on a live art piece, not a glitch.

### 8.3 Roster-scale limits that will break

Current caps, all sized for the 18-solid roster and all exceeded by F = 1082:

- `prev_used[128]` in `build_palette_mapping`, with a hard
  `HS_CHECK(handoff.prev_faces <= 128)` (`core/animation/mesh.h:638-639`) — **traps** at
  F = 1082.
- `MAX_NODE_FACES = 92`, `MAX_HANKIN_FACES = 256` (`HankinSolids.h:115-118`)
- `narrow_face_count` bounds to **UINT8_MAX = 255** (`core/animation/mesh.h:337`) — the tighter
  and more reachable of the two narrowing guards; a `kis` or `needle` on a
  high-valence orbit can hit it. (`narrow_index`, the INT16_MAX one, is at
  `core/animation/mesh.h:323`.)

**The limit that matters most is not a cap — it is a tolerance.**
`PROVENANCE_TOL_SQ = 0.15²` is documented as "under half the face-centroid
spacing of the largest node (~0.37 chord at 92 faces)" (`core/animation/mesh.h:512`). At
F = 1082 the centroid spacing is ~0.108 chord and half-spacing ~0.054, so **the
tolerance is roughly 3× the discriminating distance.** The bijection
`HS_CHECK` at `:653` will keep passing while the nearest-centroid match is simply
wrong. The failure mode is therefore *silent* — misrouted palettes, i.e. exactly
the on-screen "position jump" that `conway_morph_spec.md` overrule 2 exists to
prevent — rather than a trap.

This must be fixed before any high-F recipe morphs: scale the tolerance from the
mesh's actual centroid spacing rather than pinning it to a constant sized for
F ≤ 92, and assert the resulting match is unambiguous (nearest is meaningfully
closer than second-nearest) instead of merely within tolerance.

`build_palette_mapping`'s nearest-centroid search is O(F_prev × F_new). At
1082 × 1082 that is ~1.2 M dot products **once per leg construction** (cold path,
`HS_COLD_MEMBER`), estimated ~5–10 ms. Acceptable as a one-off; if a leg
boundary visibly hitches, bucket the centroids on a coarse spherical grid before
searching. Do not pre-optimize it — measure first.

---

## 9. Test and CI gates

Gates 1–7 and 10 are C++ and wire into `tests/CMakeLists.txt` per project
convention. **Gates 8 and 9 do not:** gate 8 is `scripts/wasm_smoke.mjs` (node,
Holosphere), and gate 9 is `tests/solid_codegen.test.js` — which lives in the
**daydream** repository, as do `tools/solids.html` and `tools/solid_codegen.js`.
Everything §4.3 and §4.4 propose therefore crosses the two-repo split, and Phase 1
is a **two-repo landing**: engine binding on Holosphere master, tool changes on
daydream master. Qualify paths accordingly when implementing.

1. **Recipe ⇔ generator bitwise equivalence.** For every entry with a non-null
   `recipe`: `build_recipe(*e.recipe)` equals `e.generate()` bitwise — vertices,
   `face_counts`, `faces`. This is the anchor gate; without it the recipe table
   is unverified duplication. Extends the existing bitwise-registry discipline
   in `tests/test_solids.h`.

2. **Per-leg topology constancy.** Already asserted at runtime
   (`core/animation/mesh.h:435`); add a host test that steps every leg of every recipe frame by
   frame and checks `compiled.face_counts.size()` never moves within a leg.

3. **Scratch and persistent high-water sweep.** Extend
   `test_edge_morph_frames_fit_scratch_budget` (`tests/test_conway_morph.h:883`) to
   every leg of every morphable recipe, asserting against `IslamicStars`'
   configured split. Include the leg-seed persistent hold from §8.2. **This gate
   decides which recipes get a recipe at all.**

4. **Persistent zero-growth soak.** Extend the `tests/test_conway_soak.h`
   `HankinWalkProbe` pattern to a full round-robin over the morphable roster,
   asserting the persistent arena high-water is flat after N cycles. Arena OOM
   traps, so a leak is a crash.

5. **Crease-pop metric per leg.** Render the frame before and after every bookend
   swap and every gate step via the framebuffer-dump harness
   (`docs/` reference; native render → PNG), and report the fraction of pixels
   changing by more than a threshold. Precedent for the criterion and its teeth:
   the ConwayMorph pyritohedron leg was **dropped** at 91.6 % crease pop.
   Proposed thresholds — **inflate-op bookends ≤ 0.5 %** (they should be
   near-exact, ±2 LSB per the `conway_morph_spec.md` §7.6 bookend pin),
   **gated swaps ≤ 2 %**, **dual swap: measure and decide** (§3.3, calibrated by gate 6).

6. **Gated-swap seam measurement** (not a pass/fail identity test). Render an
   n-gon and its `kis` fan flat-shaded in the same colour and report the pixel
   delta. §3.1 predicts a ~2 px band up to 25 % dark on every interior edge, so
   this is a **calibration** that fixes the threshold for gate 5, not a
   correctness check. If the measured delta is materially worse than predicted,
   §3.3 is dead and §3.4's crossfade is the only path. Run this **before**
   scheduling any partition-op work (§10 Phase 3).

7. **Composite lowering round-trip.** For every recipe:
   `build_recipe(expand_to_primitives(r))` equals `build_recipe(r)` bitwise. This
   is what makes §2's decomposition table an assertion rather than a claim — a
   wrong entry (say `needle = kd` transposed to `dk`) fails here rather than
   producing a subtly wrong solid at runtime.

8. **`getRecipe` reconstruction parity** (`scripts/wasm_smoke.mjs`). For every
   entry with a recipe: applying the reconstructed chain to the seed via the
   `MeshOps` op bindings yields the same V/E/F/I — and the same vertex buffer —
   as `fromSolidName(entry)`. This gates §4.4's whole premise from the JS side,
   independently of §9.1's C++ side.

9. **Codegen round-trip — mesh only, not names** (`tests/solid_codegen.test.js`).
   `generateFuncAndRecipe(item)` takes **one item object**, not `(seed, ops)`
   (`solid_codegen.js:102`). More importantly its suffix rules (`:123-172`)
   **cannot** reproduce most registry names: `_truncate50d` → `_truncate87`,
   `_truncate001` → `_truncate01`, `_relax` → `_relax8`, `_hankin59` → `_hk59`,
   `_snub` → `_snub50_tw00`, `_bevel2` → `_bevel20`. Only a handful
   (`_hk62`, `_ambo`, `_dual`, `_relax100`, `_chamfer63`, `_kis`, `_gyro`) round
   trip. **Gate on the emitted chain building an identical mesh**, and treat name
   reproduction as a separate, narrower assertion over the entries that can
   satisfy it. §4.2's reason for storing authored ops (composites would rename
   `_gyro` → `_snub_dual`) is about not *churning* names, which is a weaker and
   achievable requirement.

10. **Device profile.** Per the standing report cadence: strict per-phase colour,
   🟢 only if every frame of every phase holds its target. Report build-phase
   frames against the shape's own hold frame (§8.1's prediction), not against
   16 fps. Use the `teensy-profile` flow; pin the board with `HS_TEENSY_PORT`.

---

## 10. Phasing

**Phase 0 — the ambo-on-hankin probe.** One host test, before anything else:
step an `ambo` leg on a hankin seed frame by frame and assert constant compiled
face count and a closed manifold throughout (§2.4). Eight of the 13 Phase-2
recipes depend on this holding, including both of Phase 1's. If it fails, `T_EPS`
needs a per-mesh derivation scaled from the smallest face radius rather than a
global constant, and that fix belongs in Phase 1. Cheapest de-risk in the plan.

**Phase 1 — machinery, two recipes.** `Solids::Recipe` + table + equivalence
gate (§9.1). Refactor `ConwayMorph` → `OpLeg`, leaving room for all three kinds
but landing only CONWAY_SWEEP and HANKIN_SWEEP. Populate recipes for
`dodecahedron_hk62_ambo_hk62` and `octahedron_hk17_ambo_hk73` — both are
seed → hankin → ambo → hankin, both are small, and between them they exercise
every Phase-1 leg kind. Everything else keeps a null recipe and today's
behaviour, so the other 22 shapes are untouched.

Note both recipes are ambo-on-hankin, so **Phase 1 carries Phase 0's risk and no
other.** This is also a **two-repo landing** (§9): engine binding on Holosphere
master, tool changes on daydream master.

Land §4.3's `getRecipe` binding and §4.4's tool reconstruction **in this phase,
not later.** It is a few hundred lines against machinery Phase 1 builds anyway,
it makes the recipe table visually verifiable before a single frame of morph
exists, and its recipe-less badge is the coverage dashboard for Phases 2–3. The
tool is also the natural place to author the remaining recipes: load, inspect,
export, paste.

**Phase 2 — the inflate-only roster (13 recipes).** §2.3's pure-inflate set needs
**no partition machinery**, but it is not free: `chamfer` brought under
`clamp_param` (one recipe), the `hankin` leg kind generalized in the scheduler
(all 13), §2.2's standalone relax leg (five recipes, sites at `solids.h:892, 965,
993`), and `expand_to_primitives`' `bevel(0.5) → AMBO` special case (§2.1).
Also fix §8.3's `prev_used[128]` trap and the `PROVENANCE_TOL_SQ` scaling — both
are prerequisites for any high-F recipe. Then pass the §9.3 memory gate, tune leg
frames (§7), device profile and report.

**This is the phase that delivers the feature.** Thirteen of twenty-four recipes,
including the visually richest multi-hankin chains, morph op-by-op with nothing
downstream of §3's rasterizer constraints in the path. Everything after it is
extension.

**Phase 3a — `dual` only, measurement first (6 recipes + 13 Catalans).** Run
§9.6's seam calibration for the `dual` event **before scheduling any
implementation**. It decides between §3.3's gated swap, §3.4's local crossfade,
and not doing it at all. Only then land the `dual` leg. This is the whole value
of the partition work: 19 solids, versus 3 for `kis` (§2.3).

**Phase 3b — `kis` (3 recipes).** Strictly optional and strictly after 3a: all
three `kis` recipes contain `dual`, so `kis` unlocks nothing by itself. Its swap
frame reveals `n·F` new spoke edges at once (§3.3) and §3.1 kills the only
mechanism that could have hidden them. **Judge 3b on 3a's measured result plus
eyeball review, and be willing to drop it** — three recipes keeping today's
TerminatorSweep is a small loss, and the tool (§4.4) makes authoring `kis`-free
alternatives cheap if the roster gap matters.

§2.1's four blocked recipes stay null-recipe unless the `t > 0.5` leg kind is
separately justified.

Do not merge Phases 2 and 3. The whole point of the split is that Phase 3's
premise is unproven and Phase 2's is not.

**Phase 4 — seed affinity (optional, stretch).** The roster has heavy seed reuse:
`dodecahedron` seeds 6 recipes, `truncatedIcosahedron` 7. Order the round-robin
by seed and, between two recipes sharing one, **unbuild** the chain in reverse
(`OpLeg` already supports reverse legs) back to the shared seed, then build the
next — no fade at all. Between different seeds, traverse `ConwayGraph` as
`HankinSolids` does. The endgame is an unbroken continuous morph across the
entire roster with TerminatorSweep only at effect start. Honest cost: unbuilding
roughly doubles on-screen time per shape, so it needs the §7 product call
answered first.

---

## 11. Out of scope, open questions, risks

**Out of scope.**
- Parameterizing `dual` as a continuous sweep — §3.3's gated swap is the answer;
  see §12.
- `propeller` / `whirl` / `loft` — unimplemented (`conway.h:1087-1088` TODO).
- Replacing the generator functions with recipes (§4).
- **Parsing registry names into ops** (§4.1). Lossy and a competing source of
  truth; the recipe table is the answer.
- Persisting the tool's saved-solid collection. It is in-memory today
  (`solids.html:1111-1174`) and stays that way; recipe reconstruction removes the
  main reason to want persistence, since any shipped solid can now be re-loaded
  as an editable chain.
- Hankin-on-Catalan legs, which resonate on kis- and snub-duals (rhombic duals
  are clean). Any Catalan recipe ending in a `hankin` leg needs the per-solid
  check and the far-star guard (`STAR_FAR_BLEND_START_RATIO_SQ = 2.25` through
  `STAR_FAR_RATIO_SQ = 4.0`, gated on `plane_cross_sq`).

**Decisions (owner, 2026-07-20).** All five open questions are answered; the
sections above describe the decided design, not options.

1. **Shape duration — longer is fine.** §7's ~1.15–1.55× per-shape stretch is
   accepted. No compensating cuts to leg counts or the ripple stage. Leg lengths
   are tuned for legibility, not for a duration target.
2. **Palette convergence — spread.** Colour converges across all legs, weighted
   by leg frame count, so colour and geometry arrive together (§5.3).
3. **Seed fade-in — sweep the seed in.** The base solid arrives on
   `TerminatorSweep` as §6.1 proposes. Phase 4's fade-free continuous roster
   stays a stretch goal, not the Phase 3 target.
4. **`HankinSolids` adopting `OpLeg` — not yet.** The §5.1 refactor must leave
   `HankinSolids` working exactly as it does today. Converging the two is a
   separate change, deliberately deferred; do not let it creep into this one.
5. **Tool load semantics — seed + chain.** Reconstruction is the default, not a
   toggle (§4.4). Clearing the ops recovers the opaque view; the reverse is not
   recoverable.

**Risks, ranked.**
1. **The ambo-on-hankin assumption, §2.4.** Eight of the 13 Phase-2 recipes —
   and both Phase-1 recipes — sweep a Conway op on a hankin mesh, while the
   entire `T_EPS` / zero-area-birth analysis is characterized against simple
   solids of ≤ 92 faces. If `compile` drops a degenerate face mid-leg this traps
   at `core/animation/mesh.h:435`. Phase 0 exists to answer it first; it is
   ranked here because nothing else in the plan proceeds if it fails.
2. **Silent palette misrouting at high F, §8.3.** Highest, because it is the only
   one that fails *quietly*. `PROVENANCE_TOL_SQ` is sized for F ≤ 92 and is ~3×
   the discriminating centroid distance at F = 1082; the bijection check keeps
   passing while matches are wrong. Fix before any high-F recipe morphs.
3. **Persistent arena, §8.2.** Holding a near-final `PolyMesh` as the leg seed on
   a 298 KB device arena, where OOM traps. All three mitigations have costs the
   first draft understated. §9.3 forces the answer early.
4. **Partition ops may simply not be doable, §3.** No shading-domain trick makes
   a coplanar partition seamless (§12.1), so what remains is a bounded-seam gate
   whose acceptability is unmeasured, plus a crossfade fallback that contradicts a
   landed aesthetic decision. Phase 3's nine recipes are genuinely at risk —
   *contained* rather than fatal only because §2.3 puts 13 outside it.
5. **`chamfer` is unproven, §2.** It has never been swept and has no `T_EPS`
   characterization, yet one Phase-2 recipe needs it — and that recipe runs it on
   a hankin mesh (§2.4), so the two risks compound.
6. **Roster-scale caps, §8.3.** Mechanical but numerous; these fail as traps, so
   they surface loudly.

**Method rule.** Any claim of the form "this changes no pixels" must be checked
against `Scan::Mesh::draw`'s per-face AA (`scan.h:1238-1241`) and `SDF::Face`'s
gnomonic projection (`sdf.h:2489-2517`) specifically, and against a framebuffer
dump, before it is built on. Both constraints in §3 are invisible from the mesh
and animation layers where operator design naturally happens, and every dead end
in §12 was found there rather than reasoned out in advance.

---

## 12. Measured-dead — do not re-propose

Each of these was designed, then killed by evidence. Evidence retained so the
analysis is not redone.

1. **Flat gate as an exactness mechanism** (collapse shading `gain` to 0 so a
   coplanar partition is pixel-identical to its parent). **DEAD:**
   `Scan::Mesh::draw` antialiases each face against its own boundary and `plot`
   composites *over* (`scan.h:1238-1241`, `filter.h:151`), so every interior edge
   of a partition shows a ~2 px band up to 25 % dark at any gain including 0. The
   gradient and coverage discontinuities are independent; `gain` addresses only
   the first. **Do not re-attempt via any shading-domain trick** — the only exact
   fix is a coverage-accumulating plot path. The gate survives in §3.3 as a
   *seam-reduction* measure, never as exactness.

2. **`kis` apex-radius sweep** (`kis(mesh, target, temp, float h)`, apex swept
   from face plane to sphere). **DEAD:** `SDF::Face` projects gnomonically and
   derives edge normals from cross products (`sdf.h:2489-2517`) — both
   scale-invariant, so the silhouette never moves. Non-uniform radial motion does
   rotate the frame and rescale `size` / `fragment_edge_dist`, so the result is a
   gradient shimmer with no motion. **Do not re-attempt radial "grow it outward"
   animation for `SDF::Face` geometry**; see §3.2 for the precise scope of the
   rule.

3. **Sweeping a `truncate` leg to `t > 0.5`** (the two `truncate50d` recipes).
   **DEAD:** `t = 0.5` short-circuits to `ambo` (`conway.h:567`), halving vertex
   count and every primary degree. `clamp_param` prevents the trap but lands the
   leg on a cuboctahedron-like form that then clean-swaps to the recipe's
   self-intersecting target — a full-screen pop. §2.1.

4. **Deriving recipes by parsing registry names.** **DEAD:** the names are lossy
   — `hk35` is a rounded degree, `truncate50d` is `50·D2R = 0.873` used as a raw
   `t`, `relax` drops its iteration count. A parser would be a weaker second
   source of truth competing with the recipe table. §4.1.

5. **Folding every `relax` into the preceding leg's settle.** **DEAD** for six of
   the nine sites: the settle slerps per-vertex and requires vertex-count
   identity (`core/animation/mesh.h:423`), but an ambo leg sweeps `2E` vertices
   toward a relaxed form with `E`. §2.2.
