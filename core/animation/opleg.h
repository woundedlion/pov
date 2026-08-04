/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file opleg.h
 * @brief Animation fragment: OpLeg Conway-chain morph legs.
 */

#include "mesh/conway.h"
#include "mesh/conway_graph.h"
#include "mesh/hankin.h"

namespace Animation {

/**
 * @brief Animates one operator-sweep leg: a Conway-operator parameter sweep
 * along a graph edge (docs/conway_morph_spec.md, section 4.1) or a recipe
 * step, or a hankin contact-angle sweep on a fixed seed
 * (docs/opchain_morph_spec.md, section 5.1).
 * @details Per frame: produce the swept mesh (the edge's single op at
 * t(frame) settle-slerped toward the relaxed endpoint inside the settle
 * window, or update_hankin at theta(frame)) in scratch, compile, attach the
 * leg's hoisted classification, pre-blend the (from, to) palette ramps at
 * w(frame), and hand the mesh to the draw callback. Exactly one mesh is
 * drawn per frame. Bulk state lives in an arena-allocated Transients that no
 * destructor reclaims; the caller compacts the arena between legs.
 */
class OpLeg : public AnimationBase<OpLeg> {
public:
  static constexpr int PALETTES = BakedPaletteBank::N;
  /** Distinct (from, to) ramp pairs a leg may carry — the full pair space, so
   * no leg can overflow the table. Bounds the per-frame blended-LUT scratch
   * (PAIRS x 2 KB in scratch_arena_b); only the pairs a leg actually uses are
   * allocated, so the ceiling costs nothing until a leg needs it. */
  static constexpr int MAX_BLEND_PAIRS = PALETTES * PALETTES;

  /** Leg kind, dispatched once at construction by the chosen constructor. */
  enum class LegKind : uint8_t {
    CONWAY_SWEEP,
    HANKIN_SWEEP,
    RELAX_SLERP,
    MEDIAL_SLERP,
    GATED_SWAP
  };

  /** Partition operator a GATED_SWAP leg swaps to. */
  enum class SwapOp : uint8_t { KIS, DUAL };

  /**
   * @brief Graph-edge Conway sweep: the edge's single op swept between the
   * edge's endpoint parameters (docs/conway_morph_spec.md, section 4.1).
   */
  struct EdgeSweepSpec {
    const ConwayGraph::EdgeSpec *edge =
        nullptr;           /**< Graph edge being traversed. */
    bool reverse = false;  /**< True when traversing to_node -> from_node. */
    int sweep_frames = 0;  /**< Operator-sweep frames (N). */
    int settle_frames = 0; /**< Relax-slerp frames (S); positive on a settling
                              edge, 0 otherwise. */
  };

  /**
   * @brief Recipe-step Conway sweep: one primitive op swept t_start -> t_end
   * on a fixed seed, no graph edge (docs/opchain_morph_spec.md, section 5.1).
   */
  struct ParamSweepSpec {
    ConwayGraph::MorphOp op;  /**< Swept operator. */
    float t_start = 0.0f;     /**< Sweep parameter at frame 0; clamped to the
                                topology-constant open interval (T_EPS floor;
                                TRUNCATE capped below the ambo point). */
    float t_end = 0.0f;       /**< Arrival parameter; same clamp. */
    float twist_start = 0.0f; /**< Snub twist at frame 0. */
    float twist_end = 0.0f;   /**< Arrival snub twist. */
    int sweep_frames = 0;     /**< Operator-sweep frames (N). */
    bool bridge_provenance =
        false; /**< Dual-bridge leg: geometric provenance takes its start
                  centroids from the leg's own start (closing) or arrival
                  (opening) mesh, since the departed face blocks are
                  transposed against the handoff. */
    bool borrow_seed =
        false; /**< Sweep the caller's live mesh each frame instead of a
                  leg-local clone; legal only where the seed outlives the leg
                  unmoved. */
  };

  /** @brief Hankin contact-angle sweep on a fixed seed. */
  struct HankinSweepSpec {
    float theta_start = 0.0f; /**< Contact angle at frame 0, radians;
                                 negatives clamp to 0, and it enters only as
                                 the opening slerp fraction theta_start /
                                 theta_end, floored at K_EPS. */
    float theta_end = 0.0f;   /**< Arrival contact angle, radians; floored at
                                 THETA_EPS. */
    int sweep_frames = 0;     /**< Angle-sweep frames (N). */
  };

  /**
   * @brief Relax leg: every vertex slerps from its seed position to its
   * relaxed one (docs/opchain_morph_spec.md, section 2.2).
   */
  struct RelaxSpec {
    int iterations = 0; /**< Spring-relaxation passes of the arrival form. */
    const MeshOps::RelaxBake *bake =
        nullptr; /**< Shipped converged relaxation to land on; null instead
                    runs `iterations` live steps, and one of the two is
                    required. */
    int sweep_frames = 0; /**< Slerp frames (N). */
  };

  /** @brief Medial-slerp leg (the Conway-dual bridge's middle leg). */
  struct MedialSpec {
    int sweep_frames = 0; /**< Slerp frames (N). */
  };

  /** @brief Reconcile leg closing a smooth kis/needle path. */
  struct ReconcileSpec {
    const Vector *to_positions =
        nullptr;          /**< Authored vertex positions, one per seed vertex
                    (index-corresponded through the caller's bijection); slerp
                    endpoints. */
    int sweep_frames = 0; /**< Slerp frames (N). */
  };

  /** @brief Gated-swap leg: a partition op with no sweep. */
  struct GatedSwapSpec {
    SwapOp op = SwapOp::KIS; /**< Partition operator applied at the swap
                                frame. */
    int gate_frames = 0;     /**< Frames on each side of the swap (F_gate); the
                            leg runs 2 * gate_frames + 1 frames. */
  };

  /** Floor on a hankin leg's arrival contact angle (the T_EPS analog): the
   * arrival angle divides the opening angle into the leg's opening slerp
   * fraction. */
  static constexpr float THETA_EPS = 0.02f;

  /** Slerp-fraction floor for a hankin leg's opening frame: at fraction 0 every
   * rosette face has exactly zero area, so the floor lifts them past
   * MeshOps::compile's degenerate-face drop (which would move the compiled face
   * count mid-leg). Measured sufficient for every Phase-1 leg. */
  static constexpr float K_EPS = 0.005f;

  /** Default trailing blend window: frames over which a trailing-blend leg's
   * colour diverges from the held source palette to its target
   * (trailing_blend). Colour holds at `from` until this many frames remain,
   * then ramps to `to` at the arrival. Keyed off frames-from-the-end so the
   * window is the same few frames on any leg length: the leg's whole colour
   * movement lands at its very end instead of drifting mid-morph. */
  static constexpr int TRAILING_BLEND_FRAMES = 4;

  /**
   * @brief Per-frame shading handed to the draw callback.
   * @details The fragment path stays a single BakedPalette::get(t):
   * ramps[face_ramp[face]] is the face's pre-blended LUT. Scratch-backed,
   * valid for the current frame only.
   */
  struct Shading {
    const BakedPalette *ramps; /**< One blended LUT per (from, to) pair. */
    const uint8_t *face_ramp;  /**< Face index -> ramp index. */
    size_t faces;              /**< Face count (bounds face_ramp). */
    float gain = 1.0f; /**< Multiplier the shader applies to the edge-distance
                          gradient; OpLeg holds it at 1 on every kind. */

    /**
     * @brief The face's pre-blended ramp.
     * @param face Face index.
     * @return The blended LUT to sample; an out-of-range face falls back to
     * ramp 0 rather than reading past the table.
     */
    const BakedPalette &ramp_for(size_t face) const {
      return ramps[face < faces ? face_ramp[face] : 0];
    }
  };

  /**
   * @brief Non-owning draw callback: `void(Canvas&, MeshState&,
   * const Shading&)`. StoredFunctionRef rejects rvalue temporaries.
   */
  using MorphDrawFn =
      StoredFunctionRef<void(Canvas &, MeshState &, const Shading &)>;

  /**
   * @brief Crossfade-weight curve of a swept leg: resolved weight in [0, 1]
   * for a 1-based leg frame over the whole leg duration (sweep plus settle).
   * Supplied at construction (classic_blend default, trailing_blend for the
   * build legs), mirroring easing_fn.
   */
  using BlendWeightFn = float (*)(int frame, int duration);

  /** @brief Face-order relationship between departed and swept meshes. */
  enum class FaceCorrespondence : uint8_t {
    GEOMETRIC,
    IDENTITY,
    DUAL_CLOSING,
  };

  /**
   * @brief Palette provenance of the departed node (spec sections 2.5/2.6).
   * @details prev_face_palette describes the node base mesh the leg departs
   * from, in emission order; consumed by the constructor only. When
   * prev_face_centroid is supplied, the mapping is geometric (see
   * build_palette_mapping).
   */
  struct PaletteHandoff {
    const BakedPaletteBank *bank =
        nullptr; /**< The effect's baked source LUTs. */
    const uint8_t *prev_face_palette =
        nullptr;           /**< Per-face palette of the departed base mesh. */
    size_t prev_faces = 0; /**< Face count of the departed base mesh. */
    const Vector *prev_face_centroid =
        nullptr; /**< Unit centroid per departed base face; enables the
                    geometric provenance mapping. */
    FaceCorrespondence correspondence =
        FaceCorrespondence::GEOMETRIC; /**< Palette provenance mapping. */
  };

  /**
   * @brief Bookend grouping of the arrival node (spec sections 2.5/2.6).
   * @details topology[f] is the class the effect displays arrival face f
   * with at the closing bookend (the hankin star-face classification of the
   * arrival base mesh, which can be coarser than both the swept and the
   * clean-endpoint classifications). Color targets key on it, so every face
   * that grouping merges converges to one color by w = 1 and the swap
   * changes nothing. A null topology keys targets on the swept
   * classification instead.
   */
  struct BookendClasses {
    const uint16_t *topology = nullptr; /**< Bookend class per arrival face. */
    size_t faces = 0;                   /**< Arrival base-mesh face count. */
  };

  /**
   * @brief Leg-static arrival data the effect's completion consumes.
   * @details Arena-backed; valid until the leg arena is compacted.
   */
  struct Landing {
    const uint16_t *topology =
        nullptr; /**< Target (bookend-keyed) classification, one id per swept
                    face. */
    size_t faces = 0;         /**< Swept face count. */
    size_t primary_faces = 0; /**< Seed face count (emission-order prefix). */
    std::array<uint8_t, PALETTES>
        to_palette{}; /**< Slot -> landed palette index. */
    const uint8_t *from_palette =
        nullptr;         /**< Per-swept-face FROM palette id (arena-backed). */
    int blend_pairs = 0; /**< Distinct (from, to) ramp pairs the leg carries. */

    /**
     * @brief Palette id face f displays once the leg arrives (w = 1).
     * @param f Swept face index.
     * @details Every face — newborn or carried — lands on its target-class
     * palette.
     */
    uint8_t landed_palette(size_t f) const {
      return to_palette[wrap(static_cast<int>(topology[f]), PALETTES)];
    }
    const PolyMesh *arrival_topology =
        nullptr; /**< Fixed connectivity of a packed arrival endpoint. */
    const Snorm3 *arrival_point =
        nullptr; /**< Packed arrival vertices for a later bridge handoff. */
    size_t arrival_points = 0; /**< Length of arrival_point. */
    const CompiledHankin *hankin =
        nullptr; /**< Baked arrival topology (HANKIN_SWEEP legs only, null
                    otherwise); with star_point it rebuilds the arrival mesh
                    through arrival_mesh(). */
    const Snorm3 *star_point =
        nullptr;            /**< Arrival star points (snorm16-packed). */
    size_t star_points = 0; /**< Star-point count. */
  };

  /**
   * @brief Constructs a graph-edge sweep leg: clones the seed, computes the
   * arrival classification (relaxed form when settling), and builds the
   * palette mappings.
   * @param seed Seed mesh the op sweeps on (cloned, not borrowed).
   * @param spec Traversed edge and frame counts.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed node.
   * @param bookend Bookend grouping of the arrival node (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the sweep parameter.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const EdgeSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames + spec.settle_frames, false),
        easing_fn(easing_fn), draw_fn(draw) {
    HS_CHECK(spec.edge, "OpLeg: edge sweep carries no graph edge");
    const ConwayGraph::EdgeSpec &edge = *spec.edge;
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(spec.settle_frames >= 0 && edge.settle == (spec.settle_frames > 0),
             "OpLeg: settle frames disagree with the edge");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: edge sweep leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    MeshOps::clone(seed, tr.seed, arena);
    tr.seed_ref = &tr.seed;
    tr.seed_faces = seed.face_counts.size();
    tr.op = edge.op;
    tr.reverse = spec.reverse;
    tr.sweep_frames = spec.sweep_frames;
    tr.settle_frames = spec.settle_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;

    // Clamp both endpoints inside the topology-constant open interval; the
    // truncate upper clamp dodges the ambo short-circuit at exactly 0.5, and
    // the jitterbug clamp stops the snub bridge while its collapsing edge is
    // still a positive 0.02-chord sliver.
    auto clamp_param = [&](float t) {
      t = std::max(t, ConwayGraph::T_EPS);
      if (edge.op == ConwayGraph::MorphOp::TRUNCATE)
        t = std::min(t, 0.5f - ConwayGraph::T_EPS_AMBO);
      if (ConwayGraph::is_jitterbug_edge(edge))
        t = std::max(t, ConwayGraph::T_EPS_JITTERBUG);
      return t;
    };
    tr.t_start = clamp_param(spec.reverse ? edge.t_to : edge.t_from);
    tr.t_end = clamp_param(spec.reverse ? edge.t_from : edge.t_to);
    tr.twist_start = spec.reverse ? edge.twist_to : edge.twist_from;
    tr.twist_end = spec.reverse ? edge.twist_from : edge.twist_to;

    init_conway(handoff, bookend, arena, edge.settle,
                ConwayGraph::is_jitterbug_edge(edge),
                /*bridge_provenance=*/false);
  }

  /**
   * @brief Constructs a recipe-step Conway sweep leg: one primitive op swept
   * t_start -> t_end on a fixed seed, no graph edge
   * (docs/opchain_morph_spec.md, section 5.1).
   * @param seed Seed mesh the op sweeps on (cloned unless spec.borrow_seed).
   * @param spec Swept operator, parameter endpoints and frame count.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the sweep parameter.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const ParamSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: param sweep leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    // Borrowed seed: the swept op reads the caller's live mesh each frame
    // instead of a leg-local clone, dropping one full copy of the (tripled)
    // dual-bridge seed from the persistent arena. Only legal where the source
    // outlives the leg unmoved (the dt-bridge seed is persistent-resident and
    // untouched until the leg completes); every other caller clones.
    if (spec.borrow_seed) {
      tr.seed_ref = &seed;
    } else {
      MeshOps::clone(seed, tr.seed, arena);
      tr.seed_ref = &tr.seed;
    }
    tr.seed_faces = seed.face_counts.size();
    tr.op = spec.op;
    tr.sweep_frames = spec.sweep_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;

    // Truncate births below T_EPS when the arrival is itself below T_EPS, so a
    // 0.01 target sweeps from a smaller positive birth instead of clamping both
    // endpoints to T_EPS (a still image). Every arrival >= 0.1 keeps the T_EPS
    // birth unchanged.
    const bool truncate = spec.op == ConwayGraph::MorphOp::TRUNCATE;
    // A far-side truncate leg reaches past the ambo pinch at either endpoint
    // and sweeps through it on the constant-topology truncate branch; the
    // ambo-equivalent leg (both endpoints <= 0.5, exactly 0.5 included) keeps
    // its 0.495 cap and clean-swaps to ambo.
    const bool far_side = truncate && std::max(spec.t_start, spec.t_end) > 0.5f;
    const float trunc_floor =
        std::min(ConwayGraph::T_EPS, std::max(spec.t_start, spec.t_end) *
                                         ConwayGraph::T_EPS_TRUNCATE_FRAC);
    auto clamp_param = [&](float t) {
      t = std::max(t, truncate ? trunc_floor : ConwayGraph::T_EPS);
      if (truncate)
        t = std::min(t, far_side ? ConwayGraph::T_EPS_TRUNCATE_FAR_MAX
                                 : 0.5f - ConwayGraph::T_EPS_AMBO);
      return t;
    };
    tr.t_start = clamp_param(spec.t_start);
    tr.t_end = clamp_param(spec.t_end);
    tr.twist_start = spec.twist_start;
    tr.twist_end = spec.twist_end;

    init_conway(handoff, bookend, arena, /*settle=*/false,
                /*jitterbug=*/false, spec.bridge_provenance);
  }

  /**
   * @brief Constructs a hankin-sweep leg: bakes the angle-independent hankin
   * topology once, computes the arrival classification, and builds the palette
   * mappings.
   * @param seed Base mesh the hankin pattern sweeps on; read only here, never
   * cloned into the leg arena.
   * @param spec Contact-angle endpoints and frame count.
   * @param arena Leg arena backing the compiled hankin topology and hoisted
   * state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed node.
   * @param bookend Bookend grouping of the arrival node (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the sweep angle.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const HankinSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: hankin sweep leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    // No seed clone: the compiled hankin topology carries the base vertices
    // and every per-frame read goes through it, so the seed is needed only
    // inside this constructor.
    tr.seed_faces = seed.face_counts.size();
    tr.kind = LegKind::HANKIN_SWEEP;
    tr.sweep_frames = spec.sweep_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;

    // Hankin legs sweep the slerp fraction, not the contact angle: re-solving
    // the contact-plane intersection per frame sends star points on geodesic
    // excursions far outside their own rosette mid-sweep on ambo-of-hankin
    // seeds, tripping the far-star fallback and flipping branches, which draws
    // as lines crossing the pattern. Growing each star point out from its
    // collapsed corner is monotone and lands on the same arrival geometry.
    HS_CHECK(spec.theta_start <= spec.theta_end,
             "OpLeg: hankin leg sweeps back to a smaller contact angle");
    const float theta_hi = std::max(spec.theta_end, THETA_EPS);
    tr.t_start = std::max(std::max(spec.theta_start, 0.0f) / theta_hi, K_EPS);
    tr.t_end = 1.0f;

    {
      ScratchScope sa(scratch_arena_a);
      ScratchScope sb(scratch_arena_b);

      MeshOps::compile_hankin(seed, tr.hankin, arena, scratch_arena_a,
                              /*borrow_base_vertices=*/true);

      PolyMesh arrival;
      MeshOps::update_hankin(tr.hankin, arrival, scratch_arena_a, theta_hi);

      // Hoist classes before the snorm16 pack unless the exact bookend already
      // supplies them; quantization can split classifier angle buckets.
      hoist_arrival_topology(arrival, bookend, arena, tr.topo);

      // Only the star points are stored: the midpoint prefix is already in the
      // compiled topology, and each star point's collapsed position is its own
      // corner, reachable through the same instruction hankin_at walks.
      const size_t statics = tr.hankin.static_vertices.size();
      HS_CHECK(arrival.vertices.size() >= statics);
      const size_t dyn = arrival.vertices.size() - statics;
      tr.hk_final.bind(arena, dyn);
      for (size_t i = 0; i < dyn; ++i)
        tr.hk_final.push_back(Snorm3::encode(arrival.vertices[statics + i]));
      tr.landing.hankin = &tr.hankin;
      tr.landing.star_point = tr.hk_final.data();
      tr.landing.star_points = tr.hk_final.size();

      const Vector *start_centroid = nullptr;
      PolyMesh start_mesh;
      if (handoff.prev_face_centroid) {
        hankin_at(tr, start_mesh, scratch_arena_b, tr.t_start);
        HS_CHECK(start_mesh.face_counts.size() == tr.topo.size(),
                 "OpLeg: start face count differs from arrival");
        start_centroid = face_centroids(start_mesh, scratch_arena_b);
      }

      // Star faces are emitted first, in base-face order, so the seed's face
      // count is the emission-order prefix corresponding 1:1 to it.
      const size_t survivors = tr.seed_faces;
      build_palette_mapping(tr, arrival, handoff, bookend, arena,
                            start_centroid, survivors);
    }
  }

  /**
   * @brief Constructs a relax leg: clones the seed, relaxes it once, and
   * slerps every vertex from its seed position to its relaxed one
   * (docs/opchain_morph_spec.md, section 2.2).
   * @param seed Mesh being relaxed; its geometry is cloned, its class ids are
   * not.
   * @param spec Relaxation source (live iterations or bake) and frame count.
   * @param arena Leg arena backing the cloned seed geometry and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the slerp fraction.
   * @note relax preserves vertex count and vertex order, so the leg needs no
   * correspondence pass and its topology is constant by construction.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const RelaxSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(spec.bake || spec.iterations >= 1,
             "OpLeg: relax leg needs a positive iteration count");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: relax leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    // relax_at slerps out of the seed vertices every frame, so the seed stays
    // in the leg arena; only its per-face class ids are dead here.
    clone_geometry(seed, tr.seed, arena);
    tr.seed_faces = seed.face_counts.size();
    tr.kind = LegKind::RELAX_SLERP;
    tr.sweep_frames = spec.sweep_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;
    tr.t_start = 0.0f;
    tr.t_end = 1.0f;

    {
      ScratchScope sa(scratch_arena_a);
      ScratchScope sb(scratch_arena_b);

      // With a bake the leg lands on the shipped converged mesh (the generator
      // used relax_baked too); otherwise it runs `iterations` live steps.
      PolyMesh arrival =
          spec.bake ? MeshOps::relax_baked(tr.seed, scratch_arena_a, *spec.bake)
                    : MeshOps::relax(tr.seed, scratch_arena_a, scratch_arena_b,
                                     spec.iterations);
      HS_CHECK(arrival.vertices.size() == tr.seed.vertices.size(),
               "OpLeg: relax changed the vertex count");
      tr.relaxed.bind(arena, arrival.vertices.size());
      tr.relaxed.append_bulk(arrival.vertices.data(), arrival.vertices.size());

      hoist_arrival_topology(arrival, bookend, arena, tr.topo);

      // The opening frame is the seed verbatim, so its face centroids are the
      // provenance source.
      const Vector *start_centroid = nullptr;
      if (handoff.prev_face_centroid)
        start_centroid = face_centroids(tr.seed, scratch_arena_a);

      build_palette_mapping(tr, arrival, handoff, bookend, arena,
                            start_centroid, tr.seed_faces);
    }
  }

  /**
   * @brief Constructs the Conway-dual bridge's medial leg: builds the shared
   * rectified connectivity of P and slerps every vertex from ambo(P) to
   * ambo(dual(P)) at a fixed emission order (docs conway dual morph, leg 2).
   * @param seed Mesh whose dual bridge this leg spans; its medial is built here.
   * @param spec Slerp frame count.
   * @param arena Leg arena backing the medial connectivity and both
   * snorm16-packed vertex sets.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh (ambo(P), one face
   * per primal face + one per primal vertex).
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the slerp fraction.
   * @note The connectivity is ambo(P) exactly, so the leg holds one fixed
   * emission order and every medial face keeps its identity across the slerp;
   * only the vertex positions move. A degree-2 seed vertex makes the s=1
   * positions many-to-one (a lossy dual), which leaves the faces well-formed.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const MedialSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: medial leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    tr.kind = LegKind::MEDIAL_SLERP;
    tr.sweep_frames = spec.sweep_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;
    tr.t_start = 0.0f;
    tr.t_end = 1.0f;

    {
      ScratchScope sa(scratch_arena_a);
      ScratchScope sb(scratch_arena_b);

      // The medial is built in scratch; only its connectivity (tr.seed) and the
      // two snorm16-packed endpoint sets cross into persistent. a_e (ambo(P))
      // and b_e (ambo(dual(P))) are unit-sphere directions, so the 6-byte pack
      // halves their 12-byte resident cost (mirrors hk_final). Classification
      // and centroids run on the full-precision scratch copies: unlike hk_final
      // no arrival rebuild reads the packed points back, and the topology
      // classifier splits symmetry orbits under the quantum, so perturbing its
      // input would inflate the leg's distinct palette-pair count.
      PolyMesh med;
      ArenaVector<Vector> med_b;
      MeshOps::medial(seed, med, med_b, scratch_arena_a, scratch_arena_b);
      tr.seed_faces = med.face_counts.size();
      HS_CHECK(med_b.size() == med.vertices.size(),
               "OpLeg: medial vertex sets differ in size");

      const size_t medial_verts = med.vertices.size();
      copy_topology(tr.seed, arena, med.face_counts, med.faces);
      tr.medial_a.bind(arena, medial_verts);
      tr.medial_b.bind(arena, medial_verts);
      for (size_t i = 0; i < medial_verts; ++i) {
        tr.medial_a.push_back(Snorm3::encode(med.vertices[i]));
        tr.medial_b.push_back(Snorm3::encode(med_b[i]));
      }

      // The opening frame is ambo(P) verbatim (the packed a_e), so med's s=0
      // vertices are the geometric provenance source; snapshot their centroids
      // before the in-place overwrite below.
      const Vector *start_centroid = nullptr;
      if (handoff.prev_face_centroid)
        start_centroid = face_centroids(med, scratch_arena_a);

      // Prepare the s=1 arrival (ambo(dual(P))) in place. The exact bookend can
      // supply its classes; otherwise the full-precision mesh is classified.
      for (size_t i = 0; i < medial_verts; ++i)
        med.vertices[i] = med_b[i];

      hoist_arrival_topology(med, bookend, arena, tr.topo);
      tr.landing.arrival_topology = &tr.seed;
      tr.landing.arrival_point = tr.medial_b.data();
      tr.landing.arrival_points = tr.medial_b.size();

      // Full correspondence: every medial face lives the whole bridge, so the
      // survivor prefix is the whole face list.
      build_palette_mapping(tr, med, handoff, bookend, arena, start_centroid,
                            tr.seed_faces);
    }
  }

  /**
   * @brief Constructs the reconcile leg closing a smooth kis/needle path: slerps
   * every vertex of the topology-exact identity mesh (dt/dtd) onto its
   * counterpart in the authored kis/needle mesh, along the caller's
   * nearest-vertex bijection (docs/opchain_morph_spec.md, smooth kis/needle).
   * @param from_mesh Identity mesh (dual(truncate(...))): its connectivity and
   * vertex order are the leg's fixed emission order. Cloned, not borrowed.
   * @param spec Authored slerp endpoints and frame count.
   * @param arena Leg arena backing the connectivity and both snorm16-packed
   * vertex sets.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh (the identity mesh;
   * one face per identity face).
   * @param bookend Bookend grouping of the arrival mesh (target keying).
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the slerp fraction.
   * @note Identical connectivity is held fixed and only the positions move, so
   * the leg reuses the MEDIAL_SLERP step path verbatim; the caller owns the
   * nearest-vertex correspondence (a checked bijection between two near-identical
   * meshes) so the leg needs no per-frame matching.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &from_mesh, const ReconcileSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(spec.sweep_frames, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
    HS_CHECK(spec.to_positions, "OpLeg: reconcile leg carries no endpoints");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: reconcile leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    tr.kind = LegKind::MEDIAL_SLERP;
    tr.sweep_frames = spec.sweep_frames;
    tr.bank = handoff.bank;
    tr.blend_fn = blend_fn;
    tr.t_start = 0.0f;
    tr.t_end = 1.0f;

    {
      ScratchScope sa(scratch_arena_a);
      ScratchScope sb(scratch_arena_b);

      const size_t n = from_mesh.vertices.size();
      tr.seed_faces = from_mesh.face_counts.size();
      copy_topology(tr.seed, arena, from_mesh.face_counts, from_mesh.faces);
      // Both endpoint sets are unit-sphere directions, snorm16-packed to halve
      // their resident cost (mirrors the medial leg).
      tr.medial_a.bind(arena, n);
      tr.medial_b.bind(arena, n);
      for (size_t i = 0; i < n; ++i) {
        tr.medial_a.push_back(Snorm3::encode(from_mesh.vertices[i]));
        tr.medial_b.push_back(Snorm3::encode(spec.to_positions[i]));
      }

      // Arrival: the identity connectivity carrying the authored positions.
      PolyMesh arrival;
      arrival.vertices.bind(scratch_arena_a, n);
      arrival.vertices.append_bulk(spec.to_positions, n);
      copy_topology(arrival, scratch_arena_a, from_mesh.face_counts,
                    from_mesh.faces);

      hoist_arrival_topology(arrival, bookend, arena, tr.topo);

      // The opening frame is the identity mesh verbatim (the packed a_e), so its
      // face centroids are the geometric provenance source.
      const Vector *start_centroid = nullptr;
      if (handoff.prev_face_centroid)
        start_centroid = face_centroids(from_mesh, scratch_arena_a);

      build_palette_mapping(tr, arrival, handoff, bookend, arena,
                            start_centroid, tr.seed_faces);
    }
  }

  /**
   * @brief Constructs a gated-swap leg: a partition op with no sweep, drawn as
   * the seed then op(seed) at constant gain, the swap masked by holding the
   * inherited source colour until the late fade (docs/opchain_morph_spec.md,
   * section 3.3).
   * @param seed Mesh the partition op runs on (cloned, not borrowed).
   * @param spec Partition operator and gate length.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @param easing_fn Unused by the gate (no sweep); kept for the shared
   * constructor signature.
   * @note Radial and apex motion are invisible to SDF::Face (spec 3.2), so
   * there is no sweep segment; the leg's compiled face count is constant on
   * each side of the swap and changes exactly once, at it.
   */
  HS_COLD_MEMBER
  OpLeg(const PolyMesh &seed, const GatedSwapSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(2 * spec.gate_frames + 1, false), easing_fn(easing_fn),
        draw_fn(draw) {
    HS_CHECK(spec.gate_frames >= 1, "OpLeg needs a positive gate length");
    HS_CHECK(handoff.bank && handoff.prev_face_palette &&
                 handoff.prev_faces > 0,
             "OpLeg: gated swap leg has an incomplete palette handoff");
    Transients &tr = bind_transients(arena);

    MeshOps::clone(seed, tr.seed, arena);
    tr.seed_faces = seed.face_counts.size();
    tr.kind = LegKind::GATED_SWAP;
    tr.swap_op = spec.op;
    tr.sweep_frames = spec.gate_frames;
    tr.bank = handoff.bank;

    init_gated(handoff, bookend, arena);
  }

  /**
   * @brief Steps the sweep: the kind's swept mesh (op at t(frame) plus settle
   * slerp, update_hankin at theta(frame), or the relax slerp), then compile,
   * palette pre-blend, draw.
   * @param canvas The canvas passed through to the draw callback.
   * @details HS_COLD: once-per-frame orchestration; the hot loops live in the
   * (already cold) Conway ops and the mesh scan.
   */
  HS_COLD_MEMBER void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    draw_frame(canvas);
  }

  /** @brief Redraws the current morph frame without advancing the leg. */
  HS_COLD_MEMBER void step_paused(Canvas &canvas) override {
    draw_frame(canvas);
  }

private:
  HS_COLD_MEMBER void draw_frame(Canvas &canvas) {
    check_alive();
    Transients &tr = *buf;
    const int frame = static_cast<int>(
        std::min<uint32_t>(t, static_cast<uint32_t>(duration)));

    if (tr.kind == LegKind::GATED_SWAP) {
      step_gated(canvas, frame);
      return;
    }

    // Frame -> sweep position + settle blend. Settling legs run one eased
    // clock over the whole leg, split proportionally by frame counts: the
    // sweep/settle seam lands mid-easing instead of at the sweep's
    // zero-slope tail, so per-frame motion crosses it without a hitch while
    // the leg still starts and ends at zero slope against the static
    // bookends. Forward legs settle at the end; reverse legs un-settle over
    // the opening window, symmetrically.
    const float progress =
        easing_fn(static_cast<float>(frame) / static_cast<float>(duration));
    float k = progress;
    float settle_alpha = 0.0f;
    if (tr.settle_frames > 0) {
      float split = (tr.reverse ? tr.settle_frames : tr.sweep_frames) /
                    static_cast<float>(duration);
      if (!tr.reverse) {
        k = std::min(progress / split, 1.0f);
        settle_alpha = std::max(0.0f, (progress - split) / (1.0f - split));
      } else {
        settle_alpha = 1.0f - std::min(progress / split, 1.0f);
        k = std::max(0.0f, (progress - split) / (1.0f - split));
      }
    }
    float tp = tr.t_start + (tr.t_end - tr.t_start) * k;
    float tw = tr.twist_start + (tr.twist_end - tr.twist_start) * k;

    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    PolyMesh swept;
    {
      HS_PROFILE(hk_conway_sweep);
      if (tr.kind == LegKind::HANKIN_SWEEP) {
        hankin_at(tr, swept, scratch_arena_a, tp);
      } else if (tr.kind == LegKind::RELAX_SLERP) {
        relax_at(tr, swept, scratch_arena_a, tp);
      } else if (tr.kind == LegKind::MEDIAL_SLERP) {
        medial_at(tr, swept, scratch_arena_a, tp);
      } else {
        HS_CHECK(tr.seed_ref, "OpLeg: leg carries no seed mesh");
        swept = run_op(tr.op, *tr.seed_ref, scratch_arena_a, scratch_arena_b,
                       tp, tw);
        if (settle_alpha > 0.0f) {
          HS_CHECK(swept.vertices.size() == tr.relaxed.size());
          // Alpha 1 copies the relaxed endpoint verbatim, so the settled
          // bookend is bitwise it (mirrors slerp_vertices' k >= 1 shortcut).
          if (settle_alpha >= 1.0f)
            for (size_t i = 0; i < swept.vertices.size(); ++i)
              swept.vertices[i] = tr.relaxed[i];
          else
            for (size_t i = 0; i < swept.vertices.size(); ++i)
              swept.vertices[i] =
                  slerp(swept.vertices[i], tr.relaxed[i], settle_alpha);
        }
      }
    }

    finish_frame(canvas, swept, tr.blend_fn(frame, duration));
  }

public:
  /**
   * @brief Arrival data for the effect's completion handler.
   * @return Arena-backed Landing; stable until the leg arena is compacted.
   */
  const Landing &landing() const {
    check_alive();
    return buf->landing;
  }

  /**
   * @brief Rebuilds a hankin leg's arrival mesh from its baked topology.
   * @param landing Landing of a HANKIN_SWEEP leg.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @details The static vertices and the topology are MeshOps::hankin's, but
   * the star points are snorm16-decoded and renormalized, so the result is
   * bitwise the leg's own final swept frame rather than MeshOps::hankin
   * evaluated at the arrival angle. A chain reads its clean endpoint from the
   * finished leg instead of carrying a second copy of it through the sweep.
   */
  HS_COLD_MEMBER static void arrival_mesh(const Landing &landing, PolyMesh &out,
                                          Arena &arena) {
    HS_CHECK(landing.hankin, "OpLeg: leg carries no baked arrival");
    const CompiledHankin &hk = *landing.hankin;
    const size_t statics = hk.static_vertices.size();
    out.vertices.bind(arena, statics + landing.star_points);
    out.vertices.append_bulk(hk.static_vertices.data(), statics);
    for (size_t i = 0; i < landing.star_points; ++i)
      out.vertices.push_back(landing.star_point[i].decode().normalized());
    copy_topology(out, arena, hk.face_counts, hk.faces);
  }

  /**
   * @brief Mid-leg crossfade weight (spec 2.6), the swept ctors' default:
   * exactly 0 through the first 20% of the leg, smoothstep to exactly 1 by
   * 80%.
   * @param frame 1-based leg frame (1..duration).
   * @param duration Whole leg length in frames (sweep plus settle).
   */
  static float classic_blend(int frame, int duration) {
    return blend_weight(static_cast<float>(frame) /
                        static_cast<float>(duration));
  }

  /**
   * @brief Trailing-window crossfade weight: exactly 0 until the final @p
   * window frames of the leg, then smoothstep to exactly 1 at the arrival
   * frame.
   * @param frame 1-based leg frame (1..duration).
   * @param duration Whole leg length in frames (sweep plus settle).
   * @param window Trailing window in frames, clamped to duration - 1 so a leg
   * shorter than the window still opens on the zero plateau.
   */
  static float trailing_blend(int frame, int duration, int window) {
    const int win = std::max(1, std::min(window, duration - 1));
    const int from_end = duration - frame;
    if (from_end >= win)
      return 0.0f;
    const float u =
        static_cast<float>(win - from_end) / static_cast<float>(win);
    return u * u * (3.0f - 2.0f * u);
  }

  /**
   * @brief Trailing blend over the default TRAILING_BLEND_FRAMES window (the
   * build legs and the gate); the BlendWeightFn form.
   * @param frame 1-based leg frame (1..duration).
   * @param duration Whole leg length in frames (sweep plus settle).
   * @details Holds the inherited source palette through most of the leg, then
   * diverges to the target only over the last few frames. Reaching exactly 1 on
   * the final frame lands the leg on its target colour before the next leg
   * departs.
   */
  static float trailing_blend(int frame, int duration) {
    return trailing_blend(frame, duration, TRAILING_BLEND_FRAMES);
  }

  /**
   * @brief Sums the positions of one face's vertices — the unnormalized
   * centroid, leaving the degenerate-face policy to the caller.
   * @param vertices Vertex positions the face indices address.
   * @param faces Flat per-face vertex index list.
   * @param off Start of this face's span in `faces`.
   * @param n Side count of this face.
   * @return Sum of the face's vertex positions.
   */
  __attribute__((always_inline)) static Vector
  face_vertex_sum(const Vector *vertices, const uint16_t *faces, size_t off,
                  int n) {
    Vector c(0.0f, 0.0f, 0.0f);
    for (int k = 0; k < n; ++k)
      c = c + vertices[faces[off + k]];
    return c;
  }

  /**
   * @brief Fills a caller-allocated array with every face's unit centroid.
   * @param m Mesh whose faces are reduced.
   * @param out Receives one unit vector per face, in emission order.
   */
  HS_COLD_MEMBER static void face_centroids_into(const PolyMesh &m,
                                                 Vector *out) {
    size_t off = 0;
    for (size_t f = 0; f < m.face_counts.size(); ++f) {
      const int n = m.face_counts[f];
      out[f] = face_vertex_sum(m.vertices.data(), m.faces.data(), off, n)
                   .normalized();
      off += n;
    }
  }

private:
  /**
   * @brief Arena-allocated leg state — keeps OpLeg inline size small.
   */
  struct Transients {
    PolyMesh seed; /**< Cloned leg seed; empty on HANKIN_SWEEP legs and on a
                      borrowed CONWAY_SWEEP seed. */
    const PolyMesh *seed_ref =
        nullptr; /**< CONWAY_SWEEP swept-op source: &seed for a cloned seed, the
                    caller's live mesh for a borrowed one (dual-bridge legs whose
                    seed is already persistent-resident for the leg's whole
                    life). */
    size_t seed_faces = 0; /**< Seed face count; set by every kind, including
                              the ones that keep no seed. */
    LegKind kind = LegKind::CONWAY_SWEEP; /**< Swept-mesh production path. */
    CompiledHankin hankin; /**< Baked topology (HANKIN_SWEEP legs). */
    ArenaVector<Snorm3>
        hk_final; /**< Arrival star points, snorm16-packed (HANKIN_SWEEP). */
    ConwayGraph::MorphOp op =
        ConwayGraph::MorphOp::TRUNCATE; /**< Swept operator (CONWAY_SWEEP). */
    SwapOp swap_op = SwapOp::KIS;       /**< Partition op (GATED_SWAP). */
    bool reverse = false;               /**< Traversing to_node -> from_node. */
    int sweep_frames = 1;               /**< Operator-sweep frames. */
    int settle_frames = 0;              /**< Relax-slerp frames. */
    float t_start = 0, t_end = 0;       /**< Clamped sweep endpoints. */
    float twist_start = 0, twist_end = 0; /**< Snub twist endpoints. */
    ArenaVector<Vector>
        relaxed; /**< Relaxed endpoint vertices (settling and relax legs). */
    ArenaVector<Snorm3>
        medial_a; /**< Medial a_e (ambo(P)) endpoint, snorm16-packed
                     (MEDIAL_SLERP). */
    ArenaVector<Snorm3>
        medial_b; /**< Medial b_e (ambo(dual(P))) endpoint, snorm16-packed
                     (MEDIAL_SLERP). */
    ArenaVector<uint16_t>
        topo; /**< Hoisted arrival classification (drawn grouping). */
    ArenaVector<uint16_t>
        target_topo; /**< Per-swept-face target class (bookend grouping). */
    ArenaVector<uint8_t> face_ramp; /**< Face -> (from, to) ramp pair. */
    ArenaVector<uint16_t>
        seed_topo; /**< Seed classification (GATED_SWAP closing half). */
    ArenaVector<uint8_t>
        seed_face_ramp; /**< Seed face -> identity ramp (GATED_SWAP). */
    uint8_t seed_ramp_pal[PALETTES] = {}; /**< Departed palette per identity
                                             ramp (GATED_SWAP). */
    int seed_num_ramps = 0;               /**< Distinct departed palettes. */
    const BakedPaletteBank *bank =
        nullptr; /**< Source LUTs for the pre-blend. */
    uint8_t ramp_from[MAX_BLEND_PAIRS] = {}; /**< Per-pair from palette. */
    uint8_t ramp_to[MAX_BLEND_PAIRS] = {};   /**< Per-pair to palette. */
    int num_ramps = 0;                       /**< Distinct pair count. */
    BlendWeightFn blend_fn = classic_blend;  /**< Swept crossfade curve. */
    Landing landing; /**< Arrival data exposed to the effect. */
  };

  HS_COLD_MEMBER static const uint8_t *
  dual_closing_palettes(const PolyMesh &dual, const PaletteHandoff &handoff,
                        Arena &scratch) {
    const size_t dual_faces = dual.face_counts.size();
    const size_t primal_faces = dual.vertices.size();
    HS_CHECK(!handoff.prev_face_centroid &&
                 handoff.prev_faces == dual_faces + primal_faces,
             "OpLeg: invalid structural dual handoff");

    uint8_t *from = scratch.allocate_n<uint8_t>(handoff.prev_faces);
    for (size_t f = 0; f < dual_faces; ++f)
      from[f] = handoff.prev_face_palette[primal_faces + f];

    bool *seen = scratch.allocate_n<bool>(primal_faces);
    std::fill_n(seen, primal_faces, false);
    size_t out = dual_faces;
    for (uint16_t vertex : dual.faces) {
      HS_CHECK(vertex < primal_faces, "OpLeg: dual face index out of range");
      if (seen[vertex])
        continue;
      seen[vertex] = true;
      from[out++] = handoff.prev_face_palette[vertex];
    }
    HS_CHECK(out == handoff.prev_faces,
             "OpLeg: structural dual handoff is incomplete");
    return from;
  }

  /**
   * @brief Shared CONWAY_SWEEP construction tail: arrival classification
   * (relaxed form when settling), start centroids, and the palette mappings.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh.
   * @param arena Leg arena for the hoisted state.
   * @param settle Whether the leg settles (relax-slerps) at its arrival end.
   * @param jitterbug Whether the leg is the jitterbug bridge (vertex-orbit
   * faces survive into the node mesh).
   * @details Requires tr.op, tr.reverse, the tr.t_start/t_end and
   * tr.twist_start/twist_end endpoints and tr.settle_frames already set by the
   * calling constructor.
   */
  HS_COLD_MEMBER void init_conway(const PaletteHandoff &handoff,
                                  const BookendClasses &bookend, Arena &arena,
                                  bool settle, bool jitterbug,
                                  bool bridge_provenance) {
    Transients &tr = *buf;
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);
    HS_CHECK(tr.seed_ref, "OpLeg: leg carries no seed mesh");
    const PolyMesh &seed = *tr.seed_ref;
    const uint8_t *forced_from = nullptr;
    if (handoff.correspondence == FaceCorrespondence::DUAL_CLOSING)
      forced_from = dual_closing_palettes(seed, handoff, scratch_arena_a);

    // A closing bridge leg departs the ambo point with its face blocks
    // transposed against the handoff order (medial [P-faces][P-vertices] vs
    // truncate(dual) [D-faces][D-vertices]), so the departed centroids cannot
    // stand in for the start centroids: provenance needs the start mesh's own,
    // built here, before the arrival, so the two never co-reside in scratch.
    const Vector *start_centroid = nullptr;
    if (bridge_provenance && handoff.prev_face_centroid &&
        tr.t_start > tr.t_end) {
      Vector *cen = scratch_arena_a.allocate_n<Vector>(handoff.prev_faces);
      {
        ScratchScope ta(scratch_arena_a);
        ScratchScope tb(scratch_arena_b);
        PolyMesh start = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                                tr.t_start, tr.twist_start);
        HS_CHECK(start.face_counts.size() == handoff.prev_faces,
                 "OpLeg: closing bridge start faces differ from the handoff");
        face_centroids_into(start, cen);
      }
      start_centroid = cen;
    }

    PolyMesh arrival = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                              tr.t_end, tr.twist_end);
    // Classification is hoisted per leg, taken at arrival geometry; a
    // settling forward leg lands on the relaxed form, so classify that.
    PolyMesh *classified = &arrival;
    PolyMesh relaxed_mesh;
    if (settle) {
      if (!tr.reverse) {
        relaxed_mesh =
            MeshOps::relax(arrival, scratch_arena_b, scratch_arena_a, 50);
        classified = &relaxed_mesh;
      } else {
        // Reverse legs un-settle at the start parameter.
        PolyMesh start = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                                tr.t_start, tr.twist_start);
        relaxed_mesh =
            MeshOps::relax(start, scratch_arena_b, scratch_arena_a, 50);
      }
      tr.relaxed.bind(arena, relaxed_mesh.vertices.size());
      tr.relaxed.append_bulk(relaxed_mesh.vertices.data(),
                             relaxed_mesh.vertices.size());
      HS_CHECK(tr.relaxed.size() == arrival.vertices.size(),
               "OpLeg: relax changed the vertex count");
    }
    hoist_arrival_topology(*classified, bookend, arena, tr.topo);

    // Geometric provenance needs the mesh the first frame draws: the
    // relaxed start on a reverse settling leg (settle_alpha == 1 there),
    // the plain op at t_start otherwise.
    PolyMesh start_mesh;
    if (handoff.prev_face_centroid) {
      if (bridge_provenance) {
        // Opening bridge leg: centroids are read only for its newborn faces,
        // where the arrival's own centroids stand in for the birth mesh
        // (rebuilding it would co-reside with the arrival at the scratch
        // peak). Closing legs took theirs from the pre-arrival start build.
        if (!start_centroid)
          start_centroid = face_centroids(*classified, scratch_arena_a);
      } else {
        const PolyMesh *start = &relaxed_mesh;
        if (!(settle && tr.reverse)) {
          start_mesh = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                              tr.t_start, tr.twist_start);
          start = &start_mesh;
        }
        HS_CHECK(start->face_counts.size() == tr.topo.size(),
                 "OpLeg: start face count differs from arrival");
        start_centroid = face_centroids(*start, scratch_arena_a);
      }
    }

    // Emission-order prefix corresponding 1:1 to a node base mesh at the
    // boundary swaps: the seed's face count, plus its vertex-orbit faces on
    // the jitterbug bridge (they survive into the octahedron-end node mesh;
    // only the 12 edge-orbit faces collapse).
    const size_t survivors =
        jitterbug ? tr.seed_faces + seed.vertices.size() : tr.seed_faces;
    build_palette_mapping(tr, *classified, handoff, bookend, arena,
                          start_centroid, survivors, forced_from);
  }

  /**
   * @brief GATED_SWAP construction: both classifications, the across-swap
   * provenance, and both ramp tables.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh.
   * @param arena Leg arena for the hoisted state.
   * @details The compiled face count differs across the swap, so the seed and
   * arrival classifications are held separately and the leg's face-count
   * invariant is checked against whichever side the frame draws.
   */
  HS_COLD_MEMBER void init_gated(const PaletteHandoff &handoff,
                                 const BookendClasses &bookend, Arena &arena) {
    Transients &tr = *buf;
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    MeshOps::classify_faces_by_topology(tr.seed, scratch_arena_a,
                                        scratch_arena_b, arena);
    tr.seed_topo = std::move(tr.seed.topology);
    HS_CHECK(tr.seed_topo.size() == tr.seed.face_counts.size(),
             "OpLeg: gated seed classification is not one class per face");
    HS_CHECK(handoff.prev_faces == tr.seed.face_counts.size(),
             "OpLeg: gated swap departs from a mesh that is not its seed");

    PolyMesh arrival = swap_at(tr, scratch_arena_a, scratch_arena_b);
    MeshOps::classify_faces_by_topology(arrival, scratch_arena_a,
                                        scratch_arena_b, arena);
    tr.topo = std::move(arrival.topology);
    HS_CHECK(tr.topo.size() == arrival.face_counts.size(),
             "OpLeg: gated arrival classification is not one class per face");

    // Provenance across the swap: the children inherit their parent face's
    // (kis) or nearest orbit face's (dual) palette — colour locality.
    uint8_t *from = scratch_arena_a.allocate_n<uint8_t>(tr.topo.size());
    if (tr.swap_op == SwapOp::KIS)
      kis_provenance(tr, arrival, handoff, from);
    else
      dual_provenance(tr, arrival, handoff, from);

    build_palette_mapping(tr, arrival, handoff, bookend, arena,
                          /*start_centroid=*/nullptr,
                          tr.seed.face_counts.size(), from);

    // Seed-side table: the closing half draws at blend weight 0, so only each
    // pair's from palette is ever sampled and an identity pair suffices.
    const size_t seed_faces = tr.seed.face_counts.size();
    tr.seed_face_ramp.bind(arena, seed_faces);
    for (size_t f = 0; f < seed_faces; ++f) {
      const uint8_t pal = handoff.prev_face_palette[f];
      HS_CHECK(pal < PALETTES);
      int ramp = -1;
      for (int r = 0; r < tr.seed_num_ramps; ++r) {
        if (tr.seed_ramp_pal[r] == pal) {
          ramp = r;
          break;
        }
      }
      if (ramp < 0) {
        ramp = tr.seed_num_ramps++;
        tr.seed_ramp_pal[ramp] = pal;
      }
      tr.seed_face_ramp.push_back(static_cast<uint8_t>(ramp));
    }
  }

  /**
   * @brief Steps a gated-swap leg: the seed or the partitioned mesh, then the
   * shared frame tail.
   * @param canvas The canvas passed through to the draw callback.
   * @param frame Clamped frame index; the seed side draws frames [1, gate], the
   * swap and opening side [gate + 1, 2 * gate + 1].
   */
  HS_COLD_MEMBER void step_gated(Canvas &canvas, int frame) {
    Transients &tr = *buf;
    const bool seed_side = frame <= tr.sweep_frames;

    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    PolyMesh swapped;
    {
      HS_PROFILE(hk_gate_swap);
      if (!seed_side)
        swapped = swap_at(tr, scratch_arena_a, scratch_arena_b);
    }
    const PolyMesh &mesh = seed_side ? tr.seed : swapped;

    // Colour holds at the departed palettes across the swap and most of the
    // gate, so the children open in the colour already painted where they land,
    // and converges to the arrival targets only over the final frames.
    finish_frame(canvas, mesh, trailing_blend(frame, duration), seed_side);
  }

  /**
   * @brief Applies the leg's partition operator to its seed.
   * @param tr Leg transients holding the seed and the operator.
   * @param target Arena receiving the output mesh.
   * @param temp Arena holding the operator's transient scratch.
   * @return The partitioned mesh.
   */
  HS_COLD_MEMBER static PolyMesh swap_at(const Transients &tr, Arena &target,
                                         Arena &temp) {
    return tr.swap_op == SwapOp::KIS ? MeshOps::kis(tr.seed, target, temp)
                                     : MeshOps::dual(tr.seed, target, temp);
  }

  /**
   * @brief Per-face from-palettes across a kis swap: every child triangle
   * inherits its parent face's palette.
   * @param tr Leg transients holding the seed.
   * @param arrival Partitioned mesh.
   * @param handoff Departed-mesh provenance.
   * @param from Receives one palette id per arrival face.
   * @details kis emits its parent's children consecutively, one per parent
   * side, so the mapping is exact.
   */
  HS_COLD_MEMBER static void kis_provenance(const Transients &tr,
                                            const PolyMesh &arrival,
                                            const PaletteHandoff &handoff,
                                            uint8_t *from) {
    size_t f = 0;
    for (size_t p = 0; p < tr.seed.face_counts.size(); ++p) {
      for (int k = 0; k < tr.seed.face_counts[p]; ++k) {
        HS_CHECK(f < arrival.face_counts.size(), "OpLeg: kis child overflow");
        from[f++] = handoff.prev_face_palette[p];
      }
    }
    HS_CHECK(f == arrival.face_counts.size(),
             "OpLeg: kis child count differs from the seed's side count");
  }

  /**
   * @brief Per-face from-palettes across a dual swap: each dual face takes the
   * palette of the orbit face whose centroid is nearest the source vertex it
   * opens on.
   * @param tr Leg transients holding the seed.
   * @param arrival Dual mesh.
   * @param handoff Departed-mesh provenance.
   * @param from Receives one palette id per arrival face.
   * @details dual's vertices are its source faces' normalized centroids indexed
   * by source face, and a dual face's vertex list is exactly its source
   * vertex's face orbit, so the orbit needs no second walk. Colour locality is
   * the goal: pixel identity is not available across a partition (spec 3.1).
   */
  HS_COLD_MEMBER static void dual_provenance(const Transients &tr,
                                             const PolyMesh &arrival,
                                             const PaletteHandoff &handoff,
                                             uint8_t *from) {
    size_t off = 0;
    for (size_t f = 0; f < arrival.face_counts.size(); ++f) {
      const int n = arrival.face_counts[f];
      const Vector c =
          normalized_or(face_vertex_sum(arrival.vertices.data(),
                                        arrival.faces.data(), off, n),
                        arrival.vertices[arrival.faces[off]]);

      // The orbit's own source vertex: the dual face's centroid lies in that
      // vertex's cell, so the nearest seed vertex is it.
      Vector v = tr.seed.vertices[0];
      float best_dot = -2.0f;
      for (size_t i = 0; i < tr.seed.vertices.size(); ++i) {
        const float d = dot(tr.seed.vertices[i], c);
        if (d > best_dot) {
          best_dot = d;
          v = tr.seed.vertices[i];
        }
      }

      uint16_t best_face = arrival.faces[off];
      float best_sq = 1e9f;
      for (int k = 0; k < n; ++k) {
        const uint16_t j = arrival.faces[off + k];
        const Vector d = arrival.vertices[j] - v;
        const float dsq = dot(d, d);
        if (dsq < best_sq) {
          best_sq = dsq;
          best_face = j;
        }
      }
      HS_CHECK(best_face < handoff.prev_faces,
               "OpLeg: dual face index outside the departed face list");
      from[f] = handoff.prev_face_palette[best_face];
      off += n;
    }
  }

  /**
   * @brief Runs the edge's operator on the seed at one parameter value.
   */
  HS_COLD_MEMBER static PolyMesh run_op(ConwayGraph::MorphOp op,
                                        const PolyMesh &seed, Arena &target,
                                        Arena &temp, float t, float twist) {
    switch (op) {
    case ConwayGraph::MorphOp::TRUNCATE:
      // A far-side leg sweeps through t = 0.5; nudge that one sample off the
      // ambo short-circuit so the frame keeps the truncate topology (2E
      // vertices) instead of popping to ambo (V vertices).
      return MeshOps::truncate(seed, target, temp,
                               ConwayGraph::truncate_off_pinch(t));
    case ConwayGraph::MorphOp::EXPAND:
      return MeshOps::expand(seed, target, temp, t);
    case ConwayGraph::MorphOp::CHAMFER:
      return MeshOps::chamfer(seed, target, temp, t);
    default:
      return MeshOps::snub(seed, target, temp, t, twist);
    }
  }

  /** @brief Hoists arrival classes, reusing an exact bookend classification. */
  HS_COLD_MEMBER static void
  hoist_arrival_topology(PolyMesh &arrival, const BookendClasses &bookend,
                         Arena &arena, ArenaVector<uint16_t> &topology) {
    const size_t faces = arrival.face_counts.size();
    if (bookend.topology && bookend.faces == faces) {
      topology.bind(arena, faces);
      topology.append_bulk(bookend.topology, faces);
      return;
    }
    MeshOps::classify_faces_by_topology(arrival, scratch_arena_a,
                                        scratch_arena_b, arena);
    topology = std::move(arrival.topology);
    HS_CHECK(topology.size() == faces);
  }

  /**
   * @brief Kind-agnostic frame tail: compile the swept mesh, attach the
   * hoisted classification, pre-blend the palette ramps at the crossfade weight,
   * draw.
   * @param canvas The canvas passed through to the draw callback.
   * @param swept This frame's swept mesh (scratch-backed).
   * @param w Crossfade weight in [0, 1] the caller already resolved (the leg's
   * blend_fn for the swept kinds, trailing_blend for the gate).
   * @param seed_side Draw the gate's seed-side tables (GATED_SWAP only): the
   * seed classification and its identity ramp table, which w == 0 leaves at the
   * departed palettes.
   */
  HS_COLD_MEMBER void finish_frame(Canvas &canvas, const PolyMesh &swept,
                                   float w, bool seed_side = false) {
    Transients &tr = *buf;
    const ArenaVector<uint16_t> &topo = seed_side ? tr.seed_topo : tr.topo;
    const ArenaVector<uint8_t> &face_ramp =
        seed_side ? tr.seed_face_ramp : tr.face_ramp;

    MeshState compiled;
    {
      HS_PROFILE(hk_conway_compile);
      MeshOps::compile(swept, compiled, scratch_arena_a, scratch_arena_b);
    }
    HS_CHECK(compiled.face_counts.size() == topo.size(),
             "OpLeg: sweep changed the compiled face count");

    const int num_ramps = seed_side ? tr.seed_num_ramps : tr.num_ramps;
    BakedPalette *ramps = scratch_arena_b.allocate_n<BakedPalette>(num_ramps);
    for (int r = 0; r < num_ramps; ++r) {
      new (&ramps[r]) BakedPalette();
      if (seed_side) {
        ramps[r] = tr.bank->entries[tr.seed_ramp_pal[r]];
        continue;
      }
      const BakedPalette &from = tr.bank->entries[tr.ramp_from[r]];
      const BakedPalette &to = tr.bank->entries[tr.ramp_to[r]];
      if (tr.ramp_from[r] == tr.ramp_to[r])
        ramps[r] = to;
      else
        bake_palette_blend(ramps[r], scratch_arena_b, from, to, w);
    }

    Shading sh{ramps, face_ramp.data(), face_ramp.size()};
    draw_fn(canvas, compiled, sh);
  }

  /**
   * @brief Fills a swept mesh's vertices by slerping between two endpoints.
   * @param out Output mesh whose vertices are bound and filled.
   * @param arena Arena backing the vertex array.
   * @param from Opening vertices (indices [shared, n) only).
   * @param to Arrival vertices.
   * @param n Vertex count.
   * @param shared Prefix of @p to copied verbatim (vertices that never move).
   * @param k Slerp fraction in [0, 1]; 1 copies @p to verbatim so the closing
   * bookend swap is bitwise, not 1 ULP off.
   */
  HS_COLD_MEMBER static void slerp_vertices(PolyMesh &out, Arena &arena,
                                            const Vector *from,
                                            const Vector *to, size_t n,
                                            size_t shared, float k) {
    out.vertices.bind(arena, n);
    if (k >= 1.0f) {
      out.vertices.append_bulk(to, n);
      return;
    }
    out.vertices.append_bulk(to, shared);
    for (size_t i = shared; i < n; ++i)
      out.vertices.push_back(slerp(from[i], to[i], k));
  }

  /**
   * @brief Copies a leg's constant topology into the swept mesh.
   * @param out Output mesh whose face arrays are bound and filled.
   * @param arena Arena backing the face arrays.
   * @param face_counts Per-face side counts.
   * @param faces Flat per-face vertex index list.
   */
  HS_COLD_MEMBER static void
  copy_topology(PolyMesh &out, Arena &arena,
                const ArenaVector<uint8_t> &face_counts,
                const ArenaVector<uint16_t> &faces) {
    out.face_counts.bind(arena, face_counts.size());
    out.face_counts.append_bulk(face_counts.data(), face_counts.size());
    out.faces.bind(arena, faces.size());
    out.faces.append_bulk(faces.data(), faces.size());
  }

  /**
   * @brief Deep-copies a seed's geometry, leaving its class ids behind.
   * @param src Seed mesh.
   * @param dst Destination mesh, populated in place.
   * @param arena Arena backing the destination arrays.
   */
  HS_COLD_MEMBER static void clone_geometry(const PolyMesh &src, PolyMesh &dst,
                                            Arena &arena) {
    dst.vertices.bind(arena, src.vertices.size());
    dst.vertices.append_bulk(src.vertices.data(), src.vertices.size());
    copy_topology(dst, arena, src.face_counts, src.faces);
  }

  /**
   * @brief Builds a hankin leg's swept mesh at one slerp fraction.
   * @param tr Leg transients holding the baked topology and both endpoints.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1].
   */
  HS_COLD_MEMBER static void hankin_at(const Transients &tr, PolyMesh &out,
                                       Arena &arena, float k) {
    const CompiledHankin &hk = tr.hankin;
    const size_t statics = hk.static_vertices.size();
    const size_t dyn = tr.hk_final.size();
    out.vertices.bind(arena, statics + dyn);
    out.vertices.append_bulk(hk.static_vertices.data(), statics);
    if (k >= 1.0f) {
      for (size_t i = 0; i < dyn; ++i)
        out.vertices.push_back(tr.hk_final[i].decode().normalized());
    } else {
      for (size_t i = 0; i < dyn; ++i) {
        const Vector corner = hk.corner(hk.dynamic_instructions[i].v_corner);
        out.vertices.push_back(
            slerp(corner.normalized(), tr.hk_final[i].decode(), k));
      }
    }
    copy_topology(out, arena, hk.face_counts, hk.faces);
  }

  /**
   * @brief Builds a relax leg's swept mesh at one slerp fraction.
   * @param tr Leg transients holding the seed and its relaxed endpoint.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1].
   */
  HS_COLD_MEMBER static void relax_at(const Transients &tr, PolyMesh &out,
                                      Arena &arena, float k) {
    slerp_vertices(out, arena, tr.seed.vertices.data(), tr.relaxed.data(),
                   tr.relaxed.size(), 0, k);
    copy_topology(out, arena, tr.seed.face_counts, tr.seed.faces);
  }

  /**
   * @brief Builds a medial leg's swept mesh at one slerp fraction, decoding the
   * snorm16-packed endpoints per vertex.
   * @param tr Leg transients holding the packed a_e/b_e sets and connectivity.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1]; slerp() renormalizes the near-unit
   * decoded inputs, so the quantization is absorbed into the direction blend.
   */
  HS_COLD_MEMBER static void medial_at(const Transients &tr, PolyMesh &out,
                                       Arena &arena, float k) {
    const size_t n = tr.medial_a.size();
    out.vertices.bind(arena, n);
    if (k >= 1.0f) {
      for (size_t i = 0; i < n; ++i)
        out.vertices.push_back(tr.medial_b[i].decode().normalized());
    } else {
      for (size_t i = 0; i < n; ++i)
        out.vertices.push_back(
            slerp(tr.medial_a[i].decode(), tr.medial_b[i].decode(), k));
    }
    copy_topology(out, arena, tr.seed.face_counts, tr.seed.faces);
  }

  /** Max distance from a start-parameter face centroid to its departed
   * counterpart: bounds the T_EPS-scale swap displacement while staying under
   * half the face-centroid spacing of the largest node (~0.37 chord at 92
   * faces). */
  static constexpr float PROVENANCE_TOL_SQ = 0.15f * 0.15f;

  /**
   * @brief Unit-sphere vertex-average centroid of every face.
   * @param m Mesh whose faces are reduced.
   * @param arena Arena receiving the centroid array.
   * @return One unit vector per face, in emission order.
   */
  HS_COLD_MEMBER static const Vector *face_centroids(const PolyMesh &m,
                                                     Arena &arena) {
    Vector *out = arena.allocate_n<Vector>(m.face_counts.size());
    face_centroids_into(m, out);
    return out;
  }

  /**
   * @brief Index of the departed face nearest to a start-face centroid.
   * @param c Unit centroid of the swept face at the start parameter.
   * @param handoff Departed-node provenance (centroids non-null).
   * @return Index into the handoff arrays.
   */
  HS_COLD_MEMBER static size_t
  nearest_prev_face(const Vector &c, const PaletteHandoff &handoff) {
    size_t best = 0;
    float best_d = 1e9f;
    for (size_t j = 0; j < handoff.prev_faces; ++j) {
      const Vector d = c - handoff.prev_face_centroid[j];
      const float dsq = dot(d, d);
      if (dsq < best_d) {
        best_d = dsq;
        best = j;
      }
    }
    return best;
  }

  /**
   * @brief Crossfade weight at fraction p of the leg (classic_blend's curve).
   */
  static float blend_weight(float p) {
    constexpr float IN = 0.2f, OUT = 0.8f;
    if (p <= IN)
      return 0.0f;
    if (p >= OUT)
      return 1.0f;
    float u = (p - IN) / (OUT - IN);
    return u * u * (3.0f - 2.0f * u);
  }

  // always_inline, not a plain helper: an out-of-line copy inherits no cold
  // attribute from its HS_COLD_MEMBER caller and would land in ITCM.
  __attribute__((always_inline)) static const uint16_t *resolve_target_topology(
      Transients &tr, const PolyMesh &arrival, const PaletteHandoff &handoff,
      const BookendClasses &bookend, Arena &arena, size_t survivors) {
    const size_t total = tr.topo.size();
    HS_CHECK(!bookend.topology || bookend.faces == total ||
                 bookend.faces == survivors,
             "OpLeg: bookend face count matches neither mapping");
    if (!bookend.topology) {
      tr.target_topo.clear();
      return tr.topo.data();
    }

    const size_t landed = bookend.faces;
    const bool structural_closing =
        handoff.correspondence == FaceCorrespondence::DUAL_CLOSING;
    // A closing leg's corner births are structural: corner face k is born on
    // the k-th first-seen vertex of the dual seed's face walk (the order
    // dual_closing_palettes maps from-palettes by), so the face whose walk
    // first reaches that vertex is the corner's host — no centroid search.
    int *corner_host = nullptr;
    if (structural_closing && landed < total) {
      const PolyMesh &dual = *tr.seed_ref;
      const size_t corners = total - landed;
      HS_CHECK(landed == dual.face_counts.size() &&
                   corners == dual.vertices.size(),
               "OpLeg: closing-leg blocks differ from the dual seed");
      corner_host = scratch_arena_a.allocate_n<int>(corners);
      bool *seen = scratch_arena_a.allocate_n<bool>(corners);
      std::fill_n(seen, corners, false);
      size_t out = 0, off = 0;
      for (size_t fi = 0; fi < landed; ++fi) {
        for (int j = 0; j < dual.face_counts[fi]; ++j) {
          const uint16_t v = dual.faces[off + j];
          if (!seen[v]) {
            seen[v] = true;
            corner_host[out++] = static_cast<int>(fi);
          }
        }
        off += dual.face_counts[fi];
      }
      HS_CHECK(out == corners, "OpLeg: dual corner hosts incomplete");
    }
    const Vector *arrival_centroid =
        landed < total && !structural_closing
            ? face_centroids(arrival, scratch_arena_a)
            : nullptr;
    tr.target_topo.bind(arena, total);
    for (size_t f = 0; f < total; ++f) {
      if (f < landed) {
        tr.target_topo.push_back(bookend.topology[f]);
        continue;
      }
      if (structural_closing) {
        tr.target_topo.push_back(bookend.topology[corner_host[f - landed]]);
        continue;
      }
      size_t host = 0;
      float best_d = 1e9f;
      for (size_t j = 0; j < landed; ++j) {
        const Vector d = arrival_centroid[f] - arrival_centroid[j];
        const float dsq = dot(d, d);
        if (dsq < best_d) {
          best_d = dsq;
          host = j;
        }
      }
      tr.target_topo.push_back(bookend.topology[host]);
    }
    return tr.target_topo.data();
  }

  // always_inline, not a plain helper: an out-of-line copy inherits no cold
  // attribute from its HS_COLD_MEMBER caller and would land in ITCM.
  __attribute__((always_inline)) static uint8_t
  intern_palette_ramp(Transients &tr, uint8_t from, uint8_t to) {
    for (int r = 0; r < tr.num_ramps; ++r)
      if (tr.ramp_from[r] == from && tr.ramp_to[r] == to)
        return static_cast<uint8_t>(r);

    HS_CHECK(tr.num_ramps < MAX_BLEND_PAIRS,
             "OpLeg: distinct palette pairs exceed the blend budget");
    const int ramp = tr.num_ramps++;
    tr.ramp_from[ramp] = from;
    tr.ramp_to[ramp] = to;
    return static_cast<uint8_t>(ramp);
  }

  /**
   * @brief Builds the per-face from-palettes (leg-swap mapping, spec 2.5), the
   * shuffled target assignment (spec 2.6), and the distinct ramp-pair table.
   * @param tr Leg transients being populated.
   * @param arrival Classified arrival mesh (for face counts).
   * @param handoff Departed-node provenance.
   * @param bookend Arrival-node bookend grouping the targets key on.
   * @param arena Leg arena for the face -> ramp table.
   * @param start_centroid Unit centroid per swept face at the start parameter,
   * or nullptr for the emission-order mapping.
   * @param survivors Emission-order face-prefix length corresponding 1:1 to a
   * node base mesh at the boundary swaps (the seed face count; plus the
   * vertex-orbit faces on the jitterbug bridge, whose octahedron-end node
   * mesh keeps them).
   * @param forced_from Per-arrival-face from-palette, or nullptr to derive one.
   * A partition op has neither a centroid nor an emission-order
   * correspondence, so its leg computes provenance itself and hands it in.
   * @details With centroids, provenance is geometric: a face inherits the
   * palette of the departed face it overlies at the start parameter. On a
   * full-correspondence departure (prev_faces == total: 0.5-end swaps, dual
   * swaps, regenerated seeds) every face maps by nearest departed centroid —
   * a checked bijection that assumes nothing about emission order and keeps
   * per-side-count class splits, both of which the legacy mappings break. On
   * a node-prefix departure (prev_faces == survivors) the prefix keeps the
   * exact emission identity, and each newborn class inherits its first face's
   * nearest departed palette, so T_EPS-wide births open in the underlying
   * face's colors instead of popping in as target-colored slivers.
   */
  HS_COLD_MEMBER void build_palette_mapping(
      Transients &tr, const PolyMesh &arrival, const PaletteHandoff &handoff,
      const BookendClasses &bookend, Arena &arena, const Vector *start_centroid,
      size_t survivors, const uint8_t *forced_from = nullptr) {
    const size_t total = tr.topo.size();
    const size_t primary = tr.seed_faces;
    tr.landing.faces = total;
    tr.landing.primary_faces = primary;

    // Target classes: the bookend grouping where a face survives the swap
    // (emission-order identity). A face that collapses to a T_EPS sliver has
    // no counterpart there, so it takes the bookend class of the arrival face
    // it collapses onto — the mirror of the newborn from-palette rule below,
    // which is what makes the leg's crossfade symmetric: a sliver closes into
    // its host's color instead of freezing in an unrelated target color.
    const uint16_t *target_topo = resolve_target_topology(
        tr, arrival, handoff, bookend, arena, survivors);
    tr.landing.topology = target_topo;

    for (int i = 0; i < PALETTES; ++i)
      tr.landing.to_palette[i] = static_cast<uint8_t>(i);
    hs::shuffle(tr.landing.to_palette.begin(), tr.landing.to_palette.end());

    // Either the whole face list corresponds (departing a mid-parameter node)
    // or an emission prefix does: the survivor prefix departing a node form,
    // the seed's primaries departing the seed form (the remaining faces are
    // births).
    HS_CHECK(handoff.prev_faces == total || handoff.prev_faces == survivors ||
                 handoff.prev_faces == primary,
             "OpLeg: handoff face count matches no mapping");

    // Full-correspondence geometric mapping: nearest departed centroid,
    // pinned as a bijection within tolerance. The bijection mark is scoped to
    // that mapping: the prefix mapping a chain leg takes never reads it.
    const bool full_correspondence =
        start_centroid && handoff.prev_faces == total;
    bool *prev_used = nullptr;
    if (full_correspondence) {
      prev_used = scratch_arena_a.allocate_n<bool>(handoff.prev_faces);
      std::fill_n(prev_used, handoff.prev_faces, false);
    }
    // Newborn classes inherit one representative from-palette (first face of
    // the class), keeping the distinct-pair count at the legacy bound.
    int newborn_from[PALETTES];
    for (int i = 0; i < PALETTES; ++i)
      newborn_from[i] = -1;

    uint8_t *from_palette = arena.allocate_n<uint8_t>(total);
    tr.landing.from_palette = from_palette;

    tr.face_ramp.bind(arena, total);
    for (size_t f = 0; f < total; ++f) {
      const uint8_t to =
          tr.landing
              .to_palette[wrap(static_cast<int>(target_topo[f]), PALETTES)];
      uint8_t from = to; // fallback: newborn faces skip the crossfade
      if (forced_from) {
        from = forced_from[f];
      } else if (handoff.correspondence == FaceCorrespondence::IDENTITY) {
        HS_CHECK(handoff.prev_faces == total,
                 "OpLeg: identity handoff face count differs");
        from = handoff.prev_face_palette[f];
      } else if (full_correspondence) {
        const size_t j = nearest_prev_face(start_centroid[f], handoff);
        const Vector d = start_centroid[f] - handoff.prev_face_centroid[j];
        HS_CHECK(!prev_used[j] && dot(d, d) < PROVENANCE_TOL_SQ,
                 "OpLeg: start face has no unique departed counterpart");
        prev_used[j] = true;
        from = handoff.prev_face_palette[j];
      } else if (start_centroid) { // prev_faces == survivors
        if (f < handoff.prev_faces) {
          from = handoff.prev_face_palette[f];
        } else {
          const int slot = wrap(static_cast<int>(target_topo[f]), PALETTES);
          if (newborn_from[slot] < 0)
            newborn_from[slot] = handoff.prev_face_palette[nearest_prev_face(
                start_centroid[f], handoff)];
          from = static_cast<uint8_t>(newborn_from[slot]);
        }
      } else if (f < handoff.prev_faces) {
        from = handoff.prev_face_palette[f];
      }
      HS_CHECK(from < PALETTES && to < PALETTES,
               "OpLeg: crossfade ramp endpoint outside the palette bank");
      from_palette[f] = from;

      tr.face_ramp.push_back(intern_palette_ramp(tr, from, to));
    }
    tr.landing.blend_pairs = tr.num_ramps;
  }

  /**
   * @brief Allocates the leg's Transients and stamps what check_alive() tests.
   * @param arena Leg arena backing the leg's state.
   * @return The fresh Transients.
   */
  Transients &bind_transients(Arena &arena) {
    buf = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();
    leg_arena = &arena;
    live_end = arena.get_offset();
#ifndef NDEBUG
    birth_generation = arena.get_generation();
#endif
    return *buf;
  }

  /**
   * @brief Traps if the leg arena was reclaimed while the leg is still live.
   * @details A rewind is caught in every build by the watermark; a reset()
   * refills the arena to any offset, so only the debug generation stamp
   * (which a rewind does not bump) detects it exactly.
   */
  void check_alive() const {
    HS_CHECK(leg_arena->get_offset() >= live_end,
             "OpLeg: leg arena rewound under a live leg");
#ifndef NDEBUG
    HS_CHECK(leg_arena->get_generation() == birth_generation,
             "OpLeg: leg arena reset under a live leg");
#endif
  }

  Transients *buf;  /**< Pointer to arena-allocated leg state. */
  Arena *leg_arena; /**< Arena backing buf and the Transients vectors. */
  size_t live_end;  /**< Arena offset just past the Transients block. */
#ifndef NDEBUG
  uint32_t birth_generation = 0; /**< Leg-arena generation at construction. */
#endif
  EasingFn easing_fn;  /**< Easing applied to the sweep parameter. */
  MorphDrawFn draw_fn; /**< Per-frame draw callback. */
};

} // namespace Animation
