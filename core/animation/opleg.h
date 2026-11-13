/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstring>

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @file opleg.h
 * @brief Animation fragment: OpLeg Conway-chain morph legs.
 */

#include "core/color/palettes.h"
#include "core/mesh/conway.h"
#include "core/mesh/conway_graph.h"
#include "core/mesh/hankin.h"
#include "core/animation/opleg/specs.h"

namespace Animation {

/**
 * @brief Animates a Conway sweep, hankin sweep, relax slerp, medial slerp, or
 * gated partition swap.
 * @details Bulk state lives in an arena-allocated Transients that no
 * destructor reclaims; the caller compacts the arena between legs.
 */
class OpLeg : public AnimationBase<OpLeg> {
public:
#if HS_ENABLE_TEST_HOOKS
  /** @brief Force classification before resolving a partial bookend. */
  inline static bool test_classify_partial_bookend = false;
  /** @brief Faces classified while constructing legs. */
  inline static size_t test_classified_faces = 0;
  /** @brief Corners visited by construction-time face classification. */
  inline static size_t test_classified_corners = 0;
#endif
  /// Palette slots a leg maps between; one per baked palette.
  static constexpr int PALETTES = BakedPaletteBank::N;
  /** Capacity of the intern table: the full (from, to) pair space. Each
   * distinct non-identity pair also bakes a LUT into scratch_arena_b each
   * frame, so scratch capacity is the per-frame ceiling. */
  static constexpr int MAX_BLEND_PAIRS = PALETTES * PALETTES;
  static_assert(MAX_BLEND_PAIRS <= 256,
                "interned pair indices are uint8_t; every index below "
                "MAX_BLEND_PAIRS must fit in it");

  /** Leg kind selected by the constructor. */
  using LegKind = detail::LegKind;

  /** Partition operator a GATED_SWAP leg swaps to. */
  using SwapOp = detail::SwapOp;

  /**
   * @brief Graph-edge Conway sweep: the edge's single op swept between the
   * edge's endpoint parameters.
   */
  using EdgeSweepSpec = detail::EdgeSweepSpec;

  /** @brief Immutable classification and newborn hosts of an opening dual bridge. */
  struct BridgePreparation {
    ArenaVector<uint16_t> topology; /**< Arrival classes in sweep order. */
    std::array<uint16_t, PALETTES>
        hosts{}; /**< Departed face per newborn slot. */

    /** @brief Copies prepared bridge data between arenas.
     * @param src Prepared source.
     * @param dst Destination.
     * @param arena Destination storage. */
    HS_FLASH_MEMBER static void clone(const BridgePreparation &src,
                                      BridgePreparation &dst, Arena &arena) {
      dst.hosts = src.hosts;
      if (!src.topology.empty()) {
        dst.topology.bind(arena, src.topology.size());
        dst.topology.append_bulk(src.topology.data(), src.topology.size());
      }
    }
  };

  /** @brief Precomputes an opening dual bridge without consuming random state.
   * @param seed Immutable bridge seed.
   * @param arena Storage for the prepared classes.
   * @return Arrival classes and departed hosts for newborn palette slots. */
  inline HS_FLASH_MEMBER static BridgePreparation
  prepare_bridge(const PolyMesh &seed, Arena &arena);

  /**
   * @brief Recipe-step Conway sweep: one primitive op swept t_start -> t_end
   * on a fixed seed, no graph edge.
   */
  struct ParamSweepSpec {
    ConwayGraph::MorphOp op =
        ConwayGraph::MorphOp::TRUNCATE; /**< Swept operator. */
    float t_start = 0.0f; /**< Sweep parameter at frame 0; T_EPS floor except
                                truncate's min(T_EPS, larger endpoint *
                                TRUNCATE_BIRTH_FRAC) birth floor. */
    float t_end = 0.0f;   /**< Arrival parameter; truncate is capped below
                                ambo for near-side legs, or at
                                T_TRUNCATE_FAR_MAX for far-side legs. */
    float twist_start = 0.0f; /**< Snub twist at frame 0. */
    float twist_end = 0.0f;   /**< Arrival snub twist. */
    int sweep_frames = 0;     /**< Operator-sweep frames (N). */
    bool bridge_provenance =
        false; /**< Dual-bridge leg: geometric provenance takes its start
                  centroids from the leg's own start (closing) or arrival
                  (opening) mesh, since the departed face blocks are
                  transposed against the handoff. */
    const BridgePreparation *prepared =
        nullptr; /**< Opening bridge cache; hosts use the seed centroids. */
  };

  /** @brief Conway sweep source; cloned unless explicitly borrowed. */
  class SweepSeed {
  public:
    /**
     * @brief Seed that clones `mesh` when the leg starts.
     * @param mesh Sweep source mesh.
     */
    SweepSeed(const PolyMesh &mesh) : mesh(mesh) {}
    SweepSeed(PolyMesh &&) = delete;
    SweepSeed(const PolyMesh &&) = delete;

    /**
     * @brief Borrows a mesh that remains unmoved until the leg completes.
     * @param mesh Sweep source mesh, used in place.
     * @return Seed that borrows `mesh`.
     */
    static SweepSeed borrow(const PolyMesh &mesh) {
      return SweepSeed(mesh, true);
    }
    static SweepSeed borrow(PolyMesh &&) = delete;
    static SweepSeed borrow(const PolyMesh &&) = delete;

  private:
    friend class OpLeg;
    const PolyMesh &mesh;
    bool borrowed = false;
    SweepSeed(const PolyMesh &mesh, bool borrowed)
        : mesh(mesh), borrowed(borrowed) {}
  };

  /** @brief Hankin contact-angle sweep on a fixed seed. */
  using HankinSweepSpec = detail::HankinSweepSpec;

  /**
   * @brief Relax leg: every vertex slerps from its seed position to its
   * relaxed one.
   */
  using RelaxSpec = detail::RelaxSpec;

  /** @brief Medial-slerp leg (the Conway-dual bridge's middle leg). */
  using MedialSpec = detail::MedialSpec;

  /** @brief Reconcile leg closing a smooth kis/needle path. */
  using ReconcileSpec = detail::ReconcileSpec;

  /** @brief Gated-swap leg: a partition op with no sweep. */
  using GatedSwapSpec = detail::GatedSwapSpec;

  /** Floor on a hankin leg's arrival contact angle (the T_EPS analog): the
   * arrival angle divides the opening angle into the leg's opening slerp
   * fraction. */
  static constexpr float THETA_EPS = 0.02f;

  /** Slerp-fraction floor for a hankin leg's opening frame: at fraction 0
   * every rosette face has zero area and compile would drop it. */
  static constexpr float K_EPS = 0.005f;

  /** Squared maximum distance from a start-parameter face centroid to its
   * departed counterpart (0.15 chord); uniqueness is enforced by the used-face
   * check. */
  static constexpr float PROVENANCE_TOL_SQ = 0.15f * 0.15f;

  /** Default trailing blend window: colour holds at `from` until this many
   * frames remain, then ramps to `to` at the arrival. */
  static constexpr int TRAILING_BLEND_FRAMES = 4;

  /**
   * @brief Per-frame shading and previous geometry handed to the draw
   * callback.
   * @details Scratch-backed; valid for the current frame only.
   */
  struct Shading {
    const BakedPalette *ramps; /**< One blended LUT per (from, to) pair. */
    const uint8_t *face_ramp;  /**< Face index -> ramp index. */
    size_t faces;              /**< Face count (bounds face_ramp). */
    /** Source LUTs of each pair's endpoints; null when the frame's ramps
     * are not crossfades. */
    const BakedPaletteBank *bank = nullptr;
    const uint8_t *ramp_from = nullptr; /**< Per-pair from palette. */
    const uint8_t *ramp_to = nullptr;   /**< Per-pair to palette. */
    float weight = 0.0f;                /**< Crossfade weight `ramps` carry. */
    float previous_weight = 0.0f; /**< Crossfade weight the leg drew a frame
                                     ago, or `weight` after a reset. */
    /** The drawn mesh's vertices at the leg frame drawn last, index for
     * index, or null when the leg drew no other frame before this one or
     * moves through `leg`; the draw may overwrite them. */
    math::Vector *previous_vertices = nullptr;
    /** Leg whose drawn mesh moves vertex by vertex from the frame drawn last,
     * as vertex_at() evaluates it; null otherwise. */
    const OpLeg *leg = nullptr;
    float drawn_frame = 0.0f; /**< Leg frame drawn last, read with `leg`. */
    float frame = 0.0f;       /**< Leg frame drawn now, read with `leg`. */

    /**
     * @brief A vertex of the drawn mesh at a fraction of the frame.
     * @param i Vertex index.
     * @param u Frame fraction: 0 at the frame drawn last, 1 now.
     * @return The vertex on the leg's own path; requires `leg`.
     */
    HS_FLASH_MEMBER math::Vector vertex_at(size_t i, float u) const {
      return leg->swept_vertex(i, drawn_frame + (frame - drawn_frame) * u);
    }

    /**
     * @brief The leg's sweep parameter at a fraction of the frame.
     * @param u Frame fraction: 0 at the frame drawn last, 1 now.
     * @return The parameter vertex_at_parameter() reads; requires `leg`.
     */
    HS_FLASH_MEMBER float parameter_at(float u) const {
      return leg->sweep_parameter(drawn_frame + (frame - drawn_frame) * u);
    }

    /**
     * @brief A vertex of the drawn mesh at a sweep parameter.
     * @param i Vertex index.
     * @param k Sweep parameter, as parameter_at() gives it.
     * @return The vertex on the leg's own path; requires `leg`.
     */
    math::Vector vertex_at_parameter(size_t i, float k) const {
      return leg->swept_vertex_at(i, k);
    }

    /**
     * @brief The face's pre-blended ramp.
     * @param face Face index.
     * @return The blended LUT to sample.
     */
    const BakedPalette &ramp_for(size_t face) const {
      HS_CHECK(face < faces, "OpLeg::Shading: ramp face out of range");
      return ramps[face_ramp[face]];
    }

    /**
     * @brief Gives a face shader the face's crossfade at both frame endpoints.
     * @tparam ShaderT Shader with FacePaletteShader's `set_blend`.
     * @param shader Shader whose palette is already the face's ramp.
     * @param face Face index.
     * @details A face whose pair holds one palette, or a frame without
     * crossfades, keeps the ramp alone.
     */
    template <typename ShaderT>
    HS_FLASH_MEMBER void blend_face(ShaderT &shader, size_t face) const {
      HS_CHECK(face < faces, "OpLeg::Shading: blend face out of range");
      const uint8_t ramp = face_ramp[face];
      if (bank == nullptr || ramp_from[ramp] == ramp_to[ramp])
        return;
      shader.set_blend(&bank->entries[ramp_from[ramp]].view(),
                       &bank->entries[ramp_to[ramp]].view(), weight,
                       previous_weight);
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
   * for a leg frame in [0, duration] (0 is the paused initial state) over the
   * whole leg duration (sweep plus settle). Defaults to classic_blend.
   */
  using BlendWeightFn = float (*)(int frame, int duration);

  /** @brief Face-order relationship between departed and swept meshes. */
  using FaceCorrespondence = detail::FaceCorrespondence;

  /**
   * @brief Palette provenance of the departed node.
   * @details prev_face_palette describes the node base mesh the leg departs
   * from, in emission order; read only during construction. Supplying
   * prev_face_centroid makes the mapping geometric.
   */
  using PaletteHandoff = detail::PaletteHandoff;

  /**
   * @brief Bookend grouping of the arrival node.
   * @details topology[f] is the class the effect displays arrival face f
   * with at the closing bookend. Color targets key on it, so faces it merges
   * converge to one color by w = 1. A null topology keys targets on the swept
   * classification instead. Palette-preserving legs retain departed colors.
   */
  using BookendClasses = detail::BookendClasses;

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
        nullptr; /**< Per-swept-face FROM palette id (arena-backed). */
    bool preserve_palette = false; /**< Retain each face's departed palette. */
    int blend_pairs = 0; /**< Distinct (from, to) ramp pairs the leg carries. */

    /**
     * @brief Palette id face f displays once the leg arrives (w = 1).
     * @param f Swept face index.
     * @return Landed palette index.
     * @details Palette-preserving legs retain departed colors; other legs
     * land on their target-class palette.
     */
    uint8_t landed_palette(size_t f) const {
      HS_CHECK(topology && f < faces,
               "OpLeg::Landing: landed_palette face out of range");
      return preserve_palette
                 ? from_palette[f]
                 : to_palette[MeshPaletteBank::slot_of(topology[f])];
    }
    /**
     * @brief Copies the arrived palette ids for a prefix of swept faces.
     * @param out Destination with room for n ids.
     * @param n Number of swept faces to copy, at most faces.
     */
    void copy_landed_palettes(uint8_t *out, size_t n) const {
      for (size_t f = 0; f < n; ++f)
        out[f] = landed_palette(f);
    }
    const PolyMesh *arrival_topology =
        nullptr; /**< Fixed connectivity of a packed arrival endpoint. */
    const math::Snorm3 *arrival_point =
        nullptr;               /**< Packed arrival endpoint vertices. */
    size_t arrival_points = 0; /**< Length of arrival_point. */
    const CompiledHankin *hankin =
        nullptr; /**< Baked arrival topology (HANKIN_SWEEP legs only, null
                    otherwise); with star_point it rebuilds the arrival mesh
                    through arrival_mesh(). */
    const math::Snorm3 *star_point =
        nullptr;            /**< Arrival star points (snorm16-packed). */
    size_t star_points = 0; /**< Star-point count. */
  };

  /**
   * @brief Constructs a graph-edge sweep leg: clones the seed, computes the
   * arrival classification (relaxed form on a forward settling leg), and
   * builds the palette mappings.
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
  inline HS_FLASH_MEMBER
  OpLeg(const PolyMesh &seed, const EdgeSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  /**
   * @brief Constructs a recipe-step Conway sweep leg: one primitive op swept
   * t_start -> t_end on a fixed seed, no graph edge.
   * @param source Sweep mesh, cloned unless wrapped by SweepSeed::borrow().
   * @param spec Swept operator, parameter endpoints and frame count.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @param blend_fn Crossfade-weight curve.
   * @param easing_fn Easing applied to the sweep parameter.
   */
  inline HS_FLASH_MEMBER
  OpLeg(SweepSeed source, const ParamSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  template <typename... Args>
  OpLeg(PolyMesh &&, const ParamSweepSpec &, Args &&...) = delete;

  /**
   * @brief Constructs a hankin-sweep leg: bakes the angle-independent hankin
   * topology once, computes the arrival classification, and builds the palette
   * mappings.
   * @param seed Base mesh the hankin pattern sweeps on; never cloned into the
   * leg arena, so it must outlive the leg and must not move.
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
  inline HS_FLASH_MEMBER
  OpLeg(const PolyMesh &seed, const HankinSweepSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  template <typename... Args>
  OpLeg(PolyMesh &&, const HankinSweepSpec &, Args &&...) = delete;
  template <typename... Args>
  OpLeg(const PolyMesh &&, const HankinSweepSpec &, Args &&...) = delete;

  /**
   * @brief Constructs a relax leg: clones the seed, relaxes it once, and
   * slerps every vertex from its seed position to its relaxed one.
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
  inline HS_FLASH_MEMBER
  OpLeg(const PolyMesh &seed, const RelaxSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  /**
   * @brief Constructs the Conway-dual bridge's medial leg: builds the shared
   * rectified connectivity of P and slerps every vertex from ambo(P) to
   * ambo(dual(P)) at a fixed emission order.
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
   * @note The connectivity is ambo(P), so every medial face keeps its identity
   * across the slerp. A degree-2 seed vertex makes the s=1 positions
   * many-to-one but leaves the faces well-formed.
   */
  inline HS_FLASH_MEMBER
  OpLeg(const PolyMesh &seed, const MedialSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  /**
   * @brief Constructs the reconcile leg closing a smooth kis/needle path: slerps
   * every vertex of the topology-exact identity mesh (dt/dtd) onto its
   * counterpart in the authored kis/needle mesh, along the caller's
   * nearest-vertex bijection.
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
   * @note Connectivity is held fixed and only the positions move (the
   * MEDIAL_SLERP step path); the caller owns the nearest-vertex correspondence.
   */
  inline HS_FLASH_MEMBER
  OpLeg(const PolyMesh &from_mesh, const ReconcileSpec &spec, Arena &arena,
        MorphDrawFn draw, const PaletteHandoff &handoff,
        const BookendClasses &bookend = BookendClasses{.topology = nullptr,
                                                       .faces = 0},
        BlendWeightFn blend_fn = classic_blend,
        EasingFn easing_fn = math::ease_in_out_sin);

  /**
   * @brief Constructs a gated-swap leg: a partition op with no sweep, drawn as
   * the seed then op(seed) at constant gain, the swap masked by holding the
   * inherited source colour until the late fade.
   * @param seed Mesh the partition op runs on (cloned, not borrowed).
   * @param spec Partition operator and gate length.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh (target keying);
   * defaults to the swept-classification fallback.
   * @note The compiled face count is constant on each side of the swap and
   * changes exactly once, at it.
   */
  inline HS_FLASH_MEMBER OpLeg(const PolyMesh &seed, const GatedSwapSpec &spec,
                               Arena &arena, MorphDrawFn draw,
                               const PaletteHandoff &handoff,
                               const BookendClasses &bookend = BookendClasses{
                                   .topology = nullptr, .faces = 0});

  /**
   * @brief Draws the leg's current frame: an op sweep plus settle, hankin,
   * relax or medial slerp, or the gated swap's seed or partitioned mesh; then
   * compile, palette pre-blend, draw.
   * @param canvas The canvas passed through to the draw callback.
   */
  inline HS_FLASH_MEMBER void step(Canvas &canvas) override;

  /**
   * @brief Redraws the current morph frame without advancing the leg.
   * @note This rebuilds and recompiles the swept mesh on every paused frame;
   *       pausing freezes progression, not morph rendering cost.
   * @param canvas The canvas passed through to the draw callback.
   */
  inline HS_FLASH_MEMBER void step_paused(Canvas &canvas) override;

private:
  /**
   * @brief Initializes the members every leg constructor shares.
   * @param duration Leg length in frames.
   * @param easing_fn Easing applied to the sweep parameter.
   * @param draw Draw callback invoked once per frame.
   */
  HS_FLASH_INLINE OpLeg(int duration, EasingFn easing_fn, MorphDrawFn draw)
      : AnimationBase(duration, false), easing_fn(easing_fn), draw_fn(draw) {}

  inline HS_FLASH_MEMBER void draw_frame(Canvas &canvas);

  /// Bytes an arena keeps free beyond a sweep's own share for alignment.
  static constexpr size_t PREVIOUS_SLACK = 256;

  /**
   * @brief Whether the scratch arenas hold the frame drawn last beside this
   *        frame's sweep.
   * @param vertices Swept vertex count.
   * @param sweep_a scratch_arena_a bytes this frame's sweep holds.
   * @param sweep_b scratch_arena_b bytes this frame's sweep holds.
   * @return True when scratch_arena_a has room for the previous vertices, a
   * second sweep and the compiled mesh, which a sweep's share bounds, and
   * scratch_arena_b for a second sweep's scratch; a leg whose frames do not
   * fit draws without its previous vertices.
   */
  inline HS_FLASH_MEMBER static bool
  previous_fits(size_t vertices, size_t sweep_a, size_t sweep_b);

  /**
   * @brief One vertex of a per-vertex leg's swept mesh at a fractional leg
   * frame, as swept_at() places it.
   * @param i Vertex index.
   * @param frame Leg frame in [0, duration].
   * @return The vertex; a hankin, relax or medial leg's only.
   */
  HS_FLASH_MEMBER math::Vector swept_vertex(size_t i, float frame) const {
    return swept_vertex_at(i, sweep_parameter(frame));
  }

  /**
   * @brief Sweep parameter of a per-vertex leg's swept mesh at a fractional
   * leg frame.
   * @param frame Leg frame in [0, duration].
   * @return The parameter swept_vertex_at() reads.
   */
  inline HS_FLASH_MEMBER float sweep_parameter(float frame) const;

  /**
   * @brief One vertex of a per-vertex leg's swept mesh at a sweep parameter.
   * @param i Vertex index.
   * @param k Sweep parameter, as sweep_parameter() gives it.
   * @return The vertex; a hankin, relax or medial leg's only.
   */
  inline HS_FLASH_MEMBER math::Vector swept_vertex_at(size_t i, float k) const;

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
   * @details Star points are snorm16-decoded and renormalized, so the result
   * is bitwise the leg's own final swept frame, not MeshOps::hankin at the
   * arrival angle. @p out is reset first, dropping any prior binding and class
   * ids.
   */
  inline HS_FLASH_MEMBER static void arrival_mesh(const Landing &landing,
                                                  PolyMesh &out, Arena &arena);

  /**
   * @brief Mid-leg crossfade weight, the swept ctors' default:
   * exactly 0 through the first 20% of the leg, smoothstep to exactly 1 by
   * 80%.
   * @param frame Leg frame in [0, duration]; 0 is the paused initial state.
   * @param duration Whole leg length in frames (sweep plus settle).
   * @return Crossfade weight in [0, 1].
   */
  HS_FLASH_MEMBER static float classic_blend(int frame, int duration) {
    return blend_weight(static_cast<float>(frame) /
                        static_cast<float>(duration));
  }

  /**
   * @brief Trailing-window crossfade weight: exactly 0 until the final @p
   * window frames of the leg, then smoothstep to exactly 1 at the arrival
   * frame.
   * @param frame Leg frame in [0, duration]; 0 is the paused initial state.
   * @param duration Whole leg length in frames (sweep plus settle).
   * @param window Trailing window in frames, clamped to duration - 1 so a leg
   * shorter than the window still opens on the zero plateau.
   * @return Crossfade weight in [0, 1].
   */
  inline static float trailing_blend(int frame, int duration, int window);

  /**
   * @brief Trailing blend over the default TRAILING_BLEND_FRAMES window;
   * the BlendWeightFn form.
   * @param frame Leg frame in [0, duration]; 0 is the paused initial state.
   * @param duration Whole leg length in frames (sweep plus settle).
   * @return Crossfade weight in [0, 1].
   * @details Reaches exactly 1 on the final frame.
   */
  static float trailing_blend(int frame, int duration) {
    return trailing_blend(frame, duration, TRAILING_BLEND_FRAMES);
  }

private:
  /**
   * @brief Arena-allocated leg state — keeps OpLeg inline size small.
   */
  struct Transients {
    PolyMesh seed; /**< Cloned geometry for owned CONWAY_SWEEP, RELAX_SLERP
                      and GATED_SWAP legs; connectivity only for MEDIAL_SLERP;
                      empty for HANKIN_SWEEP and borrowed CONWAY_SWEEP. */
    const PolyMesh *seed_ref =
        nullptr; /**< CONWAY_SWEEP swept-op source: &seed for a cloned seed, the
                    caller's live mesh for a borrowed one. */
    size_t seed_faces = 0; /**< Seed face count; set by every kind, including
                              the ones that keep no seed. */
    LegKind kind = LegKind::CONWAY_SWEEP; /**< Swept-mesh production path. */
    CompiledHankin hankin; /**< Baked topology (HANKIN_SWEEP legs). */
    ArenaVector<math::Snorm3>
        hk_final; /**< Arrival star points, snorm16-packed (HANKIN_SWEEP). */
    ConwayGraph::MorphOp op =
        ConwayGraph::MorphOp::TRUNCATE; /**< Swept operator (CONWAY_SWEEP). */
    SwapOp swap_op = SwapOp::KIS;       /**< Partition op (GATED_SWAP). */
    bool reverse = false;               /**< Traversing to_node -> from_node. */
    /** Operator-sweep frames N; GATED_SWAP gate length on each side of the swap. */
    int sweep_frames = 1;
    int settle_frames = 0;                /**< Relax-slerp frames. */
    float t_start = 0, t_end = 0;         /**< Clamped sweep endpoints. */
    float twist_start = 0, twist_end = 0; /**< Snub twist endpoints. */
    ArenaVector<math::Vector>
        relaxed; /**< Relaxed endpoint vertices (settling and relax legs). */
    ArenaVector<math::Snorm3>
        medial_a; /**< Opening endpoint, snorm16-packed (MEDIAL_SLERP:
                     ambo(P) on the medial bridge, the identity mesh's vertices
                     on a reconcile leg). */
    ArenaVector<math::Snorm3>
        medial_b; /**< Arrival endpoint, snorm16-packed (MEDIAL_SLERP:
                     ambo(dual(P)) on the medial bridge, the authored positions
                     on a reconcile leg). */
    ArenaVector<uint16_t>
        topo; /**< Arrival classes, or placeholders for a partial bookend. */
    bool topo_is_bookend = false; /**< topo holds the bookend classification
                                     verbatim, so target_topo is not built. */
    ArenaVector<uint16_t>
        target_topo; /**< Per-swept-face target class (bookend grouping); empty
                        when topo already carries it. */
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
    Landing landing;           /**< Arrival data exposed to the effect. */
    float drawn_weight = 0.0f; /**< Crossfade weight of the last draw. */
    int drawn_frame = 0;       /**< Leg frame of the last draw. */
    size_t drawn_vertices = 0; /**< Swept vertex count of the last draw. */
    bool drawn = false;        /**< Whether the leg has drawn a frame. */
  };

  /** @brief Split of a swept leg's eased progress into its two stages. */
  using SweepPhase = detail::SweepPhase;

  /**
   * @brief Sweep and settle fractions at one leg frame.
   * @details One eased clock spans sweep and settle; reverse legs un-settle
   * first.
   */
  inline __attribute__((always_inline)) SweepPhase
  sweep_phase(const Transients &tr, float progress) const;

  /**
   * @brief Builds a Conway leg's swept mesh at one parameter, slerped toward
   * the relaxed endpoint by @p settle_alpha.
   */
  inline __attribute__((always_inline)) static void
  conway_at(const Transients &tr, PolyMesh &out, Arena &target, Arena &temp,
            float t, float twist, float settle_alpha);

  /**
   * @brief Builds a swept leg's mesh at one leg frame.
   * @param tr Leg transients.
   * @param frame Leg frame in [0, duration].
   * @param out Output mesh, allocated from scratch_arena_a.
   */
  inline HS_FLASH_MEMBER void swept_at(const Transients &tr, int frame,
                                       PolyMesh &out) const;

  /** @brief Slerps a swept mesh's vertices toward the relaxed endpoint. */
  inline __attribute__((always_inline)) static void
  settle_toward_relaxed(const Transients &tr, PolyMesh &swept,
                        float settle_alpha);

  /**
   * @brief Clamps a recipe-step sweep endpoint into its op's
   * topology-constant interval.
   * @param t Authored endpoint.
   * @param truncate Whether the swept op is truncate.
   * @param far_side Whether a truncate leg sweeps past ambo.
   * @param trunc_floor Truncate's birth floor for this leg.
   */
  inline __attribute__((always_inline)) static float
  clamp_sweep_param(float t, bool truncate, bool far_side, float trunc_floor);

  inline HS_FLASH_MEMBER Transients &
  init_transients(LegKind kind, int sweep_frames, Arena &arena,
                  const PaletteHandoff &handoff, BlendWeightFn blend_fn);

  inline HS_FLASH_MEMBER static void
  bind_sweep_seed(Transients &tr, SweepSeed source, Arena &arena);

  /** @brief Visits dual corners in first-seen face-walk order.
   * @tparam Emit Corner consumer.
   * @param dual Dual seed mesh.
   * @param arena Scratch storage for visited vertices.
   * @param emit Receives each vertex index and its first host face. */
  template <typename Emit>
  inline HS_FLASH_MEMBER static void
  first_seen_corners(const PolyMesh &dual, Arena &arena, Emit &&emit);

  inline HS_FLASH_MEMBER static const uint8_t *
  dual_closing_palettes(const PolyMesh &dual, const PaletteHandoff &handoff,
                        Arena &scratch);

  /** @brief Per-leg switches for init_conway. */
  using ConwayLegFlags = detail::ConwayLegFlags;

  /**
   * @brief Shared CONWAY_SWEEP construction tail: arrival classification
   * (relaxed form on a forward settling leg), start centroids, and palette
   * mappings.
   * @param handoff Palette provenance of the departed mesh.
   * @param bookend Bookend grouping of the arrival mesh.
   * @param arena Leg arena for the hoisted state.
   * @param flags Settle, jitterbug and bridge-provenance switches.
   * @details Requires tr.op, tr.reverse, the tr.t_start/t_end and
   * tr.twist_start/twist_end endpoints and tr.settle_frames already set by the
   * calling constructor.
   */
  inline HS_FLASH_MEMBER void init_conway(const PaletteHandoff &handoff,
                                          const BookendClasses &bookend,
                                          Arena &arena,
                                          const ConwayLegFlags &flags);

  /**
   * @brief Start-parameter face centroids of a closing dual-bridge leg.
   * @return Centroids in the leg's own start emission order, one per handoff
   * face, allocated from scratch_arena_a.
   */
  inline __attribute__((always_inline)) static const math::Vector *
  closing_bridge_centroids(const Transients &tr, const PolyMesh &seed,
                           const PaletteHandoff &handoff);

  /**
   * @brief Relaxes a settling leg's endpoint into tr.relaxed.
   * @param relaxed_mesh Receives the relaxed mesh (scratch-backed).
   * @return The mesh the arrival classification is taken from: the relaxed
   * form on a forward leg, which lands on it; @p arrival on a reverse leg,
   * which un-settles at the start parameter.
   */
  inline __attribute__((always_inline)) static PolyMesh *
  hoist_settle_endpoint(Transients &tr, const PolyMesh &seed, PolyMesh &arrival,
                        Arena &arena, PolyMesh &relaxed_mesh);

  /**
   * @brief Face centroids of the mesh a Conway leg's first frame draws, the
   * geometric provenance source.
   * @param relaxed_mesh Relaxed start, drawn first on a reverse settling leg
   * (settle_alpha == 1 there); otherwise the plain op at t_start is.
   * @param start_mesh Receives the t_start mesh when it is the one drawn.
   * @return Centroids allocated from scratch_arena_a.
   */
  inline __attribute__((always_inline)) static const math::Vector *
  opening_frame_centroids(const Transients &tr, const PolyMesh &seed,
                          const PolyMesh &relaxed_mesh, bool settle,
                          PolyMesh &start_mesh);

  /**
   * @brief Emission-order prefix of a Conway leg corresponding 1:1 to a node
   * base mesh at the boundary swaps: the seed's face count, plus its
   * vertex-orbit faces on the jitterbug bridge.
   */
  __attribute__((always_inline)) static size_t
  survivor_faces(const Transients &tr, const PolyMesh &seed, bool jitterbug) {
    return jitterbug ? tr.seed_faces + seed.vertices.size() : tr.seed_faces;
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
  inline HS_FLASH_MEMBER void init_gated(const PaletteHandoff &handoff,
                                         const BookendClasses &bookend,
                                         Arena &arena);

  /**
   * @brief Index of the seed-side identity ramp for a departed palette,
   * appending one on first use.
   */
  inline __attribute__((always_inline)) static uint8_t
  intern_seed_ramp(Transients &tr, uint8_t pal);

  /**
   * @brief Steps a gated-swap leg: the seed or the partitioned mesh, then the
   * shared frame tail.
   * @param canvas The canvas passed through to the draw callback.
   * @param frame Clamped frame index; the seed side draws frames [0, gate],
   * with 0 the paused initial state; the swap and opening side draws
   * [gate + 1, 2 * gate + 1].
   */
  inline HS_FLASH_MEMBER void step_gated(Canvas &canvas, int frame);

  /**
   * @brief Applies the leg's partition operator to its seed.
   * @param tr Leg transients holding the seed and the operator.
   * @param target Arena receiving the output mesh.
   * @param temp Arena holding the operator's transient scratch.
   * @return The partitioned mesh.
   */
  inline HS_FLASH_MEMBER static PolyMesh swap_at(const Transients &tr,
                                                 Arena &target, Arena &temp);

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
  inline HS_FLASH_MEMBER static void
  kis_provenance(const Transients &tr, const PolyMesh &arrival,
                 const PaletteHandoff &handoff, uint8_t *from);

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
   * vertex's face orbit.
   */
  inline HS_FLASH_MEMBER static void
  dual_provenance(const Transients &tr, const PolyMesh &arrival,
                  const PaletteHandoff &handoff, uint8_t *from);

  /**
   * @brief Index of the vertex of the face at @p off (with @p n sides)
   * nearest @p v; the first wins a tie.
   */
  inline __attribute__((always_inline)) static uint16_t
  nearest_orbit_vertex(const PolyMesh &mesh, size_t off, int n,
                       const math::Vector &v);

  /**
   * @brief Runs the edge's operator on the seed at one parameter value.
   */
  inline HS_FLASH_MEMBER static PolyMesh run_op(ConwayGraph::MorphOp op,
                                                const PolyMesh &seed,
                                                Arena &target, Arena &temp,
                                                float t, float twist);

  /**
   * @brief Hoists arrival classes or reserves a partial bookend's face slots.
   * @return Whether the bookend classification was reused verbatim.
   */
  inline HS_FLASH_MEMBER static bool
  hoist_arrival_topology(PolyMesh &arrival, const BookendClasses &bookend,
                         Arena &arena, ArenaVector<uint16_t> &topology);

  /**
   * @brief Compile the swept mesh, check its face count against the hoisted
   * classification, pre-blend ramps and pass the mesh and Shading to the callback.
   * @param canvas The canvas passed through to the draw callback.
   * @param swept This frame's swept mesh (scratch-backed).
   * @param w Crossfade weight in [0, 1] the caller already resolved.
   * @param seed_side Draw the gate's seed-side tables (GATED_SWAP only): the
   * seed classification and its identity ramp table, which w == 0 leaves at the
   * departed palettes.
   * @param drawn @p swept's vertices at the leg frame drawn last, index for
   * index, or null when the draw does not move; they reach the draw as
   * `Shading::previous_vertices` when compiling keeps every vertex in place.
   * @param previous_frame The leg frame drawn last when the leg moves vertex by
   * vertex from it, or -1; it reaches the draw through `Shading::leg` when
   * compiling keeps every vertex in place.
   * @param current_frame The leg frame being drawn, or -1 for a gated draw.
   */
  inline HS_FLASH_MEMBER void
  finish_frame(Canvas &canvas, const PolyMesh &swept, float w,
               bool seed_side = false, math::Vector *drawn = nullptr,
               int previous_frame = -1, int current_frame = -1);

  /**
   * @brief Bakes this frame's ramp per interned pair into scratch_arena_b.
   * @param tr Leg transients holding the ramp tables.
   * @param w Crossfade weight in [0, 1].
   * @param seed_side Bake the seed-side identity ramps instead.
   * @return One LUT per ramp index, scratch-backed.
   */
  inline __attribute__((always_inline)) static BakedPalette *
  blend_ramps(const Transients &tr, float w, bool seed_side);

  /**
   * @brief One dynamic Hankin vertex at an arrival or sweep fraction.
   * @tparam ARRIVAL Whether to return the normalized arrival endpoint.
   * @param tr Baked leg endpoints.
   * @param i Dynamic vertex index.
   * @param k Sweep fraction.
   * @return Unit-sphere vertex position.
   */
  template <bool ARRIVAL>
  inline __attribute__((always_inline)) static math::Vector
  hankin_vertex_at(const Transients &tr, size_t i, float k);

  /**
   * @brief One relax vertex at an arrival or sweep fraction.
   * @tparam ARRIVAL Whether to return the arrival endpoint verbatim.
   * @param tr Seed and relaxed endpoints.
   * @param i Vertex index.
   * @param k Sweep fraction.
   * @return Unit-sphere vertex position.
   */
  template <bool ARRIVAL>
  inline __attribute__((always_inline)) static math::Vector
  relax_vertex_at(const Transients &tr, size_t i, float k);

  /**
   * @brief One medial vertex at an arrival or sweep fraction.
   * @tparam ARRIVAL Whether to return the normalized arrival endpoint.
   * @param tr Packed leg endpoints.
   * @param i Vertex index.
   * @param k Sweep fraction.
   * @return Unit-sphere vertex position.
   */
  template <bool ARRIVAL>
  inline __attribute__((always_inline)) static math::Vector
  medial_vertex_at(const Transients &tr, size_t i, float k);

  /**
   * @brief Copies a leg's constant topology into the swept mesh.
   * @param out Output mesh whose face arrays are bound and filled.
   * @param arena Arena backing the face arrays.
   * @param face_counts Per-face side counts.
   * @param faces Flat per-face vertex index list.
   */
  inline HS_FLASH_MEMBER static void
  copy_topology(PolyMesh &out, Arena &arena,
                const ArenaVector<uint8_t> &face_counts,
                const ArenaVector<uint16_t> &faces);

  /**
   * @brief Deep-copies a seed's geometry, leaving its class ids behind.
   * @param src Seed mesh.
   * @param dst Destination mesh, populated in place.
   * @param arena Arena backing the destination arrays.
   */
  inline HS_FLASH_MEMBER static void
  clone_geometry(const PolyMesh &src, PolyMesh &dst, Arena &arena);

  /**
   * @brief Builds a hankin leg's swept mesh at one slerp fraction.
   * @param tr Leg transients holding the baked topology and both endpoints.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1].
   */
  inline HS_FLASH_MEMBER static void
  hankin_at(const Transients &tr, PolyMesh &out, Arena &arena, float k);

  /**
   * @brief Builds a relax leg's swept mesh at one slerp fraction.
   * @param tr Leg transients holding the seed and its relaxed endpoint.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1].
   */
  inline HS_FLASH_MEMBER static void
  relax_at(const Transients &tr, PolyMesh &out, Arena &arena, float k);

  /**
   * @brief Builds a medial leg's swept mesh at one slerp fraction, decoding the
   * snorm16-packed endpoints per vertex.
   * @param tr Leg transients holding the packed a_e/b_e sets and connectivity.
   * @param out Output mesh, allocated from @p arena.
   * @param arena Arena backing the output vectors.
   * @param k Slerp fraction in [0, 1]; slerp() renormalizes the near-unit
   * decoded inputs, so the quantization is absorbed into the direction blend.
   */
  inline HS_FLASH_MEMBER static void
  medial_at(const Transients &tr, PolyMesh &out, Arena &arena, float k);

  /**
   * @brief Unit-sphere vertex-average centroid of every face.
   * @param m Mesh whose faces are reduced.
   * @param arena Arena receiving the centroid array.
   * @return One unit vector per face, in emission order.
   */
  inline HS_FLASH_MEMBER static const math::Vector *
  face_centroids(const PolyMesh &m, Arena &arena);

  /**
   * @brief Index of the departed face nearest to a start-face centroid.
   * @param c Unit centroid of the swept face at the start parameter.
   * @param handoff Departed-node provenance (centroids non-null).
   * @return Index into the handoff arrays.
   */
  inline HS_FLASH_MEMBER static size_t
  nearest_prev_face(const math::Vector &c, const PaletteHandoff &handoff);

  /**
   * @brief Index of the point in @p points[0, n) nearest @p c by squared
   * distance; the first wins a tie.
   */
  inline __attribute__((always_inline)) static size_t
  nearest_point(const math::Vector &c, const math::Vector *points, size_t n);

  /**
   * @brief Crossfade weight at fraction p of the leg (classic_blend's curve).
   */
  inline static float blend_weight(float p);

  // always_inline: an out-of-line copy inherits no cold attribute from its
  // HS_FLASH_MEMBER caller and would land in ITCM.
  inline __attribute__((always_inline)) static const uint16_t *
  resolve_target_topology(Transients &tr, const PolyMesh &arrival,
                          const PaletteHandoff &handoff,
                          const BookendClasses &bookend, Arena &arena,
                          size_t survivors);

  /**
   * @brief Landed host face of each corner face on a closing leg.
   * @return One host face index per corner face, allocated from
   * scratch_arena_a.
   * @details Corner face k is born on the k-th first-seen vertex of the dual
   * seed's face walk; the face whose walk first reaches that vertex is the
   * corner's host.
   */
  inline __attribute__((always_inline)) static int *
  dual_corner_hosts(const Transients &tr, size_t landed, size_t total);

  /** @brief Per-face from-palette selection of one palette mapping. */
  struct DepartedPalettes {
    const PaletteHandoff &handoff; /**< Departed-node provenance. */
    const uint8_t *forced_from;    /**< Per-face override, or null to derive. */
    const math::Vector *start_centroid; /**< Start-parameter centroids, or
                                           null for the emission-order
                                           mapping. */
    const uint16_t *target_topo;        /**< Target class per swept face. */
    size_t total;                       /**< Swept face count. */
    bool full_correspondence;           /**< Every swept face maps to a
                                           departed face by centroid. */
    bool *prev_used;   /**< Claimed departed faces (full correspondence
                        only). */
    int *newborn_from; /**< Newborn palette slot -> from-palette of its first
                          face; -1 until set. */

    /**
     * @brief From-palette of swept face @p f.
     * @param to The face's landed palette, the fallback for a newborn face,
     * which skips the crossfade.
     */
    inline __attribute__((always_inline)) uint8_t for_face(size_t f,
                                                           uint8_t to) const;

    /** @brief Palette of the departed face nearest @p f's start centroid,
     * claimed so no other face maps to it. */
    inline __attribute__((always_inline)) uint8_t claim_nearest(size_t f) const;

    /** @brief Newborn face @p f's from-palette, shared across its palette
     * slot and taken from the slot's first face's nearest departed face. */
    inline __attribute__((always_inline)) uint8_t newborn(size_t f) const;
  };

  // always_inline: an out-of-line copy inherits no cold attribute from its
  // HS_FLASH_MEMBER caller and would land in ITCM.
  inline __attribute__((always_inline)) static uint8_t
  intern_palette_ramp(Transients &tr, uint8_t from, uint8_t to);

  /**
   * @brief Builds the per-face from-palettes, the
   * shuffled target assignment, and the distinct ramp-pair table.
   * @param tr Leg transients being populated.
   * @param arrival Classified arrival mesh; its face centroids pick each
   * collapsed face's host class.
   * @param handoff Departed-node provenance.
   * @param bookend Arrival-node bookend grouping the targets key on.
   * @param arena Leg arena for the face -> ramp table.
   * @param start_centroid Unit centroid per swept face at the start parameter,
   * or nullptr for the emission-order mapping.
   * @param survivors Emission-order face-prefix length corresponding 1:1 to a
   * node base mesh at the boundary swaps.
   * @param forced_from Per-arrival-face from-palette, or nullptr to derive one.
   * @details With centroids, a full-correspondence departure maps every face to
   * its nearest departed centroid as a checked bijection; on a node-prefix
   * departure each newborn palette slot inherits its first face's nearest
   * departed palette.
   */
  inline HS_FLASH_MEMBER void
  build_palette_mapping(Transients &tr, const PolyMesh &arrival,
                        const PaletteHandoff &handoff,
                        const BookendClasses &bookend, Arena &arena,
                        const math::Vector *start_centroid, size_t survivors,
                        const uint8_t *forced_from = nullptr);

  /**
   * @brief Whether the per-frame scratch peak fits scratch_arena_b: the ramp
   * array plus one baked LUT per non-identity pair.
   */
  inline __attribute__((always_inline)) static bool
  ramps_fit_scratch(const Transients &tr);

  /**
   * @brief Allocates the leg's Transients.
   * @param arena Leg arena backing the leg's state.
   * @return The fresh Transients.
   */
  inline Transients &bind_transients(Arena &arena);

  /**
   * @brief Stamps the watermark check_alive() tests past the leg's last
   * allocation.
   */
  inline void seal_transients();

  /**
   * @brief Traps if the leg arena was reclaimed while the leg is still live.
   * @details The debug block stamp also detects a rewind followed by refill.
   */
  inline void check_alive() const;

  Transients *buf = nullptr;  /**< Pointer to arena-allocated leg state. */
  Arena *leg_arena = nullptr; /**< Arena backing buf and the Transients
                                 vectors. */
  size_t live_end = 0; /**< Arena offset past all sealed leg allocations. */
#ifndef NDEBUG
  ArenaBlockStamp stamp;
  size_t live_bytes = 0;
#endif
  EasingFn easing_fn;  /**< Easing applied to the sweep parameter. */
  MorphDrawFn draw_fn; /**< Per-frame draw callback. */
};

} // namespace Animation

#include "core/animation/opleg/build.h"
#include "core/animation/opleg/frame.h"
#include "core/animation/opleg/sweep.h"
#include "core/animation/opleg/swaps.h"
#include "core/animation/opleg/ramps.h"
