/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

#include "color/composition.h"
#include "mesh/conway.h"
#include "mesh/conway_graph.h"

namespace Animation {

/**
 * @brief Animates a vertex-interpolated crossfade between two meshes.
 * @details Owns transient state (cloned meshes + SLERP buffers) via a pointer to
 * arena-allocated storage to keep inline size small for TimelineEvent. No
 * destructor reclaims the arena, so the transient bytes persist until the caller
 * compacts or resets it.
 *
 * Crossfade contract: only the incoming mesh (mesh_B, dest topology) morphs —
 * its vertices SLERP from their nearest source vertex toward their dest each
 * frame. The outgoing mesh_A (source clone) holds geometry and fades via opacity
 * (op_A = 1 - alpha); the opacities sum to 1 for constant brightness. source is
 * cloned, not borrowed, so the animation survives the caller recycling it.
 * draw_outgoing/draw_incoming shade the two halves independently.
 */
class MeshMorph : public AnimationBase<MeshMorph> {
public:
  /**
   * @brief Non-owning per-half draw callback: `void(Canvas&, const MeshState&,
   * float opacity)`. A StoredFunctionRef (held as a member, invoked across many
   * frames) rejects rvalue temporaries, so a dangling inline lambda is a compile
   * error rather than a silent use-after-free.
   */
  using MorphDrawFn = StoredFunctionRef<void(Canvas &, const MeshState &, float)>;

  /**
   * @brief Constructs a MeshMorph with separate shading for the two halves.
   * @param source The outgoing mesh (cloned, not borrowed).
   * @param dest The incoming mesh whose topology the morph targets.
   * @param arena Arena providing backing storage for cloned meshes and buffers.
   * @param draw_outgoing Draw callback for the fading-out source clone.
   * @param draw_incoming Draw callback for the fading-in morphing mesh.
   * @param duration The crossfade duration in frames.
   * @param easing_fn The easing function applied to crossfade progress.
   * @note The cloned meshes and position buffers are arena-allocated with no
   *   per-instance reclamation; the caller must compact the arena between
   *   successive morphs or it grows unbounded (see HankinSolids/MeshFeedback).
   */
  MeshMorph(const MeshState &source, const MeshState &dest, Arena &arena,
            MorphDrawFn draw_outgoing, MorphDrawFn draw_incoming, int duration,
            EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(duration, false), easing_fn(easing_fn),
        draw_outgoing(draw_outgoing), draw_incoming(draw_incoming) {
    HS_CHECK(duration >= 1, "MeshMorph duration must be a positive frame count");
    HS_CHECK(!source.vertices.is_empty());
    HS_CHECK(!dest.vertices.is_empty());
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();

    MeshOps::clone(source, buf_->mesh_A, arena);
    MeshOps::clone(dest, buf_->mesh_B, arena);

    // step() indexes mesh_B.vertices by dest's vertex count, so trap here if
    // clone() ever welds/dedups rather than writing OOB per-frame.
    HS_CHECK(buf_->mesh_B.vertices.size() == dest.vertices.size());

    buf_->start_pos.bind(arena, dest.vertices.size());
    buf_->end_pos.bind(arena, dest.vertices.size());

    // Symmetry-breaking twist to avoid degenerate nearest-vertex mapping
    Vector twist_axis = Vector(0.0f, 0.0f, 1.0f);
    auto has_pole = [](const MeshState &m) {
      for (const auto &v : m.vertices)
        if (std::abs(v.z) > 0.99f && std::abs(v.x) < 0.01f &&
            std::abs(v.y) < 0.01f)
          return true;
      return false;
    };
    bool has_poles = has_pole(source) || has_pole(dest);
    if (has_poles) {
      twist_axis = Vector(1.0f, 1.0f, 1.0f).normalized();
    }
    Quaternion twist = make_rotation(twist_axis, 0.05f);

    // Nearest-vertex matching (greatest dot) and per-frame slerp both require
    // unit-length inputs; all mesh sources sit on the unit sphere.
    for (const auto &v : source.vertices)
      HS_CHECK(std::abs(dot(v, v) - 1.0f) < 1e-3f, "MeshMorph source vertex not unit-length");
    for (const auto &v : dest.vertices)
      HS_CHECK(std::abs(dot(v, v) - 1.0f) < 1e-3f, "MeshMorph dest vertex not unit-length");

    // Build nearest-vertex correspondence: an O(V_dest * V_source) brute force,
    // run once at construction. Matched by greatest dot product against the
    // twist-biased dest vertex (the twist breaks ties on symmetric meshes).
    for (size_t i = 0; i < dest.vertices.size(); ++i) {
      Vector v_biased = rotate(dest.vertices[i], twist);
      int best_idx = 0;
      float max_dot = -9999.0f;
      for (size_t j = 0; j < source.vertices.size(); ++j) {
        float d = dot(v_biased, source.vertices[j]);
        if (d > max_dot) {
          max_dot = d;
          best_idx = static_cast<int>(j);
        }
      }
      buf_->start_pos.push_back(source.vertices[best_idx]);
      buf_->end_pos.push_back(dest.vertices[i]);
    }
  }

  // Borrow contract: the draw callbacks are non-owning StoredFunctionRefs read
  // every frame, so they must outlive the timeline; StoredFunctionRef rejects a
  // temporary at the MorphDrawFn parameter, so no `= delete` overload is needed.

  /**
   * @brief Steps the crossfade: interpolates vertices and renders both halves.
   * @param canvas The canvas buffer passed to the draw callbacks.
   */
  void step(Canvas &canvas) override {
    // Increment-first so the final frame lands exactly on the destination mesh
    // (alpha == 1); the skipped progress==0 frame is immaterial.
    AnimationBase::step(canvas);

    float progress = hs::clamp(static_cast<float>(t) / duration, 0.0f, 1.0f);
    float alpha = easing_fn(progress);

    for (size_t i = 0; i < buf_->end_pos.size(); ++i) {
      buf_->mesh_B.vertices[i] =
          slerp(buf_->start_pos[i], buf_->end_pos[i], alpha);
    }

    float op_A = 1.0f - alpha;
    if (op_A > 0.01f)
      draw_outgoing(canvas, buf_->mesh_A, op_A);
    if (alpha > 0.01f)
      draw_incoming(canvas, buf_->mesh_B, alpha);
  }

private:
  /**
   * @brief Arena-allocated transient data — keeps MeshMorph inline size small.
   */
  struct Transients {
    MeshState mesh_A;            /**< Outgoing mesh clone. */
    MeshState mesh_B;           /**< Incoming morphing mesh clone. */
    ArenaVector<Vector> start_pos; /**< Per-vertex nearest-source start points. */
    ArenaVector<Vector> end_pos;   /**< Per-vertex dest end points. */
  };

  Transients *buf_;          /**< Pointer to arena-allocated transient state. */
  EasingFn easing_fn;        /**< Easing curve applied to crossfade progress. */
  MorphDrawFn draw_outgoing; /**< Draw callback for the outgoing half. */
  MorphDrawFn draw_incoming; /**< Draw callback for the incoming half. */
};

/**
 * @brief Animates one Conway-operator parameter sweep along a graph edge
 * (docs/conway_morph_spec.md, section 4.1).
 * @details Per frame: run the edge's single op at t(frame) in scratch,
 * settle-slerp toward the relaxed endpoint inside the settle window, compile,
 * attach the leg's hoisted classification, pre-blend the (from, to) palette
 * ramps at w(frame), and hand the mesh to the draw callback. Exactly one mesh
 * is drawn per frame. Bulk state lives in an arena-allocated Transients (same
 * survival contract as MeshMorph); the caller compacts the arena between legs.
 */
class ConwayMorph : public AnimationBase<ConwayMorph> {
public:
  static constexpr int PALETTES = BakedPaletteBank::N;
  /** Distinct (from, to) ramp pairs a leg may carry; bounds the per-frame
   * blended-LUT scratch (PAIRS x 3 KB in scratch_arena_b). */
  static constexpr int MAX_BLEND_PAIRS = 8;

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
  };

  /**
   * @brief Non-owning draw callback: `void(Canvas&, const MeshState&,
   * const Shading&)`. StoredFunctionRef rejects rvalue temporaries.
   */
  using MorphDrawFn =
      StoredFunctionRef<void(Canvas &, const MeshState &, const Shading &)>;

  /**
   * @brief Palette provenance of the departed node (spec sections 2.5/2.6).
   * @details prev_face_palette/prev_face_sides describe the node base mesh the
   * leg departs from, in emission order; consumed by the constructor only.
   */
  struct PaletteHandoff {
    const BakedPaletteBank *bank = nullptr; /**< The effect's baked source LUTs. */
    const uint8_t *prev_face_palette = nullptr; /**< Per-face palette of the departed base mesh. */
    const uint8_t *prev_face_sides = nullptr;   /**< Per-face clean side counts (class-signature mapping). */
    size_t prev_faces = 0;           /**< Face count of the departed base mesh. */
    bool by_class_signature = false; /**< DUAL_SWAP departure: map by side count, not emission order. */
  };

  /**
   * @brief Leg-static arrival data the effect's completion consumes.
   * @details Arena-backed; valid until the leg arena is compacted.
   */
  struct Landing {
    const int *topology = nullptr; /**< Arrival classification, one id per swept face. */
    size_t faces = 0;              /**< Swept face count. */
    size_t primary_faces = 0;      /**< Seed face count (emission-order prefix). */
    std::array<uint8_t, PALETTES> to_palette{}; /**< Slot -> landed palette index. */
  };

  /**
   * @brief Constructs one leg: clones the seed, computes the arrival
   * classification (relaxed form when settling), and builds the palette
   * mappings.
   * @param seed Seed mesh the op sweeps on (cloned, not borrowed).
   * @param edge Graph edge being traversed.
   * @param reverse True when traversing to_node -> from_node.
   * @param arena Leg arena backing the cloned seed and hoisted state.
   * @param draw Draw callback invoked once per frame.
   * @param handoff Palette provenance of the departed node.
   * @param sweep_frames Operator-sweep frames (N).
   * @param settle_frames Relax-slerp frames (S); 0 unless the edge settles.
   * @param easing_fn Easing applied to the sweep parameter.
   */
  HS_COLD_MEMBER ConwayMorph(const PolyMesh &seed,
                             const ConwayGraph::EdgeSpec &edge, bool reverse,
                             Arena &arena, MorphDrawFn draw,
                             const PaletteHandoff &handoff, int sweep_frames,
                             int settle_frames,
                             EasingFn easing_fn = ease_in_out_sin)
      : AnimationBase(sweep_frames + settle_frames, false),
        easing_fn(easing_fn), draw_fn(draw) {
    HS_CHECK(sweep_frames >= 1, "ConwayMorph needs a positive sweep length");
    HS_CHECK(settle_frames >= 0 && (edge.settle || settle_frames == 0),
             "ConwayMorph: settle frames on a non-settling edge");
    HS_CHECK(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0);
    buf_ = new (arena.allocate(sizeof(Transients), alignof(Transients)))
        Transients();
    Transients &tr = *buf_;

    MeshOps::clone(seed, tr.seed, arena);
    tr.op = edge.op;
    tr.reverse = reverse;
    tr.sweep_frames = sweep_frames;
    tr.settle_frames = settle_frames;
    tr.bank = handoff.bank;

    // Clamp both endpoints inside the topology-constant open interval; the
    // truncate upper clamp dodges the ambo short-circuit at exactly 0.5.
    auto clamp_param = [&](float t) {
      t = std::max(t, ConwayGraph::T_EPS);
      if (edge.op == ConwayGraph::MorphOp::TRUNCATE)
        t = std::min(t, 0.5f - ConwayGraph::T_EPS);
      return t;
    };
    tr.t_start = clamp_param(reverse ? edge.t_to : edge.t_from);
    tr.t_end = clamp_param(reverse ? edge.t_from : edge.t_to);
    tr.twist_start = reverse ? edge.twist_to : edge.twist_from;
    tr.twist_end = reverse ? edge.twist_from : edge.twist_to;

    {
      ScratchScope sa(scratch_arena_a);
      ScratchScope sb(scratch_arena_b);

      PolyMesh arrival = run_op(tr.op, tr.seed, scratch_arena_a,
                                scratch_arena_b, tr.t_end, tr.twist_end);
      // Classification is hoisted per leg, taken at arrival geometry; a
      // settling forward leg lands on the relaxed form, so classify that.
      PolyMesh *classified = &arrival;
      PolyMesh relaxed_mesh;
      if (edge.settle) {
        if (!reverse) {
          relaxed_mesh =
              MeshOps::relax(arrival, scratch_arena_b, scratch_arena_a, 50);
          classified = &relaxed_mesh;
        } else {
          // Reverse legs un-settle at the start parameter.
          PolyMesh start = run_op(tr.op, tr.seed, scratch_arena_a,
                                  scratch_arena_b, tr.t_start, tr.twist_start);
          relaxed_mesh =
              MeshOps::relax(start, scratch_arena_b, scratch_arena_a, 50);
        }
        tr.relaxed.bind(arena, relaxed_mesh.vertices.size());
        tr.relaxed.append_bulk(relaxed_mesh.vertices.data(),
                               relaxed_mesh.vertices.size());
        HS_CHECK(tr.relaxed.size() == arrival.vertices.size(),
                 "ConwayMorph: relax changed the vertex count");
      }
      MeshOps::classify_faces_by_topology(*classified, scratch_arena_a,
                                          scratch_arena_b, arena);
      tr.topo = std::move(classified->topology);
      HS_CHECK(tr.topo.size() == classified->face_counts.size());

      build_palette_mapping(tr, *classified, handoff, arena);
    }
  }

  /**
   * @brief Steps the sweep: op at t(frame), settle slerp, compile, palette
   * pre-blend, draw.
   * @param canvas The canvas passed through to the draw callback.
   * @details HS_COLD: once-per-frame orchestration; the hot loops live in the
   * (already cold) Conway ops and the mesh scan.
   */
  HS_COLD_MEMBER void step(Canvas &canvas) override {
    AnimationBase::step(canvas);
    Transients &tr = *buf_;
    const int frame =
        static_cast<int>(std::min<uint32_t>(t, static_cast<uint32_t>(duration)));

    // Frame -> sweep position + settle blend (forward settles at the end,
    // reverse legs un-settle over the opening window).
    int sweep_frame = frame;
    float settle_alpha = 0.0f;
    if (tr.settle_frames > 0) {
      if (!tr.reverse) {
        sweep_frame = std::min(frame, tr.sweep_frames);
        if (frame > tr.sweep_frames)
          settle_alpha = static_cast<float>(frame - tr.sweep_frames) /
                         static_cast<float>(tr.settle_frames);
      } else {
        sweep_frame = std::max(0, frame - tr.settle_frames);
        if (frame < tr.settle_frames)
          settle_alpha = 1.0f - static_cast<float>(frame) /
                                    static_cast<float>(tr.settle_frames);
      }
    }
    float k = easing_fn(static_cast<float>(sweep_frame) /
                        static_cast<float>(tr.sweep_frames));
    float tp = tr.t_start + (tr.t_end - tr.t_start) * k;
    float tw = tr.twist_start + (tr.twist_end - tr.twist_start) * k;

    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    PolyMesh swept;
    {
      HS_PROFILE(hk_conway_op);
      swept = run_op(tr.op, tr.seed, scratch_arena_a, scratch_arena_b, tp, tw);
      if (settle_alpha > 0.0f) {
        HS_CHECK(swept.vertices.size() == tr.relaxed.size());
        for (size_t i = 0; i < swept.vertices.size(); ++i)
          swept.vertices[i] =
              slerp(swept.vertices[i], tr.relaxed[i], settle_alpha);
      }
    }

    MeshState compiled;
    {
      HS_PROFILE(hk_conway_compile);
      MeshOps::compile(swept, compiled, scratch_arena_a, scratch_arena_b);
    }
    HS_CHECK(compiled.face_counts.size() == tr.topo.size(),
             "ConwayMorph: sweep changed the compiled face count");
    compiled.topology.bind(scratch_arena_a, tr.topo.size());
    compiled.topology.append_bulk(tr.topo.data(), tr.topo.size());

    float w = blend_weight(static_cast<float>(frame) /
                           static_cast<float>(duration));
    BakedPalette *ramps =
        scratch_arena_b.allocate_n<BakedPalette>(tr.num_ramps);
    for (int r = 0; r < tr.num_ramps; ++r) {
      const BakedPalette &from = tr.bank->entries[tr.ramp_from[r]];
      const BakedPalette &to = tr.bank->entries[tr.ramp_to[r]];
      new (&ramps[r]) BakedPalette();
      if (w <= 0.0f) {
        ramps[r] = from;
      } else if (w >= 1.0f || tr.ramp_from[r] == tr.ramp_to[r]) {
        ramps[r] = to;
      } else {
        ramps[r].bake(scratch_arena_b, RampBlend{from, to, w});
      }
    }

    Shading sh{ramps, tr.face_ramp.data(), tr.face_ramp.size()};
    draw_fn(canvas, compiled, sh);
  }

  /**
   * @brief Arrival data for the effect's completion handler.
   * @return Arena-backed Landing; stable until the leg arena is compacted.
   */
  const Landing &landing() const { return buf_->landing; }

private:
  /**
   * @brief Arena-allocated leg state — keeps ConwayMorph inline size small.
   */
  struct Transients {
    PolyMesh seed;                 /**< Cloned leg seed. */
    ConwayGraph::MorphOp op = ConwayGraph::MorphOp::TRUNCATE; /**< Swept operator. */
    bool reverse = false;          /**< Traversing to_node -> from_node. */
    int sweep_frames = 1;          /**< Operator-sweep frames. */
    int settle_frames = 0;         /**< Relax-slerp frames. */
    float t_start = 0, t_end = 0;  /**< Clamped sweep endpoints. */
    float twist_start = 0, twist_end = 0; /**< Snub twist endpoints. */
    ArenaVector<Vector> relaxed;   /**< Relaxed endpoint vertices (settling legs). */
    ArenaVector<int> topo;         /**< Hoisted arrival classification. */
    ArenaVector<uint8_t> face_ramp; /**< Face -> (from, to) ramp pair. */
    const BakedPaletteBank *bank = nullptr; /**< Source LUTs for the pre-blend. */
    uint8_t ramp_from[MAX_BLEND_PAIRS] = {}; /**< Per-pair from palette. */
    uint8_t ramp_to[MAX_BLEND_PAIRS] = {};   /**< Per-pair to palette. */
    int num_ramps = 0;             /**< Distinct pair count. */
    Landing landing;               /**< Arrival data exposed to the effect. */
  };

  /** Blend source: two baked LUTs lerped at the frame's weight. */
  struct RampBlend {
    const BakedPalette &from; /**< Inherited (mapped) ramp. */
    const BakedPalette &to;   /**< Leg target ramp. */
    float w;                  /**< Blend weight in [0, 1]. */
    /**
     * @brief Samples both LUTs and lerps.
     * @param t Lookup coordinate.
     * @return The blended color at t.
     */
    Color4 get(float t) const { return from.get(t).lerp(to.get(t), w); }
  };

  /**
   * @brief Runs the edge's operator on the seed at one parameter value.
   */
  HS_COLD_MEMBER static PolyMesh run_op(ConwayGraph::MorphOp op,
                                        const PolyMesh &seed, Arena &target,
                                        Arena &temp, float t, float twist) {
    switch (op) {
    case ConwayGraph::MorphOp::TRUNCATE:
      return MeshOps::truncate(seed, target, temp, t);
    case ConwayGraph::MorphOp::EXPAND:
      return MeshOps::expand(seed, target, temp, t);
    default:
      return MeshOps::snub(seed, target, temp, t, twist);
    }
  }

  /**
   * @brief Crossfade weight (spec 2.6): exactly 0 through the first 20% of the
   * leg, smoothstep to exactly 1 by 80%.
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

  /**
   * @brief Builds the per-face from-palettes (leg-swap mapping, spec 2.5), the
   * shuffled target assignment (spec 2.6), and the distinct ramp-pair table.
   * @param tr Leg transients being populated.
   * @param arrival Classified arrival mesh (for face counts).
   * @param handoff Departed-node provenance.
   * @param arena Leg arena for the face -> ramp table.
   */
  HS_COLD_MEMBER void build_palette_mapping(Transients &tr,
                                            const PolyMesh &arrival,
                                            const PaletteHandoff &handoff,
                                            Arena &arena) {
    const size_t total = tr.topo.size();
    const size_t primary = tr.seed.face_counts.size();
    tr.landing.topology = tr.topo.data();
    tr.landing.faces = total;
    tr.landing.primary_faces = primary;

    for (int i = 0; i < PALETTES; ++i)
      tr.landing.to_palette[i] = static_cast<uint8_t>(i);
    std::shuffle(tr.landing.to_palette.begin(), tr.landing.to_palette.end(),
                 hs::random());

    if (!handoff.by_class_signature) {
      // Emission-order mapping: either the whole face list corresponds
      // (departing a mid-parameter node) or the primary prefix does
      // (departing the seed form; orbit/edge faces are zero-area births).
      HS_CHECK(handoff.prev_faces == total || handoff.prev_faces == primary,
               "ConwayMorph: handoff face count matches neither mapping");
    } else {
      HS_CHECK(handoff.prev_face_sides,
               "ConwayMorph: class-signature mapping needs prev side counts");
    }

    tr.face_ramp.bind(arena, total);
    for (size_t f = 0; f < total; ++f) {
      uint8_t to =
          tr.landing.to_palette[wrap(tr.topo[f], PALETTES)];
      uint8_t from = to; // newborn faces skip the crossfade
      if (handoff.by_class_signature) {
        // DUAL_SWAP: side count at the shared ambo point is unambiguous.
        int sides = f < primary ? tr.seed.face_counts[f]
                                : arrival.face_counts[f];
        for (size_t j = 0; j < handoff.prev_faces; ++j) {
          if (handoff.prev_face_sides[j] == sides) {
            from = handoff.prev_face_palette[j];
            break;
          }
        }
      } else if (f < handoff.prev_faces) {
        from = handoff.prev_face_palette[f];
      }
      HS_CHECK(from < PALETTES && to < PALETTES);

      int ramp = -1;
      for (int r = 0; r < tr.num_ramps; ++r) {
        if (tr.ramp_from[r] == from && tr.ramp_to[r] == to) {
          ramp = r;
          break;
        }
      }
      if (ramp < 0) {
        HS_CHECK(tr.num_ramps < MAX_BLEND_PAIRS,
                 "ConwayMorph: distinct palette pairs exceed the blend budget");
        ramp = tr.num_ramps++;
        tr.ramp_from[ramp] = from;
        tr.ramp_to[ramp] = to;
      }
      tr.face_ramp.push_back(static_cast<uint8_t>(ramp));
    }
  }

  Transients *buf_;   /**< Pointer to arena-allocated leg state. */
  EasingFn easing_fn; /**< Easing applied to the sweep parameter. */
  MorphDrawFn draw_fn; /**< Per-frame draw callback. */
};

} // namespace Animation

/**
 * @brief Compile-time segue policies for MeshCarousel: the styles by which one
 * mesh hands the sphere to the next.
 * @details A segue owns the scheduling shape of a mesh-to-mesh transition (its
 * schedule() hook, whose return value is the delay until the next transition
 * begins) and the meaning of the phase value the scheduled animation feeds the
 * draw functor: phase ramps 0→1 over the incoming window, holds 1, and falls
 * back to 0 over the outgoing window. The shading hooks translate that phase
 * into pixels:
 *
 *   bool   visible(phase)   — whether drawing at this phase is worthwhile
 *   float  opacity(phase)   — global alpha multiplier
 *   float  fill(&t, phase)  — coverage mask; may remap the edge-distance t
 *   Color4 grade(c, phase)  — color regrade after the palette lookup
 *
 * Optional hooks, detected per policy with `requires` so unused paths compile
 * out of the render loop:
 *
 *   void   retarget(v)               — re-randomize per-transition state
 *   Vector warp(v, phase)            — pre-ripple unit-sphere vertex warp
 *   float  face_offset(center, i, cls) — per-face sweep ordering in [0, 1]
 *   float  face_fade_frac(i)          — per-face fade length as a window fraction
 *   float  face_phase(phase, offset[, fade_frac]) — face-local phase from the front
 *
 * A per-face policy may also declare `static constexpr bool LOCAL_SWEEP =
 * true` to order faces by the untransformed mesh instead of world-space
 * centers: the front then rides the mesh's rotation rather than staying
 * fixed in the room.
 *
 * Policies are resolved at compile time (no virtuals); Base's identity hooks
 * inline to nothing.
 */
namespace Segue {

/**
 * @brief Schedules one sequential fade-in/fade-out sprite: consecutive sprites
 * never overlap, so a transition renders a single mesh per frame.
 * @param timeline Timeline receiving the sprite.
 * @param draw_fn Draws the mesh at the envelope phase.
 * @param duration Total frames the mesh is on screen.
 * @param window Requested transition window in frames; clamped to duration/2
 * so the in/out windows never collide.
 * @return duration — the next transition starts as this sprite ends.
 */
inline int schedule_sequential(Timeline &timeline, SpriteFn draw_fn,
                               int duration, int window) {
  int fade = std::min(window, duration / 2);
  timeline.add(0, Animation::Sprite(std::move(draw_fn), duration, fade,
                                    ease_linear, fade, ease_linear));
  return duration;
}

/**
 * @brief Soft sweep front used by Shockwave.
 * @param phase Global segue phase in [0, 1].
 * @param offset Face's sweep ordering in [0, 1]; larger extinguishes earlier.
 * @param band Softness of the front, in phase units.
 * @return The face-local phase in [0, 1]: 1 everywhere at phase 1, 0
 * everywhere at phase 0, with faces crossing the front in offset order.
 * @details The sqrt ease keeps the hand-off out of black: both meshes sit at low
 * phase around the swap, so a linear front would leave the sphere mostly dark;
 * accelerating the front through the low-phase end compresses that to a blink.
 * Endpoints stay exact (phase 1 remains the identity plateau).
 */
inline float sweep_phase(float phase, float offset, float band) {
  float p = std::sqrt(phase);
  return hs::clamp((p - offset * (1.0f - band)) / band, 0.0f, 1.0f);
}

/**
 * @brief Identity hooks every segue inherits; a policy shadows only the hooks
 * its transition uses.
 */
struct Base {
  /** @brief Default scheduling: one sequential sprite (see schedule_sequential). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    return schedule_sequential(timeline, std::move(draw_fn), duration, window);
  }
  /** @brief Whether drawing at this phase can produce visible output. */
  bool visible(float phase) const { return phase > 0.005f; }
  /** @brief Global alpha at this phase. */
  float opacity(float) const { return 1.0f; }
  /**
   * @brief Coverage mask over the face interior.
   * @param t Fragment edge distance in [0, 1] (0 at the edge, ~1 at the face
   * center); may be remapped in place for the palette lookup.
   * @param phase Transition phase; unused by the identity policy.
   * @return Coverage alpha in [0, 1]; 0 culls the fragment.
   */
  float fill(float &t, float phase) const {
    (void)t;
    (void)phase;
    return 1.0f;
  }
  /** @brief Color regrade applied after the palette lookup. */
  Color4 grade(Color4 c, float) const { return c; }
};

/**
 * @brief Opacity cross-fade between consecutive meshes.
 * @details Phase is opacity. Each transition is one fade-in/fade-out Sprite;
 * the returned delay starts the next transition one fade window before this
 * sprite ends, so the outgoing and incoming sprites overlap and both meshes
 * render during the fade — the cost of this segue is two rasterized meshes
 * per overlap frame. Every other segue is sequential (single mesh per frame).
 */
struct Crossfade : Base {
  /**
   * @brief Schedules the incoming mesh's fading sprite.
   * @param timeline Timeline receiving the sprite.
   * @param draw_fn Draws the incoming mesh at the given opacity.
   * @param duration Total frames the mesh is on screen.
   * @param window Requested fade length in frames; clamped to duration/2 so
   * the fade windows never overlap and sprites cannot pile up beyond two.
   * @return Frames after which the next transition should begin.
   */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    int fade = std::min(window, duration / 2);
    timeline.add(0, Animation::Sprite(std::move(draw_fn), duration, fade,
                                      ease_linear, fade, ease_linear));
    return duration - fade;
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief Faces contract to glowing points at their centers, then the next
 * tessellation blooms back out of the point field.
 * @details Only fragments deeper than the phase-driven inset survive, so the
 * pattern dissolves into a constellation of face-center glints at the swap.
 * The surviving core's edge distance is renormalized so it keeps the full
 * palette gradient as it shrinks.
 */
struct IrisBloom : Base {
  static constexpr float SOFT = 0.08f; /**< Soft rim width, in edge-distance units. */
  float fill(float &t, float phase) const {
    float inset = 1.0f - phase;
    if (t < inset - SOFT)
      return 0.0f;
    float cover = hs::clamp((t - (inset - SOFT)) / SOFT, 0.0f, 1.0f);
    t = hs::clamp((t - inset) / std::max(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief The fill drains until only a glowing band along the edges survives,
 * the meshes swap as lace, then the new fill floods back in.
 * @details The inverse mask of IrisBloom: fragments within the phase-driven
 * band of an edge survive. A line network changing shape reads far less
 * jarring than filled regions changing, which hides the swap.
 */
struct Lace : Base {
  static constexpr float SOFT = 0.08f; /**< Soft band-edge width, in edge-distance units. */
  float fill(float &t, float phase) const {
    if (t > phase + SOFT)
      return 0.0f;
    float cover = hs::clamp((phase + SOFT - t) / SOFT, 0.0f, 1.0f);
    t = hs::clamp(t / std::max(phase, 1e-3f), 0.0f, 1.0f);
    return cover;
  }
};

/**
 * @brief A day/night line pinned to the mesh sweeps across it; when it reaches a
 * face, that face fades over a per-face random length in [fade_frames_min,
 * fade_frames_max] frames.
 * @details LOCAL_SWEEP anchors the line to the untransformed mesh. Each face's
 * fade length is a stable per-transition hash of its index, so the front frays
 * into an irregular edge. The front crosses over the fade window minus one face
 * fade, so face phases are exactly 1 at phase 1 and 0 at phase 0 for every fade
 * length.
 */
struct TerminatorSweep : Base {
  static constexpr bool LOCAL_SWEEP = true; /**< Sweep in mesh-local space. */
  Vector axis = Y_AXIS;          /**< Mesh-local sweep axis. */
  float fade_frames_min = 4.0f;  /**< Shortest per-face fade length, in frames. */
  float fade_frames_max = 12.0f; /**< Longest per-face fade length, in frames. */
  float fade_frac_min = 0.06f; /**< fade_frames_min over the scheduled window; set by schedule(). */
  float fade_frac_max = 0.17f; /**< fade_frames_max over the scheduled window; set by schedule(). */
  uint32_t fade_seed = 0x9e3779b9u; /**< Per-transition seed for the per-face fade hash; rolled by retarget(). */
  int schedule(Timeline &timeline, SpriteFn draw_fn, int duration, int window) {
    int fade = std::min(window, duration / 2);
    float inv = 1.0f / static_cast<float>(std::max(fade, 1));
    fade_frac_min = std::min(1.0f, fade_frames_min * inv);
    fade_frac_max = std::min(1.0f, fade_frames_max * inv);
    return schedule_sequential(timeline, std::move(draw_fn), duration, window);
  }
  void retarget(const Vector &v) {
    axis = v;
    fade_seed = static_cast<uint32_t>(hs::random()());
  }
  float face_offset(const Vector &center, int, int) const {
    return 0.5f * (1.0f + dot(center, axis));
  }
  /** @brief Per-face fade length as a window fraction: a stable hash of the
   * face index into [fade_frac_min, fade_frac_max]. Computed once per face (not
   * per fragment), so it must stay a pure function of the index and seed. */
  float face_fade_frac(int i) const {
    float t = hash01(static_cast<uint32_t>(i), fade_seed);
    return fade_frac_min + (fade_frac_max - fade_frac_min) * t;
  }
  float face_phase(float phase, float offset, float fade_frac) const {
    float ff = std::max(fade_frac, 1e-4f);
    return hs::clamp((phase - offset * (1.0f - ff)) / ff, 0.0f, 1.0f);
  }
  /** @brief Squared: alpha scales linear-light color, where a linear ramp
   * reads mostly-bright almost immediately; the square spreads the perceived
   * fade across the face's fade window. */
  float opacity(float phase) const { return phase * phase; }
};

/**
 * @brief An expanding shockwave erases the pattern outward from a point; its
 * echo redraws the new one.
 * @details Faces nearest the origin extinguish first, so the wave visibly
 * expands. Pairs naturally with the effect's ripple bursts sharing the origin.
 */
struct Shockwave : Base {
  static constexpr float BAND = 0.3f; /**< Wave-front softness, in phase units. */
  Vector origin = Y_AXIS; /**< World-space wave origin. */
  void retarget(const Vector &v) { origin = v; }
  float face_offset(const Vector &center, int, int) const {
    float angle = fast_acos(hs::clamp(dot(center, origin), -1.0f, 1.0f));
    return 1.0f - angle * (1.0f / PI_F);
  }
  float face_phase(float phase, float offset) const {
    return sweep_phase(phase, offset, BAND);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The pattern breaks down one topology class at a time: all faces of a
 * class fade together, each class fully gone before the next starts, in a
 * random class order reshuffled per transition; the new tessellation
 * reassembles class by class the same way.
 * @details Faces group by palette-slot class, so each color family vanishes as
 * a unit. Class windows are abutting equal slices of the phase range (linear,
 * not sweep_phase's eased front). The BLACK_DWELL slice nearest the swap is held
 * fully black so the last class completes before the incoming mesh appears,
 * instead of popping. reorder() derives the class count from the per-face
 * classes, so it can never disagree with the mesh.
 */
struct Breakdown : Base {
  static constexpr int MAX_CLASSES = 16; /**< rank[] capacity. */
  static constexpr float BLACK_DWELL = 0.1f; /**< Phase slice held all-black at the swap end. */
  int num_classes = 1;            /**< Live class count, derived by reorder(). */
  uint8_t rank[MAX_CLASSES] = {}; /**< rank[class]: fade position; 0 vanishes first. */
  /**
   * @brief Derives the class count from the per-face classes and re-randomizes
   *        the fade order for the next transition.
   * @param face_classes Per-face class ids (dense [0, k), the same values
   *        face_offset receives). num_classes is set to max+1, so a caller can
   *        never mis-declare it. A face class at or past MAX_CLASSES folds into
   *        rank[0] (face_offset's out-of-range branch).
   */
  template <typename Classes>
  void reorder(const Classes &face_classes) {
    int detected = 1;
    for (size_t i = 0; i < face_classes.size(); ++i) {
      int c = static_cast<int>(face_classes[i]) + 1;
      if (c > detected)
        detected = c;
    }
    num_classes = hs::clamp(detected, 1, MAX_CLASSES);
    for (int i = 0; i < num_classes; ++i)
      rank[i] = static_cast<uint8_t>(i);
    std::shuffle(rank, rank + num_classes, hs::random());
  }
  float face_offset(const Vector &, int, int cls) const {
    if (num_classes <= 1)
      return 0.0f;
    int r = rank[(cls >= 0 && cls < num_classes) ? cls : 0];
    return static_cast<float>(num_classes - 1 - r) /
           static_cast<float>(num_classes - 1);
  }
  float face_phase(float phase, float offset) const {
    // Class windows tile [BLACK_DWELL, 1]; phase 1 stays the identity plateau.
    float band = (1.0f - BLACK_DWELL) / static_cast<float>(num_classes);
    return hs::clamp(
        (phase - BLACK_DWELL - offset * (1.0f - BLACK_DWELL - band)) / band,
        0.0f, 1.0f);
  }
  float opacity(float phase) const { return phase; }
};

/**
 * @brief The whole mesh spins up around an axis until the POV display smears
 * it into bands, swaps at peak speed, and spins back down — a coin flip.
 * @details The warp is rigid, so there is no fold/overdraw hazard and the mesh
 * never fades; the swap hides in the motion blur.
 */
struct SpinFlip : Base {
  static constexpr float REVS = 3.0f; /**< Extra revolutions at peak spin. */
  Vector axis = Y_AXIS; /**< Spin axis. */
  void retarget(const Vector &v) { axis = v; }
  Vector warp(const Vector &v, float phase) const {
    float wind = 1.0f - phase;
    return rotate(v, make_rotation(axis, wind * wind * REVS * 2.0f * PI_F));
  }
};

/**
 * @brief Both palettes converge to molten gold around the swap, then the new
 * mesh blooms back into color.
 * @details Purely palette-domain: geometry never moves. A mild opacity dip at
 * the swap softens the topology pop while both meshes are monochrome.
 */
struct GoldConvergence : Base {
  Pixel gold = Color4(uint8_t{255}, uint8_t{196}, uint8_t{64}).color; /**< Linear-space convergence color. */
  Color4 grade(Color4 c, float phase) const {
    return c.lerp(Color4(gold, c.alpha), 1.0f - phase);
  }
  float opacity(float phase) const { return 0.4f + 0.6f * phase; }
};

} // namespace Segue

/**
 * @brief A double-buffered pair of persistent MeshState slots, the
 *        arena-compaction primitives effects need to swap between them, and a
 *        pluggable compile-time segue.
 * @tparam SegueT Segue policy (see namespace Segue) behind schedule_segue().
 * Clients that run their own transition animations (e.g. a `MeshMorph`) keep
 * the default and simply never call it.
 * @details Holds two MeshState slots in `persistent_arena` and a front/back
 * index. Effects own generation and drawing (generate into a slot, flip the
 * front index, reclaim the old slot); the segue owns transition scheduling.
 *
 * Usage:
 *   MeshCarousel<Segue::Crossfade> carousel;  // in effect members
 *
 *   // Build the initial shape directly into the front slot:
 *   carousel.current().clear();
 *   MeshOps::compile(mesh, carousel.current(), persistent_arena, scratch_arena_a);
 *
 *   // To transition: generate into the back slot, flip, then let the segue
 *   // schedule the animation via schedule_segue (see
 *   // IslamicStars::spawn_shape for the pattern).
 */
template <typename SegueT = Segue::Crossfade> class MeshCarousel {
public:
  /**
   * @brief Constructs an empty carousel with front slot index 0.
   */
  MeshCarousel() {}

  /**
   * @brief Gets the currently visible (front) mesh.
   * @return Const reference to the front MeshState slot.
   */
  const MeshState &current() const { return slots_[front_]; }

  /**
   * @brief Gets the currently visible (front) mesh (mutable).
   * @return Mutable reference to the front MeshState slot.
   */
  MeshState &current() { return slots_[front_]; }

  /**
   * @brief Direct slot access by index (for effects that need both).
   * @param i Slot index (0 or 1).
   * @return Const reference to the requested MeshState slot.
   */
  const MeshState &slot(int i) const {
    HS_CHECK(i == 0 || i == 1, "MeshCarousel slot index out of range");
    return slots_[i];
  }

  /**
   * @brief Direct slot access by index (mutable).
   * @param i Slot index (0 or 1).
   * @return Mutable reference to the requested MeshState slot.
   */
  MeshState &slot(int i) {
    HS_CHECK(i == 0 || i == 1, "MeshCarousel slot index out of range");
    return slots_[i];
  }

  /**
   * @brief Gets the front slot index (for capture in lambdas).
   * @return The index (0 or 1) of the front slot.
   */
  int front_index() const { return front_; }

  /**
   * @brief Manually sets the front index (for effects that manage transitions
   * themselves).
   * @param idx The new front slot index (0 or 1).
   */
  void set_front(int idx) {
    HS_CHECK(idx == 0 || idx == 1, "MeshCarousel front index out of range");
    front_ = idx;
  }

  /**
   * @brief Schedules the segue's transition animation for the (already
   * front-flipped) incoming mesh.
   * @param timeline Timeline receiving the segue's animation.
   * @param draw_fn Draws the incoming mesh; the float argument is the segue's
   * phase (opacity for Segue::Crossfade).
   * @param duration Total frames the mesh is on screen.
   * @param window Transition window length in frames, segue-interpreted.
   * @return Frames after which the effect should begin the next transition.
   */
  int schedule_segue(Timeline &timeline, SpriteFn draw_fn, int duration,
                     int window) {
    return segue_.schedule(timeline, std::move(draw_fn), duration, window);
  }

  /**
   * @brief The carousel's segue policy instance (holds per-transition state
   * such as a sweep axis or wave origin).
   */
  SegueT &segue() { return segue_; }
  /** @brief Const view of the segue policy instance. */
  const SegueT &segue() const { return segue_; }

  /**
   * @brief Compacts the persistent arena, reclaiming fragmented space.
   * @details Evacuates tracked MeshStates and resets the arena. Call before
   * allocating new persistent data.
   */
  void compact() {
    // Both evacuations share scratch_arena_a, which must hold both populated slots.
    Persist<MeshState> p0(slots_[0], scratch_arena_a, persistent_arena);
    Persist<MeshState> p1(slots_[1], scratch_arena_a, persistent_arena);
    persistent_arena.reset();
  }

  /**
   * @brief Frees the back slot and compacts, preserving only the front slot.
   * @tparam AfterReset Callable type invoked as `void(Arena&)`.
   * @param after_reset Callback run immediately after the reset, while the front
   * slot is still evacuated.
   * @details Runs `after_reset(persistent_arena)` immediately after the reset —
   * while the front slot is still evacuated — so the caller can re-bake
   * effect-owned persistent data (e.g. a palette bank) into the fresh arena
   * *before* the front mesh is restored on top of it. Use when only the visible
   * (front) shape must survive a regeneration of the back slot.
   */
  template <typename AfterReset> void compact_keep_front(AfterReset after_reset) {
    int back = 1 - front_;
    slots_[back] = MeshState();
    Persist<MeshState> p(slots_[front_], scratch_arena_b, persistent_arena);
    persistent_arena.reset();
    after_reset(persistent_arena);
  }

private:
  MeshState slots_[2]; /**< Front/back double-buffered mesh slots. */
  int front_ = 0;      /**< Index (0 or 1) of the visible front slot. */
  SegueT segue_;       /**< Segue policy instance; per-transition state lives here. */
};
