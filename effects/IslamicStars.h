/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/animation/animation.h"
#include "core/engine/engine.h"
#include "core/mesh/recipe.h"

// Unit-test accessor reaching the private build-chain state (pre-init trans
// speed, build_active_, solid_idx) for the effects-module build smoke.
namespace hs_test {
namespace effects_tests {
struct IslamicBuildProbe;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Effect that displays a sequence of Islamic-geometry polyhedra,
 *        transitioning one shape into the next while ripples distort the mesh.
 *        Entries with a non-null recipe are built op by op on screen: the
 *        recipe's seed solid sweeps in, one OpLeg per lowered primitive step
 *        morphs it into the finished pattern, then the usual still/ripple/fade
 *        choreography runs (docs/opchain_morph_spec.md, section 6).
 * @tparam W Target canvas width in pixels.
 * @tparam H Target canvas height in pixels.
 */
template <int W, int H> class IslamicStars : public Effect {

public:
  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  HS_COLD_MEMBER IslamicStars()
      : Effect(W, H, {.strobe = true}), filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  void init() override {
    // Asymmetric scratch split (194 KB total): the leg-by-leg build chain
    // peaks at ~118 KB in a and ~71 KB in b, and compact_keep_front evacuates
    // the front slot (up to 63.7 KB) through b. The remainder is persistent:
    // carousel slots + BakedPaletteBank (~15 KB). Budgets enforced by the
    // test_conway_morph.h build replay and test_solids.h's high-water sweeps.
    configure_arenas(GLOBAL_ARENA_SIZE - (114 + 80) * 1024, 120 * 1024,
                     74 * 1024);

    ripple_gen.init_storage(persistent_arena);
    palette_bank_.bake_all(persistent_arena);

    // Set BEFORE registering: register_param snaps *ptr as the slider default.
    // Amplitude starts at the fold-free ceiling (see RIPPLE_AMP_MAX).
    ripple_gen.template_params.amplitude = RIPPLE_AMP_MAX;
    ripple_gen.template_params.thickness = RIPPLE_THICKNESS;
    ripple_gen.template_params.decay = 0.1f;

    register_param("Fade", &params.fade, 0.0f, 96.0f);
    // Per-face fade length range (frames): each face draws a random fade from
    // [lo, hi] as the terminator reaches it, fraying the sweep front.
    register_param("Face Fade Lo", &carousel.segue().fade_frames_min, 0.0f,
                   32.0f);
    register_param("Face Fade Hi", &carousel.segue().fade_frames_max, 0.0f,
                   32.0f);
    // Burst/Ripp Dur ranges are clamped to the ripple pool capacity invariant
    // (see the RIPPLE* constants below).
    register_param("Burst", &params.burst_size, 1.0f, (float)BURST_MAX);
    // Amplitude slider capped at the fold-free ceiling; thickness is fixed (not a
    // slider) so amplitude/thickness can never cross the self-fold onset.
    register_param("Ripp Amp", &ripple_gen.template_params.amplitude, 0.0f,
                   RIPPLE_AMP_MAX);
    register_param("Ripp Decay", &ripple_gen.template_params.decay, 0.0f, 5.0f);
    register_param("Ripp Dur", &ripple_duration, 30.0f,
                   (float)RIPPLE_DURATION_MAX);
    register_param("Trans Speed", &params.trans_speed, 1.0f, 8.0f);
    register_param("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

    // Open on a recipe entry so the op-by-op build is the first thing drawn;
    // spawn_shape pre-increments, so seed the index one before it.
    auto solids = Solids::Collections::get_islamic_solids();
    for (size_t i = 0; i < solids.size(); ++i) {
      if (solids[i].recipe) {
        solid_idx = static_cast<int>(i) - 1;
        break;
      }
    }

    spawn_shape();
  }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(is_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(is_ripple_prepare);
      ripple_gen.prepare_frame();
    }
    {
      HS_PROFILE(is_timeline_step);
      timeline.step(canvas);
    }
  }

private:
  friend struct ::hs_test::effects_tests::IslamicBuildProbe;

  // Ripple-pool sizing: a slot is held from spawn() until the staggered ripple
  // completes. Only one burst is normally live; the pool holds two so a Ripp
  // Dur/Burst slider change mid-burst cannot drop a spawn.
  static constexpr int RIPPLE_POOL_SIZE = 8;
  static constexpr int RIPPLE_STAGGER_FRAMES = 16;
  static constexpr int RIPPLE_DURATION_MAX = 143;
  static constexpr int BURST_MAX = 4;
  static constexpr int STILL_FRAMES =
      16; /**< 1 s hold (16 fps) between fade and ripple stages. */
  // Recipe-build leg lengths (docs/opchain_morph_spec.md section 7); divided
  // by the Trans Speed divisor like every other stage.
  static constexpr int HANKIN_LEG_FRAMES = 32;
  static constexpr int SWEEP_LEG_FRAMES =
      24; /**< ambo / truncate / snub / chamfer. */
  static constexpr int RELAX_LEG_FRAMES = 16;
  static constexpr int GATE_LEG_FRAMES = 13; /**< kis / dual: 6 + 1 + 6. */
  static constexpr size_t MAX_BUILD_STEPS = 8; /**< Lowered-primitive cap. */
  /** Build-chain mesh face cap. Bounds the scratch handoff arrays only; the
   * persistent budget is what actually limits which recipes ship. */
  static constexpr size_t MAX_BUILD_FACES = 1152;
  static constexpr float RIPPLE_THICKNESS =
      0.7f; /**< Fixed ripple wavelet width (radians). */
  static constexpr float RIPPLE_AMP_MAX =
      0.15f; /**< Fold-free amplitude ceiling at RIPPLE_THICKNESS (amp/thickness < ~0.2 self-fold onset). */
  static_assert(2 * BURST_MAX <= RIPPLE_POOL_SIZE,
                "IslamicStars: ripple pool must hold two overlapping bursts");

  Orientation<> orientation;
  Timeline timeline;
  Pipeline<W, H> filters;
  RippleTransformer<RIPPLE_POOL_SIZE> ripple_gen;
  FastNoiseLite noise;
  float ripple_duration = 80.0f;
  // Effective per-shape stage lengths after the Trans Speed divisor, written by
  // spawn_shape and read by the deferred ripple() callback.
  int ripple_dur_eff_ = 80;
  int ripple_stagger_eff_ = RIPPLE_STAGGER_FRAMES;
  int solid_idx = -1;
  using SegueT = Segue::TerminatorSweep;
  MeshCarousel<SegueT> carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank_;
  // Per-slot palette indices; value-init so a missed shuffle reads palette 0,
  // not garbage.
  std::array<int, NUM_PALETTES> palettes_slots[2] = {};

  // Build-chain state (entries with a non-null recipe): the shape is built op
  // by op between the fade-in and the still hold. Null-recipe entries never
  // touch any of it.
  bool build_active_ = false; /**< Legs draw; the sprite draw_fn is muted. */
  Solids::OpStep build_steps_[MAX_BUILD_STEPS]; /**< Lowered primitive chain. */
  size_t build_step_count_ = 0; /**< Lowered step count. */
  size_t build_step_ = 0;       /**< Current leg index. */
  int build_leg_frames_[MAX_BUILD_STEPS] = {}; /**< Per-leg frame budget. */
  int build_total_frames_ = 0;                 /**< Sum of leg frames. */
  PolyMesh build_seed_;       /**< Leg-k seed (persistent). */
  PolyMesh build_next_seed_;  /**< Clean endpoint seed_{k+1}: built eagerly at
                                 leg start, or from the leg's own topology at
                                 its end on a hankin step. */
  std::array<uint8_t, Animation::OpLeg::PALETTES> build_targets_ =
      {}; /**< Per-shape pinned target palettes, shuffled before leg 0. */
  const Animation::OpLeg::Landing *build_landing_ =
      nullptr; /**< Latest leg's arrival data (leg-arena backed). */
  const uint8_t *build_from_pal_ =
      nullptr; /**< Per-face palette the previous leg landed on; survives the
                  leg-boundary compaction that drops its landing. */
  size_t build_from_faces_ = 0; /**< Length of build_from_pal_. */

  /**
   * @brief Draw callback for build-leg frames.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, const MeshState &, const Animation::OpLeg::Shading &), 8>
      draw_build_fn_{[this](Canvas &c, const MeshState &m,
                            const Animation::OpLeg::Shading &sh) {
        draw_build_mesh(c, m, sh);
      }};

  /**
   * @brief Spawns one burst of burst_size ripples from a random origin,
   *        staggered RIPPLE_STAGGER_FRAMES apart, each expanding over ripple_duration
   *        frames.
   * @param canvas Unused render target for the timer callback signature.
   */
  void ripple(Canvas &) {
    Vector origin = random_vector();
    for (int i = 0; i < (int)params.burst_size; i++) {
      if (!ripple_gen.spawn(i * ripple_stagger_eff_, origin,
                            PI_F / ripple_dur_eff_, ripple_dur_eff_))
        hs::log("IslamicStars: ripple pool full, dropping spawn");
    }
  }

  /**
   * @brief Orients and ripple-distorts a mesh into scratch_arena_a.
   * @param base_state Undistorted source mesh.
   * @return The transformed mesh (scratch_arena_a-backed; the caller holds the
   * scope). Shared by the sprite and build-leg draw paths so build frames ride
   * the exact transform chain the held shape uses.
   */
  HS_O3_FN MeshState transform_shape(const MeshState &base_state) {
    MeshState transformed_state;
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
      MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                         ripple_gen, camera);
    }
    return transformed_state;
  }

  /**
   * @brief Orients, ripple-distorts, and segue-shapes base_state, then
   *        rasterizes it with a per-face palette lookup.
   * @param canvas Render target receiving the rasterized mesh.
   * @param phase Segue phase in [0, 1] from the sprite envelope: rises over
   *        the incoming window, holds 1, falls over the outgoing window.
   * @param base_state Undistorted source mesh to transform and draw.
   * @param face_indices Maps each face to its topology class.
   * @param palette_idx Assigns a palette per topology class.
   * @note Draws on the exact SDF path, not the congruence-class LUT
   * (mesh_classes.h): ripple/segue deformation makes a canonical LUT mis-shade
   * or pop. The facility is for effects whose meshes hold still.
   */
  HS_O3_FN void draw_shape(Canvas &canvas, float phase,
                           const MeshState &base_state,
                           const ArenaVector<int> &face_indices,
                           const std::array<int, NUM_PALETTES> &palette_idx) {
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase))
      return;
    HS_PROFILE(is_draw_shape);
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state = transform_shape(base_state);

    const int *raw_indices = face_indices.data();
    const int num_faces = static_cast<int>(face_indices.size());

    // Per-face segues order faces by their center, recomputed per frame: from
    // world space by default (the front stays fixed in the room while the
    // mesh rotates through it), or from the untransformed mesh for segues
    // declaring LOCAL_SWEEP (the front rides the mesh). The third argument is
    // the face's palette-slot class, mapped exactly as the fragment shader
    // maps it; class-agnostic sweeps ignore it.
    constexpr bool PER_FACE =
        requires(const Vector &c) { seg.face_offset(c, 0, 0); };
    // phase is fixed for the whole draw call, so the segue's per-face phase
    // resolves here rather than per fragment; so does the face's palette, whose
    // slot is the same class that face_offset already needs.
    ArenaVector<float> face_phases;
    ArenaVector<const BakedPalette *> face_palettes;
    if constexpr (PER_FACE) {
      HS_PROFILE(is_face_offsets);
      constexpr bool LOCAL_SWEEP = requires { requires SegueT::LOCAL_SWEEP; };
      const MeshState &sweep_state =
          LOCAL_SWEEP ? base_state : transformed_state;
      const size_t faces = sweep_state.num_faces();
      face_phases.bind(scratch_arena_a, faces);
      face_palettes.bind(scratch_arena_a, faces);
      const uint16_t *fidx = sweep_state.get_faces_data();
      const uint16_t *foff = sweep_state.get_face_offsets_data();
      const uint8_t *fcnt = sweep_state.get_face_counts_data();
      for (size_t f = 0; f < faces; ++f) {
        Vector c(0.0f, 0.0f, 0.0f);
        for (int k = 0; k < fcnt[f]; ++k)
          c = c + sweep_state.vertices[fidx[foff[f] + k]];
        int cls = (f < static_cast<size_t>(num_faces))
                      ? wrap(raw_indices[f], NUM_PALETTES)
                      : 0;
        float off =
            seg.face_offset(normalized_or(c, UP), static_cast<int>(f), cls);
        float fade = seg.face_fade_frac(static_cast<int>(f));
        face_phases.push_back(seg.face_phase(phase, off, fade));
        face_palettes.push_back(&palette_bank_[palette_idx[cls]]);
      }
    }

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      if constexpr (PER_FACE) {
        int fi = static_cast<int>(frag.v2);
        if (fi >= 0 && fi < static_cast<int>(face_phases.size())) {
          frag.color = shade_mesh_topology(frag, *face_palettes[fi], 1.0f, seg,
                                           face_phases[fi]);
          return;
        }
      }
      frag.color =
          shade_mesh_topology(frag, raw_indices, num_faces, palette_bank_,
                              palette_idx, 1.0f, seg, phase);
    };

    {
      HS_PROFILE(is_mesh_scan);
      Scan::Mesh::draw<W, H>(filters, canvas, transformed_state,
                             fragment_shader, scratch_arena_a, params.debug_bb);
    }
  }

  /**
   * @brief Draws one build-leg frame: the compiled swept mesh, shaded from its
   *        pre-blended palette ramps on the segue's phase-1 identity plateau.
   * @param canvas Render target receiving the rasterized mesh.
   * @param mesh Compiled swept mesh (scratch-backed, this frame only).
   * @param sh Per-face blended-ramp table from the OpLeg.
   */
  HS_O3_FN void draw_build_mesh(Canvas &canvas, const MeshState &mesh,
                                const Animation::OpLeg::Shading &sh) {
    if (mesh.vertices.is_empty())
      return;
    // Own scope labels: sharing the sprite's would parent two draw paths under
    // one counter, and a build-only window then prints an empty subtree while a
    // mixed window prints the child above its own parent's total.
    HS_PROFILE(is_build_draw);
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state = transform_shape(mesh);
    const SegueT &seg = carousel.segue();

    auto fragment_shader = [&](const Vector &, Fragment &frag) {
      int fi = static_cast<int>(frag.v2);
      int ramp = (fi >= 0 && fi < static_cast<int>(sh.faces))
                     ? sh.face_ramp[fi]
                     : 0;
      frag.color = shade_mesh_topology(frag, sh.ramps[ramp], sh.gain, seg, 1.0f);
    };

    {
      HS_PROFILE(is_build_scan);
      Scan::Mesh::draw<W, H>(filters, canvas, transformed_state,
                             fragment_shader, scratch_arena_a, params.debug_bb);
    }
  }

  /**
   * @brief Leg frame budget of one lowered primitive step, before the Trans
   *        Speed divisor (docs/opchain_morph_spec.md, section 7).
   * @param op Lowered primitive op.
   * @return Frames the leg runs for.
   */
  static int leg_frames(Solids::Op op) {
    switch (op) {
    case Solids::Op::HANKIN:
      return HANKIN_LEG_FRAMES;
    case Solids::Op::RELAX:
      return RELAX_LEG_FRAMES;
    case Solids::Op::KIS:
    case Solids::Op::DUAL:
      return GATE_LEG_FRAMES;
    default:
      return SWEEP_LEG_FRAMES;
    }
  }

  /**
   * @brief Whether a lowered primitive step runs as a gated swap.
   * @param op Lowered primitive op.
   * @return True for the partition ops.
   */
  static bool is_gated_step(Solids::Op op) {
    return op == Solids::Op::KIS || op == Solids::Op::DUAL;
  }

  /**
   * @brief Half-gate length of a gated leg's frame budget.
   * @param frames Budgeted leg frames after the Trans Speed divisor.
   * @return F_gate, at least 1; the leg then runs 2 * F_gate + 1 frames.
   */
  static int gate_frames(int frames) { return std::max(1, (frames - 1) / 2); }

  /**
   * @brief Log label of a lowered primitive step.
   * @param op Lowered primitive op.
   * @return Static name string.
   */
  static const char *leg_name(Solids::Op op) {
    switch (op) {
    case Solids::Op::HANKIN:
      return "hankin";
    case Solids::Op::AMBO:
      return "ambo";
    case Solids::Op::TRUNCATE:
      return "truncate";
    case Solids::Op::SNUB:
      return "snub";
    case Solids::Op::CHAMFER:
      return "chamfer";
    case Solids::Op::RELAX:
      return "relax";
    case Solids::Op::KIS:
      return "kis";
    case Solids::Op::DUAL:
      return "dual";
    default:
      return "?";
    }
  }

  /**
   * @brief Advances to the next solid, generates it into the carousel's back
   *        slot with a freshly shuffled palette, makes it the front, schedules
   *        the segue and the shape's mid-display ripple burst, and queues the
   *        next spawn_shape call.
   */
  HS_COLD_MEMBER void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
    solid_idx = (solid_idx + 1) % solids.size();
    int back = 1 - carousel.front_index();
    MeshPaletteBank::shuffle_indices(palettes_slots[back]);

    int idx = solid_idx;
    // A recipe whose lowered chain contains a step no leg kind covers falls
    // back to today's whole-generate path, seed solid and all.
    const Solids::Recipe *recipe = solids[idx].recipe;
    if (recipe) {
      build_step_count_ =
          Solids::expand_to_primitives(*recipe, build_steps_, MAX_BUILD_STEPS);
      for (size_t k = 0; k < build_step_count_; ++k) {
        if (!Solids::is_morphable_step(build_steps_[k])) {
          hs::log("IslamicStars: %s has an unsweepable step, generating whole",
                  solids[idx].name);
          recipe = nullptr;
          build_step_count_ = 0;
          break;
        }
      }
    }

    auto draw_fn = [this, back](Canvas &canvas, float phase) {
      // During the build window an OpLeg draws the mesh; exactly one mesh per
      // frame.
      if (build_active_)
        return;
      const MeshState &mesh = carousel.slot(back);
      this->draw_shape(canvas, phase, mesh, mesh.topology,
                       palettes_slots[back]);
    };

    // Compact the back slot, rebaking palettes into the fresh arena rather than
    // tracking them through the evacuation. The ripple pool re-claims first so
    // its slots re-land at their init_storage() addresses.
    //
    // A build regenerates both slots before either is drawn again — the seed
    // below, the finished solid at finish_build — so the outgoing shape is
    // dropped rather than evacuated. Keeping it costs the whole build the
    // previous shape's mesh, which at the roster's 1082-face entry is most of
    // the persistent arena.
    auto rebake = [this](Arena &arena) {
      ripple_gen.reclaim_storage(arena);
      palette_bank_.bake_all(arena);
    };
    if (recipe)
      carousel.compact_drop_all(rebake);
    else
      carousel.compact_keep_front(rebake);

    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      if (recipe) {
        // The build starts from the recipe's seed solid; the chain is swept
        // on screen leg by leg. The seed is also held as the first leg's
        // persistent PolyMesh.
        build_seed_ = Solids::finalize_solid(
            Solids::simple_registry[recipe->seed].generate(a, b), target);
        carousel.slot(back).clear();
        MeshOps::compile(build_seed_, carousel.slot(back), target,
                         scratch_arena_a);
      } else {
        PolyMesh mesh = solids[idx].generate(a, b);
        carousel.slot(back).clear();
        MeshOps::compile(mesh, carousel.slot(back), target, scratch_arena_a);
      }
    });

    // ScratchScope frees only this call's own allocations, preserving prior
    // caller allocations in these shared arenas that a bare reset() would drop.
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(carousel.slot(back), scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Flip front eagerly for the overlapping sprite.
    carousel.set_front(back);

    // Segues with a spatial anchor (sweep axis, wave origin, spin axis) get a
    // fresh random one per transition. Safe mid-carousel: those segues are
    // sequential, so the previous sprite has already finished.
    if constexpr (requires(SegueT &s, const Vector &v) { s.retarget(v); })
      carousel.segue().retarget(random_vector());

    // Per-shape choreography: segue in, hold still one second, ripple, settle
    // one second, segue out. Duration is derived from the stage lengths so the
    // stages never overlap; the segue warps are identity on the phase-1
    // plateau, so the mesh only moves during its own stage.
    // Trans Speed divides every stage length so the carousel can be sped up
    // (e.g. for profiling) without touching the shape geometry. Each stage keeps
    // a >=1-frame floor. The effective ripple duration/stagger are cached for the
    // deferred ripple() callback, which fires before the next shape spawns.
    const float sp = std::max(1.0f, params.trans_speed);
    int fade = std::max(1, static_cast<int>(params.fade / sp));
    int still = std::max(1, static_cast<int>(STILL_FRAMES / sp));
    ripple_dur_eff_ = std::max(8, static_cast<int>(ripple_duration / sp));
    ripple_stagger_eff_ =
        std::max(1, static_cast<int>(RIPPLE_STAGGER_FRAMES / sp));
    int burst_span =
        (static_cast<int>(params.burst_size) - 1) * ripple_stagger_eff_ +
        ripple_dur_eff_;

    // Recipe entries insert a build phase on the segue's phase-1 plateau:
    // duration is lengthened by the build span rather than the carousel
    // growing an asymmetric-window API (docs/opchain_morph_spec.md,
    // section 6.1).
    int build_span = 0;
    if (recipe) {
      build_step_ = 0;
      build_total_frames_ = 0;
      for (size_t k = 0; k < build_step_count_; ++k) {
        const Solids::Op op = build_steps_[k].op;
        const int frames = std::max(1, static_cast<int>(leg_frames(op) / sp));
        // A gated leg runs both half-gates plus the swap frame, so its budget
        // rounds to the odd length the leg actually takes.
        build_leg_frames_[k] =
            is_gated_step(op) ? 2 * gate_frames(frames) + 1 : frames;
        build_total_frames_ += build_leg_frames_[k];
      }
      build_span = build_total_frames_;
      // One pinned target set per shape: colour converges across the whole
      // chain, so every leg blends toward the same assignment.
      for (int i = 0; i < NUM_PALETTES; ++i)
        build_targets_[i] = static_cast<uint8_t>(i);
      hs::shuffle(build_targets_.begin(), build_targets_.end());
    }

    int duration = fade + build_span + still + burst_span + still + fade;

    int next_delay = carousel.schedule_segue(timeline, draw_fn, duration, fade);

    // Added after the sprite: on the frame this fires the sprite has already
    // drawn the seed at the phase-1 boundary, and the first leg's first draw
    // lands on the next frame — no gap, no double draw.
    if (recipe) {
      timeline.add(fade, Animation::PeriodicTimer(
                             0,
                             [this](Canvas &) {
                               build_active_ = true;
                               start_build_leg();
                             },
                             false));
    }

    timeline.add(fade + build_span + still,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &canvas) { ripple(canvas); }, false));

    // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
    // A recipe shape spawns holding its seed, so these are the seed's counts;
    // finish_build logs the finished solid's, which are what it rasterizes.
    const auto &entry = solids[solid_idx];
    const MeshState &spawned = carousel.current();
    hs::log("Spawning Shape: %s (V=%d, E=%d, F=%d, I=%d)%s", entry.name,
            (int)spawned.vertices.size(), (int)(spawned.faces.size() / 2),
            (int)spawned.face_counts.size(), (int)spawned.faces.size(),
            recipe ? " seed" : "");

    // The segue decides when the next shape starts relative to this one.
    timeline.add(next_delay,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &) { this->spawn_shape(); }, false));
  }

  /**
   * @brief Constructs and schedules build leg build_step_: eagerly builds the
   *        clean endpoint seed and its bookend classification, derives the
   *        palette handoff, and chains completion into finish_build_leg.
   */
  HS_COLD_MEMBER void start_build_leg() {
    const size_t k = build_step_;
    const Solids::OpStep &step = build_steps_[k];

    // Eager clean endpoint seed_{k+1}: the mesh the leg lands on and the next
    // leg sweeps from. Runs first — generate() resets the scratch arenas the
    // handoff arrays below live in. A hankin leg builds none: its arrival is
    // the mesh its baked topology already carries, and its bookend grouping is
    // that arrival's own classification, which the leg computes itself.
    //
    // Each leg carries its color the whole way to the palette its own arrival
    // mesh calls for, so every intermediate solid stands fully colored rather
    // than caught mid-transition. Continuity across the boundary comes from
    // the next leg departing from this leg's target (see prev_pal below): the
    // leg lands at weight 1 and its successor opens at weight 0 on the same
    // palette. The last leg's arrival is the finished solid, so its targets
    // are the final ones and the closing swap stays exact.
    Animation::OpLeg::BookendClasses bookend;
    if (step.op != Solids::Op::HANKIN) {
      generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        build_next_seed_ =
            Solids::finalize_solid(clean_endpoint(step, a, b), target);
      });
      {
        ScratchScope a_guard(scratch_arena_a);
        ScratchScope b_guard(scratch_arena_b);
        MeshOps::classify_faces_by_topology(build_next_seed_, scratch_arena_a,
                                            scratch_arena_b, persistent_arena);
      }
      const size_t bookend_faces = build_next_seed_.face_counts.size();
      HS_CHECK(bookend_faces <= MAX_BUILD_FACES);
      bookend = {build_next_seed_.topology.data(), bookend_faces};
    }

    // Handoff arrays are ctor-scoped: scratch-backed under this scope, alive
    // through the OpLeg constructor's own LIFO-stacked scratch scopes.
    ScratchScope handoff_guard(scratch_arena_a);
    const size_t prev_faces = build_seed_.face_counts.size();
    HS_CHECK(prev_faces <= MAX_BUILD_FACES);
    Vector *prev_centroid = scratch_arena_a.allocate_n<Vector>(prev_faces);
    size_t off = 0;
    for (size_t f = 0; f < prev_faces; ++f) {
      Vector c(0.0f, 0.0f, 0.0f);
      const int n = build_seed_.face_counts[f];
      for (int j = 0; j < n; ++j)
        c = c + build_seed_.vertices[build_seed_.faces[off + j]];
      prev_centroid[f] = c.normalized();
      off += n;
    }
    const uint8_t *prev_pal;
    if (k == 0) {
      // The seed slot's displayed palettes are the chain's FROM state.
      const MeshState &seed_slot = carousel.current();
      HS_CHECK(prev_faces <= seed_slot.topology.size());
      const auto &slots = palettes_slots[carousel.front_index()];
      uint8_t *pal = scratch_arena_a.allocate_n<uint8_t>(prev_faces);
      for (size_t f = 0; f < prev_faces; ++f)
        pal[f] = static_cast<uint8_t>(
            slots[wrap(seed_slot.topology[f], NUM_PALETTES)]);
      prev_pal = pal;
    } else {
      // Depart from the palette the previous leg landed on, so the boundary is
      // continuous: that leg finished at weight 1 on exactly these ids.
      HS_CHECK(build_from_pal_ && build_from_faces_ == prev_faces,
               "IslamicStars: carried palette does not cover the leg seed");
      prev_pal = build_from_pal_;
    }

    Animation::OpLeg::PaletteHandoff handoff{
        &palette_bank_.bank, prev_pal,      nullptr, prev_faces,
        false,               prev_centroid, &build_targets_};

    const int frames = build_leg_frames_[k];
    hs::log("Build leg: %s (%d frames)", leg_name(step.op), frames);

    auto schedule = [this](Animation::OpLeg &&leg) {
      build_landing_ = &leg.landing();
      timeline.add(0, std::move(leg).then([this] { finish_build_leg(); }));
    };

    // The swept operator's endpoints: every inflate leg opens at the clamped
    // zero-area birth limit and lands on the step's own parameter; ambo is a
    // truncate swept to the short-circuit point.
    switch (step.op) {
    case Solids::Op::HANKIN:
      schedule(Animation::OpLeg(build_seed_, 0.0f, step.param,
                                persistent_arena, draw_build_fn_, handoff,
                                frames, bookend));
      break;
    case Solids::Op::RELAX:
      schedule(Animation::OpLeg(build_seed_, static_cast<int>(step.param),
                                persistent_arena, draw_build_fn_, handoff,
                                frames, bookend));
      break;
    case Solids::Op::AMBO:
      schedule(Animation::OpLeg(build_seed_, ConwayGraph::MorphOp::TRUNCATE,
                                0.0f, 0.5f, 0.0f, 0.0f, persistent_arena,
                                draw_build_fn_, handoff, frames, bookend));
      break;
    case Solids::Op::TRUNCATE:
      schedule(Animation::OpLeg(build_seed_, ConwayGraph::MorphOp::TRUNCATE,
                                0.0f, step.param, 0.0f, 0.0f, persistent_arena,
                                draw_build_fn_, handoff, frames, bookend));
      break;
    case Solids::Op::SNUB:
      schedule(Animation::OpLeg(build_seed_, ConwayGraph::MorphOp::SNUB, 0.0f,
                                step.param, 0.0f, step.twist, persistent_arena,
                                draw_build_fn_, handoff, frames, bookend));
      break;
    case Solids::Op::CHAMFER:
      schedule(Animation::OpLeg(build_seed_, ConwayGraph::MorphOp::CHAMFER,
                                0.0f, step.param, 0.0f, 0.0f, persistent_arena,
                                draw_build_fn_, handoff, frames, bookend));
      break;
    case Solids::Op::KIS:
      schedule(Animation::OpLeg(build_seed_, Animation::OpLeg::SwapOp::KIS,
                                persistent_arena, draw_build_fn_, handoff,
                                gate_frames(frames), bookend));
      break;
    case Solids::Op::DUAL:
      schedule(Animation::OpLeg(build_seed_, Animation::OpLeg::SwapOp::DUAL,
                                persistent_arena, draw_build_fn_, handoff,
                                gate_frames(frames), bookend));
      break;
    default:
      HS_CHECK(false, "IslamicStars: unsweepable primitive op reached a leg");
      break;
    }
  }

  /**
   * @brief Builds the clean endpoint mesh a leg lands on.
   * @param step Lowered primitive step; a hankin step has no entry, since its
   * leg rebuilds its own arrival from the topology it swept.
   * @param a Output arena for even pipeline stages.
   * @param b Scratch arena for odd pipeline stages.
   * @return The op applied to the current build seed at its exact parameter.
   */
  HS_COLD_MEMBER PolyMesh clean_endpoint(const Solids::OpStep &step, Arena &a,
                                         Arena &b) {
    switch (step.op) {
    case Solids::Op::AMBO:
      return MeshOps::ambo(build_seed_, a, b);
    case Solids::Op::TRUNCATE:
      return MeshOps::truncate(build_seed_, a, b, step.param);
    case Solids::Op::SNUB:
      return MeshOps::snub(build_seed_, a, b, step.param, step.twist);
    case Solids::Op::CHAMFER:
      return MeshOps::chamfer(build_seed_, a, b, step.param);
    case Solids::Op::RELAX:
      return MeshOps::relax(build_seed_, a, b, static_cast<int>(step.param));
    case Solids::Op::KIS:
      return MeshOps::kis(build_seed_, a, b);
    case Solids::Op::DUAL:
      return MeshOps::dual(build_seed_, a, b);
    default:
      HS_CHECK(false, "IslamicStars: step builds no eager endpoint");
      return PolyMesh{};
    }
  }

  /**
   * @brief Leg completion: adopts the eagerly built clean endpoint as the next
   *        leg's seed, then starts the next leg or finishes the build.
   */
  HS_COLD_MEMBER void finish_build_leg() {
    // Reclaim the finished leg. Only the endpoint the next leg sweeps from
    // crosses the reset, so a build holds one leg's storage rather than the
    // whole chain's. A hankin leg's endpoint is rebuilt here from the topology
    // it swept, into the scratch the evacuation below reads; the other kinds
    // carry the endpoint start_build_leg built eagerly. The palette the leg
    // landed on is snapshotted too: the next leg departs from it and the
    // landing does not survive.
    ScratchScope a_guard(scratch_arena_a);
    if (build_steps_[build_step_].op == Solids::Op::HANKIN) {
      build_next_seed_ = PolyMesh();
      Animation::OpLeg::arrival_mesh(*build_landing_, build_next_seed_,
                                     scratch_arena_a);
    }

    if (build_step_ + 1 >= build_step_count_) {
      // finish_build consumes the last leg's landing, so that one is reclaimed
      // by the closing compaction instead.
      build_seed_ = std::move(build_next_seed_);
      ++build_step_;
      finish_build();
      return;
    }

    const size_t landed_faces = build_landing_->faces;
    uint8_t *landed = scratch_arena_a.allocate_n<uint8_t>(landed_faces);
    for (size_t f = 0; f < landed_faces; ++f)
      landed[f] = build_landing_->to_palette[wrap(build_landing_->topology[f],
                                                  NUM_PALETTES)];
    build_landing_ = nullptr;

    {
      Persist<PolyMesh> pn(build_next_seed_, scratch_arena_b, persistent_arena);
      build_seed_ = PolyMesh();
      carousel.compact_drop_all([this](Arena &arena) {
        ripple_gen.reclaim_storage(arena);
        palette_bank_.bake_all(arena);
      });
    }
    build_seed_ = std::move(build_next_seed_);

    uint8_t *from_pal = persistent_arena.allocate_n<uint8_t>(landed_faces);
    std::memcpy(from_pal, landed, landed_faces);
    build_from_pal_ = from_pal;
    build_from_faces_ = landed_faces;

    ++build_step_;
    start_build_leg();
  }

  /**
   * @brief Build completion: recompiles the finished solid into the front slot
   *        and maps its palette slots from the last leg's landing, so the
   *        sprite's next frame is pixel-equal to the leg's last one.
   */
  HS_COLD_MEMBER void finish_build() {
    // The closing compile lands on a compacted arena, so everything it is
    // checked against is snapshotted out of the last leg's storage first.
    ScratchScope landing_guard(scratch_arena_a);
    const Animation::OpLeg::Landing &landing = *build_landing_;
    const size_t landed_faces = landing.faces;
    int *landed_topo = scratch_arena_a.allocate_n<int>(landed_faces);
    std::memcpy(landed_topo, landing.topology, landed_faces * sizeof(int));
    const std::array<uint8_t, Animation::OpLeg::PALETTES> landed_to =
        landing.to_palette;
    build_landing_ = nullptr;
    build_from_pal_ = nullptr;
    build_from_faces_ = 0;

    const int front = carousel.front_index();
    MeshState &slot = carousel.slot(front);
    // Not generate(): its depth-0 reset would drop the landing snapshot above.
    {
      ScratchScope a_guard(scratch_arena_a);
      slot.clear();
      // The seed the closing compile reads is evacuated to scratch and never
      // restored: the compiled slot is the only form the effect renders, so
      // letting the PolyMesh cross the compaction would carry a second copy of
      // the built solid for the rest of the shape's life. It is dropped again
      // before the classification, which needs scratch_arena_b in full.
      {
        ScratchScope seed_guard(scratch_arena_b);
        PolyMesh built;
        MeshOps::clone(build_seed_, built, scratch_arena_b);
        build_next_seed_ = PolyMesh();
        build_seed_ = PolyMesh();
        carousel.compact_drop_all([this](Arena &arena) {
          ripple_gen.reclaim_storage(arena);
          palette_bank_.bake_all(arena);
        });
        MeshOps::compile(built, slot, persistent_arena, scratch_arena_a);
      }
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(slot, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // The leg's targets keyed on this same classification (computed at leg
    // start from the same endpoint mesh); drift here would pop the swap.
    HS_CHECK(landed_faces == slot.topology.size(),
             "IslamicStars: built mesh face count differs from the last leg");
    MeshPaletteBank::shuffle_indices(palettes_slots[front]);
    bool slot_mapped[NUM_PALETTES] = {};
    for (size_t f = 0; f < landed_faces; ++f) {
      HS_CHECK(landed_topo[f] == slot.topology[f],
               "IslamicStars: arrival classification drifted across the build");
      const int s = wrap(slot.topology[f], NUM_PALETTES);
      if (!slot_mapped[s]) {
        slot_mapped[s] = true;
        palettes_slots[front][s] = landed_to[wrap(landed_topo[f], NUM_PALETTES)];
      }
    }
    build_active_ = false;

    hs::log("Built Shape: %s (V=%d, E=%d, F=%d, I=%d)",
            Solids::Collections::get_islamic_solids()[solid_idx].name,
            (int)slot.vertices.size(), (int)(slot.faces.size() / 2),
            (int)slot.face_counts.size(), (int)slot.faces.size());
  }

  /**
   * @brief Slider-backed runtime parameters for the effect.
   */
  struct Params {
    float fade =
        72.0f; /**< Segue window length, in frames: a 64-frame (4 s) sweep crossing plus one per-face fade tail. */
    float burst_size =
        4.0f; /**< Ripples per burst; float-backed for register_param. */
    float trans_speed =
        1.0f; /**< Divides every per-shape stage length (fade, still holds, ripple span): 1 = shipping cadence, higher cycles shapes faster. */
    bool debug_bb = false; /**< Whether to draw mesh bounding boxes. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
