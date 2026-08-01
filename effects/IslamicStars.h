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
// speed, build_active, solid_idx) for the effects-module build smoke.
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
#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
  static_assert(
      HS_ISLAMICSTARS_PROFILE_SHAPE >= 0 &&
          HS_ISLAMICSTARS_PROFILE_SHAPE <
              static_cast<int>(std::size(Solids::islamic_registry)),
      "HS_ISLAMICSTARS_PROFILE_SHAPE is outside the Islamic solid registry");
#endif

  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  HS_COLD_MEMBER IslamicStars()
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(filters)::any_crosses_segments,
                .reads_outside_band = decltype(filters)::any_reads_outside_band,
                .margin = decltype(filters)::max_segment_margin}),
        filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  HS_COLD_MEMBER void init() override {
    // Asymmetric scratch split (190 KB total): the leg-by-leg build chain
    // peaks at ~114 KB in a and ~69 KB in b, and compact_keep_front evacuates
    // the front slot (up to 63.7 KB) through b. The remainder is persistent:
    // carousel slots + BakedPaletteBank (~18 KB). Budgets enforced by the
    // test_conway_morph.h build replay and test_solids.h's high-water sweeps.
    configure_arenas(GLOBAL_ARENA_SIZE - SPLIT_SCRATCH_A_DEFAULT -
                         SPLIT_SCRATCH_B_DEFAULT,
                     SPLIT_SCRATCH_A_DEFAULT, SPLIT_SCRATCH_B_DEFAULT);

    ripple_gen.init_storage(persistent_arena);
    claim_face_palettes(persistent_arena);
    palette_bank.bake_all(persistent_arena);

    // Set BEFORE registering: register_param snaps *ptr as the slider default.
    // Amplitude starts at the fold-free ceiling (see RIPPLE_AMP_MAX).
    ripple_gen.template_params.amplitude = RIPPLE_AMP_MAX;
    ripple_gen.template_params.thickness = RIPPLE_THICKNESS;
    ripple_gen.template_params.decay = 0.1f;
#ifdef HS_PROFILE_TRANS_SPEED
    params.trans_speed = static_cast<float>(HS_PROFILE_TRANS_SPEED);
#endif

    // Per-face fade length range (frames): each face draws a random fade from
    // [lo, hi] as the terminator reaches it, fraying the sweep front.
    register_param("Face Fade Lo", &carousel.segue().fade_frames_min, 0.0f,
                   32.0f);
    register_param("Face Fade Hi", &carousel.segue().fade_frames_max, 0.0f,
                   32.0f);
    register_param("Burst", &params.burst_size, 1.0f, (float)BURST_MAX);
    // Amplitude slider capped at the fold-free ceiling; thickness is fixed (not a
    // slider) so amplitude/thickness can never cross the self-fold onset.
    register_param("Ripp Amp", &ripple_gen.template_params.amplitude, 0.0f,
                   RIPPLE_AMP_MAX);
    register_param("Ripp Decay", &ripple_gen.template_params.decay, 0.0f, 5.0f);
    register_param("Ripp Dur", &params.ripple_duration, 30.0f,
                   (float)RIPPLE_DURATION_MAX);
    register_param("Trans Speed", &params.trans_speed, 1.0f, 8.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, UP, noise));

#ifndef HS_ISLAMICSTARS_PROFILE_SHAPE
    // Open on a recipe entry so the op-by-op build is the first thing drawn;
    // spawn_shape pre-increments, so seed the index one before it.
    auto solids = Solids::Collections::get_islamic_solids();
    for (size_t i = 0; i < solids.size(); ++i) {
      if (solids[i].recipe) {
        solid_idx = static_cast<int>(i) - 1;
        break;
      }
    }
#endif

    spawn_shape();
  }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  void draw_frame() override {
    Canvas canvas(*this);
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
  static constexpr int SPRITE_FADE_FRAMES = 16;
  static constexpr int STILL_FRAMES =
      16; /**< 1 s hold (16 fps) between fade and ripple stages. */
  // Recipe-build leg lengths (docs/opchain_morph_spec.md section 7); divided
  // by the Trans Speed divisor like every other stage.
  static constexpr int HANKIN_LEG_FRAMES = 32;
  static constexpr int SWEEP_LEG_FRAMES =
      24; /**< ambo / truncate / snub / chamfer. */
  static constexpr int RELAX_LEG_FRAMES = 16;
  static constexpr int RECONCILE_LEG_FRAMES =
      24; /**< identity-mesh -> authored kis/needle slerp. */
  /** Identity-mesh truncate depth of the smooth kis/needle path: the "uniform"
   * Conway depth at which dual(truncate(X)) matches kis(dual(X)) exactly on
   * regular seeds (docs/opchain_morph_spec.md, smooth kis/needle). */
  static constexpr float MACRO_TRUNCATE_T = 1.0f / 3.0f;
  /** Per-shape arena split (bytes). A smooth kis/needle macro bridges
   * truncate(seed) (~3x the seed's half-edges); rendering and classifying that
   * tripled mesh needs a scratch_a heavier than the default, while the bridge's
   * per-leg compaction keeps its persistent well under the default. So a bridge
   * shape spawns on a scratch_a-heavy split; other recipe builds trade 2 KB of
   * scratch_b for persistent; generated whole shapes retain the wider scratch_b.
   * Each split's persistent is the device-arena remainder and both scratch
   * arenas are hard-capped at their device sizes. Needle's measured peak is
   * 131,770 / 70,228 / 96,600 bytes. */
  static constexpr size_t SPLIT_SCRATCH_A_DEFAULT = 116 * 1024;      // 118,784
  static constexpr size_t SPLIT_SCRATCH_A_BRIDGE = 129 * 1024 + 512; // 132,608
  static constexpr size_t SPLIT_SCRATCH_B_DEFAULT = 74 * 1024;       // 75,776
  static constexpr size_t SPLIT_SCRATCH_B_BUILD = 72 * 1024;         // 73,728
  static constexpr size_t SPLIT_SCRATCH_B_BRIDGE = 74 * 1024;        // 75,776
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
  // Effective per-shape stage lengths after the Trans Speed divisor, written by
  // spawn_shape and read by the deferred ripple() callback.
  int ripple_dur_eff = 80;
  int ripple_stagger_eff = RIPPLE_STAGGER_FRAMES;
  int solid_idx = -1;
  using SegueT = Segue::TerminatorSweep;
  struct SpriteFaceShading {
    const uint8_t *palette;

    SpriteFaceShading(const MeshState &mesh, const uint8_t *face_palette)
        : palette(face_palette) {
      HS_CHECK(mesh.get_topology_size() == mesh.num_faces(),
               "IslamicStars: sprite shading face count mismatch");
    }
  };

  MeshCarousel<SegueT> carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank;
  /** Per-slot per-face palette ids (persistent-arena backed, MAX_BUILD_FACES
   * each). Written at spawn (class-keyed colours) and at finish_build (the
   * last leg's landed colours); every compaction re-claims the same addresses,
   * so the contents survive the reset. */
  uint8_t *slot_face_palette[2] = {};

  // Build-chain state (entries with a non-null recipe): the shape is built op
  // by op between the fade-in and the still hold. Null-recipe entries never
  // touch any of it.
  bool build_active = false; /**< Legs draw; the sprite draw_fn is muted. */
  Solids::OpStep
      build_step_chain[MAX_BUILD_STEPS];      /**< Lowered primitive chain. */
  size_t build_step_count = 0;                /**< Lowered step count. */
  size_t build_step = 0;                      /**< Current leg index. */
  int build_leg_frames[MAX_BUILD_STEPS] = {}; /**< Per-leg frame budget. */
  int build_total_frames = 0;                 /**< Sum of leg frames. */
  PolyMesh build_seed;                        /**< Leg-k seed (persistent). */
  PolyMesh build_next_seed;  /**< Clean endpoint seed_{k+1}: built eagerly at
                                 leg start, or from the leg's own topology at
                                 its end on a hankin step. */
  PolyMesh dual_bridge_ambo; /**< ambo(P) held across a DUAL bridge: leg 1's
                                 arrival grouping and leg 2's departed mesh. */
  size_t dual_bridge_ambo_faces =
      0; /**< ambo(P) face count, kept for leg 3's handoff length after the mesh
            itself is dropped at the medial leg (persistent-budget relief). */
  /** Device persistent budget of the current shape's split, set by spawn_shape
   * before any read; the host arena is over-provisioned, so gates check the
   * resident persistent high-water against this. */
  size_t device_persistent_budget = 0;
  const Animation::OpLeg::Landing *build_landing =
      nullptr; /**< Latest leg's arrival data (leg-arena backed). */
  const uint8_t *build_from_pal =
      nullptr; /**< Per-face palette the previous leg landed on; survives the
                  leg-boundary compaction that drops its landing. */
  size_t build_from_faces = 0; /**< Length of build_from_pal. */
  int dual_bridges_built = 0;  /**< DUAL bridges scheduled (test coverage). */
  int build_macro_sweep_frames =
      SWEEP_LEG_FRAMES; /**< Truncate leg of a smooth
                                                       kis/needle macro. */
  int build_reconcile_frames =
      RECONCILE_LEG_FRAMES; /**< Reconcile leg length. */
  /** Continuation the smooth dual bridge chains after its closing leg: plain
   * finish_build_leg for a lone DUAL, or a macro's next stage. Set before every
   * schedule_dual_bridge call. */
  Fn<void(), 16> dual_bridge_done{[this] { finish_build_leg(); }};

  /**
   * @brief Draw callback for build-leg frames.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, MeshState &, const Animation::OpLeg::Shading &), 8>
      draw_build_fn{
          [this](Canvas &c, MeshState &m, const Animation::OpLeg::Shading &sh) {
            draw_build_mesh(c, m, sh);
          }};

  /**
   * @brief Claims the two per-slot face-palette arrays from the arena.
   * @param arena Persistent arena the arrays live in.
   * @details Runs at init and inside every compaction rebake, directly after
   * the ripple pool's claim: the allocation order is fixed, so the arrays
   * re-land at their original addresses and their bytes survive the reset.
   */
  void claim_face_palettes(Arena &arena) {
    for (uint8_t *&pal : slot_face_palette)
      pal = arena.allocate_n<uint8_t>(MAX_BUILD_FACES);
  }

  /**
   * @brief Rebake callback for every persistent compaction: the ripple pool
   * first (its slots re-land at their init_storage addresses), then the
   * face-palette arrays, then the palette bank.
   * @param arena The freshly reset persistent arena.
   */
  HS_COLD_MEMBER void reclaim_persistent(Arena &arena) {
    ripple_gen.reclaim_storage(arena);
    claim_face_palettes(arena);
    palette_bank.bake_all(arena);
  }

  /**
   * @brief Spawns one burst of burst_size ripples from a random origin,
   *        staggered RIPPLE_STAGGER_FRAMES apart, each expanding over
   *        params.ripple_duration frames.
   * @param canvas Unused render target for the timer callback signature.
   */
  void ripple(Canvas &) {
    Vector origin = random_vector();
    for (int i = 0; i < (int)params.burst_size; i++) {
      if (!ripple_gen.spawn(i * ripple_stagger_eff, origin,
                            PI_F / ripple_dur_eff, ripple_dur_eff))
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
   * @brief Sprite draw callback: draws the held shape for one envelope frame.
   * @param canvas Render target.
   * @param phase Sprite envelope phase: rises over the incoming window, holds 1,
   *        falls over the outgoing window.
   * @param back Carousel slot the shape was spawned into.
   * @details Cold (flash): runs once per frame, so its own body stays out of
   * ITCM (phantasm sits at the granule edge); only draw_shape's per-pixel scan
   * is hot. During the build window an OpLeg draws instead (one mesh per frame).
   */
  HS_COLD_MEMBER void draw_sprite(Canvas &canvas, float phase, int back) {
    if (build_active)
      return;
    const MeshState &mesh = carousel.slot(back);
    const SpriteFaceShading shading(mesh, slot_face_palette[back]);
    draw_shape(canvas, phase, mesh, shading);
  }

  /**
   * @brief Orients, ripple-distorts, and segue-shapes base_state, then
   *        rasterizes it with a per-face palette lookup.
   * @param canvas Render target receiving the rasterized mesh.
   * @param phase Segue phase in [0, 1] from the sprite envelope: rises over
   *        the incoming window, holds 1, falls over the outgoing window.
   * @param base_state Undistorted source mesh to transform and draw; carries
   *        the per-face topology classes.
   * @param shading Per-face palette ids.
   * @note Draws on the exact SDF path, not the congruence-class LUT
   * (mesh_classes.h): ripple/segue deformation makes a canonical LUT mis-shade
   * or pop. The facility is for effects whose meshes hold still.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE void
  draw_shape(Canvas &canvas, float phase, const MeshState &base_state,
             const SpriteFaceShading &shading) {
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase))
      return;
    const uint8_t *face_palette = shading.palette;

    HS_PROFILE(is_draw_shape);
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state = transform_shape(base_state);
    // transform borrows the source's per-face classes into the transformed mesh.
    const uint16_t *face_classes = transformed_state.get_topology_data();

    // Per-face segues order faces by their center, recomputed per frame: from
    // world space by default (the front stays fixed in the room while the
    // mesh rotates through it), or from the untransformed mesh for segues
    // declaring LOCAL_SWEEP (the front rides the mesh). The third argument is
    // the face's palette-slot class, mapped exactly as the fragment shader
    // maps it; class-agnostic sweeps ignore it.
    // phase is fixed for the whole draw call, so the segue's per-face phase
    // resolves here rather than per fragment; so does the face's palette, whose
    // slot is the same class that face_offset already needs.
    ArenaVector<float> face_phases;
    ArenaVector<const BakedPalette *> face_palettes;
    {
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
        const Vector c = Animation::OpLeg::face_vertex_sum(
            sweep_state.vertices.data(), fidx, foff[f], fcnt[f]);
        const int cls = MeshPaletteBank::slot_of(face_classes[f]);
        float off =
            seg.face_offset(normalized_or(c, UP), static_cast<int>(f), cls);
        float fade = seg.face_fade_frac(static_cast<int>(f));
        face_phases.push_back(seg.face_phase(phase, off, fade));
        face_palettes.push_back(&palette_bank[face_palette[f]]);
      }
    }

    {
      HS_PROFILE(is_mesh_scan);
      FacePaletteShader fragment_shader;
      auto select_face = [&](size_t fi, float size) {
        HS_CHECK(fi < face_phases.size(),
                 "IslamicStars: sprite shading face mismatch");
        fragment_shader.set_palette(face_palettes[fi]);
        fragment_shader.alpha = seg.opacity(face_phases[fi]);
        fragment_shader.scale = size > math::TOLERANCE ? 1.0f / size : 0.0f;
      };
      Scan::Mesh::draw_specialized<W, H>(filters, canvas, transformed_state,
                                         fragment_shader, scratch_arena_a,
                                         nullptr, select_face);
    }
  }

  /**
   * @brief Draws one build-leg frame: the compiled swept mesh, shaded from its
   *        per-face palettes on the segue's phase-1 identity plateau.
   * @param canvas Render target receiving the rasterized mesh.
   * @param mesh Compiled swept mesh (scratch-backed, this frame only).
   * @param sh Per-face palette table from the OpLeg.
   */
  HS_O3_FN void draw_build_mesh(Canvas &canvas, MeshState &mesh,
                                const Animation::OpLeg::Shading &sh) {
    if (mesh.vertices.is_empty())
      return;
    // Own scope labels: sharing the sprite's would parent two draw paths under
    // one counter, and a build-only window then prints an empty subtree while a
    // mixed window prints the child above its own parent's total.
    HS_PROFILE(is_build_draw);
    // The frame-local mesh and its source fill scratch_a at the 1082-face peak.
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
      MeshOps::transform_in_place(mesh, ripple_gen, camera);
    }
    const SegueT &seg = carousel.segue();
    FacePaletteShader fragment_shader;
    fragment_shader.alpha = seg.opacity(1.0f);

    auto select_face = [&](size_t fi, float size) {
      HS_CHECK(fi < sh.faces, "IslamicStars: build shading face mismatch");
      fragment_shader.set_palette(&sh.ramp_for(fi));
      fragment_shader.scale = size > math::TOLERANCE ? sh.gain / size : 0.0f;
    };

    {
      HS_PROFILE(is_build_scan);
      // Rasterize from scratch_b: the swept+compiled mesh fills scratch_a to
      // ~120.9 KB during a build leg, leaving no room for the scan's per-face
      // SDF::FaceScratchBuffer, while scratch_b is near-empty here (the Conway
      // op and compile temps have unwound). The sprite path scans from
      // scratch_a, where its transformed copy already lives.
      Scan::Mesh::draw_specialized<W, H>(filters, canvas, mesh, fragment_shader,
                                         scratch_arena_b, nullptr, select_face);
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
    case Solids::Op::DUAL:
      // The smooth dual is a three-leg bridge (truncate to ambo, medial slerp,
      // truncate down to the dual), each a normal single-mesh sweep.
      return 3 * SWEEP_LEG_FRAMES;
    default:
      return SWEEP_LEG_FRAMES;
    }
  }

  /**
   * @brief Whether a lowered DUAL step pairs with a trailing KIS (needle = kd =
   *        dt): the pair builds as the smooth dt macro over both steps.
   * @param k Lowered step index.
   */
  bool dt_pair_at(size_t k) const {
    return k + 1 < build_step_count &&
           build_step_chain[k].op == Solids::Op::DUAL &&
           build_step_chain[k + 1].op == Solids::Op::KIS;
  }

  /**
   * @brief Whether a lowered KIS step stands alone (not the tail of a dt pair);
   *        such a kis builds as the dtd macro (kis = dtd).
   * @param k Lowered step index.
   */
  bool standalone_kis_at(size_t k) const {
    return build_step_chain[k].op == Solids::Op::KIS &&
           !(k > 0 && build_step_chain[k - 1].op == Solids::Op::DUAL);
  }

  /**
   * @brief Whether the lowered chain runs a smooth kis/needle bridge (a dt pair
   *        or a standalone kis), which spawns on the scratch_a-heavy split.
   * @details Scanned at spawn to pick the per-shape arena split before the build
   * grows the arenas. Every such shape fits the bridge split (needle is the
   * largest).
   */
  bool build_uses_smooth_bridge() const {
    for (size_t k = 0; k < build_step_count; ++k)
      if (dt_pair_at(k) || standalone_kis_at(k))
        return true;
    return false;
  }

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
#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
    solid_idx = HS_ISLAMICSTARS_PROFILE_SHAPE;
#else
    solid_idx = (solid_idx + 1) % solids.size();
#endif
    int back = 1 - carousel.front_index();
    // The spawned mesh's shuffled palette order, consumed by class ordinal.
    std::array<int, NUM_PALETTES> palette_slots;
    MeshPaletteBank::shuffle_indices(palette_slots);

    int idx = solid_idx;
    // A recipe whose lowered chain contains a step no leg kind covers falls
    // back to today's whole-generate path, seed solid and all.
    const Solids::Recipe *recipe = solids[idx].recipe;
    if (recipe) {
      build_step_count = Solids::expand_to_primitives(*recipe, build_step_chain,
                                                      MAX_BUILD_STEPS);
      for (size_t k = 0; k < build_step_count; ++k) {
        if (!Solids::is_morphable_step(build_step_chain[k])) {
          hs::log("IslamicStars: %s has an unsweepable step, generating whole",
                  solids[idx].name);
          recipe = nullptr;
          build_step_count = 0;
          break;
        }
      }
    }

    auto draw_fn = [this, back](Canvas &canvas, float phase) {
      this->draw_sprite(canvas, phase, back);
    };

    // Compact the back slot, rebaking palettes into the fresh arena rather than
    // tracking them through the evacuation.
    //
    // A build regenerates both slots before either is drawn again — the seed
    // below, the finished solid at finish_build — so the outgoing shape is
    // dropped rather than evacuated. Keeping it costs the whole build the
    // previous shape's mesh, which at the roster's 1082-face entry is most of
    // the persistent arena.
    auto rebake = [this](Arena &arena) { reclaim_persistent(arena); };
    if (recipe)
      carousel.compact_drop_all(rebake);
    else
      carousel.compact_keep_front(rebake);

    // Per-shape arena re-split. Safe here: the compact above left persistent at
    // its ~baseline (carousel slots + palette bank, low addresses) and both
    // scratch arenas idle, so resplit_arenas moves only the boundaries -- the
    // long-lived content survives (it never resets the persistent offset). A
    // smooth kis/needle bridge shape gets the scratch_a-heavy split; other
    // recipes trade unused scratch_b for persistent; whole-generated shapes
    // keep the full generation scratch_b. Persistent takes the remainder.
    const bool bridge_split = recipe && build_uses_smooth_bridge();
    const size_t split_a =
        bridge_split ? SPLIT_SCRATCH_A_BRIDGE : SPLIT_SCRATCH_A_DEFAULT;
    const size_t split_b = bridge_split ? SPLIT_SCRATCH_B_BRIDGE
                                        : (recipe ? SPLIT_SCRATCH_B_BUILD
                                                  : SPLIT_SCRATCH_B_DEFAULT);
    device_persistent_budget = DEVICE_GLOBAL_ARENA_SIZE - split_a - split_b;
    resplit_arenas(GLOBAL_ARENA_SIZE - split_a - split_b, split_a, split_b);

    generate(persistent_arena, [&](Arena &target, Arena &a,
                                   Arena &b) HS_COLD_MEMBER {
      if (recipe) {
        // The build starts from the recipe's seed solid; the chain is swept
        // on screen leg by leg. The seed is also held as the first leg's
        // persistent PolyMesh.
        build_seed = Solids::finalize_solid(
            Solids::simple_registry[recipe->seed].generate(a, b), target);
        carousel.slot(back).clear();
        MeshOps::compile(build_seed, carousel.slot(back), target,
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

    // Colours of the spawned mesh (a build's seed, or the whole solid): the
    // shuffled palette order consumed by class ordinal — class ids are dense,
    // so class c is the c-th cohort.
    {
      const MeshState &spawned_mesh = carousel.slot(back);
      const size_t spawned_faces = spawned_mesh.topology.size();
      HS_CHECK(spawned_faces <= MAX_BUILD_FACES);
      MeshPaletteBank::assign_by_class(spawned_mesh.topology.data(),
                                       spawned_faces, palette_slots,
                                       slot_face_palette[back]);
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
    int fade = std::max(1, static_cast<int>(SPRITE_FADE_FRAMES / sp));
    int still = std::max(1, static_cast<int>(STILL_FRAMES / sp));
    ripple_dur_eff = std::max(8, static_cast<int>(params.ripple_duration / sp));
    ripple_stagger_eff =
        std::max(1, static_cast<int>(RIPPLE_STAGGER_FRAMES / sp));
    int burst_span =
        (static_cast<int>(params.burst_size) - 1) * ripple_stagger_eff +
        ripple_dur_eff;

    // Recipe entries insert a build phase on the segue's phase-1 plateau:
    // duration is lengthened by the build span rather than the carousel
    // growing an asymmetric-window API (docs/opchain_morph_spec.md,
    // section 6.1).
    int build_span = 0;
    if (recipe) {
      build_step = 0;
      build_from_pal = nullptr; // first leg departs the carousel seed slot
      build_from_faces = 0;
      build_total_frames = 0;
      build_macro_sweep_frames =
          std::max(1, static_cast<int>(SWEEP_LEG_FRAMES / sp));
      build_reconcile_frames =
          std::max(1, static_cast<int>(RECONCILE_LEG_FRAMES / sp));
      const int bridge_frames =
          std::max(1, static_cast<int>(leg_frames(Solids::Op::DUAL) / sp));
      // A smooth kis/needle macro spans more legs than its lowered step count:
      // a trailing dual,kis (dt) is truncate + dual bridge + reconcile; a
      // standalone kis (dtd) is dual bridge + truncate + dual bridge + reconcile.
      // build_leg_frames[k] carries the dual-bridge budget the bridge splits by
      // three; the truncate and reconcile legs draw fixed member budgets.
      for (size_t k = 0; k < build_step_count; ++k) {
        const Solids::Op op = build_step_chain[k].op;
        if (dt_pair_at(k)) {
          build_leg_frames[k] = bridge_frames; // dual bridge (DUAL step)
          build_total_frames +=
              build_macro_sweep_frames + bridge_frames + build_reconcile_frames;
          ++k; // the trailing kis is consumed by the dt macro
          continue;
        }
        if (standalone_kis_at(k)) {
          build_leg_frames[k] = bridge_frames; // both dtd bridges split this
          build_total_frames += bridge_frames + build_macro_sweep_frames +
                                bridge_frames + build_reconcile_frames;
          continue;
        }
        const int frames = std::max(1, static_cast<int>(leg_frames(op) / sp));
        build_leg_frames[k] = frames;
        build_total_frames += frames;
      }
      build_span = build_total_frames;
    }

    int duration = fade + build_span + still + burst_span + still + fade;

    int next_delay =
        carousel.schedule_segue(timeline, back, draw_fn, duration, fade);

    // Added after the sprite: on the frame this fires the sprite has already
    // drawn the seed at the phase-1 boundary, and the first leg's first draw
    // lands on the next frame — no gap, no double draw.
    if (recipe) {
      timeline.add(fade, Animation::PeriodicTimer(
                             0,
                             [this](Canvas &) {
                               build_active = true;
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
   * @brief Constructs and schedules build leg build_step: eagerly builds the
   *        clean endpoint seed and its bookend classification, derives the
   *        palette handoff, and chains completion into finish_build_leg.
   */
  HS_COLD_MEMBER void start_build_leg() {
    const size_t k = build_step;
    const Solids::OpStep &step = build_step_chain[k];

    // Smooth kis/needle macros (docs/opchain_morph_spec.md, smooth kis/needle):
    // a trailing dual,kis is the dt macro (spanning both steps), a standalone
    // kis is the dtd macro; each ends on a reconcile leg onto the exact authored
    // mesh. The recipe/expand_to_primitives still lower needle to {DUAL,KIS}.
    // Every dt/dtd chain takes the smooth path -- spawn_shape sizes the arena
    // split so even the needle's tripled bridge fits (see build_uses_smooth_bridge).
    if (dt_pair_at(k)) {
      schedule_dt_macro();
      return;
    }
    if (standalone_kis_at(k)) {
      schedule_dtd_macro();
      return;
    }

    // A lone DUAL (or a macro-ineligible dt pair's DUAL) is the smooth
    // three-leg bridge; it builds its own endpoints and chains its legs, then
    // rejoins at finish_build_leg.
    if (step.op == Solids::Op::DUAL) {
      dual_bridge_done = [this] { finish_build_leg(); };
      schedule_dual_bridge();
      return;
    }

    // Eager clean endpoint seed_{k+1}: the mesh the leg lands on and the next
    // leg sweeps from. Runs first — generate() resets the scratch arenas the
    // handoff arrays below live in. A hankin leg builds none: its arrival is
    // the mesh its baked topology already carries, and its bookend grouping is
    // that arrival's own classification, which the leg computes itself.
    //
    // Colours re-key per leg: the arrival's own classification maps to a
    // freshly shuffled palette set, and every face crossfades from the palette
    // the previous leg landed on to its new class target over the leg.
    Animation::OpLeg::BookendClasses bookend;
    if (step.op != Solids::Op::HANKIN) {
      generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        build_next_seed =
            Solids::finalize_solid(clean_endpoint(step, a, b), target);
      });
      {
        ScratchScope a_guard(scratch_arena_a);
        ScratchScope b_guard(scratch_arena_b);
        MeshOps::classify_faces_by_topology(build_next_seed, scratch_arena_a,
                                            scratch_arena_b, persistent_arena);
      }
      const size_t bookend_faces = build_next_seed.face_counts.size();
      HS_CHECK(bookend_faces <= MAX_BUILD_FACES);
      bookend = {.topology = build_next_seed.topology.data(),
                 .faces = bookend_faces};
    }

    // Handoff arrays are ctor-scoped: scratch-backed under this scope, alive
    // through the OpLeg constructor's own LIFO-stacked scratch scopes.
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);

    const int frames = build_leg_frames[k];
    hs::log("Build leg: %s (%d frames)", leg_name(step.op), frames);

    // The swept operator's endpoints: every inflate leg opens at the clamped
    // zero-area birth limit and lands on the step's own parameter; ambo is a
    // truncate swept to the short-circuit point.
    switch (step.op) {
    case Solids::Op::HANKIN:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::HankinSweepSpec{.theta_start = 0.0f,
                                            .theta_end = step.param,
                                            .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::RELAX:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::RelaxSpec{.iterations =
                                          static_cast<int>(step.param),
                                      .bake = step.bake,
                                      .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::AMBO:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                           .t_start = 0.0f,
                                           .t_end = 0.5f,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::TRUNCATE:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::SNUB:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::SNUB,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .twist_end = step.twist,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    case Solids::Op::CHAMFER:
      schedule_build_leg(Animation::OpLeg(
          build_seed,
          Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::CHAMFER,
                                           .t_start = 0.0f,
                                           .t_end = step.param,
                                           .sweep_frames = frames},
          persistent_arena, draw_build_fn, handoff, bookend));
      break;
    default:
      // Neither DUAL nor KIS reaches here: both route through the smooth
      // bridges above.
      HS_CHECK(false, "IslamicStars: unsweepable primitive op reached a leg");
      break;
    }
  }

  HS_COLD_MEMBER void schedule_build_leg(Animation::OpLeg &&leg) {
    build_landing = &leg.landing();
    timeline.add(0, std::move(leg).then([this] { finish_build_leg(); }));
  }

  /**
   * @brief Captures palette provenance for a leg departing the current seed.
   * @param scratch Arena the centroid and palette arrays live in; must outlive
   * the OpLeg constructor that reads them.
   * @param correspondence Departed-to-swept face-order relationship.
   */
  HS_COLD_MEMBER Animation::OpLeg::PaletteHandoff
  seed_handoff(Arena &scratch,
               Animation::OpLeg::FaceCorrespondence correspondence =
                   Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
    const size_t prev_faces = build_seed.face_counts.size();
    HS_CHECK(prev_faces <= MAX_BUILD_FACES);
    Vector *prev_centroid = nullptr;
    if (correspondence == Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
      prev_centroid = scratch.allocate_n<Vector>(prev_faces);
      Animation::OpLeg::face_centroids_into(build_seed, prev_centroid);
    }
    const uint8_t *prev_pal;
    if (!build_from_pal) {
      // The build's first leg departs the carousel seed slot: its per-face
      // spawn colours are the chain's FROM state. Keyed on build_from_pal
      // (reset to null per build) rather than build_step == 0, so a smooth
      // kis/needle macro whose later sub-legs still sit on step 0
      // (icosahedron_kis_gyro's leading kis) departs from the carried palette,
      // not the seed slot's.
      HS_CHECK(prev_faces <= carousel.current().topology.size());
      prev_pal = slot_face_palette[carousel.front_index()];
    } else {
      // Depart from the palette the previous leg landed on.
      HS_CHECK(build_from_pal && build_from_faces == prev_faces,
               "IslamicStars: carried palette does not cover the leg seed");
      prev_pal = build_from_pal;
    }
    return {.bank = &palette_bank.bank,
            .prev_face_palette = prev_pal,
            .prev_faces = prev_faces,
            .prev_face_centroid = prev_centroid,
            .correspondence = correspondence};
  }

  /**
   * @brief Captures palette provenance for a leg departing a prior landing.
   * @param departed Mesh the next leg departs from, in the previous leg's
   * landing face order (so its face f carries build_landing face f's palette).
   * @param scratch Arena the arrays live in.
   * @param correspondence Departed-to-swept face-order relationship.
   */
  HS_COLD_MEMBER Animation::OpLeg::PaletteHandoff
  landing_handoff(const PolyMesh &departed, Arena &scratch,
                  Animation::OpLeg::FaceCorrespondence correspondence =
                      Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
    const size_t nf = departed.face_counts.size();
    HS_CHECK(nf <= MAX_BUILD_FACES && build_landing &&
             build_landing->faces >= nf);
    Vector *cen = nullptr;
    uint8_t *pal = scratch.allocate_n<uint8_t>(nf);
    if (correspondence == Animation::OpLeg::FaceCorrespondence::GEOMETRIC) {
      cen = scratch.allocate_n<Vector>(nf);
      Animation::OpLeg::face_centroids_into(departed, cen);
    }
    for (size_t f = 0; f < nf; ++f)
      pal[f] = build_landing->landed_palette(f);
    return {.bank = &palette_bank.bank,
            .prev_face_palette = pal,
            .prev_faces = nf,
            .prev_face_centroid = cen,
            .correspondence = correspondence};
  }

  /**
   * @brief Frame count of DUAL bridge sub-leg `sub` (0/1/2), summing to the
   * step's budget.
   */
  int dual_sub_frames(int sub) const {
    const int total = build_leg_frames[build_step];
    const int third = std::max(1, total / 3);
    return sub < 2 ? third : std::max(1, total - 2 * third);
  }

  /**
   * @brief Schedules the smooth dual as three legs: truncate P -> ambo(P), a
   * medial slerp to ambo(dual(P)), and truncate dual(P) back down to dual(P).
   * @details Only ambo(P) (leg 1's arrival) is built up front; dual(P) is
   * deferred to leg 3 so it never co-resides with the medial leg's peak. Each
   * leg's scheduler compacts the arena before it runs, reclaiming the finished
   * legs and the endpoints they no longer need -- the heaviest gyro and
   * ambo_dual seeds run the whole bridge co-resident ~21 KB over budget.
   */
  HS_COLD_MEMBER void schedule_dual_bridge() {
    ++dual_bridges_built;
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      dual_bridge_ambo =
          Solids::finalize_solid(MeshOps::ambo(build_seed, a, b), target);
    });
    HS_CHECK(dual_bridge_ambo.face_counts.size() <= MAX_BUILD_FACES);

    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);
    const int frames = dual_sub_frames(0);
    hs::log("Build leg: dual bridge 1/3 truncate->ambo (%d frames)", frames);
    Animation::OpLeg leg(
        build_seed,
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.0f,
                                         .t_end = 0.5f,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true,
                                         .borrow_seed = true},
        persistent_arena, draw_build_fn, handoff);
    build_landing = &leg.landing();
    timeline.add(0, std::move(leg).then([this] { schedule_dual_medial(); }));
  }

  /**
   * @brief Schedules the dual bridge's medial-slerp leg (ambo(P) ->
   * ambo(dual(P))), departing from leg 1's landing. No bookend: the medial
   * arrival's own classification is the target leg 3 departs from.
   */
  HS_COLD_MEMBER void schedule_dual_medial() {
    ScratchScope handoff_guard(scratch_arena_a);
    // Snapshot the handoff off leg 1's landing (centroids from ambo(P), palette
    // from the leg-1 landing) before compacting away that finished leg. The
    // arrays live in scratch_a, which the persistent reset below leaves intact.
    Animation::OpLeg::PaletteHandoff handoff =
        landing_handoff(dual_bridge_ambo, scratch_arena_a,
                        Animation::OpLeg::FaceCorrespondence::IDENTITY);
    const size_t medial_faces = dual_bridge_ambo.face_counts.size();
    HS_CHECK(build_landing && build_landing->faces == medial_faces);
    uint16_t *medial_topology =
        scratch_arena_a.allocate_n<uint16_t>(medial_faces);
    std::copy_n(build_landing->topology, medial_faces, medial_topology);
    build_landing = nullptr;

    // Only ambo(P)'s face count survives to leg 3 (its handoff length); the mesh
    // is dead once the handoff centroids above are snapshotted, so drop it here
    // rather than carrying a full ambo copy through the medial leg's whole render
    // life -- the persistent spike that kept the needle over budget.
    dual_bridge_ambo_faces = dual_bridge_ambo.face_counts.size();
    dual_bridge_ambo = PolyMesh();

    // Compact: keep only the seed P (leg 2 rebuilds its medial from it, leg 3 its
    // dual), drop leg 1 and the ambo endpoint.
    {
      Persist<PolyMesh> ps(build_seed, scratch_arena_b, persistent_arena);
      carousel.compact_drop_all(
          [this](Arena &arena) { reclaim_persistent(arena); });
    }

    const int frames = dual_sub_frames(1);
    hs::log("Build leg: dual bridge 2/3 medial (%d frames)", frames);
    Animation::OpLeg leg(build_seed,
                         Animation::OpLeg::MedialSpec{.sweep_frames = frames},
                         persistent_arena, draw_build_fn, handoff,
                         Animation::OpLeg::BookendClasses{
                             .topology = medial_topology,
                             .faces = medial_faces});
    build_landing = &leg.landing();
    timeline.add(0,
                 std::move(leg).then([this] { schedule_dual_untruncate(); }));
  }

  /**
   * @brief Schedules the dual bridge's closing leg: truncate dual(P) from the
   * ambo point down to dual(P), landing on the macro's clean endpoint. The
   * departed centroids come from the medial leg's cached arrival, whose face
   * order matches leg 2's landing palettes.
   */
  HS_COLD_MEMBER void schedule_dual_untruncate() {
    ScratchScope a_guard(scratch_arena_a);
    const size_t nf = dual_bridge_ambo_faces;
    HS_CHECK(nf <= MAX_BUILD_FACES && build_landing &&
             build_landing->faces >= nf);
    uint8_t *pal = scratch_arena_a.allocate_n<uint8_t>(nf);
    for (size_t f = 0; f < nf; ++f)
      pal[f] = build_landing->landed_palette(f);
    build_landing = nullptr;

    // Compact: keep the seed P (leg 3 builds dual(P) from it), drop leg 2 (the
    // ambo(P) endpoint was already dropped at the medial leg).
    {
      Persist<PolyMesh> ps(build_seed, scratch_arena_b, persistent_arena);
      carousel.compact_drop_all(
          [this](Arena &arena) { reclaim_persistent(arena); });
    }

    // Build the deferred leg-3 seed dual(P) into the compacted arena. Not
    // generate(): its depth-0 scratch_a reset would drop the cen/pal snapshot
    // above, so scope the pipeline off the live frame instead.
    {
      ScratchScope da(scratch_arena_a);
      ScratchScope db(scratch_arena_b);
      build_next_seed = Solids::finalize_solid(
          MeshOps::dual(build_seed, scratch_arena_a, scratch_arena_b),
          persistent_arena);
    }
    {
      ScratchScope ca(scratch_arena_a);
      ScratchScope cb(scratch_arena_b);
      MeshOps::classify_faces_by_topology(build_next_seed, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    HS_CHECK(build_next_seed.face_counts.size() <= MAX_BUILD_FACES);

    Animation::OpLeg::PaletteHandoff handoff{
        .bank = &palette_bank.bank,
        .prev_face_palette = pal,
        .prev_faces = nf,
        .correspondence = Animation::OpLeg::FaceCorrespondence::DUAL_CLOSING};
    Animation::OpLeg::BookendClasses bookend{
        .topology = build_next_seed.topology.data(),
        .faces = build_next_seed.face_counts.size()};
    const int frames = dual_sub_frames(2);
    hs::log("Build leg: dual bridge 3/3 truncate->dual (%d frames)", frames);
    Animation::OpLeg leg(
        build_next_seed,
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.5f,
                                         .t_end = 0.0f,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true,
                                         .borrow_seed = true},
        persistent_arena, draw_build_fn, handoff, bookend);
    build_landing = &leg.landing();
    // Rejoin the caller's continuation: finish_build_leg for a lone DUAL, or the
    // next stage of a smooth kis/needle macro.
    timeline.add(0, std::move(leg).then([this] { dual_bridge_done(); }));
  }

  /**
   * @brief Adopts the eagerly built endpoint (build_next_seed) as the next
   * leg's seed, snapshots the palette the finished leg landed on so its
   * successor departs continuously, and compacts the finished leg's storage.
   * @details The shared middle of finish_build_leg's non-last path, reused by
   * the smooth kis/needle macro stages, which chain their own next leg rather
   * than advancing build_step.
   */
  HS_COLD_MEMBER void carry_landing_to_seed() {
    const size_t landed_faces = build_next_seed.face_counts.size();
    HS_CHECK(landed_faces <= build_landing->faces,
             "IslamicStars: next seed larger than the leg landing");
    HS_CHECK(landed_faces <= MAX_BUILD_FACES,
             "IslamicStars: next seed exceeds the slot palette capacity");
    // The carry lives in the idle slot's face-palette array: the outgoing
    // shape was dropped at the recipe spawn, and the array's same-address
    // re-claim keeps the bytes across every boundary compaction.
    uint8_t *carry = slot_face_palette[1 - carousel.front_index()];
    for (size_t f = 0; f < landed_faces; ++f)
      carry[f] = build_landing->landed_palette(f);
    build_landing = nullptr;

    {
      Persist<PolyMesh> pn(build_next_seed, scratch_arena_b, persistent_arena);
      build_seed = PolyMesh();
      carousel.compact_drop_all(
          [this](Arena &arena) { reclaim_persistent(arena); });
    }
    build_seed = std::move(build_next_seed);

    build_from_pal = carry;
    build_from_faces = landed_faces;
  }

  /**
   * @brief Schedules a plain truncate sweep of the current build seed to
   * MACRO_TRUNCATE_T, landing on build_next_seed = truncate(seed, 1/3).
   * @param log Log label ("dt truncate" / "dtd truncate").
   * @param then Completion chained after the leg.
   */
  template <typename Then>
  HS_COLD_MEMBER void schedule_macro_truncate(const char *log, Then &&then) {
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      build_next_seed = Solids::finalize_solid(
          MeshOps::truncate(build_seed, a, b, MACRO_TRUNCATE_T), target);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(build_next_seed, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    HS_CHECK(build_next_seed.face_counts.size() <= MAX_BUILD_FACES);
    Animation::OpLeg::BookendClasses bookend{
        .topology = build_next_seed.topology.data(),
        .faces = build_next_seed.face_counts.size()};
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(scratch_arena_a);
    const int frames = build_macro_sweep_frames;
    hs::log("Build leg: %s (%d frames)", log, frames);
    Animation::OpLeg leg(
        build_seed,
        Animation::OpLeg::ParamSweepSpec{.op = ConwayGraph::MorphOp::TRUNCATE,
                                         .t_start = 0.0f,
                                         .t_end = MACRO_TRUNCATE_T,
                                         .sweep_frames = frames,
                                         .bridge_provenance = true,
                                         .borrow_seed = true},
        persistent_arena, draw_build_fn, handoff, bookend);
    build_landing = &leg.landing();
    timeline.add(0, std::move(leg).then(std::forward<Then>(then)));
  }

  /**
   * @brief Schedules the smooth dt macro for a trailing dual,kis (needle = kd =
   * dt): truncate(X, 1/3) sweep, dual bridge on it, then reconcile onto the
   * exact kis(dual(X)) mesh. Covers both the DUAL step and its trailing KIS.
   */
  HS_COLD_MEMBER void schedule_dt_macro() {
    schedule_macro_truncate("dt truncate", [this] { dt_after_truncate(); });
  }
  HS_COLD_MEMBER void dt_after_truncate() {
    carry_landing_to_seed(); // build_seed = truncate(X, 1/3)
    dual_bridge_done = [this] { dt_after_bridge(); };
    schedule_dual_bridge();
  }
  HS_COLD_MEMBER void dt_after_bridge() {
    carry_landing_to_seed(); // build_seed = dual(truncate(X, 1/3)) (identity)
    ++build_step; // advance onto the KIS index the reconcile finishes at
    schedule_reconcile(build_step - 1, /*kis_of_dual=*/true);
  }

  /**
   * @brief Schedules the smooth dtd macro for a standalone kis (kis = dtd):
   * dual bridge on X, truncate(dual(X), 1/3) sweep, dual bridge on it, then
   * reconcile onto the exact kis(X) mesh. Runs entirely on the KIS step.
   */
  HS_COLD_MEMBER void schedule_dtd_macro() {
    dual_bridge_done = [this] { dtd_after_bridge1(); };
    schedule_dual_bridge();
  }
  HS_COLD_MEMBER void dtd_after_bridge1() {
    carry_landing_to_seed(); // build_seed = dual(X)
    schedule_macro_truncate("dtd truncate", [this] { dtd_after_truncate(); });
  }
  HS_COLD_MEMBER void dtd_after_truncate() {
    carry_landing_to_seed(); // build_seed = truncate(dual(X), 1/3)
    dual_bridge_done = [this] { dtd_after_bridge2(); };
    schedule_dual_bridge();
  }
  HS_COLD_MEMBER void dtd_after_bridge2() {
    carry_landing_to_seed(); // build_seed = dual(truncate(dual(X), 1/3))
    schedule_reconcile(build_step, /*kis_of_dual=*/false);
  }

  /**
   * @brief Relocates the authored kis/needle mesh's vertices onto the identity
   * mesh's connectivity through the nearest-vertex bijection.
   * @param identity Identity mesh (dt/dtd result): its connectivity and vertex
   * order are kept; each vertex is matched to the authored vertex nearest it.
   * @param authored Exact kis/needle mesh (same V as identity).
   * @param out Receives identity's connectivity carrying the matched authored
   * positions, in identity's vertex order.
   * @param target Arena backing @p out.
   * @param scratch Arena for the injectivity bookkeeping.
   * @details The correspondence is a proven bijection (identical topology, ~3.4%
   * position gap); a non-injective match traps rather than silently folding a
   * face.
   */
  HS_COLD_MEMBER void build_reconcile_endpoint(const PolyMesh &identity,
                                               const PolyMesh &authored,
                                               PolyMesh &out, Arena &target,
                                               Arena &scratch) {
    const size_t V = identity.vertices.size();
    HS_CHECK(authored.vertices.size() == V,
             "IslamicStars: reconcile endpoints differ in vertex count");
    ScratchScope guard(scratch);
    bool *used = scratch.allocate_n<bool>(V);
    std::fill_n(used, V, false);
    out.vertices.bind(target, V);
    for (size_t i = 0; i < V; ++i) {
      int best = -1;
      float best_dot = -2.0f;
      for (size_t j = 0; j < V; ++j) {
        const float d = dot(identity.vertices[i], authored.vertices[j]);
        if (d > best_dot) {
          best_dot = d;
          best = static_cast<int>(j);
        }
      }
      HS_CHECK(best >= 0 && !used[best],
               "IslamicStars: reconcile nearest-vertex map is not a bijection");
      used[best] = true;
      out.vertices.push_back(authored.vertices[best]);
    }
    out.face_counts.bind(target, identity.face_counts.size());
    out.face_counts.append_bulk(identity.face_counts.data(),
                                identity.face_counts.size());
    out.faces.bind(target, identity.faces.size());
    out.faces.append_bulk(identity.faces.data(), identity.faces.size());
  }

  /**
   * @brief Schedules the reconcile leg closing a smooth kis/needle macro: a
   * per-vertex great-circle slerp from the identity mesh (build_seed) onto the
   * exact authored kis/needle mesh, landing on the generator's mesh.
   * @param x_prefix Lowered-step count replayed to rebuild X, the mesh the macro
   * departed from (the generator's exact intermediate).
   * @param kis_of_dual True for a dt macro (authored = kis(dual(X)) = needle(X)),
   * false for a dtd macro (authored = kis(X)).
   */
  HS_COLD_MEMBER void schedule_reconcile(size_t x_prefix, bool kis_of_dual) {
    const uint8_t seed =
        Solids::Collections::get_islamic_solids()[solid_idx].recipe->seed;
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      PolyMesh X = Solids::build_steps(seed, build_step_chain, x_prefix,
                                       scratch_arena_a, scratch_arena_b);
      PolyMesh Xc;
      MeshOps::clone(X, Xc, scratch_arena_a);
      PolyMesh authored =
          kis_of_dual ? MeshOps::needle(Xc, scratch_arena_a, scratch_arena_b)
                      : MeshOps::kis(Xc, scratch_arena_a, scratch_arena_b);
      build_reconcile_endpoint(build_seed, authored, build_next_seed,
                               persistent_arena, scratch_arena_a);
    }
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(build_next_seed, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }
    HS_CHECK(build_next_seed.face_counts.size() <= MAX_BUILD_FACES);
    Animation::OpLeg::BookendClasses bookend{
        .topology = build_next_seed.topology.data(),
        .faces = build_next_seed.face_counts.size()};
    ScratchScope handoff_guard(scratch_arena_a);
    Animation::OpLeg::PaletteHandoff handoff = seed_handoff(
        scratch_arena_a, Animation::OpLeg::FaceCorrespondence::IDENTITY);
    const int frames = build_reconcile_frames;
    hs::log("Build leg: reconcile (%d frames)", frames);
    Animation::OpLeg leg(build_seed,
                         Animation::OpLeg::ReconcileSpec{
                             .to_positions = build_next_seed.vertices.data(),
                             .sweep_frames = frames},
                         persistent_arena, draw_build_fn, handoff, bookend);
    build_landing = &leg.landing();
    timeline.add(0, std::move(leg).then([this] { finish_build_leg(); }));
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
      return MeshOps::ambo(build_seed, a, b);
    case Solids::Op::TRUNCATE:
      return MeshOps::truncate(build_seed, a, b, step.param);
    case Solids::Op::SNUB:
      return MeshOps::snub(build_seed, a, b, step.param, step.twist);
    case Solids::Op::CHAMFER:
      return MeshOps::chamfer(build_seed, a, b, step.param);
    case Solids::Op::RELAX:
      if (step.bake)
        return MeshOps::relax_baked(build_seed, a, *step.bake);
      return MeshOps::relax(build_seed, a, b, static_cast<int>(step.param));
    case Solids::Op::KIS:
      return MeshOps::kis(build_seed, a, b);
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
    if (build_step_chain[build_step].op == Solids::Op::HANKIN) {
      build_next_seed = PolyMesh();
      Animation::OpLeg::arrival_mesh(*build_landing, build_next_seed,
                                     scratch_arena_a);
    }

    if (build_step + 1 >= build_step_count) {
      // finish_build consumes the last leg's landing, so that one is reclaimed
      // by the closing compaction instead.
      build_seed = std::move(build_next_seed);
      ++build_step;
      finish_build();
      return;
    }

    // Carry the emission-order prefix the next leg departs from (the whole
    // landing for a face-count-preserving leg, the survivor prefix where a leg's
    // arrival has more faces than its clean endpoint -- the DUAL bridge's
    // closing truncate lands V+F faces but hands off the V dual faces), then
    // advance to the next lowered step.
    carry_landing_to_seed();
    ++build_step;
    start_build_leg();
  }

  /**
   * @brief Build completion: recompiles the finished solid into the front slot
   *        and hands the last leg's per-face colours to the sprite, so its
   *        next frame is pixel-equal to the leg's last one.
   */
  HS_COLD_MEMBER void finish_build() {
    // The finished solid is build_seed (the last leg's clean endpoint): its
    // face count is the emission-order prefix of the landing (the whole landing
    // for a normal leg, the surviving dual faces for the DUAL bridge's closing
    // truncate, whose zero-area corner births drop at the compile below).
    const size_t landed_faces = build_seed.face_counts.size();
    HS_CHECK(landed_faces <= build_landing->faces,
             "IslamicStars: finished solid larger than the leg landing");
    HS_CHECK(landed_faces <= MAX_BUILD_FACES,
             "IslamicStars: finished solid exceeds the slot palette capacity");
    HS_CHECK(build_landing->topology,
             "IslamicStars: finished leg has no topology");
    ScratchScope topology_guard(scratch_arena_b);
    uint16_t *landed_topology =
        scratch_arena_b.allocate_n<uint16_t>(landed_faces);
    std::copy_n(build_landing->topology, landed_faces, landed_topology);
    const int front = carousel.front_index();
    // Per-face sprite handoff: copied before the compaction below, whose
    // same-address re-claim keeps the array's bytes.
    for (size_t f = 0; f < landed_faces; ++f)
      slot_face_palette[front][f] = build_landing->landed_palette(f);
    build_landing = nullptr;
    build_from_pal = nullptr;
    build_from_faces = 0;

    MeshState &slot = carousel.slot(front);
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
        MeshOps::clone(build_seed, built, scratch_arena_b);
        build_next_seed = PolyMesh();
        build_seed = PolyMesh();
        carousel.compact_drop_all(
            [this](Arena &arena) { reclaim_persistent(arena); });
        MeshOps::compile(built, slot, persistent_arena, scratch_arena_a);
      }
      slot.topology.bind(persistent_arena, landed_faces);
      for (size_t f = 0; f < landed_faces; ++f)
        slot.topology.push_back(landed_topology[f]);
    }

    // The per-face colours index the compiled slot by emission order, and
    // compile() strips faces with fewer than 3 sides, so the counts must agree
    // exactly.
    HS_CHECK(landed_faces == slot.num_faces(),
             "IslamicStars: compiled face count differs from the leg landing");
    build_active = false;

    hs::log("Built Shape: %s (V=%d, E=%d, F=%d, I=%d)",
            Solids::Collections::get_islamic_solids()[solid_idx].name,
            (int)slot.vertices.size(), (int)(slot.faces.size() / 2),
            (int)slot.face_counts.size(), (int)slot.faces.size());
  }

  /**
   * @brief Slider-backed runtime parameters for the effect.
   */
  struct Params {
    float burst_size =
        4.0f; /**< Ripples per burst; float-backed for register_param. */
    float ripple_duration =
        80.0f; /**< Frames each ripple takes to expand across the sphere. */
    float trans_speed =
        1.0f; /**< Divides every per-shape stage length (fade, still holds, ripple span): 1 = shipping cadence, higher cycles shapes faster. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(IslamicStars)
