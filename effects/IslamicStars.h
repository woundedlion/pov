/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file IslamicStars.h
 * @brief Sequence of Islamic-geometry polyhedra, built op by op and morphed
 *        one into the next.
 */

#include "core/animation/orientation.h"
#include "core/animation/animation.h"
#include "core/engine/engine.h"
#include "core/animation/recipe_build.h"

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
 *        choreography runs
 *        (docs/specs/opchain_morph_spec.md, "Recipe model").
 * @tparam W Target canvas width in pixels.
 * @tparam H Target canvas height in pixels.
 */
template <int W, int H>
class IslamicStars
    : public Effect,
      private Animation::RecipeBuild<IslamicStars<W, H>, 8, 1152> {

public:
#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
  static_assert(
      HS_ISLAMICSTARS_PROFILE_SHAPE >= 0 &&
          HS_ISLAMICSTARS_PROFILE_SHAPE <
              static_cast<int>(std::size(Solids::islamic_registry)),
      "HS_ISLAMICSTARS_PROFILE_SHAPE is outside the Islamic solid registry");
#endif

  struct ArenaBudget {
    size_t scratch_a;
    size_t scratch_b;

    constexpr size_t persistent(size_t total = DEVICE_GLOBAL_ARENA_SIZE) const {
      return total - scratch_a - scratch_b;
    }
  };

  static constexpr ArenaBudget GENERATED_BUDGET{116 * 1024, 74 * 1024};
  static constexpr ArenaBudget RECIPE_BUDGET{116 * 1024, 72 * 1024};
  static constexpr ArenaBudget BRIDGE_BUDGET{129 * 1024 + 512, 74 * 1024};
  static_assert(BRIDGE_BUDGET.scratch_a + BRIDGE_BUDGET.scratch_b <
                DEVICE_GLOBAL_ARENA_SIZE);

  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  HS_COLD_MEMBER IslamicStars()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  HS_COLD_MEMBER void init() override {
    configure_arenas(GENERATED_BUDGET.persistent(GLOBAL_ARENA_SIZE),
                     GENERATED_BUDGET.scratch_a, GENERATED_BUDGET.scratch_b);
    device_persistent_budget = GENERATED_BUDGET.persistent();

    ripple_gen.init_storage(persistent_arena);
    claim_face_palettes(persistent_arena);
    palette_bank.bake_all(persistent_arena);

    // Set BEFORE registering: register_param snaps *ptr as the slider default.
    // Amplitude starts at the slider ceiling (see RIPPLE_AMP_MAX).
    ripple_gen.template_params.amplitude = RIPPLE_AMP_MAX;
    ripple_gen.template_params.thickness = RIPPLE_THICKNESS;
    ripple_gen.template_params.decay = 0.1f;
#ifdef HS_PROFILE_TRANS_SPEED
    static_assert(HS_PROFILE_TRANS_SPEED >= 1 && HS_PROFILE_TRANS_SPEED <= 8,
                  "HS_PROFILE_TRANS_SPEED must be in [1, 8]");
    params.trans_speed = static_cast<float>(HS_PROFILE_TRANS_SPEED);
#endif

    // Per-face fade length range (frames): each face draws a random fade from
    // [lo, hi] as the terminator reaches it, fraying the sweep front.
    register_param("Face Fade Lo", &carousel.segue().fade_frames_min, 0.0f,
                   32.0f);
    register_param("Face Fade Hi", &carousel.segue().fade_frames_max, 0.0f,
                   32.0f);
    register_int_param("Burst", &params.burst_size, 1, BURST_MAX);
    // Amplitude slider capped at RIPPLE_AMP_MAX; thickness is fixed (not a
    // slider), so no setting exceeds the ratio RIPPLE_AMP_MAX is sized for.
    register_param("Ripp Amp", &ripple_gen.template_params.amplitude, 0.0f,
                   RIPPLE_AMP_MAX);
    register_param("Ripp Decay", &ripple_gen.template_params.decay, 0.0f, 5.0f);
    register_param("Ripp Dur", &params.ripple_duration, 30.0f,
                   (float)RIPPLE_DURATION_MAX);
    register_param("Trans Speed", &params.trans_speed, 1.0f, 8.0f);

    timeline.add(0, Animation::RandomWalk<W>(orientation, math::UP, noise));

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
  using Builder = Animation::RecipeBuild<IslamicStars<W, H>, 8, 1152>;
  friend Builder;
  using Builder::build_active;
  using Builder::build_entry;
  using Builder::build_step_chain;
  using Builder::build_step_count;
  using Builder::build_step;
  using Builder::build_seed;
  using Builder::build_total_frames;
  using Builder::device_persistent_budget;
  using Builder::dual_bridges_built;
  using typename Builder::BuildContinuation;
  using Builder::check_build_budget;
  using Builder::schedule_dual_bridge;
  using Builder::clean_endpoint;
  using Builder::build_uses_smooth_bridge;
  using Builder::plan_build_legs;
  using Builder::start_build_leg;
  static constexpr size_t MAX_BUILD_STEPS = 8;
  static constexpr size_t MAX_BUILD_FACES = 1152;
  static_assert(Solids::max_lowered_step_count(Solids::islamic_registry) <=
                MAX_BUILD_STEPS);
  friend struct ::hs_test::effects_tests::IslamicBuildProbe;

  // Ripple-pool sizing: a slot is held from spawn() until the staggered ripple
  // completes. Only one burst is normally live; the pool holds two so a Ripp
  // Dur/Burst slider change mid-burst cannot drop a spawn.
  static constexpr int RIPPLE_POOL_SIZE = 8;
  static constexpr int RIPPLE_STAGGER_FRAMES = 16;
  /** Ripp Dur slider ceiling. Nothing downstream bounds it: the per-shape
   * display window is derived from the burst span, so the burst always
   * completes before the next shape spawns. */
  static constexpr int RIPPLE_DURATION_MAX = 143;
  static constexpr int BURST_MAX = 4;
  static constexpr int SPRITE_FADE_FRAMES = 16;
  static constexpr int STILL_FRAMES =
      16; /**< 1 s hold (16 fps) between fade and ripple stages. */
  static constexpr float RIPPLE_THICKNESS =
      0.7f; /**< Fixed ripple wavelet width (radians). */
  /** Amplitude ceiling. Equals RIPPLE_SMALL_ANGLE_MAX, so every ripple
   * rotation takes the series-form quaternion. The displacement map
   * d -> d + theta(d) is injective only below amp/thickness = 0.181; at
   * 0.15/0.7 it folds across a <= 0.1 rad band by <= 0.012 rad (about half
   * a pixel at W=288), and at the default decay only within ~1.7 rad of the
   * origin. */
  static constexpr float RIPPLE_AMP_MAX = 0.15f;
  static_assert(2 * BURST_MAX <= RIPPLE_POOL_SIZE,
                "IslamicStars: ripple pool must hold two overlapping bursts");

  // orientation and noise are borrowed by the timeline-resident RandomWalk, so
  // they are declared before the Timeline to outlive it; ripple_gen must stay
  // after it, since a TransformerPool drops its clear hook through the
  // reference it holds.
  math::Orientation<> orientation;
  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H> filters;
  RippleTransformer<RIPPLE_POOL_SIZE> ripple_gen;
  // Effective per-shape stage lengths after the Trans Speed divisor, written by
  // spawn_shape and read by the deferred ripple() callback.
  int ripple_dur_eff = 80;
  int ripple_stagger_eff = RIPPLE_STAGGER_FRAMES;
  int burst_size_eff = 4;
  int solid_idx = -1;
  using SegueT = Segue::TerminatorSweep;

  MeshCarousel<SegueT> carousel;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank;
  /** Per-slot per-face palette ids (persistent-arena backed, MAX_BUILD_FACES
   * each). Written at spawn (class-keyed colours) and at finish_build (the
   * last leg's landed colours); every compaction re-claims the same addresses,
   * so the contents survive the reset. */
  uint8_t *slot_face_palette[2] = {};

  /**
   * @brief Claims the two per-slot face-palette arrays from the arena.
   * @param arena Persistent arena the arrays live in.
   * @details Runs at init and inside every compaction rebake, directly after
   * the ripple pool's claim: the allocation order is fixed, so the arrays
   * re-land at their original addresses (asserted) and their bytes survive the
   * reset.
   */
  void claim_face_palettes(Arena &arena) {
    for (uint8_t *&pal : slot_face_palette) {
      uint8_t *prev = pal;
      pal = arena.allocate_n<uint8_t>(MAX_BUILD_FACES);
      HS_CHECK(prev == nullptr || prev == pal,
               "IslamicStars: face-palette array moved across the compaction");
    }
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
   *        staggered ripple_stagger_eff frames apart, each expanding over
   *        ripple_dur_eff frames — the Trans-Speed-divided values spawn_entry
   *        caches, which bottom out at 2 and max(8, ripple_duration / 8).
   * @param canvas Unused render target for the timer callback signature.
   */
  void ripple(Canvas &) {
    math::Vector origin = math::random_vector();
    for (int i = 0; i < burst_size_eff; i++) {
      if (!ripple_gen.spawn(i * ripple_stagger_eff, origin,
                            math::PI_F / ripple_dur_eff, ripple_dur_eff))
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
    draw_shape(canvas, phase, mesh, slot_face_palette[back]);
  }

  /**
   * @brief Orients, ripple-distorts, and segue-shapes base_state, then
   *        rasterizes it with a per-face palette lookup.
   * @param canvas Render target receiving the rasterized mesh.
   * @param phase Segue phase in [0, 1] from the sprite envelope: rises over
   *        the incoming window, holds 1, falls over the outgoing window.
   * @param base_state Undistorted source mesh to transform and draw; carries
   *        the per-face topology classes.
   * @param face_palette Per-face palette ids.
   * @note Draws on the exact SDF path, not the congruence-class LUT
   * (face_class_bake.h): ripple/segue deformation makes a canonical LUT mis-shade
   * or pop. The facility is for effects whose meshes hold still.
   */
  HS_O3_FN HS_NOINLINE_NOCLONE void draw_shape(Canvas &canvas, float phase,
                                               const MeshState &base_state,
                                               const uint8_t *face_palette) {
    HS_CHECK(base_state.get_topology_size() == base_state.num_faces(),
             "IslamicStars: sprite shading face count mismatch");
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase))
      return;

    HS_PROFILE(is_draw_shape);
    ScratchScope a_guard(scratch_arena_a);
    MeshState transformed_state = transform_shape(base_state);
    // transform borrows the source's per-face classes into the transformed mesh.
    const uint16_t *face_classes = transformed_state.get_topology_data();

    // Per-face segues order faces by their center, recomputed per frame: from
    // world space by default, or from the untransformed mesh for segues
    // declaring LOCAL_SWEEP. The third argument is the face's palette-slot
    // class, mapped exactly as the fragment shader maps it; class-agnostic
    // sweeps ignore it.
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
        const math::Vector c = Animation::OpLeg::face_vertex_sum(
            sweep_state.vertices.data(), fidx, foff[f], fcnt[f]);
        const int cls = MeshPaletteBank::slot_of(face_classes[f]);
        float off = seg.face_offset(math::normalized_or(c, math::UP),
                                    static_cast<int>(f), cls);
        float fade = seg.face_fade_frac(static_cast<int>(f));
        face_phases.push_back(seg.face_phase(phase, off, fade));
        face_palettes.push_back(&palette_bank[face_palette[f]].view());
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
    // Opened after the OpLeg's blended ramps, which the shader below still
    // reads, so only the scan's own scratch_b allocations unwind here.
    ScratchScope b_guard(scratch_arena_b);
    // The frame-local mesh and its source fill scratch_a at the 1082-face peak.
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
      MeshOps::transform_in_place(mesh, ripple_gen, camera);
    }
    const SegueT &seg = carousel.segue();

    {
      HS_PROFILE(is_build_scan);
      // Rasterize from scratch_b: the swept+compiled mesh fills scratch_a to
      // ~120.9 KB during a build leg, leaving no room for the scan's per-face
      // SDF::FaceScratchBuffer. The sprite path scans from scratch_a, where its
      // transformed copy already lives.
      Scan::Mesh::draw_opleg_shading<W, H>(filters, canvas, mesh, sh, 1.0f,
                                           seg.opacity(1.0f), scratch_arena_b);
    }
  }

  /**
   * @brief Advances to the next registry solid and spawns it.
   */
  HS_COLD_MEMBER void spawn_shape() {
    auto solids = Solids::Collections::get_islamic_solids();
#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
    solid_idx = HS_ISLAMICSTARS_PROFILE_SHAPE;
#else
    solid_idx = (solid_idx + 1) % solids.size();
#endif
    spawn_entry(solids[solid_idx]);
  }

  /**
   * @brief Re-splits the arenas for the shape about to be generated.
   * @param has_recipe Whether the shape builds through a swept recipe chain.
   * @details Valid only with persistent at its ~baseline and both scratch
   * arenas idle, as the caller's compact leaves them. A smooth kis/needle
   * bridge shape gets the scratch_a-heavy split; other recipes trade unused
   * scratch_b for persistent; whole-generated shapes keep the full generation
   * scratch_b. Persistent takes the remainder.
   */
  HS_COLD_MEMBER void resplit_for_spawn(bool has_recipe) {
    const bool bridge_split = has_recipe && build_uses_smooth_bridge();
    const ArenaBudget budget = bridge_split ? BRIDGE_BUDGET
                               : has_recipe ? RECIPE_BUDGET
                                            : GENERATED_BUDGET;
    device_persistent_budget = budget.persistent();
    resplit_arenas(budget.persistent(GLOBAL_ARENA_SIZE), budget.scratch_a,
                   budget.scratch_b);
  }

  /**
   * @brief Generates @p entry into the carousel's back slot with a freshly
   *        shuffled palette, makes it the front, schedules the segue and the
   *        shape's mid-display ripple burst, and queues the next spawn_shape
   *        call.
   * @param entry Solid spawned; its recipe, if any, drives the build chain.
   */
  HS_COLD_MEMBER void spawn_entry(const Solids::Entry &entry) {
    build_entry = &entry;
    int back = 1 - carousel.front_index();
    // The spawned mesh's shuffled palette order, consumed by class ordinal.
    std::array<int, NUM_PALETTES> palette_slots;
    MeshPaletteBank::shuffle_indices(palette_slots);

    // A recipe whose lowered chain contains a step no leg kind covers falls
    // back to today's whole-generate path, seed solid and all.
    const Solids::Recipe *recipe = entry.recipe;
    build_step_count = 0;
    if (recipe) {
      build_step_count = Solids::expand_to_primitives(*recipe, build_step_chain,
                                                      MAX_BUILD_STEPS);
      for (size_t k = 0; k < build_step_count; ++k) {
        if (!Solids::is_morphable_step(build_step_chain[k])) {
          hs::log("IslamicStars: %s has an unsweepable step, generating whole",
                  entry.name);
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
    // tracking them through the evacuation. A build regenerates both slots
    // before either is drawn again, so the outgoing shape is dropped.
    auto rebake = [this](Arena &arena) { reclaim_persistent(arena); };
    if (recipe)
      carousel.compact_drop_all(rebake);
    else
      carousel.compact_keep_front(back, rebake);

    resplit_for_spawn(recipe != nullptr);

    hs::generate(persistent_arena, [&](Arena &target, Arena &a,
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
        PolyMesh mesh = entry.generate(a, b);
        carousel.slot(back).clear();
        MeshOps::compile(mesh, carousel.slot(back), target, scratch_arena_a);
      }
    });

    MeshOps::classify_faces_by_topology(carousel.slot(back), scratch_arena_a,
                                        scratch_arena_b, persistent_arena);

    // Colours of the spawned mesh (a build's seed, or the whole solid): the
    // shuffled palette order consumed by class ordinal — class ids are dense,
    // so class c is the c-th cohort.
    {
      const MeshState &spawned_mesh = carousel.slot(back);
      const size_t spawned_faces = spawned_mesh.topology.size();
      HS_CHECK(spawned_faces <= MAX_BUILD_FACES,
               "IslamicStars: spawned mesh exceeds MAX_BUILD_FACES");
      MeshPaletteBank::assign_by_class(spawned_mesh.topology.data(),
                                       spawned_faces, palette_slots,
                                       slot_face_palette[back]);
    }

    // Flip front eagerly for the overlapping sprite.
    carousel.set_front(back);

    // Segues with a spatial anchor (sweep axis, wave origin, spin axis) get a
    // fresh random one per transition. Safe mid-carousel: those segues are
    // sequential, so the previous sprite has already finished.
    if constexpr (requires(SegueT &s, const math::Vector &v) { s.retarget(v); })
      carousel.segue().retarget(math::random_vector());

    // Per-shape choreography: segue in, hold still one second, ripple, settle
    // one second, segue out. Duration is derived from the stage lengths so the
    // stages never overlap. Trans Speed divides every stage length, each with a
    // >=1-frame floor. The effective ripple duration/stagger are cached for the
    // deferred ripple() callback, which fires before the next shape spawns.
    const float sp = std::max(1.0f, params.trans_speed);
    int fade = std::max(1, static_cast<int>(SPRITE_FADE_FRAMES / sp));
    int still = std::max(1, static_cast<int>(STILL_FRAMES / sp));
    ripple_dur_eff = std::max(8, static_cast<int>(params.ripple_duration / sp));
    ripple_stagger_eff =
        std::max(1, static_cast<int>(RIPPLE_STAGGER_FRAMES / sp));
    burst_size_eff = params.burst_size;
    int burst_span = (burst_size_eff - 1) * ripple_stagger_eff + ripple_dur_eff;

    // Recipe entries insert a build phase on the segue's phase-1 plateau:
    // duration is lengthened by the build span rather than the carousel
    // growing an asymmetric-window API.
    const int build_span = recipe ? plan_build_legs(sp) : 0;

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
   * @brief Slider-backed runtime parameters for the effect.
   */
  struct Params {
    uint8_t burst_size = 4; /**< Ripples per burst. */
    float ripple_duration =
        80.0f; /**< Frames each ripple takes to expand across the sphere. */
    float trans_speed =
        1.0f; /**< Divides every per-shape stage length (fade, still holds, ripple span) and every build-leg budget: 1 = shipping cadence, higher cycles shapes faster. */
  } params;
};

#include "core/control/registry.h"
