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

#include "core/render/scan/mesh_motion.h"
#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "core/animation/animation.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/animation/recipe_build.h"

namespace IslamicStarsDetail {
/// Playback roster in preset order.
inline constexpr std::string_view SHAPE_NAMES[] = {
    "dodecahedron_hk62_ambo_hk62",
    "truncatedIcosahedron_hk58_chamfer63",
    "dodecahedron_bevel2_relax_gyro",
    "truncatedCube_hk90_meta",
    "truncatedIcosahedron_ambo_relax_truncate001_hankin59",
    "icosahedron_ambo_truncate033_hankin59",
    "dodecahedron_hk35_ambo_hk62_ambo_relax_hk42",
    "octahedron_hk17_ambo_hk73",
    "truncatedIcosidodecahedron_truncate50d_ambo_dual",
    "icosidodecahedron_truncate5d_ambo_dual",
    "snubDodecahedron_truncate5d_ambo_dual",
    "rhombicuboctahedron_hk63_ambo_hk63",
    "dodecahedron_hk54_ambo_hk72",
    "dodecahedron_hk72_ambo_dual_hk20",
    "icosahedron_snub_relax_truncate033_hankin62",
    "icosahedron_truncate25_hk60_ambo_relax100_hk60"};
/// Lowered primitive capacity of an on-screen recipe build.
inline constexpr size_t MAX_BUILD_OPS = 8;
/// Mesh face capacity of an on-screen recipe build.
inline constexpr size_t MAX_BUILD_FACES = 1152;
} // namespace IslamicStarsDetail

namespace hs_test {
namespace effects_tests {
struct IslamicBuildProbe;
} // namespace effects_tests
} // namespace hs_test

/** @brief IslamicStars' slider-bound parameters, shared by every preset. */
struct IslamicStarsParams {
  uint8_t burst_size = 4; /**< Ripples per burst. */
  float ripple_duration =
      80.0f; /**< Frames each ripple takes to expand across the sphere. */
  float trans_speed =
      1.0f; /**< Divides every per-shape stage length and build-leg budget. */
};

/**
 * @brief Effect that displays a sequence of Islamic-geometry polyhedra,
 *        transitioning one shape into the next while ripples distort the mesh.
 * @details Each roster shape is one preset. Entries with a non-null recipe
 *          are built op by op on screen: the seed solid sweeps in, then OpLegs
 *          morph the lowered chain into the finished pattern. Each shape's
 *          segue schedules the next automatic advance; a manual or synchronized
 *          change cuts the current shape and spawns the selected one. Pause
 *          holds a shape at the end of its plateau with repeating ripple bursts,
 *          before its segue out.
 * @tparam W Target canvas width in pixels.
 * @tparam H Target canvas height in pixels.
 */
template <int W, int H>
class IslamicStars
    : public ChoreographedEffect<IslamicStars<W, H>, IslamicStarsParams>,
      private Animation::RecipeBuild<IslamicStars<W, H>,
                                     IslamicStarsDetail::MAX_BUILD_OPS,
                                     IslamicStarsDetail::MAX_BUILD_FACES> {
  using Choreography =
      ChoreographedEffect<IslamicStars<W, H>, IslamicStarsParams>;
  friend Choreography;

public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::IslamicStars::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /// Live parameter set.
  using Params = IslamicStarsParams;
  /// Snapshot schema version; changes with the `Params` layout.
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** Bookkeeping only; each shape's segue schedules the next advance. */
  static constexpr uint16_t PRESET_DWELL_FRAMES = 192;

  /** @brief Whether every parameter lies inside its slider range.
   *  @param p Parameter set to check.
   *  @return True if all fields are in range. */
  static constexpr bool valid_params(const Params &p) {
    return p.burst_size >= 1 && p.burst_size <= BURST_MAX &&
           p.ripple_duration >= RIPPLE_DURATION_MIN &&
           p.ripple_duration <= RIPPLE_DURATION_MAX &&
           p.trans_speed >= TRANS_SPEED_MIN && p.trans_speed <= TRANS_SPEED_MAX;
  }

#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
  static_assert(
      HS_ISLAMICSTARS_PROFILE_SHAPE >= 0 &&
          HS_ISLAMICSTARS_PROFILE_SHAPE <
              static_cast<int>(std::size(IslamicStarsDetail::SHAPE_NAMES)),
      "HS_ISLAMICSTARS_PROFILE_SHAPE is outside the IslamicStars roster");
#endif

  /// Scratch split, in bytes, for a shape without a recipe.
  static constexpr ArenaSplit GENERATED_BUDGET{116 * 1024, 74 * 1024};
  /// Scratch split, in bytes, for a recipe build.
  static constexpr ArenaSplit RECIPE_BUDGET{116 * 1024, 72 * 1024};
  /// Scratch split, in bytes, for a recipe build with a smooth bridge.
  static constexpr ArenaSplit BRIDGE_BUDGET{129 * 1024 + 512, 74 * 1024};
  static_assert(BRIDGE_BUDGET.scratch_a + BRIDGE_BUDGET.scratch_b <
                DEVICE_GLOBAL_ARENA_SIZE);

  /**
   * @brief Constructs the effect, binding the ripple generator to the timeline.
   */
  HS_FLASH_MEMBER IslamicStars()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  HS_FLASH_MEMBER void init() override {
    begin_choreography();
    GENERATED_BUDGET.configure();
    device_persistent_budget = GENERATED_BUDGET.device_persistent();

    ripple_gen.init_storage(persistent_arena);
    claim_face_palettes(persistent_arena);
    palette_bank.bake_all(persistent_arena);

    ripple_gen.template_params.amplitude = RIPPLE_AMP_MAX;
    ripple_gen.template_params.thickness = RIPPLE_THICKNESS;
    ripple_gen.template_params.decay = 0.1f;
#ifdef HS_PROFILE_TRANS_SPEED
    static_assert(HS_PROFILE_TRANS_SPEED >= TRANS_SPEED_MIN &&
                      HS_PROFILE_TRANS_SPEED <= TRANS_SPEED_MAX,
                  "HS_PROFILE_TRANS_SPEED must be in [1, 8]");
    params.trans_speed = static_cast<float>(HS_PROFILE_TRANS_SPEED);
#endif

    this->template register_described_params<&IslamicStars::segue_fields>(
        carousel.segue());
    this->template register_described_params<&IslamicStars::burst_fields>(
        params);
    this->template register_described_params<&IslamicStars::ripple_fields>(
        ripple_gen.template_params);
    this->template register_described_params<&IslamicStars::timing_fields>(
        params);

    add_orientation_walk();

#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
    preset_index = HS_ISLAMICSTARS_PROFILE_SHAPE;
#ifdef HS_PROFILE_ENABLE
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(preset_index),
            static_cast<unsigned>(PRESETS.size()));
#endif
#else
    // Open on a recipe entry so the op-by-op build is the first thing drawn.
    for (size_t i = 0; i < PRESETS.size(); ++i) {
      if (PRESETS[i].params->recipe) {
        preset_index = i;
        break;
      }
    }
#endif

    spawn_shape();
  }

  /**
   * @brief Advances ripple state once and runs the timeline for this frame.
   */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    plateau_hold.step(anims_paused);
    if (plateau_hold.held() && ripple_gen.active_count() == 0)
      ripple(canvas);
    {
      HS_PROFILE(is_ripple_prepare);
      ripple_gen.prepare_frame();
    }
    {
      HS_PROFILE(is_timeline_step);
      timeline.step(canvas);
    }
    this->prepare_build();
  }

private:
  using Builder = Animation::RecipeBuild<IslamicStars<W, H>,
                                         IslamicStarsDetail::MAX_BUILD_OPS,
                                         IslamicStarsDetail::MAX_BUILD_FACES>;
  friend Builder;
  using Choreography::anims_paused;
  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::preset_index;
  using Choreography::timeline;
  using Builder::abandon_build;
  using Builder::build_active;
  using Builder::build_entry;
  using Builder::build_step_chain;
  using Builder::build_step_count;
  using Builder::build_seed;
  using Builder::device_persistent_budget;
#if HS_ENABLE_TEST_HOOKS
  using Builder::dual_bridges_built;
#endif
  using typename Builder::BuildContinuation;
  using Builder::check_build_budget;
  using Builder::schedule_dual_bridge;
  using Builder::clean_endpoint;
  using Builder::build_uses_smooth_bridge;
  using Builder::plan_build_legs;
  using Builder::start_build_leg;
  static constexpr size_t MAX_BUILD_STEPS = IslamicStarsDetail::MAX_BUILD_OPS;
  static constexpr size_t MAX_BUILD_FACES = IslamicStarsDetail::MAX_BUILD_FACES;
  static_assert(Solids::max_lowered_step_count(Solids::islamic_registry) <=
                MAX_BUILD_STEPS);
  friend struct ::hs_test::effects_tests::IslamicBuildProbe;

  // Two-burst capacity reserves headroom for the previous shape's live ripples.
  static constexpr int RIPPLE_POOL_SIZE = 8;
  static constexpr int RIPPLE_STAGGER_FRAMES = 16;
  static constexpr float RIPPLE_DURATION_MIN = 30.0f;
  static constexpr float RIPPLE_DURATION_MAX = 143.0f;
  static constexpr int BURST_MAX = 4;
  static constexpr float TRANS_SPEED_MIN = 1.0f;
  static constexpr float TRANS_SPEED_MAX = 8.0f;
  static constexpr int SPRITE_FADE_FRAMES = 16;
  static constexpr int STILL_FRAMES =
      16; /**< Hold frames between fade and ripple, before the Trans Speed divisor. */
  static constexpr float RIPPLE_THICKNESS =
      0.7f; /**< Fixed ripple wavelet width (radians). */
  /** Amplitude ceiling; every ripple rotation takes the series-form
   * quaternion. The displacement map d -> d + theta(d) is injective only below
   * amp/thickness = 0.181, so it folds slightly at this ceiling. */
  static constexpr float RIPPLE_AMP_MAX = RIPPLE_SMALL_ANGLE_MAX;
  static_assert(
      2 * BURST_MAX <= RIPPLE_POOL_SIZE,
      "IslamicStars: ripple pool must reserve capacity for two bursts");

  math::Orientation<> orientation;
  FastNoiseLite noise;
  Pipeline<W, H> filters;
  RippleTransformer<RIPPLE_POOL_SIZE> ripple_gen;
  // Effective per-shape stage lengths after the Trans Speed divisor.
  int ripple_dur_eff = 80;
  int ripple_stagger_eff = RIPPLE_STAGGER_FRAMES;
  int burst_size_eff = 4;
  /** Departure gate at the sprite plateau boundary. */
  Segue::PlateauHold plateau_hold;
  using SegueT = Segue::TerminatorSweep;

  MeshCarousel<SegueT> carousel;

  /** @brief Descriptions of the global per-face fade range parameters.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto segue_fields() {
    return std::tuple{
        // Per-face fade length range (frames): each face draws a random fade
        // from [lo, hi] as the terminator reaches it, fraying the sweep front.
        Control::Field<SegueT, float>{
            .id = "fade_frames_min",
            .member = &SegueT::fade_frames_min,
            .name = "Face Fade Lo",
            .spec = {.min = 0.0f, .max = 32.0f, .preset = false}},
        Control::Field<SegueT, float>{
            .id = "fade_frames_max",
            .member = &SegueT::fade_frames_max,
            .name = "Face Fade Hi",
            .spec = {.min = 0.0f, .max = 32.0f, .preset = false}}};
  }

  /** @brief Description of the ripple burst size parameter.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto burst_fields() {
    return std::tuple{Control::Field<IslamicStarsParams, uint8_t>{
        .id = "burst_size",
        .member = &IslamicStarsParams::burst_size,
        .name = "Burst",
        .spec = {.min = 1, .max = BURST_MAX}}};
  }

  /** @brief Descriptions of the global ripple template parameters.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto ripple_fields() {
    return std::tuple{
        // Thickness is fixed, so RIPPLE_AMP_MAX bounds amp/thickness.
        Control::Field<Animation::RippleParams, float>{
            .id = "amplitude",
            .member = &Animation::RippleParams::amplitude,
            .name = "Ripp Amp",
            .spec = {.min = 0.0f, .max = RIPPLE_AMP_MAX, .preset = false}},
        Control::Field<Animation::RippleParams, float>{
            .id = "decay",
            .member = &Animation::RippleParams::decay,
            .name = "Ripp Decay",
            .spec = {.min = 0.0f, .max = 5.0f, .preset = false}}};
  }

  /** @brief Descriptions of the ripple duration and transition speed
   *         parameters.
   *  @return Tuple of `Control::Field` descriptors, one per parameter. */
  static constexpr auto timing_fields() {
    return std::tuple{
        Control::Field<IslamicStarsParams, float>{
            .id = "ripple_duration",
            .member = &IslamicStarsParams::ripple_duration,
            .name = "Ripp Dur",
            .spec = {.min = RIPPLE_DURATION_MIN, .max = RIPPLE_DURATION_MAX}},
        Control::Field<IslamicStarsParams, float>{
            .id = "trans_speed",
            .member = &IslamicStarsParams::trans_speed,
            .name = "Trans Speed",
            .spec = {.min = TRANS_SPEED_MIN, .max = TRANS_SPEED_MAX}}};
  }

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank;
  /** Per-slot per-face palette ids (persistent-arena backed, MAX_BUILD_FACES
   * each); compaction re-claims the same addresses, so the contents survive
   * the reset. */
  uint8_t *slot_face_palette[2] = {};
  /** Per-slot segue phase of the last sprite callback. */
  Spacetime::Retained<float> phase_history[2];

  /**
   * @brief Claims the two per-slot face-palette arrays from the arena.
   * @param arena Persistent arena the arrays live in.
   * @details Must directly follow the ripple pool's claim so the arrays
   * re-land at their original addresses (asserted).
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
  HS_FLASH_MEMBER void reclaim_persistent(Arena &arena) {
    ripple_gen.reclaim_storage(arena);
    claim_face_palettes(arena);
    palette_bank.bake_all(arena);
  }

  /**
   * @brief Spawns one burst of burst_size ripples from a random origin,
   *        staggered ripple_stagger_eff frames apart, each expanding over
   *        ripple_dur_eff frames.
   */
  void ripple(Canvas &) {
    math::Vector origin = math::random_vector();
    for (int i = 0; i < burst_size_eff; i++) {
      if (!ripple_gen.spawn(i * ripple_stagger_eff, origin,
                            math::PI_F / ripple_dur_eff, ripple_dur_eff))
        hs::log("IslamicStars: ripple pool full, dropping spawn");
    }
  }

  /** @brief Consecutive endpoint ripples sharing one rejection band. */
  using RippleGroup = RippleKernel::EndpointGroup<Animation::RippleParams>;

  /** @brief Endpoint ripple chain with shared-center group rejection. */
  struct EndpointRipples {
    RippleGroup *waves = nullptr; ///< Arena-backed waves and group bounds.
    int count = 0;                ///< Number of nonzero waves.

    /** @brief Applies the endpoint chain in spawn order. */
    HS_O3_FN math::Vector operator()(math::Vector v) const {
      return RippleKernel::sequential_groups(
          v, waves, count,
          []() HS_HOT_INLINE {
#if HS_ENABLE_TEST_HOOKS
            ++test_endpoint_groups;
#endif
          },
          []() HS_HOT_INLINE {
#if HS_ENABLE_TEST_HOOKS
            ++test_endpoint_bands;
#endif
          });
    }
  };

  /** @brief Prepares endpoint rejection groups in the draw's scratch scope. */
  HS_HOT_FLASH_MEMBER EndpointRipples endpoint_ripples(Arena &arena) const {
    EndpointRipples result;
    if (ripple_gen.active_count() == 0)
      return result;
    result.waves = arena.allocate_n<RippleGroup>(ripple_gen.active_count());
    int first = 0;
    for (int k = 0; k < ripple_gen.active_count(); ++k) {
      const auto &p = ripple_gen.active_params(k);
      if (p.amplitude <= RippleKernel::NEGLIGIBLE_AMPLITUDE)
        continue;
      RippleKernel::append_endpoint(result.waves, result.count, first, p);
    }
    return result;
  }
#if HS_ENABLE_TEST_HOOKS
  inline static bool test_reference_endpoints = false; ///< Ungrouped chain.
  inline static size_t test_endpoint_groups = 0; ///< Group rejection tests.
  inline static size_t test_endpoint_bands = 0;  ///< Individual band tests.
#endif

  /**
   * @brief Orients and ripple-distorts a mesh into scratch_arena_a.
   * @param base_state Undistorted source mesh.
   * @param visible Per-vertex board visibility bits, or null.
   * @return The transformed mesh (scratch_arena_a-backed; the caller holds the
   * scope).
   */
  HS_O3_FN MeshState transform_shape(const MeshState &base_state,
                                     const uint8_t *visible = nullptr) {
    MeshState transformed_state;
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
#if HS_ENABLE_TEST_HOOKS
      if (test_reference_endpoints) {
        MeshOps::transform(base_state, transformed_state, scratch_arena_a,
                           ripple_gen, camera);
        for (int k = 0; k < ripple_gen.active_count(); ++k)
          if (ripple_gen.active_params(k).amplitude >
              RippleKernel::NEGLIGIBLE_AMPLITUDE)
            test_endpoint_bands += base_state.vertices.size();
      } else
#endif
      {
        const auto ripples = endpoint_ripples(scratch_arena_a);
        MeshOps::transform(base_state, transformed_state, scratch_arena_a);
        MeshOps::transform_selected_in_place(transformed_state, visible,
                                             ripples, camera);
      }
    }
    return transformed_state;
  }

  // Scan::plan_motion()'s cap radii are freed before the scan's face scratch and
  // face selection take their place.
  static_assert(MAX_BUILD_FACES * sizeof(float) <=
                MAX_BUILD_FACES + sizeof(SDF::FaceStorage));

  /**
   * @brief Sprite draw callback: draws the held shape for one envelope frame.
   * @param canvas Render target.
   * @param phase Sprite envelope phase: rises over the incoming window, holds 1,
   *        falls over the outgoing window.
   * @param back Carousel slot the shape was spawned into.
   * @details Draws nothing during the build window.
   */
  HS_FLASH_MEMBER void draw_sprite(Canvas &canvas, float phase, int back) {
    const float previous_phase = record_phase(back, phase);
    if (build_active)
      return;
    const MeshState &mesh = carousel.slot(back);
    draw_shape(canvas, phase, previous_phase, mesh, slot_face_palette[back]);
  }

  /**
   * @brief Records a slot's sprite phase for this frame.
   * @param slot Carousel slot whose sprite draws.
   * @param phase Envelope phase now.
   * @return The phase the slot's sprite drew last, or @p phase for its first
   * draw since spawning and with temporal sampling off.
   */
  HS_FLASH_MEMBER float record_phase(int slot, float phase) {
    const float previous = phase_history[slot].start(phase);
    phase_history[slot].capture(phase);
    return previous;
  }

  /**
   * @brief Per-face segue opacities at the previous frame endpoint.
   * @param base_state Undistorted mesh the sweep orders faces on.
   * @param transformed_state Displayed mesh carrying the face classes.
   * @param previous_phase Sprite envelope phase a frame ago.
   * @return One opacity per face, backed by scratch_arena_a for the caller's
   * scope on it.
   */
  HS_FLASH_MEMBER const float *
  previous_opacities(const MeshState &base_state,
                     const MeshState &transformed_state, float previous_phase) {
    ArenaVector<float> phases;
    phases.bind(scratch_arena_a, transformed_state.num_faces());
    carousel.fill_face_phases(base_state, transformed_state, previous_phase,
                              phases);
    const SegueT &seg = carousel.segue();
    for (float &value : phases)
      value = seg.opacity(value);
    return phases.data();
  }

  /**
   * @brief Orients, ripple-distorts, and segue-shapes base_state, then
   *        rasterizes it with a per-face palette lookup, each column seeing
   *        the orientation at its sweep time and, while ripples move, the
   *        ripples at their capture times.
   * @param canvas Render target receiving the rasterized mesh.
   * @param phase Segue phase in [0, 1] from the sprite envelope: rises over
   *        the incoming window, holds 1, falls over the outgoing window.
   * @param previous_phase Segue phase a frame ago; each column shades each
   *        face's opacity at its sweep time between the two.
   * @param base_state Undistorted source mesh to transform and draw; carries
   *        the per-face topology classes.
   * @param face_palette Per-face palette ids.
   */
  HS_HOT_FLASH_MEMBER void draw_shape(Canvas &canvas, float phase,
                                      float previous_phase,
                                      const MeshState &base_state,
                                      const uint8_t *face_palette) {
    HS_CHECK(base_state.get_topology_size() == base_state.num_faces(),
             "IslamicStars: sprite shading face count mismatch");
    const SegueT &seg = carousel.segue();
    if (!seg.visible(phase) && !seg.visible(previous_phase))
      return;

    HS_PROFILE(is_draw_shape);
    ScratchScope a_guard(scratch_arena_a);
    const Scan::ModelPath path{base_state.vertices.data()};
    const float bound =
        Scan::motion_bound(path, base_state.vertices.size(), ripple_gen);
    Scan::ColumnBasis column_basis = Scan::prepare_column_basis<W, H>(
        scratch_arena_a, orientation,
        previous_phase != phase || Scan::mesh_motion_active<W>(bound));
    Scan::select_mesh_visibility<W, H>(
        canvas.clip(), base_state, OrientTransformer(orientation),
        Scan::visibility_reach(ripple_gen, bound), scratch_arena_a,
        column_basis);
    Scan::ColumnBasisScope column_scope(canvas, column_basis);
    MeshState transformed_state =
        transform_shape(base_state, column_basis.mesh.visible_vertices);
    ScratchScope b_guard(scratch_arena_b);
    const Scan::EndpointCamera camera(orientation);
    const Scan::MotionPlan plan =
        Scan::plan_motion<W>(canvas, transformed_state, path, bound, ripple_gen,
                             camera, scratch_arena_a, scratch_arena_b);
    ArenaVector<float> face_phases;
    {
      HS_PROFILE(is_face_offsets);
      const size_t faces = transformed_state.num_faces();
      face_phases.bind(scratch_arena_a, faces);
      carousel.fill_face_phases(base_state, transformed_state, phase,
                                face_phases);
    }
    const float *previous_opacity =
        previous_phase == phase
            ? nullptr
            : previous_opacities(base_state, transformed_state, previous_phase);

    {
      HS_PROFILE(is_mesh_scan);
      FacePaletteShader fragment_shader;
      auto select_face = [&](size_t fi, float size) {
        HS_CHECK(fi < face_phases.size(),
                 "IslamicStars: sprite shading face mismatch");
        fragment_shader.set_palette(&palette_bank[face_palette[fi]].view());
        const float opacity = seg.opacity(face_phases[fi]);
        fragment_shader.set_alpha(
            opacity, previous_opacity ? previous_opacity[fi] : opacity);
        fragment_shader.scale = size > math::TOLERANCE ? 1.0f / size : 0.0f;
      };
      Scan::deform_planned<W, H>(
          canvas, transformed_state, plan, &path, ripple_gen, camera,
          scratch_arena_a, bound, [&]() HS_HOT_INLINE {
            Scan::Mesh::draw_specialized<W, H>(
                filters, canvas, transformed_state, fragment_shader,
                scratch_arena_a, nullptr, select_face);
          });
    }
  }

  /**
   * @brief Draws one build-leg frame: the compiled swept mesh, shaded from its
   *        per-face palettes on the segue's phase-1 identity plateau, each
   *        column seeing the orientation at its sweep time and, while the
   *        mesh or a ripple moves, the mesh at their capture times.
   * @param canvas Render target receiving the rasterized mesh.
   * @param mesh Compiled swept mesh (scratch-backed, this frame only).
   * @param sh Per-face palette table and previous vertices from the OpLeg.
   */
  HS_HOT_FLASH_MEMBER void
  draw_build_mesh(Canvas &canvas, MeshState &mesh,
                  const Animation::OpLeg::Shading &sh) {
    if (mesh.vertices.is_empty())
      return;
    HS_PROFILE(is_build_draw);
    // Opened after the OpLeg's blended ramps, which the shader reads.
    ScratchScope b_guard(scratch_arena_b);
    const size_t count = mesh.vertices.size();
    Scan::ModelPath path = Scan::ModelPath::of(mesh.vertices.data(), sh);
    const float bound = Scan::motion_bound(path, count, ripple_gen);
    const bool active = Scan::mesh_motion_active<W>(bound);
    // A leg moving vertex by vertex places from its own path; otherwise the
    // scan keeps the current endpoint's model-space vertices, packed.
    if (active && !path.leg) {
      math::Snorm3 *model = scratch_arena_b.allocate_n<math::Snorm3>(count);
      for (size_t i = 0; i < count; ++i)
        model[i] = math::Snorm3::encode(mesh.vertices[i]);
      path.packed = model;
      path.current = nullptr;
    }
    Scan::ColumnBasis column_basis = Scan::prepare_column_basis<W, H>(
        scratch_arena_b, orientation,
        sh.previous_weight != sh.weight || active);
    Scan::select_mesh_visibility<W, H>(
        canvas.clip(), mesh, OrientTransformer(orientation),
        Scan::visibility_reach(ripple_gen, bound), scratch_arena_b,
        column_basis);
    Scan::ColumnBasisScope column_scope(canvas, column_basis);
    // The frame-local mesh and its source occupy scratch_a.
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
#if HS_ENABLE_TEST_HOOKS
      if (test_reference_endpoints) {
        MeshOps::transform_in_place(mesh, ripple_gen, camera);
        for (int k = 0; k < ripple_gen.active_count(); ++k)
          if (ripple_gen.active_params(k).amplitude >
              RippleKernel::NEGLIGIBLE_AMPLITUDE)
            test_endpoint_bands += mesh.vertices.size();
      } else
#endif
      {
        ScratchScope ripple_scope(scratch_arena_a);
        const auto ripples = endpoint_ripples(scratch_arena_a);
        MeshOps::transform_selected_in_place(
            mesh, column_basis.mesh.visible_vertices, ripples, camera);
      }
    }
    HS_PROFILE(is_build_scan);
    scan_build_mesh(canvas, mesh, sh, active ? &path : nullptr, bound);
  }

  /**
   * @brief Scans a build-leg frame's placed mesh, in vertex deformation while it
   *        moves.
   * @param canvas Render target with the frame's column-time state.
   * @param mesh Compiled swept mesh, placed at the current endpoint.
   * @param sh Per-face palette table and previous vertices from the OpLeg.
   * @param path The mesh's model path, or null when the frame cannot need
   *        deformation.
   * @param bound Scan::motion_bound() of the draw.
   */
  HS_FLASH_MEMBER void scan_build_mesh(Canvas &canvas, MeshState &mesh,
                                       const Animation::OpLeg::Shading &sh,
                                       const Scan::ModelPath *path,
                                       float bound) {
    const Scan::EndpointCamera camera(orientation);
    const Scan::MotionPlan plan =
        path ? Scan::plan_motion<W>(canvas, mesh, *path, bound, ripple_gen,
                                    camera, scratch_arena_b, scratch_arena_b)
             : Scan::MotionPlan{};
    const SegueT &seg = carousel.segue();
    // Per-face scan scratch must remain disjoint from the mesh in scratch_a.
    Scan::deform_planned<W, H>(canvas, mesh, plan, path, ripple_gen, camera,
                               scratch_arena_b, bound, [&]() HS_HOT_INLINE {
                                 Scan::Mesh::draw_opleg_shading<W, H>(
                                     filters, canvas, mesh, sh, 1.0f,
                                     seg.opacity(1.0f), scratch_arena_b);
                               });
  }

  /** @brief One snapping preset per playback roster shape. */
  static constexpr auto PRESETS = [] {
    std::array<PresetEntry<const Solids::Entry *>,
               std::size(IslamicStarsDetail::SHAPE_NAMES)>
        rows{};
    for (size_t i = 0; i < rows.size(); ++i)
      for (const auto &entry : Solids::islamic_registry)
        if (IslamicStarsDetail::SHAPE_NAMES[i] == entry.name)
          rows[i] = {&entry, Segue::Preset::Snap{}};
    return rows;
  }();
  static_assert(
      [] {
        for (const auto &preset : PRESETS)
          if (!preset.params)
            return false;
        return true;
      }(),
      "IslamicStars roster names must identify registered solids");

  /** @brief Startup parameters: the slider defaults. */
  static Params initial_params() { return {}; }

  /** @brief Every preset keeps the live slider values. */
  Params preset_params(size_t) const { return params; }

  /** @brief Restarts the shape sequence on a manual, synchronized, or restored change.
   *  @param change The committed preset change. */
  HS_FLASH_MEMBER void
  preset_changed(const Effect::PresetChange &change) override {
    if (change.origin != Effect::PresetChangeOrigin::AUTOMATIC)
      restart_shape();
  }

  /** @brief Adds the camera's random walk to the timeline. */
  HS_FLASH_MEMBER void add_orientation_walk() {
    timeline.add(0, Animation::RandomWalk<W>(orientation, math::UP, noise));
  }

  /**
   * @brief Cuts the resident shape and spawns the committed preset.
   * @details Clears every scheduled event, including any in-flight build leg
   * and ripple; the orientation walk resumes from the held orientation.
   */
  HS_FLASH_MEMBER void restart_shape() {
    timeline.clear();
    abandon_build();
    add_orientation_walk();
    spawn_shape();
  }

  /** @brief Spawns the committed preset's solid. */
  HS_FLASH_MEMBER void spawn_shape() {
    spawn_entry(*PRESETS[preset_index].params);
  }

  /** @brief Advances to the next preset and spawns it. */
  HS_FLASH_MEMBER void advance_shape() {
#ifndef HS_ISLAMICSTARS_PROFILE_SHAPE
    const bool advanced = this->advance_preset();
    HS_CHECK(advanced, "IslamicStars: automatic preset advance must succeed");
#endif
    spawn_shape();
  }

  /**
   * @brief Re-splits the arenas for the shape about to be generated.
   * @param has_recipe Whether the shape builds through a swept recipe chain.
   * @details Valid only after the caller's carousel compaction: both scratch
   * arenas empty and the live persistent offset within the chosen split's
   * persistent capacity.
   */
  HS_FLASH_MEMBER void resplit_for_spawn(bool has_recipe) {
    const bool bridge_split = has_recipe && build_uses_smooth_bridge();
    const ArenaSplit budget = bridge_split ? BRIDGE_BUDGET
                              : has_recipe ? RECIPE_BUDGET
                                           : GENERATED_BUDGET;
    device_persistent_budget = budget.device_persistent();
    resplit_arenas(budget.persistent(GLOBAL_ARENA_SIZE), budget.scratch_a,
                   budget.scratch_b);
  }

  /**
   * @brief Generates @p entry into the carousel's back slot, makes it the
   *        front, and schedules its segue, ripple burst and the next spawn.
   * @param entry Solid spawned; its recipe, if any, drives the build chain.
   */
  HS_FLASH_MEMBER void spawn_entry(const Solids::Entry &entry) {
    build_entry = &entry;
    int back = 1 - carousel.front_index();
    // The spawned mesh's shuffled palette order, consumed by class ordinal.
    std::array<int, NUM_PALETTES> palette_slots;
    MeshPaletteBank::shuffle_indices(palette_slots);

    // A recipe with a lowered step outside admitted recipe-morph coverage
    // falls back to the whole-generate path, seed solid and all.
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

    // Sequential segues leave no active sprite; the build rebuilds and selects
    // the back slot before drawing resumes.
    auto rebake = [this](Arena &arena) { reclaim_persistent(arena); };
    if (recipe)
      carousel.compact_drop_all(rebake);
    else
      carousel.compact_keep_front(back, rebake);

    resplit_for_spawn(recipe != nullptr);

    hs::generate(persistent_arena, [&](Arena &target, Arena &a,
                                       Arena &b) HS_FLASH_MEMBER {
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
    phase_history[back].reset();

    // Segues with a spatial anchor get a fresh random one per transition.
    if constexpr (requires(SegueT &s, const math::Vector &v) { s.retarget(v); })
      carousel.segue().retarget(math::random_vector());

    // Per-shape choreography: segue in, hold, ripple, settle, segue out. Trans
    // Speed divides every stage length, floored at 1 frame (ripple duration
    // at 8).
    const float sp = std::max(1.0f, params.trans_speed);
    int fade = std::max(1, static_cast<int>(SPRITE_FADE_FRAMES / sp));
    int still = std::max(1, static_cast<int>(STILL_FRAMES / sp));
    ripple_dur_eff = std::max(8, static_cast<int>(params.ripple_duration / sp));
    ripple_stagger_eff =
        std::max(1, static_cast<int>(RIPPLE_STAGGER_FRAMES / sp));
    burst_size_eff = params.burst_size;
    int burst_span = (burst_size_eff - 1) * ripple_stagger_eff + ripple_dur_eff;

    // Recipe entries insert a build phase on the segue's phase-1 plateau.
    const int build_span = recipe ? plan_build_legs(sp) : 0;

    int duration = fade + build_span + still + burst_span + still + fade;

    plateau_hold.reset();
    int next_delay = carousel.schedule_segue(timeline, back, draw_fn, duration,
                                             fade, plateau_hold.gate());

    // Added after the sprite, so the first leg draws the frame after the
    // sprite's last seed draw: no gap, no double draw.
    if (recipe) {
      timeline.add(fade, Animation::PeriodicTimer(0, [this](Canvas &) {
                     build_active = true;
                     start_build_leg();
                   }));
    }

    timeline.add(fade + build_span + still,
                 Animation::PeriodicTimer(
                     0, [this](Canvas &canvas) { ripple(canvas); }));

    plateau_hold.arm_plateau_end(timeline, duration - fade);

    // On a closed 2-manifold faces.size() (Σ face degrees) is exactly 2·E.
    // A recipe shape spawns holding its seed, so these are the seed's counts.
    const MeshState &spawned = carousel.current();
    hs::log("Spawning Shape: %s (V=%d, E=%d, F=%d, I=%d)%s", entry.name,
            (int)spawned.vertices.size(), (int)(spawned.faces.size() / 2),
            (int)spawned.face_counts.size(), (int)spawned.faces.size(),
            recipe ? " seed" : "");

    // The segue decides when the next shape starts relative to this one.
    timeline.add_pausable(next_delay,
                          Animation::PeriodicTimer(
                              0, [this](Canvas &) { this->advance_shape(); }),
                          plateau_hold.gate());
  }

public:
  /** @return Frames in one automatic roster cycle at shipping defaults. */
  static consteval size_t cycle_frames() {
    constexpr Params DEFAULTS{};
    const int BASE = 2 * SPRITE_FADE_FRAMES + 2 * STILL_FRAMES +
                     (DEFAULTS.burst_size - 1) * RIPPLE_STAGGER_FRAMES +
                     static_cast<int>(DEFAULTS.ripple_duration);
    size_t frames = 0;
    for (const auto &preset : PRESETS)
      frames += BASE + Builder::recipe_frame_budget(preset.params->recipe,
                                                    DEFAULTS.trans_speed);
    return frames;
  }

  /** @brief Explicit choreography customization bindings.
   * @return Type-only hook list. */
  static consteval auto choreography_hooks() {
    return Control::HookSet<
        IslamicStars,
        Control::HookBinding<Control::ChoreographyHook::VALID_PARAMS,
                             &IslamicStars::valid_params>,
        Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                             &IslamicStars::initial_params>,
        Control::HookBinding<Control::ChoreographyHook::PRESET_PARAMS,
                             &IslamicStars::preset_params>>{};
  }
};
