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
  HS_COLD_MEMBER IslamicStars()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        filters(), ripple_gen(timeline) {}

  /**
   * @brief Bakes palettes, registers the UI sliders, and seeds the timeline
   *        with the orientation walk and the first shape.
   */
  HS_COLD_MEMBER void init() override {
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

    // Per-face fade length range (frames): each face draws a random fade from
    // [lo, hi] as the terminator reaches it, fraying the sweep front.
    this->register_param("Face Fade Lo", &carousel.segue().fade_frames_min,
                         0.0f, 32.0f);
    this->register_param("Face Fade Hi", &carousel.segue().fade_frames_max,
                         0.0f, 32.0f);
    this->register_int_param("Burst", &params.burst_size, 1, BURST_MAX);
    // Thickness is fixed, so RIPPLE_AMP_MAX bounds amp/thickness.
    this->register_param("Ripp Amp", &ripple_gen.template_params.amplitude,
                         0.0f, RIPPLE_AMP_MAX);
    this->register_param("Ripp Decay", &ripple_gen.template_params.decay, 0.0f,
                         5.0f);
    for (const char *name :
         {"Face Fade Lo", "Face Fade Hi", "Ripp Amp", "Ripp Decay"})
      this->mark_global(name);
    this->register_param("Ripp Dur", &params.ripple_duration,
                         RIPPLE_DURATION_MIN, RIPPLE_DURATION_MAX);
    this->register_param("Trans Speed", &params.trans_speed, TRANS_SPEED_MIN,
                         TRANS_SPEED_MAX);

    add_orientation_walk();

#ifdef HS_ISLAMICSTARS_PROFILE_SHAPE
    preset_index = HS_ISLAMICSTARS_PROFILE_SHAPE;
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
  using Builder::build_step;
  using Builder::build_seed;
  using Builder::build_total_frames;
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

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  MeshPaletteBank palette_bank;
  /** Per-slot per-face palette ids (persistent-arena backed, MAX_BUILD_FACES
   * each); compaction re-claims the same addresses, so the contents survive
   * the reset. */
  uint8_t *slot_face_palette[2] = {};
  /** Per-slot segue phase of the sprite's last draw. */
  float drawn_phase[2] = {};
  /** Per-slot whether the slot's sprite has drawn since it spawned. */
  bool phase_drawn[2] = {};

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
  HS_COLD_MEMBER void reclaim_persistent(Arena &arena) {
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

  /**
   * @brief Orients and ripple-distorts a mesh into scratch_arena_a.
   * @param base_state Undistorted source mesh.
   * @return The transformed mesh (scratch_arena_a-backed; the caller holds the
   * scope).
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
   * @brief Bound on how far the active ripples move any point over the frame.
   * @return Radians.
   */
  HS_FLASH_MEMBER float ripple_arc() const {
    float arc = 0.0f;
    for (int k = 0; k < ripple_gen.active_count(); ++k)
      arc += ripple_gen.active_params(k).step_arc();
    return arc;
  }

  /** @brief Where a draw's model-space vertices stand over the frame. */
  struct ModelPath {
    const math::Vector *current = nullptr;  ///< At the current endpoint.
    const math::Vector *previous = nullptr; ///< At the previous endpoint.
    /** Shading whose leg moves the model vertex by vertex, or null. */
    const Animation::OpLeg::Shading *leg = nullptr;

    /** @return Whether the model moves over the frame. */
    bool moves() const { return leg || previous; }

    /** @brief One fraction of the frame, with the leg's parameter there. */
    struct Moment {
      float u = 0.0f; ///< Frame fraction; 0 and 1 are the endpoints.
      float k = 0.0f; ///< The leg's sweep parameter at `u`; 0 without a leg.
    };

    /**
     * @brief A fraction of the frame, resolved for at().
     * @param u Frame fraction; past either endpoint the model continues
     *        along its path, holding at a leg's own ends.
     * @return The moment.
     */
    Moment moment(float u) const {
      return {u, leg ? leg->parameter_at(u) : 0.0f};
    }

    /**
     * @brief A model vertex at a moment of the frame.
     * @param i Vertex index.
     * @param m The moment, as moment() gives it.
     * @return The vertex on its leg's path, normalized-lerped from
     * `previous`, or `current` when the model holds.
     */
    HS_HOT_INLINE math::Vector at(size_t i, const Moment &m) const {
      if (leg)
        return leg->vertex_at_parameter(i, m.k);
      if (previous)
        return (previous[i] + (current[i] - previous[i]) * m.u).normalized();
      return current[i];
    }

    /**
     * @brief A model vertex at its capture fraction.
     * @param i Vertex index.
     * @param u Frame fraction.
     * @return Vertex on the model path.
     */
    HS_HOT_INLINE math::Vector at(size_t i, float u) const {
      return at(i, moment(u));
    }
  };

  /**
   * @brief Bound on how far a mesh draw's vertices move over the frame.
   * @param path The draw's model path.
   * @param count Number of vertices.
   * @return Radians; the ripples' bound plus the model's largest motion.
   */
  HS_FLASH_MEMBER float motion_bound(const ModelPath &path,
                                     size_t count) const {
    float min_cos = 1.0f;
    if (path.moves()) {
      const auto start = path.moment(0.0f);
      const auto end = path.moment(1.0f);
      for (size_t i = 0; i < count; ++i)
        min_cos = fminf(min_cos, math::dot(path.at(i, start), path.at(i, end)));
    }
    return ripple_arc() + acosf(hs::clamp(min_cos, -1.0f, 1.0f));
  }

  /**
   * @brief Flags the vertices the model or an active ripple's last step can
   *        move.
   * @param path The draw's model path.
   * @param count Number of vertices.
   * @param leg_arc The model's largest motion (radians).
   * @param moving Receives 1 for each vertex the model moves more than
   *        Scan::HELD_VERTEX_ARC, or that stands within @p leg_arc of a
   *        moving ripple's band at either end of its step.
   * @param arena Arena holding the bands for the call.
   * @details A vertex outside every band is left in place by each ripple in
   * turn.
   */
  HS_FLASH_MEMBER void flag_moving(const ModelPath &path, size_t count,
                                   float leg_arc, uint8_t *moving,
                                   Arena &arena) const {
    struct Band {
      math::Vector center; ///< Ripple center.
      float cos_near;      ///< Cosine of the band's inner edge.
      float cos_far;       ///< Cosine of the band's outer edge.
    };
    ScratchScope scope(arena);
    Band *bands = arena.allocate_n<Band>(2 * ripple_gen.active_count() + 1);
    int band_count = 0;
    for (int k = 0; k < ripple_gen.active_count(); ++k) {
      const Animation::RippleParams &p = ripple_gen.active_params(k);
      if (p.step_arc() == 0.0f)
        continue;
      for (const float phase : {p.previous_phase, p.phase}) {
        const float reach = 2.0f * p.half_width() + leg_arc;
        bands[band_count++] = {
            p.center, cosf(hs::clamp(phase - reach, 0.0f, math::PI_F)),
            cosf(hs::clamp(phase + reach, 0.0f, math::PI_F))};
      }
    }
    const float held_cos = cosf(Scan::HELD_VERTEX_ARC<W>);
    const auto start = path.moment(0.0f);
    const auto end = path.moment(1.0f);
    for (size_t i = 0; i < count; ++i) {
      const math::Vector model = path.at(i, end);
      bool moves =
          path.moves() && math::dot(path.at(i, start), model) < held_cos;
      for (int b = 0; b < band_count && !moves; ++b) {
        const float c = math::dot(model, bands[b].center);
        moves = c <= bands[b].cos_near && c >= bands[b].cos_far;
      }
      moving[i] = moves;
    }
  }

  // plan_motion()'s cap radii are freed before the scan's face scratch and
  // face selection take their place.
  static_assert(MAX_BUILD_FACES * sizeof(float) <=
                MAX_BUILD_FACES + sizeof(SDF::FaceScratchBuffer));

  /**
   * @brief Whether scratch_arena_b holds a build draw's vertex deformation beside
   *        the rest of its draw.
   * @param vertices Mesh vertex count.
   * @param faces Mesh face count.
   * @param copies_model Whether the scan keeps a copy of the model.
   * @return True when the motion tables, the ripple tables, the model copy and
   * the draw's column basis, face scratch and pads fit.
   */
  HS_FLASH_MEMBER bool build_motion_fits(size_t vertices, size_t faces,
                                         bool copies_model) const {
    const size_t need =
        (copies_model ? vertices * sizeof(math::Vector) : 0) +
        Scan::mesh_motion_bytes(vertices, faces) +
        (2 * static_cast<size_t>(ripple_gen.active_count()) + 2) *
            sizeof(Animation::RippleParams) +
        sizeof(Scan::ColumnBasis<W, H>) + sizeof(SDF::FaceScratchBuffer) +
        H * sizeof(float);
    return scratch_arena_b.get_capacity() - scratch_arena_b.get_offset() >=
           need;
  }

  /** @brief Vertex deformation state of one mesh draw. */
  struct MotionPlan {
    bool active = false;                       ///< Whether vertices deform.
    uint8_t *moving = nullptr;                 ///< Nonzero per moving vertex.
    uint16_t *columns = nullptr;               ///< Per-face seeing columns.
    Animation::RippleParams *during = nullptr; ///< One per active ripple.
  };

  /**
   * @brief Fills one params per active ripple at a fraction of its last step.
   * @param during Receives the params.
   * @param u Frame fraction in [0, 1].
   */
  HS_FLASH_MEMBER void ripples_at(Animation::RippleParams *during,
                                  float u) const {
    for (int k = 0; k < ripple_gen.active_count(); ++k)
      during[k] = ripple_gen.active_params(k).during_step(u);
  }

  /**
   * @brief Plans the vertex deformation of one mesh draw.
   * @param canvas Canvas with the frame's column-time state installed.
   * @param mesh The mesh, placed at the current endpoint.
   * @param path The draw's model path.
   * @param bound motion_bound() of the draw.
   * @param arena Arena backing the plan's tables for the caller's scope.
   * @return Inactive with temporal sampling off or when every visible vertex holds.
   * @details The faces the clip can see with a moving vertex are picked; the
   * largest vertex motion activates deformation. A picked face
   * whose vertices all move at most Scan::HELD_VERTEX_ARC draws once. A
   * placed vertex strays from the current endpoint by at most its motion
   * plus the bend of its path between the endpoints, which a column of arc
   * covers.
   */
  HS_FLASH_MEMBER MotionPlan plan_motion(const Canvas &canvas,
                                         const MeshState &mesh,
                                         const ModelPath &path, float bound,
                                         Arena &arena) {
    MotionPlan plan;
    if (!Scan::mesh_motion_active<W>(bound))
      return plan;
    const size_t count = mesh.vertices.size();
    plan.moving = arena.allocate_n<uint8_t>(count);
    flag_moving(path, count, bound - ripple_arc(), plan.moving, arena);
    plan.columns = arena.allocate_n<uint16_t>(2 * mesh.num_faces());
    plan.during = arena.allocate_n<Animation::RippleParams>(
        ripple_gen.active_count() + 1);
    ScratchScope radii_scope(arena);
    float *radii = arena.allocate_n<float>(mesh.num_faces());
    if (Scan::pick_moving_faces<W>(canvas, mesh, plan.moving, bound,
                                   plan.columns, radii) == 0)
      return plan;
    ripples_at(plan.during, 0.0f);
    OrientTransformer camera(orientation);
    const float held_cos = cosf(Scan::HELD_VERTEX_ARC<W>);
    const auto start = path.moment(0.0f);
    float min_cos = 1.0f;
    for (size_t i = 0; i < count; ++i) {
      if (!plan.moving[i])
        continue;
      const float c = math::dot(
          camera(ripple_gen.transform(path.at(i, start), plan.during)),
          mesh.vertices[i]);
      plan.moving[i] = c < held_cos;
      min_cos = fminf(min_cos, c);
    }
    const float arc = acosf(hs::clamp(min_cos, -1.0f, 1.0f));
    if (Scan::drop_held_faces(mesh, plan.moving, plan.columns) == 0)
      return plan;
    plan.active = Scan::mesh_motion_active<W>(arc);
    if (plan.active)
      Scan::moving_face_columns<W>(canvas, mesh, arc + math::TWO_PI_F / W,
                                   plan.columns, radii);
    for (int k = 0; k < ripple_gen.active_count(); ++k) {
      plan.during[k].cos_threshold_min =
          fmaxf(plan.during[k].cos_threshold_min,
                ripple_gen.active_params(k).cos_threshold_min);
      plan.during[k].cos_threshold_max =
          fminf(plan.during[k].cos_threshold_max,
                ripple_gen.active_params(k).cos_threshold_max);
    }
    return plan;
  }

  /**
   * @brief Evaluates one vertex on its model and ripple paths in camera space.
   * @param plan Frame's ripple cache and motion selection.
   * @param path Model-space vertex path.
   * @param i Vertex index.
   * @param u Frame fraction.
   * @return Vertex at the requested frame fraction.
   */
  HS_HOT_INLINE math::Vector capture_position(const MotionPlan &plan,
                                              const ModelPath &path, size_t i,
                                              float u) const {
    for (int k = 0; k < ripple_gen.active_count(); ++k) {
      const auto &p = ripple_gen.active_params(k);
      plan.during[k].phase =
          p.previous_phase + (p.phase - p.previous_phase) * u;
      plan.during[k].amplitude =
          p.previous_amplitude + (p.amplitude - p.previous_amplitude) * u;
    }
    return OrientTransformer(orientation)(
        ripple_gen.transform(path.at(i, u), plan.during));
  }

  /**
   * @brief Sprite draw callback: draws the held shape for one envelope frame.
   * @param canvas Render target.
   * @param phase Sprite envelope phase: rises over the incoming window, holds 1,
   *        falls over the outgoing window.
   * @param back Carousel slot the shape was spawned into.
   * @details Draws nothing during the build window.
   */
  HS_COLD_MEMBER void draw_sprite(Canvas &canvas, float phase, int back) {
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
    const float previous =
        phase_drawn[slot] && Temporal::temporal_sampling_enabled()
            ? drawn_phase[slot]
            : phase;
    drawn_phase[slot] = phase;
    phase_drawn[slot] = true;
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
  HS_O3_FN HS_NOINLINE_NOCLONE void draw_shape(Canvas &canvas, float phase,
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
    const ModelPath path{base_state.vertices.data()};
    const float bound = motion_bound(path, base_state.vertices.size());
    const Temporal::ColumnTime column_time = Scan::prepare_column_basis<W, H>(
        scratch_arena_a, orientation,
        previous_phase != phase || Scan::mesh_motion_active<W>(bound));
    Scan::ColumnTimeScope column_scope(canvas, column_time);
    MeshState transformed_state = transform_shape(base_state);
    const MotionPlan plan =
        plan_motion(canvas, transformed_state, path, bound, scratch_arena_a);
    ArenaVector<float> face_phases;
    ArenaVector<const BakedPalette *> face_palettes;
    {
      HS_PROFILE(is_face_offsets);
      const size_t faces = transformed_state.num_faces();
      face_phases.bind(scratch_arena_a, faces);
      face_palettes.bind(scratch_arena_a, faces);
      carousel.fill_face_phases(base_state, transformed_state, phase,
                                face_phases);
      for (size_t face = 0; face < faces; ++face)
        face_palettes.push_back(&palette_bank[face_palette[face]].view());
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
        fragment_shader.set_palette(face_palettes[fi]);
        const float opacity = seg.opacity(face_phases[fi]);
        fragment_shader.set_alpha(
            opacity, previous_opacity ? previous_opacity[fi] : opacity);
        fragment_shader.scale = size > math::TOLERANCE ? 1.0f / size : 0.0f;
      };
      Scan::deform_mesh<W, H>(
          canvas, transformed_state, plan.active, plan.columns, plan.moving,
          scratch_arena_a, bound,
          [&](size_t i, float u)
              HS_HOT_INLINE { return capture_position(plan, path, i, u); },
          [&]() HS_HOT_INLINE {
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
   *        mesh or a ripple moves, the mesh at their capture times; a frame
   *        whose motion scratch_arena_b cannot hold draws at its current
   *        endpoint.
   * @param canvas Render target receiving the rasterized mesh.
   * @param mesh Compiled swept mesh (scratch-backed, this frame only).
   * @param sh Per-face palette table and previous vertices from the OpLeg.
   */
  HS_O3_FN void draw_build_mesh(Canvas &canvas, MeshState &mesh,
                                const Animation::OpLeg::Shading &sh) {
    if (mesh.vertices.is_empty())
      return;
    HS_PROFILE(is_build_draw);
    // Opened after the OpLeg's blended ramps, which the shader reads.
    ScratchScope b_guard(scratch_arena_b);
    const size_t count = mesh.vertices.size();
    ModelPath path{mesh.vertices.data(), sh.previous_vertices,
                   sh.leg ? &sh : nullptr};
    const float bound = motion_bound(path, count);
    const bool active =
        Scan::mesh_motion_active<W>(bound) &&
        build_motion_fits(count, mesh.num_faces(), path.leg == nullptr);
    // A leg moving vertex by vertex places from its own path; otherwise the
    // The scan keeps the current endpoint's model-space vertices.
    if (active && !path.leg) {
      math::Vector *model = scratch_arena_b.allocate_n<math::Vector>(count);
      std::copy_n(mesh.vertices.data(), count, model);
      path.current = model;
    }
    const Temporal::ColumnTime column_time = Scan::prepare_column_basis<W, H>(
        scratch_arena_b, orientation,
        sh.previous_weight != sh.weight || active);
    Scan::ColumnTimeScope column_scope(canvas, column_time);
    // The frame-local mesh and its source occupy scratch_a.
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(is_mesh_transform);
      MeshOps::transform_in_place(mesh, ripple_gen, camera);
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
   * @param bound motion_bound() of the draw.
   */
  HS_FLASH_MEMBER void scan_build_mesh(Canvas &canvas, MeshState &mesh,
                                       const Animation::OpLeg::Shading &sh,
                                       const ModelPath *path, float bound) {
    const MotionPlan plan =
        path ? plan_motion(canvas, mesh, *path, bound, scratch_arena_b)
             : MotionPlan{};
    const SegueT &seg = carousel.segue();
    // Per-face scan scratch must remain disjoint from the mesh in scratch_a.
    Scan::deform_mesh<W, H>(
        canvas, mesh, plan.active, plan.columns, plan.moving, scratch_arena_b,
        bound,
        [&](size_t i, float u)
            HS_HOT_INLINE { return capture_position(plan, *path, i, u); },
        [&]() HS_HOT_INLINE {
          Scan::Mesh::draw_opleg_shading<W, H>(filters, canvas, mesh, sh, 1.0f,
                                               seg.opacity(1.0f),
                                               scratch_arena_b);
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
  HS_COLD_MEMBER void
  preset_changed(const Effect::PresetChange &change) override {
    if (change.origin != Effect::PresetChangeOrigin::AUTOMATIC)
      restart_shape();
  }

  /** @brief Adds the camera's random walk to the timeline. */
  HS_COLD_MEMBER void add_orientation_walk() {
    timeline.add(0, Animation::RandomWalk<W>(orientation, math::UP, noise));
  }

  /**
   * @brief Cuts the resident shape and spawns the committed preset.
   * @details Clears every scheduled event, including any in-flight build leg
   * and ripple; the orientation walk resumes from the held orientation.
   */
  HS_COLD_MEMBER void restart_shape() {
    timeline.clear();
    abandon_build();
    add_orientation_walk();
    spawn_shape();
  }

  /** @brief Spawns the committed preset's solid. */
  HS_COLD_MEMBER void spawn_shape() {
    spawn_entry(*PRESETS[preset_index].params);
  }

  /** @brief Advances to the next preset and spawns it. */
  HS_COLD_MEMBER void advance_shape() {
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
  HS_COLD_MEMBER void resplit_for_spawn(bool has_recipe) {
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
  HS_COLD_MEMBER void spawn_entry(const Solids::Entry &entry) {
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
    phase_drawn[back] = false;

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
