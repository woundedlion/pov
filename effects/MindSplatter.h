/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file MindSplatter.h
 * @brief Particles sprayed from Platonic-solid emitters toward dual-solid
 *        attractors through a Mobius warp.
 */

#include "core/control/choreography.h"
#include "core/engine/engine.h"
// 256 x 256 Pixels = 393,216 B of flash, about a fifth of the Teensy budget;
// the trade is stated in tools/mindsplatter_palette_gen.cpp.
#include "core/color/triadic_palette_luts.h"

// Unit-test accessor for emitter and hole-kernel invariants.
namespace hs_test {
namespace effects_tests {
struct MindSplatterWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Animated MindSplatter parameters snapshot.
 * @details active_count is engine-written (read-only); the remaining fields
 *          are driven by the preset Lerp or by user input when paused.
 */
struct MindSplatterParams {
  Solids::BaseMesh base_mesh = Solids::BaseMesh::CUBE;
  float friction = 0.85f;       /**< Velocity retention per step in [0.5, 1]. */
  float well_strength = 1.0f;   /**< Attractor pull strength in [0, 20]. */
  float initial_speed = 0.025f; /**< Spawn speed in [0, 0.5] (units/step). */
  float angular_speed = 0.2f; /**< Emission phase rate in [0, 1] (rad/emit). */
  float warp_scale = 0.6f;    /**< Mobius warp magnitude in [0, 5]. */
  float active_count = 0.0f;  /**< Live particle count (engine-written). */

  /**
   * @brief Linearly interpolates each field between two preset snapshots.
   * @param start Source snapshot (interpolation parameter t = 0).
   * @param target Destination snapshot (interpolation parameter t = 1).
   * @param t Interpolation factor in [0, 1].
   */
  void lerp(const MindSplatterParams &start, const MindSplatterParams &target,
            float t) {
    // Trips if the field set changes, so a new preset float can't silently go
    // un-interpolated (engine-written active_count is excluded on purpose).
    static_assert(sizeof(MindSplatterParams) == 7 * sizeof(float),
                  "MindSplatter::Params field set changed — update lerp");
    base_mesh = t < 0.5f ? start.base_mesh : target.base_mesh;
    friction = start.friction + (target.friction - start.friction) * t;
    well_strength =
        start.well_strength + (target.well_strength - start.well_strength) * t;
    initial_speed =
        start.initial_speed + (target.initial_speed - start.initial_speed) * t;
    angular_speed =
        start.angular_speed + (target.angular_speed - start.angular_speed) * t;
    warp_scale = start.warp_scale + (target.warp_scale - start.warp_scale) * t;
  }
};

/**
 * @brief Particle effect spraying from Platonic-solid emitters toward
 *        dual-solid attractors through a Mobius warp.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Presets cyclically lerp friction/well-strength/speed params, and
 *          the whole field is randomly re-warped on a timer.
 */
template <int W, int H>
class MindSplatter
    : public ChoreographedEffect<MindSplatter<W, H>, MindSplatterParams> {
  using Choreography =
      ChoreographedEffect<MindSplatter<W, H>, MindSplatterParams>;
  friend Choreography;

public:
  using BaseMesh = Solids::BaseMesh;
  using Params = MindSplatterParams;

  /** Preset policy: an automatic change crossfades the live parameters over
      48 frames; pause freezes an in-flight crossfade. */
  static constexpr Segue::Lerp PRESET_SEGUE{48, ease_linear, /*pausable=*/true};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  /** Dwell + blend = the 160-frame preset cadence. */
  static constexpr uint16_t PRESET_DWELL_FRAMES = 112;

  /** @brief Initial live parameters: the defaults, not a preset entry. */
  static Params initial_params() { return Params{}; }

  static constexpr size_t authored_preset_count() { return PRESETS.size(); }

  /**
   * @brief Constructs the effect, seeding the filters.
   */
  HS_COLD_MEMBER MindSplatter()
      : Choreography(W, H,
                     pipeline_config<decltype(filters)>({.strobe = true})),
        particle_system() {}

  /**
   * @brief Registers params, builds the particle system, bakes the palette,
   *        and kicks off the warp scheduler.
   */
  void init() override {
    configure_presets(PRESETS.size());
    static constexpr size_t SCRATCH_BYTES = 6 * 1024;
    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_BYTES, SCRATCH_BYTES, 0);

    // Compile-time device-budget guard: GLOBAL_ARENA_SIZE is inflated on the
    // host test build, so check against the real device arena literal.
    static constexpr size_t DEVICE_ARENA_BYTES = DEVICE_GLOBAL_ARENA_SIZE;
    static constexpr size_t POOL_BYTES =
        sizeof(Animation::Particle<TRAIL_LEN>) * NUM_PARTICLES;
    static constexpr size_t AUX_RESERVE_BYTES = 6 * 1024;
    static_assert(
        POOL_BYTES + AUX_RESERVE_BYTES <= DEVICE_ARENA_BYTES - SCRATCH_BYTES,
        "MindSplatter particle pool + palette/aux overflow the device "
        "persistent arena");

    // One trail at a time stages its fragments, the pre-shader positions the
    // deferred hole shader reads, and the gate arrays in scratch_a, all live
    // across rasterize. A live tip is appended past the stored history.
    static constexpr size_t TRAIL_POINTS = TRAIL_LEN + 1;
    static_assert(TRAIL_POINTS * (sizeof(Fragment) + sizeof(Vector)) +
                          Plot::rasterize_scratch_a_bytes<W>(0, TRAIL_POINTS) <=
                      SCRATCH_BYTES,
                  "MindSplatter trail staging exceeds its scratch_a split; "
                  "retune TRAIL_LEN or enlarge the split");

    register_animated_param(
        "Base Mesh", &params.base_mesh, Solids::BASE_MESH_OPTIONS,
        Solids::BASE_MESH_EXPORT_OPTIONS, Solids::PLATONIC_BASE_MESH_COUNT);
    register_animated_param("Friction", &params.friction, 0.5f, 1.0f);
    register_animated_param("Well Str", &params.well_strength, 0.0f, 20.0f);
    register_animated_param("Init Spd", &params.initial_speed,
                            INITIAL_SPEED_MIN, INITIAL_SPEED_MAX);
    register_animated_param("Ang Spd", &params.angular_speed, 0.0f, 1.0f);
    register_animated_param("Warp", &params.warp_scale, WARP_SCALE_MIN,
                            WARP_SCALE_MAX);
    register_readonly_param("Particles", &params.active_count, 0.0f,
                            (float)NUM_PARTICLES);

    timeline.add(0,
                 Animation::RandomWalk<W, 4, true>(orientation, Y_AXIS, noise));

    // First dwell spans a full cadence period, so the opening preset holds as
    // long as every later one (dwell + blend).
    hold_initial_preset(PRESET_DWELL_FRAMES + PRESET_SEGUE.frames);

    build_particle_system();
    schedule_warp();
  }

  /**
   * @brief Steps the timeline, pushes live params into the particle system,
   *        advances it, then renders the particles.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(msp_timeline_step);
      timeline.step(canvas);
    }
    begin_automatic_transition();

    if (!particle_geometry_ready || params.base_mesh != active_base_mesh)
      configure_particle_geometry(params.base_mesh);

    particle_system.friction = params.friction;
    // All attractors share the one live (preset-animated) Well Str slider; the
    // strength passed at add_attractor time is just a seed overwritten here.
    for (size_t i = 0; i < particle_system.attractors.size(); ++i)
      particle_system.attractors[i].strength = params.well_strength;
    {
      HS_PROFILE(msp_particle_step);
      particle_system.step(canvas);
    }
    params.active_count = (float)particle_system.active();

    draw_particles(canvas);
  }

private:
  using Choreography::begin_automatic_transition;
  using Choreography::configure_presets;
  using Choreography::hold_initial_preset;
  using Choreography::params;
  using Choreography::register_animated_param;
  using Choreography::register_readonly_param;
  using Choreography::timeline;
  using Choreography::transition;

  /** @brief Writes the interpolated parameters of an in-flight transition. */
  void blend_params(float progress) {
    params.lerp(transition.from, transition.to, progress);
  }

  // Test seam for emitter and attractor invariants.
  friend struct ::hs_test::effects_tests::MindSplatterWhiteBox;

  /** @brief Per-particle trail length (feeds the pool footprint below). */
  static constexpr int TRAIL_LEN = 8;

  /** @brief Frames between stored trail anchors. */
  static constexpr int TRAIL_SAMPLE_STRIDE = 3;

  /**
   * @brief Fixed particle pool capacity.
   * @details Footprint is host/device-identical (fixed-width trail storage:
   *          snorm16 entries, uint32_t ring indices): 88 B/particle × 1672 =
   *          143.7 KiB of the 298 KiB device arena. init()'s static_assert
   *          enforces the budget at compile time against the real device arena
   *          literal.
   */
  static constexpr int NUM_PARTICLES = 1672;

  static constexpr int MAX_EMITTERS = Solids::Dodecahedron::NUM_VERTS;
  static constexpr int MAX_ATTRACTORS = Solids::Dodecahedron::NUM_VERTS;

  typedef Animation::ParticleSystem<W, NUM_PARTICLES, TRAIL_LEN, MAX_EMITTERS,
                                    MAX_ATTRACTORS, true, TRAIL_SAMPLE_STRIDE>
      ParticleSystem;

  static constexpr float FRICTION_MIN = 0.5f, FRICTION_MAX = 1.0f;
  static constexpr float WELL_STRENGTH_MIN = 0.0f, WELL_STRENGTH_MAX = 20.0f;
  static constexpr float INITIAL_SPEED_MIN = 0.0f, INITIAL_SPEED_MAX = 0.5f;
  static constexpr float ANGULAR_SPEED_MIN = 0.0f, ANGULAR_SPEED_MAX = 1.0f;
  static constexpr float WARP_SCALE_MIN = 0.0f, WARP_SCALE_MAX = 5.0f;
  static constexpr float EVENT_HORIZON = 0.2f;
  static constexpr float GRAVITY = 0.001f;
  static constexpr float PARTICLE_LIFETIME_FRAMES = 160.0f;
  static constexpr float ATTRACTOR_KILL_RADIUS = 0.003f;

  template <size_t N>
  static consteval bool
  event_horizons_do_not_overlap(const std::array<Vector, N> &vertices) {
    constexpr float MIN_CHORD_SQUARED = 4.0f * EVENT_HORIZON * EVENT_HORIZON;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = i + 1; j < N; ++j) {
        const float dx = vertices[i].x - vertices[j].x;
        const float dy = vertices[i].y - vertices[j].y;
        const float dz = vertices[i].z - vertices[j].z;
        if (dx * dx + dy * dy + dz * dz <= MIN_CHORD_SQUARED)
          return false;
      }
    }
    return true;
  }
  static_assert(
      event_horizons_do_not_overlap(Solids::Tetrahedron::vertices) &&
          event_horizons_do_not_overlap(Solids::Cube::vertices) &&
          event_horizons_do_not_overlap(Solids::Octahedron::vertices) &&
          event_horizons_do_not_overlap(Solids::Dodecahedron::vertices) &&
          event_horizons_do_not_overlap(Solids::Icosahedron::vertices),
      "MindSplatter event-horizon caps must not overlap");

  static consteval bool attractors_are_signed_axes() {
    constexpr std::array<Vector, 6> AXES = {X_AXIS,  -X_AXIS, Y_AXIS,
                                            -Y_AXIS, Z_AXIS,  -Z_AXIS};
    if (Solids::Octahedron::vertices.size() != AXES.size())
      return false;
    for (size_t i = 0; i < AXES.size(); ++i) {
      const Vector &v = Solids::Octahedron::vertices[i];
      if (v.x != AXES[i].x || v.y != AXES[i].y || v.z != AXES[i].z)
        return false;
    }
    return true;
  }
  // ParticleSystem's signed-axis physics pairs attractors[2p]/[2p+1] as the
  // plus/minus poles of axis p, so the emission order is part of the contract.
  static_assert(attractors_are_signed_axes(),
                "MindSplatter attractors must be the six signed axes in "
                "+X,-X,+Y,-Y,+Z,-Z order");

  static inline float octahedral_hole_alpha(const Vector &p,
                                            float cos_event_horizon) {
    const float m =
        std::max(std::abs(p.x), std::max(std::abs(p.y), std::abs(p.z)));
    if (m < cos_event_horizon) {
      HS_MSP_COUNT(hole_early_outs);
      return 1.0f;
    }
    const float d = fast_acos(hs::clamp(m, -1.0f, 1.0f));
    return quintic_kernel(d / EVENT_HORIZON);
  }

  float attractor_hole_alpha(const Vector &p, float cos_event_horizon) const {
    float alpha = 1.0f;
    for (const auto &attractor : particle_system.attractors) {
      const float cos_distance = dot(p, attractor.position);
      if (cos_distance < cos_event_horizon)
        continue;
      const float distance = fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f));
      alpha *= quintic_kernel(distance / EVENT_HORIZON);
      if (alpha != 1.0f)
        return alpha;
    }
    return alpha;
  }

#if HS_ENABLE_TEST_ORACLES
  float reference_attractor_hole_alpha(const Vector &p,
                                       float cos_event_horizon) const {
    float alpha = 1.0f;
    for (const auto &attractor : particle_system.attractors) {
      const float cos_distance = dot(p, attractor.position);
      if (cos_distance < cos_event_horizon)
        continue;
      const float distance = fast_acos(hs::clamp(cos_distance, -1.0f, 1.0f));
      alpha *= quintic_kernel(distance / EVENT_HORIZON);
    }
    return alpha;
  }
#endif

  /** @brief True iff every preset-driven field of @p p lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Params &p) {
    return static_cast<size_t>(p.base_mesh) <
               Solids::PLATONIC_BASE_MESH_COUNT &&
           p.friction >= FRICTION_MIN && p.friction <= FRICTION_MAX &&
           p.well_strength >= WELL_STRENGTH_MIN &&
           p.well_strength <= WELL_STRENGTH_MAX &&
           p.initial_speed >= INITIAL_SPEED_MIN &&
           p.initial_speed <= INITIAL_SPEED_MAX &&
           p.angular_speed >= ANGULAR_SPEED_MIN &&
           p.angular_speed <= ANGULAR_SPEED_MAX &&
           p.warp_scale >= WARP_SCALE_MIN && p.warp_scale <= WARP_SCALE_MAX;
  }

  /** @brief Whether a parameter set is admissible for a snapshot restore. */
  static bool valid_params(const Params &p) { return preset_in_ranges(p); }

  // well_strength is pre-scaled by the preset's own friction because the
  // integrator applies v <- friction*v + impulse, dragging velocity before the
  // attractor impulse.
  static constexpr size_t PRESET_COUNT = 8;
  static constexpr std::array<PresetEntry<Params>, PRESET_COUNT> PRESETS{{
      {{.base_mesh = BaseMesh::CUBE,
        .friction = 0.85f,
        .well_strength = 0.85f,
        .initial_speed = 0.025f,
        .angular_speed = 0.2f}},
      {{.base_mesh = BaseMesh::DODECAHEDRON,
        .friction = 1.0f,
        .well_strength = 9.06f,
        .initial_speed = 0.5f,
        .angular_speed = 0.069f}},
      {{.base_mesh = BaseMesh::ICOSAHEDRON,
        .friction = 0.9645f,
        .well_strength = 16.280001f,
        .initial_speed = 0.1f,
        .angular_speed = 1.0f}},
      {{.base_mesh = BaseMesh::OCTAHEDRON,
        .friction = 0.93f,
        .well_strength = 1.74f,
        .initial_speed = 0.1f,
        .angular_speed = 1.0f}},
      {{.base_mesh = BaseMesh::TETRAHEDRON,
        .friction = 1.0f,
        .well_strength = 4.6f,
        .initial_speed = 0.5f,
        .angular_speed = 0.055f,
        .warp_scale = 0.0f}},
      {{.base_mesh = BaseMesh::DODECAHEDRON,
        .friction = 1.0f,
        .well_strength = 15.32f,
        .initial_speed = 0.5f,
        .angular_speed = 0.361f}},
      {{.base_mesh = BaseMesh::ICOSAHEDRON,
        .friction = 0.7465f,
        .well_strength = 1.54f,
        .initial_speed = 0.5f,
        .angular_speed = 0.164f}},
      {{.base_mesh = BaseMesh::CUBE,
        .friction = 1.0f,
        .well_strength = 4.6f,
        .initial_speed = 0.5f,
        .angular_speed = 0.164f,
        .warp_scale = 0.0f}},
  }};
  static_assert(all_presets_in_ranges(PRESETS, preset_in_ranges),
                "a MindSplatter preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset (the "
                "range exposes the presets, it does not clamp them)");

  /** @brief Params for the preset at @p index. */
  static Params preset_params(size_t index) { return PRESETS[index].params; }

  // orientation/noise/mobius are borrowed by timeline-resident animations.
  Orientation<> orientation;
  FastNoiseLite noise;
  MobiusParams mobius; /**< Current Mobius warp parameters. */
  Filter::Screen::DirectAntiAliasSink<W, H> filters;
  ParticleSystem particle_system;
  /**
   * @brief Per-emitter accumulated emission angle (radians, wrapped to
   *        [0, 2pi)).
   * @details Integrates Ang Spd on every emitter tick, including ticks whose
   *          spawn is dropped at pool capacity.
   */
  std::array<float, MAX_EMITTERS> emit_phases;
  uint8_t palette_sequence = 0;
  /**
   * @brief Per-emitter tangent-plane basis, built once in init().
   * @details The emitter callback is stored in a 32-byte EmitterFn, too small
   *          to also capture a 36-byte Basis, so it indexes this array by the
   *          captured i.
   */
  std::array<Basis, MAX_EMITTERS> emitter_basis;
  std::array<Vector, MAX_EMITTERS> emitter_positions;
  BaseMesh active_base_mesh = BaseMesh::CUBE;
  bool particle_geometry_ready = false;

#if HS_ENABLE_TEST_ORACLES
  bool reference_orientation = false;
  bool reference_vertex_pass = false;
  bool reference_hole_kernel = false;
#endif

  __attribute__((always_inline)) static Pixel
  sample_trail_palette(const Pixel *colors, float t) {
    const float index =
        t * static_cast<float>(MINDSPLATTER_PALETTE_LUT_SIZE - 1);
    if (index <= 0.0f)
      return colors[0];
    return lut_sample_pixel(colors, MINDSPLATTER_PALETTE_LUT_SIZE, index);
  }

  static std::span<const Vector> emitter_vertices(BaseMesh base_mesh) {
    switch (base_mesh) {
    case BaseMesh::TETRAHEDRON:
      return Solids::Tetrahedron::vertices;
    case BaseMesh::CUBE:
      return Solids::Cube::vertices;
    case BaseMesh::OCTAHEDRON:
      return Solids::Octahedron::vertices;
    case BaseMesh::DODECAHEDRON:
      return Solids::Dodecahedron::vertices;
    case BaseMesh::ICOSAHEDRON:
      return Solids::Icosahedron::vertices;
    default:
      HS_CHECK(false, "MindSplatter base mesh must be Platonic");
      return {};
    }
  }

  static std::span<const Vector> attractor_vertices(BaseMesh base_mesh) {
    switch (base_mesh) {
    case BaseMesh::TETRAHEDRON:
      return Solids::Tetrahedron::vertices;
    case BaseMesh::CUBE:
      return Solids::Octahedron::vertices;
    case BaseMesh::OCTAHEDRON:
      return Solids::Cube::vertices;
    case BaseMesh::DODECAHEDRON:
      return Solids::Icosahedron::vertices;
    case BaseMesh::ICOSAHEDRON:
      return Solids::Dodecahedron::vertices;
    default:
      HS_CHECK(false, "MindSplatter base mesh must be Platonic");
      return {};
    }
  }

  HS_FLASH_MEMBER void configure_particle_geometry(BaseMesh base_mesh) {
    const std::span<const Vector> emitters = emitter_vertices(base_mesh);
    const std::span<const Vector> attractors = attractor_vertices(base_mesh);

    particle_system.emitters.clear();
    particle_system.attractors.clear();
    particle_system.set_signed_axis_attractors(base_mesh == BaseMesh::CUBE);
    emit_phases.fill(0.0f);

    for (const Vector &v : attractors) {
      const Vector position = base_mesh == BaseMesh::TETRAHEDRON ? -v : v;
      particle_system.add_attractor(position, params.well_strength,
                                    ATTRACTOR_KILL_RADIUS, EVENT_HORIZON);
    }

    for (size_t i = 0; i < emitters.size(); ++i) {
      emitter_positions[i] = emitters[i];
      emitter_basis[i] = make_basis(Quaternion(), emitter_positions[i]);
      particle_system.add_emitter([this, i](ParticleSystem &) {
        float angle = emit_phases[i];
        emit_phases[i] =
            fmodf(emit_phases[i] + params.angular_speed, 2.0f * PI_F);

        const Basis &basis = emitter_basis[i];
        Vector vel = (basis.u * fast_cosf(angle) + basis.w * fast_sinf(angle)) *
                     params.initial_speed;

        const uint16_t color_seed = static_cast<uint16_t>(palette_sequence++)
                                    << 8;
        if (particle_system.active() < particle_system.pool.capacity()) {
          particle_system.spawn(emitter_positions[i], vel, color_seed);
        }
      });
    }

    active_base_mesh = base_mesh;
    particle_geometry_ready = true;
  }

  /**
   * @brief Builds the particle system.
   * @details Inits the pool and installs the selected Platonic emitter solid
   *          with its dual as the attractor solid. Single-shot: the arena has
   *          no per-allocation free, so ParticleSystem::init traps on a second
   *          call.
   */
  HS_COLD_MEMBER void build_particle_system() {
    particle_system.init(persistent_arena, params.friction, GRAVITY,
                         PARTICLE_LIFETIME_FRAMES);
    palette_sequence = 0;
    configure_particle_geometry(params.base_mesh);
  }

  /**
   * @brief Renders all particles through the Mobius warp, dimming each fragment
   *        by the attractor event-horizon kernels and coloring from its
   *        per-particle flash palette.
   * @param canvas Target canvas to draw the particle system into.
   * @param opacity Global opacity multiplier in [0, 1] applied to each fragment.
   */
  void draw_particles(Canvas &canvas, float opacity = 1.0f) {
    HS_PROFILE(msp_draw_particles);
    draw_particles_with(filters, canvas, opacity);
  }

  template <typename Sink>
  void draw_particles_with(Sink &sink, Canvas &canvas, float opacity = 1.0f) {
    if constexpr (requires { sink.prepare(canvas); })
      sink.prepare(canvas);

    const float cos_event_horizon = fast_cosf(EVENT_HORIZON);
    const RotationMatrix rotation(orientation.get());
    const Pixel *trail_palette = nullptr;

    // Position pass: Mobius warp + orientation (decides cullability).
    auto vertex_shader = [&](Fragment &f) {
      f.pos = mobius_transform(f.pos, mobius);
#if HS_ENABLE_TEST_ORACLES
      if (reference_orientation) {
        f.pos = orientation.orient(f.pos);
        return;
      }
#endif
      f.pos = rotation.apply(f.pos);
    };

    // Signed-axis event-horizon falloff from the pre-warp position.
    auto hole_shader = [&](FragmentRegisters f, const Vector &original_pos) {
#if HS_ENABLE_TEST_ORACLES
      if (reference_hole_kernel) {
        f.v3 *= reference_attractor_hole_alpha(original_pos, cos_event_horizon);
        return;
      }
#endif
      if (active_base_mesh == BaseMesh::CUBE)
        f.v3 *= octahedral_hole_alpha(original_pos, cos_event_horizon);
      else
        f.v3 *= attractor_hole_alpha(original_pos, cos_event_horizon);
    };

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      const float alpha = std::max(0.0f, std::min(f.v0, f.v3));
      const float palette_t = 1.0f - f.v0;
      if (f.v0 <= 0.0f || f.v0 >= 1.0f)
        HS_MSP_COUNT(palette_endpoints);
      else
        HS_MSP_COUNT(palette_interpolated);
      f.color = Color4(sample_trail_palette(trail_palette, palette_t),
                       alpha * opacity);
    };

    // v2 mapper used as a per-particle palette bind: ParticleSystem::draw calls
    // it before the particle's fragments are shaded, so fragment_shader always
    // reads this particle's palette. v2 itself is unused by the shader.
    auto prepare_trail_palette = [&](const auto &p, int) {
      trail_palette = MINDSPLATTER_PALETTES[p.color_seed >> 8];
      return 0.0f;
    };

    HS_CHECK(particle_system.active() <= particle_system.pool.capacity(),
             "MindSplatter particle index space exceeds pool capacity");
    {
      HS_PROFILE(msp_particle_scan);
#if HS_ENABLE_TEST_ORACLES
      if (reference_vertex_pass) {
        Plot::ParticleSystem::draw<W, H>(sink, canvas, particle_system,
                                         fragment_shader, vertex_shader,
                                         hole_shader, prepare_trail_palette);
        return;
      }
#endif
      using DirectSink = Filter::Screen::DirectAntiAliasSink<W, H>;
      if constexpr (std::same_as<std::remove_cvref_t<Sink>, DirectSink>) {
        Plot::ParticleSystem::draw_fused_vertex<W, H, true>(
            sink, canvas, particle_system, fragment_shader, vertex_shader,
            hole_shader, prepare_trail_palette);
      } else {
        Plot::ParticleSystem::draw_fused_vertex<W, H>(
            sink, canvas, particle_system, fragment_shader, vertex_shader,
            hole_shader, prepare_trail_palette);
      }
    }
  }

  /**
   * @brief Arms a one-shot timer (180-300 steps) that triggers the next warp.
   */
  void schedule_warp() {
    auto timer =
        Animation::RandomTimer(180, 300, [this](Canvas &) { perform_warp(); });
    timeline.add(0, timer);
  }

  /**
   * @brief Runs one Mobius warp animation, then re-arms the timer for the next.
   */
  void perform_warp() {
    auto warp = Animation::MobiusWarp(mobius, params.warp_scale, 160, false);
    warp.bind_scale(params.warp_scale);
    warp.then([this]() { schedule_warp(); });
    timeline.add(0, warp);
  }
};

#include "core/control/registry.h"
REGISTER_EFFECT(MindSplatter)
