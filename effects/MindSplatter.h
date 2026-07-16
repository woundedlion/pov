/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

// Unit-test accessor for the per-emitter emit-phase / hue wrap invariants.
namespace hs_test {
namespace effects_tests {
struct MindSplatterWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Particle effect spraying from cube-vertex emitters toward
 *        octahedron-vertex attractors through a Mobius warp.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Presets cyclically lerp friction/well-strength/speed params, and
 *          the whole field is randomly re-warped on a timer.
 */
template <int W, int H> class MindSplatter : public Effect {
public:
  /**
   * @brief Constructs the effect, seeding the preset table and filters.
   */
  HS_COLD_MEMBER MindSplatter()
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(filters)::any_crosses_segments}),
        presets{PRESETS},
        filters(Filter::Screen::AntiAlias<W, H>()),
        particle_system() {}

  /**
   * @brief Registers params, builds the particle system and presets timeline,
   *        bakes the palette, and kicks off the warp scheduler.
   */
  void init() override {
    static constexpr size_t SCRATCH_BYTES = 6 * 1024;
    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_BYTES, SCRATCH_BYTES, 0);

    // Compile-time device-budget guard: GLOBAL_ARENA_SIZE is inflated on the
    // host test build, so check against the real device arena literal.
    static constexpr size_t DEVICE_ARENA_BYTES = DEVICE_GLOBAL_ARENA_SIZE;
    static constexpr size_t POOL_BYTES =
        sizeof(Animation::Particle<TRAIL_LEN>) * NUM_PARTICLES;
    static constexpr size_t AUX_RESERVE_BYTES = 6 * 1024;
    static_assert(POOL_BYTES + AUX_RESERVE_BYTES <=
                      DEVICE_ARENA_BYTES - SCRATCH_BYTES,
                  "MindSplatter particle pool + palette/aux overflow the device "
                  "persistent arena");

    register_animated_param("Friction", &params.friction, 0.5f, 1.0f);
    register_animated_param("Well Str", &params.well_strength, 0.0f, 20.0f);
    register_animated_param("Init Spd", &params.initial_speed, 0.0f, 0.1f);
    register_animated_param("Ang Spd", &params.angular_speed, 0.0f, 1.0f);
    register_readonly_param("Particles", &params.active_count, 0.0f,
                          (float)NUM_PARTICLES);

    timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise));

    auto preset_timer = Animation::PeriodicTimer(
        160,
        [this](Canvas &) {
          if (animations_paused())
            return;
          presets.next();
          timeline.add(0, Animation::Lerp(params, presets.prev_get(),
                                          presets.get(), 48, ease_linear,
                                          &anims_paused_));
        },
        true);
    timeline.add(0, preset_timer);

    rebuild();
    baked_palette_.bake(persistent_arena, base_palette);
    start_warp();
  }

  /**
   * @brief Steps the timeline, pushes live params into the particle system,
   *        advances it, then renders the particles.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(msp_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(msp_timeline_step);
      timeline.step(canvas);
    }

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
  // Test seam: reaches the per-emitter emit-phase / hue wrap invariants.
  friend struct ::hs_test::effects_tests::MindSplatterWhiteBox;

  /** @brief Per-particle trail length (feeds the pool footprint below). */
  static constexpr int TRAIL_LEN = 23;

  /**
   * @brief Fixed particle pool capacity.
   * @details Footprint is host/device-identical (fixed-width trail storage:
   *          snorm16 entries, uint32_t ring indices): ~180 B/particle × 1024 =
   *          180 KiB of the 330 KiB device arena. init()'s static_assert
   *          enforces the budget at compile time against the real device arena
   *          literal.
   */
  static const int NUM_PARTICLES = 1024;

  typedef Solids::Cube EmitSolid;
  typedef Solids::Octahedron AttractSolid;

  typedef Animation::ParticleSystem<W, NUM_PARTICLES, TRAIL_LEN,
                                    EmitSolid::NUM_VERTS, AttractSolid::NUM_VERTS>
      ParticleSystem;

  /**
   * @brief Animated effect parameters snapshot.
   * @details active_count is engine-written (read-only); the remaining fields
   *          are driven by the preset Lerp or by user input when paused.
   */
  struct Params {
    float friction = 0.85f;       /**< Velocity retention per step in [0.5, 1]. */
    float well_strength = 1.0f;   /**< Attractor pull strength in [0, 20]. */
    float initial_speed = 0.025f; /**< Spawn speed in [0, 0.1] (units/step). */
    float angular_speed = 0.2f;   /**< Emission phase rate in [0, 1] (rad/emit). */
    float active_count = 0.0f;    /**< Live particle count (engine-written). */

    /**
     * @brief Linearly interpolates each field between two preset snapshots.
     * @param start Source snapshot (interpolation parameter t = 0).
     * @param target Destination snapshot (interpolation parameter t = 1).
     * @param t Interpolation factor in [0, 1].
     */
    void lerp(const Params &start, const Params &target, float t) {
      // Trips if the field set changes, so a new preset float can't silently go
      // un-interpolated (engine-written active_count is excluded on purpose).
      static_assert(sizeof(Params) == 5 * sizeof(float),
                    "MindSplatter::Params field set changed — update lerp");
      friction = start.friction + (target.friction - start.friction) * t;
      well_strength = start.well_strength +
                      (target.well_strength - start.well_strength) * t;
      initial_speed = start.initial_speed +
                      (target.initial_speed - start.initial_speed) * t;
      angular_speed = start.angular_speed +
                      (target.angular_speed - start.angular_speed) * t;
    }
  } params;

  static constexpr float FRICTION_MIN = 0.5f, FRICTION_MAX = 1.0f;
  static constexpr float WELL_STRENGTH_MIN = 0.0f, WELL_STRENGTH_MAX = 20.0f;
  static constexpr float INITIAL_SPEED_MIN = 0.0f, INITIAL_SPEED_MAX = 0.1f;
  static constexpr float ANGULAR_SPEED_MIN = 0.0f, ANGULAR_SPEED_MAX = 1.0f;

  /** @brief True iff every preset-driven field of @p p lies within its
   *  registered slider range (see the range constants above). */
  static constexpr bool preset_in_ranges(const Params &p) {
    return p.friction >= FRICTION_MIN && p.friction <= FRICTION_MAX &&
           p.well_strength >= WELL_STRENGTH_MIN &&
           p.well_strength <= WELL_STRENGTH_MAX &&
           p.initial_speed >= INITIAL_SPEED_MIN &&
           p.initial_speed <= INITIAL_SPEED_MAX &&
           p.angular_speed >= ANGULAR_SPEED_MIN &&
           p.angular_speed <= ANGULAR_SPEED_MAX;
  }

  // well_strength is pre-scaled by friction (0.85) because the integrator
  // applies v <- friction*v + impulse, dragging velocity before the attractor
  // impulse.
  static constexpr std::array<PresetEntry<Params>, 5> PRESETS{{
      {{.friction = 0.85f, .well_strength = 0.85f, .initial_speed = 0.025f, .angular_speed = 0.2f}},
      {{.friction = 0.85f, .well_strength = 0.85f, .initial_speed = 0.025f, .angular_speed = 0.52f}},
      {{.friction = 0.85f, .well_strength = 0.85f, .initial_speed = 0.025f, .angular_speed = 0.92f}},
      {{.friction = 0.85f, .well_strength = 1.7f, .initial_speed = 0.094f, .angular_speed = 0.2f}},
      {{.friction = 0.85f, .well_strength = 2.55f, .initial_speed = 0.035f, .angular_speed = 1.0f}},
  }};
  static_assert(preset_in_ranges(PRESETS[0].params) &&
                    preset_in_ranges(PRESETS[1].params) &&
                    preset_in_ranges(PRESETS[2].params) &&
                    preset_in_ranges(PRESETS[3].params) &&
                    preset_in_ranges(PRESETS[4].params),
                "a MindSplatter preset drives a param outside its registered "
                "slider range; widen the range to accommodate the preset (the "
                "range exposes the presets, it does not clamp them)");

  Presets<Params, 5> presets;

  Orientation<> orientation;
  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
  ParticleSystem particle_system;
  GenerativePalette base_palette{
      GradientShape::CIRCULAR, HarmonyType::COMPLEMENTARY,
      BrightnessProfile::FLAT, SaturationProfile::MID};
  BakedPalette baked_palette_;
  std::array<float, EmitSolid::NUM_VERTS> emitter_hues;
  /**
   * @brief Per-emitter accumulated emission angle (radians, wrapped to
   *        [0, 2pi)).
   * @details Driven by integrating Ang Spd each emission, so a live speed
   *          change alters the emission rate going forward only.
   */
  std::array<float, EmitSolid::NUM_VERTS> emit_phases;
  /**
   * @brief Per-emitter tangent-plane basis, built once in init().
   * @details Each basis depends only on its fixed emitter axis, so it is cached
   *          here rather than rebuilt inside the per-frame emitter callback. The
   *          callback is stored in a 32-byte EmitterFn, too small to also
   *          capture a 36-byte Basis, so the array is indexed by the captured i.
   */
  std::array<Basis, EmitSolid::NUM_VERTS> emitter_basis_;

  MobiusParams mobius;      /**< Current Mobius warp parameters. */
  float warp_scale = 0.6f;  /**< Magnitude of each warp animation. */

  /**
   * @brief (Re)builds the particle system from scratch.
   * @details Inits the pool, places attractors on the octahedron vertices,
   *          seeds emitter hues/phases, and installs an emitter at each cube
   *          vertex.
   */
  HS_COLD_MEMBER void rebuild() {
    particle_system.init(persistent_arena, params.friction, 0.001f, 160.0f);

    for (const auto &v : AttractSolid::vertices) {
      particle_system.add_attractor(v, params.well_strength, 0.003f, 0.2f);
    }
    HS_CHECK(particle_system.attractors.size() <= AttractSolid::NUM_VERTS,
             "attractor count exceeds draw_particles cos_eh capacity");

    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      emitter_hues[i] = hs::rand_f();
    }

    emit_phases.fill(0.0f);

    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      emitter_basis_[i] = make_basis(Quaternion(), EmitSolid::vertices[i]);

      particle_system.add_emitter([this, i](ParticleSystem &) {
        float angle = emit_phases[i];
        emit_phases[i] =
            fmodf(emit_phases[i] + params.angular_speed, 2.0f * PI_F);

        const Basis &basis = emitter_basis_[i];
        Vector vel = (basis.u * fast_cosf(angle) + basis.w * fast_sinf(angle)) *
                     params.initial_speed;

        emitter_hues[i] = fmodf(emitter_hues[i] + INV_PHI * 0.1f, 1.0f);

        if (particle_system.active() < particle_system.pool.capacity()) {
          uint16_t seed_u16 = static_cast<uint16_t>(emitter_hues[i] * 65535.0f);
          particle_system.spawn(EmitSolid::vertices[i], vel, seed_u16);
        }
      });
    }
  }

  /**
   * @brief Renders all particles through the Mobius warp, dimming each fragment
   *        by the attractor event-horizon kernels and coloring from the baked
   *        palette.
   * @param canvas Target canvas to draw the particle system into.
   * @param opacity Global opacity multiplier in [0, 1] applied to each fragment.
   */
  void draw_particles(Canvas &canvas, float opacity = 1.0f) {
    HS_PROFILE(msp_draw_particles);
    // cos(event_horizon) per attractor for the dot-product fast-reject below.
    std::array<float, AttractSolid::NUM_VERTS> cos_eh{};
    for (size_t i = 0; i < particle_system.attractors.size(); ++i) {
      cos_eh[i] = fast_cosf(particle_system.attractors[i].event_horizon);
    }

    // Position pass: Mobius warp + orientation (decides cullability).
    auto vertex_shader = [&](Fragment &f) {
      f.pos = mobius_transform(f.pos, mobius);
      f.pos = orientation.orient(f.pos);
    };

    // Deferred pass, run only for trails the segment cull keeps: per-attractor
    // alpha falloff from the pre-warp position.
    auto hole_shader = [&](Fragment &f, const Vector &original_pos) {
      float hole_alpha = 1.0f;
      for (size_t ai = 0; ai < particle_system.attractors.size(); ++ai) {
        const auto &attr = particle_system.attractors[ai];
        float cos_d = dot(original_pos, attr.position);
        if (cos_d < cos_eh[ai])
          continue;
        float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
        float t = d / attr.event_horizon;
        hole_alpha *= quintic_kernel(t);
      }
      f.v3 *= hole_alpha;
    };

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      float alpha = std::min(f.v0, f.v3);
      size_t p_idx = static_cast<size_t>(f.v2 + 0.5f);
      // Fragments only exist for live particles, so active() >= 1 here; the
      // clamp keeps a float-rounding overshoot in range. The active() guard
      // avoids an unsigned underflow to SIZE_MAX if that precondition is ever
      // violated.
      if (particle_system.active())
        p_idx = std::min<size_t>(p_idx, particle_system.active() - 1);

      const auto &p = particle_system.pool[p_idx];
      float seed_f = static_cast<float>(p.color_seed) / 65535.0f;
      float t_shifted = wrap_t(f.v0 + seed_f);
      Color4 c = baked_palette_.get(t_shifted);
      c.alpha = c.alpha * alpha * alpha * opacity;
      f.color = c;
    };

    HS_CHECK(particle_system.active() <= particle_system.pool.capacity(),
             "MindSplatter particle index space exceeds pool capacity");
    {
      HS_PROFILE(msp_particle_scan);
      Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                       fragment_shader, vertex_shader,
                                       hole_shader);
    }
  }

  /**
   * @brief Begins the self-rescheduling warp cycle.
   */
  void start_warp() { schedule_warp(); }

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
    auto warp = Animation::MobiusWarp(mobius, warp_scale, 160, false);
    warp.then([this]() { schedule_warp(); });
    timeline.add(0, warp);
  }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(MindSplatter)
