/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine.h"

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
  FLASHMEM MindSplatter()
      : Effect(W, H), presets{{{{{0.85f, 1.0f, 0.025f, 0.2f}},
                                {{0.85f, 1.0f, 0.025f, 0.52f}},
                                {{0.85f, 1.0f, 0.025f, 0.92f}},
                                {{0.85f, 2.0f, 0.094f, 0.2f}},
                                {{0.85f, 3.0f, 0.035f, 1.0f}}}}},
        filters(Filter::Screen::AntiAlias<W, H>()),
        particle_system() {}

  /**
   * @brief Registers params, builds the particle system and presets timeline,
   *        bakes the palette, and kicks off the warp scheduler.
   */
  void init() override {
    // The particle pool needs the bulk of the arena, so leave the persistent
    // arena large: carve only a small scratch_a (11 KiB) and no scratch_b.
    static constexpr size_t SCRATCH_BYTES = 11 * 1024;
    configure_arenas(GLOBAL_ARENA_SIZE - SCRATCH_BYTES, SCRATCH_BYTES, 0);

    registerParam("Friction", &params.friction, 0.5f, 1.0f);
    registerParam("Well Str", &params.well_strength, 0.0f, 20.0f);
    registerParam("Init Spd", &params.initial_speed, 0.0f, 0.1f);
    registerParam("Ang Spd", &params.angular_speed, 0.0f, 1.0f);
    registerParam("Particles", &params.active_count, 0.0f,
                  (float)NUM_PARTICLES);
    // The preset Lerp drives these four every cycle; flag them animated so the
    // "Pause Animation" toggle governs them (touching one pauses the presets
    // and hands the value to the user).
    markAnimated("Friction");
    markAnimated("Well Str");
    markAnimated("Init Spd");
    markAnimated("Ang Spd");
    // Particles is an engine-written active count, not an input.
    markReadonly("Particles");

    timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS, noise));

    auto preset_timer = Animation::PeriodicTimer(
        160,
        [this](Canvas &) {
          if (animationsPaused())
            return;
          presets.next();
          timeline.add(0, Animation::Lerp(params, presets.prev_get(),
                                          presets.get(), 48, ease_mid,
                                          &anims_paused_));
        },
        true);
    timeline.add(0, preset_timer);

    rebuild();
    baked_palette_.bake(persistent_arena, base_palette);
    start_warp();
  }

  /**
   * @brief Reports whether the engine should clear to a background each frame.
   * @return Always false; this effect manages its own canvas contents.
   */
  bool show_bg() const override { return false; }

  /**
   * @brief Steps the timeline, pushes live params into the particle system,
   *        advances it, then renders the particles.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);

    particle_system.friction = params.friction;
    // Apply Well Str live (attractors otherwise captured it once at rebuild()).
    for (size_t i = 0; i < particle_system.attractors.size(); ++i)
      particle_system.attractors[i].strength = params.well_strength;
    particle_system.step(canvas);
    params.active_count = (float)particle_system.active_count;

    draw_particles(canvas, 1.0f);
  }

private:
  /**
   * @brief Fixed particle pool capacity.
   * @details Pool footprint is identical on the 32-bit device and the 64-bit
   *          native build, since VectorTrail's circular-buffer indices are
   *          uint32_t. Particle<23> = pos(12) + vel(12) + seed(2) + life(2) +
   *          VectorTrail<23>(288) = 316 B; pool = 1024 * 316 = 316 KB. That
   *          leaves ~19 KB under the 335 KB arena to cover the 11 KB scratch
   *          carve plus the baked palette and attractors/emitters, so the
   *          device fits with headroom.
   */
  static const int NUM_PARTICLES = 1024;

  typedef Solids::Cube EmitSolid;
  typedef Solids::Octahedron AttractSolid;

  typedef Animation::ParticleSystem<W, NUM_PARTICLES, 23, EmitSolid::NUM_VERTS,
                                    AttractSolid::NUM_VERTS>
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
      friction = start.friction + (target.friction - start.friction) * t;
      well_strength = start.well_strength +
                      (target.well_strength - start.well_strength) * t;
      initial_speed = start.initial_speed +
                      (target.initial_speed - start.initial_speed) * t;
      angular_speed = start.angular_speed +
                      (target.angular_speed - start.angular_speed) * t;
    }
  } params;

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

  MobiusParams mobius;      /**< Current Mobius warp parameters. */
  float warp_scale = 0.6f;  /**< Magnitude of each warp animation. */

  /**
   * @brief (Re)builds the particle system from scratch.
   * @details Inits the pool, places attractors on the octahedron vertices,
   *          seeds emitter hues/phases, and installs an emitter at each cube
   *          vertex.
   */
  FLASHMEM void rebuild() {
    particle_system.init(persistent_arena, params.friction, 0.001f, 160.0f);

    for (const auto &v : AttractSolid::vertices) {
      particle_system.add_attractor(v, params.well_strength, 0.003f, 0.2f);
    }

    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      emitter_hues[i] = hs::rand_f();
    }

    emit_phases.fill(0.0f);

    for (size_t i = 0; i < EmitSolid::NUM_VERTS; ++i) {
      Vector axis = EmitSolid::vertices[i];

      particle_system.add_emitter([this, i, axis](ParticleSystem &) {
        // Integrate Ang Spd into a wrapped per-emitter phase so a live speed
        // change affects only the rate going forward; the [0,2pi) wrap keeps
        // cos/sin precise over long runs.
        float angle = emit_phases[i];
        emit_phases[i] =
            fmodf(emit_phases[i] + params.angular_speed, 2.0f * PI_F);

        // Velocity in the emitter's tangent plane, rotated by the phase angle.
        auto basis = make_basis(Quaternion(), axis);
        Vector vel = (basis.u * cosf(angle) + basis.w * sinf(angle)) *
                     params.initial_speed;

        // Advance hue by the golden ratio for an even spread across emissions.
        emitter_hues[i] = fmodf(emitter_hues[i] + G * 0.1f, 1.0f);

        if (particle_system.active_count < particle_system.pool.capacity()) {
          uint16_t seed_u16 = static_cast<uint16_t>(emitter_hues[i] * 65535.0f);
          particle_system.spawn(axis, vel, seed_u16);
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
    // Precompute cos(event_horizon) per attractor for a dot-product fast-reject:
    // points with dot < this threshold lie beyond the horizon (no influence)
    // and are skipped below.
    std::array<float, AttractSolid::NUM_VERTS> cos_eh;
    for (size_t i = 0; i < particle_system.attractors.size(); ++i) {
      cos_eh[i] = fast_cosf(particle_system.attractors[i].event_horizon);
    }

    // Accumulate per-attractor alpha falloff from the pre-warp position, then
    // apply the Mobius warp and orientation to the fragment position.
    //
    // This inner loop runs per emitted fragment, so the shader is O(fragments x
    // attractors). That is accepted, not overlooked: AttractSolid::NUM_VERTS is a
    // small compile-time constant (6 for the octahedron) and the cos_eh dot-
    // product reject above skips most of those iterations, so the constant stays
    // tiny — no spatial cull is warranted at this attractor count. If the
    // attractor solid ever grows, bucket the pre-warp position (already in hand
    // here, before the warp) to cull distant attractors before the loop.
    auto vertex_shader = [&](Fragment &f) {
      Vector original_pos = f.pos;
      float holeAlpha = 1.0f;
      for (size_t ai = 0; ai < particle_system.attractors.size(); ++ai) {
        const auto &attr = particle_system.attractors[ai];
        float cos_d = dot(original_pos, attr.position);
        if (cos_d < cos_eh[ai])
          continue;
        float d = fast_acos(hs::clamp(cos_d, -1.0f, 1.0f));
        float t = d / attr.event_horizon;
        holeAlpha *= quintic_kernel(t);
      }

      f.pos = mobius_transform(f.pos, mobius);
      f.pos = orientation.orient(f.pos);
      f.v3 *= holeAlpha;
    };

    // Color the fragment from the per-particle hue seed (skipping inactive
    // particles), squaring alpha for a softer falloff.
    auto fragment_shader = [&](const Vector &, Fragment &f) {
      float alpha = std::min(f.v0, f.v3);
      size_t p_idx = static_cast<size_t>(f.v2 + 0.5f);

      if (p_idx >= particle_system.active_count) {
        f.color = Color4(CRGB(0, 0, 0), 0.0f);
        return;
      }

      const auto &p = particle_system.pool[p_idx];
      float seed_f = static_cast<float>(p.color_seed) / 65535.0f;
      float t_shifted = wrap_t(f.v0 + seed_f);
      Color4 c = baked_palette_.get(t_shifted);
      c.alpha = c.alpha * alpha * alpha * opacity;
      f.color = c;
    };

    Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system,
                                     fragment_shader, vertex_shader);
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

#include "core/effect_registry.h"
REGISTER_EFFECT(MindSplatter)
