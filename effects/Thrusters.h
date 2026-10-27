/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file Thrusters.h
 * @brief Rotating distorted ring that periodically fires, spawning opposed
 *        expanding rings.
 */

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "core/engine/engine.h"

// Unit-test accessor for warp_decay endpoints and on_fire_thruster().
namespace hs_test {
namespace effects_tests {
struct ThrustersWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Rotating, palette-shaded distorted ring that periodically "fires".
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each fire warps the ring, spins the global orientation about a
 *          derived axis, and spawns a pair of opposed expanding thruster rings
 *          that grow and fade.
 */
template <int W, int H> class Thrusters : public Effect {
public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::Thrusters::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /**
   * @brief Constructs the effect, seeding the filters and warp state.
   */
  HS_COLD_MEMBER Thrusters()
      : Effect(W, H, pipeline_config<decltype(filters)>({.strobe = true})),
        filters(Filter::Screen::AntiAlias<W, H>()), ring_vec(0.5f, 0.5f, 0.5f),
        amplitude(0), warp_phase(0), t_global(0),
        warp_anim(
            amplitude, [](float) { return 0.0f; }, 0, math::ease_linear) {}

  /**
   * @brief Registers tunable params, bakes the ring palette, and seeds the
   *        timeline.
   * @details Adds a persistent ring sprite plus a random thruster-fire timer.
   */
  HS_COLD_MEMBER void init() override {
    // Radius 0 and 2 both collapse the ring to a point (ring_vec and its
    // antipode), so the slider stops short of each end.
    register_param("Radius", &params.radius, 0.1f, 1.9f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    // Inigo Quilez cosine palette: color(t) = bias + amp*cos(2π(freq*t + phase)).
    palette.bake(persistent_arena,
                 dot_keyed(ProceduralPalette(/*bias*/ {0.5f, 0.5f, 0.5f},
                                             /*amp*/ {0.5f, 0.5f, 0.5f},
                                             /*freq*/ {0.3f, 0.3f, 0.3f},
                                             /*phase*/ {0.0f, 0.2f, 0.6f})));

    ring_vec.normalize();

    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) { draw_ring(c, opacity); }, -1,
               {.fade_in = {16, math::ease_in_sin}}));

    timeline.add(
        0, Animation::RandomTimer({.min = 16, .max = 48, .repeat = true},
                                  [this](Canvas &) { on_fire_thruster(); }));
  }

  /**
   * @brief Advances one frame of the effect.
   * @details Steps the timeline and warp animation, expires dead thrusters,
   *          then advances and draws each live thruster.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // ring_fn completes three temporal cycles per 32 frames.
    t_global = (t_global + 1) % 32;

    {
      HS_PROFILE(th_timeline_step);
      timeline.step(canvas);
    }
    if (!warp_anim.done()) {
      HS_PROFILE(th_warp_step);
      warp_anim.step(canvas);
    }

    {
      HS_PROFILE(th_thrusters);
      // Expire finished thrusters from the front (FIFO: all share the same LIFE,
      // so the front always expires first).
      while (!thrusters.is_empty() && thrusters.front().expired()) {
        thrusters.pop_front();
      }

      for (size_t i = 0; i < thrusters.size(); ++i) {
        ThrusterContext &ctx = thrusters[i];
        float progress = static_cast<float>(ctx.age) / ThrusterContext::LIFE;
        float opacity =
            1.0f - math::ease_out_expo(hs::clamp(progress, 0.0f, 1.0f));
        draw_thruster(canvas, ctx, ctx.radius_at(), opacity);
        ctx.previous_age = ctx.age;
        ctx.has_history = true;
        ++ctx.age;
      }
    }
  }

private:
  /**
   * @brief One spawned thruster ring aged frame-by-frame to drive growth/fade.
   * @details A point on the unit sphere with its own orientation snapshot.
   */
  struct ThrusterContext {
    static constexpr int LIFE = 16; /**< Visible lifetime in frames. */
    /**
     * @brief Frames over which the ring radius grows from 0 to RADIUS_MAX.
     * @details After this many frames the radius holds.
     */
    static constexpr int RADIUS_GROW_FRAMES = 8;
    static constexpr float RADIUS_MAX = 0.3f; /**< Peak ring radius. */

    int age = 0;              /**< Frames elapsed since (re)spawn. */
    int previous_age = 0;     /**< Age at the preceding draw. */
    bool has_history = false; /**< A preceding draw exists for this spawn. */
    math::Orientation<> orientation; /**< Orientation snapshot at spawn time. */
    math::Vector point;              /**< Thrust point on the unit sphere. */

    /**
     * @brief Reinitializes this slot for a freshly spawned thruster.
     * @param o Orientation snapshot to copy in.
     * @param p Thrust point on the unit sphere.
     */
    void reset(const math::Orientation<> &o, const math::Vector &p) {
      orientation = o;
      point = p;
      age = 0;
      has_history = false;
    }

    /**
     * @brief Linear ring radius for the frame being drawn.
     * @return Radius as a fraction of a hemisphere, in [0, RADIUS_MAX].
     * @details Uses age + 1 (not age) so the first draw renders one linear step
     *          in rather than radius 0, reaching RADIUS_MAX after
     *          RADIUS_GROW_FRAMES frames.
     */
    float radius_at() const {
      float t = hs::clamp(static_cast<float>(age + 1) / RADIUS_GROW_FRAMES,
                          0.0f, 1.0f);
      return RADIUS_MAX * t;
    }

    /**
     * @brief Whether the slot has outlived LIFE and is recyclable.
     */
    bool expired() const { return age >= LIFE; }
  };

  StaticCircularBuffer<ThrusterContext, 16>
      thrusters; /**< Live thruster ring slots (FIFO). */

  // Test seam: reaches warp_decay and the fire path.
  friend struct ::hs_test::effects_tests::ThrustersWhiteBox;

  /**
   * @brief Warp amplitude decay curve over a fire's life, t in [0, 1].
   * @param t Normalized progress through the warp Mutation, 0 at fire, 1 at end.
   * @return Warp amplitude: exactly 0.7 at t=0, decaying to exactly 0 at t=1.
   * @details 0.7*exp(-2t), shifted and renormalized to reach zero at t=1.
   */
  static float warp_decay(float t) {
    constexpr float FLOOR = 0.1353352832f; // expf(-2)
    return 0.7f * (expf(-2.0f * t) - FLOOR) / (1.0f - FLOOR);
  }

  /**
   * @brief Handles a fire event.
   * @details Computes an opposed pair of thrust points on the ring, picks a
   *          random warp phase, restarts the warp decay, spins the orientation
   *          about the axis derived from the thrust point, and spawns both
   *          thrusters.
   */
  HS_COLD_MEMBER void on_fire_thruster() {
    // The whole warp state is read before the fire mutates any of it, so the
    // thrust points sit on the ring as it is currently displayed rather than on
    // this fire's about-to-start 0.7 warp.
    const float phase = warp_phase;
    const float amp = amplitude;
    const int frame = t_global;

    warp_phase = hs::rand_f() * 2 * math::PI_F;
    warp_reset = true;

    auto r_fn = [phase, amp, frame](float t) {
      return ring_fn(t, phase, amp, frame);
    };
    math::Basis basis = math::make_basis(math::Quaternion(), ring_vec);
    // Use params.radius (matching the visible ring) so the thrust pairs and the
    // derived spin axis track the ring under the Radius slider.
    math::Vector thrust_point =
        Plot::DistortedRing::fn_point(r_fn, basis, params.radius, phase);
    math::Vector thrust_opp = Plot::DistortedRing::fn_point(
        r_fn, basis, params.radius, phase + math::PI_F);

    warp_anim =
        Animation::Mutation(amplitude, warp_decay, 32, math::ease_linear);

    // The warp can carry thrust_point onto ring_vec, vanishing their cross
    // product; the spin axis is arbitrary when parallel, so fall back.
    math::Vector thrust_axis =
        math::normalized_or(math::cross(orientation.orient(thrust_point),
                                        orientation.orient(ring_vec)),
                            math::Y_AXIS);
    timeline.add(0, Animation::Rotation<W>(orientation, thrust_axis,
                                           2 * math::PI_F, 8 * 16,
                                           math::ease_out_expo));

    spawn_thruster(thrust_point);
    spawn_thruster(thrust_opp);
  }

  /**
   * @brief Pushes a fresh thruster at `point`, evicting the oldest if full.
   * @param point Thrust point on the unit sphere where the ring spawns.
   */
  HS_COLD_MEMBER void spawn_thruster(const math::Vector &point) {
    thrusters.push_back(ThrusterContext());
    thrusters.back().reset(orientation, point);
  }

  /**
   * @brief Radial distortion of the ring at parameter t.
   * @param t Ring parameter in [0, 1) around the circumference.
   * @param phase Spatial warp phase in radians.
   * @param amp Warp amplitude in radians.
   * @param frame Frame counter driving the temporal wave.
   * @return Signed angular radius offset in radians.
   * @details Product of a spatial warp wave and a temporal wave, scaled by
   *          `amp`. Pure in its arguments.
   */
  static float ring_fn(float t, float phase, float amp, int frame) {
    // phase is radians; sin_wave's phase is cycles.
    return math::sin_wave(-1, 1, 2, phase / (2 * math::PI_F))(t) *
           math::sin_wave(-1, 1, 3, 0)(static_cast<float>(frame) / 32.0f) * amp;
  }

  /**
   * @brief Draws one thruster as a white ring at the given radius.
   * @param c Canvas to render into.
   * @param ctx Thruster slot supplying orientation and thrust point.
   * @param radius Ring radius as a fraction of a hemisphere.
   * @param opacity Fade factor in [0, 1]; multiplied by the global alpha param.
   */
  void draw_thruster(Canvas &c, const ThrusterContext &ctx, float radius,
                     float opacity) {
    // Master alpha below one slider LSB: skip rasterizing. Thrusters keep aging in
    // draw_frame(), so they reappear mid-life when alpha rises.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;
    HS_PROFILE(th_thruster_draw);
    math::Basis basis = math::make_basis(ctx.orientation.get(), ctx.point);
    const float start_age =
        ctx.has_history && ctx.age + 1 < ThrusterContext::LIFE
            ? float(ctx.previous_age)
            : float(ctx.age);
    auto age_at = [&](float t) {
      return start_age + (float(ctx.age) - start_age) * t;
    };
    auto fragment_shader = [&, this](const math::Vector &, Fragment &f,
                                     float t) {
      const float fade =
          t == 1.0f
              ? opacity
              : 1.0f - math::ease_out_expo(hs::clamp(
                           age_at(t) / ThrusterContext::LIFE, 0.0f, 1.0f));
      f.color = Color4(CRGB(255, 255, 255));
      // Opacity scales both color and alpha, so the lifetime fade is quadratic.
      f.color.color = f.color.color * fade;
      f.color.alpha = fade * params.alpha;
    };
    auto vertex = [&](size_t i, float t) {
      const float r =
          ThrusterContext::RADIUS_MAX *
          hs::clamp((age_at(t) + 1.0f) / ThrusterContext::RADIUS_GROW_FRAMES,
                    0.0f, 1.0f);
      const float polar = r * (math::PI_F * 0.5f);
      return Plot::ring_vertex<W>(basis, polar, i);
    };
    const float start_radius =
        ThrusterContext::RADIUS_MAX *
        hs::clamp((start_age + 1.0f) / ThrusterContext::RADIUS_GROW_FRAMES,
                  0.0f, 1.0f);
    Plot::GeometryCallback geometry{
        vertex, Spacetime::PathBound::rotation(fabsf(radius - start_radius) *
                                               (math::PI_F * 0.5f))};
    Plot::Ring::draw_timed<W, H>(filters, c, basis, radius, fragment_shader,
                                 geometry, start_radius);
  }

  /**
   * @brief Draws the main warped ring, palette-shaded by fragment angle.
   * @param c Canvas to render into.
   * @param opacity Fade factor in [0, 1]; multiplied by the global alpha param.
   */
  void draw_ring(Canvas &c, float opacity) {
    const RingEndpoint current{orientation.get(), amplitude, warp_phase,
                               t_global};
    RingEndpoint start = has_ring_history ? previous_ring : current;
    if (warp_reset) {
      start.amplitude = current.amplitude;
      start.phase = current.phase;
      start.frame = current.frame;
    }
    previous_ring = current;
    has_ring_history = true;
    warp_reset = false;
    // Master alpha below one slider LSB: skip rasterizing. The warp and spin still
    // advance in draw_frame(), so the ring reappears in place when alpha rises.
    if (params.alpha < MIN_VISIBLE_ALPHA)
      return;
    HS_PROFILE(th_ring_draw);
    math::Basis basis = math::make_basis(orientation.get(), ring_vec);

    auto fragment_shader = [this, opacity](const math::Vector &v, Fragment &f) {
      // v is the world-space fragment direction and is unit (the rasterizer
      // renormalizes every shaded position), so dot(world X, v) is just v.x;
      // the LUT is baked in that cos domain (dot_keyed).
      f.color = Color4(palette.get_color(dot_key(v.x)), params.alpha * opacity);
    };

    const float phase = warp_phase;
    const float amp = amplitude;
    const int frame = t_global;
    auto shift = [phase, amp, frame](float t) {
      return ring_fn(t, phase, amp, frame);
    };
    if (Spacetime::flat_surface()) {
      Plot::DistortedRing::draw<W, H>(filters, c, basis, params.radius, shift,
                                      fragment_shader);
      return;
    }
    ScratchScope pose_storage(scratch_arena_b);
    math::Orientation<2> drawn_pose(start.pose);
    drawn_pose.push(current.pose);
    Spacetime::PoseTable poses;
    poses.bind(scratch_arena_b, W);
    poses.build(drawn_pose);
    const math::Basis local_basis =
        math::make_basis(math::Quaternion(), ring_vec);
    const auto local_frame = Plot::ring_frame(local_basis, params.radius);
    auto vertex = [&](size_t i, float t) {
      const float coordinate = float(i % (W / 2)) / W;
      auto distortion = [&](float) {
        if (t == 0.0f)
          return ring_fn(coordinate, start.phase, start.amplitude, start.frame);
        if (t == 1.0f)
          return shift(coordinate);
        const float first =
            ring_fn(coordinate, start.phase, start.amplitude, start.frame);
        return first + (shift(coordinate) - first) * t;
      };
      return Plot::ring_vertex<W>(local_frame.basis,
                                  local_frame.theta_eq + distortion(0.0f), i);
    };
    Plot::GeometryCallback geometry{
        vertex, Spacetime::PathBound::rotation(fabsf(start.amplitude) +
                                               fabsf(current.amplitude))};
    const float spin = Spacetime::rotation_chord_bound(
        math::RotationMatrix(start.pose), math::RotationMatrix(current.pose));
    Plot::DistortedRing::draw_timed<W, H>(filters, c, basis, params.radius,
                                          shift, fragment_shader, geometry,
                                          &poses, spin);
  }

  /** @brief State sampled by the main ring's sprite draw. */
  struct RingEndpoint {
    math::Quaternion pose;  /**< Drawn world orientation. */
    float amplitude = 0.0f; /**< Drawn warp amplitude. */
    float phase = 0.0f;     /**< Drawn spatial phase, radians. */
    int frame = 0;          /**< Drawn modulation frame. */
  };
  RingEndpoint previous_ring;    /**< Last sprite draw's endpoint. */
  bool has_ring_history = false; /**< The sprite has drawn an endpoint. */
  bool warp_reset = false; /**< A fire changed the warp since the last draw. */

  BakedPaletteStorage palette; /**< Ring-shading LUT, keyed by dot(X, v). */

  // One palette LUT in the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "Thrusters persistent footprint exceeds the default partition");

  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>>
      filters; /**< Anti-aliasing render pipeline. */

  math::Vector ring_vec; /**< Unit normal of the main ring's plane. */
  float amplitude;       /**< Current warp amplitude driven by warp_anim. */
  float warp_phase; /**< Spatial warp phase in radians, randomized per fire. */
  int t_global; /**< Frame counter, wrapped to [0, 32) — ring_fn's modulation period; never overflows. */
  Animation::Mutation
      warp_anim; /**< Restartable warp-amplitude decay animation. */

  /**
   * @brief Global orientation, spun by each fire.
   * @details Declared before `timeline` so it outlives the Rotations that point
   * here, which ~Timeline clears on teardown.
   */
  math::Orientation<> orientation;
  Timeline timeline; /**< Animation timeline for sprite/timer/spin. */

  /**
   * @brief User-tunable parameters exposed via register_param.
   */
  struct Params {
    float radius = 1.0f; /**< Ring radius as a fraction of a hemisphere. */
    float alpha = 0.2f;  /**< Global opacity multiplier in [0, 1]. */
  } params;
};
