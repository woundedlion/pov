/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once
#include "core/engine/engine.h"

// Forward declaration of the unit-test accessor (tests/test_effects.h) that
// pins warp_decay's endpoint invariants; the smoke harness only proves the
// effect renders, not that the warp relaxes to exactly 0 at t=1.
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
  /**
   * @brief Constructs the effect, seeding the palette, filters and warp state.
   * @details The Inigo Quilez cosine palette is color(t) = bias +
   *          amp*cos(2π(freq*t + phase)).
   */
  HS_COLD_MEMBER Thrusters()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}), palette(/*bias*/ {0.5f, 0.5f, 0.5f},
                              /*amp*/ {0.5f, 0.5f, 0.5f},
                              /*freq*/ {0.3f, 0.3f, 0.3f},
                              /*phase*/ {0.0f, 0.2f, 0.6f}),
        filters(Filter::Screen::AntiAlias<W, H>()), ring_vec(0.5f, 0.5f, 0.5f),
        amplitude(0), warp_phase(0), t_global(0),
        warp_anim(amplitude, [](float) { return 0.0f; }, 0, ease_linear) {}

  /**
   * @brief Registers tunable params and seeds the timeline.
   * @details Adds a persistent ring sprite plus a random timer that fires
   *          thrusters every 16-48 frames.
   */
  void init() override {
    register_param("Radius", &params.radius, 0.1f, 2.0f);
    register_param("Alpha", &params.alpha, 0.0f, 1.0f);

    ring_vec.normalize();

    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) { draw_ring(c, opacity); }, -1,
               16, ease_in_sin, 16, ease_out_sin));

    timeline.add(0, Animation::RandomTimer(
                        16, 48, [this](Canvas &) { on_fire_thruster(); }, true));
  }

  /**
   * @brief Advances one frame of the effect.
   * @details Steps the timeline and warp animation, expires dead thrusters,
   *          then advances and draws each live thruster.
   */
  void draw_frame() override {
    Canvas canvas(*this);

    // Wrap at 32, ring_fn's modulation period.
    t_global = (t_global + 1) % 32;

    timeline.step(canvas);
    if (!warp_anim.done()) {
      warp_anim.step(canvas);
    }

    // Expire finished thrusters from the front (FIFO: all share the same LIFE,
    // so the front always expires first).
    while (!thrusters.is_empty() && thrusters.front().expired()) {
      thrusters.pop_front();
    }

    // Advance and draw each live thruster; radius and opacity are pure
    // functions of `age`.
    for (size_t i = 0; i < thrusters.size(); ++i) {
      ThrusterContext &ctx = thrusters[i];
      float progress = static_cast<float>(ctx.age) / ThrusterContext::LIFE;
      float opacity = 1.0f - ease_out_expo(hs::clamp(progress, 0.0f, 1.0f));
      draw_thruster(canvas, ctx, ctx.radius_at(), opacity);
      ++ctx.age;
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
     * @details After this many frames the radius holds. Derived from `age` in
     *          draw_frame() (see radius_at()).
     */
    static constexpr int RADIUS_GROW_FRAMES = 8;
    static constexpr float RADIUS_MAX = 0.3f; /**< Peak ring radius. */

    int age = 0;               /**< Frames elapsed since (re)spawn. */
    Orientation<> orientation; /**< Orientation snapshot at spawn time. */
    Vector point;              /**< Thrust point on the unit sphere. */

    /**
     * @brief Reinitializes this slot for a freshly spawned thruster.
     * @param o Orientation snapshot to copy in.
     * @param p Thrust point on the unit sphere.
     * @details Trivially copyable: the circular buffer relocates slots by plain
     *          memberwise copy. No member holds a reference into this slot, so a
     *          copy never has to be rebound.
     */
    void reset(const Orientation<> &o, const Vector &p) {
      orientation = o;
      point = p;
      age = 0;
    }

    /**
     * @brief Linear ring radius for the frame being drawn.
     * @return Radius in unit-sphere units, in [0, RADIUS_MAX].
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

  StaticCircularBuffer<ThrusterContext, 16> thrusters; /**< Live thruster ring slots (FIFO). */

  // Test seam: reaches the private warp_decay endpoint invariants.
  friend struct ::hs_test::effects_tests::ThrustersWhiteBox;

  /**
   * @brief Warp amplitude decay curve over a fire's life, t in [0, 1].
   * @param t Normalized progress through the warp Mutation, 0 at fire, 1 at end.
   * @return Warp amplitude: exactly 0.7 at t=0, decaying to exactly 0 at t=1.
   * @details A bare 0.7*exp(-2t) bottoms out at 0.7*e^-2 ~= 0.095 and, once
   *          done() freezes it, leaves a residual wobble; the shift-and-
   *          renormalize below lands it on zero at t=1 so the ring fully relaxes
   *          between fires.
   */
  static float warp_decay(float t) {
    constexpr float FLOOR = 0.1353352832f; // expf(-2)
    return 0.7f * (expf(-2.0f * t) - FLOOR) / (1.0f - FLOOR);
  }

  /**
   * @brief Handles a fire event.
   * @details Picks a random warp phase, computes an opposed pair of thrust
   *          points on the ring, restarts the warp decay, spins the orientation
   *          about the axis derived from the thrust point, and spawns both
   *          thrusters.
   */
  HS_COLD_MEMBER void on_fire_thruster() {
    warp_phase = hs::rand_f() * 2 * PI_F;

    // Snapshot the warp state into the closure so the thrust-point geometry
    // can't depend on member-mutation order. amp is the residual amplitude
    // before warp_anim restarts below, so the thrust points sit on the ring as
    // it is currently displayed, not on this fire's about-to-start 0.7 warp.
    const float phase = warp_phase;
    const float amp = amplitude;
    const int frame = t_global;
    auto r_fn = [phase, amp, frame](float t) { return ring_fn(t, phase, amp, frame); };
    Basis basis = make_basis(Quaternion(), ring_vec);
    // Use params.radius (matching the visible ring) so the thrust pairs and the
    // derived spin axis track the ring under the Radius slider.
    Vector thrust_point =
        Plot::DistortedRing::fn_point(r_fn, basis, params.radius, phase);
    Vector thrust_opp =
        Plot::DistortedRing::fn_point(r_fn, basis, params.radius,
                                      phase + PI_F);

    warp_anim = Animation::Mutation(amplitude, warp_decay, 32, ease_linear);

    // Under a large warp the two oriented vectors can become near-parallel, so
    // their cross product collapses toward zero; fall back to a fixed axis on the
    // degenerate case (the spin axis is arbitrary when the pair is parallel).
    Vector thrust_axis = normalized_or(
        cross(orientation.orient(thrust_point), orientation.orient(ring_vec)),
        Y_AXIS);
    timeline.add(0, Animation::Rotation<W>(orientation, thrust_axis, 2 * PI_F,
                                           8 * 16, ease_out_expo));

    // spawn
    spawn_thruster(thrust_point);
    spawn_thruster(thrust_opp);
  }

  /**
   * @brief Pushes a fresh thruster at `point`, evicting the oldest if full.
   * @param point Thrust point on the unit sphere where the ring spawns.
   */
  HS_COLD_MEMBER void spawn_thruster(const Vector &point) {
    if (thrusters.is_full())
      thrusters.pop_front();
    thrusters.push_back(ThrusterContext());
    thrusters.back().reset(orientation, point);
  }

  /**
   * @brief Radial distortion of the ring at parameter t.
   * @param t Ring parameter in [0, 1) around the circumference.
   * @param phase Spatial warp phase in radians.
   * @param amp Warp amplitude (unit-sphere units).
   * @param frame Frame counter driving the slow temporal wave.
   * @return Signed radial offset, in unit-sphere units.
   * @details Product of a spatial warp wave (from `phase`) and a slow temporal
   *          wave (period 32 frames, from `frame`), scaled by `amp`. Pure in its
   *          arguments: callers pass an explicit snapshot of the warp state so
   *          the result can never depend on the order in which warp_phase /
   *          amplitude / t_global are mutated relative to the call.
   */
  static float ring_fn(float t, float phase, float amp, int frame) {
    // phase is radians; sin_wave's phase is cycles.
    return sin_wave(-1, 1, 2, phase / (2 * PI_F))(t) *
           sin_wave(-1, 1, 3, 0)(static_cast<float>(frame) / 32.0f) *
           amp;
  }

  /**
   * @brief Draws one thruster as a white ring at the given radius.
   * @param c Canvas to render into.
   * @param ctx Thruster slot supplying orientation and thrust point.
   * @param radius Ring radius in unit-sphere units.
   * @param opacity Fade factor in [0, 1]; multiplied by the global alpha param.
   */
  void draw_thruster(Canvas &c, const ThrusterContext &ctx, float radius,
                     float opacity) {
    Basis basis = make_basis(ctx.orientation.get(), ctx.point);
    auto fragment_shader = [this, opacity](const Vector &, Fragment &f) {
      f.color = Color4(CRGB(255, 255, 255));
      // Opacity drives both color and alpha (a quadratic edge falloff) for the
      // alpha-blended thruster fade.
      f.color.color = f.color.color * opacity;
      f.color.alpha = opacity * params.alpha;
    };
    Plot::Ring::draw<W, H>(filters, c, basis, radius, fragment_shader);
  }

  /**
   * @brief Draws the main warped ring, palette-shaded by fragment angle.
   * @param c Canvas to render into.
   * @param opacity Fade factor in [0, 1]; multiplied by the global alpha param.
   */
  void draw_ring(Canvas &c, float opacity) {
    Basis basis = make_basis(orientation.get(), ring_vec);

    auto fragment_shader = [this, opacity](const Vector &v, Fragment &f) {
      // v is the world-space fragment direction (the basis bakes orientation),
      // so banding is anchored to the fixed world X axis: the ring sweeps through
      // static world-space color bands.
      float angle = angle_between(X_AXIS, v);
      f.color = palette.get(angle / PI_F);
      f.color.alpha *= params.alpha * opacity;
    };

    const float phase = warp_phase;
    const float amp = amplitude;
    const int frame = t_global;
    Plot::DistortedRing::draw<W, H>(
        filters, c, basis, params.radius,
        [phase, amp, frame](float t) { return ring_fn(t, phase, amp, frame); },
        fragment_shader);
  }

  ProceduralPalette palette;                          /**< Cosine palette for ring shading. */
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters; /**< Anti-aliasing render pipeline. */

  Vector ring_vec;   /**< Unit normal of the main ring's plane. */
  float amplitude;   /**< Current warp amplitude driven by warp_anim. */
  float warp_phase;  /**< Spatial warp phase in radians, randomized per fire. */
  int t_global;      /**< Frame counter, wrapped to [0, 32) — ring_fn's modulation period; never overflows. */
  Animation::Mutation warp_anim;    /**< Restartable warp-amplitude decay animation. */

  Timeline timeline;                /**< Animation timeline for sprite/timer/spin. */
  Orientation<> orientation;        /**< Global orientation, spun by each fire. */

  /**
   * @brief User-tunable parameters exposed via register_param.
   */
  struct Params {
    float radius = 1.0f; /**< Ring radius in unit-sphere units. */
    float alpha = 0.2f;  /**< Global opacity multiplier in [0, 1]. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Thrusters)
