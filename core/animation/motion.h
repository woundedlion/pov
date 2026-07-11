/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#ifndef HS_ANIMATION_INTERNAL
#error internal fragment of animation.h; include "animation.h" instead
#endif

/**
 * @brief Represents a customizable path.
 * @tparam W Display width (carried for downstream sizing).
 * @tparam RESOLUTION Capacity of the internal point ring buffer.
 * @details Retains an internal buffer for state, but draws to the pipeline.
 */
template <int W, int RESOLUTION = 1024> class Path {
public:
  /**
   * @brief Constructs an empty path.
   */
  Path() {}

  /**
   * @brief Appends a procedurally generated segment to the path.
   * @param plot The function generating points.
   * @param domain The input domain scale.
   * @param samples The number of sample intervals to take (>= 1); the loop
   *        emits samples + 1 points so the endpoint easing(1.0) is always hit.
   * @param easing The easing function for sample distribution.
   * @return Reference to self for chaining.
   */
  Path &append_segment(PlotFn plot, float domain, int samples,
                       ScalarFn easing) {
    // samples >= 1 also keeps the t / samples divide below non-zero.
    HS_CHECK(samples >= 1);
    // Account for the pop_back below: a non-empty path drops its last point
    // before appending samples + 1, so the final size is (size - 1) + samples + 1.
    size_t retained = points.is_empty() ? points.size() : points.size() - 1;
    HS_CHECK(retained + static_cast<size_t>(samples) + 1 <= RESOLUTION);
    if (!points.is_empty())
      points.pop_back();
    for (int t = 0; t <= samples; t++) {
      points.push_back(plot(easing(static_cast<float>(t) / samples) * domain));
    }
    return *this;
  }

  /**
   * @brief Retrieves a point along the path by interpolation.
   * @param t Progress along the path (0.0 to 1.0).
   * @return The linear interpolation of the two adjacent samples. Because it is
   *         a chord between unit-sphere samples, the result has magnitude < 1
   *         (and non-uniform angular speed); callers needing a unit direction
   *         must normalize.
   */
  Vector get_point(float t) const {
    if (points.is_empty())
      return Vector(0, 0, 0);
    // Clamp: a negative t makes raw_index negative, and casting that to size_t
    // is UB (t > 1 is caught by the i >= size-1 guard below, t < 0 is not).
    t = hs::clamp(t, 0.0f, 1.0f);
    float raw_index = t * (points.size() - 1);
    size_t i = static_cast<size_t>(raw_index);
    float f = raw_index - i;
    if (i >= points.size() - 1)
      return points.back();
    const Vector &p1 = points[i];
    const Vector &p2 = points[i + 1];
    return p1 * (1.0f - f) + p2 * f;
  }

  /**
   * @brief Collapses the path to its newest point only.
   */
  void collapse() {
    // Clear in place rather than assigning a fresh buffer: a full
    // StaticCircularBuffer temporary (~12.3 KB) overflows the WASM 8 KB stack.
    if (points.size() > 1) {
      Vector last = points.back();
      points.clear();
      points.push_back(last);
    }
  }

private:
  StaticCircularBuffer<Vector, RESOLUTION> points;
};

/**
 * @brief Represents a path defined by a single procedural function.
 * @details Matches the interface of Path for use in Motion animations.
 */
struct ProceduralPath {
  PlotFn f; /**< The procedural function mapping t to a point. */

  /**
   * @brief Constructs an empty procedural path.
   */
  ProceduralPath() = default;

  /**
   * @brief Constructs a procedural path from a plotting function.
   * @param path_fn Function mapping a parameter to a point on the path.
   */
  ProceduralPath(PlotFn path_fn) : f(std::move(path_fn)) {}

  /**
   * @brief Evaluates the path at a given parameter.
   * @param t Path parameter.
   * @return The point produced by the procedural function at t.
   */
  Vector get_point(float t) const { return f(t); }
};

namespace Animation {

/**
 * @brief Sub-rotation count needed to keep each interpolated step within
 * `max_angle` over a total sweep of `angle` radians.
 * @param angle Total sweep angle in radians.
 * @param max_angle Per-column smoothness threshold in radians.
 * @return Sub-step count, always >= 1 (tight ceil, no extra subdivision).
 * @details The caller upsamples the orientation trail to (result + 1) frames,
 * so the loop's (len-1) sub-intervals each stay <= max_angle. Shared by Motion
 * and Rotation so both size the trail identically for the same angle.
 */
inline int rotation_substeps(float angle, float max_angle) {
  return static_cast<int>(std::ceil(std::max(1.0f, angle / max_angle)));
}

/**
 * @brief An animation that moves an Orientation along a defined Path.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class Motion : public AnimationBase<Motion<W, CAP>> {
  // Interpolates across (len-1) sub-intervals, len <= CAP. CAP == 1 collapses to
  // len == 1: zero sub-intervals (Motion's loop never runs, Rotation divides by
  // zero).
  static_assert(CAP >= 2, "Motion requires an Orientation capacity >= 2");

public:
  /**
   * @brief Constructs a Motion animation.
   * @tparam P Path type exposing get_point(float t).
   * @param orientation The Orientation object to update.
   * @param path_obj The path to follow (must have get_point(float t)).
   * @param duration The duration in frames.
   * @param repeat If true, the motion repeats.
   * @param space Coordinate space for the applied rotations.
   * @note `path_obj` is borrowed by reference and read every frame, so it must
   *       outlive this Motion; a temporary is rejected at compile time (the
   *       rvalue overload is deleted).
   */
  // No ctor emptiness guard (P is generic): an empty/origin-crossing path traps
  // downstream when step() feeds the origin into angle_between()/normalized().
  template <typename P>
  Motion(Orientation<CAP> &orientation, const P &path_obj, int duration,
         bool repeat = false, Space space = Space::World)
      : AnimationBase<Motion<W, CAP>>(duration, repeat),
        orientation(orientation),
        path_fn([&path_obj](float t) { return path_obj.get_point(t); }),
        space(space) {
    // Reject the perpetual -1 the base permits: step() samples path_fn(t/duration),
    // so a -1 duration walks the path at negative, decreasing parameters.
    HS_CHECK(duration >= 0, "Motion duration must be >= 0");
  }

  // Borrow contract: path_fn captures &path_obj and reads it every frame, so the
  // path must outlive the timeline — reject a temporary at compile time.
  template <typename P,
            typename = std::enable_if_t<!std::is_lvalue_reference_v<P>>>
  Motion(Orientation<CAP> &orientation, P &&path_obj, int duration,
         bool repeat = false, Space space = Space::World) = delete;

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return orientation.get(); }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation.
   */
  const void *orientation_id() const override { return &get_orientation(); }

  /**
   * @brief Live-updates the traversal duration (frames per path loop).
   * @param frames New duration in frames; any value below 1 is clamped to 1.
   * @details Clamps to 1 to avoid a divide-by-zero (0) or backward walk
   * (negative). A changed duration reanchors the baseline (see reanchor()): the
   * carried prev_frame_ was sampled under the old parameterization, so the next
   * delta would otherwise teleport the Orientation.
   */
  void set_duration(int frames) {
    int clamped = (frames < 1 ? 1 : frames);
    if (clamped != this->duration)
      have_prev_frame_ = false;
    this->duration = clamped;
  }

  /**
   * @brief Discards the carried baseline frame so the next step re-seeds it from
   * the current path.
   * @details step() integrates the path as RELATIVE deltas off `prev_frame_`.
   * After the borrowed path is swapped, prev_frame_ still holds the old
   * function's frame, so call this to re-seed the baseline and keep the first
   * delta incremental.
   */
  void reanchor() { have_prev_frame_ = false; }

  /**
   * @brief Steps the animation, calculates intermediate rotation steps along
   * the path, and pushes them to the Orientation.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    AnimationBase<Motion<W, CAP>>::step(canvas);

    // Motion advances the Orientation by RELATIVE deltas only, so it composes
    // with any co-driver of the same Orientation. Each delta's frame is a pure
    // function of the path parameter (see path_frame), so deltas telescope
    // drift-free.
    float t_prev = static_cast<float>(this->t - 1);
    Vector current_v = path_fn(t_prev / this->duration);
    float t_curr = static_cast<float>(this->t);
    Vector target_v = path_fn(t_curr / this->duration);
    float total_angle = angle_between(current_v, target_v);
    int num_steps = rotation_substeps(total_angle, MAX_ANGLE);

    orientation.get().upsample(num_steps + 1);
    int len = orientation.get().length();

    // Baseline = the frame Motion last drove to, carried across frames. At a
    // repeat seam the carried baseline makes the first delta a relative jump-back
    // a co-driver rides along. The i==0 delta is identity, so the loop starts at
    // i==1.
    if (!have_prev_frame_) {
      prev_frame_ = path_frame(t_prev / this->duration);
      have_prev_frame_ = true;
    }
    Quaternion base_inv = prev_frame_.conjugate();
    Quaternion frame = prev_frame_;

    for (int i = 1; i < len; ++i) {
      // i maps the sub-interval [t-1, t]: i=0 is t-1, i=len-1 is t.
      float sub_t = t_prev + (static_cast<float>(i) / (len - 1));
      frame = path_frame(sub_t / this->duration);

      // Relative delta from the baseline frame to this substep's frame: World
      // pre-multiplies, Local post-multiplies.
      Quaternion &current_q = orientation.get().at(i);
      if (space == Space::Local) {
        current_q = current_q * (base_inv * frame);
      } else {
        current_q = (frame * base_inv) * current_q;
      }
      current_q.normalize();
    }
    // Carry the final substep's frame as the next frame's baseline.
    prev_frame_ = frame;
  }

  /**
   * @brief The path's own orientation frame at path parameter s.
   * @param s Normalized path parameter (the value passed to path_fn).
   * @return A unit quaternion mapping body +X to the point path(s), body +Y to
   * the travel direction tangent to the sphere, and body +Z to their cross.
   * @details A pure function of s, so consecutive-frame deltas telescope (no
   * drift). Tangent is a fixed-step forward difference; on a vanishing forward
   * tangent it falls back to a backward difference, then to a deterministic
   * perpendicular of the point.
   */
  Quaternion path_frame(float s) const {
    Vector point = path_fn(s).normalized();
    Vector ahead = path_fn(s + FRAME_TANGENT_H).normalized();
    Vector tangent = ahead - point;
    tangent = tangent - dot(tangent, point) * point;
    Vector seed = least_parallel_axis(point);
    Vector b1;
    if (dot(tangent, tangent) < math::EPS_NORMALIZE_SQ) {
      // Forward tangent degenerate (clamped endpoint): use the backward one.
      Vector behind = path_fn(s - FRAME_TANGENT_H).normalized();
      Vector back = point - behind;
      back = back - dot(back, point) * point;
      b1 = normalized_or(back, cross(seed, point).normalized());
    } else {
      b1 = tangent.normalized();
    }
    Vector b2 = cross(point, b1);
    return quaternion_from_basis(point, b1, b2);
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  /** Forward step (in path-parameter space) used to finite-difference the
   * travel tangent in path_frame(). Fixed so the frame is a pure function of s. */
  static constexpr float FRAME_TANGENT_H = 1e-3f;
  std::reference_wrapper<Orientation<CAP>>
      orientation;               /**< Reference to the Orientation state. */
  Fn<Vector(float), 16> path_fn; /**< Function to retrieve path points. */
  Space space;                   /**< The coordinate space for rotation. */
  Quaternion prev_frame_;        /**< Path frame Motion last drove to; the
                                      baseline carried across the repeat seam. */
  bool have_prev_frame_ = false; /**< Whether prev_frame_ has been seeded. */
};

/**
 * @brief An animation that applies a fixed, time-eased rotation.
 * @tparam W The width of the LED display (used for calculating maximum rotation
 * step).
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class Rotation : public AnimationBase<Rotation<W, CAP>> {
  // CAP == 1 collapses len to 1 and makes step_angle = delta / (len - 1) a
  // divide-by-zero; at least two sub-frames are required (see Motion).
  static_assert(CAP >= 2, "Rotation requires an Orientation capacity >= 2");

public:
  /**
   * @brief Default constructor. Creates an inactive/identity rotation.
   */
  Rotation()
      : AnimationBase<Rotation<W, CAP>>(0, false), orientation(nullptr),
        axis(X_AXIS), total_angle(0), easing_fn(ease_linear), last_angle(0),
        space(Space::World) {}

  /**
   * @brief Constructs a Rotation animation.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis; normalized here so a non-unit axis cannot
   *        silently skew the angle (make_rotation's trailing normalize folds a
   *        length-L axis into tan(φ/2) = L·tan(θ/2), rotating by the wrong
   *        amount with no visual tell). A zero axis traps in normalized().
   * @param angle The total rotation angle in radians.
   * @param duration The duration in frames.
   * @param easing_fn The easing function to use.
   * @param repeat If true, the rotation repeats.
   * @param space The coordinate space for rotation ("World" or "Local").
   */
  Rotation(Orientation<CAP> &orientation, const Vector &axis, float angle,
           int duration, EasingFn easing_fn, bool repeat = false,
           Space space = Space::World)
      : AnimationBase<Rotation<W, CAP>>(duration, repeat),
        orientation(&orientation), axis(axis.normalized()), total_angle(angle),
        easing_fn(std::move(easing_fn)), last_angle(0), space(space) {
    HS_CHECK(duration >= 0, "Rotation duration must be >= 0");
  }

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return *orientation; }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation (nullptr if unbound).
   */
  const void *orientation_id() const override { return orientation; }

  /**
   * @brief Checks if the rotation has a valid orientation bound.
   * @return True if an Orientation is bound.
   */
  bool has_orientation() const { return orientation != nullptr; }

  /**
   * @brief Steps the animation, calculates the incremental rotation delta, and
   * pushes it to the Orientation.
   * @param canvas The canvas buffer (forwarded to the base step).
   */
  void step(Canvas &canvas) override {
    HS_CHECK(orientation != nullptr);
    if (this->t == 0) {
      last_angle = 0;
    }
    AnimationBase<Rotation<W, CAP>>::step(canvas);
    float target_angle =
        easing_fn(static_cast<float>(this->t) / this->duration) * total_angle;
    float delta = target_angle - last_angle;

    // Sub-threshold increments accumulate across a multi-frame sweep: last_angle
    // persists until the sum crosses TOLERANCE. A one-shot (duration 1, the
    // animate() path) has no next frame to accumulate into, so it must apply
    // unconditionally or a very slow per-frame driver freezes forever.
    if (this->duration > 1 && std::abs(delta) < TOLERANCE) {
      return;
    }

    int num_steps = rotation_substeps(std::abs(delta), MAX_ANGLE);
    orientation->upsample(num_steps + 1);
    int len = orientation->length();

    float step_angle = delta / (len - 1);

    auto apply_rotation = [&](Quaternion &target, const Quaternion &source) {
      if (space == Space::Local) {
        target = target * source;
      } else {
        target = source * target;
      }
    };

    for (int i = 1; i < len; ++i) {
      float angle = step_angle * i;
      Quaternion q = make_rotation(axis, angle);

      Quaternion &current_q = orientation->at(i);
      apply_rotation(current_q, q);
      current_q.normalize();
    }
    last_angle = target_angle;
  }

  /**
   * @brief Static helper to apply a rotation immediately (duration = 1 frame).
   * @param canvas The canvas buffer.
   * @param orientation The Orientation object to update.
   * @param axis The rotation axis (unit vector).
   * @param angle The rotation angle in radians.
   * @param easing_fn The easing function to use.
   * @param space The coordinate space for rotation.
   */
  static void animate(Canvas &canvas, Orientation<CAP> &orientation,
                      const Vector &axis, float_t angle, EasingFn easing_fn,
                      Space space = Space::World) {
    Rotation<W, CAP> r(orientation, axis, angle, 1, easing_fn, false, space);
    r.step(canvas);
  }

private:
  static constexpr float MAX_ANGLE =
      2 * PI_F /
      W; /**< Maximum rotation angle per step to ensure smoothness. */
  Orientation<CAP> *orientation; /**< Pointer to the Orientation state. */
  Vector axis;                      /**< The axis of rotation. */
  float total_angle;                /**< The total angle to sweep. */
  EasingFn easing_fn;               /**< Easing curve. */
  float last_angle; /**< The angle reached in the previous frame. */
  Space space;      /**< The coordinate space for rotation. */
};

/**
 * @brief An animation that simulates a particle/camera performing a random walk
 * across the sphere's surface.
 * @details Uses Perlin noise to create continuous, turbulent pivoting motion.
 *
 * PERPETUAL (duration -1, no repeat or self-reset): never reaches done(), so a
 * `.then()` callback NEVER fires and a repeating spawn leaks recycled slots.
 * Drive follow-on behavior from a finite animation or cancel() the walk
 * explicitly.
 * @tparam W The width of the LED display.
 * @tparam CAP Orientation sub-frame capacity.
 */
template <int W, int CAP = 4>
class RandomWalk : public AnimationBase<RandomWalk<W, CAP>> {
public:
  /**
   * @brief Tunable parameters for the random walk.
   */
  struct Options {
    float speed = 0.02f; /**< Movement speed per frame. */
    float pivot_strength =
        0.1f; /**< Strength of the direction change (noise amplitude). */
    float noise_scale = 0.02f;  /**< Frequency of the Perlin noise. */
    float smoothing = 0.85f;    /**< Angular momentum (0 = none, 0.95 = very sluggish). */
    float drift = 0.5f;         /**< Temporal drift speed for spatial noise. */

    /**
     * @brief Preset for slow, gentle motion.
     * @return Options tuned for a languid walk.
     */
    static Options Languid() { return {0.02f, 0.1f, 0.02f, 0.85f, 0.5f}; }

    /**
     * @brief Preset for fast, lively motion.
     * @return Options tuned for an energetic walk.
     */
    static Options Energetic() { return {0.05f, 0.4f, 0.08f, 0.7f, 1.0f}; }
  };

  /**
   * @brief Constructs a RandomWalk animation.
   * @param orientation The Orientation object to update.
   * @param v_start The starting direction vector.
   * @param noise External noise generator (caller owns lifetime).
   * @param options Configuration options.
   * @param seed Noise seed; 0 selects a random seed.
   */
  RandomWalk(Orientation<CAP> &orientation, const Vector &v_start,
             FastNoiseLite &noise, Options options = Options(),
             int seed = 0)
      : AnimationBase<RandomWalk<W, CAP>>(-1, false), orientation(orientation),
        v(Vector(v_start).normalized()), options(options),
        noise_generator(noise) {
    Vector u = least_parallel_axis(v);
    direction = cross(v, u).normalized();
    noise_generator.get().SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise_generator.get().SetFrequency(options.noise_scale);
    if (seed == 0) {
      noise_generator.get().SetSeed(static_cast<int>(hs::random()()));
    } else {
      noise_generator.get().SetSeed(seed);
    }
  }

  /**
   * @brief Live-updates the per-frame movement speed (e.g. from a GUI slider).
   * @param new_speed The new per-frame movement speed.
   */
  void set_speed(float new_speed) { options.speed = new_speed; }

  /**
   * @brief Accesses the associated Orientation.
   * @return Reference to the bound Orientation.
   */
  Orientation<CAP> &get_orientation() const { return orientation.get(); }

  /**
   * @brief Collapses the bound Orientation's motion-blur history.
   */
  void collapse_orientation() override { get_orientation().collapse(); }

  /**
   * @brief Identity of the bound Orientation.
   * @return Pointer identifying the bound Orientation.
   */
  const void *orientation_id() const override { return &get_orientation(); }

  /**
   * @brief Steps the walk: pivots direction based on noise, then rotates the
   * view along the calculated axis.
   * @param canvas The canvas buffer (forwarded to the base step and rotation).
   */
  void step(Canvas &canvas) override {
    AnimationBase<RandomWalk<W, CAP>>::step(canvas);
    // noise_scale is applied once via SetFrequency() in the ctor; the 100x is a
    // fixed base sample scale (scaling coords by noise_scale here too would make
    // the spatial frequency quadratic in it).
    // Accepted limit: past t == 2^24 (~77 h at 60 fps, sooner at higher drift)
    // float can't represent consecutive frames and the drift coordinate freezes.
    float target_pivot =
        noise_generator.get().GetNoise(
            v.x * 100.0f, v.y * 100.0f,
            v.z * 100.0f + static_cast<float>(this->t) * options.drift) *
        options.pivot_strength;
    angular_velocity = angular_velocity * options.smoothing +
                       target_pivot * (1.0f - options.smoothing);
    direction = rotate(direction, make_rotation(v, angular_velocity)).normalized();
    // If v and direction drift near-parallel the cross collapses to zero; fall
    // back to a deterministic perpendicular of v.
    const Vector axis_seed = least_parallel_axis(v);
    Vector walk_axis =
        normalized_or(cross(v, direction), cross(v, axis_seed).normalized());
    v = rotate(v, make_rotation(walk_axis, options.speed)).normalized();
    direction =
        rotate(direction, make_rotation(walk_axis, options.speed)).normalized();
    Rotation<W, CAP>::animate(canvas, orientation, walk_axis, options.speed,
                              ease_linear);
  }

private:
  std::reference_wrapper<Orientation<CAP>>
      orientation;  /**< Reference to the global Orientation state. */
  Vector v;         /**< Current forward direction vector. */
  Vector direction; /**< Current pivoting direction (orthogonal to v). */
  Options options;  /**< Configuration options. */
  float angular_velocity = 0.0f; /**< Smoothed pivot rate (angular momentum). */
  std::reference_wrapper<FastNoiseLite>
      noise_generator; /**< External noise generator. */
};


} // namespace Animation
