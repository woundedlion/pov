/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Column tolerance a timed route holds interpolants to, radians. */
inline constexpr float MOTION_TOLERANCE = 0.2f * math::TWO_PI_F / 288.0f;

/**
 * @brief Largest distance of a path's interpolant from the exact motion.
 * @tparam Exact Callable `math::Vector(float t)`.
 * @param path Interpolated path.
 * @param exact The motion it stands for.
 * @return The largest distance over 1,001 evenly spaced fractions.
 */
template <typename Exact>
inline float interpolant_departure(const Temporal::VertexPath &path,
                                   Exact &&exact) {
  float worst = 0.0f;
  for (int i = 0; i <= 1000; ++i) {
    const float t = static_cast<float>(i) / 1000.0f;
    worst = fmaxf(worst, (path.at(t) - exact(t)).magnitude());
  }
  return worst;
}

/**
 * @brief Largest distance of a path from its chord, point for point.
 * @param path Interpolated path.
 * @return The largest `|at(t) - chord(t)|` over 1,001 fractions.
 */
inline float chord_departure(const Temporal::VertexPath &path) {
  float worst = 0.0f;
  for (int i = 0; i <= 1000; ++i) {
    const float t = static_cast<float>(i) / 1000.0f;
    const math::Vector chord = path.start + (path.end - path.start) * t;
    worst = fmaxf(worst, (path.at(t) - chord).magnitude());
  }
  return worst;
}

/**
 * @brief Checks each interpolant of one exact motion against its bound.
 * @tparam Exact Callable `math::Vector(float t)`.
 * @param exact The motion.
 * @param bound Its derivative bounds.
 * @param axis Its circle's axis, or null for a great circle.
 * @param on_arc Whether the motion is a circle at a constant rate.
 */
template <typename Exact>
inline void check_interpolants(Exact &&exact, const Temporal::PathBound &bound,
                               const math::Vector *axis, bool on_arc) {
  const math::Vector a = exact(0.0f), m = exact(0.5f), b = exact(1.0f);
  const auto chord = Temporal::VertexPath::chord(a, b);
  const float chord_error = interpolant_departure(chord, exact);
  HS_EXPECT_LE(chord_error, bound.chord_error() + 1e-6f);
  const auto parabola = Temporal::VertexPath::quadratic(a, m, b);
  HS_EXPECT_LE(interpolant_departure(parabola, exact),
               bound.quadratic_error() + 1e-6f);
  HS_EXPECT_LE(chord_departure(parabola), parabola.deviation() + 1e-6f);
  if (on_arc) {
    const auto arc = Temporal::VertexPath::arc(a, b, axis);
    HS_EXPECT_TRUE(arc.kind == Temporal::Interpolant::ARC);
    HS_EXPECT_LE(interpolant_departure(arc, exact), 2e-6f);
    HS_EXPECT_LE(chord_departure(arc), arc.deviation() + 1e-6f);
  }
  const Temporal::Interpolant kind =
      Temporal::choose_interpolant(bound, MOTION_TOLERANCE, on_arc);
  if (kind == Temporal::Interpolant::CHORD)
    HS_EXPECT_LE(chord_error, MOTION_TOLERANCE);
  else
    HS_EXPECT_GT(bound.chord_error(), MOTION_TOLERANCE);
  HS_EXPECT_TRUE((kind == Temporal::Interpolant::ARC) ==
                 (on_arc && bound.chord_error() > MOTION_TOLERANCE));
  HS_EXPECT_LE(Temporal::interpolant_error(kind, bound),
               kind == Temporal::Interpolant::CHORD ? MOTION_TOLERANCE
                                                    : bound.quadratic_error());
}

/** @brief Rigid rotation: the chord bound holds and is tight, and the arc
 *         about the rotation axis is exact. */
inline void test_interpolants_bound_rotation() {
  const math::Vector axis = math::Vector(0.2f, 1.0f, -0.3f).normalized();
  for (const float angle : {0.01f, 0.1f, 0.4f, -0.7f}) {
    for (const float lift : {0.0f, 0.5f, 0.9f}) {
      HS_CONTEXT("angle, lift x100", std::llround(angle * 100.0f),
                 std::llround(lift * 100.0f));
      const math::Vector across = math::cross(axis, math::X_AXIS).normalized();
      const math::Vector p0 =
          (axis * lift + across * sqrtf(1.0f - lift * lift)).normalized();
      const float radius = sqrtf(1.0f - lift * lift);
      const auto exact = [&](float t) {
        return math::rotate(p0, math::make_rotation(axis, angle * t));
      };
      const auto bound = Temporal::PathBound::rotation(angle, radius);
      check_interpolants(exact, bound, &axis, true);
      const auto chord = Temporal::VertexPath::chord(exact(0.0f), exact(1.0f));
      // The sagitta r (1 - cos(angle / 2)) approaches r angle^2 / 8.
      HS_EXPECT_GE(interpolant_departure(chord, exact),
                   0.97f * bound.chord_error());
    }
  }
}

/** @brief A meridian at a constant rate is a great arc the ARC interpolant
 *         reproduces without an axis. */
inline void test_interpolants_bound_meridian() {
  for (const float angle : {0.05f, 0.3f, 1.2f}) {
    HS_CONTEXT("angle x100", std::llround(angle * 100.0f));
    const auto exact = [&](float t) {
      const float colatitude = 0.4f + angle * t;
      return math::Vector(sinf(colatitude) * 0.6f, cosf(colatitude),
                          sinf(colatitude) * 0.8f);
    };
    check_interpolants(exact, Temporal::PathBound::rotation(angle), nullptr,
                       true);
  }
}

/** @brief An orbit on a small circle about a tilted axis. */
inline void test_interpolants_bound_orbit() {
  const math::Vector axis = math::Vector(0.5f, 0.4f, 0.7f).normalized();
  const math::Vector centre = axis * 0.8f;
  const math::Vector u = math::cross(axis, math::Y_AXIS).normalized();
  const math::Vector w = math::cross(axis, u);
  for (const float angle : {0.08f, 0.6f, 2.5f}) {
    HS_CONTEXT("angle x100", std::llround(angle * 100.0f));
    const auto exact = [&](float t) {
      const float turn = 0.3f + angle * t;
      return centre + (u * cosf(turn) + w * sinf(turn)) * 0.6f;
    };
    check_interpolants(exact, Temporal::PathBound::rotation(angle, 0.6f), &axis,
                       true);
  }
}

/** @brief A curved non-rigid path: the parabola through the middle capture
 *         meets the third-derivative bound, which is tight for a cubic. */
inline void test_interpolants_bound_curved_path() {
  const math::Vector a(0.1f, 0.7f, 0.6f), v(0.04f, -0.02f, 0.03f),
      w(-0.05f, 0.03f, 0.02f), z(0.03f, 0.05f, -0.04f);
  for (const float scale : {0.2f, 1.0f, 3.0f}) {
    HS_CONTEXT("scale x10", std::llround(scale * 10.0f));
    const auto exact = [&](float t) {
      return a + (v * t + w * (t * t) + z * (t * t * t)) * scale;
    };
    const Temporal::PathBound bound{
        (2.0f * w.magnitude() + 6.0f * z.magnitude()) * scale,
        6.0f * z.magnitude() * scale};
    check_interpolants(exact, bound, nullptr, false);
    const auto parabola =
        Temporal::VertexPath::quadratic(exact(0.0f), exact(0.5f), exact(1.0f));
    HS_EXPECT_GE(interpolant_departure(parabola, exact),
                 0.99f * bound.quadratic_error());
  }
}

/** @brief A path linear in Mercator ordinate and longitude stays within its
 *         loxodrome bounds, which it meets crossing the equator. */
inline void test_interpolants_bound_loxodrome() {
  struct Motion {
    double ordinate, longitude, start;
  };
  const auto position = [](const Motion &motion, double t, int axis) {
    const double u = motion.start + motion.ordinate * t;
    const double theta = 0.4 + motion.longitude * t;
    const double horizontal = 1.0 / cosh(u);
    return axis == 0   ? horizontal * cos(theta)
           : axis == 1 ? horizontal * sin(theta)
                       : tanh(u);
  };
  for (const Motion motion :
       {Motion{0.3, 0.0, -0.15}, Motion{0.0, 0.5, 0.0}, Motion{0.2, -0.3, -0.1},
        Motion{0.6, 0.4, 1.1}, Motion{-0.25, 0.9, -1.7}}) {
    HS_CONTEXT("ordinate, longitude x100", std::llround(motion.ordinate * 100),
               std::llround(motion.longitude * 100));
    const auto bound =
        Temporal::PathBound::loxodrome(static_cast<float>(motion.ordinate),
                                       static_cast<float>(motion.longitude));
    const auto exact = [&](float t) {
      return math::Vector(static_cast<float>(position(motion, t, 0)),
                          static_cast<float>(position(motion, t, 1)),
                          static_cast<float>(position(motion, t, 2)));
    };
    check_interpolants(exact, bound, nullptr, false);
    constexpr double STEP = 1e-3;
    double second = 0.0, third = 0.0;
    for (int i = 0; i <= 100; ++i) {
      const double t = i / 100.0;
      double d2 = 0.0, d3 = 0.0;
      for (int axis = 0; axis < 3; ++axis) {
        const auto at = [&](int k) {
          return position(motion, t + k * STEP, axis);
        };
        const double c2 = (at(1) - 2.0 * at(0) + at(-1)) / (STEP * STEP);
        const double c3 = (at(2) - 2.0 * at(1) + 2.0 * at(-1) - at(-2)) /
                          (2.0 * STEP * STEP * STEP);
        d2 += c2 * c2;
        d3 += c3 * c3;
      }
      second = std::max(second, sqrt(d2));
      third = std::max(third, sqrt(d3));
    }
    HS_EXPECT_LE(second, bound.second * 1.0001 + 1e-6);
    HS_EXPECT_LE(third, bound.third * 1.0001 + 1e-6);
    if (motion.start + 0.5 * motion.ordinate == 0.0) {
      HS_EXPECT_GE(second, 0.999 * bound.second);
      if (motion.longitude == 0.0)
        HS_EXPECT_GE(third, 0.999 * bound.third);
    }
  }
}

/** @brief A point turning along a great circle under a moving plane
 *         translation stays within its warped bounds, which a fixed
 *         translation of a turning point and a moving translation of a still
 *         point meet. */
inline void test_interpolants_bound_warped_great_circle() {
  struct Motion {
    double x, y, z;  ///< Frame-end source position.
    double hx, hz;   ///< Turn direction, completed to a tangent.
    double turn;     ///< Turn over the frame, radians.
    double b0x, b0z; ///< Translation at the frame's start.
    double b1x, b1z; ///< Translation at the frame's end.
  };
  const auto position = [](const Motion &m, double t, double out[3]) {
    double p[3] = {m.x, m.y, m.z};
    const double n = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    for (double &c : p)
      c /= n;
    double h[3] = {m.hx, 0.37, m.hz};
    const double along = h[0] * p[0] + h[1] * p[1] + h[2] * p[2];
    for (int i = 0; i < 3; ++i)
      h[i] -= along * p[i];
    const double hn = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    const double back = m.turn * (1.0 - t);
    double q[3];
    for (int i = 0; i < 3; ++i)
      q[i] = p[i] * cos(back) - h[i] / hn * sin(back);
    // Stereographic from +y, 1 - y kept exact near the pole.
    const double s =
        q[1] > 0.0 ? (q[0] * q[0] + q[2] * q[2]) / (1.0 + q[1]) : 1.0 - q[1];
    const double bx = m.b0x + (m.b1x - m.b0x) * t;
    const double bz = m.b0z + (m.b1z - m.b0z) * t;
    const double nx = q[0] + bx * s, nz = q[2] + bz * s;
    const double n2 = nx * nx + nz * nz, s2 = s * s, den = n2 + s2;
    out[0] = 2.0 * nx * s / den;
    out[1] = (n2 - s2) / den;
    out[2] = 2.0 * nz * s / den;
  };
  uint32_t state = 29u;
  auto uniform = [&](double lo, double hi) {
    state = state * 1664525u + 1013904223u;
    return lo + (hi - lo) * double(state >> 8) * (1.0 / 16777216.0);
  };
  double tight_second = 0.0, tight_third = 0.0;
  for (int trial = 0; trial < 3000; ++trial) {
    Motion m;
    const double polar =
        trial % 3 == 0 ? uniform(1e-3, 0.05) : uniform(0.0, 3.1);
    const double azimuth = uniform(0.0, 6.283);
    m.x = sin(polar) * cos(azimuth);
    m.y = cos(polar);
    m.z = sin(polar) * sin(azimuth);
    m.hx = uniform(-1.0, 1.0);
    m.hz = uniform(-1.0, 1.0);
    const int kind = trial % 4;
    m.turn = kind == 1 ? 0.0 : uniform(0.001, 0.0218);
    const double scale = uniform(0.0, 5.0);
    const double start = uniform(0.0, 6.283);
    const double sweep = kind == 0 ? 0.0 : uniform(0.0, 0.08);
    m.b0x = scale * (cos(start) - 1.0);
    m.b0z = scale * sin(start);
    m.b1x = scale * (cos(start + sweep) - 1.0);
    m.b1z = scale * sin(start + sweep);
    const double shift = hypot(m.b1x - m.b0x, m.b1z - m.b0z);
    const double reach = std::max(hypot(m.b0x, m.b0z), hypot(m.b1x, m.b1z));
    const auto bound = Temporal::PathBound::warped_great_circle(
        float(m.turn), float(shift), float(reach));
    constexpr double STEP = 1e-2;
    double second = 0.0, third = 0.0;
    for (int i = 0; i <= 20; ++i) {
      const double t = i / 20.0;
      double at[5][3];
      for (int k = -2; k <= 2; ++k)
        position(m, t + k * STEP, at[k + 2]);
      double d2 = 0.0, d3 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double c2 =
            (at[3][c] - 2.0 * at[2][c] + at[1][c]) / (STEP * STEP);
        const double c3 =
            (at[4][c] - 2.0 * at[3][c] + 2.0 * at[1][c] - at[0][c]) /
            (2.0 * STEP * STEP * STEP);
        d2 += c2 * c2;
        d3 += c3 * c3;
      }
      second = std::max(second, sqrt(d2));
      third = std::max(third, sqrt(d3));
    }
    HS_CONTEXT("trial", trial);
    HS_EXPECT_LE(second, double(bound.second) * 1.001 + 1e-9);
    HS_EXPECT_LE(third, double(bound.third) * 1.001 + 1e-9);
    if (bound.second > 0.0f)
      tight_second = std::max(tight_second, second / double(bound.second));
    if (bound.third > 0.0f)
      tight_third = std::max(tight_third, third / double(bound.third));
  }
  HS_EXPECT_GE(tight_second, 0.9);
  HS_EXPECT_GE(tight_third, 0.5);
}

/** @brief Frame-end captures carry into the next frame's start, so each
 *         vertex costs one new evaluation, plus one at the middle for a
 *         QUADRATIC path. */
inline void test_vertex_motion_carries_frame_ends() {
  constexpr size_t VERTICES = 5;
  ScratchScope scope(scratch_arena_a);
  Temporal::VertexMotion<math::Vector> exact_motion;
  Temporal::VertexMotion<math::Snorm3> packed_motion;
  exact_motion.bind(scratch_arena_a, VERTICES);
  packed_motion.bind(scratch_arena_a, VERTICES);
  HS_EXPECT_EQ(exact_motion.capacity(), VERTICES);
  int evaluations = 0;
  const auto position = [&](size_t i, float time) {
    ++evaluations;
    const float turn = 0.3f * time * static_cast<float>(i + 1);
    return math::Vector(cosf(turn) * 0.8f, 0.6f, sinf(turn) * 0.8f);
  };
  for (int frame = 0; frame < 4; ++frame) {
    HS_CONTEXT("frame", frame);
    std::array<math::Vector, VERTICES> previous_end{};
    for (size_t i = 0; frame > 0 && i < VERTICES; ++i)
      previous_end[i] = exact_motion.end(i);
    exact_motion.begin_frame(VERTICES, frame > 0);
    packed_motion.begin_frame(VERTICES, frame > 0);
    evaluations = 0;
    int middles = 0;
    for (size_t i = 0; i < VERTICES; ++i) {
      const math::Vector end = position(i, static_cast<float>(frame + 1));
      exact_motion.capture_end(i, end);
      packed_motion.capture_end(i, end);
      const float angle = 0.3f * static_cast<float>(i + 1);
      const auto kind =
          exact_motion.choose(i, Temporal::PathBound::rotation(angle, 0.8f),
                              MOTION_TOLERANCE, false);
      if (kind == Temporal::Interpolant::QUADRATIC) {
        ++middles;
        exact_motion.capture_middle(
            i, position(i, static_cast<float>(frame) + 0.5f));
      }
      HS_EXPECT_TRUE(exact_motion.interpolant(i) == kind);
    }
    HS_EXPECT_EQ(evaluations, static_cast<int>(VERTICES) + middles);
    HS_EXPECT_GT(middles, 0);
    for (size_t i = 0; i < VERTICES; ++i) {
      const math::Vector start = exact_motion.start(i);
      const math::Vector want =
          frame > 0 ? previous_end[i] : exact_motion.end(i);
      HS_EXPECT_EQ(start.x, want.x);
      HS_EXPECT_EQ(start.y, want.y);
      HS_EXPECT_EQ(start.z, want.z);
      HS_EXPECT_LE((packed_motion.end(i) - exact_motion.end(i)).magnitude(),
                   3e-5f);
      HS_EXPECT_LE((packed_motion.start(i) - exact_motion.start(i)).magnitude(),
                   3e-5f);
      const auto path = exact_motion.path(i);
      const auto truth = [&](float t) {
        const float turn =
            0.3f * (static_cast<float>(frame) + t) * static_cast<float>(i + 1);
        return math::Vector(cosf(turn) * 0.8f, 0.6f, sinf(turn) * 0.8f);
      };
      if (frame > 0)
        HS_EXPECT_LE(
            interpolant_departure(path, truth),
            Temporal::interpolant_error(
                path.kind, Temporal::PathBound::rotation(
                               0.3f * static_cast<float>(i + 1), 0.8f)) +
                1e-6f);
    }
    for (size_t i = 0; i < VERTICES; ++i) {
      exact_motion.hold(i);
      const math::Vector held = exact_motion.start(i);
      HS_EXPECT_EQ(held.x, exact_motion.end(i).x);
    }
  }
}

/** @brief The chord bound between two rotations bounds every unit vector's
 *         move. */
inline void test_rotation_chord_bound_covers_unit_vectors() {
  const math::RotationMatrix a(math::make_rotation(math::Y_AXIS, 0.4f));
  const math::RotationMatrix b(
      math::make_rotation(math::Vector(0.3f, 1.0f, 0.1f).normalized(), 0.55f));
  const float bound = Temporal::rotation_chord_bound(a, b);
  float worst = 0.0f;
  for (int i = 0; i < 400; ++i) {
    const float u = static_cast<float>(i) * 0.61803f;
    const float z = 1.0f - 2.0f * (static_cast<float>(i) + 0.5f) / 400.0f;
    const float r = sqrtf(1.0f - z * z);
    const math::Vector v(r * cosf(u * math::TWO_PI_F), z,
                         r * sinf(u * math::TWO_PI_F));
    worst = fmaxf(worst, (a.apply(v) - b.apply(v)).magnitude());
  }
  HS_EXPECT_LE(worst, bound);
  HS_EXPECT_GT(worst, 0.5f * bound);
}
