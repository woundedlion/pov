/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file rings_temporal.h
 * @brief Frozen-column oracles for growing and warped rings. */

/** @brief Requires lit, matching intensity on opposite arms. */
inline void expect_arm_symmetric(std::span<const float> pixels) {
  float energy = 0;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W / 2; ++x) {
      energy += pixels[y * DEFAULT_W + x];
      HS_EXPECT_NEAR(pixels[y * DEFAULT_W + x],
                     pixels[y * DEFAULT_W + x + DEFAULT_W / 2], 3e-5f);
    }
  HS_EXPECT_GT(energy, 0.01f);
}

/** @brief Frozen ring geometry and luminance comparisons. */
struct RingTemporalOracle {
  /** @brief Direction of a physical reference displacement. */
  enum class Direction { NORMAL, AZIMUTH, MERIDIAN };
  /** @brief Signed physical reference displacement. */
  struct Offset {
    float pixels;        /**< Row-pitch pixel distance. */
    Direction direction; /**< Unit tangent direction. */
    /** @brief Sets a reference displacement.
     * @param pixels Signed pixel distance.
     * @param direction Tangent direction. */
    Offset(float pixels, Direction direction = Direction::NORMAL)
        : pixels(pixels), direction(direction) {}
  };
  /** @brief Selects a focused scenario when requested by the native runner.
   * @param scenario Scenario number.
   * @return Whether this scenario is selected. */
  static bool selected(int scenario) {
    const char *filter = hs_test::env_var("HS_RING_TEMPORAL_SCENARIO");
    return !filter || std::atoi(filter) == scenario;
  }

  /** @brief Displaces geometry by a physical pixel distance normal to its stroke.
   * @param point Unit direction.
   * @param stroke Tangent along the stroke.
   * @param offset Signed physical displacement.
   * @return Direction turned along the unit normal tangent. */
  static math::Vector displace(const math::Vector &point,
                               const math::Vector &stroke, Offset offset) {
    if (offset.pixels == 0.0f)
      return point;
    const double r = std::hypot(double(point.x), double(point.z));
    const math::Vector azimuth =
        r > 0.0 ? math::Vector(float(-point.z / r), 0, float(point.x / r))
                : math::X_AXIS;
    math::Vector tangent = math::normalized_or(math::cross(stroke, point),
                                               math::cross(azimuth, point));
    if (offset.direction == Direction::AZIMUTH)
      tangent = azimuth;
    if (offset.direction == Direction::MERIDIAN)
      tangent = math::cross(azimuth, point);
    const double angle =
        double(offset.pixels) * math::RADIANS_PER_ROW<DEFAULT_H>;
    return (point * float(std::cos(angle)) + tangent * float(std::sin(angle)))
        .normalized();
  }

  /** @brief Draws independently evaluated vertices at one frozen instant.
   * @param canvas Destination canvas.
   * @param vertex Independent vertex evaluator.
   * @param shader Frozen fragment shader.
   * @param rows Reference displacement in pixels.
   * @param fraction Frozen frame fraction.
   * @param stride Material-grid spacing. */
  template <typename Vertex, typename Shader>
  static void frozen(Canvas &canvas, Vertex vertex, Shader shader, Offset rows,
                     float fraction, int stride) {
    ScratchScope storage(scratch_arena_a);
    const int count = DEFAULT_W / stride;
    auto *points = scratch_arena_a.allocate_n<math::Vector>(DEFAULT_W);
    auto *unshifted = scratch_arena_a.allocate_n<math::Vector>(DEFAULT_W);
    auto *starts = scratch_arena_a.allocate_n<math::Vector>(DEFAULT_W);
    auto *ends = scratch_arena_a.allocate_n<math::Vector>(DEFAULT_W);
    ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_a, DEFAULT_W);
    for (int i = 0; i < count; ++i) {
      unshifted[i] = vertex(i * stride, fraction);
      starts[i] = vertex(i * stride, 0.0f);
      ends[i] = vertex(i * stride, 1.0f);
    }
    for (int i = 0; i < count; ++i)
      points[i] = displace(unshifted[i],
                           unshifted[(i + 1) % count] -
                               unshifted[(i + count - 1) % count],
                           rows);
    Pipeline<DEFAULT_W, DEFAULT_H,
             Filter::Screen::AntiAlias<DEFAULT_W, DEFAULT_H>>
        pipeline;
    Plot::PointRoute<DEFAULT_W, DEFAULT_H> route;
    route.prepare(canvas, 0.0f);
    for (int i = 0; i < count; ++i) {
      const int next = (i + 1) % count;
      const math::Vector first[2] = {starts[i], starts[next]};
      const math::Vector last[2] = {ends[i], ends[next]};
      auto shade = [&](float) {
        Fragment f;
        f.pos = points[i];
        shader(f.pos, f);
        return f.color;
      };
      Plot::GeodesicEdgeSampler arcs[2];
      if (Plot::EdgeWorldsheet::steady_lattice<DEFAULT_W, DEFAULT_H>(
              first, last, 1, arcs) == 1 &&
          Plot::EdgeWorldsheet::pieces_of<DEFAULT_W, DEFAULT_H>(arcs[0]) <=
              0.5f &&
          Plot::EdgeWorldsheet::pieces_of<DEFAULT_W, DEFAULT_H>(arcs[1]) <=
              0.5f)
        route.plot(
            pipeline, canvas, points[i], [&] { return points[i]; }, shade);
      else
        edges.push_back({uint16_t(i), uint16_t(next)});
    }
    Plot::GeometryCallback geometry;
    geometry.chained = true;
    geometry.rigid_lattice = true;
    Plot::Mesh::draw_columns<DEFAULT_W, DEFAULT_H>(
        pipeline, canvas, Plot::Mesh::VertexSpan{{points, size_t(count)}},
        edges, shader, geometry, 0);
  }

  /** @brief Reads displayed linear luminance.
   * @param effect Displayed effect.
   * @return Row-major luminance. */
  template <typename FX> static std::vector<float> pixels(FX &effect) {
    effect.advance_display();
    return capture_luminance<DEFAULT_W, DEFAULT_H>(effect);
  }

  /** @brief Compares full and segmented draws with frozen columns.
   * @param render Timed or independently frozen renderer.
   * @param name Diagnostic label.
   * @param moving Whether the frame has continuous motion. */
  template <typename Render>
  static void compare(Render render, const char *name, bool moving) {
    const auto actual = render(true, 1.0f, 0.0f, -1);
    hs_test::SeamCompareOptions options;
    options.bidirectional = true;
    options.require_motion = moving;
    options.widen_budget = [&](int x, float t, std::span<const float> expected,
                               float budget) {
      for (Direction direction : {Direction::AZIMUTH, Direction::MERIDIAN})
        for (float sign : {-1.0f, 1.0f}) {
          const auto displaced =
              render(false, t, Offset(sign * 0.2f, direction), -1);
          float axis_budget = 0.0f;
          for (int y = 0; y < DEFAULT_H; ++y) {
            const size_t i = y * DEFAULT_W + x;
            axis_budget += fabsf(displaced[i] - expected[i]);
          }
          budget = std::max(budget, axis_budget);
        }
      return budget;
    };
    const auto result = hs_test::seam_compare_render<DEFAULT_W, DEFAULT_H>(
        [&](bool timed, float fraction, float displacement) {
          return render(timed, fraction, displacement, -1);
        },
        name, options);
    HS_EXPECT_GT(result.energy, 0.01f);
    for (int segment = 0; segment < 8; ++segment) {
      const auto clipped = render(true, 1.0f, 0.0f, segment);
      const auto bounds = segment_bounds(segment);
      for (int y = bounds.y0; y < bounds.y1; ++y)
        for (int x = bounds.x0; x < bounds.x1; ++x)
          HS_EXPECT_NEAR(clipped[y * DEFAULT_W + x], actual[y * DEFAULT_W + x],
                         1e-5f);
    }
  }

  /** @brief Computes one device segment's pixel bounds.
   * @param segment Segment number.
   * @return Half-open row and column bounds. */
  static pov::SegmentClip segment_bounds(int segment) {
    const auto map = pov::segment_map(segment % 4, DEFAULT_H * 2, 4);
    return pov::segment_clip(map, segment < 4, DEFAULT_H * 2, 4, DEFAULT_W);
  }

  /** @brief Applies the device's segment clip.
   * @param effect Destination effect.
   * @param segment Segment number, or -1 for a full frame. */
  template <typename FX> static void clip(FX &effect, int segment) {
    if (segment < 0)
      return;
    const auto bounds = segment_bounds(segment);
    effect.set_clip(bounds.y0, bounds.y1, bounds.x0, bounds.x1);
  }

  /** @brief Evaluates a ring vertex independently in double precision.
   * @param basis Ring basis.
   * @param i Material vertex index.
   * @param polar Colatitude in radians.
   * @return Unit direction. */
  static math::Vector vertex(const math::Basis &basis, int i, double polar) {
    const double angle = 2.0 * double(math::PI_F) * i / DEFAULT_W;
    return (basis.v * float(std::cos(polar)) +
            (basis.u * float(std::cos(angle)) +
             basis.w * float(std::sin(angle))) *
                float(std::sin(polar)))
        .normalized();
  }

  /** @brief Captures use two endpoint evaluations per distinct vertex.
   * @param vertices Distinct material vertices. */
  static void counts(size_t vertices) {
    HS_EXPECT_EQ(Spacetime::geometry_counters.endpoint_evaluations,
                 2u * vertices);
    HS_EXPECT_EQ(Spacetime::geometry_counters.midpoint_evaluations, 0u);
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_draw_setups, 1u);
  }
};

inline void RingShowerWhiteBox::check_temporal_rings() {
  enum class Scenario {
    YOUNG_GROWTH = 0,
    NO_HISTORY = 1,
    RESPAWN = 2,
    FINAL_FRAME_HOLD = 3,
    MATURE_GROWTH = 4,
    SHORT_LIFE_GROWTH = 5,
  };
  using FX = RingShower<DEFAULT_W, DEFAULT_H>;
  const ScopedTemporalSampling temporal_sampling(true);
  {
    reset_globals();
    FX effect;
    effect.init();
    auto &ring = effect.rings[0];
    ring.normal = math::Y_AXIS;
    ring.life = 64;
    ring.age = 30;
    ring.previous_age = 29;
    ring.has_history = true;
    ring.palette.rebake(dot_keyed(
        ProceduralPalette({1, 1, 1}, {0, 0, 0}, {1, 1, 1}, {0, 0, 0})));
    {
      Canvas canvas(effect);
      effect.draw_ring(canvas, ring);
    }
    const auto pixels = RingTemporalOracle::pixels(effect);
    expect_arm_symmetric(pixels);
  }
  for (const auto &[scenario, label] : {
           std::pair{Scenario::YOUNG_GROWTH, "young growth"},
           std::pair{Scenario::NO_HISTORY, "no history"},
           std::pair{Scenario::RESPAWN, "respawn"},
           std::pair{Scenario::FINAL_FRAME_HOLD, "final frame hold"},
           std::pair{Scenario::MATURE_GROWTH, "mature growth"},
           std::pair{Scenario::SHORT_LIFE_GROWTH, "short life growth"},
       })
    for (float direction : {-1.0f, 1.0f}) {
      if (!RingTemporalOracle::selected(static_cast<int>(scenario)))
        continue;
      auto render = [=](bool timed, float t, RingTemporalOracle::Offset rows,
                        int segment) {
        reset_globals();
        hs::random().seed(137);
        FX effect;
        effect.init();
        RingTemporalOracle::clip(effect, segment);
        auto &ring = effect.rings[0];
        ring.life = 64;
        ring.age = scenario == Scenario::YOUNG_GROWTH ? 3 : 32;
        ring.previous_age = ring.age - 1;
        if (scenario == Scenario::SHORT_LIFE_GROWTH) {
          ring.life = 8;
          ring.age = 3;
          ring.previous_age = 2;
        }
        ring.has_history = scenario != Scenario::NO_HISTORY;
        if (scenario == Scenario::RESPAWN) {
          ring.age = ring.life;
          effect.spawn_ring();
          ring.life = 64;
        }
        if (scenario == Scenario::FINAL_FRAME_HOLD)
          ring.age = ring.life - 1;
        ring.normal = math::Vector(direction, 0.7f, 0.02f).normalized();
        const float start = scenario == Scenario::NO_HISTORY ||
                                    scenario == Scenario::RESPAWN ||
                                    scenario == Scenario::FINAL_FRAME_HOLD
                                ? float(ring.age)
                                : float(ring.previous_age);
        const float age = start + (float(ring.age) - start) * t;
        const int stride = Plot::Ring::timed_stride<DEFAULT_W>(
            2.0f * (start + 1.0f) / ring.life,
            2.0f * (ring.age + 1.0f) / ring.life);
        const auto basis = math::make_basis(math::Quaternion(), ring.normal);
        {
          Canvas canvas(effect);
          if (timed) {
            Spacetime::geometry_counters = {};
            effect.draw_ring(canvas, ring);
            RingTemporalOracle::counts(DEFAULT_W / stride);
          } else {
            auto vertex = [&](int i, float time) {
              const float at = start + (float(ring.age) - start) * time;
              return RingTemporalOracle::vertex(
                  basis, i, double(math::PI_F) * (at + 1.0) / ring.life);
            };
            auto shade = [&](const math::Vector &v, Fragment &f) {
              f.color = ring.palette.get(dot_key(v.x));
              f.color.alpha *=
                  hs::clamp((age + 1.0f) / 4, 0.0f, 1.0f) * effect.params.alpha;
            };
            RingTemporalOracle::frozen(canvas, vertex, shade, rows, t, stride);
          }
        }
        return RingTemporalOracle::pixels(effect);
      };
      char name[80];
      std::snprintf(name, sizeof(name),
                    "RingShower scenario=%s(%d) direction=%.0f", label,
                    static_cast<int>(scenario), direction);
      RingTemporalOracle::compare(render, name,
                                  scenario == Scenario::YOUNG_GROWTH ||
                                      scenario == Scenario::MATURE_GROWTH ||
                                      scenario == Scenario::SHORT_LIFE_GROWTH);
    }
}

inline void ThrustersWhiteBox::check_temporal_rings() {
  enum class Scenario {
    MAIN_INNER_RADIUS = 0,
    MAIN_NO_HISTORY = 1,
    MAIN_FIRE_HOLD = 2,
    MAIN_FRAME_WRAP = 3,
    THRUSTER_MOTION = 4,
    THRUSTER_BIRTH = 5,
  };
  const ScopedTemporalSampling temporal_sampling(true);
  for (bool main_ring : {false, true}) {
    reset_globals();
    FX effect;
    effect.init();
    effect.ring_vec = math::Y_AXIS;
    effect.params.radius = 0.9f;
    effect.palette.rebake(dot_keyed(
        ProceduralPalette({1, 1, 1}, {0, 0, 0}, {1, 1, 1}, {0, 0, 0})));
    Slot slot;
    slot.reset(effect.orientation, math::Y_AXIS);
    slot.age = 2;
    slot.previous_age = 1;
    slot.has_history = true;
    {
      Canvas canvas(effect);
      if (main_ring)
        effect.draw_ring(canvas, 1.0f);
      else
        effect.draw_thruster(canvas, slot, slot.radius_at(), 0.5f);
    }
    const auto pixels = RingTemporalOracle::pixels(effect);
    expect_arm_symmetric(pixels);
  }
  {
    reset_globals();
    FX effect;
    effect.init();
    effect.on_fire_thruster();
    for (int frame = 0; frame < 3; ++frame) {
      const auto pose = effect.orientation.get();
      const float amplitude = effect.amplitude;
      const float phase = effect.warp_phase;
      effect.draw_frame();
      effect.advance_display();
      HS_EXPECT_EQ(effect.previous_ring.pose, pose);
      HS_EXPECT_EQ(effect.previous_ring.amplitude, amplitude);
      HS_EXPECT_EQ(effect.previous_ring.phase, phase);
      HS_EXPECT_EQ(effect.previous_ring.frame, frame + 1);
    }
  }
  for (const auto &[scenario, label] : {
           std::pair{Scenario::MAIN_INNER_RADIUS, "main inner radius"},
           std::pair{Scenario::MAIN_NO_HISTORY, "main no history"},
           std::pair{Scenario::MAIN_FIRE_HOLD, "main fire hold"},
           std::pair{Scenario::MAIN_FRAME_WRAP, "main frame wrap"},
           std::pair{Scenario::THRUSTER_MOTION, "thruster motion"},
           std::pair{Scenario::THRUSTER_BIRTH, "thruster birth"},
       })
    for (float direction : {-1.0f, 1.0f}) {
      if (!RingTemporalOracle::selected(static_cast<int>(scenario)))
        continue;
      const bool main_ring = scenario != Scenario::THRUSTER_MOTION &&
                             scenario != Scenario::THRUSTER_BIRTH;
      auto render = [=](bool timed, float t, RingTemporalOracle::Offset rows,
                        int segment) {
        reset_globals();
        hs::random().seed(137);
        FX effect;
        effect.init();
        effect.params.radius = scenario == Scenario::MAIN_INNER_RADIUS ? 0.63f
                               : scenario == Scenario::MAIN_NO_HISTORY ? 1.37f
                                                                       : 1.0f;
        RingTemporalOracle::clip(effect, segment);
        const auto first = math::make_rotation(math::Y_AXIS, 0.13f);
        const auto last =
            math::make_rotation(math::Y_AXIS, direction * 0.04f) * first;
        effect.orientation.set(last);
        effect.amplitude = 0.45f;
        effect.warp_phase = 0.27f;
        effect.t_global = scenario == Scenario::MAIN_FRAME_WRAP ? 0 : 5;
        effect.previous_ring = {first, 0.48f, 0.27f,
                                scenario == Scenario::MAIN_FRAME_WRAP ? 31 : 4};
        effect.has_ring_history = scenario != Scenario::MAIN_NO_HISTORY;
        if (scenario == Scenario::MAIN_FIRE_HOLD)
          effect.on_fire_thruster();
        const auto current = typename FX::RingEndpoint{
            last, effect.amplitude, effect.warp_phase, effect.t_global};
        auto start = effect.has_ring_history ? effect.previous_ring : current;
        if (scenario == Scenario::MAIN_FIRE_HOLD) {
          start.amplitude = current.amplitude;
          start.phase = current.phase;
          start.frame = current.frame;
        }
        Slot slot;
        slot.reset(effect.orientation,
                   math::Vector(direction, 0.7f, 0.02f).normalized());
        slot.age = scenario == Scenario::THRUSTER_BIRTH ? 0 : 2;
        slot.previous_age = 1;
        slot.has_history = scenario != Scenario::THRUSTER_BIRTH;
        const float age =
            scenario == Scenario::THRUSTER_BIRTH ? 0.0f : 1.0f + t;
        const int stride =
            main_ring
                ? 1
                : Plot::Ring::timed_stride<DEFAULT_W>(
                      0.3f *
                          (scenario == Scenario::THRUSTER_BIRTH ? 1.0f : 2.0f) /
                          8.0f,
                      0.3f * (slot.age + 1.0f) / 8.0f);
        {
          Canvas canvas(effect);
          if (timed) {
            Spacetime::geometry_counters = {};
            if (main_ring)
              effect.draw_ring(canvas, 1.0f);
            else
              effect.draw_thruster(
                  canvas, slot, slot.radius_at(),
                  1.0f - math::ease_out_expo(float(slot.age) / Slot::LIFE));
            RingTemporalOracle::counts(DEFAULT_W / stride);
          } else if (main_ring) {
            auto basis = math::make_basis(math::Quaternion(), effect.ring_vec);
            if (effect.params.radius > 1.0f) {
              basis.u = -basis.u;
              basis.v = -basis.v;
            }
            auto warp = [](int i, const auto &state) {
              const double spatial =
                  std::sin(4.0 * double(math::PI_F) * i / DEFAULT_W -
                           double(math::PI_F) / 2 + state.phase);
              const double temporal =
                  std::sin(6.0 * double(math::PI_F) * state.frame / 32 -
                           double(math::PI_F) / 2);
              return spatial * temporal * state.amplitude;
            };
            auto vertex = [&](int i, float time) {
              const double polar =
                  std::min(double(effect.params.radius),
                           2.0 - double(effect.params.radius)) *
                      double(math::PI_F) / 2 +
                  (1 - time) * warp(i, start) + time * warp(i, current);
              return math::rotate(RingTemporalOracle::vertex(basis, i, polar),
                                  math::slerp(start.pose, last, time));
            };
            auto shade = [&](const math::Vector &v, Fragment &f) {
              f.color = Color4(effect.palette.get_color(dot_key(v.x)),
                               effect.params.alpha);
            };
            RingTemporalOracle::frozen(canvas, vertex, shade, rows, t, stride);
          } else {
            const auto basis =
                math::make_basis(slot.orientation.get(), slot.point);
            auto vertex = [&](int i, float time) {
              const float at =
                  scenario == Scenario::THRUSTER_BIRTH ? 0.0f : 1.0f + time;
              return RingTemporalOracle::vertex(
                  basis, i, 0.3 * (at + 1.0) / 8 * double(math::PI_F) / 2);
            };
            auto shade = [&](const math::Vector &, Fragment &f) {
              const float opacity = 1.0f - math::ease_out_expo(age / 16);
              f.color = Color4(CRGB(255, 255, 255));
              f.color.color = f.color.color * opacity;
              f.color.alpha = opacity * effect.params.alpha;
            };
            RingTemporalOracle::frozen(canvas, vertex, shade, rows, t, stride);
          }
        }
        return RingTemporalOracle::pixels(effect);
      };
      char name[80];
      std::snprintf(name, sizeof(name),
                    "Thrusters scenario=%s(%d) direction=%.0f", label,
                    static_cast<int>(scenario), direction);
      RingTemporalOracle::compare(render, name,
                                  scenario != Scenario::MAIN_NO_HISTORY &&
                                      scenario != Scenario::THRUSTER_BIRTH);
    }
}
