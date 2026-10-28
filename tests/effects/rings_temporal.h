/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file rings_temporal.h
 * @brief Frozen-column oracles for growing and warped rings. */

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
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    const char *filter = std::getenv("HS_RING_TEMPORAL_SCENARIO");
#pragma clang diagnostic pop
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
                     float fraction, int stride = 1) {
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
    struct View {
      std::span<const math::Vector> vertices;
    };
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
        pipeline, canvas, View{{points, size_t(count)}}, edges, shader,
        geometry, 0);
  }

  /** @brief Reads displayed linear luminance.
   * @param effect Displayed effect.
   * @return Row-major luminance. */
  template <typename FX> static std::vector<float> pixels(FX &effect) {
    effect.advance_display();
    std::vector<float> result(DEFAULT_W * DEFAULT_H);
    for (int y = 0; y < DEFAULT_H; ++y)
      for (int x = 0; x < DEFAULT_W; ++x) {
        const Pixel p = effect.get_pixel(x, y);
        result[y * DEFAULT_W + x] =
            (float(p.r) + float(p.g) + float(p.b)) / (3.0f * 65535.0f);
      }
    return result;
  }

  /** @brief Compares full and segmented draws with frozen columns.
   * @param render Timed or independently frozen renderer.
   * @param name Diagnostic label.
   * @param moving Whether the frame has continuous motion. */
  template <typename Render>
  static void compare(Render render, const char *name, bool moving) {
    const auto actual = render(true, 1.0f, 0.0f, -1);
    const auto endpoint = render(false, 1.0f, 0.0f, -1);
    float total = 0, endpoint_total = 0, energy = 0;
    int failures = 0;
    for (int c = 0; c < DEFAULT_W / 2; ++c) {
      const float t = float(c) / (DEFAULT_W / 2);
      const auto expected = render(false, t, 0.0f, -1);
      const auto shifted = render(false, t, 0.2f, -1);
      const auto reversed = render(false, t, -0.2f, -1);
      for (int x : {c, c + DEFAULT_W / 2}) {
        float error = 0, budget = 0, reverse_budget = 0;
        for (int y = 0; y < DEFAULT_H; ++y) {
          const size_t i = y * DEFAULT_W + x;
          error += fabsf(actual[i] - expected[i]);
          budget += fabsf(shifted[i] - expected[i]);
          reverse_budget += fabsf(reversed[i] - expected[i]);
          endpoint_total += fabsf(endpoint[i] - expected[i]);
          energy += expected[i];
        }
        total += error;
        budget = std::max(budget, reverse_budget);
        if (error > budget + 1e-5f) {
          for (Direction direction : {Direction::AZIMUTH, Direction::MERIDIAN})
            for (float sign : {-1.0f, 1.0f}) {
              const auto displaced =
                  render(false, t, Offset(sign * 0.2f, direction), -1);
              float axis_budget = 0;
              for (int y = 0; y < DEFAULT_H; ++y) {
                const size_t i = y * DEFAULT_W + x;
                axis_budget += fabsf(displaced[i] - expected[i]);
              }
              budget = std::max(budget, axis_budget);
            }
        }
        if (error > budget + 1e-5f) {
          if (failures < 4)
            std::printf("  %s col=%d error=%g budget=%g\n", name, x, error,
                        budget);
          ++failures;
        }
      }
    }
    std::printf("  %s error=%g endpoint=%g limit=%g bad_columns=%d\n", name,
                total, endpoint_total, endpoint_total / 20.0f, failures);
    HS_EXPECT_GT(energy, 0.01f);
    HS_EXPECT_EQ(failures, 0);
    if (moving) {
      HS_EXPECT_GT(endpoint_total, 0.01f);
      HS_EXPECT(
          hs_test::temporal_frame_error_within_budget(total, endpoint_total),
          "ring frame meets one twentieth endpoint error");
    }
    for (int segment = 0; segment < 8; ++segment) {
      const auto clipped = render(true, 1.0f, 0.0f, segment);
      const auto map = pov::segment_map(segment % 4, DEFAULT_H * 2, 4);
      const auto bounds =
          pov::segment_clip(map, segment < 4, DEFAULT_H * 2, 4, DEFAULT_W);
      for (int y = bounds.y0; y < bounds.y1; ++y)
        for (int x = bounds.x0; x < bounds.x1; ++x)
          HS_EXPECT_NEAR(clipped[y * DEFAULT_W + x], actual[y * DEFAULT_W + x],
                         1e-5f);
    }
  }

  /** @brief Applies the device's segment clip.
   * @param effect Destination effect.
   * @param segment Segment number, or -1 for a full frame. */
  template <typename FX> static void clip(FX &effect, int segment) {
    if (segment < 0)
      return;
    const auto map = pov::segment_map(segment % 4, DEFAULT_H * 2, 4);
    const auto bounds =
        pov::segment_clip(map, segment < 4, DEFAULT_H * 2, 4, DEFAULT_W);
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
  static void counts(size_t vertices = DEFAULT_W) {
    HS_EXPECT_EQ(Spacetime::geometry_counters.endpoint_evaluations,
                 2u * vertices);
    HS_EXPECT_EQ(Spacetime::geometry_counters.midpoint_evaluations, 0u);
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_draw_setups, 1u);
  }
};

inline void RingShowerWhiteBox::check_temporal_rings() {
  using FX = RingShower<DEFAULT_W, DEFAULT_H>;
  const ScopedTemporalSampling temporal_sampling(true);
  {
    reset_effect_globals();
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
    float energy = 0;
    for (int y = 0; y < DEFAULT_H; ++y)
      for (int x = 0; x < DEFAULT_W / 2; ++x) {
        energy += pixels[y * DEFAULT_W + x];
        HS_EXPECT_NEAR(pixels[y * DEFAULT_W + x],
                       pixels[y * DEFAULT_W + x + DEFAULT_W / 2], 3e-5f);
      }
    HS_EXPECT_GT(energy, 0.01f);
  }
  for (int scenario = 0; scenario < 6; ++scenario)
    for (float direction : {-1.0f, 1.0f}) {
      if (!RingTemporalOracle::selected(scenario))
        continue;
      auto render = [=](bool timed, float t, RingTemporalOracle::Offset rows,
                        int segment) {
        reset_effect_globals();
        hs::random().seed(137);
        FX effect;
        effect.init();
        RingTemporalOracle::clip(effect, segment);
        auto &ring = effect.rings[0];
        ring.life = 64;
        ring.age = scenario == 0 ? 3 : 32;
        ring.previous_age = ring.age - 1;
        if (scenario == 5) {
          ring.life = 8;
          ring.age = 3;
          ring.previous_age = 2;
        }
        ring.has_history = scenario != 1;
        if (scenario == 2) {
          ring.age = ring.life;
          effect.spawn_ring();
          ring.life = 64;
        }
        if (scenario == 3)
          ring.age = ring.life - 1;
        ring.normal = math::Vector(direction, 0.7f, 0.02f).normalized();
        const float start = scenario == 1 || scenario == 2 || scenario == 3
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
      std::snprintf(name, sizeof(name), "RingShower scenario=%d direction=%.0f",
                    scenario, direction);
      RingTemporalOracle::compare(render, name, scenario == 0 || scenario >= 4);
    }
}

inline void ThrustersWhiteBox::check_temporal_rings() {
  const ScopedTemporalSampling temporal_sampling(true);
  for (bool main : {false, true}) {
    reset_effect_globals();
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
      if (main)
        effect.draw_ring(canvas, 1.0f);
      else
        effect.draw_thruster(canvas, slot, slot.radius_at(), 0.5f);
    }
    const auto pixels = RingTemporalOracle::pixels(effect);
    float energy = 0;
    for (int y = 0; y < DEFAULT_H; ++y)
      for (int x = 0; x < DEFAULT_W / 2; ++x) {
        energy += pixels[y * DEFAULT_W + x];
        HS_EXPECT_NEAR(pixels[y * DEFAULT_W + x],
                       pixels[y * DEFAULT_W + x + DEFAULT_W / 2], 3e-5f);
      }
    HS_EXPECT_GT(energy, 0.01f);
  }
  {
    reset_effect_globals();
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
  for (int scenario = 0; scenario < 6; ++scenario)
    for (float direction : {-1.0f, 1.0f}) {
      if (!RingTemporalOracle::selected(scenario))
        continue;
      auto render = [=](bool timed, float t, RingTemporalOracle::Offset rows,
                        int segment) {
        reset_effect_globals();
        hs::random().seed(137);
        FX effect;
        effect.init();
        effect.params.radius = scenario == 0   ? 0.63f
                               : scenario == 1 ? 1.37f
                                               : 1.0f;
        RingTemporalOracle::clip(effect, segment);
        const auto first = math::make_rotation(math::Y_AXIS, 0.13f);
        const auto last =
            math::make_rotation(math::Y_AXIS, direction * 0.04f) * first;
        effect.orientation.set(last);
        effect.amplitude = 0.45f;
        effect.warp_phase = 0.27f;
        effect.t_global = scenario == 3 ? 0 : 5;
        effect.previous_ring = {first, 0.48f, 0.27f, scenario == 3 ? 31 : 4};
        effect.has_ring_history = scenario != 1;
        if (scenario == 2)
          effect.on_fire_thruster();
        const auto current = typename FX::RingEndpoint{
            last, effect.amplitude, effect.warp_phase, effect.t_global};
        auto start = effect.has_ring_history ? effect.previous_ring : current;
        if (scenario == 2) {
          start.amplitude = current.amplitude;
          start.phase = current.phase;
          start.frame = current.frame;
        }
        Slot slot;
        slot.reset(effect.orientation,
                   math::Vector(direction, 0.7f, 0.02f).normalized());
        slot.age = scenario == 5 ? 0 : 2;
        slot.previous_age = 1;
        slot.has_history = scenario != 5;
        const float age = scenario == 5 ? 0.0f : 1.0f + t;
        const int stride =
            scenario < 4 ? 1
                         : Plot::Ring::timed_stride<DEFAULT_W>(
                               0.3f * (scenario == 5 ? 1.0f : 2.0f) / 8.0f,
                               0.3f * (slot.age + 1.0f) / 8.0f);
        {
          Canvas canvas(effect);
          if (timed) {
            Spacetime::geometry_counters = {};
            if (scenario < 4)
              effect.draw_ring(canvas, 1.0f);
            else
              effect.draw_thruster(
                  canvas, slot, slot.radius_at(),
                  1.0f - math::ease_out_expo(float(slot.age) / Slot::LIFE));
            RingTemporalOracle::counts(DEFAULT_W / stride);
          } else if (scenario < 4) {
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
              const float at = scenario == 5 ? 0.0f : 1.0f + time;
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
      std::snprintf(name, sizeof(name), "Thrusters scenario=%d direction=%.0f",
                    scenario, direction);
      RingTemporalOracle::compare(render, name, scenario != 1 && scenario != 5);
    }
}
