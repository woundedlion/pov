/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief MobiusRings white-box cases. */

/** @brief White-box accessor for MobiusRings' singular numeric branches. */
struct MobiusRingsWhiteBox {
  using MR = MobiusRings<DEFAULT_W, DEFAULT_H>;
  /** @brief Independent double-precision homogeneous plane translation.
   * @param v Source unit vector.
   * @param real,imag Translation coefficients.
   * @return Warped direction. */
  static math::Vector translated(const math::Vector &v, double real,
                                 double imag) {
    double x = v.x, z = v.z, scale = 1.0 - v.y;
    if (x * x + z * z < math::MOBIUS_POLE_EPS &&
        scale < math::MOBIUS_POLE_EPS) {
      x = 1;
      z = 0;
      scale = 0;
    }
    const double nr = x + real * scale, ni = z + imag * scale;
    const double n2 = nr * nr + ni * ni, m2 = scale * scale;
    return {float(2 * nr * scale / (n2 + m2)), float((n2 - m2) / (n2 + m2)),
            float(2 * ni * scale / (n2 + m2))};
  }

  /** @brief Renders the production warp or independently frozen vertices.
   * @tparam W,H Canvas dimensions.
   * @param timed Whether to invoke the production moving draw.
   * @param t Frozen frame fraction.
   * @param displacement Reference shift in pixels.
   * @param direction Signed motion direction.
   * @param clipped Whether to use the device segment clip.
   * @param hold Whether history starts at the current endpoint.
   * @param wrap Whether latitude scrolling crosses its periodic cut.
   * @param count_change Positive births or negative retires latitude curves.
   * @return Linear framebuffer luminance. */
  template <int W, int H>
  static std::vector<float> render(bool timed, float t, float displacement,
                                   float direction, bool clipped = false,
                                   bool hold = false, bool wrap = false,
                                   int count_change = 0) {
    reset_effect_globals();
    Spacetime::geometry_counters = {};
    hs::random().seed(1337u);
    pin_frame_clock(0);
    using FX = MobiusRings<W, H>;
    FX fx;
    fx.init();
    fx.params.num_rings = 2.5f;
    fx.params.num_lines = 2.5f;
    fx.params.alpha = 0.45f;
    typename FX::Endpoint start, end;
    start.warp.b = math::Complex(0.7f, 0.15f);
    end.warp.b =
        math::Complex(0.7f + direction * 0.025f, 0.15f - direction * 0.025f);
    start.rotation = math::make_rotation(math::Z_AXIS, -0.13f);
    end.rotation =
        math::make_rotation(math::Z_AXIS, -0.13f + direction * 0.025f);
    start.rings = end.rings = 2.5f;
    start.lines = end.lines = 2.5f;
    start.phase = wrap ? (direction > 0 ? 0.995f : 0.005f) : 0.13f;
    end.phase = wrap ? (direction > 0 ? 0.005f : 0.995f) : 0.14f;
    start.alpha = 0.4f;
    end.alpha = 0.45f;
    if (count_change != 0) {
      start.phase = end.phase;
      if (count_change > 0)
        start.rings = 0.0f;
      else
        end.rings = 0.0f;
    }
    if (hold)
      start = end;
    const auto first_pose = math::make_rotation(math::Y_AXIS, 0.017f);
    const auto last_pose =
        hold ? first_pose
             : math::make_rotation(math::Y_AXIS, 0.017f + direction * 0.02f);
    if (clipped) {
      const auto map = pov::segment_map(1, H * 2, 4);
      const auto clip = pov::segment_clip(map, true, H * 2, 4, W);
      fx.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
    }
    if (timed) {
      fx.orientation.set(first_pose);
      fx.orientation.push(last_pose);
      Canvas canvas(fx);
      fx.draw_timed_grid(canvas, start, end);
    } else {
      displacement *= 0.70710678118f;
      int column = int(t * (W / 2) + 0.5f) % (W / 2);
      const auto meridian =
          math::pixel_to_vector<W, H>(float(column), H / 2.0f);
      const auto across = math::cross(math::Y_AXIS, meridian).normalized();
      fx.orientation.set(
          math::make_rotation(across, displacement * math::RADIANS_PER_ROW<H>) *
          math::make_rotation(math::Y_AXIS,
                              displacement * math::RADIANS_PER_COLUMN<W>) *
          math::slerp(first_pose, last_pose, t));
      const double real =
          double(start.warp.b.re) + t * double(end.warp.b.re - start.warp.b.re);
      const double imag =
          double(start.warp.b.im) + t * double(end.warp.b.im - start.warp.b.im);
      const auto counter = math::slerp(start.rotation, end.rotation, t);
      const auto north =
          math::rotate(translated(math::Y_AXIS, real, imag), counter)
              .normalized();
      const auto south =
          math::rotate(translated(-math::Y_AXIS, real, imag), counter)
              .normalized();
      const auto mask = [&](const math::Vector &v, const math::Vector &origin) {
        const double distance = std::acos(std::clamp(
            double(math::dot(v, origin)) / std::sqrt(double(math::dot(v, v))),
            -1.0, 1.0));
        const double x = std::clamp(distance / 1.2, 0.0, 1.0);
        return float(x * x * x * (x * (x * 6 - 15) + 10));
      };
      Canvas canvas(fx);
      for (int axis = 0; axis < 2; ++axis)
        for (int i = 0; i < 3; ++i) {
          ScratchScope scope(scratch_arena_a);
          math::Basis basis;
          float radius = 1;
          if (axis == 0) {
            const float coordinate = math::wrap(
                i / 2.5f + start.phase +
                    std::remainder(end.phase - start.phase, 1.0f) * t,
                1.0f);
            radius = (4 / math::PI_F) * atanf(expf(2.5f - coordinate * 5));
            basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
          } else {
            const float angle = i / 2.5f * math::PI_F;
            const math::Vector normal(cosf(angle), 0, -sinf(angle));
            basis = {math::cross(normal, math::Y_AXIS), normal, math::Y_AXIS};
          }
          Fragments source;
          source.bind(scratch_arena_a, FX::CURVE_SIDES + 2);
          Plot::Polygon<Plot::GeodesicProjection>::sample(source, basis, radius,
                                                          FX::CURVE_SIDES);
          auto *vertices =
              scratch_arena_a.allocate_n<math::Vector>(FX::CURVE_SIDES);
          ArenaVector<Plot::Mesh::Edge> edges;
          edges.bind(scratch_arena_a, FX::CURVE_SIDES);
          for (int k = 0; k < FX::CURVE_SIDES; ++k) {
            vertices[k] =
                math::rotate(translated(source[k].pos, real, imag), counter)
                    .normalized();
            edges.push_back({uint16_t(k), uint16_t((k + 1) % FX::CURVE_SIDES)});
          }
          struct View {
            std::span<const math::Vector> vertices;
          };
          auto shader = [&](const math::Vector &v, Fragment &f) {
            const float phase =
                start.phase + std::remainder(end.phase - start.phase, 1.0f) * t;
            const float coordinate =
                axis == 0
                    ? i / 2.5f
                    : FX::conformal_coord(
                          math::fast_sinf(((f.v2 + f.v0) / FX::CURVE_SIDES +
                                           0.5f / FX::CURVE_SIDES) *
                                          2 * math::PI_F),
                          phase);
            f.color = fx.baked_palette.get(coordinate);
            f.color.alpha *= std::min(1.0f, 2.5f - i) *
                             hs::lerp(start.alpha, end.alpha, t) *
                             mask(v, north) * mask(v, south);
            if (axis == 0 && count_change != 0)
              f.color.alpha *= count_change > 0 ? t : 1.0f - t;
          };
          Plot::GeometryCallback geometry;
          geometry.chained = true;
          geometry.rigid_lattice = true;
          Plot::Mesh::draw_column_edges<W, H>(fx.filters.next.next, canvas,
                                              View{{vertices, FX::CURVE_SIDES}},
                                              edges, shader, geometry, 0);
        }
    }
    fx.advance_display();
    return capture_luminance<W, H>(fx);
  }

  /** @brief Warps, pole fades and colour follow both arm cuts in both directions. */
  static void check_column_time_oracle() {
    for (float direction : {-1.0f, 1.0f}) {
      seam_compare_render<96, 20>(
          [&](bool timed, float t, float shift) {
            return render<96, 20>(timed, t, shift, direction);
          },
          "MobiusRings 96x20");
      seam_compare_render<288, 144>(
          [&](bool timed, float t, float shift) {
            return render<288, 144>(timed, t, shift, direction);
          },
          "MobiusRings 288x144");
    }
    for (float direction : {-1.0f, 1.0f}) {
      seam_compare_render<96, 20>(
          [&](bool timed, float t, float shift) {
            return render<96, 20>(timed, t, shift, direction, false, false,
                                  true);
          },
          "MobiusRings wrap 96x20");
      seam_compare_render<288, 144>(
          [&](bool timed, float t, float shift) {
            return render<288, 144>(timed, t, shift, direction, false, false,
                                    true);
          },
          "MobiusRings wrap 288x144");
    }
    for (int count_change : {-1, 1}) {
      seam_compare_render<96, 20>(
          [&](bool timed, float t, float shift) {
            return render<96, 20>(timed, t, shift, 1.0f, false, false, false,
                                  count_change);
          },
          "MobiusRings count change 96x20");
      seam_compare_render<288, 144>(
          [&](bool timed, float t, float shift) {
            return render<288, 144>(timed, t, shift, 1.0f, false, false, false,
                                    count_change);
          },
          "MobiusRings count change 288x144");
    }
    const auto full = render<288, 144>(true, 1.0f, 0.0f, 1.0f);
    const auto clipped = render<288, 144>(true, 1.0f, 0.0f, 1.0f, true);
    const auto work = Spacetime::geometry_counters;
    constexpr unsigned VERTICES = 6 * (288 / 4);
    std::printf(
        "  MobiusRings device clip: callback evaluations=%u vertices=%u ratio=%.3f columns=%u culled=%u\n",
        VERTICES + work.endpoint_evaluations + work.midpoint_evaluations,
        VERTICES,
        float(VERTICES + work.endpoint_evaluations +
              work.midpoint_evaluations) /
            VERTICES,
        work.columns_visited, work.column_edges_culled);
    HS_EXPECT_LE(work.endpoint_evaluations + work.midpoint_evaluations,
                 2 * VERTICES);
    const auto map = pov::segment_map(1, 288, 4);
    const auto clip = pov::segment_clip(map, true, 288, 4, 288);
    for (int y = clip.y0; y < clip.y1; ++y)
      for (int x = clip.x0; x < clip.x1; ++x)
        HS_EXPECT_NEAR(full[y * 288 + x], clipped[y * 288 + x], 1e-5f);
    const auto first = render<96, 20>(true, 1.0f, 0.0f, 1.0f, false, true);
    const auto frozen = render<96, 20>(false, 0.0f, 0.0f, 1.0f, false, true);
    float error = 0;
    int lit = 0;
    for (size_t i = 0; i < first.size(); ++i) {
      error += fabsf(first[i] - frozen[i]);
      lit += first[i] > 0.0f;
    }
    std::printf("  MobiusRings first frame lit=%d\n", lit);
    HS_EXPECT_GT(lit, 350);
    HS_EXPECT_LE(error, 0.02f);
  }
  static float conformal_coord(float z, float phase) {
    return MR::conformal_coord(z, phase);
  }
  static math::Quaternion counter_rotation(const math::Vector &mid) {
    return MR::counter_rotation(mid);
  }
};

inline void test_mobius_rings_conformal_and_counter_rotation() {
  using WB = MobiusRingsWhiteBox;

  const float zs[] = {-1.0f, -0.999999f, -0.5f, 0.0f, 0.5f, 0.999999f, 1.0f};
  for (float z : zs) {
    for (float phase : {0.0f, 0.25f, 0.6f, 0.95f}) {
      float coord = WB::conformal_coord(z, phase);
      HS_EXPECT_TRUE(std::isfinite(coord));
      HS_EXPECT_GE(coord, 0.0f);
      HS_EXPECT_LE(coord, 1.0f);
    }
  }
  HS_EXPECT_NEAR(WB::conformal_coord(1.0f, 0.3f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::conformal_coord(-1.0f, 0.7f), 1.0f, 1e-6f);

  HS_EXPECT_EQ(WB::counter_rotation(math::Vector(0.0f, 0.0f, 0.0f)),
               math::Quaternion());

  math::Vector mid(0.3f, -0.7f, 0.4f);
  math::Vector r = math::rotate(mid.normalized(), WB::counter_rotation(mid));
  HS_EXPECT_VEC(r, math::Z_AXIS, 1e-3f);
}
