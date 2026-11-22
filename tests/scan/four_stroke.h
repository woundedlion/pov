/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke.h
 * @brief Forward boundary-register and material conformance tests. */

#include "tests/scan/context.h"

namespace hs_test::scan_tests {

template <typename Arrow>
concept BindableFourStrokeArrow = requires(Arrow &&arrow) {
  Scan::four_stroke_tangent_direction<96, 20>(std::forward<Arrow>(arrow));
};
using FourStrokeTestArrow =
    decltype([](const math::Vector &, uint64_t) { return math::X_AXIS; });
static_assert(BindableFourStrokeArrow<FourStrokeTestArrow &>);
static_assert(!BindableFourStrokeArrow<FourStrokeTestArrow>);

/** @brief Both public rasterizers expose signed contour data after the base shader. */
inline void test_four_stroke_fragment_registers() {
  constexpr int W = 96, H = 20;
  const math::Basis basis = math::make_basis(math::Quaternion{}, math::X_AXIS);
  for (bool plot : {false, true}) {
    int base_calls = 0, edge_calls = 0;
    bool leading = false, trailing = false, interior = false,
         registers_valid = true;
    bool north = false, south = false;
    const auto base = [&](const math::Vector &, Fragment &f) {
      ++base_calls;
      f.v0 = 0.375f;
      f.age = 0.25f;
      f.color = Color4(Pixel(0, 65535, 0), 0.5f);
    };
    const auto inspect = [&](const math::Vector &p, Fragment &f) {
      ++edge_calls;
      registers_valid &= f.v0 == 0.375f && f.v1 <= 0 && f.v2 == 1 &&
                         std::abs(f.v3) <= 1.0001f && f.age == 0.25f &&
                         f.color.alpha == 0.5f && f.color.color.g == 65535;
      if (f.v1 <= -1)
        interior = true;
      else {
        leading |= f.v3 > 0.25f;
        trailing |= f.v3 < -0.25f;
      }
      if (plot && f.v1 > -0.8f && p.x < -0.93f) {
        if (p.y > 0.15f) {
          north = true;
          registers_valid &= f.v3 < -0.8f;
        }
        if (p.y < -0.15f) {
          south = true;
          registers_valid &= f.v3 > 0.8f;
        }
      }
      f.color.color = Pixel(65535, 0, 0);
    };
    hs_test::StubEffect fx(W, H);
    {
      Canvas canvas(fx);
      Pipeline<W, H> filters;

      Scan::FourStrokeMaterial material{
          {.heading = plot ? 90.f : 37.f}, base, inspect};
      if (plot) {
        Fragment a, b;
        a.pos =
            math::Vector(math::Spherical(0.75f * math::PI_F, math::PI_F / 2));
        b.pos =
            math::Vector(math::Spherical(1.25f * math::PI_F, math::PI_F / 2));
        Plot::Line::draw<W, H>(filters, canvas, a, b, 0.36f, material);
      } else {
        Scan::Star::draw<W, H>(filters, canvas, basis, 0.46f, 5, material);
      }
    }
    fx.advance_display();
    HS_EXPECT_TRUE(base_calls > 0 && edge_calls == base_calls);
    HS_EXPECT_TRUE(registers_valid);
    HS_EXPECT_TRUE(!plot || (north && south));
    HS_EXPECT_TRUE(leading && trailing && interior);
    uint16_t peak = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const auto c = fx.get_pixel(x, y);
        peak = std::max(peak, c.r);
        HS_EXPECT_EQ(c.g, 0);
      }
    HS_EXPECT_NEAR(float(peak), 32768.f, 2.f);
  }
}

/** @brief A geodesic's end cap keeps its endpoint material across the angle seam. */
inline void test_four_stroke_line_cap_material_across_seam() {
  constexpr int W = 96, H = 48;
  constexpr float DEG = math::PI_F / 180.0f;
  StubEffect effect(W, H);
  Canvas canvas(effect);
  Pipeline<W, H> pipeline;
  Fragment start, end;
  start.pos = math::X_AXIS;
  end.pos = {cosf(175.0f * DEG), 0.0f, sinf(175.0f * DEG)};
  start.color = Color4(Pixel(65535, 0, 0), 1.0f);
  end.color = Color4(Pixel(0, 0, 65535), 1.0f);
  int sampled = 0;
  const auto base = [&](const math::Vector &p, Fragment &f) {
    if (p.x < -0.98f && p.z < -0.06f && fabsf(p.y) < 0.1f) {
      ++sampled;
      HS_EXPECT_EQ(f.v0, 1.0f);
      HS_EXPECT_EQ(f.color.color.b, 65535);
    }
  };
  Plot::Line::draw<W, H>(pipeline, canvas, start, end, 15.0f * DEG,
                         Scan::FourStrokeMaterial{{.strength = 0.0f}, base});
  HS_EXPECT_GT(sampled, 0);
}

/** @brief Disabled four-stroke stars render the ordinary scan bit for bit. */
inline void test_four_stroke_disabled_star() {
  constexpr int W = 96, H = 20;
  const auto base = [](const math::Vector &, Fragment &f) {
    f.color = Color4(179, 62, 93, 0.65f);
  };
  const auto draw = [&](bool wrapped) {
    return render_stub_frame<W, H>(
        [](auto &) {},
        [&](auto &sink, Canvas &canvas) {
          const auto basis = math::make_basis(math::Quaternion{}, math::X_AXIS);
          if (wrapped)
            Scan::Star::draw<W, H>(
                sink, canvas, basis, 0.38f, 5,
                Scan::FourStrokeMaterial{{.strength = 0}, base});
          else
            Scan::Star::draw<W, H>(sink, canvas, basis, 0.38f, 5, base);
        });
  };
  const auto expected = draw(false), actual = draw(true);
  bool lit = false;
  for (size_t i = 0; i < expected.size(); ++i) {
    lit |= expected[i].r != 0;
    HS_EXPECT_EQ(expected[i].r, actual[i].r);
    HS_EXPECT_EQ(expected[i].g, actual[i].g);
    HS_EXPECT_EQ(expected[i].b, actual[i].b);
  }
  HS_EXPECT_TRUE(lit);
}
/** @brief Projected chords track the exact curved star perimeter in output pixels. */
inline void test_four_stroke_star_contour_error() {
  constexpr int W = 96, H = 20;
  float worst = 0;
  int projected = 0;
  for (float latitude : {0.f, 0.5f, 1.f, math::PI_F / 2}) {
    const math::Vector axis{std::cos(latitude), std::sin(latitude), 0};
    const auto basis = math::make_basis(math::Quaternion{}, axis);
    constexpr float RADIUS = 0.31f, PHASE = 0.37f;
    const SDF::Star star(basis, RADIUS, 5, PHASE);
    Scan::ProjectedStar<W, H> contour;
    const bool ok = contour.prepare(star);
    if (latitude == math::PI_F / 2)
      HS_EXPECT_FALSE(ok);
    projected += ok;
    if (!ok)
      continue;
    for (int edge = 0; edge < 10; ++edge) {
      const auto vertex = [](int i) {
        const float a = i * math::PI_F / 5 - PHASE;
        const float r =
            RADIUS * (math::PI_F / 2) * (i % 2 ? Render::STAR_INNER_RATIO : 1);
        return math::Vector{r * std::cos(a), 0, r * std::sin(a)};
      };
      const auto a = vertex(edge), b = vertex((edge + 1) % 10);
      for (int k = 0; k <= 64; ++k) {
        const auto t = a + (b - a) * (k / 64.f);
        const float radius = t.length();
        const auto direction =
            basis.v * std::cos(radius) +
            (basis.u * t.x + basis.w * t.z) * (std::sin(radius) / radius);
        const float x =
            std::atan2(direction.z, direction.x) * (W / math::TWO_PI_F);
        const float y = math::DisplayGeometry<H>::phi_to_row(
            std::acos(std::clamp(direction.y, -1.f, 1.f)));
        const auto boundary = contour.query(
            contour.unwrap(x), y, [](FourStroke::Direction n) { return n.x; });
        worst = std::max(worst, boundary.depth);
      }
    }
  }
  std::printf("  projected star contour max error %.6f pixels\n",
              double(worst));
  HS_EXPECT_GE(projected, 2);
  HS_EXPECT_TRUE(worst < 0.025f);
}
/** @brief Retained contours shade only covered star taps and invalidate on a pole-crossing outline. */
inline void test_four_stroke_prepared_star() {
  constexpr int W = 96, H = 20;
  Scan::ProjectedStar<W, H> outline;
  const auto equator_basis = math::make_basis(math::Quaternion{}, math::X_AXIS);
  const auto pole_basis = math::make_basis(math::Quaternion{}, math::Y_AXIS);
  const SDF::Star star(equator_basis, 0.31f, 5, 0.37f);
  const SDF::Star pole(pole_basis, 0.31f, 5, 0.37f);
  int calls = 0;
  const auto base = [&](const auto &, Fragment &f) {
    ++calls;
    f.color = Color4(243, 173, 95, 0.7f);
  };
  const Scan::FourStrokeMaterial material{{}, base};
  const auto draw = [&](const SDF::Star &shape, bool cached) {
    return render_stub_frame<W, H>(
        [](auto &) {},
        [&](auto &sink, Canvas &canvas) {
          if (cached)
            Scan::Star::draw<W, H>(sink, canvas, shape, outline, material);
          else
            Scan::rasterize<W, H>(sink, canvas, shape, material);
        });
  };
  for (const auto *shape : {&star, &pole}) {
    const bool expected = shape == &star;
    HS_EXPECT_EQ(outline.prepare(*shape), expected);
    HS_EXPECT_EQ(outline.ready(), expected);
    calls = 0;
    const auto cached = draw(*shape, true);
    const int cached_calls = calls;
    HS_EXPECT_TRUE(cached_calls > 0 && cached_calls < W * H * 16 / 4);
    if (expected)
      HS_EXPECT_TRUE(cached_calls <
                     16 * (outline.end_column() - outline.begin_column()) *
                         (outline.end_row() - outline.begin_row()));
    calls = 0;
    const auto fresh = draw(*shape, false);
    HS_EXPECT_EQ(cached_calls, calls);
    for (size_t i = 0; i < fresh.size(); ++i) {
      HS_EXPECT_EQ(fresh[i].r, cached[i].r);
      HS_EXPECT_EQ(fresh[i].g, cached[i].g);
      HS_EXPECT_EQ(fresh[i].b, cached[i].b);
    }
  }
}
/** @brief Counts SDF probes in the sampled fallback. */
struct CountedStrokeStar : SDF::Star {
  using SDF::Star::Star;
  mutable int probes = 0;
  template <bool UV>
  void distance(const math::Vector &at, SDF::DistanceResult &hit) const {
    if constexpr (!UV)
      ++probes;
    SDF::Star::distance<UV>(at, hit);
  }
};
/** @brief Same field with conservative bounds covering the entire frame. */
struct UnboundedStrokeStar : CountedStrokeStar {
  using CountedStrokeStar::CountedStrokeStar;
  void swept(float) const = delete;
  template <int H> auto get_vertical_bounds() const {
    auto bounds = SDF::Star::get_vertical_bounds<H>();
    bounds.y_min = 0;
    bounds.y_max = H - 1;
    return bounds;
  }
};
/** @brief The sampled scan preserves pole and seam coverage with fewer field probes. */
inline void test_four_stroke_sampled_bounds() {
  constexpr int W = 96, H = 20;
  auto base = [](const auto &, Fragment &f) {
    f.color = Color4(241, 173, 95, 0.8f);
  };
  const Scan::FourStrokeMaterial material{{}, base};
  const auto draw = [&](const auto &shape) {
    return render_stub_frame<W, H>([](auto &) {},
                                   [&](auto &sink, Canvas &canvas) {
                                     Scan::rasterize<W, H>(sink, canvas, shape,
                                                           material);
                                   });
  };
  for (const auto axis : {math::X_AXIS, math::Y_AXIS}) {
    const auto basis = math::make_basis(math::Quaternion{}, axis);
    CountedStrokeStar bounded(basis, 0.31f, 5, 0.37f);
    UnboundedStrokeStar brute(basis, 0.31f, 5, 0.37f);
    const auto expected = draw(brute), actual = draw(bounded);
    std::printf("  sampled star probes %d bounded / %d full\n", bounded.probes,
                brute.probes);
    HS_EXPECT_TRUE(bounded.probes > 0 && bounded.probes < brute.probes);
    for (size_t i = 0; i < actual.size(); ++i) {
      HS_EXPECT_EQ(expected[i].r, actual[i].r);
      HS_EXPECT_EQ(expected[i].g, actual[i].g);
      HS_EXPECT_EQ(expected[i].b, actual[i].b);
    }
  }
}
} // namespace hs_test::scan_tests
