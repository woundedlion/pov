/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke.h
 * @brief Forward boundary-register and material conformance tests. */

#include "tests/scan/context.h"

namespace hs_test::scan_tests {

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
} // namespace hs_test::scan_tests
