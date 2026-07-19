/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Searches (sides, rho, phi) for near-polar faces whose AA fringe the constant
 * azimuth pad clips, to seed the polar regression case. Not a CTest.
 */
#include <cstdio>
#include <vector>

#include "core/engine/memory.h"
#include "core/render/sdf.h"
#include "tests/aa_audit.h"
#include "tests/test_fixture.h"

namespace {
constexpr int W = 256, H = 128;
constexpr int HV = H + hs::H_OFFSET;

/** @brief Counts fringe pixels the emitted intervals never visit. */
int missed_for(int sides, float rho, float phi, float spin) {
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();
  Vector axis(sinf(phi), cosf(phi), 0.0f);
  Basis basis = make_basis(Quaternion(), axis);
  Vector v3[16];
  uint16_t idx[16];
  for (int i = 0; i < sides; ++i) {
    float a = (2.0f * PI_F * i) / sides + spin;
    v3[i] = (basis.v * cosf(rho) +
             (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(rho))
                .normalized();
    idx[i] = static_cast<uint16_t>(i);
  }
  SDF::FaceScratchBuffer scratch;
  SDF::Face face(std::span<const Vector>(v3, sides),
                 std::span<const uint16_t>(idx, sides), 0.0f, scratch, HV, H);
  if (face.y_min > face.y_max)
    return -1;

  std::vector<uint8_t> vis(static_cast<size_t>(W) * H, 0);
  int y_lo = face.y_min < 0 ? 0 : face.y_min;
  int y_hi = face.y_max > H - 1 ? H - 1 : face.y_max;
  bool handled = face.get_horizontal_intervals<W, H>(y_lo, [&](float f1,
                                                               float f2) {
    int x1 = static_cast<int>(f1), x2 = static_cast<int>(f2);
    for (int x = x1; x <= x2; ++x)
      for (int y = y_lo; y <= y_hi; ++y)
        vis[static_cast<size_t>(y) * W + ((x % W) + W) % W] = 1;
  });
  if (!handled)
    for (int y = y_lo; y <= y_hi; ++y)
      for (int x = 0; x < W; ++x)
        vis[static_cast<size_t>(y) * W + x] = 1;

  const float *ct = TrigLUT<W, H>::sin_theta.data() + W / 4;
  const float *st = TrigLUT<W, H>::sin_theta.data();
  const float pw = 2.0f * PI_F / W;
  int missed = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      Vector p(TrigLUT<W, H>::sin_phi[y] * ct[x], TrigLUT<W, H>::cos_phi[y],
               TrigLUT<W, H>::sin_phi[y] * st[x]);
      if (face.distance(p).dist < pw && !vis[static_cast<size_t>(y) * W + x])
        ++missed;
    }
  return missed;
}
} // namespace

int main() {
  hs_test::reset_globals();
  hs_aa::g_audit.legacy_pad = true;
  int best = 0;
  for (int sides = 3; sides <= 6; ++sides)
    for (int ri = 1; ri <= 20; ++ri)
      for (int pi = 1; pi <= 30; ++pi)
        for (int si = 0; si < 4; ++si) {
          float rho = 0.03f * ri, phi = 0.03f * pi;
          if (rho >= phi * 0.95f)
            continue;
          int m = missed_for(sides, rho, phi, 0.37f * si);
          if (m > 0) {
            best = m;
            std::printf("sides %d rho %.2f phi %.2f spin %.2f -> missed %d\n",
                        sides, rho, phi, 0.37f * si, m);
          }
        }
  return 0;
}
