/*
 * Airocean constants/transforms and Peirce coefficients derived from PROJ
 * commit 71eb06bfa597ed1e6d673723fc509e7604fe72af.
 *
 * Copyright (c) 2000, Frank Warmerdam
 * Copyright (c) 2019, 2020, 2021, 2022, 2023, 2024 PROJ contributors
 * Copyright (c) 2005, 2006, 2009 Gerald I. Evenden
 * Copyright (c) 2020 Kristian Evers
 * Copyright (c) 2021 Toby C Wilkinson
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
#pragma once

#include "core/math/3dmath.h"

namespace shadierball {

enum class ProjectionTrait : uint8_t {
  NONE = 0,
  CUT = 1U << 0,
  GLUED = 1U << 1,
  PERIODIC = 1U << 2,
  FOLDED = 1U << 3,
  SINGULAR = 1U << 4
};

constexpr uint8_t projection_traits(ProjectionTrait a) {
  return static_cast<uint8_t>(a);
}

constexpr uint8_t projection_traits(ProjectionTrait a, ProjectionTrait b,
                                    ProjectionTrait c = ProjectionTrait::NONE) {
  return projection_traits(a) | projection_traits(b) | projection_traits(c);
}

struct ProjectionKernelResult {
  Complex coords;
  uint8_t region_id;
  uint8_t component_id;
  uint8_t boundary_flags;
  float fade_edge_distance;
  uint8_t flags;
  uint8_t traits;
  uint8_t edge_class;

  constexpr ProjectionKernelResult(Complex coords, uint8_t region_id,
                                   uint8_t component_id, uint8_t boundary_flags,
                                   float fade_edge_distance, uint8_t flags,
                                   uint8_t traits = 0, uint8_t edge_class = 0)
      : coords(coords), region_id(region_id), component_id(component_id),
        boundary_flags(boundary_flags), fade_edge_distance(fade_edge_distance),
        flags(flags), traits(traits), edge_class(edge_class) {}
};

inline float wrap_longitude(float longitude) {
  float wrapped = fmodf(longitude + PI_F, TWO_PI_F);
  if (wrapped < 0.0f)
    wrapped += TWO_PI_F;
  return wrapped - PI_F;
}

inline ProjectionKernelResult bonne_projection(const Vector &v,
                                               float central_meridian,
                                               float standard_parallel) {
  const float longitude = wrap_longitude(atan2f(v.z, v.x) - central_meridian);
  const float latitude = asinf(hs::clamp(v.y, -1.0f, 1.0f));
  const float cut_distance = cosf(latitude) * (PI_F - fabsf(longitude));
  const float c = fabsf(standard_parallel) == 0.5f * PI_F
                      ? 0.0f
                      : cosf(standard_parallel) / sinf(standard_parallel);
  const float rho = c + standard_parallel - latitude;
  float e = 0.0f;
  if (fabsf(rho) > 1e-7f)
    e = longitude * cosf(latitude) / rho;
  const Complex coords(rho * sinf(e),
                       latitude - standard_parallel +
                           2.0f * rho * sinf(0.5f * e) * sinf(0.5f * e));
  return {coords,
          static_cast<uint8_t>(longitude < 0.0f),
          0,
          1,
          std::max(0.0f, cut_distance),
          0,
          projection_traits(ProjectionTrait::CUT),
          0};
}

inline float peirce_elliptic_integral(float phi) {
  constexpr float C0 = 2.19174570831038f;
  constexpr float C[] = {-8.58691003636495e-7f, 2.02692115653689e-7f,
                         3.12960480765314e-5f,  5.30394739921063e-5f,
                         -0.0012804644680613f,  -0.00575574836830288f,
                         0.0914203033408211f};
  float y = phi * (2.0f / PI_F);
  y = 2.0f * y * y - 1.0f;
  const float y2 = 2.0f * y;
  float d1 = 0.0f;
  float d2 = 0.0f;
  for (float coefficient : C) {
    const float previous = d1;
    d1 = y2 * d1 - d2 + coefficient;
    d2 = previous;
  }
  return phi * (y * d1 - d2 + 0.5f * C0);
}

inline float peirce_sector_longitude(const Vector &v, float central_meridian) {
  if (v.x == 0.0f && v.z == 0.0f)
    return 0.5f * PI_F;
  float longitude = wrap_longitude(atan2f(v.z, v.x) - central_meridian);
  constexpr float BOUNDARIES[] = {-0.75f * PI_F, -0.25f * PI_F, 0.25f * PI_F,
                                  0.75f * PI_F};
  constexpr float TIE_EPSILON = 2e-6f;
  for (float boundary : BOUNDARIES)
    if (fabsf(longitude - boundary) <= TIE_EPSILON)
      return boundary;
  return longitude;
}

inline ProjectionKernelResult peirce_projection(const Vector &v,
                                                float central_meridian,
                                                uint8_t layout, float scroll) {
  constexpr float INV_SQRT_TWO = 0.7071067811865475f;
  constexpr float K = 1.8540746773013719f;
  constexpr float SHIFT = 2.0f * K;
  const float longitude = peirce_sector_longitude(v, central_meridian);
  const float latitude = asinf(hs::clamp(v.y, -1.0f, 1.0f));
  const float sl = sinf(longitude);
  const float cl = cosf(longitude);
  const float cp = cosf(latitude);
  const float a = acosf(hs::clamp(cp * (sl + cl) * INV_SQRT_TWO, -1.0f, 1.0f));
  const float b = acosf(hs::clamp(cp * (sl - cl) * INV_SQRT_TWO, -1.0f, 1.0f));
  float m = asinf(sqrtf(std::max(0.0f, 1.0f + std::min(0.0f, cosf(a + b)))));
  float n = asinf(sqrtf(fabsf(1.0f - std::max(0.0f, cosf(a - b)))));
  if (sl < 0.0f)
    m = -m;
  if (cl > 0.0f)
    n = -n;
  float x = peirce_elliptic_integral(m);
  float y = peirce_elliptic_integral(n);
  uint8_t region = 0;
  uint8_t flags = 0;
  if (latitude < 0.0f) {
    flags = 1;
    if (longitude < -0.75f * PI_F || longitude >= 0.75f * PI_F)
      region = 1;
    else if (longitude < -0.25f * PI_F)
      region = 2;
    else if (longitude < 0.25f * PI_F)
      region = 3;
    else
      region = 4;
    if (layout <= 1) {
      if (longitude < -0.75f * PI_F) {
        y = SHIFT - y;
      } else if (longitude < -0.25f * PI_F) {
        x = -SHIFT - x;
      } else if (longitude < 0.25f * PI_F) {
        y = -SHIFT - y;
      } else if (longitude < 0.75f * PI_F) {
        x = SHIFT - x;
      } else {
        y = SHIFT - y;
      }
    }
  }
  if (layout == 1) {
    const float old_x = x;
    x = INV_SQRT_TWO * (x - y);
    y = INV_SQRT_TWO * (old_x + y);
  } else if (layout == 2) {
    if (latitude < 0.0f)
      x = SHIFT - x;
    x -= K;
    x += scroll * 2.0f * SHIFT;
    x = fmodf(x + SHIFT, 2.0f * SHIFT);
    if (x < 0.0f)
      x += 2.0f * SHIFT;
    x -= SHIFT;
  } else if (layout == 3) {
    if (latitude < 0.0f)
      y = SHIFT - y;
    y -= K;
    y += scroll * 2.0f * SHIFT;
    y = fmodf(y + SHIFT, 2.0f * SHIFT);
    if (y < 0.0f)
      y += 2.0f * SHIFT;
    y -= SHIFT;
  }
  const float rotated_x = cp * cosf(longitude);
  const float rotated_z = cp * sinf(longitude);
  const float edge = acosf(
      hs::clamp(std::max(fabsf(rotated_x), fabsf(rotated_z)), 0.0f, 1.0f));
  uint8_t edge_class = fabsf(rotated_x) >= fabsf(rotated_z)
                           ? static_cast<uint8_t>(rotated_x < 0.0f)
                           : static_cast<uint8_t>(2 + (rotated_z < 0.0f));
  if (layout == 2)
    edge_class = 4;
  else if (layout == 3)
    edge_class = 5;
  return {Complex(x, y),
          region,
          0,
          2,
          edge,
          flags,
          projection_traits(ProjectionTrait::GLUED, ProjectionTrait::FOLDED,
                            ProjectionTrait::PERIODIC) |
              projection_traits(ProjectionTrait::SINGULAR),
          edge_class};
}

struct AiroceanVector {
  float x;
  float y;
  float z;
};

struct AiroceanPoint {
  float x;
  float y;
};

static constexpr AiroceanVector AIROCEAN_FACES[23][3] = {
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {-0.5154559599f, -0.3817168983f, 0.7672009925f}},
    {{0.4201524267f, 0.0781452494f, 0.9040825506f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f}},
    {{0.9950094394f, -0.09134779528f, 0.04014717588f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{0.5154559599f, 0.3817168983f, -0.7672009925f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {0.5188367303f, 0.8354203804f, 0.1813318376f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f}},
    {{-0.5154559599f, -0.3817168983f, 0.7672009925f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f}},
    {{-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f}},
    {{0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.3557814025f, -0.8435800025f, 0.4022342266f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{0.5154559599f, 0.3817168983f, -0.7672009925f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.9950094394f, -0.09134779528f, 0.04014717588f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {-0.5188367303f, -0.8354203804f, -0.1813318376f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {0.4146822253f, -0.6559624054f, -0.6306758079f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.3879669146f, 0.3827173765f, -0.6531583886f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f}},
    {{-0.4201524267f, -0.0781452494f, -0.9040825506f},
     {0.5154559599f, 0.3817168983f, -0.7672009925f},
     {-0.3879669146f, 0.3827173765f, -0.6531583886f}},
    {{-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f}},
    {{-0.3557814025f, 0.8435800025f, -0.4022342266f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f}},
    {{-0.9950094394f, 0.09134779528f, -0.04014717588f},
     {-0.5884910224f, 0.5302967344f, 0.0627648018f},
     {-0.4146822253f, 0.6559624054f, 0.6306758079f}}};
static constexpr AiroceanVector AIROCEAN_CENTERS[23] = {
    {0.6446661988f, 0.2740726115f, 0.375187188f},
    {0.1747689772f, 0.5231760117f, 0.5720300654f},
    {-0.1699952529f, 0.1174635855f, 0.7673197837f},
    {0.08682595643f, -0.3823838838f, 0.6911725899f},
    {0.5903144229f, -0.2855941828f, 0.4488213177f},
    {0.6764340432f, 0.3752631611f, -0.1819073264f},
    {0.2261704292f, 0.6869057604f, -0.3293677939f},
    {-0.08387563251f, 0.7783209294f, 0.1365911396f},
    {-0.6417158749f, 0.1218644341f, 0.4525765415f},
    {-0.6764340432f, -0.3752631611f, 0.1819073264f},
    {-0.2261704292f, -0.6869057604f, 0.3293677939f},
    {0.08387563251f, -0.7783209294f, -0.1365911396f},
    {0.5884910224f, -0.5302967344f, -0.0627648018f},
    {0.6417158749f, -0.1218644341f, -0.4525765415f},
    {-0.5903144229f, 0.2855941828f, -0.4488213177f},
    {-0.6446661988f, -0.2740726115f, -0.375187188f},
    {-0.1747689772f, -0.5231760117f, -0.5720300654f},
    {0.1699952529f, -0.1174635855f, -0.7673197837f},
    {-0.0760974524f, 0.5360047591f, -0.6075312026f},
    {-0.09755446046f, 0.2287630085f, -0.7748139772f},
    {-0.6464272881f, 0.4884081774f, -0.1265388669f},
    {-0.4529848834f, 0.6766130474f, 0.09706879436f},
    {-0.6660608957f, 0.4258689784f, 0.2177644779f}};
static constexpr AiroceanVector AIROCEAN_NORMALS[23] = {
    {0.8112534709f, 0.3448953238f, 0.4721387736f},
    {0.2199307791f, 0.658369178f, 0.7198475379f},
    {-0.2139234835f, 0.147817183f, 0.9656017935f},
    {0.1092625279f, -0.4811951573f, 0.8697775121f},
    {0.7428567302f, -0.3593941678f, 0.5648005937f},
    {0.8512303986f, 0.4722343789f, -0.2289137389f},
    {0.284614807f, 0.8644080973f, -0.4144792552f},
    {-0.105549815f, 0.9794457296f, 0.171887461f},
    {-0.807540758f, 0.1533552486f, 0.5695261995f},
    {-0.8512303986f, -0.4722343789f, 0.2289137389f},
    {-0.284614807f, -0.8644080973f, 0.4144792552f},
    {0.105549815f, -0.9794457296f, -0.171887461f},
    {0.7405621474f, -0.6673299565f, -0.07898376463f},
    {0.807540758f, -0.1533552486f, -0.5695261995f},
    {-0.7428567302f, 0.3593941678f, -0.5648005937f},
    {-0.8112534709f, -0.3448953238f, -0.4721387736f},
    {-0.2199307791f, -0.658369178f, -0.7198475379f},
    {0.2139234835f, -0.147817183f, -0.9656017935f},
    {-0.1092625279f, 0.4811951573f, -0.8697775121f},
    {-0.1092625279f, 0.4811951573f, -0.8697775121f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f},
    {-0.7405621474f, 0.6673299565f, 0.07898376463f}};
static constexpr AiroceanPoint AIROCEAN_PLANAR_FACES[23][3] = {
    {{1.821185995f, 3.154386673f},
     {1.821185995f, 4.205848897f},
     {2.731778992f, 3.680117785f}},
    {{1.821185995f, 3.154386673f},
     {0.9105929973f, 3.680117785f},
     {1.821185995f, 4.205848897f}},
    {{1.821185995f, 3.154386673f},
     {0.9105929973f, 2.628655561f},
     {0.9105929973f, 3.680117785f}},
    {{1.821185995f, 3.154386673f},
     {1.821185995f, 2.102924448f},
     {0.9105929973f, 2.628655561f}},
    {{1.821185995f, 3.154386673f},
     {2.731778992f, 3.680117785f},
     {2.731778992f, 2.628655561f}},
    {{2.731778992f, 3.680117785f},
     {1.821185995f, 4.205848897f},
     {2.731778992f, 4.731580009f}},
    {{1.821185995f, 5.257311121f},
     {1.821185995f, 4.205848897f},
     {0.9105929973f, 4.731580009f}},
    {{0.9105929973f, 4.731580009f},
     {1.821185995f, 4.205848897f},
     {0.9105929973f, 3.680117785f}},
    {{0.9105929973f, 2.628655561f},
     {0.0f, 3.154386673f},
     {0.9105929973f, 3.680117785f}},
    {{0.9105929973f, 2.628655561f},
     {0.9105929973f, 1.577193336f},
     {0.0f, 2.102924448f}},
    {{0.9105929973f, 2.628655561f},
     {1.821185995f, 2.102924448f},
     {0.9105929973f, 1.577193336f}},
    {{0.9105929973f, 1.577193336f},
     {1.821185995f, 2.102924448f},
     {1.821185995f, 1.051462224f}},
    {{1.821185995f, 1.051462224f},
     {1.821185995f, 2.102924448f},
     {2.731778992f, 1.577193336f}},
    {{1.821185995f, 0.0f},
     {1.821185995f, 1.051462224f},
     {2.731778992f, 0.5257311121f}},
    {{0.0f, 5.257311121f}, {0.9105929973f, 4.731580009f}, {0.0f, 4.205848897f}},
    {{0.0f, 1.051462224f}, {0.0f, 2.102924448f}, {0.9105929973f, 1.577193336f}},
    {{0.9105929973f, 0.5257311121f},
     {0.9105929973f, 1.577193336f},
     {1.821185995f, 1.051462224f}},
    {{0.9105929973f, 0.5257311121f},
     {1.821185995f, 1.051462224f},
     {1.821185995f, 0.0f}},
    {{0.9105929973f, 4.731580009f},
     {0.4552964987f, 4.994445565f},
     {0.9105929973f, 5.783042233f}},
    {{0.9105929973f, 0.5257311121f},
     {1.821185995f, 0.0f},
     {0.9105929973f, 0.0f}},
    {{0.0f, 4.205848897f},
     {0.9105929973f, 4.731580009f},
     {0.6070619982f, 4.205848897f}},
    {{0.9105929973f, 4.731580009f},
     {0.9105929973f, 3.680117785f},
     {0.6070619982f, 4.205848897f}},
    {{0.0f, 3.154386673f},
     {0.3035309991f, 3.680117785f},
     {0.9105929973f, 3.680117785f}}};
static constexpr float AIROCEAN_TRANSFORMS[23][2][4] = {
    {{0.5771127853f, -0.6019490725f, -0.5519041105f, 2.124716994f},
     {0.09385435001f, 0.7202114479f, -0.6873767753f, 3.680117785f}},
    {{0.9709901201f, -0.2187361325f, -0.09660585362f, 1.517654996f},
     {0.09385435001f, 0.7202114479f, -0.6873767753f, 3.680117785f}},
    {{0.9721374115f, -0.06476823822f, 0.2252863255f, 1.214123996f},
     {0.09584151699f, 0.9868916636f, -0.1298431665f, 3.154386673f}},
    {{0.9921258754f, -0.001098710673f, -0.1252399307f, 1.517654996f},
     {0.061220482f, 0.876612807f, 0.4772861187f, 2.628655561f}},
    {{0.280304148f, -0.5991800397f, -0.7499419075f, 2.428247993f},
     {0.6079419899f, 0.7154153424f, -0.3443652491f, 3.154386673f}},
    {{0.2596066191f, -0.7580069591f, -0.5983559587f, 2.428247993f},
     {-0.4560824616f, 0.4499112594f, -0.7678337365f, 4.205848897f}},
    {{0.9586365701f, -0.2580860646f, 0.1200312867f, 1.517654996f},
     {-0.003215303703f, -0.4314976531f, -0.902108329f, 4.731580009f}},
    {{0.9928349404f, 0.09405868119f, 0.07370037669f, 1.214123996f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.5819727896f, 0.05026939416f, 0.8116530418f, 0.6070619982f},
     {0.09584151699f, 0.9868916636f, -0.1298431665f, 3.154386673f}},
    {{0.5247823075f, -0.7686380597f, 0.3657855423f, 0.6070619982f},
     {0.003215303703f, 0.4314976531f, 0.902108329f, 2.102924448f}},
    {{0.9586365701f, -0.2580860646f, 0.1200312867f, 1.214123996f},
     {0.003215303703f, 0.4314976531f, 0.902108329f, 2.102924448f}},
    {{0.9928349404f, 0.09405868119f, 0.07370037669f, 1.517654996f},
     {-0.05601801133f, -0.1784349382f, 0.9823558191f, 1.577193336f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 2.124716994f},
     {-0.05601801133f, -0.1784349382f, 0.9823558191f, 1.577193336f}},
    {{0.5819727896f, 0.05026939416f, 0.8116530418f, 2.124716994f},
     {-0.09584151699f, -0.9868916636f, 0.1298431665f, 0.5257311121f}},
    {{0.3863411333f, 0.9191578806f, 0.07674189989f, 0.3035309991f},
     {0.5467215079f, -0.1611974646f, -0.8216513678f, 4.731580009f}},
    {{0.2072761413f, -0.9246959463f, 0.3193336941f, 0.3035309991f},
     {-0.5467215079f, 0.1611974646f, 0.8216513678f, 1.577193336f}},
    {{0.9709901201f, -0.2187361325f, -0.09660585362f, 1.214123996f},
     {-0.09385435001f, -0.7202114479f, 0.6873767753f, 1.051462224f}},
    {{0.9721374115f, -0.06476823822f, 0.2252863255f, 1.517654996f},
     {-0.09584151699f, -0.9868916636f, 0.1298431665f, 0.5257311121f}},
    {{0.5490814303f, 0.7586196048f, 0.3507219383f, 0.6070619982f},
     {0.8285959708f, -0.4392579149f, -0.347104021f, 5.257311121f}},
    {{0.9921258754f, -0.001098710673f, -0.1252399307f, 1.214123996f},
     {-0.061220482f, -0.876612807f, -0.4772861187f, 0.0f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 0.6070619982f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.6696489291f, 0.7230710214f, 0.1695246581f, 0.6070619982f},
     {0.05601801133f, 0.1784349382f, -0.9823558191f, 4.205848897f}},
    {{0.2863114437f, 0.2070063213f, 0.9355074239f, 0.3035309991f},
     {0.6079419899f, 0.7154153424f, -0.3443652491f, 3.680117785f}}};

static constexpr uint8_t AIROCEAN_CUT_FACES[] = {
    3,  4,  4,  5,  5,  6,  6,  8,  9,  12, 12, 13, 13,
    14, 14, 15, 15, 16, 18, 18, 19, 19, 20, 21, 22, 22};
static constexpr uint8_t AIROCEAN_CUT_EDGES[] = {0, 1, 2, 1, 2, 0, 2, 0, 2,
                                                 1, 2, 1, 2, 0, 2, 0, 2, 0,
                                                 1, 2, 1, 2, 2, 1, 0, 1};

inline bool airocean_edge_is_cut(uint8_t face, uint8_t edge) {
  for (size_t index = 0; index < std::size(AIROCEAN_CUT_FACES); ++index)
    if (AIROCEAN_CUT_FACES[index] == face && AIROCEAN_CUT_EDGES[index] == edge)
      return true;
  return false;
}

inline bool airocean_same_point(const AiroceanPoint &a,
                                const AiroceanPoint &b) {
  constexpr float EPSILON = 1e-6f;
  return fabsf(a.x - b.x) <= EPSILON && fabsf(a.y - b.y) <= EPSILON;
}

inline uint8_t airocean_edge_identity(uint8_t face, uint8_t edge) {
  const AiroceanPoint &a = AIROCEAN_PLANAR_FACES[face][edge];
  const AiroceanPoint &b = AIROCEAN_PLANAR_FACES[face][(edge + 1U) % 3U];
  for (uint8_t candidate_face = 0; candidate_face < 23; ++candidate_face) {
    for (uint8_t candidate_edge = 0; candidate_edge < 3; ++candidate_edge) {
      const AiroceanPoint &c =
          AIROCEAN_PLANAR_FACES[candidate_face][candidate_edge];
      const AiroceanPoint &d =
          AIROCEAN_PLANAR_FACES[candidate_face][(candidate_edge + 1U) % 3U];
      if ((airocean_same_point(a, c) && airocean_same_point(b, d)) ||
          (airocean_same_point(a, d) && airocean_same_point(b, c)))
        return static_cast<uint8_t>(candidate_face * 3U + candidate_edge);
    }
  }
  return static_cast<uint8_t>(face * 3U + edge);
}

inline float airocean_det(const AiroceanVector &u, const AiroceanVector &v,
                          const AiroceanVector &w) {
  return u.x * (v.y * w.z - v.z * w.y) - v.x * (u.y * w.z - u.z * w.y) +
         w.x * (u.y * v.z - u.z * v.y);
}

inline bool airocean_contains(const AiroceanVector &p,
                              const AiroceanVector (&face)[3]) {
  return airocean_det(p, face[1], face[2]) <= 0.0f &&
         airocean_det(face[0], p, face[2]) <= 0.0f &&
         airocean_det(face[0], face[1], p) <= 0.0f;
}

inline float airocean_outside_score(const AiroceanVector &p,
                                    const AiroceanVector (&face)[3]) {
  return std::max(0.0f, std::max(airocean_det(p, face[1], face[2]),
                                 std::max(airocean_det(face[0], p, face[2]),
                                          airocean_det(face[0], face[1], p))));
}

inline float point_segment_distance(const AiroceanPoint &p,
                                    const AiroceanPoint &a,
                                    const AiroceanPoint &b) {
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float denom = dx * dx + dy * dy;
  const float t = denom == 0.0f
                      ? 0.0f
                      : hs::clamp(((p.x - a.x) * dx + (p.y - a.y) * dy) / denom,
                                  0.0f, 1.0f);
  const float ex = p.x - (a.x + t * dx);
  const float ey = p.y - (a.y + t * dy);
  return sqrtf(ex * ex + ey * ey);
}

inline ProjectionKernelResult
airocean_projection(const Vector &v, float central_meridian, bool horizontal) {
  const float c = cosf(central_meridian);
  const float s = sinf(central_meridian);
  const AiroceanVector p{v.x * c + v.z * s, v.z * c - v.x * s, v.y};
  uint8_t face = 0;
  for (; face < 23; ++face)
    if (airocean_contains(p, AIROCEAN_FACES[face]))
      break;
  if (face == 23) {
    float best_score = 65536.0f;
    for (uint8_t candidate = 0; candidate < 23; ++candidate) {
      const float score = airocean_outside_score(p, AIROCEAN_FACES[candidate]);
      if (score < best_score) {
        best_score = score;
        face = candidate;
      }
    }
  }
  const AiroceanVector &center = AIROCEAN_CENTERS[face];
  const AiroceanVector &normal = AIROCEAN_NORMALS[face];
  const float plane =
      center.x * normal.x + center.y * normal.y + center.z * normal.z;
  const float ray = p.x * normal.x + p.y * normal.y + p.z * normal.z;
  const float scale = plane / ray;
  const AiroceanVector q{p.x * scale, p.y * scale, p.z * scale};
  const float(&transform)[2][4] = AIROCEAN_TRANSFORMS[face];
  AiroceanPoint output{transform[0][0] * q.x + transform[0][1] * q.y +
                           transform[0][2] * q.z + transform[0][3],
                       transform[1][0] * q.x + transform[1][1] * q.y +
                           transform[1][2] * q.z + transform[1][3]};
  float edge = 65536.0f;
  for (size_t index = 0; index < std::size(AIROCEAN_CUT_FACES); ++index) {
    const uint8_t cut_face = AIROCEAN_CUT_FACES[index];
    const uint8_t cut_edge = AIROCEAN_CUT_EDGES[index];
    const AiroceanPoint &a = AIROCEAN_PLANAR_FACES[cut_face][cut_edge];
    const AiroceanPoint &b =
        cut_face == 14 && cut_edge == 0
            ? AIROCEAN_PLANAR_FACES[18][1]
            : AIROCEAN_PLANAR_FACES[cut_face][(cut_edge + 1) % 3];
    edge = std::min(edge, point_segment_distance(output, a, b));
  }
  uint8_t edge_class = 0;
  float face_edge_distance = 65536.0f;
  for (uint8_t candidate = 0; candidate < 3; ++candidate) {
    const float distance = point_segment_distance(
        output, AIROCEAN_PLANAR_FACES[face][candidate],
        AIROCEAN_PLANAR_FACES[face][(candidate + 1) % 3]);
    if (distance < face_edge_distance) {
      face_edge_distance = distance;
      edge_class = candidate;
    }
  }
  bool cut_edge = airocean_edge_is_cut(face, edge_class);
  uint8_t edge_identity = airocean_edge_identity(face, edge_class);
  if (face == 14 && edge_class == 0) {
    const AiroceanPoint &a = AIROCEAN_PLANAR_FACES[14][0];
    const AiroceanPoint &b = AIROCEAN_PLANAR_FACES[14][1];
    const float dx = b.x - a.x;
    const float dy = b.y - a.y;
    const float t =
        ((output.x - a.x) * dx + (output.y - a.y) * dy) / (dx * dx + dy * dy);
    cut_edge = t <= 0.5f;
    if (!cut_edge)
      edge_identity = airocean_edge_identity(18, 0);
  }
  if (horizontal)
    output = {5.78304223331047f - output.y, output.x};
  return {Complex(output.x, output.y),
          face,
          0,
          static_cast<uint8_t>(cut_edge),
          edge,
          0,
          projection_traits(cut_edge ? ProjectionTrait::CUT
                                     : ProjectionTrait::GLUED),
          edge_identity};
}

} // namespace shadierball
