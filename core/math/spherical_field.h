/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <utility>

#include "math/geometry.h"

namespace hs {

/**
 * @brief Allocation-free layout for a latitude-ring field on a sphere.
 * @tparam W Longitude-domain width.
 * @tparam H Number of rendered latitude rows.
 * @tparam HOffset Virtual rows below the rendered domain.
 * @details Rings are uniformly spaced in latitude. Each ring's periodic
 * longitude count follows sin(phi), producing approximately uniform physical
 * sample spacing and compact contiguous storage.
 */
template <int W, int H, int HOffset = H_OFFSET> class SphericalFieldLayout {
public:
  struct Ring {
    int y;
    int samples;
    int offset;
  };

  struct Row {
    Ring lower;
    Ring upper;
    float mix;
  };

  struct Longitude {
    int left;
    int right;
    float mix;
  };

  struct Coordinates {
    float x;
    float y;
  };

  constexpr explicit SphericalFieldLayout(int spacing, int north_infill = 0,
                                          int south_infill = 0,
                                          int equator_samples = 0)
      : spacing(spacing), north_infill(north_infill),
        south_infill(south_infill), equator_samples(equator_samples) {}

  constexpr int latitude_spacing() const { return spacing; }

  constexpr int ring_count() const {
    if (spacing <= 0)
      return 0;
    int count = 1;
    for (int y = 0; y < H - 1; ++count)
      y = next_ring_y(y);
    return count;
  }

  constexpr int sample_count() const {
    if (spacing <= 0)
      return 0;
    int count = 0;
    for (int ring_index = 0; ring_index < ring_count(); ++ring_index)
      count += samples_on_ring(ring_y(ring_index));
    return count;
  }

  constexpr Ring ring(int ring_index) const {
    int offset = 0;
    for (int i = 0; i < ring_index; ++i)
      offset += samples_on_ring(ring_y(i));
    const int y = ring_y(ring_index);
    return {y, samples_on_ring(y), offset};
  }

  constexpr Coordinates sample_coordinates(const Ring &ring,
                                           int sample_index) const {
    return {static_cast<float>(sample_index * W) / ring.samples,
            static_cast<float>(ring.y)};
  }

  Vector sample_vector(const Ring &ring, int sample_index) const {
    const Coordinates point = sample_coordinates(ring, sample_index);
    const float theta = point.x * (2.0f * PI_F) / W;
    const float phi = point.y * PI_F / (H + HOffset - 1);
    return Vector(Spherical(theta, phi));
  }

  Coordinates project(const Vector &value) const {
    const Spherical spherical(value);
    return {(spherical.theta * W) / (2.0f * PI_F),
            (spherical.phi * (H + HOffset - 1)) / PI_F};
  }

  constexpr Row row(float y) const {
    const float bounded_y = std::clamp(y, 0.0f, static_cast<float>(H - 1));
    const int lower_index = ring_index_at_or_before(bounded_y);
    int upper_index = std::min(lower_index + 1, ring_count() - 1);
    Ring lower = ring(lower_index);
    Ring upper = ring(upper_index);
    const int height = upper.y - lower.y;
    const float mix =
        height > 0 ? (bounded_y - lower.y) / static_cast<float>(height) : 0.0f;
    return {lower, upper, mix};
  }

  constexpr int ring_index_at_or_before(float y) const {
    const float bounded_y = std::clamp(y, 0.0f, static_cast<float>(H - 1));
    int index = 0;
    while (index + 1 < ring_count() && ring_y(index + 1) <= bounded_y)
      ++index;
    return index;
  }

  constexpr int ring_index_at_or_after(float y) const {
    const int lower = ring_index_at_or_before(y);
    return ring(lower).y < y ? std::min(lower + 1, ring_count() - 1) : lower;
  }

  constexpr Longitude longitude(const Ring &ring, float x) const {
    float wrapped_x = std::fmod(x, static_cast<float>(W));
    if (wrapped_x < 0.0f)
      wrapped_x += W;
    const float position = wrapped_x * ring.samples / W;
    const int left = static_cast<int>(position);
    const int right = left + 1 < ring.samples ? left + 1 : 0;
    return {ring.offset + left, ring.offset + right, position - left};
  }

  /**
   * @brief Locates an integer longitude known to be inside the domain.
   * @param ring Target latitude ring.
   * @param x Longitude coordinate in [0, W).
   */
  constexpr Longitude longitude_bounded(const Ring &ring, int x) const {
    assert(x >= 0 && x < W);
    const int position = x * ring.samples;
    const int left = position / W;
    const int phase = position - left * W;
    const int right = left + 1 < ring.samples ? left + 1 : 0;
    return {ring.offset + left, ring.offset + right,
            static_cast<float>(phase) / W};
  }

  template <typename Value, typename Lerp>
  constexpr Value interpolate(const Value *values, float x, float y,
                              Lerp &&lerp) const {
    const Row latitude = row(y);
    const Longitude lower = longitude(latitude.lower, x);
    const Longitude upper = longitude(latitude.upper, x);
    const Value a = lerp(values[lower.left], values[lower.right], lower.mix);
    const Value b = lerp(values[upper.left], values[upper.right], upper.mix);
    return lerp(a, b, latitude.mix);
  }

private:
  static_assert(W > 0 && H > 1 && H + HOffset > 1);

  constexpr int ring_y(int ring_index) const {
    int y = 0;
    for (int i = 0; i < ring_index; ++i)
      y = next_ring_y(y);
    return y;
  }

  constexpr int maximum_longitude_samples() const {
    return equator_samples > 0
               ? equator_samples
               : (2 * (H + HOffset - 1) + spacing / 2) / spacing;
  }

  static constexpr float latitude_sine(int y) {
    float phi = (static_cast<float>(y) * PI_F) / (H + HOffset - 1);
    if (phi > PI_F * 0.5f)
      phi = PI_F - phi;
    const float phi2 = phi * phi;
    return phi *
           (1.0f + phi2 * (-1.0f / 6.0f +
                           phi2 * (1.0f / 120.0f + phi2 * (-1.0f / 5040.0f +
                                                           phi2 / 362880.0f))));
  }

  constexpr int samples_on_ring(int y) const {
    if (y < north_infill || y >= H - south_infill)
      return maximum_longitude_samples();
    const int count =
        static_cast<int>(maximum_longitude_samples() * latitude_sine(y) + 0.5f);
    return std::max(count, 1);
  }

  constexpr int next_ring_y(int y) const {
    if (y < north_infill - 1)
      return y + 1;
    const int south_begin = H - south_infill;
    if (y >= south_begin)
      return std::min(y + 1, H - 1);
    const int next_regular = ((y / spacing) + 1) * spacing;
    if (next_regular >= south_begin)
      return south_infill > 0 ? south_begin : H - 1;
    return std::min(next_regular, H - 1);
  }

  int spacing;
  int north_infill;
  int south_infill;
  int equator_samples;
};

/**
 * @brief Non-owning value field stored on a SphericalFieldLayout.
 * @tparam Value Stored sample type.
 */
template <typename Value, int W, int H, int HOffset = H_OFFSET>
class SphericalField {
public:
  using Layout = SphericalFieldLayout<W, H, HOffset>;
  using Ring = typename Layout::Ring;

  constexpr SphericalField(Value *values, const Layout &layout)
      : layout(layout), values(values) {}

  constexpr int sample_count() const { return layout.sample_count(); }

  template <typename Populate>
  void populate(int ring_begin, int ring_end, Populate &&populate_sample) {
    for (int ring_index = ring_begin; ring_index <= ring_end; ++ring_index) {
      const Ring ring = layout.ring(ring_index);
      const Vector meridian = pixel_to_vector<W, H>(0, ring.y);
      const float theta_step = 2.0f * PI_F / ring.samples;
      const float step_cos = cosf(theta_step);
      const float step_sin = sinf(theta_step);
      float theta_cos = 1.0f;
      float theta_sin = 0.0f;
      for (int sample = 0; sample < ring.samples; ++sample) {
        const Vector position(meridian.x * theta_cos, meridian.y,
                              meridian.x * theta_sin);
        values[ring.offset + sample] =
            populate_sample(position, layout.sample_coordinates(ring, sample));
        const float next_cos = theta_cos * step_cos - theta_sin * step_sin;
        theta_sin = theta_sin * step_cos + theta_cos * step_sin;
        theta_cos = next_cos;
      }
    }
  }

  template <typename Lerp>
  Value interpolate(float x, float y, Lerp &&lerp) const {
    return layout.interpolate(values, x, y, std::forward<Lerp>(lerp));
  }

  template <typename Lerp>
  Value interpolate(const Ring &ring, int x, Lerp &&lerp) const {
    const auto longitude = layout.longitude_bounded(ring, x);
    return lerp(values[longitude.left], values[longitude.right], longitude.mix);
  }

  Value *data() { return values; }
  const Value *data() const { return values; }

  const Layout layout;

private:
  Value *values;
};

} // namespace hs
