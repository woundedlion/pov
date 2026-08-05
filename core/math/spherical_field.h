/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file spherical_field.h
 * @brief SphericalFieldLayout and SphericalField: allocation-free fields of
 *        latitude rings over the sphere.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

#include "engine/util.h"
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
  /**
   * @brief One latitude ring.
   * @details y is its row, samples its periodic longitude count, and offset the
   * index of its first sample in the field's contiguous storage.
   */
  struct Ring {
    int y;
    int samples;
    int offset;
  };

  /**
   * @brief A fractional latitude bracketed by two rings.
   * @details mix is the weight of upper, in [0, 1].
   */
  struct Row {
    Ring lower;
    Ring upper;
    float mix;
  };

  /**
   * @brief A fractional longitude bracketed by two absolute sample indices.
   * @details mix is the weight of right, in [0, 1]; it reaches exactly 1 at the
   *   seam, where left is clamped to the ring's last sample.
   */
  struct Longitude {
    int left;
    int right;
    float mix;
  };

  /** @brief Field coordinates: longitude x in [0, W), latitude row y. */
  struct Coordinates {
    float x;
    float y;
  };

  /**
   * @brief Builds the ring chain over the rendered domain.
   * @param spacing Latitude rows between consecutive rings; must be > 0.
   * @param north_infill Leading rows given one ring each at full longitude
   *   resolution.
   * @param south_infill Trailing rows given one ring each at full longitude
   *   resolution.
   * @param equator_samples Longitude samples on the widest ring; 0 derives the
   *   count from spacing.
   */
  constexpr explicit SphericalFieldLayout(int spacing, int north_infill = 0,
                                          int south_infill = 0,
                                          int equator_samples = 0)
      : spacing(spacing), north_infill(north_infill),
        south_infill(south_infill), equator_samples(equator_samples) {
    HS_CHECK(spacing > 0, "SphericalFieldLayout: spacing must be > 0");
  }

  /** @brief Rows that alias a whole longitude range onto one physical point:
   *  row 0 always, plus row H-1 when no virtual rows separate it from the
   *  south pole. */
  static constexpr int POLE_COUNT = HOffset == 0 ? 2 : 1;

  /** @brief Rings in the chain, counting both endpoint rows. */
  constexpr int ring_count() const {
    int count = 1;
    for (int y = 0; y < H - 1; ++count)
      y = next_ring_y(y);
    return count;
  }

  /** @brief Samples across every ring, i.e. the storage a field needs. */
  constexpr int sample_count() const {
    int count = 0;
    for (int y = 0;; y = next_ring_y(y)) {
      count += samples_on_ring(y);
      if (y >= H - 1)
        return count;
    }
  }

  /**
   * @brief Walks the ring chain to one ring.
   * @param ring_index Ring position in [0, ring_count()); O(ring_index).
   * @details Prefer next_ring() when iterating: an index loop over the chain is
   * quadratic.
   */
  constexpr Ring ring(int ring_index) const {
    HS_CHECK(ring_index >= 0, "SphericalFieldLayout: negative ring index %d",
             ring_index);
    int offset = 0;
    int y = 0;
    for (int i = 0; i < ring_index; ++i) {
      HS_CHECK(y < H - 1, "SphericalFieldLayout: ring index %d out of range",
               ring_index);
      offset += samples_on_ring(y);
      y = next_ring_y(y);
    }
    return {y, samples_on_ring(y), offset};
  }

  /**
   * @brief Advances one ring down the chain in constant time.
   * @param ring Current ring.
   * @return The next ring; the last ring is its own successor.
   */
  constexpr Ring next_ring(const Ring &ring) const {
    if (ring.y >= H - 1)
      return ring;
    const int y = next_ring_y(ring.y);
    return {y, samples_on_ring(y), ring.offset + ring.samples};
  }

  /**
   * @brief Field coordinates of one sample on a ring.
   * @param ring Target latitude ring.
   * @param sample_index Sample position in [0, ring.samples).
   * @return Coordinates with x in [0, W) and y at the ring's row.
   */
  constexpr Coordinates sample_coordinates(const Ring &ring,
                                           int sample_index) const {
    return {static_cast<float>(sample_index * W) / ring.samples,
            static_cast<float>(ring.y)};
  }

  /**
   * @brief Reconstructs the unit vector for one sample on a ring.
   * @param ring Target latitude ring.
   * @param sample_index Sample position in [0, ring.samples).
   * @return Unit vector for sample_coordinates(), whose x lies in [0, W).
   */
  Vector sample_vector(const Ring &ring, int sample_index) const {
    const Coordinates point = sample_coordinates(ring, sample_index);
    const float theta = point.x * (2.0f * PI_F) / W;
    const float phi = y_to_phi_virtual(point.y, H + HOffset);
    return Vector(Spherical(theta, phi));
  }

  /**
   * @brief Maps a unit vector back to fractional field coordinates.
   * @param value Unit vector on the sphere.
   * @return Coordinates with x in [-W/2, W/2] and y in [0, H + HOffset - 1].
   *   x is signed because Spherical::theta is atan2's [-pi, pi], unlike
   *   sample_coordinates()' [0, W); longitude() accepts either convention.
   */
  Coordinates project(const Vector &value) const {
    const Spherical spherical(value);
    return {(spherical.theta * W) / (2.0f * PI_F),
            phi_to_y_virtual(spherical.phi, H + HOffset)};
  }

  /**
   * @brief Returns an odd longitude footprint with equatorial pixel width.
   * @param y Latitude row in the rendered domain.
   * @return Odd width; pole rows, where every longitude collapses onto one
   *   point, saturate to the widest odd footprint the row admits.
   */
  int longitude_filter_width(int y) const {
    const int maximum_odd_width = (W & 1) ? W : W - 1;
    const float sine = latitude_sine(y);
    if (sine < 1.0f / W)
      return maximum_odd_width;
    int width = static_cast<int>(ceilf(1.0f / sine)) | 1;
    return hs::clamp(width, 1, maximum_odd_width);
  }

  /**
   * @brief Reconstructs one equirectangular row at its spherical width.
   * @tparam Accumulator Default-constructible rolling accumulator with
   *   add(), remove(), and average().
   * @param source Dense row containing W values.
   * @param y Latitude row controlling the longitude footprint.
   * @param emit Receives each destination column and filtered value.
   */
  template <typename Accumulator, typename Value, typename Emit>
  void reconstruct_longitude_row(const Value *source, int y,
                                 Emit &&emit) const {
    const int width = longitude_filter_width(y);
    const int radius = width / 2;
    Accumulator accumulator;
    for (int dx = -radius; dx <= radius; ++dx) {
      int x = dx;
      if (x < 0)
        x += W;
      accumulator.add(source[x]);
    }

    for (int x = 0; x < W; ++x) {
      emit(x, accumulator.average(width));
      int remove = x - radius;
      if (remove < 0)
        remove += W;
      int add = x + radius + 1;
      if (add >= W)
        add -= W;
      accumulator.remove(source[remove]);
      accumulator.add(source[add]);
    }
  }

  /** @brief Reflects one lattice coordinate across spherical seams and poles. */
  static bool wrap_sample(int &x, int &y) {
    return ::pole_wrap<W, H, HOffset>(x, y);
  }

  /**
   * @brief Bilinearly samples a dense equirectangular field across seams and
   * poles.
   * @param x Fractional longitude in [-W, 2W).
   * @param y Fractional latitude row.
   * @param poles POLE_COUNT shared values for the longitude-aliased pole rows:
   *   [0] for row 0, [1] for row H-1 where that row is the south pole.
   * @param outside Value returned where the rendered domain has no sample.
   * @param load Loads an in-domain, non-pole lattice sample.
   * @param combine Combines four topology-correct taps and fractional
   *   coordinates into the result.
   */
  template <typename Value, typename Load, typename Combine>
  __attribute__((always_inline)) decltype(auto)
  sample_bilinear(float x, float y, const Value *poles, const Value &outside,
                  Load &&load, Combine &&combine) const {
    const float floor_x = std::floor(x);
    const float floor_y = std::floor(y);
    int x0 = static_cast<int>(floor_x);
    const int y0 = static_cast<int>(floor_y);
    const float fx = x - floor_x;
    const float fy = y - floor_y;

    x0 = ::fast_wrap(x0, W);
    const int x1 = x0 + 1 < W ? x0 + 1 : 0;

    // Row H-1 is the south pole only when HOffset is 0; it then leaves the
    // direct-load band so the pole substitution reaches it.
    constexpr int LAST_DIRECT_ROW = POLE_COUNT == 2 ? H - 3 : H - 2;
    if (y0 > 0 && y0 <= LAST_DIRECT_ROW) {
      return std::forward<Combine>(combine)(load(x0, y0), load(x1, y0),
                                            load(x0, y0 + 1), load(x1, y0 + 1),
                                            fx, fy);
    }

    auto tap = [&](int sample_x, int sample_y) {
      if (!wrap_sample(sample_x, sample_y))
        return outside;
      if (sample_y == 0)
        return poles[0];
      if constexpr (POLE_COUNT == 2) {
        if (sample_y == H - 1)
          return poles[1];
      }
      return load(sample_x, sample_y);
    };
    return std::forward<Combine>(combine)(
        tap(x0, y0), tap(x1, y0), tap(x0, y0 + 1), tap(x1, y0 + 1), fx, fy);
  }

  /**
   * @brief Bilinearly samples a row-major three-channel field.
   * @param source Dense W-by-H source field.
   * @param poles POLE_COUNT shared values for the longitude-aliased pole rows:
   *   [0] for row 0, [1] for row H-1 where that row is the south pole.
   * @param x Fractional longitude in [-W, 2W).
   * @param y Fractional latitude row.
   * @param r Out: interpolated red channel.
   * @param g Out: interpolated green channel.
   * @param b Out: interpolated blue channel.
   */
  template <typename Pixel>
  __attribute__((always_inline)) void
  sample_bilinear_rgb(const Pixel *source, const Pixel *poles, float x, float y,
                      float &r, float &g, float &b) const {
    static constexpr Pixel OUTSIDE{};
    sample_bilinear(
        x, y, poles, OUTSIDE,
        [source](int sample_x, int sample_y) {
          return source[sample_y * W + sample_x];
        },
        [&](const Pixel &p00, const Pixel &p10, const Pixel &p01,
            const Pixel &p11, float fx, float fy) {
          const float w00 = (1.0f - fx) * (1.0f - fy);
          const float w10 = fx * (1.0f - fy);
          const float w01 = (1.0f - fx) * fy;
          const float w11 = fx * fy;
          r = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
          g = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
          b = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
        });
  }

  /**
   * @brief Brackets a fractional latitude between two rings.
   * @param y Latitude row, clamped to [0, H-1].
   * @details Walks the ring chain; O(ring index).
   */
  constexpr Row row(float y) const {
    const float bounded_y = hs::clamp(y, 0.0f, static_cast<float>(H - 1));
    const Ring lower = ring(ring_index_at_or_before(bounded_y));
    const Ring upper = next_ring(lower);
    const int height = upper.y - lower.y;
    const float mix =
        height > 0 ? (bounded_y - lower.y) / static_cast<float>(height) : 0.0f;
    return {lower, upper, mix};
  }

  /**
   * @brief Index of the last ring whose row is at or above y.
   * @param y Latitude row, clamped to [0, H-1].
   */
  constexpr int ring_index_at_or_before(float y) const {
    const float bounded_y = hs::clamp(y, 0.0f, static_cast<float>(H - 1));
    int index = 0;
    int current_y = 0;
    while (current_y < H - 1 && next_ring_y(current_y) <= bounded_y) {
      current_y = next_ring_y(current_y);
      ++index;
    }
    return index;
  }

  /**
   * @brief Index of the first ring whose row is at or below y, saturating at
   * the last ring.
   * @param y Latitude row, clamped to [0, H-1].
   */
  constexpr int ring_index_at_or_after(float y) const {
    const int lower = ring_index_at_or_before(y);
    return ring(lower).y < y ? std::min(lower + 1, ring_count() - 1) : lower;
  }

  /**
   * @brief Locates the samples bracketing a fractional longitude.
   * @param ring Target latitude ring.
   * @param x Longitude in any range; wrapped into [0, W), so both project()'s
   *   signed x and sample_coordinates()' [0, W) are accepted. A non-finite x
   *   saturates to the ring's seam.
   * @return Absolute sample indices (ring.offset applied) and their mix.
   */
  constexpr Longitude longitude(const Ring &ring, float x) const {
    // fmod turns either infinity into NaN, which wrap() passes through; the
    // clamp saturates it before the cast, which would otherwise be UB.
    const float wrapped_x = wrap(x, static_cast<float>(W));
    const float position = hs::clamp(wrapped_x * ring.samples / W, 0.0f,
                                     static_cast<float>(ring.samples));
    // The scale rounds up to exactly ring.samples near the seam, so the
    // truncation would otherwise land one sample past the ring.
    const int left = std::min(static_cast<int>(position), ring.samples - 1);
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

private:
  static_assert(W > 0 && H > 1 && H + HOffset > 1);

  constexpr int maximum_longitude_samples() const {
    return equator_samples > 0
               ? equator_samples
               : (2 * (H + HOffset - 1) + spacing / 2) / spacing;
  }

  /** @brief sin(phi) at row y, as a Taylor series because sinf is not
   *  constexpr. */
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
    Ring ring = layout.ring(ring_begin);
    for (int ring_index = ring_begin; ring_index <= ring_end;
         ++ring_index, ring = layout.next_ring(ring)) {
      const Vector meridian = layout.sample_vector(ring, 0);
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

  Value *data() { return values; }
  const Value *data() const { return values; }

  const Layout layout;

private:
  Value *values;
};

} // namespace hs
