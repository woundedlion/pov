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
#include <cstdint>
#include <type_traits>
#include <utility>

#include "math/geometry.h"

namespace hs {

/**
 * @brief Allocation-free layout for a latitude-ring field on a sphere.
 * @tparam W Longitude-domain width.
 * @tparam H Number of rendered latitude rows.
 * @tparam HOffset South offset; -1 selects the active display profile.
 * @details Rings follow the requested latitude-row spacing, with endpoint and
 *          optional infill rings added. Each ring's periodic
 * longitude count follows sin(phi), producing approximately uniform physical
 * sample spacing and compact contiguous storage.
 */
template <int W, int H, int HOffset = -1> class SphericalFieldLayout {
public:
  using Geometry = math::DisplayGeometry<H, HOffset>;
#if HS_RUNTIME_DISPLAY_GEOMETRY
  inline static const bool &HAS_NORTH_POLE = Geometry::HAS_NORTH_POLE;
  inline static const bool &HAS_SOUTH_POLE = Geometry::HAS_SOUTH_POLE;
#else
  static constexpr bool HAS_NORTH_POLE = Geometry::HAS_NORTH_POLE;
  static constexpr bool HAS_SOUTH_POLE = Geometry::HAS_SOUTH_POLE;
#endif

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

  /** @brief Field coordinates: longitude x in columns (possibly signed), latitude row y. */
  struct Coordinates {
    float x;
    float y;
  };

  /**
   * @brief Builds the ring chain over the rendered domain.
   * @param spacing Latitude rows between consecutive rings; must be > 0.
   * @param north_infill Leading rows given one ring each at full longitude
   *   resolution; the two infill bands must be non-negative and must not
   *   overlap.
   * @param south_infill Trailing rows given one ring each at full longitude
   *   resolution.
   * @param equator_samples Longitude samples on the widest ring; must be
   *   non-negative, and 0 derives the count from spacing.
   */
  constexpr explicit SphericalFieldLayout(int spacing, int north_infill = 0,
                                          int south_infill = 0,
                                          int equator_samples = 0)
      : spacing(spacing), north_infill(north_infill),
        south_infill(south_infill), equator_samples(equator_samples) {
    HS_CHECK(spacing > 0, "SphericalFieldLayout: spacing must be > 0");
    HS_CHECK(north_infill >= 0 && south_infill >= 0 &&
                 north_infill + south_infill <= H,
             "SphericalFieldLayout: infills %d + %d must be non-negative and "
             "fit within H = %d",
             north_infill, south_infill, H);
    HS_CHECK(equator_samples >= 0,
             "SphericalFieldLayout: equator_samples %d must be >= 0",
             equator_samples);
  }

  /** @brief Pole-value storage capacity; runtime geometry reserves both poles. */
#if HS_RUNTIME_DISPLAY_GEOMETRY
  static constexpr int POLE_STORAGE_COUNT = 2;
#else
  static constexpr int POLE_COUNT = int(HAS_NORTH_POLE) + int(HAS_SOUTH_POLE);
  static constexpr int POLE_STORAGE_COUNT = POLE_COUNT > 0 ? POLE_COUNT : 1;
#endif

  /** @brief Bound a sampler's fractional row must stay under in absolute
   *  value: NaN, an infinity, or a row this far out makes the truncation to
   *  the lattice index UB. */
  static constexpr float ROW_LIMIT = 1e9f;

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
  math::Vector sample_vector(const Ring &ring, int sample_index) const {
    const Coordinates point = sample_coordinates(ring, sample_index);
    const float theta = point.x * (2.0f * math::PI_F) / W;
    const float phi = Geometry::row_to_phi(point.y);
    return math::Vector(math::Spherical(theta, phi));
  }

  /**
   * @brief Maps a unit vector back to fractional field coordinates.
   * @note Derives `theta`/`phi` with the approximate `fast_atan2`/`fast_acos`,
   *   so the coordinates are sub-sample inexact and vector → coordinates →
   *   vector does not bit-exactly invert the exact-trig `sample_vector()`.
   * @param value Unit vector on the sphere.
   * @return Coordinates with x in [-W/2, W/2]; missing caps project outside [0,H-1].
   */
  Coordinates project(const math::Vector &value) const {
    const math::Spherical spherical(value);
    return {(spherical.theta * W) / (2.0f * math::PI_F),
            Geometry::phi_to_row(spherical.phi)};
  }

  /**
   * @brief sin(phi) at row y.
   * @param y Latitude row in the rendered domain.
   * @return The row's latitude sine (constexpr Taylor series); the longitude
   *   density of the row relative to the equator.
   */
  static constexpr float latitude_sine(int y) {
    float phi = Geometry::row_to_phi(static_cast<float>(y));
    if (phi > math::PI_F * 0.5f)
      phi = math::PI_F - phi;
    const float phi2 = phi * phi;
    return phi *
           (1.0f + phi2 * (-1.0f / 6.0f +
                           phi2 * (1.0f / 120.0f + phi2 * (-1.0f / 5040.0f +
                                                           phi2 / 362880.0f))));
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
   * @note A footprint that spans the whole row emits the row mean, constant in
   *   longitude.
   */
  template <typename Accumulator, typename Value, typename Emit>
  void reconstruct_longitude_row(const Value *source, int y,
                                 Emit &&emit) const {
    const int width = longitude_filter_width(y);
    Accumulator accumulator;
    if (width >= W - 1) {
      for (int x = 0; x < W; ++x)
        accumulator.add(source[x]);
      const auto mean = accumulator.average(W);
      for (int x = 0; x < W; ++x)
        emit(x, mean);
      return;
    }

    const int radius = width / 2;
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

  /** @brief Reflects one lattice coordinate across poles.
   * @pre The column is already in [0, W).
   */
  static bool wrap_sample(int &x, int &y) {
    return ::math::pole_wrap<W, H, HOffset>(x, y);
  }

  /**
   * @brief Bilinearly samples a dense equirectangular field across seams and
   * poles.
   * @param x Fractional longitude in [-W, 2W).
   * @param y Fractional latitude row.
   * @param poles Shared values for true pole rows, north then south when present.
   *   The south pole uses slot [0] when no north pole is present.
   * @param outside Value returned where the rendered domain has no sample.
   * @param load Loads an in-domain, non-pole lattice sample.
   * @param combine Combines four topology-correct taps and fractional
   *   coordinates into the result.
   */
  template <typename Value, typename Load, typename Combine>
  __attribute__((always_inline)) decltype(auto)
  sample_bilinear(float x, float y, const Value (&poles)[POLE_STORAGE_COUNT],
                  const Value &outside, Load &&load, Combine &&combine) const {
    const Footprint tap = bilinear_footprint(x, y);
    const int x0 = tap.x0;
    const int x1 = tap.x1;
    const int y0 = tap.y0;

    if (in_direct_band(y0)) {
      return std::forward<Combine>(combine)(load(x0, y0), load(x1, y0),
                                            load(x0, y0 + 1), load(x1, y0 + 1),
                                            tap.fx, tap.fy);
    }
    return std::forward<Combine>(combine)(
        pole_tap(x0, y0, poles, outside, load),
        pole_tap(x1, y0, poles, outside, load),
        pole_tap(x0, y0 + 1, poles, outside, load),
        pole_tap(x1, y0 + 1, poles, outside, load), tap.fx, tap.fy);
  }

  /**
   * @brief Bilinearly samples a row-major three-channel field.
   * @param source Dense W-by-H source field.
   * @param poles Shared values for true pole rows, north then south when present.
   *   The south pole uses slot [0] when no north pole is present.
   * @param x Fractional longitude in [-W, 2W).
   * @param y Fractional latitude row.
   * @param r Out: interpolated red channel.
   * @param g Out: interpolated green channel.
   * @param b Out: interpolated blue channel.
   * @note Out-of-domain taps are black; use sample_bilinear() to supply a
   *   different outside value.
   * @details Inside the direct band, 16-bit unsigned channels blend through
   *   combine_rgb_q15; other channel types, and seam or pole rows, blend in
   *   float.
   */
  template <typename Pixel>
  __attribute__((always_inline)) void
  sample_bilinear_rgb(const Pixel *source,
                      const Pixel (&poles)[POLE_STORAGE_COUNT], float x,
                      float y, float &r, float &g, float &b) const {
    const Footprint tap = bilinear_footprint(x, y);

    if (!in_direct_band(tap.y0)) {
      sample_bilinear_rgb_poles(source, poles, tap.x0, tap.x1, tap.y0, tap.fx,
                                tap.fy, r, g, b);
      return;
    }
    const int row = tap.y0 * W;
    const int next_row = row + W;
    if constexpr (sizeof(source->r) == 2 &&
                  std::is_unsigned_v<decltype(source->r)>)
      combine_rgb_q15(source[row + tap.x0], source[row + tap.x1],
                      source[next_row + tap.x0], source[next_row + tap.x1],
                      tap.fx, tap.fy, r, g, b);
    else
      combine_rgb(source[row + tap.x0], source[row + tap.x1],
                  source[next_row + tap.x0], source[next_row + tap.x1], tap.fx,
                  tap.fy, r, g, b);
  }

  /**
   * @brief Brackets a fractional latitude between two rings.
   * @param y Latitude row, clamped to [0, H-1].
   * @details Walks the ring chain; O(ring index).
   */
  constexpr Row row(float y) const {
    const float bounded_y = hs::clamp(y, 0.0f, static_cast<float>(H - 1));
    const Ring lower = ring_at_or_before(bounded_y).ring;
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
    return ring_at_or_before(hs::clamp(y, 0.0f, static_cast<float>(H - 1)))
        .index;
  }

  /**
   * @brief Index of the first ring whose row is at or below y, saturating at
   * the last ring.
   * @param y Latitude row, clamped to [0, H-1].
   */
  constexpr int ring_index_at_or_after(float y) const {
    const IndexedRing lower =
        ring_at_or_before(hs::clamp(y, 0.0f, static_cast<float>(H - 1)));
    // The last ring, the only one with row H-1, is its own successor.
    return lower.ring.y < y && lower.ring.y < H - 1 ? lower.index + 1
                                                    : lower.index;
  }

  /**
   * @brief Locates the samples bracketing a fractional longitude.
   * @param ring Target latitude ring.
   * @param x Longitude in any range; wrapped into [0, W). A non-finite x
   *   saturates to the ring's seam.
   * @return Absolute sample indices (ring.offset applied) and their mix.
   */
  constexpr Longitude longitude(const Ring &ring, float x) const {
    // fmod turns either infinity into NaN, which wrap() passes through; the
    // clamp saturates it before the cast, which would otherwise be UB.
    const float wrapped_x = math::wrap(x, static_cast<float>(W));
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
  static_assert(W > 0 && H > 1);

  /** @brief A ring paired with its position in the chain. */
  struct IndexedRing {
    Ring ring;
    int index;
  };

  /**
   * @brief Walks the chain once to the last ring whose row is at or above y.
   * @param bounded_y Latitude row already clamped to [0, H-1].
   */
  constexpr IndexedRing ring_at_or_before(float bounded_y) const {
    int index = 0;
    int offset = 0;
    int y = 0;
    while (y < H - 1) {
      const int next = next_ring_y(y);
      if (next > bounded_y)
        break;
      offset += samples_on_ring(y);
      y = next;
      ++index;
    }
    return {{y, samples_on_ring(y), offset}, index};
  }

  /** @brief True when row y0's bilinear footprint needs no seam or pole
   *  substitution. */
  static constexpr bool in_direct_band(int y0) {
    return y0 >= (HAS_NORTH_POLE ? 1 : 0) &&
           y0 <= (HAS_SOUTH_POLE ? H - 3 : H - 2);
  }

  /**
   * @brief The two lattice columns, upper row, and fractional weights of a
   * bilinear footprint.
   */
  struct Footprint {
    int x0;
    int x1;
    int y0;
    float fx;
    float fy;
  };

  /**
   * @brief Resolves a fractional coordinate to its bilinear footprint.
   * @param x Fractional longitude in [-W, 2W).
   * @param y Fractional latitude row; rows outside [0, H) are allowed but must
   *   stay under ROW_LIMIT in absolute value.
   * @details Wraps the columns across the seam; the row is left for the
   * caller's pole policy.
   */
  __attribute__((always_inline)) static Footprint bilinear_footprint(float x,
                                                                     float y) {
    assert(x >= -static_cast<float>(W) && x < 2.0f * static_cast<float>(W));
    assert(std::fabs(y) < ROW_LIMIT);
    if (Geometry::OFFSET < 0 && !(HAS_NORTH_POLE && HAS_SOUTH_POLE)) {
      if (y < Geometry::NORTH_POLE_ROW || y > Geometry::SOUTH_POLE_ROW)
        math::pole_wrap<W, H, HOffset>(x, y);
    }
    const float floor_x = std::floor(x);
    const float floor_y = std::floor(y);
    const int x0 = ::math::fast_wrap(static_cast<int>(floor_x), W);
    return {x0, x0 + 1 < W ? x0 + 1 : 0, static_cast<int>(floor_y), x - floor_x,
            y - floor_y};
  }

  /**
   * @brief One topology-correct tap for a footprint outside the direct band.
   * @param sample_x Lattice column, before seam reflection.
   * @param sample_y Lattice row, before pole reflection.
   * @param poles Shared pole-row values.
   * @param outside Value for a coordinate the rendered domain has no sample
   *   for.
   * @param load Loads an in-domain, non-pole lattice sample.
   */
  template <typename Value, typename Load>
  __attribute__((always_inline)) Value
  pole_tap(int sample_x, int sample_y, const Value (&poles)[POLE_STORAGE_COUNT],
           const Value &outside, Load &&load) const {
    if (Geometry::OFFSET < 0 && !(HAS_NORTH_POLE && HAS_SOUTH_POLE)) {
      if (sample_y < 0 || sample_y >= H)
        return outside;
    } else if (!wrap_sample(sample_x, sample_y)) {
      return outside;
    }
    if (HAS_NORTH_POLE) {
      if (sample_y == 0)
        return poles[0];
    }
    if (HAS_SOUTH_POLE) {
      if (sample_y == H - 1) {
        if constexpr (POLE_STORAGE_COUNT > 1)
          return poles[HAS_NORTH_POLE ? 1 : 0];
        else
          return poles[0];
      }
    }
    return load(sample_x, sample_y);
  }

  /**
   * @brief Bilinearly blends four 16-bit taps with Q15 integer weights.
   * @details The truncated Q15 weights sum to one exactly and each channel
   * accumulates in 31 bits. The result differs from combine_rgb only by the
   * weight truncation.
   */
  template <typename Pixel>
  __attribute__((always_inline)) static void
  combine_rgb_q15(const Pixel &p00, const Pixel &p10, const Pixel &p01,
                  const Pixel &p11, float fx, float fy, float &r, float &g,
                  float &b) {
    static_assert(sizeof(p00.r) == 2 && std::is_unsigned_v<decltype(p00.r)>,
                  "combine_rgb_q15 blends 16-bit unsigned channels");
    static_assert(sizeof(p00.g) == 2 && sizeof(p00.b) == 2);
    constexpr uint32_t ONE = 1u << 15;
    const uint32_t wx = static_cast<uint32_t>(fx * static_cast<float>(ONE));
    const uint32_t wy = static_cast<uint32_t>(fy * static_cast<float>(ONE));
    const uint32_t w11 = (wx * wy) >> 15;
    const uint32_t w10 = wx - w11;
    const uint32_t w01 = wy - w11;
    const uint32_t w00 = ONE - wx - wy + w11;
    const uint32_t rs = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
    const uint32_t gs = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
    const uint32_t bs = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
    constexpr float INVERSE_ONE = 1.0f / static_cast<float>(ONE);
    r = static_cast<float>(rs) * INVERSE_ONE;
    g = static_cast<float>(gs) * INVERSE_ONE;
    b = static_cast<float>(bs) * INVERSE_ONE;
  }

  /** @brief Bilinearly blends four taps into unclamped float channels. */
  template <typename Pixel>
  __attribute__((always_inline)) static void
  combine_rgb(const Pixel &p00, const Pixel &p10, const Pixel &p01,
              const Pixel &p11, float fx, float fy, float &r, float &g,
              float &b) {
    const float w00 = (1.0f - fx) * (1.0f - fy);
    const float w10 = fx * (1.0f - fy);
    const float w01 = (1.0f - fx) * fy;
    const float w11 = fx * fy;
    r = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
    g = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
    b = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
  }

  /**
   * @brief sample_bilinear_rgb()'s footprint outside the direct band.
   * @details Handles footprints that touch pole rows or leave the rendered
   * domain.
   */
  template <typename Pixel>
  HS_NOINLINE_NOCLONE void sample_bilinear_rgb_poles(
      const Pixel *source, const Pixel (&poles)[POLE_STORAGE_COUNT], int x0,
      int x1, int y0, float fx, float fy, float &r, float &g, float &b) const {
    static constexpr Pixel OUTSIDE{};
    auto load = [source](int sample_x, int sample_y) {
      return source[sample_y * W + sample_x];
    };
    combine_rgb(pole_tap(x0, y0, poles, OUTSIDE, load),
                pole_tap(x1, y0, poles, OUTSIDE, load),
                pole_tap(x0, y0 + 1, poles, OUTSIDE, load),
                pole_tap(x1, y0 + 1, poles, OUTSIDE, load), fx, fy, r, g, b);
  }

  constexpr int maximum_longitude_samples() const {
    return equator_samples > 0
               ? equator_samples
               : std::max(static_cast<int>(2.0f * math::PI_F *
                                               Geometry::ROWS_PER_RADIAN /
                                               spacing +
                                           0.5f),
                          1);
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
template <typename Value, int W, int H, int HOffset = -1> class SphericalField {
public:
  using Layout = SphericalFieldLayout<W, H, HOffset>;
  using Ring = typename Layout::Ring;

  constexpr SphericalField(Value *values, const Layout &layout)
      : layout(layout), values(values) {}

  constexpr int sample_count() const { return layout.sample_count(); }

  /**
   * @brief Fills an inclusive range of rings from a per-sample callback.
   * @tparam Populate Callable
   *   (const Vector &position, const Layout::Coordinates &), or the same with
   *   a trailing int receiving the sample's absolute index into data().
   * @param ring_begin First ring index to fill.
   * @param ring_end Last ring index to fill; must be below ring_count().
   * @param populate_sample Returns the Value stored at each sample.
   * @details Walks each ring by an incremental rotation, so `position` is only
   *   approximately unit and drifts from the exact-trig sample_vector(); a
   *   callback that needs the exact vector calls sample_vector() itself.
   */
  template <typename Populate>
  void populate(int ring_begin, int ring_end, Populate &&populate_sample) {
    HS_CHECK(ring_end < layout.ring_count(),
             "SphericalField::populate: ring_end %d past the last ring %d",
             ring_end, layout.ring_count() - 1);
    Ring ring = layout.ring(ring_begin);
    for (int ring_index = ring_begin; ring_index <= ring_end;
         ++ring_index, ring = layout.next_ring(ring)) {
      const math::Vector meridian = layout.sample_vector(ring, 0);
      const float theta_step = 2.0f * math::PI_F / ring.samples;
      const float step_cos = cosf(theta_step);
      const float step_sin = sinf(theta_step);
      float theta_cos = 1.0f;
      float theta_sin = 0.0f;
      for (int sample = 0; sample < ring.samples; ++sample) {
        const math::Vector position(meridian.x * theta_cos, meridian.y,
                                    meridian.x * theta_sin);
        const int index = ring.offset + sample;
        if constexpr (std::is_invocable_v<Populate &, const math::Vector &,
                                          const typename Layout::Coordinates &,
                                          int>)
          values[index] = populate_sample(
              position, layout.sample_coordinates(ring, sample), index);
        else
          values[index] = populate_sample(
              position, layout.sample_coordinates(ring, sample));
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
