/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file spider_web.h
 * @brief Orb-web thread field in a gnomonic plane centred on the hub. */

#include "core/math/fast_math.h"
#include "core/math/scalar_kernels.h"
#include "core/math/scalar_noise.h"
#include <array>
#include <cfloat>
#include <cmath>

#include "core/pattern/field.h"

namespace Pattern {

/** @brief Passive orb-web controls. */
struct SpiderWebParams {
  float web_radius = 1.2f; ///< Outermost capture turn, plane units.
  float spokes = 24.0f;    ///< Rounded radial thread count.
  float turns = 16.0f;     ///< Capture spiral turns from hub to rim.
  float hub_radius = 0.2f; ///< Capture-free zone as a fraction of web radius.
  float sag = 0.1f;        ///< Capture-chord bow toward the hub, chord lengths.
  float jitter = 0.3f;     ///< Spoke-angle and turn-spacing irregularity.
  float dew = 0.4f;        ///< Fraction of capture-chord slots holding a drop.
  float seed = 0.0f;       ///< Rounded hash stream for jitter and dew.
  float line_thickness = 0.006f; ///< Thread half-width, radians.
  float line_softness = 0.01f;   ///< Thread edge fade width, radians.
};

/**
 * @brief Orb web of straight geodesic threads in a gnomonic plane.
 * @details The plane origin is the hub. Great circles project to straight
 * lines, so spokes are rays and capture threads are chords between adjacent
 * spokes. Widths are angular under a gnomonic projection.
 */
struct PlaneSpiderWeb {
  using Coordinate = PlanePosition;     ///< Plane position; origin at the hub.
  using Params = SpiderWebParams;       ///< Web layout and thread controls.
  using State = EmptyState;             ///< Stationary source state.
  static constexpr int MAX_SPOKES = 48; ///< Maximum radial thread count.
  static constexpr float MAX_TURNS = 64.0f; ///< Maximum capture turns.

  /** @brief Spoke table and derived web radii. */
  struct Prepared {
    std::array<float, MAX_SPOKES> dir_x{}; ///< Unit spoke direction x.
    std::array<float, MAX_SPOKES> dir_y{}; ///< Unit spoke direction y.
    /// Counterclockwise spoke angles; entry N repeats entry 0 plus a turn.
    std::array<float, MAX_SPOKES + 1> angle{};
    int spokes = 0;              ///< Active spokes.
    int crossings = 0;           ///< Capture-spiral spoke crossings.
    float inner_radius = 0.0f;   ///< Capture spiral start radius.
    float step = 0.0f;           ///< Radius gained per spoke crossing.
    float radius_jitter = 0.0f;  ///< Peak crossing-radius offset.
    float frame_radius = 0.0f;   ///< Radius of the taut outer frame.
    float hub_core = 0.0f;       ///< Outer hub-ring radius.
    float reach = 0.0f;          ///< Radius beyond which nothing lights.
    float sag = 0.0f;            ///< Capture-chord bow in chord lengths.
    float dew = 0.0f;            ///< Drop occupancy.
    float line_thickness = 0.0f; ///< Thread half-width, radians.
    float line_softness = 0.0f;  ///< Thread edge fade width, radians.
    uint32_t seed = 0;           ///< Hash stream.
  };

  /** @brief Checks the finite bounded web recipe.
   * @param params Web controls.
   * @return Whether preparation admits the recipe. */
  static bool valid(const Params &params) {
    return parameter_in_range(params.web_radius, 0.05f, 16.0f) &&
           parameter_in_range(params.spokes, 3.0f,
                              static_cast<float>(MAX_SPOKES)) &&
           parameter_in_range(params.turns, 1.0f, MAX_TURNS) &&
           parameter_in_range(params.hub_radius, 0.0f, 0.9f) &&
           parameter_in_range(params.sag, 0.0f, 0.5f) &&
           parameter_in_range(params.jitter, 0.0f, 1.0f) &&
           parameter_in_range(params.dew, 0.0f, 1.0f) &&
           parameter_in_range(params.seed, 0.0f, 65535.0f) &&
           parameter_in_range(params.line_thickness, 0.0f, 0.1f) &&
           parameter_in_range(params.line_softness, 1.0f / 1024.0f, 0.1f);
  }

  /** @brief Builds the jittered spoke table and spiral radii.
   * @param params Web controls.
   * @param state Empty source state.
   * @return Immutable web. */
  static Prepared prepare(const Params &params, const State &state) {
    (void)state;
    HS_CHECK(valid(params), "pattern spider web: invalid parameters");
    Prepared out;
    out.seed = static_cast<uint32_t>(lroundf(params.seed));
    out.spokes = static_cast<int>(lroundf(params.spokes));
    out.crossings = static_cast<int>(
        lroundf(params.turns * static_cast<float>(out.spokes)));
    const float sector = math::TWO_PI_F / static_cast<float>(out.spokes);
    for (int k = 0; k < out.spokes; ++k) {
      const float offset =
          SPOKE_JITTER * params.jitter *
          (math::hash01(static_cast<uint32_t>(k), out.seed) - 0.5f);
      out.angle[k] = sector * (static_cast<float>(k) + offset);
      out.dir_x[k] = cosf(out.angle[k]);
      out.dir_y[k] = sinf(out.angle[k]);
    }
    out.angle[out.spokes] = out.angle[0] + math::TWO_PI_F;
    out.inner_radius = params.hub_radius * params.web_radius;
    out.step = (params.web_radius - out.inner_radius) /
               static_cast<float>(out.crossings);
    const float turn_spacing = out.step * static_cast<float>(out.spokes);
    out.radius_jitter = RADIUS_JITTER * params.jitter * turn_spacing;
    out.frame_radius =
        params.web_radius + out.radius_jitter + FRAME_GAP * turn_spacing;
    out.hub_core = HUB_CORE * out.inner_radius;
    out.sag = params.sag;
    for (int k = 0; k < out.spokes; ++k) {
      const int b = k + 1 == out.spokes ? 0 : k + 1;
      const float cosine =
          out.dir_x[k] * out.dir_x[b] + out.dir_y[k] * out.dir_y[b];
      const float sine =
          out.dir_x[k] * out.dir_y[b] - out.dir_y[k] * out.dir_x[b];
      for (int c = k; c < out.crossings; c += out.spokes) {
        const float ra = crossing_radius(out, c);
        const float rb = crossing_radius(out, c + 1);
        const float near = fminf(ra, rb), far = fmaxf(ra, rb);
        const float denominator = 4.0f * (far - near * cosine);
        if (denominator > 0.0f)
          out.sag = fminf(out.sag, near * sine / denominator);
      }
    }
    out.dew = params.dew;
    out.line_thickness = params.line_thickness;
    out.line_softness = params.line_softness;
    const float cutoff =
        2.0f * (params.line_thickness * FRAME_WEIGHT + params.line_softness);
    out.reach =
        cutoff * out.frame_radius < 1.0f
            ? (out.frame_radius + cutoff) / (1.0f - cutoff * out.frame_radius)
            : FLT_MAX;
    return out;
  }

  /** @brief Samples threads and dew drops.
   * @param coordinate Plane position; the hub sits at the origin.
   * @param prepared Immutable web.
   * @return 1 on a thread or drop core, -1 off the web, ramped between. */
  HS_HOT_FLASH_MEMBER static float sample(const Coordinate &coordinate,
                                          const Prepared &prepared) {
    const float x = coordinate.value.re, y = coordinate.value.im;
    const float rho_sq = x * x + y * y;
    if (!(rho_sq < prepared.reach * prepared.reach))
      return -1.0f;
    const float rho = sqrtf(rho_sq);
    const Scales scales{sqrtf(1.0f + rho_sq), rho};
    const int n = prepared.spokes;

    float phi = math::fast_atan2(y, x);
    if (phi < prepared.angle[0])
      phi += math::TWO_PI_F;
    int k = static_cast<int>(phi * static_cast<float>(n) / math::TWO_PI_F);
    k = k < 0 ? 0 : (k >= n ? n - 1 : k);
    for (int pass = 0; pass < 2; ++pass) {
      if (cross(prepared, k, x, y) < 0.0f)
        k = k == 0 ? n - 1 : k - 1;
      else if (cross(prepared, k + 1 == n ? 0 : k + 1, x, y) >= 0.0f)
        k = k + 1 == n ? 0 : k + 1;
      else
        break;
    }
    const int b = k + 1 == n ? 0 : k + 1;
    const float cross_a = cross(prepared, k, x, y);
    const float cross_b = -cross(prepared, b, x, y);

    const float overshoot = fmaxf(rho - prepared.frame_radius, 0.0f);
    const float spoke_across = fminf(fabsf(cross_a), fabsf(cross_b));
    float thread = hypotf(spoke_across / scales.tangential,
                          overshoot / (1.0f + rho * prepared.frame_radius));

    const float frame = chord(prepared, k, b, prepared.frame_radius,
                              prepared.frame_radius, 0.0f, x, y, scales)
                            .distance;
    thread =
        fminf(thread, frame - prepared.line_thickness * (FRAME_WEIGHT - 1.0f));

    const float margin = prepared.reach - prepared.frame_radius;
    if (prepared.hub_core > 0.0f && rho < prepared.hub_core + margin) {
      for (int ring = 1; ring <= HUB_RINGS; ++ring) {
        const float radius = prepared.hub_core * static_cast<float>(ring) /
                             static_cast<float>(HUB_RINGS);
        thread = fminf(
            thread,
            chord(prepared, k, b, radius, radius, 0.0f, x, y, scales).distance);
      }
    }

    const float width = prepared.angle[k + 1] - prepared.angle[k];
    const float u_est = cross_a / fmaxf(cross_a + cross_b, 1e-12f);
    const float bow = 4.0f * u_est * (1.0f - u_est) *
                      (width * width * 0.125f + prepared.sag * width);
    const float linear = rho / fmaxf(1.0f - bow, 0.25f);
    const float turn_est = ((linear - prepared.inner_radius) / prepared.step -
                            static_cast<float>(k) - u_est) /
                           static_cast<float>(n);
    const int first = static_cast<int>(floorf(turn_est)) - 1;
    Chord best{1e30f, 0.0f, 0.0f, 0.0f};
    int best_crossing = -1;
    for (int turn = first; turn <= first + 3; ++turn) {
      const int c = turn * n + k;
      if (c < 0 || c + 1 > prepared.crossings)
        continue;
      const Chord capture =
          chord(prepared, k, b, crossing_radius(prepared, c),
                crossing_radius(prepared, c + 1), prepared.sag, x, y, scales);
      if (capture.distance < best.distance) {
        best = capture;
        best_crossing = c;
      }
    }
    thread = fminf(thread, best.distance);

    const float edge = prepared.line_thickness;
    float value =
        1.0f - math::smooth_ramp(edge, edge + prepared.line_softness, thread);
    if (best_crossing >= 0 && prepared.dew > 0.0f)
      value = fmaxf(value, drop(prepared, best_crossing, best));
    return 2.0f * value - 1.0f;
  }

private:
  static constexpr float SPOKE_JITTER = 0.3f;  ///< Peak offset, sectors.
  static constexpr float RADIUS_JITTER = 0.2f; ///< Peak offset, turn spacings.
  static constexpr float FRAME_GAP =
      0.75f; ///< Rim-to-frame gap, turn spacings.
  static constexpr float FRAME_WEIGHT =
      1.5f;                               ///< Frame width over thread width.
  static constexpr float HUB_CORE = 0.4f; ///< Hub rings' extent in hub radii.
  static constexpr int HUB_RINGS = 3;     ///< Concentric hub rings.
  static constexpr float DROP_SIZE =
      3.5f; ///< Largest drop radius, thread widths.
  static constexpr float DROP_PITCH =
      3.0f; ///< Drop slot length, largest radii.
  static constexpr uint32_t MAX_DROP_SLOTS = 4096u; ///< Drop slots per chord.
  static constexpr uint32_t CROSSING_STREAM = 0x51EDu; ///< Radius jitter.
  static constexpr uint32_t DROP_STREAM = 0xD3E5u;     ///< Drop occupancy.
  static constexpr uint32_t SIZE_STREAM = 0x2C1Bu;     ///< Drop radius.

  /** @brief Gnomonic plane-to-angle divisors at the sample. */
  struct Scales {
    float tangential; ///< Plane length per radian across a spoke.
    float radius;     ///< Sample radius in plane units.
  };

  /** @brief Nearest point on one sagging chord. */
  struct Chord {
    float distance; ///< Angular distance to the chord, radians.
    float u;        ///< Clamped chord parameter in [0, 1].
    float across;   ///< Angular offset across the chord, radians.
    float length;   ///< Angular chord length at its midpoint, radians.
  };

  /** @brief Signed plane distance from a spoke's line.
   * @param prepared Immutable web.
   * @param k Spoke index.
   * @param x Plane x.
   * @param y Plane y.
   * @return Distance, positive counterclockwise of the spoke. */
  static float cross(const Prepared &prepared, int k, float x, float y) {
    return prepared.dir_x[k] * y - prepared.dir_y[k] * x;
  }

  /** @brief Jittered radius of one capture-spiral spoke crossing.
   * @param prepared Immutable web.
   * @param c Crossing index from the hub outward.
   * @return Nonnegative plane radius. */
  static float crossing_radius(const Prepared &prepared, int c) {
    const float offset = math::hash01(static_cast<uint32_t>(c),
                                      prepared.seed + CROSSING_STREAM) -
                         0.5f;
    return fmaxf(prepared.inner_radius + prepared.step * static_cast<float>(c) +
                     2.0f * prepared.radius_jitter * offset,
                 0.0f);
  }

  /** @brief Nearest point on a chord between two spokes, bowed toward the hub.
   * @param prepared Immutable web.
   * @param a Start spoke.
   * @param b End spoke, counterclockwise of @p a.
   * @param ra Start radius on @p a.
   * @param rb End radius on @p b.
   * @param sag Bow depth at the chord middle, chord lengths.
   * @param x Plane x.
   * @param y Plane y.
   * @param scales Plane-to-angle divisors at the sample.
   * @return Angular distance, parameter, offset and length. */
  static Chord chord(const Prepared &prepared, int a, int b, float ra, float rb,
                     float sag, float x, float y, const Scales &scales) {
    const float ax = ra * prepared.dir_x[a], ay = ra * prepared.dir_y[a];
    const float dx = rb * prepared.dir_x[b] - ax;
    const float dy = rb * prepared.dir_y[b] - ay;
    const float length = sqrtf(dx * dx + dy * dy);
    const bool degenerate = !(length > 1e-9f);
    const float inv = degenerate ? 0.0f : 1.0f / length;
    const float tx = degenerate ? 1.0f : dx * inv;
    const float ty = degenerate ? 0.0f : dy * inv;
    const float qx = x - ax, qy = y - ay;
    const float along = qx * tx + qy * ty;
    const float u = hs::clamp(along * inv, 0.0f, 1.0f);
    const float bow = sag * length * 4.0f * u * (1.0f - u);
    const float offset = qx * ty - qy * tx + bow;
    const float across =
        offset / (1.0f + scales.radius * (scales.radius - offset));
    const float over = (along - u * length) / scales.tangential;
    const float mid = 0.5f * (ra + rb);
    return {hypotf(across, over), u, across, length / sqrtf(1.0f + mid * mid)};
  }

  /** @brief Dew drop coverage around the nearest capture chord.
   * @param prepared Immutable web.
   * @param crossing Crossing index the chord starts at.
   * @param near Nearest point on that chord.
   * @return Drop coverage in [0, 1]. */
  static float drop(const Prepared &prepared, int crossing, const Chord &near) {
    const float largest =
        DROP_SIZE * (prepared.line_thickness + prepared.line_softness);
    const float slots_f =
        hs::clamp(floorf(near.length / (DROP_PITCH * largest)), 1.0f,
                  static_cast<float>(MAX_DROP_SLOTS));
    const int slot = static_cast<int>(fminf(near.u * slots_f, slots_f - 1.0f));
    const uint32_t id = static_cast<uint32_t>(crossing) * MAX_DROP_SLOTS +
                        static_cast<uint32_t>(slot);
    if (math::hash01(id, prepared.seed + DROP_STREAM) >= prepared.dew)
      return 0.0f;
    const float radius =
        largest * (0.4f + 0.6f * math::hash01(id, prepared.seed + SIZE_STREAM));
    const float center = (static_cast<float>(slot) + 0.5f) / slots_f;
    const float along = (near.u - center) * near.length;
    const float distance = hypotf(along, near.across);
    return 1.0f - math::smooth_ramp(0.35f * radius,
                                    radius + prepared.line_softness, distance);
  }
};

} // namespace Pattern
