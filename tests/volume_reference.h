/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cfloat>
#include "core/math/3dmath.h"

namespace hs_test::scan_tests {

/** @brief Dense fixed-step oracle, independent of sphere-tracing policy. */
struct VolumeReference {
  struct Sample {
    math::Vector position;
    float distance;
  };

  template <typename Shape, typename Visit>
  static void march(const Shape &shape, const math::Vector &origin,
                    const math::Vector &direction, float radius, Visit visit) {
    constexpr int STEPS = 16384;
    const float END = radius - math::dot(origin, direction);
    for (int i = 0; i <= STEPS; ++i) {
      const float T = END * static_cast<float>(i) / STEPS;
      const math::Vector P(origin.x + direction.x * T,
                           origin.y + direction.y * T,
                           origin.z + direction.z * T);
      if (!visit(Sample{P, shape.distance(P)}))
        break;
    }
  }

  static float coverage(float distance, float threshold, float width) {
    const double X =
        std::clamp((width - distance) / (width - threshold), 0.0f, 1.0f);
    return static_cast<float>(X * X * X * (10.0 + X * (-15.0 + 6.0 * X)));
  }

  template <typename Shape>
  static Sample trace(const Shape &shape, const math::Vector &origin,
                      const math::Vector &direction, float radius, float aa) {
    Sample closest{origin, FLT_MAX};
    march(shape, origin, direction, radius, [&](Sample sample) {
      if (sample.distance < closest.distance)
        closest = sample;
      else if (closest.distance < aa)
        return false;
      return closest.distance > 0.02f * aa;
    });
    return closest;
  }

  template <typename Shape>
  static Sample behind(const Shape &shape, const math::Vector &origin,
                       const math::Vector &direction, float radius) {
    Sample closest{origin, FLT_MAX};
    float previous = FLT_MAX;
    bool climbing = false;
    march(shape, origin, direction, radius, [&](Sample sample) {
      if (sample.distance > previous)
        climbing = true;
      if (climbing && sample.distance < previous &&
          sample.distance < closest.distance)
        closest = sample;
      previous = sample.distance;
      return true;
    });
    return closest;
  }
};

// Vector-accumulating baseline for scalar ray-state regression tests.
struct VolumeScalarRegression {
  static constexpr float OVERRELAX_OMEGA = 1.3f;
  static constexpr int PROBE_STEPS = 24;
  static constexpr int PROBE_NEAR_STEPS = 6;
  static constexpr float PROBE_FLOOR_NEAR = 0.04f;
  static constexpr float PROBE_FLOOR_FAR = 0.12f;
  template <typename Shape>
  static __attribute__((always_inline)) float
  trace_closest(const Shape &shape, const math::Vector &local_ro,
                const math::Vector &local_vd, float bounds_radius,
                int max_steps, float aa_width, math::Vector &closest_local) {
    HS_PROFILE_DEEP(vol_trace);
    math::Vector local_p = local_ro;
    closest_local = local_ro;
    float closest_d = FLT_MAX;
    float omega = OVERRELAX_OMEGA;
    float prev_r = 0.0f;
    float step_len = 0.0f;
    for (int i = 0; i < max_steps; ++i) {
      if (local_p.x * local_vd.x + local_p.y * local_vd.y +
              local_p.z * local_vd.z >
          bounds_radius)
        break;
      float d = shape.distance(local_p);
      float r = d < 0.0f ? -d : d;
      if (omega > 1.0f && r + prev_r < step_len) {
        float back = prev_r - step_len;
        local_p = math::Vector(local_p.x + local_vd.x * back,
                               local_p.y + local_vd.y * back,
                               local_p.z + local_vd.z * back);
        omega = 1.0f;
        prev_r = 0.0f;
        step_len = 0.0f;
        continue;
      }
      prev_r = r;
      if (d < closest_d) {
        closest_d = d;
        closest_local = local_p;
        if (closest_d <= aa_width * 0.02f)
          break;
      } else if (closest_d < aa_width) {
        break;
      }
      if (d < -aa_width)
        break;
      step_len = std::max(d * 0.9f * omega, 1e-5f);
      local_p = math::Vector(local_p.x + local_vd.x * step_len,
                             local_p.y + local_vd.y * step_len,
                             local_p.z + local_vd.z * step_len);
    }
    return closest_d;
  }
  struct Occluder {
    bool solid;
    math::Vector behind;
    float distance;
    float soft;
  };
  template <typename Shape>
  static __attribute__((always_inline)) Occluder
  probe_occluder(const Shape &shape, const math::Vector &closest_local,
                 const math::Vector &local_vd, float bounds_radius,
                 float hit_threshold, float aa_width) {
    HS_PROFILE_DEEP(vol_probe);
    math::Vector probe = closest_local;
    float prev = FLT_MAX;
    bool climbing = false;
    float min_behind = FLT_MAX;
    math::Vector min_pos = closest_local;
    float s = 0.0f, prev_s = 0.0f, min_s = 0.0f;
    float bef_s = 0.0f, bef_pd = FLT_MAX;
    float aft_s = 0.0f, aft_pd = FLT_MAX;
    bool need_aft = false;
    for (int i = 0; i < PROBE_STEPS; ++i) {
      if (probe.x * local_vd.x + probe.y * local_vd.y + probe.z * local_vd.z >
          bounds_radius)
        break;
      float pd = shape.distance(probe);
      if (pd < hit_threshold)
        return {true, probe, pd, 0.0f};
      if (need_aft) {
        aft_s = s;
        aft_pd = pd;
        need_aft = false;
      }
      if (pd > prev)
        climbing = true;
      else if (climbing && pd < min_behind) {
        min_behind = pd;
        min_pos = probe;
        min_s = s;
        bef_s = prev_s;
        bef_pd = prev;
        need_aft = true;
        aft_pd = FLT_MAX;
      }
      prev = pd;
      prev_s = s;
      float floor = bounds_radius *
                    (i < PROBE_NEAR_STEPS ? PROBE_FLOOR_NEAR : PROBE_FLOOR_FAR);
      float step = std::max(pd * 0.9f, floor);
      probe =
          math::Vector(probe.x + local_vd.x * step, probe.y + local_vd.y * step,
                       probe.z + local_vd.z * step);
      s += step;
    }
    if (min_behind < 2.0f * aa_width && bef_pd != FLT_MAX &&
        aft_pd != FLT_MAX) {
      float p = min_s - bef_s;
      float q = min_s - aft_s;
      float yb = bef_pd - min_behind;
      float ya = aft_pd - min_behind;
      float den = q * yb - p * ya;
      if (den < -1e-12f) {
        float ds = -0.5f * (q * q * yb - p * p * ya) / den;
        if (ds > -p && ds < -q) {
          math::Vector rp(min_pos.x + local_vd.x * ds,
                          min_pos.y + local_vd.y * ds,
                          min_pos.z + local_vd.z * ds);
          float rpd = shape.distance(rp);
          if (rpd < min_behind) {
            min_behind = rpd;
            min_pos = rp;
          }
          if (min_behind < hit_threshold)
            return {true, min_pos, min_behind, 0.0f};
        }
      }
    }
    float soft =
        (min_behind < aa_width)
            ? VolumeReference::coverage(min_behind, hit_threshold, aa_width)
            : 0.0f;
    return {false, min_pos, min_behind, soft};
  }
};

} // namespace hs_test::scan_tests
