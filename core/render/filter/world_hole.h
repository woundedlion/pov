/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "render/filter/pipeline.h"
#include "math/geometry.h"
#include "color/color.h"

/**
 * @file world_hole.h
 * @brief Filter::World::Hole: attenuates alpha inside an angular radius of a
 * center direction, leaving world positions unmoved.
 */

namespace Filter {
namespace World {

/**
 * @brief Creates a spherical hole by masking points within a radius.
 */
class Hole : public Is3D {
public:
  static constexpr bool requires_unit_world_input = true;
  /** @brief Attenuates alpha only; world points pass through unmoved. */
  static constexpr bool world_transform_is_identity = true;
  /**
   * @brief Constructs a hole centered at @p origin with angular @p radius.
   * @param origin Center of the hole (unit vector).
   * @param radius Angular radius of the hole in radians.
   */
  Hole(const Vector &origin, float radius) : origin(origin), radius(radius) {}

  /** @brief Moves the hole center to a new unit vector. */
  void set_origin(const Vector &new_origin) { origin = new_origin; }

  /**
   * @brief Changes the hole's angular radius in radians.
   * @param new_radius Angular radius; values <= 0 disable the hole.
   */
  void set_radius(float new_radius) { radius = new_radius; }
  /**
   * @brief Attenuates points near the hole center, leaving others unchanged.
   * @param v World-space point to test.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled by a quintic falloff inside the
   * radius, and the tap is dropped entirely once the falloff reaches zero.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    float d = angle_between(v, origin);
    if (d >= radius)
      pass(v, color, age, alpha);
    else {
      float mask = quintic_kernel(d / radius);
      if (mask > 1e-8f)
        pass(v, color, age, alpha * mask);
    }
  }

private:
  Vector origin; /**< Center of the hole (unit vector). */
  float radius;  /**< Angular radius of the hole in radians. */
};

} // namespace World
} // namespace Filter
