/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file HyperLattice.h
 * @brief Reflective flight through cubic and four-dimensional hypercubic
 *        lattices.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <tuple>

#include "core/color/effect_palette_recipes.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/render/pullback.h"

namespace hs_test {
namespace hyper_lattice_tests {
struct HyperLatticeWhiteBox;
} // namespace hyper_lattice_tests
} // namespace hs_test

namespace HyperLatticeDetail {

constexpr int DIMENSIONS = 4;
constexpr int MAX_SHELLS = 3;
constexpr float DIRECTION_EPSILON = 1.0e-4f;

enum class ReflectionMode : uint8_t { CHROME, RADIAL };
enum class ColorMode : uint8_t { DEPTH, AXIS };
enum class ShellCount : uint8_t { ONE, TWO, THREE };

struct Vec4 {
  float v[DIMENSIONS]{};

  constexpr float &operator[](int index) { return v[index]; }
  constexpr float operator[](int index) const { return v[index]; }
};

struct Mat4 {
  float m[DIMENSIONS][DIMENSIONS]{};

  static constexpr Mat4 identity() {
    Mat4 result;
    for (int i = 0; i < DIMENSIONS; ++i)
      result.m[i][i] = 1.0f;
    return result;
  }

  Vec4 apply(const Vec4 &input) const {
    Vec4 result;
    for (int row = 0; row < DIMENSIONS; ++row)
      result[row] = m[row][0] * input[0] + m[row][1] * input[1] +
                    m[row][2] * input[2] + m[row][3] * input[3];
    return result;
  }
};

inline void rotate_plane(Mat4 &matrix, int a, int b, float angle) {
  const float c = fast_cosf(angle);
  const float s = fast_sinf(angle);
  for (int column = 0; column < DIMENSIONS; ++column) {
    const float av = matrix.m[a][column];
    const float bv = matrix.m[b][column];
    matrix.m[a][column] = c * av - s * bv;
    matrix.m[b][column] = s * av + c * bv;
  }
}

inline float periodic_distance(float value) {
  return fabsf(value - floorf(value + 0.5f));
}

struct EdgeMetric {
  float distance_sq;
  uint8_t free_axis;
};

inline EdgeMetric edge_metric_3d(const Vec4 &point, int plane_axis) {
  float best_fixed = 1.0f;
  float farthest = -1.0f;
  uint8_t free_axis = 0;
  for (int axis = 0; axis < 3; ++axis) {
    if (axis == plane_axis)
      continue;
    const float distance = periodic_distance(point[axis]);
    const float distance_sq = distance * distance;
    best_fixed = std::min(best_fixed, distance_sq);
    if (distance_sq > farthest) {
      farthest = distance_sq;
      free_axis = static_cast<uint8_t>(axis);
    }
  }
  return {best_fixed, free_axis};
}

inline EdgeMetric edge_metric_4d(const Vec4 &point, int plane_axis) {
  float sum = 0.0f;
  float largest = -1.0f;
  uint8_t free_axis = 0;
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    if (axis == plane_axis)
      continue;
    const float distance = periodic_distance(point[axis]);
    const float distance_sq = distance * distance;
    sum += distance_sq;
    if (distance_sq > largest) {
      largest = distance_sq;
      free_axis = static_cast<uint8_t>(axis);
    }
  }
  return {sum - largest, free_axis};
}

struct Params {
  float dimension = 0.0f;
  float sphere_radius = 0.4f;
  float wire_radius = 0.055f;
  float softness = 0.012f;
  float far_cells = 7.0f;
  float aa_strength = 1.0f;
  float speed = 0.018f;
  float spin_3d = 0.0024f;
  float spin_4d = 0.0f;
  ReflectionMode reflection = ReflectionMode::CHROME;
  ColorMode color = ColorMode::DEPTH;
  ShellCount shells = ShellCount::TWO;

  void lerp(const Params &start, const Params &target, float amount) {
    dimension = hs::lerp(start.dimension, target.dimension, amount);
    sphere_radius = hs::lerp(start.sphere_radius, target.sphere_radius, amount);
    wire_radius = hs::lerp(start.wire_radius, target.wire_radius, amount);
    softness = hs::lerp(start.softness, target.softness, amount);
    far_cells = hs::lerp(start.far_cells, target.far_cells, amount);
    aa_strength = hs::lerp(start.aa_strength, target.aa_strength, amount);
    speed = hs::lerp(start.speed, target.speed, amount);
    spin_3d = hs::lerp(start.spin_3d, target.spin_3d, amount);
    spin_4d = hs::lerp(start.spin_4d, target.spin_4d, amount);
    reflection = amount < 0.5f ? start.reflection : target.reflection;
    color = amount < 0.5f ? start.color : target.color;
    shells = amount < 0.5f ? start.shells : target.shells;
  }
};

struct FrameState {
  Params params;
  Vec4 origin;
  std::array<float, 6> rotation_phase;
  float pixel_half_angle;
  const BakedPalette *palette;
};

struct Binding {
  using FrameState = HyperLatticeDetail::FrameState;
  using Instrumentation = Pullback::NoInstrumentation;
};

struct PreparedTrace {
  Params params;
  Vec4 origin;
  Mat4 world_to_lattice;
  float inv_far;
  float pixel_half_angle;
};

template <int W, int H> constexpr float pixel_half_angle() {
  constexpr float HORIZONTAL = TWO_PI_F / static_cast<float>(W);
  constexpr float VERTICAL = PI_F / static_cast<float>(H + hs::H_OFFSET - 1);
  return 0.5f * (HORIZONTAL > VERTICAL ? HORIZONTAL : VERTICAL);
}

inline float wire_coverage(float metric_sq, float radius, float half_width) {
  const float signed_distance = sqrtf(metric_sq) - radius;
  return 1.0f - smooth_ramp(-half_width, half_width, signed_distance);
}

inline float projected_half_width(float distance, float plane_component,
                                  float pixel_half_angle, float aa_strength,
                                  float softness) {
  return softness +
         aa_strength * pixel_half_angle * distance / fabsf(plane_component);
}

inline float near_field_coverage(float distance, float wire_radius) {
  return smooth_ramp(1.5f * wire_radius, 4.0f * wire_radius, distance);
}

inline float next_plane_offset(float origin, bool positive) {
  const float fraction = wrap_t(origin);
  if (fraction == 0.0f)
    return 1.0f;
  return positive ? 1.0f - fraction : fraction;
}

inline PreparedTrace prepare_trace(const FrameState &frame) {
  PreparedTrace prepared;
  prepared.params = frame.params;
  prepared.origin = frame.origin;
  prepared.world_to_lattice = Mat4::identity();
  rotate_plane(prepared.world_to_lattice, 0, 1, frame.rotation_phase[0]);
  rotate_plane(prepared.world_to_lattice, 0, 2, frame.rotation_phase[1]);
  rotate_plane(prepared.world_to_lattice, 1, 2, frame.rotation_phase[2]);
  const float dimension = frame.params.dimension;
  rotate_plane(prepared.world_to_lattice, 0, 3,
               dimension * frame.rotation_phase[3]);
  rotate_plane(prepared.world_to_lattice, 1, 3,
               dimension * frame.rotation_phase[4]);
  rotate_plane(prepared.world_to_lattice, 2, 3,
               dimension * frame.rotation_phase[5]);
  prepared.inv_far = 1.0f / frame.params.far_cells;
  prepared.pixel_half_angle = frame.pixel_half_angle;
  return prepared;
}

inline Vec4 reflected_direction(const Vector &normal, ReflectionMode mode) {
  if (mode == ReflectionMode::RADIAL)
    return {{normal.x, normal.y, normal.z, 0.0f}};
  const float twice_x = 2.0f * normal.x;
  return {{twice_x * normal.x - 1.0f, twice_x * normal.y, twice_x * normal.z,
           0.0f}};
}

struct TraceHit {
  float coverage = 0.0f;
  float distance = 0.0f;
  uint8_t free_axis = 0;
};

inline TraceHit trace_plane(const Vec4 &ray_origin, const Vec4 &direction,
                            int plane_axis, float distance,
                            const PreparedTrace &prepared) {
  Vec4 point;
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    point[axis] = ray_origin[axis] + distance * direction[axis];

  float metric_sq;
  uint8_t free_axis;
  float dimensional_coverage = 1.0f;
  if (prepared.params.dimension == 0.0f) {
    const EdgeMetric metric_3d = edge_metric_3d(point, plane_axis);
    metric_sq = metric_3d.distance_sq;
    free_axis = metric_3d.free_axis;
  } else {
    const EdgeMetric metric_4d = edge_metric_4d(point, plane_axis);
    metric_sq = metric_4d.distance_sq;
    free_axis = metric_4d.free_axis;
    if (plane_axis < 3 && prepared.params.dimension < 1.0f) {
      const EdgeMetric metric_3d = edge_metric_3d(point, plane_axis);
      metric_sq = hs::lerp(metric_3d.distance_sq, metric_4d.distance_sq,
                           prepared.params.dimension);
      if (prepared.params.dimension < 0.5f)
        free_axis = metric_3d.free_axis;
    } else if (plane_axis == 3) {
      dimensional_coverage = prepared.params.dimension;
    }
  }

  const float half_width = projected_half_width(
      distance, direction[plane_axis], prepared.pixel_half_angle,
      prepared.params.aa_strength, prepared.params.softness);
  const float edge =
      wire_coverage(metric_sq, prepared.params.wire_radius, half_width);
  const float fog = std::max(0.0f, 1.0f - distance * prepared.inv_far);
  const float near = near_field_coverage(distance, prepared.params.wire_radius);
  return {edge * fog * fog * near * dimensional_coverage, distance, free_axis};
}

struct TraceCursor {
  float distance = 0.0f;
  float step = 0.0f;
  uint8_t shell = 0;
  bool active = false;
};

template <typename ConsumeFn>
inline void trace_layers(const Vector &normal, const PreparedTrace &prepared,
                         ConsumeFn consume) {
  constexpr float GROUP_EPSILON = 1.0e-4f;
  const Vec4 surface_normal =
      prepared.world_to_lattice.apply({{normal.x, normal.y, normal.z, 0.0f}});
  Vec4 ray_origin = prepared.origin;
  for (int axis = 0; axis < DIMENSIONS; ++axis)
    ray_origin[axis] += prepared.params.sphere_radius * surface_normal[axis];
  const Vec4 direction = prepared.world_to_lattice.apply(
      reflected_direction(normal, prepared.params.reflection));
  const uint8_t shell_count = static_cast<uint8_t>(prepared.params.shells) + 1;
  TraceCursor cursors[DIMENSIONS];
  for (int axis = 0; axis < DIMENSIONS; ++axis) {
    if (axis == 3 && prepared.params.dimension == 0.0f)
      continue;
    const float component = direction[axis];
    const float magnitude = fabsf(component);
    if (magnitude < DIRECTION_EPSILON)
      continue;
    TraceCursor &cursor = cursors[axis];
    cursor.step = 1.0f / magnitude;
    cursor.distance =
        next_plane_offset(ray_origin[axis], component > 0.0f) * cursor.step;
    cursor.active = cursor.distance < prepared.params.far_cells;
  }

  for (int event = 0; event < DIMENSIONS * MAX_SHELLS; ++event) {
    float nearest = prepared.params.far_cells;
    for (const TraceCursor &cursor : cursors)
      if (cursor.active && cursor.distance < nearest)
        nearest = cursor.distance;
    if (nearest >= prepared.params.far_cells)
      break;

    const float tolerance = GROUP_EPSILON * std::max(1.0f, nearest);
    TraceHit layer;
    for (int axis = 0; axis < DIMENSIONS; ++axis) {
      TraceCursor &cursor = cursors[axis];
      if (!cursor.active || fabsf(cursor.distance - nearest) > tolerance)
        continue;
      const TraceHit candidate =
          trace_plane(ray_origin, direction, axis, cursor.distance, prepared);
      if (candidate.coverage > layer.coverage)
        layer = candidate;
      ++cursor.shell;
      cursor.distance += cursor.step;
      cursor.active = cursor.shell < shell_count &&
                      cursor.distance < prepared.params.far_cells;
    }
    if (layer.coverage > 0.0f && !consume(layer))
      return;
  }
}

inline TraceHit trace(const Vector &normal, const PreparedTrace &prepared) {
  TraceHit nearest;
  trace_layers(normal, prepared, [&](const TraceHit &hit) {
    nearest = hit;
    return false;
  });
  return nearest;
}

struct LayerComposite {
  float red = 0.0f;
  float green = 0.0f;
  float blue = 0.0f;
  float remaining = 1.0f;

  void add(const Pixel &color, float coverage) {
    const float layer_coverage = hs::clamp(coverage, 0.0f, 1.0f);
    const float weight = remaining * layer_coverage;
    red += static_cast<float>(color.r) * weight;
    green += static_cast<float>(color.g) * weight;
    blue += static_cast<float>(color.b) * weight;
    remaining *= 1.0f - layer_coverage;
  }

  Color4 finish() const {
    const float alpha = 1.0f - remaining;
    if (alpha <= 0.0f)
      return {};
    const float inverse_alpha = 1.0f / alpha;
    const Pixel color{static_cast<uint16_t>(hs::clamp(
                          red * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          green * inverse_alpha + 0.5f, 0.0f, 65535.0f)),
                      static_cast<uint16_t>(hs::clamp(
                          blue * inverse_alpha + 0.5f, 0.0f, 65535.0f))};
    return {color, alpha};
  }
};

inline Color4 shade(const Pullback::SphereSample &input,
                    const FrameState &frame, const PreparedTrace &prepared) {
  LayerComposite composite;
  trace_layers(input.dir, prepared, [&](const TraceHit &hit) {
    const float path_length = input.path_length + hit.distance;
    const float depth = hs::clamp(path_length * prepared.inv_far, 0.0f, 1.0f);
    const float value =
        prepared.params.color == ColorMode::DEPTH
            ? 1.0f - depth
            : (static_cast<float>(hit.free_axis) + 0.75f * depth) / 4.0f;
    Color4 color = frame.palette->get(value);
    color.color = color.color * (0.45f + 0.55f * (1.0f - depth));
    composite.add(color.color, hit.coverage * color.alpha);
    return composite.remaining > MIN_ENCODABLE_ALPHA;
  });
  return composite.finish();
}

struct ShadeStage
    : Pullback::Stage::Contract<ShadeStage, Pullback::SphereSample, Color4> {
  using Policies = std::tuple<>;

  template <typename PipelineBinding>
  static PreparedTrace
  prepare(const typename PipelineBinding::FrameState &frame) {
    return prepare_trace(frame);
  }

  template <typename PipelineBinding>
  __attribute__((always_inline)) static Color4
  run(const Pullback::SphereSample &input,
      const typename PipelineBinding::FrameState &frame,
      const PreparedTrace &prepared) {
    return shade(input, frame, prepared);
  }
};

using RenderPipeline = Pullback::Pipeline<Binding, ShadeStage>;

} // namespace HyperLatticeDetail

template <int W, int H>
class HyperLattice : public ChoreographedEffect<HyperLattice<W, H>,
                                                HyperLatticeDetail::Params> {
  using Choreography =
      ChoreographedEffect<HyperLattice<W, H>, HyperLatticeDetail::Params>;
  friend Choreography;

public:
  using Params = HyperLatticeDetail::Params;
  using ReflectionMode = HyperLatticeDetail::ReflectionMode;
  using ColorMode = HyperLatticeDetail::ColorMode;
  using ShellCount = HyperLatticeDetail::ShellCount;

  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "cubic-flight", "deep-grid", "dimensional-rift", "hypercube-flight"};
  static constexpr Segue::Lerp PRESET_SEGUE{240, ease_in_out_sin,
                                            /*pausable=*/true};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;

  static constexpr Params initial_params() { return {}; }

  static constexpr Params preset_params(size_t index) {
    Params value;
    switch (index) {
    case 1:
      value.wire_radius = 0.035f;
      value.softness = 0.008f;
      value.far_cells = 12.0f;
      value.speed = 0.03f;
      value.reflection = ReflectionMode::RADIAL;
      value.shells = ShellCount::THREE;
      break;
    case 2:
      value.dimension = 0.55f;
      value.wire_radius = 0.085f;
      value.softness = 0.015f;
      value.far_cells = 9.0f;
      value.spin_4d = 0.004f;
      value.color = ColorMode::AXIS;
      break;
    case 3:
      value.dimension = 1.0f;
      value.wire_radius = 0.11f;
      value.softness = 0.018f;
      value.far_cells = 10.0f;
      value.speed = 0.024f;
      value.spin_4d = 0.0041f;
      value.color = ColorMode::AXIS;
      value.shells = ShellCount::THREE;
      break;
    default:
      break;
    }
    return value;
  }

  static bool valid_params(const Params &value) {
    return value.dimension >= 0.0f && value.dimension <= 1.0f &&
           value.sphere_radius >= 0.0f && value.sphere_radius <= 1.5f &&
           value.wire_radius >= 0.015f && value.wire_radius <= 0.18f &&
           value.softness >= 0.002f && value.softness <= 0.08f &&
           value.far_cells >= 2.0f && value.far_cells <= 16.0f &&
           value.aa_strength >= 0.0f && value.aa_strength <= 2.0f &&
           value.speed >= 0.0f && value.speed <= 0.05f &&
           value.spin_3d >= 0.0f && value.spin_3d <= 0.015f &&
           value.spin_4d >= 0.0f && value.spin_4d <= 0.015f &&
           static_cast<uint8_t>(value.reflection) <=
               static_cast<uint8_t>(ReflectionMode::RADIAL) &&
           static_cast<uint8_t>(value.color) <=
               static_cast<uint8_t>(ColorMode::AXIS) &&
           static_cast<uint8_t>(value.shells) <=
               static_cast<uint8_t>(ShellCount::THREE);
  }

  HS_COLD_MEMBER HyperLattice() : Choreography(W, H, {.strobe = true}) {}

  void init() override {
    configure_presets(PRESET_IDS.size());
    register_animated_param("Dimension", &params.dimension, 0.0f, 1.0f);
    register_animated_param("Sphere Radius", &params.sphere_radius, 0.0f, 1.5f);
    register_animated_param("Wire Radius", &params.wire_radius, 0.015f, 0.18f);
    register_animated_param("Softness", &params.softness, 0.002f, 0.08f);
    register_animated_param("Far Cells", &params.far_cells, 2.0f, 16.0f);
    register_animated_param("AA Strength", &params.aa_strength, 0.0f, 2.0f);
    register_animated_param("Speed", &params.speed, 0.0f, 0.05f);
    register_animated_param("3D Spin", &params.spin_3d, 0.0f, 0.015f);
    register_animated_param("4D Spin", &params.spin_4d, 0.0f, 0.015f);
    register_animated_param("Reflection", &params.reflection,
                            REFLECTION_OPTIONS, REFLECTION_EXPORT_OPTIONS,
                            std::size(REFLECTION_OPTIONS));
    register_animated_param("Color", &params.color, COLOR_OPTIONS,
                            COLOR_EXPORT_OPTIONS, std::size(COLOR_OPTIONS));
    register_animated_param("Shells", &params.shells, SHELL_OPTIONS,
                            SHELL_EXPORT_OPTIONS, std::size(SHELL_OPTIONS));
    const GenerativePalette palette{EffectPaletteRecipes::raymarch()};
    baked_palette.bake(persistent_arena, palette);
  }

  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(hl_timeline_step);
      timeline.step(canvas);
    }
    begin_automatic_transition();
    advance_state();
    const HyperLatticeDetail::FrameState context{
        params, origin, rotation_phase,
        HyperLatticeDetail::pixel_half_angle<W, H>(), &baked_palette};
    const auto frame = HyperLatticeDetail::RenderPipeline::prepare(context);
    {
      HS_PROFILE(hl_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas, [&frame](const Vector &view) {
        return HyperLatticeDetail::RenderPipeline::evaluate(view, frame.ctx,
                                                            frame.prepared);
      });
    }
  }

private:
  using Choreography::begin_automatic_transition;
  using Choreography::configure_presets;
  using Choreography::params;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;

  void blend_params(float progress) {
    params.lerp(transition.from, transition.to, progress);
  }

  void advance_state() {
    static constexpr float VELOCITY[HyperLatticeDetail::DIMENSIONS] = {
        0.4815434f, 0.2993373f, 0.4034555f, 0.7223151f};
    static constexpr float RATE[6] = {1.0f, 0.731f, 0.517f,
                                      1.0f, 0.707f, 0.419f};
    for (int axis = 0; axis < HyperLatticeDetail::DIMENSIONS; ++axis)
      origin[axis] = wrap_t(origin[axis] + params.speed * VELOCITY[axis]);
    for (int plane = 0; plane < 3; ++plane)
      rotation_phase[plane] =
          wrap(rotation_phase[plane] + params.spin_3d * RATE[plane], TWO_PI_F);
    for (int plane = 3; plane < 6; ++plane)
      rotation_phase[plane] =
          wrap(rotation_phase[plane] + params.spin_4d * RATE[plane], TWO_PI_F);
  }

  static constexpr const char *REFLECTION_OPTIONS[] = {"Chrome", "Radial"};
  static constexpr const char *REFLECTION_EXPORT_OPTIONS[] = {
      "ReflectionMode::CHROME", "ReflectionMode::RADIAL"};
  static constexpr const char *COLOR_OPTIONS[] = {"Depth", "Axis"};
  static constexpr const char *COLOR_EXPORT_OPTIONS[] = {"ColorMode::DEPTH",
                                                         "ColorMode::AXIS"};
  static constexpr const char *SHELL_OPTIONS[] = {"1", "2", "3"};
  static constexpr const char *SHELL_EXPORT_OPTIONS[] = {
      "ShellCount::ONE", "ShellCount::TWO", "ShellCount::THREE"};

  HyperLatticeDetail::Vec4 origin{{0.17f, 0.31f, 0.43f, 0.59f}};
  std::array<float, 6> rotation_phase{};
  BakedPalette baked_palette;

  friend struct hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;

  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HyperLattice persistent footprint exceeds the default "
                "partition");
};

#include "core/control/registry.h"
REGISTER_EFFECT(HyperLattice)
