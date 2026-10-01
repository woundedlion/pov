/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file HyperLattice.h
 * @brief Analytic flight through periodic lattices and shells in 3D and 4D.
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <tuple>

#include "core/color/effect_palette_recipes.h"
#include "core/render/sdf/lattice.h"
#include "core/render/ray/shade.h"
#include "core/render/pullback/ray.h"
#include "core/control/choreography.h"
#include "core/engine/engine.h"
#include "core/math/4dmath.h"
#include "core/render/pullback.h"
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
#include "core/render/sdf/lattice_trace.h"
#endif

namespace hs_test {
namespace hyper_lattice_tests {
struct HyperLatticeWhiteBox;
} // namespace hyper_lattice_tests
} // namespace hs_test

namespace HyperLatticeDetail {
constexpr int DIMENSIONS = math::VEC4_DIMENSIONS;
using LatticeMode = SDF::Lattice::Domain;
using ShellCount = SDF::Lattice::ShellCount;
enum class Pattern : uint8_t {
  CUBIC_WIRE,
  OCTET,
  DIAMOND,
  HEXAGONAL,
  RHOMBIC,
  AFFINE_CUBIC,
  SHELLS
};
enum class ConfigurationId : uint8_t {
  CUBIC_3D,
  CUBIC_4D,
  OCTET_3D,
  OCTET_4D,
  DIAMOND_3D,
  HEXAGONAL_3D,
  RHOMBIC_3D,
  AFFINE_3D,
  AFFINE_4D,
  SHELLS_3D,
  SHELLS_4D,
  INVALID
};
struct Params {
  LatticeMode mode = LatticeMode::THREE_D;
  Pattern pattern = Pattern::CUBIC_WIRE;
  float sphere_radius = 1.0f;
  float cell_size = 1.0f;
  float wire_radius = 0.055f;
  float softness = 0.012f;
  float near_fade = 0.5f;
  float far_distance = 7.0f;
  float aa_strength = 1.0f;
  float speed = 0.018f;
  float spin_3d = 0.0024f;
  float spin_4d = 0.0f;
  ShellCount shells = ShellCount::TWO;
  float shear = .55f;
  float stretch = 1.4f;
  float shell_radius = .30f;

  /**
   * @brief Interpolates every continuous field; the pattern, view and shell
   *        count switch together at the midpoint.
   */
  void lerp(const Params &start, const Params &target, float amount) {
    const Params &NEAREST = amount < 0.5f ? start : target;
    mode = NEAREST.mode;
    pattern = NEAREST.pattern;
    sphere_radius = hs::lerp(start.sphere_radius, target.sphere_radius, amount);
    cell_size = hs::lerp(start.cell_size, target.cell_size, amount);
    wire_radius = hs::lerp(start.wire_radius, target.wire_radius, amount);
    softness = hs::lerp(start.softness, target.softness, amount);
    near_fade = hs::lerp(start.near_fade, target.near_fade, amount);
    far_distance = hs::lerp(start.far_distance, target.far_distance, amount);
    aa_strength = hs::lerp(start.aa_strength, target.aa_strength, amount);
    speed = hs::lerp(start.speed, target.speed, amount);
    spin_3d = hs::lerp(start.spin_3d, target.spin_3d, amount);
    spin_4d = hs::lerp(start.spin_4d, target.spin_4d, amount);
    shells = NEAREST.shells;
    shear = hs::lerp(start.shear, target.shear, amount);
    stretch = hs::lerp(start.stretch, target.stretch, amount);
    shell_radius = hs::lerp(start.shell_radius, target.shell_radius, amount);
  }

  /**
   * @brief Compile-time field-set pin for lerp() and valid_params(); never
   *        called.
   * @details The binding names every member. Interior padding after pattern
   *          and shells can absorb small fields without changing sizeof().
   */
  static void pin_field_set(const Params &p) {
    const auto &[mode, pattern, sphere_radius, cell_size, wire_radius, softness,
                 near_fade, far_distance, aa_strength, speed, spin_3d, spin_4d,
                 shells, shear, stretch, shell_radius] = p;
    (void)mode, (void)pattern, (void)sphere_radius, (void)cell_size,
        (void)wire_radius, (void)softness, (void)near_fade, (void)far_distance,
        (void)aa_strength, (void)speed, (void)spin_3d, (void)spin_4d,
        (void)shells, (void)shear, (void)stretch, (void)shell_radius;
  }
};

// Width pin; pin_field_set() is what catches an added or removed field. Every
// enum has a fixed uint8_t base, so the size holds under ARM -fshort-enums.
static_assert(sizeof(Params) == 60, "HyperLattice::Params width changed");

using SDF::Lattice::CrossingList;

struct FrameState {
  Params params;
  math::Vec4 origin;
  std::array<float, 6> rotation_phase;
  float pixel_half_angle;
  const BakedPalette *depth_palette;
  float gain = 1.0f; /**< Brightness scale of the whole frame. */
  CrossingList *crossings =
      nullptr; /**< Arena scratch the cubic trace sorts. */
};

struct Binding {
  using FrameState = HyperLatticeDetail::FrameState;
  using Instrumentation = Pullback::NoInstrumentation;
};

using PreparedTrace = SDF::Lattice::PreparedShading;
using SDF::Lattice::composite_crossings;

template <int W, int H> constexpr float pixel_half_angle() {
  return .5f * math::coarse_pixel_pitch<W, H>();
}
HS_FLASH_INLINE inline math::Mat4 view_embedding(const FrameState &frame) {
  math::Mat4 embedding = math::Mat4::identity();
  math::rotate_plane(embedding, 0, 1, frame.rotation_phase[0]);
  math::rotate_plane(embedding, 0, 2, frame.rotation_phase[1]);
  math::rotate_plane(embedding, 1, 2, frame.rotation_phase[2]);
  if (frame.params.mode == LatticeMode::FOUR_D_SLICE) {
    math::rotate_plane(embedding, 0, 3, frame.rotation_phase[3]);
    math::rotate_plane(embedding, 1, 3, frame.rotation_phase[4]);
    math::rotate_plane(embedding, 2, 3, frame.rotation_phase[5]);
  }
  return embedding;
}
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
namespace Experiment = SDF::LatticeTrace;

static_assert(static_cast<uint8_t>(Pattern::OCTET) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::OCTET));
static_assert(static_cast<uint8_t>(Pattern::DIAMOND) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::DIAMOND));
static_assert(static_cast<uint8_t>(Pattern::HEXAGONAL) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::HEXAGONAL));
static_assert(static_cast<uint8_t>(Pattern::RHOMBIC) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::RHOMBIC));
static_assert(static_cast<uint8_t>(Pattern::AFFINE_CUBIC) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::AFFINE_CUBIC));
static_assert(static_cast<uint8_t>(Pattern::SHELLS) - 1 ==
              static_cast<uint8_t>(Experiment::Geometry::SHELLS));

HS_FLASH_INLINE inline Experiment::Settings
experimental_settings(const FrameState &frame, math::Vec4 center) {
  const auto &p = frame.params;
  if (p.mode == LatticeMode::THREE_D)
    center[3] = 0;
  return {
      p.mode == LatticeMode::FOUR_D_SLICE ? Raycast::SamplingDomain::SLICE_4D
                                          : Raycast::SamplingDomain::SPATIAL_3D,
      p.cell_size,
      p.wire_radius * p.cell_size,
      p.sphere_radius,
      p.far_distance,
      p.near_fade,
      p.aa_strength,
      center,
      view_embedding(frame),
      frame.pixel_half_angle,
      frame.depth_palette,
      static_cast<Experiment::Geometry>(static_cast<uint8_t>(p.pattern) - 1),
      p.shear,
      p.stretch,
      p.shell_radius,
      frame.gain};
}
#endif
inline PreparedTrace prepare_trace(const FrameState &frame) {
  HS_CHECK(frame.crossings, "HyperLattice: frame has no crossing list");
  const auto embedding = view_embedding(frame);
  const auto &p = frame.params;
  const SDF::Lattice::Settings settings{
      p.mode,     p.sphere_radius, p.cell_size, p.wire_radius,
      p.softness, p.aa_strength,   p.shells};
  const float near_scale = p.cell_size * (1 + p.sphere_radius);
  return {SDF::Lattice::prepare(settings, frame.origin, embedding,
                                p.far_distance, frame.pixel_half_angle),
          {1.0f / p.far_distance, 1.5f * p.wire_radius * near_scale,
           1.0f / (p.near_fade * near_scale), frame.depth_palette, frame.gain},
          frame.crossings};
}

template <bool SLICE_4D = false, uint8_t SHELLS = 0> struct Renderer {
  static PreparedTrace prepare(const FrameState &frame) {
    return prepare_trace(frame);
  }
  __attribute__((always_inline)) static Color4
  shade(const math::Vector &normal, const FrameState &,
        const PreparedTrace &prepared) {
    return composite_crossings<SLICE_4D, SHELLS>(normal, prepared).finish();
  }
  /** @brief shade() premultiplied by its alpha, without the round trip. */
  __attribute__((always_inline)) static Pixel
  shade_premultiplied(const math::Vector &normal,
                      const PreparedTrace &prepared) {
    return composite_crossings<SLICE_4D, SHELLS>(normal, prepared)
        .premultiplied();
  }
};
using RenderPipeline =
    Pullback::Pipeline<Binding, Pullback::RayStage<Renderer<>>>;
template <uint8_t SHELL_COUNT>
using SpecializedRenderPipeline =
    Pullback::Pipeline<Binding,
                       Pullback::RayStage<Renderer<true, SHELL_COUNT>>>;
} // namespace HyperLatticeDetail

/**
 * @brief Flights through periodic lattices, shells, and their 4D slices.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H>
class HyperLattice : public ChoreographedEffect<HyperLattice<W, H>,
                                                HyperLatticeDetail::Params> {
  using Choreography =
      ChoreographedEffect<HyperLattice<W, H>, HyperLatticeDetail::Params>;
  friend Choreography;

public:
  static constexpr const char *EFFECT_ID = "HyperLattice";

  using Params = HyperLatticeDetail::Params;
  using LatticeMode = HyperLatticeDetail::LatticeMode;
  using ShellCount = HyperLatticeDetail::ShellCount;
  using Pattern = HyperLatticeDetail::Pattern;
  using ConfigurationId = HyperLatticeDetail::ConfigurationId;

  static constexpr auto PRESET_IDS = std::to_array<std::string_view>({
      "cubic-flight",
      "cubic-wide-flight",
      "hypercube-flight",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "experimental-octet-flight",
      "experimental-octet-wide-flight",
      "experimental-octet-4d-flight",
      "experimental-shell-flight",
      "experimental-shell-close-flight",
      "experimental-shell-4d-flight",
#endif
  });
  static constexpr size_t CUBIC_PRESET_INDEX = 0;
  static constexpr size_t WIDE_PRESET_INDEX = 1;
  static constexpr size_t HYPERCUBE_PRESET_INDEX = 2;
  static constexpr size_t OCTET_PRESET_INDEX = 3;
  static constexpr size_t OCTET_WIDE_PRESET_INDEX = 4;
  static constexpr size_t OCTET_4D_PRESET_INDEX = 5;
  static constexpr size_t SHELL_PRESET_INDEX = 6;
  static constexpr size_t SHELL_CLOSE_PRESET_INDEX = 7;
  static constexpr size_t SHELL_4D_PRESET_INDEX = 8;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 14;

  /**
   * @brief The preset at @p index and how it departs.
   * @details A preset morphs through its parameters into the next preset of
   * its pattern and view; any other departure fades through black, since the
   * geometry cannot interpolate across a pattern or view.
   */
  HS_COLD_MEMBER static constexpr PresetEntry<Params> preset(size_t index) {
    constexpr Segue::Preset::Lerp MORPH{240, math::ease_in_out_sin,
                                        /*pausable=*/true};
    constexpr Segue::Preset::Fade FADE{240};
    const bool MORPHS = index == CUBIC_PRESET_INDEX ||
                        index == OCTET_PRESET_INDEX ||
                        index == SHELL_PRESET_INDEX;
    Params value;
    switch (index) {
    case CUBIC_PRESET_INDEX:
    case WIDE_PRESET_INDEX:
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = index == CUBIC_PRESET_INDEX ? 1.0f : 0.0f;
      value.cell_size = index == CUBIC_PRESET_INDEX ? 1.0f : 2.38525f;
      value.wire_radius = 0.055f;
      value.softness = 0.08f;
      value.near_fade = index == CUBIC_PRESET_INDEX ? 0.5f : 2.0f;
      value.far_distance = index == CUBIC_PRESET_INDEX ? 4.198f : 11.66f;
      value.aa_strength = 1.0f;
      value.speed = 0.05f;
      value.spin_3d = 0.015f;
      value.spin_4d = 0.0f;
      value.shells = ShellCount::TWO;
      break;
    case HYPERCUBE_PRESET_INDEX:
      value.mode = LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0.0f;
      value.cell_size = 1.0f;
      value.wire_radius = 0.03546f;
      value.softness = 0.029612f;
      value.far_distance = 8.0f;
      value.aa_strength = 1.0f;
      value.speed = 0.03f;
      value.spin_3d = 0.01089f;
      value.spin_4d = 0.015f;
      value.shells = ShellCount::TWO;
      break;
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    case OCTET_PRESET_INDEX:
      value.pattern = Pattern::OCTET;
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = 0;
      value.cell_size = 1.74175f;
      value.wire_radius = .015f;
      value.softness = .012f;
      value.far_distance = 4.5f;
      value.near_fade = 2.0f;
      value.aa_strength = 2.0f;
      value.speed = .078f;
      value.spin_3d = .008265f;
      break;
    case OCTET_WIDE_PRESET_INDEX:
      value.pattern = Pattern::OCTET;
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = 1.0f;
      value.cell_size = 3.82825f;
      value.wire_radius = .015f;
      value.near_fade = 2.0f;
      value.far_distance = 10.666f;
      value.aa_strength = 2.0f;
      value.speed = .12750001f;
      value.spin_3d = .010155f;
      break;
    case OCTET_4D_PRESET_INDEX:
      value.pattern = Pattern::OCTET;
      value.mode = LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0;
      value.cell_size = 3.448f;
      value.wire_radius = .02919f;
      value.near_fade = 2.0f;
      value.far_distance = 16.0f;
      value.aa_strength = 2.0f;
      value.speed = .008f;
      value.spin_3d = .01146f;
      value.spin_4d = .010875f;
      break;
    case SHELL_PRESET_INDEX:
    case SHELL_CLOSE_PRESET_INDEX: {
      const bool FLIGHT = index == SHELL_PRESET_INDEX;
      value.pattern = Pattern::SHELLS;
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = 0;
      value.cell_size = FLIGHT ? .78625f : .4645f;
      value.near_fade = FLIGHT ? 2.0f : .6f;
      value.far_distance = FLIGHT ? 10.736f : 5.836f;
      value.aa_strength = 2.0f;
      value.speed = .025f;
      value.spin_3d = FLIGHT ? .015f : .003f;
      value.stretch = 1.0f;
      value.shell_radius = FLIGHT ? .1f : .15f;
      break;
    }
    case SHELL_4D_PRESET_INDEX:
      value.pattern = Pattern::SHELLS;
      value.mode = LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0;
      value.cell_size = 1.0f;
      value.near_fade = 2.0f;
      value.far_distance = 16.0f;
      value.aa_strength = 2.0f;
      value.speed = .05f;
      value.spin_3d = .005f;
      value.spin_4d = .005f;
      value.stretch = 1.0f;
      value.shell_radius = .15f;
      break;
#endif
    default:
      break;
    }
    return {value, MORPHS ? Segue::Preset::Departure{MORPH}
                          : Segue::Preset::Departure{FADE}};
  }

  static constexpr Params pattern_defaults(Pattern pattern, LatticeMode mode) {
    const bool SLICE = mode == LatticeMode::FOUR_D_SLICE;
    if (pattern == Pattern::CUBIC_WIRE)
      return preset(SLICE ? HYPERCUBE_PRESET_INDEX : CUBIC_PRESET_INDEX).params;
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    if (pattern == Pattern::OCTET)
      return preset(SLICE ? OCTET_4D_PRESET_INDEX : OCTET_PRESET_INDEX).params;
    if (pattern == Pattern::SHELLS)
      return preset(SLICE ? SHELL_4D_PRESET_INDEX : SHELL_PRESET_INDEX).params;
#endif
    Params value;
    value.pattern = pattern;
    value.mode = mode;
    value.sphere_radius = 0;
    value.cell_size = pattern == Pattern::HEXAGONAL ? 1.4f : 2.0f;
    value.wire_radius = .025f;
    value.far_distance = 6.0f;
    value.near_fade = .6f;
    value.speed = .025f;
    value.spin_3d = .003f;
    value.spin_4d = mode == LatticeMode::FOUR_D_SLICE ? .004f : 0;
    value.stretch = 1.4f;
    return value;
  }

  struct Configuration {
    Pattern pattern;
    LatticeMode domain;
    uint8_t max_candidates;
    uint8_t max_layers;
  };
  static constexpr auto CONFIGURATIONS = std::to_array<Configuration>({
      {Pattern::CUBIC_WIRE, LatticeMode::THREE_D, 9, 9},
      {Pattern::CUBIC_WIRE, LatticeMode::FOUR_D_SLICE, 12, 12},
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      {Pattern::OCTET, LatticeMode::THREE_D, 64, 32},
      {Pattern::OCTET, LatticeMode::FOUR_D_SLICE, 64, 32},
      {Pattern::DIAMOND, LatticeMode::THREE_D, 192, 32},
      {Pattern::HEXAGONAL, LatticeMode::THREE_D, 192, 32},
      {Pattern::RHOMBIC, LatticeMode::THREE_D, 192, 32},
      {Pattern::AFFINE_CUBIC, LatticeMode::THREE_D, 96, 32},
      {Pattern::AFFINE_CUBIC, LatticeMode::FOUR_D_SLICE, 96, 32},
      {Pattern::SHELLS, LatticeMode::THREE_D, 64, 32},
      {Pattern::SHELLS, LatticeMode::FOUR_D_SLICE, 64, 32},
#endif
  });

  static constexpr ConfigurationId configuration_id(const Params &value) {
    for (size_t i = 0; i < CONFIGURATIONS.size(); ++i)
      if (CONFIGURATIONS[i].pattern == value.pattern &&
          CONFIGURATIONS[i].domain == value.mode)
        return static_cast<ConfigurationId>(i);
    return ConfigurationId::INVALID;
  }

  static constexpr bool supported_combination(const Params &value) {
    return configuration_id(value) != ConfigurationId::INVALID;
  }

  static constexpr float SHEAR_MIN = -1.0f;
  static constexpr float SHEAR_MAX = 1.0f;
  static constexpr float STRETCH_MIN = 1.0f;
  static constexpr float STRETCH_MAX = 1.5f;
  static constexpr float SHELL_RADIUS_MIN = .10f;
  static constexpr float SHELL_RADIUS_MAX = .32f;

  static constexpr bool valid_params(const Params &value) {
    return supported_combination(value) &&
           value.sphere_radius >= SPHERE_RADIUS_MIN &&
           value.sphere_radius <= SPHERE_RADIUS_MAX &&
           value.cell_size >= CELL_SIZE_MIN &&
           value.cell_size <= CELL_SIZE_MAX &&
           value.wire_radius >= WIRE_RADIUS_MIN &&
           value.wire_radius <= WIRE_RADIUS_MAX &&
           value.softness >= SOFTNESS_MIN && value.softness <= SOFTNESS_MAX &&
           value.near_fade >= NEAR_FADE_MIN &&
           value.near_fade <= NEAR_FADE_MAX &&
           value.far_distance >= FAR_DISTANCE_MIN &&
           value.far_distance <= FAR_DISTANCE_MAX &&
           value.aa_strength >= AA_STRENGTH_MIN &&
           value.aa_strength <= AA_STRENGTH_MAX && value.speed >= SPEED_MIN &&
           value.speed <= SPEED_MAX && value.spin_3d >= SPIN_3D_MIN &&
           value.spin_3d <= SPIN_3D_MAX && value.spin_4d >= SPIN_4D_MIN &&
           value.spin_4d <= SPIN_4D_MAX &&
           static_cast<uint8_t>(value.shells) <=
               static_cast<uint8_t>(ShellCount::THREE) &&
           value.shear >= SHEAR_MIN && value.shear <= SHEAR_MAX &&
           value.stretch >= STRETCH_MIN && value.stretch <= STRETCH_MAX &&
           value.shell_radius >= SHELL_RADIUS_MIN &&
           value.shell_radius <= SHELL_RADIUS_MAX;
  }

  /**
   * @brief Whether a parameter set matches a specialized slice pipeline.
   * @param value Parameters to test.
   * @return true when the parameters match the specialized trace assumptions.
   * @details The shape is the hypercube preset's; the assert in draw_frame()
   *          ties the two, so retuning that preset cannot leave the gate
   *          behind.
   */
  static constexpr bool uses_specialized_slice(const Params &value) {
    return value.pattern == Pattern::CUBIC_WIRE &&
           value.mode == LatticeMode::FOUR_D_SLICE &&
           (value.shells == ShellCount::TWO ||
            value.shells == ShellCount::THREE);
  }

  HS_COLD_MEMBER HyperLattice() : Choreography(W, H, {.strobe = true}) {}

  HS_COLD_MEMBER void init() override {
    begin_choreography();
    register_animated_param("Pattern", &params.pattern, PATTERN_OPTIONS,
                            PATTERN_EXPORT_OPTIONS, std::size(PATTERN_OPTIONS));
    register_animated_param("View", &params.mode, VIEW_OPTIONS,
                            VIEW_EXPORT_OPTIONS, std::size(VIEW_OPTIONS));
    register_animated_param("Sphere Radius", &params.sphere_radius,
                            SPHERE_RADIUS_MIN, SPHERE_RADIUS_MAX);
    register_animated_param("Cell Size", &params.cell_size, CELL_SIZE_MIN,
                            CELL_SIZE_MAX);
    register_animated_param("Wire Radius", &params.wire_radius, WIRE_RADIUS_MIN,
                            WIRE_RADIUS_MAX);
    register_animated_param("Softness", &params.softness, SOFTNESS_MIN,
                            SOFTNESS_MAX);
    register_animated_param("Near Fade", &params.near_fade, NEAR_FADE_MIN,
                            NEAR_FADE_MAX);
    register_animated_param("Far Distance", &params.far_distance,
                            FAR_DISTANCE_MIN, FAR_DISTANCE_MAX);
    register_animated_param("AA Strength", &params.aa_strength, AA_STRENGTH_MIN,
                            AA_STRENGTH_MAX);
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("3D Spin", &params.spin_3d, SPIN_3D_MIN,
                            SPIN_3D_MAX);
    register_animated_param("4D Spin", &params.spin_4d, SPIN_4D_MIN,
                            SPIN_4D_MAX);
    register_animated_param("Lattice Planes", &params.shells, SHELL_OPTIONS,
                            SHELL_EXPORT_OPTIONS, std::size(SHELL_OPTIONS));
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    register_animated_param("Shear", &params.shear, SHEAR_MIN, SHEAR_MAX);
    register_animated_param("Stretch", &params.stretch, STRETCH_MIN,
                            STRETCH_MAX);
    register_animated_param("Shell Radius", &params.shell_radius,
                            SHELL_RADIUS_MIN, SHELL_RADIUS_MAX);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    this->register_readonly_param("Unfinished Rays", &unfinished_rays, 0,
                                  (W + 2) * (H + 2));
#endif
#endif
    refresh_configuration_schema();
    depth_palette.init_generated(persistent_arena, next_depth_palette, nullptr,
                                 0, PALETTE_FADE_FRAMES, math::ease_in_out_sin);
    crossing_list =
        persistent_arena.allocate_n<HyperLatticeDetail::CrossingList>(1);
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    crossing_storage =
        persistent_arena.allocate_n<SDF::OctetTrace::CrossingStorage>(1);
    shell_layers = persistent_arena.allocate_n<SDF::ShellLayerStorage>(1);
    cellular_hits =
        persistent_arena.allocate_n<SDF::CellularWire::HitStorage>(1);
#endif
  }

  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(hl_timeline_step);
      timeline.step(canvas);
    }
    step_choreography();
    advance_state();
    depth_palette.step();
    const HyperLatticeDetail::FrameState context{
        params,
        origin,
        rotation_phase,
        HyperLatticeDetail::pixel_half_angle<W, H>(),
        &depth_palette.palette(),
        preset_gain,
        crossing_list};
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    unfinished_rays = 0;
    if (params.pattern != Pattern::CUBIC_WIRE) {
      HS_PROFILE(hl_shader_draw);
      draw_experimental(canvas, context);
      return;
    }
#endif
    const auto prepared = HyperLatticeDetail::prepare_trace(context);
    {
      HS_PROFILE(hl_shader_draw);
      static_assert(
          uses_specialized_slice(preset(HYPERCUBE_PRESET_INDEX).params),
          "the hypercube preset no longer selects the specialized slice trace");
      if (uses_specialized_slice(params) && params.shells == ShellCount::TWO) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&prepared](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::Renderer<true, 2>::shade_premultiplied(
                  view, prepared);
            });
      } else if (uses_specialized_slice(params)) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&prepared](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::Renderer<true>::shade_premultiplied(
                  view, prepared);
            });
      } else if (params.shells == ShellCount::TWO) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&prepared](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::Renderer<
                  false, 2>::shade_premultiplied(view, prepared);
            });
      } else {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&prepared](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::Renderer<>::shade_premultiplied(
                  view, prepared);
            });
      }
    }
  }

#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    HS_CHECK(index < PRESET_IDS.size(),
             "HyperLattice profile preset index out of range");
    HS_CHECK(this->selectPreset(index),
             "HyperLattice profile preset selection failed");
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESET_IDS.size()));
  }
#endif

private:
  using Choreography::begin_choreography;
  using Choreography::params;
  using Choreography::step_choreography;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;

#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_COLD_MEMBER bool parameter_write_admitted(const ParamDef &parameter,
                                               float value) override {
    Params candidate = params;
    if (parameter.target == &params.pattern) {
      candidate.pattern = static_cast<Pattern>(value);
      if (!supported_combination(candidate))
        candidate.mode = LatticeMode::THREE_D;
    } else if (parameter.target == &params.mode)
      candidate.mode = static_cast<LatticeMode>(value);
    else
      return true;
    return supported_combination(candidate);
  }
#endif

  HS_COLD_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    if (!supported_combination(params))
      params.mode = LatticeMode::THREE_D;
    const auto configuration = configuration_id(params);
    if (selected_configuration != configuration) {
      const float near_fade = params.near_fade;
      params = pattern_defaults(params.pattern, params.mode);
      params.near_fade = near_fade;
      selected_configuration = configuration;
      refresh_configuration_schema();
    }
  }
  void adopt_params(const Params &target) {
    params = target;
    if (configuration_id(params) != selected_configuration) {
      selected_configuration = configuration_id(params);
      refresh_configuration_schema();
    }
  }
  HS_COLD_MEMBER void blend_params(float progress) {
    params.lerp(transition.from, transition.to, progress);
    if (configuration_id(params) != selected_configuration) {
      selected_configuration = configuration_id(params);
      refresh_configuration_schema();
    }
  }

  HS_COLD_MEMBER void refresh_configuration_schema() {
    if (auto *parameter = this->getParameters().find("4D Spin")) {
      const bool readonly = params.mode == LatticeMode::THREE_D;
      if (parameter->readonly != readonly) {
        this->mark_readonly("4D Spin", readonly);
      }
    }
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    const bool CUBIC = params.pattern == Pattern::CUBIC_WIRE;
    this->mark_readonly("Lattice Planes", !CUBIC);
    this->mark_readonly("Softness", !CUBIC);
    const bool AFFINE = params.pattern == Pattern::AFFINE_CUBIC;
    const bool SHELLS = params.pattern == Pattern::SHELLS;
    this->mark_readonly("Shear", !AFFINE);
    this->mark_readonly("Stretch", !AFFINE);
    this->mark_readonly("Shell Radius", !SHELLS);
    this->mark_readonly("Wire Radius", SHELLS);
    this->mark_readonly("View", params.pattern == Pattern::DIAMOND ||
                                    params.pattern == Pattern::HEXAGONAL ||
                                    params.pattern == Pattern::RHOMBIC);
#endif
  }

  HS_FLASH_MEMBER void advance_state() {
    static constexpr float VELOCITY[HyperLatticeDetail::DIMENSIONS] = {
        0.4815434f, 0.2993373f, 0.4034555f, 0.7223151f};
    static constexpr float RATE[6] = {1.0f, 0.731f, 0.517f,
                                      1.0f, 0.707f, 0.419f};
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    if (params.pattern != Pattern::CUBIC_WIRE) {
      for (int axis = 0; axis < HyperLatticeDetail::DIMENSIONS; ++axis)
        experimental_center[axis] += params.speed * VELOCITY[axis];
      if (params.pattern == Pattern::AFFINE_CUBIC) {
        experimental_center =
            SDF::AffineLattice{params.cell_size, params.shear, params.stretch}
                .wrap(experimental_center);
      } else {
        for (int axis = 0; axis < HyperLatticeDetail::DIMENSIONS; ++axis) {
          float period = params.cell_size;
          if (params.pattern == Pattern::OCTET)
            period *= 1.4142135623730951f;
          else if (params.pattern == Pattern::HEXAGONAL && axis < 2)
            period *= axis == 0 ? 3.0f : 1.7320508075688772f;
          experimental_center[axis] =
              math::wrap(experimental_center[axis], period);
        }
      }
    } else
#endif
      for (int axis = 0; axis < HyperLatticeDetail::DIMENSIONS; ++axis)
        origin[axis] =
            math::wrap_t(origin[axis] + params.speed * VELOCITY[axis]);
    for (int plane = 0; plane < 3; ++plane)
      rotation_phase[plane] = math::wrap(
          rotation_phase[plane] + params.spin_3d * RATE[plane], math::TWO_PI_F);
    const float SPIN_4D_STEP =
        (params.mode == LatticeMode::FOUR_D_SLICE ? params.spin_4d : 0.0f);
    for (int plane = 3; plane < 6; ++plane)
      rotation_phase[plane] = math::wrap(
          rotation_phase[plane] + SPIN_4D_STEP * RATE[plane], math::TWO_PI_F);
  }

#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
  HS_FLASH_MEMBER void
  draw_experimental(Canvas &canvas,
                    const HyperLatticeDetail::FrameState &context) {
    using namespace HyperLatticeDetail::Experiment;
    auto settings =
        HyperLatticeDetail::experimental_settings(context, experimental_center);
    settings.crossings = crossing_storage;
    settings.shell_layers = shell_layers;
    settings.cellular_hits = cellular_hits;
    auto prepared = prepare(settings);
    const auto &configuration =
        CONFIGURATIONS[static_cast<size_t>(configuration_id(params))];
    prepared.limits.max_candidates = configuration.max_candidates;
    prepared.limits.max_layers = configuration.max_layers;
    if (prepared.valid && prepared.geometry == Geometry::SHELLS &&
        prepared.periodic_shells.single_owner) {
      Scan::Shader::draw_cached<W, H, 1>(
          canvas,
          [&prepared, this](const math::Vector &view) HS_HOT_FLASH_MEMBER {
            const auto result = SDF::trace_periodic_shells_3d(
                prepared.periodic_shells, prepared.camera, view,
                prepared.limits, prepared.appearance);
            if (result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED)
              unfinished_rays += 1;
            return result.color;
          });
      return;
    }
    const auto draw_march = [&]<int DIMENSIONS>() {
      Scan::Shader::draw_cached<W, H, 1>(
          canvas,
          [&prepared, this](const math::Vector &view) HS_HOT_FLASH_MEMBER {
            const auto result = SDF::trace_periodic_shells_march<DIMENSIONS>(
                prepared.periodic_shells, prepared.camera, view,
                prepared.limits, prepared.appearance, *prepared.shell_layers);
            if (result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED)
              unfinished_rays += 1;
            return result.color;
          });
    };
    if (prepared.valid && prepared.geometry == Geometry::SHELLS &&
        prepared.periodic_shells.march) {
      if (params.mode == LatticeMode::FOUR_D_SLICE)
        draw_march.template operator()<4>();
      else
        draw_march.template operator()<3>();
      return;
    }
    const auto draw_octet = [&]<bool SLICE_4D>() {
      Scan::Shader::draw_cached<W, H, 1>(
          canvas,
          [&prepared, this](const math::Vector &view) HS_HOT_FLASH_MEMBER {
            const auto result = shade_octet<SLICE_4D>(view, prepared);
            if (result.status == Raycast::TraceStatus::BUDGET_EXHAUSTED)
              unfinished_rays += 1;
            return result.color;
          });
    };
    if (prepared.valid && prepared.geometry == Geometry::OCTET) {
      if (params.mode == LatticeMode::FOUR_D_SLICE)
        draw_octet.template operator()<true>();
      else
        draw_octet.template operator()<false>();
      return;
    }
    using Shade = Sample (*)(const math::Vector &, const Prepared &);
    const Shade shade_ray =
        params.mode == LatticeMode::FOUR_D_SLICE ? &shade<true> : &shade<false>;
    Scan::Shader::draw_cached<W, H, 1>(
        canvas, [&prepared, shade_ray, this](const math::Vector &view)
                    HS_HOT_FLASH_MEMBER {
                      const auto result = shade_ray(view, prepared);
                      const auto STATUS = result.status;
                      if (STATUS != Raycast::TraceStatus::SURFACE &&
                          STATUS != Raycast::TraceStatus::RANGE_COMPLETE &&
                          STATUS != Raycast::TraceStatus::SATURATED)
                        unfinished_rays += 1;
                      return result.color;
                    });
  }

  math::Vec4 experimental_center{{.255f, .465f, .645f, .375f}};
  float unfinished_rays = 0;
#endif

  static void next_depth_palette(void *, uint32_t sequence,
                                 GenerativePalette &out) {
    static constexpr std::array<float, 3> HUE_OFFSETS{0.0f, -8.0f / 360.0f,
                                                      8.0f / 360.0f};
    out = GenerativePalette{EffectPaletteRecipes::hyper_lattice(
        HUE_OFFSETS[sequence % HUE_OFFSETS.size()])};
  }

  static constexpr float SPHERE_RADIUS_MIN = 0.0f, SPHERE_RADIUS_MAX = 2.0f;
  static constexpr float CELL_SIZE_MIN = 0.25f, CELL_SIZE_MAX = 10.0f;
  static constexpr float WIRE_RADIUS_MIN = 0.015f, WIRE_RADIUS_MAX = 0.18f;
  static constexpr float SOFTNESS_MIN = 0.002f, SOFTNESS_MAX = 0.08f;
  static constexpr float NEAR_FADE_MIN = 0.01f, NEAR_FADE_MAX = 2.0f;
  static constexpr float FAR_DISTANCE_MIN = 2.0f, FAR_DISTANCE_MAX = 16.0f;
  static constexpr float AA_STRENGTH_MIN = 0.0f, AA_STRENGTH_MAX = 2.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 0.3f;
  static constexpr float SPIN_3D_MIN = 0.0f, SPIN_3D_MAX = 0.015f;
  static constexpr float SPIN_4D_MIN = 0.0f, SPIN_4D_MAX = 0.015f;

  static constexpr int PALETTE_FADE_FRAMES = 960;

  static constexpr const char *PATTERN_OPTIONS[] = {
      "Cubic",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "Experimental / Octet Truss",
      "Experimental / Diamond",
      "Experimental / Honeycomb",
      "Experimental / Rhombic Cells",
      "Experimental / Sheared Cubic",
      "Experimental / Shells",
#endif
  };
  static constexpr const char *PATTERN_EXPORT_OPTIONS[] = {
      "Pattern::CUBIC_WIRE",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "Pattern::OCTET",      "Pattern::DIAMOND",      "Pattern::HEXAGONAL",
      "Pattern::RHOMBIC",    "Pattern::AFFINE_CUBIC", "Pattern::SHELLS",
#endif
  };
  static constexpr const char *VIEW_OPTIONS[] = {"3D perspective", "4D slice"};
  static constexpr const char *VIEW_EXPORT_OPTIONS[] = {
      "LatticeMode::THREE_D", "LatticeMode::FOUR_D_SLICE"};
  static constexpr const char *SHELL_OPTIONS[] = {"1", "2", "3"};
  static constexpr const char *SHELL_EXPORT_OPTIONS[] = {
      "ShellCount::ONE", "ShellCount::TWO", "ShellCount::THREE"};

  ConfigurationId selected_configuration = ConfigurationId::CUBIC_3D;
  float preset_gain = 1.0f;

  /** @brief Receives a fading departure's opacity as the frame's gain. */
  void set_preset_opacity(float value) { preset_gain = value; }
  math::Vec4 origin{{0.17f, 0.31f, 0.43f, 0.59f}};
  std::array<float, 6> rotation_phase{};
  PaletteCycler depth_palette;
  HyperLatticeDetail::CrossingList *crossing_list = nullptr;
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
  SDF::OctetTrace::CrossingStorage *crossing_storage = nullptr;
  SDF::ShellLayerStorage *shell_layers = nullptr;
  SDF::CellularWire::HitStorage *cellular_hits = nullptr;
#endif

  friend struct hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;

  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes() +
      sizeof(HyperLatticeDetail::CrossingList) +
      alignof(HyperLatticeDetail::CrossingList)
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      + sizeof(SDF::OctetTrace::CrossingStorage) +
      alignof(SDF::OctetTrace::CrossingStorage) +
      sizeof(SDF::ShellLayerStorage) + alignof(SDF::ShellLayerStorage) +
      sizeof(SDF::CellularWire::HitStorage) +
      alignof(SDF::CellularWire::HitStorage)
#endif
      ;
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HyperLattice persistent footprint exceeds the default "
                "partition");
};

static_assert(
    [] {
      using Effect = HyperLattice<1, 1>;
      for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
        if (!Effect::valid_params(Effect::preset(index).params))
          return false;
      return true;
    }(),
    "HyperLattice preset is outside a registered slider range");
