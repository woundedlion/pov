/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file HyperLattice.h
 * @brief Analytic flight through cubic and octet lattices in 3D and 4D.
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
#include "effects/HyperLatticeExperimental.h"
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
using ColorMode = Raycast::ColorMode;
enum class Pattern : uint8_t { CUBIC_WIRE, OCTET };
enum class ConfigurationId : uint8_t { CUBIC_3D, CUBIC_4D, OCTET_3D, OCTET_4D };
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
  ColorMode color = ColorMode::DEPTH;
  ShellCount shells = ShellCount::TWO;

  void lerp(const Params &start, const Params &target, float amount) {
    if (start.mode != target.mode || start.pattern != target.pattern) {
      *this = amount < 0.5f ? start : target;
      near_fade = hs::lerp(start.near_fade, target.near_fade, amount);
      return;
    }
    mode = start.mode;
    pattern = start.pattern;
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
    color = amount < 0.5f ? start.color : target.color;
    shells = amount < 0.5f ? start.shells : target.shells;
  }

  /**
   * @brief Compile-time field-set pin for lerp() and valid_params(); never
   *        called.
   * @details The binding names every member, so adding or removing a field is a
   *          build error here. sizeof() cannot stand in: the trailing enums
   *          leave tail padding that absorbs an added small field and leaves
   *          the size unchanged, after which the new field holds at the
   *          departing preset's value across every crossfade.
   */
  static void pin_field_set(const Params &p) {
    const auto &[mode, pattern, sphere_radius, cell_size, wire_radius, softness,
                 near_fade, far_distance, aa_strength, speed, spin_3d, spin_4d,
                 color, shells] = p;
    (void)mode, (void)pattern, (void)sphere_radius, (void)cell_size,
        (void)wire_radius, (void)softness, (void)near_fade, (void)far_distance,
        (void)aa_strength, (void)speed, (void)spin_3d, (void)spin_4d,
        (void)color, (void)shells;
  }
};

// Width pin; pin_field_set() is what catches an added or removed field. Every
// enum has a fixed uint8_t base, so the size holds under ARM -fshort-enums.
static_assert(sizeof(Params) == 48, "HyperLattice::Params width changed");

struct FrameState {
  Params params;
  math::Vec4 origin;
  std::array<float, 6> rotation_phase;
  float pixel_half_angle;
  const BakedPalette *depth_palette;
  const BakedPalette *axis_palette;
};

struct Binding {
  using FrameState = HyperLatticeDetail::FrameState;
  using Instrumentation = Pullback::NoInstrumentation;
};

struct PreparedTrace {
  SDF::Lattice::PreparedTrace lattice;
  Raycast::Appearance appearance;
};
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
HS_FLASH_INLINE inline HyperLatticeExperimental::Settings
experimental_settings(const FrameState &frame, float phase) {
  const auto &p = frame.params;
  const math::Vec4 CENTER{
      {.255f + .3f * sinf(phase), .465f + .225f * sinf(2 * phase),
       .645f + .27f * sinf(3 * phase),
       p.mode == LatticeMode::FOUR_D_SLICE ? .375f + .21f * sinf(5 * phase)
                                           : 0.0f}};
  return {p.mode == LatticeMode::FOUR_D_SLICE
              ? Raycast::SamplingDomain::SLICE_4D
              : Raycast::SamplingDomain::SPATIAL_3D,
          p.cell_size,
          p.wire_radius * p.cell_size,
          p.sphere_radius,
          p.far_distance,
          p.near_fade,
          p.aa_strength,
          CENTER,
          view_embedding(frame),
          frame.pixel_half_angle,
          frame.depth_palette,
          frame.axis_palette,
          p.color};
}
#endif
inline PreparedTrace prepare_trace(const FrameState &frame) {
  const auto embedding = view_embedding(frame);
  const auto &p = frame.params;
  const SDF::Lattice::Settings settings{
      p.mode,     p.sphere_radius, p.cell_size, p.wire_radius,
      p.softness, p.aa_strength,   p.shells};
  const float near_scale = p.cell_size * (1 + p.sphere_radius);
  return {SDF::Lattice::prepare(settings, frame.origin, embedding,
                                p.far_distance, frame.pixel_half_angle),
          {1.0f / p.far_distance, 1.5f * p.wire_radius * near_scale,
           1.0f / (p.near_fade * near_scale), p.color, frame.depth_palette,
           frame.axis_palette}};
}
template <bool SLICE_4D = false, uint8_t SHELLS = 0> struct Renderer {
  static PreparedTrace prepare(const FrameState &frame) {
    return prepare_trace(frame);
  }
  __attribute__((always_inline)) static Color4
  shade(const math::Vector &normal, const FrameState &,
        const PreparedTrace &prepared) {
    HS_PROFILE_DEEP(hl_shade);
    SDF::Lattice::Events<SLICE_4D, SHELLS> events(normal, prepared.lattice);
    Raycast::TraceLimits limits;
    limits.max_candidates =
        DIMENSIONS * (SHELLS ? SHELLS : SDF::Lattice::MAX_SHELLS);
    return Raycast::shade_events<SLICE_4D>(events,
                                           {0, prepared.lattice.far_distance},
                                           limits, prepared.appearance)
        .color;
  }
};
using RenderPipeline =
    Pullback::Pipeline<Binding, Pullback::RayStage<Renderer<>>>;
template <uint8_t SHELL_COUNT>
using SpecializedRenderPipeline =
    Pullback::Pipeline<Binding, Pullback::RayStage<Renderer<true>>>;
} // namespace HyperLatticeDetail

/**
 * @brief Flights through cubic and octet lattices and their 4D slices.
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
  using Params = HyperLatticeDetail::Params;
  using LatticeMode = HyperLatticeDetail::LatticeMode;
  using ColorMode = HyperLatticeDetail::ColorMode;
  using ShellCount = HyperLatticeDetail::ShellCount;
  using Pattern = HyperLatticeDetail::Pattern;
  using ConfigurationId = HyperLatticeDetail::ConfigurationId;

  static constexpr auto PRESET_IDS = std::to_array<std::string_view>({
      "cubic-flight",
      "hypercube-flight",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "experimental-octet-flight",
      "experimental-octet-4d-slice",
#endif
      "cubic-wide-flight",
  });
  static constexpr size_t WIDE_PRESET_INDEX = PRESET_IDS.size() - 1;
  static constexpr Segue::Preset::Lerp PRESET_SEGUE{240, math::ease_in_out_sin,
                                                    /*pausable=*/true};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 320;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 12;

  HS_COLD_MEMBER static constexpr Params preset_params(size_t index) {
    Params value;
    switch (index) {
    case 0:
    case WIDE_PRESET_INDEX:
      value.mode = LatticeMode::THREE_D;
      value.sphere_radius = index == 0 ? 1.0f : 0.0f;
      value.cell_size = index == 0 ? 1.0f : 2.38525f;
      value.wire_radius = 0.055f;
      value.softness = 0.08f;
      value.near_fade = index == 0 ? 0.5f : 2.0f;
      value.far_distance = index == 0 ? 4.198f : 11.66f;
      value.aa_strength = 1.0f;
      value.speed = 0.05f;
      value.spin_3d = 0.015f;
      value.spin_4d = 0.0f;
      value.color = ColorMode::DEPTH;
      value.shells = ShellCount::TWO;
      break;
    case 1:
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
      value.color = ColorMode::DEPTH;
      value.shells = ShellCount::TWO;
      break;
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    case 2:
    case 3:
      value.pattern = Pattern::OCTET;
      value.mode =
          index == 2 ? LatticeMode::THREE_D : LatticeMode::FOUR_D_SLICE;
      value.sphere_radius = 0;
      value.cell_size = 1.5f;
      value.wire_radius = .055f;
      value.softness = .012f;
      value.far_distance = 4.5f;
      value.near_fade = .08f;
      value.speed = .008f;
      value.spin_3d = .0024f;
      value.spin_4d = index == 3 ? .0024f : 0.0f;
      value.color = ColorMode::DEPTH;
      break;
#endif
    default:
      break;
    }
    return value;
  }

  enum class Backend : uint8_t { ANALYTIC_EVENTS };
  enum class Policy : uint8_t { LEGACY_COVERAGE, EXPERIMENTAL };
  struct Configuration {
    Pattern pattern;
    LatticeMode domain;
    Backend backend;
    Policy policy;
    uint8_t default_preset;
    uint8_t max_candidates;
    uint8_t max_layers;
  };
  static constexpr auto CONFIGURATIONS = std::to_array<Configuration>({
      {Pattern::CUBIC_WIRE, LatticeMode::THREE_D, Backend::ANALYTIC_EVENTS,
       Policy::LEGACY_COVERAGE, 0, 9, 9},
      {Pattern::CUBIC_WIRE, LatticeMode::FOUR_D_SLICE, Backend::ANALYTIC_EVENTS,
       Policy::LEGACY_COVERAGE, 1, 12, 12},
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      {Pattern::OCTET, LatticeMode::THREE_D, Backend::ANALYTIC_EVENTS,
       Policy::EXPERIMENTAL, 2, 64, 32},
      {Pattern::OCTET, LatticeMode::FOUR_D_SLICE, Backend::ANALYTIC_EVENTS,
       Policy::EXPERIMENTAL, 3, 64, 32},
#endif
  });

  static constexpr ConfigurationId configuration_id(const Params &value) {
    return static_cast<ConfigurationId>(
        2 * static_cast<uint8_t>(value.pattern) +
        static_cast<uint8_t>(value.mode));
  }

  static constexpr bool supported_combination(const Params &value) {
    const auto index = static_cast<size_t>(configuration_id(value));
    return static_cast<uint8_t>(value.pattern) <=
               static_cast<uint8_t>(Pattern::OCTET) &&
           static_cast<uint8_t>(value.mode) <=
               static_cast<uint8_t>(LatticeMode::FOUR_D_SLICE) &&
           index < CONFIGURATIONS.size() &&
           CONFIGURATIONS[index].pattern == value.pattern &&
           CONFIGURATIONS[index].domain == value.mode;
  }

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
           static_cast<uint8_t>(value.color) <=
               static_cast<uint8_t>(ColorMode::AXIS) &&
           static_cast<uint8_t>(value.shells) <=
               static_cast<uint8_t>(ShellCount::THREE);
  }

  /**
   * @brief Whether a parameter set matches a specialized slice pipeline.
   * @param value Parameters to test.
   * @return true when the parameters match the specialized trace assumptions.
   * @details The shape is preset 1's; the assert in draw_frame() ties the two,
   *          so retuning that preset cannot leave the gate behind.
   */
  static constexpr bool uses_specialized_slice(const Params &value) {
    return value.pattern == Pattern::CUBIC_WIRE &&
           value.mode == LatticeMode::FOUR_D_SLICE &&
           value.color == ColorMode::DEPTH &&
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
    register_animated_param("Color", &params.color, COLOR_OPTIONS,
                            COLOR_EXPORT_OPTIONS, std::size(COLOR_OPTIONS));
    register_animated_param("Lattice Planes", &params.shells, SHELL_OPTIONS,
                            SHELL_EXPORT_OPTIONS, std::size(SHELL_OPTIONS));
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    this->register_readonly_param("Unfinished Rays", &unfinished_rays, 0,
                                  (W + 2) * (H + 2));
#endif
    refresh_configuration_schema();
    depth_palette.init_generated(persistent_arena, next_depth_palette, nullptr,
                                 0, PALETTE_FADE_FRAMES, math::ease_in_out_sin);
    const GenerativePalette fixed_axis_palette{
        EffectPaletteRecipes::hyper_lattice()};
    axis_palette.bake(persistent_arena, fixed_axis_palette);
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
        &axis_palette.view()};
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    unfinished_rays = 0;
    if (params.pattern != Pattern::CUBIC_WIRE) {
      HS_PROFILE(hl_shader_draw);
      draw_experimental(canvas, context);
      return;
    }
#endif
    const auto frame = HyperLatticeDetail::RenderPipeline::prepare(context);
    {
      HS_PROFILE(hl_shader_draw);
      static_assert(uses_specialized_slice(preset_params(1)),
                    "preset 1 no longer selects the specialized slice trace");
      if (uses_specialized_slice(frame.ctx.params)) {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::SpecializedRenderPipeline<2>::evaluate(
                  view, frame.ctx, frame.prepared);
            });
      } else {
        Scan::Shader::draw_cached<W, H, 1>(
            canvas, [&frame](const math::Vector &view) HS_HOT_FLASH_MEMBER {
              return HyperLatticeDetail::RenderPipeline::evaluate(
                  view, frame.ctx, frame.prepared);
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
    if (parameter.target == &params.pattern)
      candidate.pattern = static_cast<Pattern>(value);
    else if (parameter.target == &params.mode)
      candidate.mode = static_cast<LatticeMode>(value);
    else
      return true;
    return supported_combination(candidate);
  }
#endif

  HS_COLD_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    const auto configuration = configuration_id(params);
    if (selected_configuration != configuration) {
      const auto color = params.color;
      const float near_fade = params.near_fade;
      params = preset_params(
          CONFIGURATIONS[static_cast<size_t>(configuration)].default_preset);
      params.color = color;
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
  void blend_params(float progress) {
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
#endif
  }

  HS_FLASH_MEMBER void advance_state() {
    static constexpr float VELOCITY[HyperLatticeDetail::DIMENSIONS] = {
        0.4815434f, 0.2993373f, 0.4034555f, 0.7223151f};
    static constexpr float RATE[6] = {1.0f, 0.731f, 0.517f,
                                      1.0f, 0.707f, 0.419f};
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
    if (params.pattern != Pattern::CUBIC_WIRE)
      experimental_phase =
          math::wrap(experimental_phase + params.speed, math::TWO_PI_F);
    else
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
    using namespace HyperLatticeExperimental;
    const auto SETTINGS =
        HyperLatticeDetail::experimental_settings(context, experimental_phase);
    const auto prepared = prepare(SETTINGS);
    using Shade =
        Raycast::ShadedTrace (*)(const math::Vector &, const Prepared &);
    const Shade shade_ray =
        params.mode == LatticeMode::FOUR_D_SLICE ? &shade<true> : &shade<false>;
    Scan::Shader::draw_cached<W, H, 1>(
        canvas, [&prepared, shade_ray, this](const math::Vector &view)
                    HS_HOT_FLASH_MEMBER {
                      const auto result = shade_ray(view, prepared);
                      const auto STATUS = result.trace.status;
                      if (STATUS != Raycast::TraceStatus::SURFACE &&
                          STATUS != Raycast::TraceStatus::RANGE_COMPLETE &&
                          STATUS != Raycast::TraceStatus::SATURATED)
                        unfinished_rays += 1;
                      return result.color;
                    });
  }

  float experimental_phase = 0;
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
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 0.05f;
  static constexpr float SPIN_3D_MIN = 0.0f, SPIN_3D_MAX = 0.015f;
  static constexpr float SPIN_4D_MIN = 0.0f, SPIN_4D_MAX = 0.015f;

  static constexpr int PALETTE_FADE_FRAMES = 960;

  static constexpr const char *PATTERN_OPTIONS[] = {
      "Cubic",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "Experimental / Octet Truss",
#endif
  };
  static constexpr const char *PATTERN_EXPORT_OPTIONS[] = {
      "Pattern::CUBIC_WIRE",
#if HS_ENABLE_HYPERLATTICE_EXPERIMENTS
      "Pattern::OCTET",
#endif
  };
  static constexpr const char *VIEW_OPTIONS[] = {"3D perspective", "4D slice"};
  static constexpr const char *VIEW_EXPORT_OPTIONS[] = {
      "LatticeMode::THREE_D", "LatticeMode::FOUR_D_SLICE"};
  static constexpr const char *COLOR_OPTIONS[] = {"Depth", "Axis"};
  static constexpr const char *COLOR_EXPORT_OPTIONS[] = {"ColorMode::DEPTH",
                                                         "ColorMode::AXIS"};
  static constexpr const char *SHELL_OPTIONS[] = {"1", "2", "3"};
  static constexpr const char *SHELL_EXPORT_OPTIONS[] = {
      "ShellCount::ONE", "ShellCount::TWO", "ShellCount::THREE"};

  ConfigurationId selected_configuration = ConfigurationId::CUBIC_3D;
  math::Vec4 origin{{0.17f, 0.31f, 0.43f, 0.59f}};
  std::array<float, 6> rotation_phase{};
  PaletteCycler depth_palette;
  BakedPaletteStorage axis_palette;

  friend struct hs_test::hyper_lattice_tests::HyperLatticeWhiteBox;

  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes() +
      BakedPalette::required_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "HyperLattice persistent footprint exceeds the default "
                "partition");
};

static_assert(
    [] {
      using Effect = HyperLattice<1, 1>;
      for (size_t index = 0; index < Effect::PRESET_IDS.size(); ++index)
        if (!Effect::valid_params(Effect::preset_params(index)))
          return false;
      return true;
    }(),
    "HyperLattice preset is outside a registered slider range");
