/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <type_traits>

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/engine/fixed_pipeline.h"
#include "core/math/noise_field.h"
#include "core/render/pullback.h"

namespace FixedLook {

enum class HueMode : uint8_t { NONE, NOISE, PATH_LENGTH };

struct GridSourceParams {
  float pattern_freq = 1.0f;
  float speed = 0.0f;
  float complexity = 0.0f;
  float pattern_mix = 0.0f;
  float secondary_rate = 0.0f;
  float angle_rate = 0.0f;
};

struct TwinWaveSourceParams {
  float pattern_freq = 1.0f;
  float speed = 0.0f;
  float secondary_rate = 0.0f;
  float angle_rate = 0.0f;
};

struct NoiseSourceParams {
  float noise_scale = 1.0f;
  float noise_contrast = 0.0f;
  float noise_time_rate = 0.0f;
  int32_t noise_seed = 2927;
};

struct LatticeSourceParams {
  float lattice_cell_scale = 1.0f;
  float lattice_shape_blend = 0.0f;
  float lattice_softness = 0.05f;
  float lattice_radius = 0.25f;
};

struct ProjectionParams {
  float pole_fade = 1.0f;
  float spin_rate = 0.0f;
  float wander = 0.0f;
  float camera_wander = 0.0f;
  float central_meridian = 0.0f;
};

struct NoWarpParams {
  float speed = 0.0f;
};

struct MirrorParams {
  float speed = 0.0f;
  float rotation = 0.0f;
  float cell_x = 1.0f;
  float cell_y = 1.0f;
  float offset_x = 0.0f;
  float offset_y = 0.0f;
};

struct WaveShearParams {
  float speed = 0.0f;
  float strength = 0.0f;
  float frequency = 1.0f;
  float field_angle = 0.0f;
  float edge_width = 0.1f;
};

struct VectorNoiseParams {
  float speed = 0.0f;
  float strength = 0.0f;
  float scale = 1.0f;
  float vector_angle = 0.0f;
  float edge_width = 0.1f;
  int32_t seed = 1337;
};

struct AffineParams {
  float speed = 0.0f;
  float rotation_rate = 0.0f;
  float translation_x = 0.0f;
  float translation_y = 0.0f;
  float scale_x = 1.0f;
  float scale_y = 1.0f;
  float shear = 0.0f;
};

struct PolarParams {
  float speed = 0.0f;
  float radial_scale = 1.0f;
  float radial_phase = 0.0f;
  float angular_phase = 0.0f;
};

struct NoLensParams {};
struct MobiusLensParams {
  MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                      0.7071067811865475f, 0.0f};
};

struct LinearValueParams {};
struct EdgeValueParams {
  float edge_width = 0.1f;
};
struct IsoValueParams {
  float iso_level = 0.5f;
  float iso_width = 0.05f;
};

struct ColorParams {
  float hue_shift_amount = 0.0f;
  float hue_noise_scale = 1.0f;
  float hue_noise_speed = 0.0f;
  float palette_chroma = 0.62f;
  float mapping_frequency = 1.0f;
  float mapping_phase = 0.0f;
  float phase_oscillation_depth = 0.0f;
  float phase_oscillation_speed = 0.0f;
  float brightness_depth = 1.0f;
  float opacity_low = 1.0f;
  float opacity_high = 1.0f;
  Pullback::Color::PaletteMapping palette_mapping =
      Pullback::Color::PaletteMapping::LINEAR;
};

template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT = NoLensParams, typename ValueT = LinearValueParams>
struct Params {
  using source_type = SourceT;
  using outer_warp_type = OuterWarpT;
  using inner_warp_type = InnerWarpT;
  using lens_type = LensT;
  using value_type = ValueT;
  SourceT source;
  ProjectionParams projection;
  OuterWarpT outer_warp;
  InnerWarpT inner_warp;
  LensT lens;
  ValueT value;
  ColorParams color;
};

struct PreparedSource {
  float primary;
  float secondary;
  float angle;
  float angle_cos;
  float angle_sin;
};

struct PreparedAffine {
  float translation_x;
  float translation_y;
  float scale_x;
  float scale_y;
  float shear;
};

struct PreparedMirror {
  float offset_x;
  float offset_y;
};

struct PreparedNoiseLoop {
  float diagonal;
  float z;
};

union PreparedWarpTransform {
  PreparedAffine affine;
  PreparedMirror mirror;
  PreparedNoiseLoop noise_loop;
};

struct PreparedWarp {
  float rotation_cos;
  float rotation_sin;
  PreparedWarpTransform transform;
};

template <typename ParamsT> struct FrameState {
  Quaternion projection_conjugate;
  Quaternion outer_conjugate;
  PreparedSource prepared_source;
  PreparedWarp prepared_outer;
  PreparedWarp prepared_inner;
  const FastNoiseLite *outer_noise;
  const FastNoiseLite *source_noise;
  const BakedPalette *palette;
  const Pixel *hue_rotation_lut;
  const int8_t *hue_noise_lut;
  ParamsT params;
  Pullback::Color::PaletteMappingWeights palette_mapping;
  float outer_phase;
  float inner_phase;
  float source_noise_time;
  float palette_oscillation_phase;
};

template <typename FrameT> struct Binding {
  using FrameState = FrameT;
  using Instrumentation = Pullback::NoInstrumentation;
};

template <typename BindingT> struct OuterCameraProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const Quaternion &conjugate(const FrameState &frame) {
    return frame.outer_conjugate;
  }
};

template <typename BindingT> struct ProjectionProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const Quaternion &conjugate(const FrameState &frame) {
    return frame.projection_conjugate;
  }
  static float pole_fade(const FrameState &frame) {
    return frame.params.projection.pole_fade;
  }
  static float central_meridian(const FrameState &frame) {
    return frame.params.projection.central_meridian;
  }
};

template <typename BindingT> struct LensProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const MobiusParams &params(const FrameState &frame) {
    return frame.params.lens.mobius;
  }
};

template <typename BindingT, bool Outer, bool TrackPath = false>
struct WarpProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const auto &params(const FrameState &frame) {
    if constexpr (Outer)
      return frame.params.outer_warp;
    else
      return frame.params.inner_warp;
  }
  static const PreparedWarp &prepared(const FrameState &frame) {
    if constexpr (Outer)
      return frame.prepared_outer;
    else
      return frame.prepared_inner;
  }
  static float phase(const FrameState &frame) {
    if constexpr (Outer)
      return frame.outer_phase;
    else
      return frame.inner_phase;
  }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.outer_noise;
  }
  static bool path_length_required(const FrameState &) { return TrackPath; }
};

template <typename BindingT> struct SourceProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static const auto &params(const FrameState &frame) {
    return frame.params.source;
  }
  static const PreparedSource &prepared(const FrameState &frame) {
    return frame.prepared_source;
  }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.source_noise;
  }
  static float noise_scale(const FrameState &frame) {
    return frame.params.source.noise_scale;
  }
  static float noise_time(const FrameState &frame) {
    return frame.source_noise_time;
  }
  static float noise_contrast(const FrameState &frame) {
    return frame.params.source.noise_contrast;
  }
};

template <typename BindingT> struct ValueProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static float iso_level(const FrameState &frame) {
    return frame.params.value.iso_level;
  }
  static float iso_width(const FrameState &frame) {
    return frame.params.value.iso_width;
  }
  static float edge_width(const FrameState &frame) {
    return frame.params.value.edge_width;
  }
};

template <typename BindingT, HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV>
struct ColorProvider {
  using Binding = BindingT;
  using FrameState = typename Binding::FrameState;
  static Pullback::Color::PaletteMappingWeights
  mapping_weights(const FrameState &frame) {
    return frame.palette_mapping;
  }
  static float mapping_frequency(const FrameState &frame) {
    return frame.params.color.mapping_frequency;
  }
  static float mapping_phase(const FrameState &frame) {
    return frame.params.color.mapping_phase;
  }
  static float oscillation_depth(const FrameState &frame) {
    return frame.params.color.phase_oscillation_depth;
  }
  static float oscillation_phase(const FrameState &frame) {
    return frame.palette_oscillation_phase;
  }
  static const BakedPalette &palette(const FrameState &frame) {
    return *frame.palette;
  }
  static Pullback::Color::HueMode hue_mode(const FrameState &) {
    if constexpr (HueV == HueMode::NOISE)
      return Pullback::Color::HueMode::NOISE;
    if constexpr (HueV == HueMode::PATH_LENGTH)
      return Pullback::Color::HueMode::PATH_LENGTH;
    return Pullback::Color::HueMode::NONE;
  }
  static float hue_shift_amount(const FrameState &frame) {
    return frame.params.color.hue_shift_amount;
  }
  static Pullback::Color::HueRotationLutView
  hue_rotation(const FrameState &frame) {
    return {frame.hue_rotation_lut,
            HueV != HueMode::NONE &&
                frame.params.color.hue_shift_amount != 0.0f};
  }
  static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
    return {frame.hue_noise_lut,
            HueV == HueMode::NOISE &&
                frame.params.color.hue_shift_amount != 0.0f};
  }
  static Pullback::Color::BrightnessEnvelope
  brightness_envelope(const FrameState &) {
    return BrightnessV;
  }
  static float brightness_depth(const FrameState &frame) {
    return frame.params.color.brightness_depth;
  }
  static float opacity_low(const FrameState &frame) {
    return frame.params.color.opacity_low;
  }
  static float opacity_high(const FrameState &frame) {
    return frame.params.color.opacity_high;
  }
};

inline float lerp(float from, float to, float progress) {
  return FixedPipeline::linear(from, to, progress);
}

inline GridSourceParams interpolate(const GridSourceParams &a,
                                    const GridSourceParams &b, float t) {
  return {lerp(a.pattern_freq, b.pattern_freq, t),
          lerp(a.speed, b.speed, t),
          lerp(a.complexity, b.complexity, t),
          lerp(a.pattern_mix, b.pattern_mix, t),
          lerp(a.secondary_rate, b.secondary_rate, t),
          lerp(a.angle_rate, b.angle_rate, t)};
}

inline TwinWaveSourceParams interpolate(const TwinWaveSourceParams &a,
                                        const TwinWaveSourceParams &b,
                                        float t) {
  return {lerp(a.pattern_freq, b.pattern_freq, t), lerp(a.speed, b.speed, t),
          lerp(a.secondary_rate, b.secondary_rate, t),
          lerp(a.angle_rate, b.angle_rate, t)};
}

inline NoiseSourceParams interpolate(const NoiseSourceParams &a,
                                     const NoiseSourceParams &b, float t) {
  return {lerp(a.noise_scale, b.noise_scale, t),
          lerp(a.noise_contrast, b.noise_contrast, t),
          lerp(a.noise_time_rate, b.noise_time_rate, t),
          t < 1.0f ? a.noise_seed : b.noise_seed};
}

inline LatticeSourceParams interpolate(const LatticeSourceParams &a,
                                       const LatticeSourceParams &b, float t) {
  return {
      FixedPipeline::log_positive(a.lattice_cell_scale, b.lattice_cell_scale,
                                  t),
      lerp(a.lattice_shape_blend, b.lattice_shape_blend, t),
      FixedPipeline::log_positive(a.lattice_softness, b.lattice_softness, t),
      lerp(a.lattice_radius, b.lattice_radius, t)};
}

inline ProjectionParams interpolate(const ProjectionParams &a,
                                    const ProjectionParams &b, float t) {
  return {lerp(a.pole_fade, b.pole_fade, t), lerp(a.spin_rate, b.spin_rate, t),
          lerp(a.wander, b.wander, t),
          lerp(a.camera_wander, b.camera_wander, t),
          FixedPipeline::shortest_periodic(a.central_meridian,
                                           b.central_meridian, t, TWO_PI_F)};
}

inline NoWarpParams interpolate(const NoWarpParams &a, const NoWarpParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t)};
}

inline MirrorParams interpolate(const MirrorParams &a, const MirrorParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t),
          FixedPipeline::shortest_periodic(a.rotation, b.rotation, t, TWO_PI_F),
          FixedPipeline::log_positive(a.cell_x, b.cell_x, t),
          FixedPipeline::log_positive(a.cell_y, b.cell_y, t),
          lerp(a.offset_x, b.offset_x, t),
          lerp(a.offset_y, b.offset_y, t)};
}

inline WaveShearParams interpolate(const WaveShearParams &a,
                                   const WaveShearParams &b, float t) {
  return {lerp(a.speed, b.speed, t), lerp(a.strength, b.strength, t),
          FixedPipeline::log_positive(a.frequency, b.frequency, t),
          FixedPipeline::shortest_periodic(a.field_angle, b.field_angle, t,
                                           TWO_PI_F),
          lerp(a.edge_width, b.edge_width, t)};
}

inline VectorNoiseParams interpolate(const VectorNoiseParams &a,
                                     const VectorNoiseParams &b, float t) {
  return {lerp(a.speed, b.speed, t),
          lerp(a.strength, b.strength, t),
          FixedPipeline::log_positive(a.scale, b.scale, t),
          FixedPipeline::shortest_periodic(a.vector_angle, b.vector_angle, t,
                                           TWO_PI_F),
          lerp(a.edge_width, b.edge_width, t),
          t < 1.0f ? a.seed : b.seed};
}

inline AffineParams interpolate(const AffineParams &a, const AffineParams &b,
                                float t) {
  return {lerp(a.speed, b.speed, t),
          lerp(a.rotation_rate, b.rotation_rate, t),
          lerp(a.translation_x, b.translation_x, t),
          lerp(a.translation_y, b.translation_y, t),
          FixedPipeline::log_positive(a.scale_x, b.scale_x, t),
          FixedPipeline::log_positive(a.scale_y, b.scale_y, t),
          lerp(a.shear, b.shear, t)};
}

inline PolarParams interpolate(const PolarParams &a, const PolarParams &b,
                               float t) {
  return {lerp(a.speed, b.speed, t),
          FixedPipeline::log_positive(a.radial_scale, b.radial_scale, t),
          FixedPipeline::shortest_periodic(a.radial_phase, b.radial_phase, t,
                                           TWO_PI_F),
          FixedPipeline::shortest_periodic(a.angular_phase, b.angular_phase, t,
                                           TWO_PI_F)};
}

inline NoLensParams interpolate(const NoLensParams &, const NoLensParams &,
                                float) {
  return {};
}

inline MobiusLensParams interpolate(const MobiusLensParams &a,
                                    const MobiusLensParams &b, float t) {
  MobiusLensParams value;
  value.mobius = t < 1.0f ? a.mobius : b.mobius;
  return value;
}

inline LinearValueParams interpolate(const LinearValueParams &,
                                     const LinearValueParams &, float) {
  return {};
}

inline EdgeValueParams interpolate(const EdgeValueParams &a,
                                   const EdgeValueParams &b, float t) {
  return {lerp(a.edge_width, b.edge_width, t)};
}

inline IsoValueParams interpolate(const IsoValueParams &a,
                                  const IsoValueParams &b, float t) {
  return {lerp(a.iso_level, b.iso_level, t),
          FixedPipeline::log_positive(a.iso_width, b.iso_width, t)};
}

inline ColorParams interpolate(const ColorParams &a, const ColorParams &b,
                               float t) {
  return {
      lerp(a.hue_shift_amount, b.hue_shift_amount, t),
      FixedPipeline::log_positive(a.hue_noise_scale, b.hue_noise_scale, t),
      lerp(a.hue_noise_speed, b.hue_noise_speed, t),
      lerp(a.palette_chroma, b.palette_chroma, t),
      FixedPipeline::log_positive(a.mapping_frequency, b.mapping_frequency, t),
      lerp(a.mapping_phase, b.mapping_phase, t),
      lerp(a.phase_oscillation_depth, b.phase_oscillation_depth, t),
      lerp(a.phase_oscillation_speed, b.phase_oscillation_speed, t),
      lerp(a.brightness_depth, b.brightness_depth, t),
      lerp(a.opacity_low, b.opacity_low, t),
      lerp(a.opacity_high, b.opacity_high, t),
      t < 1.0f ? a.palette_mapping : b.palette_mapping};
}

template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT, typename ValueT>
inline Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT>
interpolate(const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT> &from,
            const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT> &to,
            float progress) {
  return {interpolate(from.source, to.source, progress),
          interpolate(from.projection, to.projection, progress),
          interpolate(from.outer_warp, to.outer_warp, progress),
          interpolate(from.inner_warp, to.inner_warp, progress),
          interpolate(from.lens, to.lens, progress),
          interpolate(from.value, to.value, progress),
          interpolate(from.color, to.color, progress)};
}

inline bool finite_range(float value, float minimum, float maximum) {
  return std::isfinite(value) && value >= minimum && value <= maximum;
}

inline bool valid(const GridSourceParams &p) {
  return finite_range(p.pattern_freq, 0.1f, 20.0f) &&
         finite_range(p.speed, 0.0f, 0.5f) &&
         finite_range(p.complexity, 0.0f, 3.0f) &&
         finite_range(p.pattern_mix, 0.0f, 1.0f) &&
         finite_range(p.secondary_rate, 0.0f, 1.25f) &&
         finite_range(p.angle_rate, 0.0f, 0.05f);
}

inline bool valid(const TwinWaveSourceParams &p) {
  return finite_range(p.pattern_freq, 0.1f, 20.0f) &&
         finite_range(p.speed, 0.0f, 0.5f) &&
         finite_range(p.secondary_rate, 0.0f, 1.25f) &&
         finite_range(p.angle_rate, 0.0f, 0.05f);
}

inline bool valid(const NoiseSourceParams &p) {
  return finite_range(p.noise_scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.noise_contrast, 0.0f, 8.0f) &&
         finite_range(p.noise_time_rate, -1.0f / 64.0f, 1.0f / 64.0f);
}

inline bool valid(const LatticeSourceParams &p) {
  return finite_range(p.lattice_cell_scale, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.lattice_shape_blend, 0.0f, 1.0f) &&
         finite_range(p.lattice_softness, 1.0f / 1024.0f, 1.0f) &&
         finite_range(p.lattice_radius, 1.0f / 64.0f, 0.49f);
}

inline bool valid(const ProjectionParams &p) {
  return finite_range(p.pole_fade, 1.0f, 20.0f) &&
         finite_range(p.spin_rate, 0.0f, 0.05f) &&
         finite_range(p.wander, 0.0f, 1.0f) &&
         finite_range(p.camera_wander, 0.0f, 1.0f) &&
         finite_range(p.central_meridian, 0.0f, TWO_PI_F);
}

inline bool valid(const NoWarpParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f);
}

inline bool valid(const MirrorParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.rotation, 0.0f, TWO_PI_F) &&
         finite_range(p.cell_x, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.cell_y, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.offset_x, -8.0f, 8.0f) &&
         finite_range(p.offset_y, -8.0f, 8.0f);
}

inline bool valid(const WaveShearParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.strength, -30.0f, 30.0f) &&
         finite_range(p.frequency, 0.01f, 32.0f) &&
         finite_range(p.field_angle, 0.0f, TWO_PI_F) &&
         finite_range(p.edge_width, 0.0f, 1.0f);
}

inline bool valid(const VectorNoiseParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.strength, -30.0f, 30.0f) &&
         finite_range(p.scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.vector_angle, 0.0f, TWO_PI_F) &&
         finite_range(p.edge_width, 0.0f, 1.0f);
}

inline bool valid(const AffineParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.rotation_rate, -TWO_PI_F, TWO_PI_F) &&
         finite_range(p.translation_x, -4.0f, 4.0f) &&
         finite_range(p.translation_y, -4.0f, 4.0f) &&
         finite_range(p.scale_x, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.scale_y, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.shear, -4.0f, 4.0f);
}

inline bool valid(const PolarParams &p) {
  return finite_range(p.speed, -0.02f, 0.02f) &&
         finite_range(p.radial_scale, 1.0f / 64.0f, 64.0f) &&
         finite_range(p.radial_phase, -TWO_PI_F, TWO_PI_F) &&
         finite_range(p.angular_phase, -TWO_PI_F, TWO_PI_F);
}

inline bool valid(const NoLensParams &) { return true; }
inline bool valid(const MobiusLensParams &p) {
  const float values[] = {p.mobius.a.re, p.mobius.a.im, p.mobius.b.re,
                          p.mobius.b.im, p.mobius.c.re, p.mobius.c.im,
                          p.mobius.d.re, p.mobius.d.im};
  for (float value : values)
    if (!std::isfinite(value))
      return false;
  return true;
}
inline bool valid(const LinearValueParams &) { return true; }
inline bool valid(const EdgeValueParams &p) {
  return finite_range(p.edge_width, 0.0f, 1.0f);
}
inline bool valid(const IsoValueParams &p) {
  return finite_range(p.iso_level, 0.0f, 1.0f) &&
         finite_range(p.iso_width, 1.0f / 1024.0f, 1.0f);
}
inline bool valid(const ColorParams &p) {
  return finite_range(p.hue_shift_amount, -4.0f, 4.0f) &&
         finite_range(p.hue_noise_scale, 1.0f / 64.0f, 8.0f) &&
         finite_range(p.hue_noise_speed, -0.001f, 0.001f) &&
         finite_range(p.palette_chroma, 0.0f, 1.0f) &&
         finite_range(p.mapping_frequency, 1.0f, 32.0f) &&
         finite_range(p.mapping_phase, -1.0f, 1.0f) &&
         finite_range(p.phase_oscillation_depth, 0.0f, 1.0f) &&
         finite_range(p.phase_oscillation_speed, -0.01f, 0.01f) &&
         finite_range(p.brightness_depth, 0.0f, 1.0f) &&
         finite_range(p.opacity_low, 0.0f, 1.0f) &&
         finite_range(p.opacity_high, 0.0f, 1.0f) &&
         static_cast<uint8_t>(p.palette_mapping) <=
             static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE);
}

template <typename SourceT, typename OuterWarpT, typename InnerWarpT,
          typename LensT, typename ValueT>
inline bool
valid(const Params<SourceT, OuterWarpT, InnerWarpT, LensT, ValueT> &params) {
  return valid(params.source) && valid(params.projection) &&
         valid(params.outer_warp) && valid(params.inner_warp) &&
         valid(params.lens) && valid(params.value) && valid(params.color);
}

template <bool Enabled> struct OptionalNoise {};
template <> struct OptionalNoise<true> {
  FastNoiseLite noise;
};

template <int W, int H, typename Derived, typename ParamsT,
          PaletteHarmony Harmony, HueMode HueV,
          Pullback::Color::BrightnessEnvelope BrightnessV,
          bool HasOuterNoise = false, bool HasSourceNoise = false>
class Runtime : public Effect {
public:
  using Params = ParamsT;
  using Frame = FrameState<Params>;
  using PipelineBinding = Binding<Frame>;
  static constexpr uint16_t TRANSITION_DURATION = 480;

  HS_COLD_MEMBER Runtime() : Effect(W, H, {.strobe = true}) {}

  HS_COLD_MEMBER void init() override {
    configure_presets(Derived::PRESET_IDS.size());
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    configure_noise(state->color_noise, 6047);
    if constexpr (HasOuterNoise)
      configure_noise(state->outer.noise, Derived::OUTER_NOISE_SEED);
    if constexpr (HasSourceNoise)
      configure_noise(state->source.noise, Derived::SOURCE_NOISE_SEED);
    init_gamut_lut(persistent_arena, GAMUT_LUT_ANGLE_STEPS, GAMUT_LUT_L_STEPS);
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  ease_in_out_sin);
    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));
    register_parameters();
    if constexpr (requires(Derived &effect) { effect.after_fixed_init(); })
      static_cast<Derived &>(*this).after_fixed_init();
  }

  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(fx_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(fx_advance);
      begin_automatic_transition();
      prepare_transition_value();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const Frame frame = prepare_frame();
    {
      HS_PROFILE(fx_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas, [&frame](const Vector &view) {
        return Derived::shade(view, frame);
      });
    }
    finish_transition_evaluation();
  }

  struct ParameterSnapshot {
    uint32_t schema_version;
    Params params;
  };

  ParameterSnapshot serialize_parameters() const {
    return {Derived::PARAMETER_SCHEMA_VERSION, params};
  }

  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != Derived::PARAMETER_SCHEMA_VERSION ||
        !FixedLook::valid(snapshot.params))
      return false;
    transition.active = false;
    params = snapshot.params;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.color.palette_mapping);
    return true;
  }

protected:
  struct Preset {
    Params params;
  };

  HS_COLD_MEMBER void start_mobius_animation(float scale, int duration) {
    if constexpr (requires { params.lens.mobius; })
      timeline.add_pausable(0,
                            Animation::MobiusWarpCircular(
                                params.lens.mobius, scale, duration, true),
                            &anims_paused);
  }

  HS_COLD_MEMBER void hold_initial_preset(uint16_t frames) {
    preset_dwell_remaining = frames;
  }

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const Params target = Derived::preset_params(change.to);
    if (change.origin != PresetChangeOrigin::AUTOMATIC) {
      transition.active = false;
      params = target;
      palette_mapping = Pullback::Color::PaletteMappingWeights::single(
          target.color.palette_mapping);
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      return true;
    }
    transition = {params,
                  target,
                  palette_mapping,
                  Pullback::Color::PaletteMappingWeights::single(
                      target.color.palette_mapping),
                  0,
                  TRANSITION_DURATION,
                  true};
    return true;
  }

  Params params = Derived::initial_params();

private:
  static constexpr size_t PARAM_CAPACITY = 48;

  struct Transition {
    Params from{};
    Params to{};
    Pullback::Color::PaletteMappingWeights mapping_from{};
    Pullback::Color::PaletteMappingWeights mapping_to{};
    uint16_t evaluation = 0;
    uint16_t duration = TRANSITION_DURATION;
    bool active = false;
  };

  struct State {
    std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE>
        hue_rotation_lut;
    std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
    FastNoiseLite color_noise;
    OptionalNoise<HasOuterNoise> outer;
    OptionalNoise<HasSourceNoise> source;
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    float hue_noise_lut_scale = -1.0f;
    float hue_noise_lut_phase = -1.0f;
  };

  static void configure_noise(FastNoiseLite &noise, int32_t seed) {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(seed);
    noise.SetFrequency(1.0f);
  }

  template <typename T> void register_source(T &source) {
    if constexpr (requires { source.pattern_freq; })
      register_animated_param("Pattern Freq", &source.pattern_freq, 0.1f,
                              20.0f);
    if constexpr (requires { source.speed; })
      register_animated_param("Speed", &source.speed, 0.0f, 0.5f);
    if constexpr (requires { source.complexity; })
      register_animated_param("Complexity", &source.complexity, 0.0f, 3.0f);
    if constexpr (requires { source.pattern_mix; })
      register_animated_param("Pattern Mix", &source.pattern_mix, 0.0f, 1.0f);
    if constexpr (requires { source.secondary_rate; })
      register_animated_param("Drift", &source.secondary_rate, 0.0f, 1.25f);
    if constexpr (requires { source.angle_rate; })
      register_animated_param("Source Angle Speed", &source.angle_rate, 0.0f,
                              0.05f);
    if constexpr (requires { source.noise_scale; })
      register_animated_param("Source Noise Scale", &source.noise_scale,
                              1.0f / 64.0f, 64.0f);
    if constexpr (requires { source.noise_contrast; })
      register_animated_param("Source Noise Contrast", &source.noise_contrast,
                              0.0f, 8.0f);
    if constexpr (requires { source.noise_time_rate; })
      register_animated_param("Source Noise Speed", &source.noise_time_rate,
                              -1.0f / 64.0f, 1.0f / 64.0f);
    if constexpr (requires { source.lattice_cell_scale; }) {
      register_animated_param("Lattice Cell Scale", &source.lattice_cell_scale,
                              1.0f / 64.0f, 8.0f);
      register_animated_param("Lattice Shape", &source.lattice_shape_blend,
                              0.0f, 1.0f);
      register_animated_param("Lattice Softness", &source.lattice_softness,
                              1.0f / 1024.0f, 1.0f);
      register_animated_param("Lattice Radius", &source.lattice_radius,
                              1.0f / 64.0f, 0.49f);
    }
  }

  template <typename T> void register_warp(T &warp, const char *prefix) {
    if constexpr (!std::is_same_v<T, NoWarpParams>) {
      register_animated_param(prefix, &warp.speed, -0.02f, 0.02f);
      if constexpr (requires { warp.strength; })
        register_animated_param("Warp Strength", &warp.strength, -30.0f, 30.0f);
      if constexpr (requires { warp.frequency; })
        register_animated_param("Warp Frequency", &warp.frequency, 0.01f,
                                32.0f);
      if constexpr (requires { warp.field_angle; })
        register_animated_param("Warp Field Angle", &warp.field_angle, 0.0f,
                                TWO_PI_F);
      if constexpr (requires { warp.scale; })
        register_animated_param("Warp Scale", &warp.scale, 1.0f / 64.0f, 64.0f);
      if constexpr (requires { warp.vector_angle; })
        register_animated_param("Warp Vector Angle", &warp.vector_angle, 0.0f,
                                TWO_PI_F);
      if constexpr (requires { warp.rotation; }) {
        register_animated_param("Mirror Rotation", &warp.rotation, 0.0f,
                                TWO_PI_F);
        register_animated_param("Mirror Cell X", &warp.cell_x, 1.0f / 64.0f,
                                8.0f);
        register_animated_param("Mirror Cell Y", &warp.cell_y, 1.0f / 64.0f,
                                8.0f);
        register_animated_param("Mirror Offset X", &warp.offset_x, -8.0f, 8.0f);
        register_animated_param("Mirror Offset Y", &warp.offset_y, -8.0f, 8.0f);
      }
      if constexpr (requires { warp.rotation_rate; }) {
        register_animated_param("Affine Rotation Rate", &warp.rotation_rate,
                                -TWO_PI_F, TWO_PI_F);
        register_animated_param("Affine Translation X", &warp.translation_x,
                                -4.0f, 4.0f);
        register_animated_param("Affine Translation Y", &warp.translation_y,
                                -4.0f, 4.0f);
        register_animated_param("Affine Scale X", &warp.scale_x, 1.0f / 64.0f,
                                64.0f);
        register_animated_param("Affine Scale Y", &warp.scale_y, 1.0f / 64.0f,
                                64.0f);
        register_animated_param("Affine Shear", &warp.shear, -4.0f, 4.0f);
      }
      if constexpr (requires { warp.radial_scale; }) {
        register_animated_param("Polar Radial Scale", &warp.radial_scale,
                                1.0f / 64.0f, 64.0f);
        register_animated_param("Polar Radial Phase", &warp.radial_phase,
                                -TWO_PI_F, TWO_PI_F);
        register_animated_param("Polar Angular Phase", &warp.angular_phase,
                                -TWO_PI_F, TWO_PI_F);
      }
    }
  }

  void register_parameters() {
    register_source(params.source);
    register_animated_param("Pole Fade", &params.projection.pole_fade, 1.0f,
                            20.0f);
    if constexpr (Derived::ANIMATED_PROJECTION) {
      register_animated_param("Projection Spin Speed",
                              &params.projection.spin_rate, 0.0f, 0.05f);
      register_animated_param("Projection Wander", &params.projection.wander,
                              0.0f, 1.0f);
    }
    register_animated_param("Camera Wander", &params.projection.camera_wander,
                            0.0f, 1.0f);
    if constexpr (requires { Derived::USES_CENTRAL_MERIDIAN; })
      if constexpr (Derived::USES_CENTRAL_MERIDIAN)
        register_animated_param("Central Meridian",
                                &params.projection.central_meridian, 0.0f,
                                TWO_PI_F);
    register_warp(params.outer_warp, "Planar Warp 1 Speed");
    register_warp(params.inner_warp, "Planar Warp 2 Speed");
    if constexpr (requires { params.value.edge_width; })
      register_animated_param("Edge Width", &params.value.edge_width, 0.0f,
                              1.0f);
    if constexpr (requires { params.value.iso_level; }) {
      register_animated_param("Iso Level", &params.value.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.value.iso_width,
                              1.0f / 1024.0f, 1.0f);
    }
    if constexpr (requires { params.lens.mobius; }) {
      register_animated_param("Mobius A Re", &params.lens.mobius.a.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius A Im", &params.lens.mobius.a.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius B Re", &params.lens.mobius.b.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius B Im", &params.lens.mobius.b.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius C Re", &params.lens.mobius.c.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius C Im", &params.lens.mobius.c.im, -4.0f,
                              4.0f);
      register_animated_param("Mobius D Re", &params.lens.mobius.d.re, -4.0f,
                              4.0f);
      register_animated_param("Mobius D Im", &params.lens.mobius.d.im, -4.0f,
                              4.0f);
    }
    register_animated_param("Palette Chroma", &params.color.palette_chroma,
                            0.0f, 1.0f);
    register_animated_param("Palette Mapping", &params.color.palette_mapping,
                            PALETTE_MAPPING_OPTIONS,
                            PALETTE_MAPPING_EXPORT_OPTIONS,
                            std::size(PALETTE_MAPPING_OPTIONS));
    register_animated_param("Mapping Frequency",
                            &params.color.mapping_frequency, 1.0f, 32.0f);
    register_animated_param("Mapping Phase", &params.color.mapping_phase, -1.0f,
                            1.0f);
    register_animated_param("Phase Oscillation Depth",
                            &params.color.phase_oscillation_depth, 0.0f, 1.0f);
    register_animated_param("Phase Oscillation Speed",
                            &params.color.phase_oscillation_speed, -0.01f,
                            0.01f);
    if constexpr (BrightnessV != Pullback::Color::BrightnessEnvelope::NONE)
      register_animated_param("Brightness Depth",
                              &params.color.brightness_depth, 0.0f, 1.0f);
    register_animated_param("Value Opacity Low", &params.color.opacity_low,
                            0.0f, 1.0f);
    register_animated_param("Value Opacity High", &params.color.opacity_high,
                            0.0f, 1.0f);
    register_animated_param("Hue Shift Amount", &params.color.hue_shift_amount,
                            HueV == HueMode::PATH_LENGTH ? -4.0f : -1.0f,
                            HueV == HueMode::PATH_LENGTH ? 4.0f : 1.0f);
    if constexpr (HueV == HueMode::NOISE) {
      register_animated_param("Hue Noise Scale", &params.color.hue_noise_scale,
                              1.0f / 64.0f, 8.0f);
      register_animated_param("Hue Noise Speed", &params.color.hue_noise_speed,
                              -0.001f, 0.001f);
    }
  }

  HS_COLD_MEMBER void begin_automatic_transition() {
    if constexpr (Derived::PRESET_IDS.size() == 1)
      return;
    if (anims_paused || transition.active)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    if (advancePreset()) {
#ifdef HS_PROFILE_ENABLE
      hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
              static_cast<unsigned>(getPresetCount()));
#endif
    }
  }

  HS_COLD_MEMBER void prepare_transition_value() {
    if (!transition.active)
      return;
    MobiusParams animated_mobius;
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        animated_mobius = params.lens.mobius;
    const FixedPipeline::EdgeProgress progress =
        FixedPipeline::edge_progress(transition.evaluation, transition.duration,
                                     FixedPipeline::Easing::EASE_IN_OUT_SIN);
    params =
        FixedLook::interpolate(transition.from, transition.to, progress.eased);
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        params.lens.mobius = animated_mobius;
    palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
        transition.mapping_from, transition.mapping_to, progress.eased);
  }

  HS_COLD_MEMBER void finish_transition_evaluation() {
    if (!transition.active || anims_paused)
      return;
    if (transition.evaluation == transition.duration) {
      MobiusParams animated_mobius;
      if constexpr (requires { Derived::ANIMATED_MOBIUS; })
        if constexpr (Derived::ANIMATED_MOBIUS)
          animated_mobius = params.lens.mobius;
      params = transition.to;
      if constexpr (requires { Derived::ANIMATED_MOBIUS; })
        if constexpr (Derived::ANIMATED_MOBIUS)
          params.lens.mobius = animated_mobius;
      palette_mapping = transition.mapping_to;
      transition.active = false;
      preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
      return;
    }
    ++transition.evaluation;
  }

  HS_COLD_MEMBER void advance_runtime() {
    if constexpr (requires { params.source.speed; }) {
      source_primary = fmodf(source_primary + params.source.speed, TWO_PI_F);
      if constexpr (requires { params.source.secondary_rate; })
        source_secondary =
            fmodf(source_secondary +
                      params.source.speed * params.source.secondary_rate,
                  TWO_PI_F);
      if constexpr (requires { params.source.angle_rate; })
        source_angle = fmodf(source_angle + params.source.angle_rate, TWO_PI_F);
    }
    if constexpr (requires { params.source.noise_time_rate; })
      source_noise_time =
          wrap_t(source_noise_time + params.source.noise_time_rate);
    projection_spin =
        fmodf(projection_spin + params.projection.spin_rate, TWO_PI_F);
    hue_noise_phase = wrap_t(hue_noise_phase + params.color.hue_noise_speed);
    if constexpr (std::is_same_v<typename Params::outer_warp_type,
                                 AffineParams>)
      outer_rotation =
          TWO_PI_F *
          wrap_t((outer_rotation +
                  params.outer_warp.speed * params.outer_warp.rotation_rate) /
                 TWO_PI_F);
    outer_phase = wrap_t(outer_phase + params.outer_warp.speed);
    inner_phase = wrap_t(inner_phase + params.inner_warp.speed);
    palette_oscillation_phase = wrap_t(palette_oscillation_phase +
                                       params.color.phase_oscillation_speed);
  }

  HS_COLD_MEMBER void update_spatial_frames() {
    const Quaternion projection = projection_walk.get();
    const Quaternion projection_delta =
        projection * projection_walk_previous.conjugate();
    projection_walk_previous = projection;
    projection_wander =
        (FixedPipeline::scaled_rotation_delta(projection_delta.normalized(),
                                              params.projection.wander) *
         projection_wander)
            .normalized();
    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander =
        (FixedPipeline::scaled_rotation_delta(outer_delta.normalized(),
                                              params.projection.camera_wander) *
         outer_wander)
            .normalized();
    projection_conjugate = (make_rotation(Y_AXIS, projection_spin) *
                            base_orientation * projection_wander)
                               .conjugate();
    outer_conjugate = outer_wander.conjugate();
  }

  template <typename WarpT>
  HS_COLD_MEMBER PreparedWarp prepare_warp(const WarpT &warp, float phase,
                                           bool outer) const {
    PreparedWarp prepared{};
    float rotation = 0.0f;
    if constexpr (requires { warp.rotation; }) {
      rotation = warp.rotation;
      prepared.transform.mirror = {
          wrap_t(warp.offset_x / warp.cell_x + phase) * warp.cell_x,
          wrap_t(warp.offset_y / warp.cell_y) * warp.cell_y};
    } else if constexpr (requires { warp.field_angle; }) {
      rotation = warp.field_angle;
    } else if constexpr (requires { warp.vector_angle; }) {
      rotation = warp.vector_angle;
      const float angle = TWO_PI_F * wrap_t(phase);
      prepared.transform.noise_loop = {NOISE_LOOP_RADIUS * sinf(angle) *
                                           0.7071067811865475f,
                                       NOISE_LOOP_RADIUS * cosf(angle)};
    } else if constexpr (requires { warp.rotation_rate; }) {
      const float cycle = TWO_PI_F * wrap_t(phase);
      const float cycle_cos = cosf(cycle);
      const float period = 1.0f / params.source.lattice_cell_scale;
      rotation = outer ? outer_rotation : 0.0f;
      prepared.transform.affine = {wrap_t(phase) * warp.translation_x * period,
                                   wrap_t(phase) * warp.translation_y * period,
                                   powf(warp.scale_x, cycle_cos),
                                   powf(warp.scale_y, cycle_cos),
                                   warp.shear * cycle_cos};
    }
    prepared.rotation_cos = cosf(rotation);
    prepared.rotation_sin = sinf(rotation);
    return prepared;
  }

  HS_COLD_MEMBER Frame prepare_frame() {
    HS_PROFILE(fx_prepare_frame);
    if constexpr (HueV == HueMode::NOISE) {
      if (state->hue_noise_lut_scale != params.color.hue_noise_scale ||
          state->hue_noise_lut_phase != hue_noise_phase) {
        Pullback::Color::prepare_hue_noise_lut(
            std::span<int8_t, Pullback::Color::HueNoiseLutView::SIZE>(
                state->hue_noise_lut),
            state->color_noise, params.color.hue_noise_scale, hue_noise_phase);
        state->hue_noise_lut_scale = params.color.hue_noise_scale;
        state->hue_noise_lut_phase = hue_noise_phase;
      }
    }
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            state->hue_rotation_lut),
        palette_cycler.palette());
    const FastNoiseLite *outer_noise = nullptr;
    const FastNoiseLite *source_noise = nullptr;
    if constexpr (HasOuterNoise)
      outer_noise = &state->outer.noise;
    if constexpr (HasSourceNoise)
      source_noise = &state->source.noise;
    return {Derived::ANIMATED_PROJECTION ? projection_conjugate : Quaternion(),
            outer_conjugate,
            {source_primary, source_secondary, source_angle,
             fast_cosf(source_angle), fast_sinf(source_angle)},
            prepare_warp(params.outer_warp, outer_phase, true),
            prepare_warp(params.inner_warp, inner_phase, false),
            outer_noise,
            source_noise,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            params,
            palette_mapping,
            outer_phase,
            inner_phase,
            source_noise_time,
            palette_oscillation_phase};
  }

  HS_COLD_MEMBER void update_palette_chroma() {
    if (palette_chroma == params.color.palette_chroma)
      return;
    palette_chroma = params.color.palette_chroma;
    palette_cycler.set_generated_chroma(palette_chroma);
  }

  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    Runtime &effect = *static_cast<Runtime *>(context);
    if (sequence > 0)
      effect.palette_hue += 159;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, Harmony, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(effect.palette_hue),
        effect.params.color.palette_chroma)};
  }

  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "Pullback::Color::PaletteMapping::CUP",
      "Pullback::Color::PaletteMapping::BELL",
      "Pullback::Color::PaletteMapping::LINEAR",
      "Pullback::Color::PaletteMapping::REVERSE"};

  State *state = nullptr;
  Pullback::Color::PaletteMappingWeights palette_mapping =
      Pullback::Color::PaletteMappingWeights::single(
          params.color.palette_mapping);
  Transition transition;
  uint16_t preset_dwell_remaining = Derived::PRESET_DWELL_FRAMES;
  Timeline timeline;
  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Quaternion projection_walk_previous;
  Quaternion outer_walk_previous;
  Quaternion projection_wander;
  Quaternion outer_wander;
  Quaternion projection_conjugate;
  Quaternion outer_conjugate;
  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float source_noise_time = 0.0f;
  float projection_spin = 0.0f;
  float hue_noise_phase = 0.0f;
  float outer_phase = 0.0f;
  float inner_phase = 0.0f;
  float outer_rotation = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

} // namespace FixedLook
