/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ShadierBall.h
 * @brief Slot-based sphere shader: discrete function, projection, and lens
 *        slots with continuous per-slot params, sequenced as presets.
 */

#include "core/engine/engine.h"

// Unit-test accessor for the accumulators, pipeline stages, and palette walk.
namespace hs_test {
namespace shadierball_tests {
struct ShadierBallWhiteBox;
} // namespace shadierball_tests
} // namespace hs_test

/**
 * @brief Sphere shader assembled from preset-selected slots.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Each preset names a pattern Function, a sphere-to-plane
 * Projection, and a sphere Lens — discrete slot tags dispatched per pixel on
 * frame-constant copies — plus the continuous params those slots consume.
 * Color comes from a PaletteCycler walking the 256 prebaked triadic profiles
 * on a golden-ratio hue step.
 */
template <int W, int H> class ShadierBall : public Effect {
public:
  HS_COLD_MEMBER ShadierBall() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief Registers slot and param controls and starts the palette walk.
   */
  HS_COLD_MEMBER void init() override {
    // register_param captures *ptr as the GUI default, so the live preset has
    // to be loaded before any registration below.
    active = PRESETS[0];
    // Slot tags register animated like every preset-driven dropdown: writes
    // engage the pause, and the future preset choreography owns them.
    register_animated_param("Function", &active.slots.function,
                            FUNCTION_OPTIONS, FUNCTION_EXPORT_OPTIONS,
                            NUM_FUNCTIONS);
    register_animated_param("Projection", &active.slots.projection,
                            PROJECTION_OPTIONS, PROJECTION_EXPORT_OPTIONS,
                            NUM_PROJECTIONS);
    register_animated_param("Lens", &active.slots.lens, LENS_OPTIONS,
                            LENS_EXPORT_OPTIONS, NUM_LENSES);
    Params &params = active.params;
    register_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_param("Pattern Freq", &params.pattern_freq, PATTERN_FREQ_MIN,
                   PATTERN_FREQ_MAX);
    register_param("Wave Spin", &params.wave_spin, WAVE_SPIN_MIN,
                   WAVE_SPIN_MAX);
    register_param("Pole Fade", &params.pole_fade, POLE_FADE_MIN,
                   POLE_FADE_MAX);
    register_param("Lens Mix", &params.lens_mix, LENS_MIX_MIN, LENS_MIX_MAX);

    palette_cycler.init_generated(persistent_arena, next_triadic_palette,
                                  &palette_hue, PALETTE_DWELL_FRAMES,
                                  PALETTE_FADE_FRAMES, ease_in_out_sin,
                                  &anims_paused);
  }

  /**
   * @brief Advances the wave phases and shades every pixel for one frame.
   * @details Wraps the phase accumulators, steps the palette walk, hoists the
   * frame-constant slot tags and wave state, then shades every pixel: lens,
   * stereographic projection, pattern function, pole normalize, palette.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    const Params &P = active.params;
    // Wrap the trig phases to 2pi so fast_sinf/fast_cosf keep precise range
    // reduction.
    phase = fmodf(phase + P.speed, TWO_PI_F);
    wave_angle = fmodf(wave_angle + P.wave_spin, TWO_PI_F);

    palette_cycler.step();

    // Frame-constant hoists; the per-pixel slot switches test these copies, so
    // the branch predictor makes the dispatch nearly free.
    const Slots slots = active.slots;
    const TwinWaveState wave{phase, fast_cosf(wave_angle),
                             fast_sinf(wave_angle)};
    const float pattern_freq = P.pattern_freq;
    const float pole_fade = P.pole_fade;
    const float lens_mix = P.lens_mix;

    auto shader = [&](const Vector &v) -> Color4 {
      const Sample s =
          sample_pipeline(v, slots, wave, pattern_freq, pole_fade, lens_mix);
      return palette_cycler.palette().get(s.value);
    };
    {
      HS_PROFILE(sdb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
  }

private:
  // Test seam: reaches the accumulators, pipeline stages, and palette walk.
  friend struct ::hs_test::shadierball_tests::ShadierBallWhiteBox;

  /** @brief Pattern-field slot tag; one body per enumerator. */
  enum class Function : uint8_t { TWIN_WAVE };
  /** @brief Sphere-to-plane projection slot tag. */
  enum class Projection : uint8_t { EQUIRECTANGULAR, STEREOGRAPHIC, GNOMONIC };
  /** @brief Sphere-lens slot tag applied before projection. */
  enum class Lens : uint8_t { NONE, GLITCH };

  /** @brief Slot tags one preset selects. */
  struct Slots {
    Function function;     /**< Pattern field evaluated per pixel. */
    Projection projection; /**< Sphere-to-plane map feeding the pattern. */
    Lens lens;             /**< Sphere lens ahead of the projection. */
  };

  /** @brief Continuous params the slots consume, one snapshot per preset. */
  struct Params {
    float speed;        /**< Wave phase advance, rad/frame. */
    float pattern_freq; /**< Spatial frequency of the pattern field. */
    float wave_spin;    /**< Second-wave rotation rate, rad/frame. */
    float pole_fade;    /**< Pole attenuation radius. */
    float lens_mix;     /**< Lens blend in [0, 1]; 0 and 1 skip a projection. */
  };

  /** @brief One sequencable look: slot tags plus their params. */
  struct Preset {
    Slots slots;
    Params params;
  };

  /** @brief Frame-constant state of the twin-wave function. */
  struct TwinWaveState {
    float phase;    /**< Shared travel phase of both waves. */
    float wave_cos; /**< Cosine of the second wave's rotation angle. */
    float wave_sin; /**< Sine of the second wave's rotation angle. */
  };

  /** @brief Inter-stage record the function stage hands the palette stage. */
  struct Sample {
    float value; /**< Pole-normalized pattern value in [0, 1]. */
  };

  /** @brief Full per-pixel pipeline: lens, projection, function, normalize. */
  static Sample sample_pipeline(const Vector &v, const Slots &slots,
                                const TwinWaveState &wave, float pattern_freq,
                                float pole_fade, float lens_mix) {
    const Complex z = project(v, slots.projection, slots.lens, lens_mix);
    const float r_sq = z.re * z.re + z.im * z.im;
    const float pattern = sample_function(
        slots.function, stereo_pattern_args(z, pattern_freq), wave);
    return {pole_normalize_pattern(pattern, r_sq, pole_fade)};
  }

  /** @brief Lens dispatch and projection to the pattern plane. */
  static Complex project(const Vector &v, Projection projection, Lens lens,
                         float lens_mix) {
    switch (lens) {
    case Lens::NONE:
      return project_point(v, projection);
    case Lens::GLITCH: {
      if (lens_mix == 0.0f)
        return project_point(v, projection);
      const Complex lensed = project_point(glitch_lens(v), projection);
      if (lens_mix == 1.0f)
        return lensed;
      const Complex direct = project_point(v, projection);
      return {hs::lerp(direct.re, lensed.re, lens_mix),
              hs::lerp(direct.im, lensed.im, lens_mix)};
    }
    }
    __builtin_unreachable();
  }

  /** @brief Projection-slot dispatch: unit direction to plane coordinates. */
  static Complex project_point(const Vector &v, Projection projection) {
    switch (projection) {
    case Projection::EQUIRECTANGULAR:
      return equirectangular(v);
    case Projection::STEREOGRAPHIC:
      return stereo(v);
    case Projection::GNOMONIC:
      return gnomonic(v);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Equirectangular map of a unit direction.
   * @param v Unit direction vector on the sphere.
   * @return re = azimuth in (-pi, pi], im = latitude in [-pi/2, pi/2].
   */
  static Complex equirectangular(const Vector &v) {
    return {fast_atan2(v.z, v.x), 0.5f * PI_F - fast_acos(v.y)};
  }

  /**
   * @brief Gnomonic map through the sphere center onto the y = 1 plane.
   * @param v Unit direction vector on the sphere.
   * @return Plane coordinates; antipodal directions map to the same point.
   */
  static Complex gnomonic(const Vector &v) {
    // Floor |y| so the equator ring divides by the epsilon instead of blowing
    // up to inf/NaN; the resulting ~1e3 coordinates stay inside the
    // pattern-arg clamp.
    float y = v.y;
    if (std::fabs(y) < GNOMONIC_AXIS_EPS)
      y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
    return {v.x / y, v.z / y};
  }

  /** @brief Function-slot dispatch on a frame-constant tag. */
  static float sample_function(Function function, const Complex &p,
                               const TwinWaveState &wave) {
    switch (function) {
    case Function::TWIN_WAVE:
      return twin_wave(p, wave);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Two traveling plane waves, the second rotated by the spin angle.
   * @param p Frequency-scaled, bounded pattern arguments.
   * @param wave Frame-constant phase and rotation of the pair.
   * @return Mean of the two waves in [-1, 1].
   */
  static float twin_wave(const Complex &p, const TwinWaveState &wave) {
    const float rotated = p.re * wave.wave_cos + p.im * wave.wave_sin;
    return 0.5f *
           (fast_sinf(p.re + wave.phase) + fast_sinf(rotated + wave.phase));
  }

  /**
   * @brief Trig-free glitch lens on a sphere direction.
   * @param v Unit direction vector on the sphere.
   * @return Direction after latitude doubling and azimuth tripling; returns
   * the up vector near the lens axis.
   */
  // Mirror of ShaderBall's apply_glitch_lens (effects/ShaderBall.h).
  static Vector glitch_lens(const Vector &v) {
    const float x2 = v.x * v.x;
    const float z2 = v.z * v.z;
    const float radius2 = x2 + z2;
    constexpr float MIN_AXIS_RADIUS2 = 1e-6f;
    if (radius2 < MIN_AXIS_RADIUS2)
      return Vector(0.0f, 1.0f, 0.0f);

    const float inverse_radius2 = 1.0f / radius2;
    const float double_y = 2.0f * v.y;
    return Vector(double_y * v.x * (4.0f * x2 * inverse_radius2 - 3.0f),
                  2.0f * v.y * v.y - 1.0f,
                  double_y * v.z * (3.0f - 4.0f * z2 * inverse_radius2));
  }

  /** @brief PaletteCycler provider: the triadic profile walking the prebaked
   *  hue wheel; sequence 0 returns the starting hue unchanged. */
  HS_COLD_MEMBER static void next_triadic_palette(void *context,
                                                  uint32_t sequence,
                                                  GenerativePalette &out) {
    uint32_t &hue = *static_cast<uint32_t *>(context);
    if (sequence > 0)
      hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(hue))};
  }

  /** @brief Zero dwell: the palette shifts continuously, fade after fade. */
  static constexpr int PALETTE_DWELL_FRAMES = 0;
  /** @brief Frames each palette fade spans. */
  static constexpr int PALETTE_FADE_FRAMES = 600;
  /** @brief Hue-wheel steps per palette. Odd, so the walk visits all 256
   *  prebaked hues before repeating; 159/256 approximates the golden-ratio
   *  conjugate. */
  static constexpr uint32_t HUE_STEP = 159;

  /** @brief |v.y| floor for the gnomonic division. */
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;

  static constexpr const char *FUNCTION_OPTIONS[] = {"Twin Wave"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE"};
  static constexpr int NUM_FUNCTIONS =
      static_cast<int>(std::size(FUNCTION_OPTIONS));
  static constexpr const char *PROJECTION_OPTIONS[] = {
      "Equirectangular", "Stereographic", "Gnomonic"};
  static constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
      "Projection::EQUIRECTANGULAR", "Projection::STEREOGRAPHIC",
      "Projection::GNOMONIC"};
  static constexpr int NUM_PROJECTIONS =
      static_cast<int>(std::size(PROJECTION_OPTIONS));
  static constexpr const char *LENS_OPTIONS[] = {"None", "Glitch"};
  static constexpr const char *LENS_EXPORT_OPTIONS[] = {"Lens::NONE",
                                                        "Lens::GLITCH"};
  static constexpr int NUM_LENSES = static_cast<int>(std::size(LENS_OPTIONS));

  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 0.5f;
  static constexpr float PATTERN_FREQ_MIN = 1.0f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float WAVE_SPIN_MIN = 0.0f, WAVE_SPIN_MAX = 0.05f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float LENS_MIX_MIN = 0.0f, LENS_MIX_MAX = 1.0f;

  /** @brief True iff every param of @p p lies within its registered range. */
  static constexpr bool preset_in_ranges(const Preset &p) {
    return p.params.speed >= SPEED_MIN && p.params.speed <= SPEED_MAX &&
           p.params.pattern_freq >= PATTERN_FREQ_MIN &&
           p.params.pattern_freq <= PATTERN_FREQ_MAX &&
           p.params.wave_spin >= WAVE_SPIN_MIN &&
           p.params.wave_spin <= WAVE_SPIN_MAX &&
           p.params.pole_fade >= POLE_FADE_MIN &&
           p.params.pole_fade <= POLE_FADE_MAX &&
           p.params.lens_mix >= LENS_MIX_MIN &&
           p.params.lens_mix <= LENS_MIX_MAX;
  }

  static constexpr Preset PRESETS[] = {
      // Rotating twin-wave interference through the full glitch lens.
      {{Function::TWIN_WAVE, Projection::STEREOGRAPHIC, Lens::GLITCH},
       {0.05f, 6.0f, 0.006f, 2.0f, 1.0f}},
  };
  static_assert(
      [] {
        for (const Preset &p : PRESETS)
          if (!preset_in_ranges(p))
            return false;
        return true;
      }(),
      "a ShadierBall preset drives a param outside its registered slider "
      "range; widen the range to accommodate the preset");

  Preset active = PRESETS[0]; /**< Live slots and params; the GUI edits these
                                 directly. */
  float phase = 0.0f;         /**< Wrapped to [0, 2pi): shared wave travel. */
  float wave_angle = 0.0f;    /**< Wrapped to [0, 2pi): second-wave rotation. */
  uint32_t palette_hue = 0;   /**< Hue-wheel index of the current triadic. */
  PaletteCycler palette_cycler;

  // init() buys the palette cycler's display LUT and morph slots from the
  // persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes();
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "ShadierBall persistent footprint exceeds the default "
                "partition");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShadierBall)
