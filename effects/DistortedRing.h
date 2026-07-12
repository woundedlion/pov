/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

namespace hs_test {
namespace effects_tests {
struct DistortedRingWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Concentric rings warped by an animated sinusoidal radial displacement.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Rings are drawn in ring-space and warped by a sum of sine waves
 * with different frequencies, each with its own drifting azimuth phase and
 * independently animated amplitude. Orientation random-walks over time and
 * each ring is shaded from a circular split-complementary palette, viewed
 * through a GUI-selectable palette-modifier composition.
 */
template <int W, int H> class DistortedRing : public Effect {
public:
  /**
   * @brief Builds the effect with its palette and ring normal.
   * @details Configures a circular split-complementary palette and an X-axis
   * ring normal.
   */
  HS_COLD_MEMBER DistortedRing()
      : Effect(W, H, {.strobe = true, .full_frame = decltype(filters)::any_crosses_segments}),
        ring_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY,
                    BrightnessProfile::FLAT),
        normal(X_AXIS) {}

  /**
   * @brief Registers params and builds the timeline.
   * @details Adds the ring sprite, orientation random walk, and one amplitude
   * mutation per distortion wave to the timeline.
   */
  void init() override {
    // One pixel of azimuth in ring-space.
    const float px = 2.0f * PI_F / W;
    params.thickness = 4.0f * px;

    ring_baked.bake(persistent_arena, ring_palette);

    warp_palette.bind(&ring_baked, &warp_mod);
    drift_palette.bind(&ring_baked, &drift_mod);
    spin_palette.bind(&ring_baked, &spin_shade);
    wobble_palette.bind(&ring_baked, &wobble_shade);
    sparkle_palette.bind(&ring_baked, &sparkle_shade);
    pulse_palette.bind(&ring_baked, &pulse_shade);
    grain_palette.bind(&ring_baked, &grain_shade);
    sheen_palette.bind(&ring_baked, &sheen_shade);

    register_param("Alpha", &params.alpha, 0.0f, 1.0f);
    register_param("MaxAmplitude", &params.max_amplitude, 0.0f, 2.0f);
    register_param("Thickness", &params.thickness, 2.0f * px, 12.0f * px);
    register_param("Rings", &params.num_rings, 1.0f, 10.0f);
    register_param("Palette Mod", &params.palette_mod, PALETTE_MODS,
                   static_cast<int>(sizeof(PALETTE_MODS) /
                                    sizeof(PALETTE_MODS[0])));
    register_param("Mod Depth", &params.mod_depth, 0.0f, 1.0f);
    register_param("Show Bounding", &params.debug_bb);

    timeline.add(0, Animation::Sprite(
                        [this](Canvas &canvas, float opacity) {
                          this->draw_fn(canvas, opacity);
                        },
                        -1, 48, ease_linear, 0, ease_linear));

    timeline.add(0, Animation::RandomWalk<W>(orientation, normal, noise));

    for (int i = 0; i < NUM_WAVES; ++i) {
      timeline.add(0, Animation::Mutation(
                          wave_amp[i],
                          [this, i](float t) {
                            return sin_wave(-params.max_amplitude,
                                            params.max_amplitude, 1.0f,
                                            WAVES[i].env_phase)(t);
                          },
                          WAVES[i].env_frames, ease_linear, true));
    }
  }

  /**
   * @brief Advances and renders the timeline for one frame.
   */
  void draw_frame() override {
    step_mod_drivers();
    for (int i = 0; i < NUM_WAVES; ++i)
      wave_phase[i] = wrap_t(wave_phase[i] + WAVES[i].drift);
    Canvas canvas(*this);
    timeline.step(canvas);
  }

  /**
   * @brief Draws all rings for this frame.
   * @param canvas Render target for the ring fragments.
   * @param opacity Sprite's animated fade in [0, 1], multiplied into each
   * fragment's alpha.
   * @details Ring radii are evenly spaced and the displacement is the sum of
   * the WAVES table's sines, baked into shift_lut once per frame and sampled
   * by azimuth with linear interpolation.
   */
  void draw_fn(Canvas &canvas, float opacity) {
    int n_rings = static_cast<int>(params.num_rings);

    for (int x = 0; x <= W; ++x) {
      float t = static_cast<float>(x) / W;
      float s = 0.0f;
      for (int i = 0; i < NUM_WAVES; ++i)
        s += wave_amp[i] * WAVES[i].weight *
             sinf(2.0f * PI_F * (WAVES[i].freq * t + wave_phase[i]));
      shift_lut[x] = s;
    }

    // The slope out-param feeds the rasterizer's slope compensation.
    auto shift_fn = [this](float t, float &slope) {
      float x = wrap_t(t) * W;
      int i = static_cast<int>(x);
      float d = shift_lut[i + 1] - shift_lut[i];
      slope = d * W;
      return shift_lut[i] + (x - i) * d;
    };

    // True upper bound on |shift_fn|: an underestimate silently culls arcs.
    float max_shift = 0.0f;
    for (int i = 0; i < NUM_WAVES; ++i)
      max_shift += std::abs(wave_amp[i]) * WAVES[i].weight;

    Basis basis = make_basis(orientation.get(), normal);
    for (int i = 0; i < n_rings; ++i) {
      float radius = 2.0f / (n_rings + 1) * (i + 1);
      auto fragment_shader = [&](const Vector &, Fragment &f) {
        // f.v0 = normalized azimuth (0..1), f.v1 = distance from center line,
        // f.size = thickness.
        f.color = sample_palette(f.v0);

        float norm_dist = hs::clamp(f.v1 / f.size, 0.0f, 1.0f);
        float falloff = quintic_kernel(1.0f - norm_dist);

        f.color.alpha = f.color.alpha * opacity * params.alpha * falloff;
      };

      Scan::DistortedRing::draw<W, H>(
          filters, canvas, basis, radius, params.thickness, shift_fn,
          max_shift, fragment_shader, 0.0f, params.debug_bb);
    }
  }

private:
  friend struct ::hs_test::effects_tests::DistortedRingWhiteBox;

  /**
   * @brief Samples the ring palette through the selected modifier composition.
   * @param t Palette coordinate (normalized ring azimuth).
   * @return The sample of the composition selected by params.palette_mod.
   * @details Case order matches PALETTE_MODS.
   */
  Color4 sample_palette(float t) const {
    switch (static_cast<int>(params.palette_mod)) {
    case 1: return warp_palette.get(t);
    case 2: return drift_palette.get(t);
    case 3: return spin_palette.get(t);
    case 4: return wobble_palette.get(t);
    case 5: return sparkle_palette.get(t);
    case 6: return pulse_palette.get(t);
    case 7: return grain_palette.get(t);
    case 8: return sheen_palette.get(t);
    default: return ring_baked.get(t);
    }
  }

  /**
   * @brief Advances the palette-modifier drivers one frame and applies the
   * Mod Depth slider to every modifier's strength knob.
   */
  void step_mod_drivers() {
    mod_time += MOD_TIME_STEP;
    mod_phase += MOD_PHASE_STEP;
    if (mod_phase >= 2.0f * PI_F)
      mod_phase -= 2.0f * PI_F;
    spin_amount += SPIN_RATE * params.mod_depth;
    if (spin_amount >= 1.0f)
      spin_amount -= 1.0f;

    const float d = params.mod_depth;
    warp_mod.amplitude = WARP_AMP * d;
    drift_mod.amplitude = DRIFT_AMP * d;
    wobble_shade.depth = WOBBLE_DEPTH * d;
    sparkle_shade.threshold = 1.0f - SPARKLE_SPAN * d;
    pulse_shade.depth = PULSE_DEPTH * d;
    grain_shade.amplitude = GRAIN_AMP * d;
    sheen_shade.weight = SHEEN_WEIGHT * d;
  }

  FastNoiseLite noise;
  Timeline timeline;
  Pipeline<W, H> filters;

  /**
   * @brief Slider-backed parameters.
   * @details Defaults are pre-registration starting values.
   */
  struct Params {
    float alpha = 0.3f;         /**< Overall ring opacity multiplier in [0, 1]. */
    float max_amplitude = 0.3f; /**< Peak radial displacement; the animated wave sweeps [-max_amplitude, +max_amplitude]. */
    float thickness = 1.0f;     /**< Ring stroke width, in radians of azimuth; init() reseeds it to 4 px scaled by resolution. */
    float num_rings = 1.0f;     /**< Number of evenly spaced concentric rings (truncated to int when drawn). */
    float palette_mod = 0.0f;   /**< Selected palette-modifier option index (see PALETTE_MODS). */
    float mod_depth = 0.75f;    /**< Strength of the selected palette modifier in [0, 1]. */
    bool debug_bb = false;      /**< When true, draws each fragment's bounding box for debugging. */
  } params;

  /** @brief One component of the radial displacement sum. */
  struct Wave {
    float freq;      /**< Cycles around the ring; integer keeps the wrap seamless. */
    float weight;    /**< Share of MaxAmplitude; weights sum to 1 so the sum's bound holds. */
    float drift;     /**< Azimuth phase advance per frame (cycles). */
    float env_phase; /**< Amplitude-envelope start offset (cycles). */
    int env_frames;  /**< Amplitude-envelope period (frames). */
  };

  static constexpr Wave WAVES[] = {
      {3.0f, 0.5f, 0.008f, 0.0f, 32},
      {5.0f, 0.3f, -0.013f, 0.33f, 44},
      {7.0f, 0.2f, 0.021f, 0.67f, 60},
  };
  static constexpr int NUM_WAVES =
      static_cast<int>(sizeof(WAVES) / sizeof(WAVES[0]));

  float wave_amp[NUM_WAVES] = {};   /**< Animated amplitude per wave (radians). */
  float wave_phase[NUM_WAVES] = {}; /**< Azimuth phase per wave (cycles, [0,1)). */
  float shift_lut[W + 1] = {};      /**< Per-frame displacement sum at pixel resolution; entry W repeats entry 0 for seamless lerp. */

  GenerativePalette ring_palette;
  BakedPalette ring_baked;
  Vector normal;
  Orientation<> orientation;

  /**
   * @brief Dropdown labels for the palette-modifier selection, in
   * sample_palette() case order.
   */
  static constexpr const char *PALETTE_MODS[] = {
      "None",        "Noise Warp", "Drift",           "Hue Spin",
      "Hue Wobble",  "Sparkle",    "Chroma Pulse",    "Lightness Grain",
      "Iridescent"};

  static constexpr float MOD_TIME_STEP = 0.05f;  /**< Noise-time advance per frame. */
  static constexpr float MOD_PHASE_STEP = 0.04f; /**< Phase advance per frame (radians). */

  // Per-modifier strength at Mod Depth = 1; step_mod_drivers() scales each by
  // the slider every frame.
  static constexpr float SPIN_RATE = 0.006f;   /**< Hue-spin turns per frame. */
  static constexpr float WARP_AMP = 0.5f;      /**< Noise-warp peak displacement. */
  static constexpr float DRIFT_AMP = 0.5f;     /**< Drift-walk peak offset. */
  static constexpr float WOBBLE_DEPTH = 0.35f; /**< Hue-wobble peak rotation (turns). */
  static constexpr float SPARKLE_SPAN = 0.45f; /**< Sparkle threshold drop below 1. */
  static constexpr float PULSE_DEPTH = 1.0f;   /**< Chroma-pulse swing. */
  static constexpr float GRAIN_AMP = 0.7f;     /**< Lightness-grain gain swing. */
  static constexpr float SHEEN_WEIGHT = 0.6f;  /**< Iridescent overlay strength. */

  // Drivers must precede the modifiers below, whose constructors take their
  // addresses. Strength knobs are placeholders here; step_mod_drivers()
  // overwrites them from Mod Depth before the first frame draws.
  float mod_time = 0.0f;
  float mod_phase = 0.0f;
  float spin_amount = 0.0f;

  NoiseWarpModifier warp_mod{&mod_time, 5.0f, 0.0f};
  DriftModifier drift_mod{&mod_time, 0.25f, 0.0f, 1u};
  HueSpinShade spin_shade{&spin_amount};
  HueWobbleShade wobble_shade{&mod_phase, 2.0f, 0.0f};
  SparkleShade sparkle_shade{&mod_time, 32.0f, 0.99f, 2u};
  ChromaPulseShade pulse_shade{&mod_phase, 0.0f};
  LightnessGrainShade grain_shade{&mod_time, 20.0f, 0.0f, 3u};
  IridescentShade sheen_shade{&mod_phase, 6.0f, 0.0f};

  StaticPalette<BakedPalette, Coords<NoiseWarpModifier>> warp_palette;
  StaticPalette<BakedPalette, Coords<DriftModifier>> drift_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<HueSpinShade>> spin_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<HueWobbleShade>> wobble_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<SparkleShade>> sparkle_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<ChromaPulseShade>> pulse_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<LightnessGrainShade>>
      grain_palette;
  StaticPalette<BakedPalette, Coords<>, Colors<IridescentShade>> sheen_palette;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(DistortedRing)
