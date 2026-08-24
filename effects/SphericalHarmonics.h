/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file SphericalHarmonics.h
 * @brief Real spherical harmonics painted on the unit sphere, continuously
 *        morphing between modes.
 */

#include "core/engine/engine.h"
#include "core/math/spherical_harmonics.h"
#include "core/render/pullback.h"
#include <array>
#include <cmath>
#include <string_view>
#include <utility>

// Unit-test accessor reaching the private morph-chain state and the baked
// palette the shader reads.
namespace hs_test {
namespace effects_tests {
struct SphericalHarmonicsWhiteBox;
} // namespace effects_tests
} // namespace hs_test

/**
 * @brief Visualizer that paints real spherical harmonics on the unit sphere.
 * @tparam W Display width in pixels.
 * @tparam H Display height in pixels.
 * @details Continuously morphs between randomly chosen modes.
 */
template <int W, int H> class SphericalHarmonics : public Effect {
public:
  static constexpr std::array<std::string_view, 24> PRESET_IDS{
      "sh-l2-m0",  "sh-l1-m-1", "sh-l1-m0",  "sh-l1-m1",  "sh-l2-m-2",
      "sh-l2-m-1", "sh-l2-m1",  "sh-l2-m2",  "sh-l3-m-3", "sh-l3-m-2",
      "sh-l3-m-1", "sh-l3-m0",  "sh-l3-m1",  "sh-l3-m2",  "sh-l3-m3",
      "sh-l4-m-4", "sh-l4-m-3", "sh-l4-m-2", "sh-l4-m-1", "sh-l4-m0",
      "sh-l4-m1",  "sh-l4-m2",  "sh-l4-m3",  "sh-l4-m4"};

  /**
   * @brief Field sampler that evaluates the (blended) harmonic at a world point.
   */
  struct HarmonicField {
    int l1, m1;
    int l2, m2;
    float blend;
    RotationMatrix orientation_conj; /**< World->local rotation (conjugate of
                                        the shape orientation). */
    float N1, N2; /**< Per-mode harmonic scales, precomputed once per shape. */

    /**
     * @brief Construct a field for blending mode (l1, m1) into (l2, m2).
     * @param l1 Degree of the first harmonic.
     * @param m1 Order of the first harmonic.
     * @param l2 Degree of the second harmonic.
     * @param m2 Order of the second harmonic.
     * @param blend Morph fraction in [0, 1] from the first toward the second.
     * @param q Orientation quaternion of the shape.
     */
    HarmonicField(int l1, int m1, int l2, int m2, float blend, Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2), blend(blend),
          orientation_conj(q.conjugate()), N1(SHMath::harmonic_scale(l1, m1)),
          N2(SHMath::harmonic_scale(l2, m2)) {}

    /**
     * @brief Sample the (possibly blended) harmonic at world point p.
     * @param p World-space sample point.
     * @return Signed field value; zero-crossings are the lobe boundaries.
     * @details Rotates p into the shape's local frame and evaluates both modes
     * there.
     */
    float sample(const Vector &p) const {
      Vector local = orientation_conj.apply(p);

      // The shape spins about an arbitrary axis, so the local frame varies
      // across a screen row even though the WORLD latitude is row-constant.
      // Hoisting any part of the harmonic to a per-row precompute is invalid.
      float val = SHMath::spherical_harmonic(l1, m1, local, N1);
      if (blend > 0.001f) {
        float val2 = SHMath::spherical_harmonic(l2, m2, local, N2);
        val += (val2 - val) * blend;
      }
      return val;
    }
  };

  /**
   * @brief Construct the visualizer with the display dimensions.
   */
  HS_COLD_MEMBER SphericalHarmonics() : Effect(W, H, {.strobe = true}) {}

  /**
   * @brief One-time setup of params, palette, shape, spin, and first morph.
   * @details Registers params, bakes the palette, seeds the shape and the
   * continuous spin, and kicks off the first morph.
   */
  void init() override {
    configure_presets(PRESET_IDS.size());
    register_param("Amplitude", &params.amplitude, 0.1f, 10.0f);

    baked_palette.bake(persistent_arena, Palettes::RICH_SUNSET);

    current_idx = SEED_MODE_IDX;

    Vector axis = Vector(0.5f, 1.0f, 0.2f).normalized();
    timeline.add(0, Animation::Rotation<W>(orientation, axis,
                                           2 * PI_F * SPIN_TURNS, SPIN_FRAMES,
                                           ease_linear, true));

    start_morph();
  }

  /**
   * @brief Render one frame of the morphing harmonic.
   * @details Decodes the current and target modes, builds the field for this
   * frame's morph state, and shades the sphere with the harmonic-coloring
   * shader.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sh_timeline_step);
      timeline.step(canvas);
    }

    auto [l1, m1] = SHMath::decode_lm(current_idx);
    auto [l2, m2] = SHMath::decode_lm(next_idx);
    const typename RenderPipeline::Frame frame = RenderPipeline::prepare(
        {{l1, m1, l2, m2, morph_alpha, orientation.get()},
         {&baked_palette, params.amplitude}});

    {
      HS_PROFILE(sh_rasterize);
      Scan::Shader::draw<W, H>(canvas, [&frame](const Vector &view) {
        return RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
      });
    }
  }

private:
  static constexpr int mode_index_for_preset(size_t preset_index) {
    if (preset_index == 0)
      return SEED_MODE_IDX;
    return static_cast<int>(preset_index) +
           (preset_index < static_cast<size_t>(SEED_MODE_IDX) ? 0 : 1);
  }

  static constexpr size_t preset_index_for_mode(int mode_index) {
    if (mode_index == SEED_MODE_IDX)
      return 0;
    return static_cast<size_t>(mode_index - (mode_index > SEED_MODE_IDX));
  }

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    ++morph_generation;
    current_idx = mode_index_for_preset(change.to);
    next_idx = current_idx;
    morph_alpha = 0.0f;
    return true;
  }

  friend struct ::hs_test::effects_tests::SphericalHarmonicsWhiteBox;

  /**
   * @brief Ambient-occlusion shaping constants for the draw_frame shader.
   * @details AO_FALLOFF is the field magnitude at which shading saturates,
   * AO_AMBIENT is the ambient floor, and AO_RANGE is the range above it
   * (AMBIENT + RANGE = 1.0 → full brightness at saturation).
   */
  static constexpr float AO_FALLOFF = 0.4f;
  static constexpr float AO_AMBIENT = 0.15f;
  static constexpr float AO_RANGE = 0.85f;

  struct HarmonicSourceState {
    int l1;
    int m1;
    int l2;
    int m2;
    float blend;
    Quaternion orientation;
  };

  struct HarmonicMaterialState {
    const BakedPalette *palette;
    float amplitude;
  };

  struct HarmonicFrameState {
    HarmonicSourceState source;
    HarmonicMaterialState material;
  };

  struct HarmonicBinding {
    using FrameState = HarmonicFrameState;
    using Instrumentation = Pullback::NoInstrumentation;
  };

  struct HarmonicSourcePolicy : Pullback::ApproximationDefaults {
    using Prepared = HarmonicField;

    HS_FLASH_INLINE static Prepared prepare(const HarmonicFrameState &frame) {
      const HarmonicSourceState &source = frame.source;
      return {source.l1, source.m1,    source.l2,
              source.m2, source.blend, source.orientation};
    }

    __attribute__((always_inline)) static float
    sample(const Pullback::SphereSample &input, const HarmonicFrameState &,
           const Prepared &field) {
      return field.sample(input.dir);
    }
  };

  __attribute__((always_inline)) static Color4
  colorize_harmonic(float val, const HarmonicMaterialState &material) {
    const float abs_val = std::abs(val);
    Color4 pos =
        material.palette->get(std::min(1.0f, abs_val * material.amplitude));

    constexpr float NEG_LOBE_GREEN_SCALE = 0.8f;
    const Color4 neg = Color4(
        Pixel(pos.color.b,
              static_cast<uint16_t>(pos.color.g * NEG_LOBE_GREEN_SCALE + 0.5f),
              pos.color.r),
        pos.alpha);

    constexpr float TRANSITION = 0.03f;
    const float blend_t =
        quintic_kernel(hs::clamp(val / TRANSITION * 0.5f + 0.5f, 0.0f, 1.0f));
    Color4 base = pos.lerp(neg, 1.0f - blend_t);

    const float shadow =
        hs::clamp((abs_val * material.amplitude) / AO_FALLOFF, 0.0f, 1.0f);
    base.color = base.color * (AO_AMBIENT + AO_RANGE * shadow);
    return base;
  }

  struct HarmonicColorPolicy : Pullback::ApproximationDefaults {
    __attribute__((always_inline)) static Color4
    apply(const Pullback::FieldSample &sample,
          const HarmonicFrameState &frame) {
      return colorize_harmonic(sample.value * 2.0f - 1.0f, frame.material);
    }
  };

  using RenderPipeline =
      Pullback::Pipeline<HarmonicBinding,
                         Pullback::Stage::SampleSphere<HarmonicSourcePolicy>,
                         Pullback::Stage::Colorize<HarmonicColorPolicy>>;

  /**
   * @brief Choose the next harmonic and animate the morph toward it.
   * @details On completion commits it as current and recurses, yielding an
   * endless morph chain. Pausable: the chain is this effect's look
   * choreography, so "Pause Animation" holds the current blend.
   */
  HS_COLD_MEMBER void start_morph() {
    const uint32_t generation = morph_generation;
#ifdef HS_PROFILE_ORDERED_CYCLE
    next_idx = (current_idx % MAX_MODE_IDX) + 1;
#else
    // Re-roll on a match with current_idx: blending a basis function into
    // itself would freeze the sphere for the whole transition.
    do {
      next_idx = hs::rand_int(1, MAX_MODE_IDX + 1);
    } while (next_idx == current_idx);
#endif
#if defined(HS_PROFILE_ENABLE)
    hs::log("Mode: %d/%d", next_idx, MAX_MODE_IDX);
#endif

    timeline.add_pausable(
        0,
        Animation::Transition(morph_alpha, 1.0f, 64, ease_linear, false, false)
            .then([this, generation]() {
              if (generation != morph_generation) {
                morph_alpha = 0.0f;
                start_morph();
                return;
              }
              current_idx = next_idx;
              morph_alpha = 0.0f;
              HS_CHECK(synchronizePreset(preset_index_for_mode(current_idx)));
              start_morph();
            }),
        &anims_paused);
  }

  Orientation<> orientation;  /**< Current sphere orientation. */
  Timeline timeline;          /**< Drives spin and morph animations. */
  BakedPalette baked_palette; /**< Precomputed color LUT for the shader. */

  // init() bakes one palette LUT into the persistent arena. Effect keeps the
  // default arena split, so the total must fit the device persistent partition.
  static constexpr size_t FOOTPRINT_BYTES =
      BakedPalette::required_arena_bytes();
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "SphericalHarmonics persistent footprint exceeds the default partition");

  // Highest harmonic degree the morph visits. Modes are flat-indexed
  // idx = l*l + l + m, so degrees [0, MAX_DEGREE] occupy idx [0, MAX_MODE_IDX].
  static constexpr int MAX_DEGREE = 4;
  static_assert(MAX_DEGREE <= 5,
                "SampleSphere requires normalized harmonics bounded by one");
  // normalization() divides factorials of up to 2*MAX_DEGREE; float holds exact
  // integers only up to 2^24, first exceeded at 11!.
  static_assert(2 * MAX_DEGREE <= 10,
                "MAX_DEGREE too large: factorial(2*MAX_DEGREE) exceeds float "
                "precision in normalization()");
  // Top flat index over those degrees: idx peaks at l = MAX_DEGREE, m = +MAX_DEGREE.
  static constexpr int MAX_MODE_IDX =
      (MAX_DEGREE + 1) * (MAX_DEGREE + 1) - 1; // 24
  static_assert(PRESET_IDS.size() == MAX_MODE_IDX);
  // Initial mode (l=2, m=0); the constant mode (idx 0) is excluded from the roll.
  static constexpr int SEED_MODE_IDX = 6;
  static_assert(SEED_MODE_IDX > 0 && SEED_MODE_IDX <= MAX_MODE_IDX,
                "seed mode must be a valid, non-constant harmonic index");
  // Looped continuous spin: SPIN_TURNS revolutions over SPIN_FRAMES frames,
  // i.e. 0.01 turn (3.6 degrees) per frame.
  static constexpr float SPIN_TURNS = 100.0f;
  static constexpr int SPIN_FRAMES = 10000;

  int current_idx = 0; /**< Flat index of the currently displayed mode. */
  int next_idx = 0;    /**< Flat index of the mode being morphed toward. */
  float morph_alpha =
      0.0f; /**< Morph progress in [0, 1] from current to next. */
  uint32_t morph_generation = 0;

  /**
   * @brief User-tunable parameters for the visualizer.
   */
  struct Params {
    float amplitude = 3.2f; /**< Field-value gain applied before coloring. */
  } params;
};

#include "core/control/registry.h"
REGISTER_EFFECT(SphericalHarmonics)
