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

#include "core/platform/effect_identity.h"
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "core/math/spherical_harmonics.h"
#include "core/render/pullback.h"
#include "core/render/spacetime/field.h"
#include "core/render/spacetime.h"
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
 * @details Continuously morphs through every nonconstant mode in order.
 */
template <int W, int H> class SphericalHarmonics : public Effect {
public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::SphericalHarmonics::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /// Preset identities: one per (l, m) mode.
  static constexpr std::array<std::string_view, 24> PRESET_IDS{
      "sh-l2-m0",  "sh-l1-m-1", "sh-l1-m0",  "sh-l1-m1",  "sh-l2-m-2",
      "sh-l2-m-1", "sh-l2-m1",  "sh-l2-m2",  "sh-l3-m-3", "sh-l3-m-2",
      "sh-l3-m-1", "sh-l3-m0",  "sh-l3-m1",  "sh-l3-m2",  "sh-l3-m3",
      "sh-l4-m-4", "sh-l4-m-3", "sh-l4-m-2", "sh-l4-m-1", "sh-l4-m0",
      "sh-l4-m1",  "sh-l4-m2",  "sh-l4-m3",  "sh-l4-m4"};
  /** Frames for each harmonic morph. */
  static constexpr int MORPH_FRAMES = 64;
  /** Frames covering every mode and returning to the starting mode. */
  static constexpr size_t FULL_CYCLE_FRAMES = PRESET_IDS.size() * MORPH_FRAMES;
  /** @brief Number of selectable harmonic modes.
   * @return Authored preset count. */
  static consteval size_t authored_preset_count() { return PRESET_IDS.size(); }
  /** @brief Frames covering the complete harmonic mode cycle.
   * @return Full cycle frame budget. */
  static consteval size_t cycle_frames() { return FULL_CYCLE_FRAMES; }

  /**
   * @brief Field sampler that evaluates the (blended) harmonic at a world point.
   */
  struct HarmonicField {
    /// Degree of the first harmonic.
    int l1, m1; ///< Order of the first harmonic.
    /// Degree of the second harmonic.
    int l2, m2;  ///< Order of the second harmonic.
    float blend; ///< Morph fraction in [0, 1] toward (l2, m2).
    math::RotationMatrix
        orientation_conj;  /**< World->local rotation (conjugate of
                                        the shape orientation). */
    float normalization_1; ///< Normalization of the first harmonic.
    float normalization_2; ///< Normalization of the second harmonic.

    /**
     * @brief Construct a field for blending mode (l1, m1) into (l2, m2).
     * @param l1 Degree of the first harmonic.
     * @param m1 Order of the first harmonic.
     * @param l2 Degree of the second harmonic.
     * @param m2 Order of the second harmonic.
     * @param blend Morph fraction in [0, 1] from the first toward the second;
     * folded to 0 when both modes are the same.
     * @param q Orientation quaternion of the shape.
     */
    HarmonicField(int l1, int m1, int l2, int m2, float blend,
                  math::Quaternion q)
        : l1(l1), m1(m1), l2(l2), m2(m2),
          blend(l1 == l2 && m1 == m2 ? 0.0f : blend),
          orientation_conj(q.conjugate()),
          normalization_1(SHMath::harmonic_scale(l1, m1)),
          normalization_2(SHMath::harmonic_scale(l2, m2)) {}

    /**
     * @brief Sample the (possibly blended) harmonic at world point p.
     * @param p World-space sample point.
     * @return Signed field value; zero-crossings are the lobe boundaries.
     */
    float sample(const math::Vector &p) const {
      math::Vector local = orientation_conj.apply(p);

      // The local frame varies across a screen row, so no part of the
      // harmonic can be hoisted to a per-row precompute.
      float val = SHMath::spherical_harmonic(l1, m1, local, normalization_1);
      if (blend > 0.001f) {
        float val2 = SHMath::spherical_harmonic(l2, m2, local, normalization_2);
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
   */
  HS_COLD_MEMBER void init() override {
    configure_presets(PRESET_IDS.size());
    register_param("Amplitude", &params.amplitude, 0.1f, 10.0f);

    baked_palette.bake(persistent_arena, Palettes::RICH_SUNSET);

    current_idx = SEED_MODE_IDX;

    math::Vector axis = math::Vector(0.5f, 1.0f, 0.2f).normalized();
    timeline.add(0, Animation::Rotation<W>(
                        orientation, axis, 2 * math::PI_F * SPIN_TURNS,
                        SPIN_FRAMES, math::ease_linear, {.repeat = true}));

    start_morph();
  }

  /**
   * @brief Render one frame of the morphing harmonic.
   */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sh_timeline_step);
      timeline.step(canvas);
    }

    const HarmonicFrameState current = frame_state();
    if (!history.initialized)
      history.capture(current);
    if (!Spacetime::flat_surface()) {
      draw_columns(canvas, history.previous, current);
    } else {
      const typename RenderPipeline::Frame frame =
          RenderPipeline::prepare(current);
      HS_PROFILE(sh_rasterize);
      Scan::Shader::draw<W, H>(canvas, [&frame](const math::Vector &view) {
        return RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
      });
    }
    history.capture(current);
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

  /**
   * @brief Starts at the mode of the target preset with no morph pending.
   * @param change Preset change; only its target index is read.
   * @return Always true.
   */
  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    ++morph_generation;
    if (change.origin != PresetChangeOrigin::SYNCHRONIZED)
      ++morph_epoch;
    current_idx = mode_index_for_preset(change.to);
    next_idx = current_idx;
    morph_alpha = 0.0f;
    return true;
  }

  friend struct ::hs_test::effects_tests::SphericalHarmonicsWhiteBox;

  /**
   * @brief Ambient-occlusion shaping constants for colorize_harmonic.
   * @details AO_FALLOFF is the amplitude-scaled field magnitude
   * (|val| * Amplitude) at which shading saturates, AO_AMBIENT is the ambient
   * floor, and AO_RANGE is the range above it (AMBIENT + RANGE = 1.0 → full
   * brightness at saturation).
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
    math::Quaternion orientation;
    uint32_t morph_epoch; ///< Changes when a preset snap breaks the morph.

    /** @brief Whether both morph endpoints are one mode.
     * @return True for a pure harmonic. */
    bool pure() const { return l1 == l2 && m1 == m2; }

    /** @brief Morph fraction as the field sees it.
     * @return Zero for a pure harmonic, else the blend. */
    float effective_blend() const { return pure() ? 0.0f : blend; }
  };

  struct HarmonicMaterialState {
    const BakedPalette *palette;
    float amplitude;
  };

  struct HarmonicFrameState {
    HarmonicSourceState source;
    HarmonicMaterialState material;

    /**
     * @brief Samples the spin, morph chain and amplitude between two frames.
     * @details The morph is re-expressed on the current leg. A leg rollover
     * continues the outgoing leg at full blend into the incoming leg at zero;
     * harmonic indices are never interpolated. Unrelated legs and a changed
     * morph epoch take the current morph.
     * @param from Previous frame.
     * @param to Current frame.
     * @param fraction Frame fraction in [0, 1).
     * @param rotation Prepared orientation interval.
     * @return The frame at @p fraction.
     */
    HS_FLASH_INLINE static HarmonicFrameState
    sample(const HarmonicFrameState &from, const HarmonicFrameState &to,
           float fraction, const Spacetime::TemporalRotation &rotation) {
      HarmonicFrameState frame = to;
      frame.source.orientation = rotation.sample(fraction);
      frame.material.amplitude =
          from.material.amplitude +
          (to.material.amplitude - from.material.amplitude) * fraction;
      const HarmonicSourceState &a = from.source;
      const HarmonicSourceState &b = to.source;
      if (a.morph_epoch != b.morph_epoch)
        return frame;
      const bool same_leg = a.l1 == b.l1 && a.m1 == b.m1 &&
                            (a.pure() || (a.l2 == b.l2 && a.m2 == b.m2));
      const bool rollover = !a.pure() && a.l2 == b.l1 && a.m2 == b.m1;
      if (!same_leg && !rollover)
        return frame;
      const float start = same_leg ? a.effective_blend() : a.blend - 1.0f;
      const float position = start + (b.effective_blend() - start) * fraction;
      if (position >= 0.0f) {
        frame.source.blend = position;
      } else {
        frame.source.l1 = a.l1;
        frame.source.m1 = a.m1;
        frame.source.l2 = a.l2;
        frame.source.m2 = a.m2;
        frame.source.blend = position + 1.0f;
      }
      return frame;
    }
    /** @brief Samples one frame fraction from endpoint orientations.
     * @param from Previous frame.
     * @param to Current frame.
     * @param fraction Frame fraction in [0, 1).
     * @return Interpolated frame. */
    HS_FLASH_INLINE static HarmonicFrameState
    sample(const HarmonicFrameState &from, const HarmonicFrameState &to,
           float fraction) {
      return sample(from, to, fraction,
                    Spacetime::TemporalRotation(from.source.orientation,
                                                to.source.orientation));
    }
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
        material.palette->get(fminf(1.0f, abs_val * material.amplitude));

    constexpr float NEG_LOBE_GREEN_SCALE = 0.8f;
    const Color4 neg = Color4(
        Pixel(pos.color.b,
              static_cast<uint16_t>(pos.color.g * NEG_LOBE_GREEN_SCALE + 0.5f),
              pos.color.r),
        pos.alpha);

    constexpr float TRANSITION = 0.03f;
    const float blend_t = math::quintic_kernel(
        hs::clamp(val / TRANSITION * 0.5f + 0.5f, 0.0f, 1.0f));
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

  /** @brief The frame the timeline has just advanced to.
   * @return Current harmonic pair, spin, palette and amplitude. */
  HarmonicFrameState frame_state() const {
    auto [l1, m1] = SHMath::decode_lm(current_idx);
    auto [l2, m2] = SHMath::decode_lm(next_idx);
    return {{l1, m1, l2, m2, morph_alpha, orientation.get(), morph_epoch},
            {&baked_palette.view(), params.amplitude}};
  }

  /**
   * @brief Shades each column at its sweep time between two frames.
   * @param canvas Destination canvas.
   * @param from Previous frame.
   * @param to Current frame.
   */
  HS_FLASH_MEMBER static void draw_columns(Canvas &canvas,
                                           const HarmonicFrameState &from,
                                           const HarmonicFrameState &to) {
    HS_PROFILE(sh_rasterize);
    const Spacetime::TemporalRotation rotation(from.source.orientation,
                                               to.source.orientation);
    Scan::Shader::draw_column_times<W, H>(
        canvas,
        [&](int c) {
          return RenderPipeline::prepare(HarmonicFrameState::sample(
              from, to, Spacetime::SamplingSurface<W>::time(c), rotation));
        },
        [](const auto &frame, int) {
          return [&frame](const math::Vector &view) {
            return RenderPipeline::evaluate(view, frame.ctx, frame.prepared);
          };
        });
  }

  /**
   * @brief Advances to the next harmonic and animates the morph toward it.
   * @details On completion commits it as current and recurses, yielding an
   * endless morph chain. Pausable: the chain is this effect's preset
   * choreography, so "Pause Presets" holds the current blend.
   */
  HS_COLD_MEMBER void start_morph() {
    const uint32_t generation = morph_generation;
    next_idx = (current_idx % MAX_MODE_IDX) + 1;
#if defined(HS_PROFILE_ENABLE)
    hs::log("Mode: %d/%d", next_idx, MAX_MODE_IDX);
#endif

    timeline.add_pausable(
        0,
        Animation::Transition(morph_alpha, 1.0f, MORPH_FRAMES,
                              math::ease_linear)
            .then([this, generation]() {
              if (generation != morph_generation) {
                morph_alpha = 0.0f;
                start_morph();
                return;
              }
              current_idx = next_idx;
              morph_alpha = 0.0f;
              const bool synchronized =
                  synchronizePreset(preset_index_for_mode(current_idx));
              HS_CHECK(synchronized,
                       "SphericalHarmonics preset synchronization failed");
              start_morph();
            }),
        &anims_paused);
  }

  math::Orientation<> orientation; /**< Current sphere orientation. */
  Timeline timeline;               /**< Drives spin and morph animations. */
  BakedPaletteStorage
      baked_palette; /**< Precomputed color LUT for the shader. */

  // One palette LUT in the persistent arena.
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
  static constexpr int MAX_MODE_IDX = (MAX_DEGREE + 1) * (MAX_DEGREE + 1) - 1;
  static_assert(PRESET_IDS.size() == MAX_MODE_IDX);
  // Initial mode (l=2, m=0); the constant mode (idx 0) is excluded from the roll.
  static constexpr int SEED_MODE_IDX = 6;
  static_assert(SEED_MODE_IDX > 0 && SEED_MODE_IDX <= MAX_MODE_IDX,
                "seed mode must be a valid, non-constant harmonic index");
  // Looped continuous spin: SPIN_TURNS revolutions over SPIN_FRAMES frames.
  static constexpr float SPIN_TURNS = 100.0f;
  static constexpr int SPIN_FRAMES = 10000;

  int current_idx = 0; /**< Flat index of the currently displayed mode. */
  int next_idx = 0;    /**< Flat index of the mode being morphed toward. */
  float morph_alpha =
      0.0f; /**< Morph progress in [0, 1] from current to next. */
  uint32_t morph_generation = 0;
  uint32_t morph_epoch = 0; /**< Bumped by preset snaps outside the chain. */
  /** Previous frame for column-time sampling. */
  Spacetime::FrameHistory<HarmonicFrameState> history;

  /**
   * @brief User-tunable parameters for the visualizer.
   */
  struct Params {
    float amplitude = 3.2f; /**< Field-value gain applied before coloring. */
  } params;
};
