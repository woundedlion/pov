/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_temporal.h
 * @brief Retained frame endpoints and column-time sampling for composed effects.
 */

/** @brief Column time within a two-arm half-turn.
 * @param column Absolute canvas column.
 * @param width Even canvas width.
 * @return Fraction of one effect frame, in [0, 1). */
constexpr float sweep_fraction(int column, int width) {
  return static_cast<float>(column % (width / 2)) / (width / 2);
}

/** @brief Interpolates the shortest arc of a periodic clock.
 * @param from Previous clock value.
 * @param to Current clock value.
 * @param fraction Fraction of one frame, in [0, 1].
 * @param period Clock period; the frame advance must be less than half of it.
 * @return Interpolated clock, equivalent modulo the period. */
inline float sample_clock(float from, float to, float fraction, float period) {
  return from + std::remainder(to - from, period) * fraction;
}

/** @brief Samples spatial state and clocks between consecutive frame endpoints.
 * @tparam Frame Composed temporal frame context.
 * @param from Previous frame, with retained palette and hue tables.
 * @param to Current frame.
 * @param fraction Fraction of one frame, in [0, 1].
 * @return Column context with resources borrowed from the endpoints. */
template <typename Frame>
HS_FLASH_INLINE Frame sample_frame(const Frame &from, const Frame &to,
                                   float fraction) {
  if (fraction <= 0.0f)
    return from;
  if (fraction >= 1.0f)
    return to;
  Frame frame = to;
  frame.projection_conjugate =
      math::slerp(from.projection_conjugate, to.projection_conjugate, fraction);
  frame.outer_conjugate =
      math::slerp(from.outer_conjugate, to.outer_conjugate, fraction);
  frame.params = interpolate(from.params, to.params, fraction);
  frame.palette_mapping = Color::PaletteMappingWeights::lerp(
      from.palette_mapping, to.palette_mapping, fraction);
  frame.palette_oscillation_phase =
      sample_clock(from.palette_oscillation_phase, to.palette_oscillation_phase,
                   fraction, 1.0f);
  frame.params.visit([&]<typename Resource>(const auto &family) {
    auto &clock = frame.resources.template get<Resource::KEY>();
    const auto &a = from.resources.template get<Resource::KEY>();
    const auto &b = to.resources.template get<Resource::KEY>();
    if constexpr (Resource::KIND == ResourceKind::SOURCE) {
      clock.primary =
          sample_clock(a.primary, b.primary, fraction, math::TWO_PI_F);
      clock.secondary =
          sample_clock(a.secondary, b.secondary, fraction, math::TWO_PI_F);
      clock.angle = sample_clock(a.angle, b.angle, fraction, math::TWO_PI_F);
      clock.noise_time =
          sample_clock(a.noise_time, b.noise_time, fraction, 1.0f);
    } else if constexpr (Resource::KIND == ResourceKind::WARP) {
      clock.phase = sample_clock(a.phase, b.phase, fraction, 1.0f);
      if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
        clock.rotation =
            sample_clock(a.rotation, b.rotation, fraction, math::TWO_PI_F);
    } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
      if constexpr (std::is_same_v<typename Resource::Family,
                                   PeriodicRippleParams>)
        clock.phase = sample_clock(a.phase, b.phase, fraction, family.period);
      else
        clock.phase = sample_clock(a.phase, b.phase, fraction, 1.0f);
    }
  });
  frame.temporal = {from.palette, from.hue_rotation_lut, from.hue_noise_lut,
                    fraction};
  return frame;
}

/** @brief Optional retained endpoint storage. */
template <typename Frame, HueMode Hue, bool Enabled> struct TemporalHistory {};

/** @brief Owned palette and hue tables for the previous frame endpoint. */
template <typename Frame, HueMode Hue>
struct TemporalHistory<Frame, Hue, true> {
  Frame previous{}; ///< Previous endpoint with owned color resource views.
  BakedPaletteStorage palette; ///< Retained palette table.
  std::array<Pixel, Hue != HueMode::NONE ? Color::HueRotationLutView::SIZE : 0>
      hue_rotation; ///< Retained hue-rotation table.
  std::array<int8_t, Hue == HueMode::NOISE ? Color::HueNoiseLutView::SIZE : 0>
      hue_noise;            ///< Retained hue-noise table.
  bool initialized = false; ///< Whether a previous endpoint exists.

  /** @brief Retains the current endpoint before the next animation step.
   * @param frame Current frame.
   * @param arena Palette storage arena, used only on the first capture. */
  HS_FLASH_MEMBER void capture(const Frame &frame, Arena &arena) {
    if (!initialized)
      palette.clone_from(*frame.palette, arena);
    else
      palette.rebake_copy(*frame.palette);
    previous = frame;
    previous.palette = &palette.view();
    if constexpr (Hue != HueMode::NONE) {
      std::copy_n(frame.hue_rotation_lut, hue_rotation.size(),
                  hue_rotation.data());
      previous.hue_rotation_lut = hue_rotation.data();
    }
    if constexpr (Hue == HueMode::NOISE) {
      std::copy_n(frame.hue_noise_lut, hue_noise.size(), hue_noise.data());
      previous.hue_noise_lut = hue_noise.data();
    }
    previous.temporal = {previous.palette, previous.hue_rotation_lut,
                         previous.hue_noise_lut, 1.0f};
    initialized = true;
  }
};
