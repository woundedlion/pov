/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_temporal.h
 * @brief Retained frame endpoints and column-time sampling for composed effects.
 */

/** @brief Gives a frame the current lens wherever a Mobius lens cannot
 *  follow its generator phase between the endpoints.
 * @tparam Frame Composed temporal frame context.
 * @param frame Frame to correct.
 * @param from Previous frame.
 * @param to Current frame. */
template <typename Frame>
HS_FLASH_INLINE void hold_discrete_lenses(Frame &frame, const Frame &from,
                                          const Frame &to) {
  to.params.visit([&]<typename Resource>(const auto &lens) {
    if constexpr (std::is_same_v<typename Resource::Family, MobiusLensParams>)
      if (!lens_phase_continuous(from.params.template get<Resource::KEY>(),
                                 lens,
                                 from.resources.template get<Resource::KEY>(),
                                 to.resources.template get<Resource::KEY>()))
        frame.params.template get<Resource::KEY>() = lens;
  });
}

/** @brief Samples spatial state and clocks between consecutive frame endpoints.
 * @tparam Frame Composed temporal frame context.
 * @param from Previous frame, with retained palette and hue tables.
 * @param to Current frame.
 * @param fraction Fraction of one frame, in [0, 1].
 * @param constant_params Interior parameter sample when both endpoints agree.
 * @return Column context with resources borrowed from the endpoints. */
template <typename Frame>
HS_FLASH_INLINE Frame
sample_frame(const Frame &from, const Frame &to, float fraction,
             const decltype(Frame::params) *constant_params = nullptr) {
  if (fraction <= 0.0f) {
    Frame frame = from;
    hold_discrete_lenses(frame, from, to);
    return frame;
  }
  if (fraction >= 1.0f)
    return to;
  Frame frame = to;
  frame.projection_conjugate =
      math::slerp(from.projection_conjugate, to.projection_conjugate, fraction);
  frame.outer_conjugate =
      math::slerp(from.outer_conjugate, to.outer_conjugate, fraction);
  frame.params = constant_params
                     ? *constant_params
                     : interpolate(from.params, to.params, fraction);
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
    } else if constexpr (std::is_same_v<typename Resource::Family,
                                        MobiusLensParams>) {
      const auto &start = from.params.template get<Resource::KEY>();
      const auto &end = to.params.template get<Resource::KEY>();
      frame.params.template get<Resource::KEY>() =
          lens_phase_continuous(start, end, a, b)
              ? sample_lens_phase(
                    end,
                    sample_clock(a.angle, b.angle, fraction, math::TWO_PI_F),
                    a.magnitude + (b.magnitude - a.magnitude) * fraction)
              : end;
    }
  });
  frame.temporal = {from.palette,
                    from.hue_rotation_lut,
                    from.hue_noise_lut,
                    fraction,
                    from.hue_rotation_linear,
                    from.palette_alpha,
                    from.hue_noise_linear};
  return frame;
}

/** @brief Shared preparation for column samples of one frame interval. */
template <typename Frame> class TemporalFrameSampler {
  using Params = decltype(Frame::params);
  const Frame &from;
  const Frame &to;
  Temporal::TemporalRotation projection;
  Temporal::TemporalRotation outer;
  decltype(Frame::resources) deltas{};
  Params interior_params;
  bool constant_params;
  float oscillation_delta;
  bool same_noise;
  bool same_palette;

public:
  /** @brief Prepares parameter identity and periodic clock deltas.
   * @param a Previous frame endpoint.
   * @param b Current frame endpoint. */
  HS_FLASH_INLINE TemporalFrameSampler(const Frame &a, const Frame &b)
      : from(a), to(b),
        projection(a.projection_conjugate, b.projection_conjugate),
        outer(a.outer_conjugate, b.outer_conjugate), interior_params(b.params),
        constant_params(std::memcmp(&a.params, &b.params, sizeof(Params)) == 0),
        oscillation_delta(std::remainder(
            b.palette_oscillation_phase - a.palette_oscillation_phase, 1.0f)),
        same_noise(a.hue_noise_lut && b.hue_noise_lut &&
                   std::memcmp(a.hue_noise_lut, b.hue_noise_lut,
                               Color::HueNoiseLutView::SIZE) == 0),
        same_palette(a.palette && b.palette &&
                     a.palette->same_entries(*b.palette) &&
                     (a.hue_rotation_lut == b.hue_rotation_lut ||
                      (a.hue_rotation_lut && b.hue_rotation_lut &&
                       std::memcmp(a.hue_rotation_lut, b.hue_rotation_lut,
                                   Color::HueRotationLutView::SIZE *
                                       sizeof(Pixel)) == 0))) {
    static_assert(std::is_trivially_copyable_v<Params>);
    if (constant_params)
      interior_params = interpolate(a.params, b.params, 0.5f);
    b.params.visit([&]<typename Resource>(const auto &) {
      auto &delta = deltas.template get<Resource::KEY>();
      const auto &start = a.resources.template get<Resource::KEY>();
      const auto &end = b.resources.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        delta.primary =
            std::remainder(end.primary - start.primary, math::TWO_PI_F);
        delta.secondary =
            std::remainder(end.secondary - start.secondary, math::TWO_PI_F);
        delta.angle = std::remainder(end.angle - start.angle, math::TWO_PI_F);
        delta.noise_time =
            std::remainder(end.noise_time - start.noise_time, 1.0f);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        delta.phase = std::remainder(end.phase - start.phase, 1.0f);
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          delta.rotation =
              std::remainder(end.rotation - start.rotation, math::TWO_PI_F);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE &&
                           !std::is_same_v<typename Resource::Family,
                                           PeriodicRippleParams>) {
        delta.phase = std::remainder(end.phase - start.phase, 1.0f);
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusLensParams>) {
        delta.angle = std::remainder(end.angle - start.angle, math::TWO_PI_F);
        delta.magnitude = end.magnitude - start.magnitude;
        delta.generated = lens_phase_continuous(
            a.params.template get<Resource::KEY>(),
            b.params.template get<Resource::KEY>(), start, end);
      }
    });
  }

  /** @brief Samples the prepared interval, preserving both endpoints.
   * @param fraction Frame progress in [0, 1].
   * @return Interpolated context borrowing the endpoint resources. */
  HS_FLASH_INLINE Frame sample(float fraction) const {
    if (fraction <= 0.0f) {
      Frame frame = from;
      hold_discrete_lenses(frame, from, to);
      return frame;
    }
    if (fraction >= 1.0f)
      return to;
    Frame frame = to;
    frame.projection_conjugate = projection.sample(fraction);
    frame.outer_conjugate = outer.sample(fraction);
    frame.params = constant_params
                       ? interior_params
                       : interpolate(from.params, to.params, fraction);
    frame.palette_mapping = Color::PaletteMappingWeights::lerp(
        from.palette_mapping, to.palette_mapping, fraction);
    frame.palette_oscillation_phase =
        from.palette_oscillation_phase + oscillation_delta * fraction;
    frame.params.visit([&]<typename Resource>(const auto &family) {
      auto &clock = frame.resources.template get<Resource::KEY>();
      const auto &start = from.resources.template get<Resource::KEY>();
      const auto &delta = deltas.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        clock.primary = start.primary + delta.primary * fraction;
        clock.secondary = start.secondary + delta.secondary * fraction;
        clock.angle = start.angle + delta.angle * fraction;
        clock.noise_time = start.noise_time + delta.noise_time * fraction;
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        clock.phase = start.phase + delta.phase * fraction;
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          clock.rotation = start.rotation + delta.rotation * fraction;
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     PeriodicRippleParams>)
          clock.phase = sample_clock(
              start.phase, to.resources.template get<Resource::KEY>().phase,
              fraction, family.period);
        else
          clock.phase = start.phase + delta.phase * fraction;
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusLensParams>) {
        const auto &end = to.params.template get<Resource::KEY>();
        frame.params.template get<Resource::KEY>() =
            delta.generated ? sample_lens_phase(
                                  end, start.angle + delta.angle * fraction,
                                  start.magnitude + delta.magnitude * fraction)
                            : end;
      }
    });
    frame.temporal = {from.palette,
                      from.hue_rotation_lut,
                      from.hue_noise_lut,
                      fraction,
                      from.hue_rotation_linear,
                      from.palette_alpha,
                      from.hue_noise_linear,
                      same_noise,
                      same_palette};
    return frame;
  }
};

/** @brief Owned palette and hue tables for the previous frame endpoint. */
template <typename Frame, HueMode Hue, bool FloatHue>
  requires BorrowsColorResources<Frame>
struct TemporalHistory<Frame, Hue, true, FloatHue> {
  Frame previous{}; ///< Previous endpoint with owned color resource views.
  BakedPaletteStorage palette; ///< Retained palette table.
  std::array<Pixel, Hue != HueMode::NONE ? Color::HueRotationLutView::SIZE : 0>
      hue_rotation; ///< Retained hue-rotation table.
  std::array<int8_t, Hue == HueMode::NOISE ? Color::HueNoiseLutView::SIZE : 0>
      hue_noise; ///< Retained hue-noise table.
  std::array<math::Vector, FloatHue ? Color::HueRotationLutView::SIZE : 0>
      hue_rotation_linear; ///< Retained float hue table.
  std::array<float, FloatHue ? Color::HueNoiseLutView::SIZE : 0>
      hue_noise_linear;     ///< Retained float noise table.
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
    previous.hue_rotation_linear = nullptr;
    previous.hue_noise_linear = nullptr;
    if constexpr (Hue != HueMode::NONE) {
      std::copy_n(frame.hue_rotation_lut, hue_rotation.size(),
                  hue_rotation.data());
      previous.hue_rotation_lut = hue_rotation.data();
    }
    if constexpr (Hue == HueMode::NOISE) {
      std::copy_n(frame.hue_noise_lut, hue_noise.size(), hue_noise.data());
      previous.hue_noise_lut = hue_noise.data();
    }
    if constexpr (FloatHue) {
      std::copy_n(frame.hue_rotation_linear, hue_rotation_linear.size(),
                  hue_rotation_linear.data());
      previous.hue_rotation_linear = hue_rotation_linear.data();
      std::copy_n(frame.hue_noise_linear, hue_noise_linear.size(),
                  hue_noise_linear.data());
      previous.hue_noise_linear = hue_noise_linear.data();
    }
    previous.temporal = {previous.palette,
                         previous.hue_rotation_lut,
                         previous.hue_noise_lut,
                         1.0f,
                         previous.hue_rotation_linear,
                         previous.palette_alpha,
                         previous.hue_noise_linear};
    initialized = true;
  }
};
