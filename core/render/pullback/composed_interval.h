/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_interval.h
 * @brief Retained frame endpoints and column-time sampling for composed effects.
 */

namespace Pullback {

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

/** @brief Shared preparation for column samples of one frame interval. */
template <typename Frame> class FrameInterval {
  using Params = decltype(Frame::params);
  const Frame &from;
  const Frame &to;
  Spacetime::TemporalRotation projection;
  Spacetime::TemporalRotation outer;
  decltype(Frame::resources) deltas{};
  Params interior_params;
  bool constant_params;
  uint32_t constant_families = 0;
  float oscillation_delta;
  bool same_palette;
  bool phase_compatible = true;

public:
  /** @brief Prepares parameter identity and periodic clock deltas.
   * @param a Previous frame endpoint.
   * @param b Current frame endpoint. */
  HS_FLASH_INLINE FrameInterval(const Frame &a, const Frame &b)
      : from(a), to(b),
        projection(a.projection_conjugate, b.projection_conjugate),
        outer(a.outer_conjugate, b.outer_conjugate), interior_params(b.params),
        constant_params(std::memcmp(&a.params, &b.params, sizeof(Params)) == 0),
        oscillation_delta(Spacetime::ClockInterval::prepare_delta(
            a.palette_oscillation_phase, b.palette_oscillation_phase, 1.0f)),
        same_palette(a.palette && b.palette &&
                     a.palette->same_entries(*b.palette) &&
                     (a.hue_rotation_lut == b.hue_rotation_lut ||
                      (a.hue_rotation_lut && b.hue_rotation_lut &&
                       std::memcmp(a.hue_rotation_lut, b.hue_rotation_lut,
                                   Color::HueRotationLutView::SIZE *
                                       sizeof(Pixel)) == 0))) {
    static_assert(std::is_trivially_copyable_v<Params>);
    static_assert(
        []<typename... Resources>(
            ComposedDetail::ResourceList<Resources...> *) {
          return sizeof...(Resources);
        }(static_cast<typename Params::ResourceTypes *>(nullptr)) <= 32,
        "constant_families holds one bit per parameter family");
    if (constant_params) {
      interior_params = interpolate(a.params, b.params, 0.5f);
    } else {
      // A family equal at both endpoints takes its interior value throughout.
      int family = 0;
      interior_params.visit(
          [&]<typename Resource>(auto &value) __attribute__((always_inline)) {
            const auto &start = a.params.template get<Resource::KEY>();
            const auto &end = b.params.template get<Resource::KEY>();
            if (std::memcmp(&start, &end, sizeof(end)) == 0) {
              value = interpolate(start, end, 0.5f);
              constant_families |= 1u << family;
            }
            ++family;
          });
    }
    b.params.visit([&]<typename Resource>(const auto &) {
      auto &delta = deltas.template get<Resource::KEY>();
      const auto &start = a.resources.template get<Resource::KEY>();
      const auto &end = b.resources.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        delta.primary = Spacetime::ClockInterval::prepare_delta(
            start.primary, end.primary, math::TWO_PI_F);
        delta.secondary = Spacetime::ClockInterval::prepare_delta(
            start.secondary, end.secondary, math::TWO_PI_F);
        delta.angle = Spacetime::ClockInterval::prepare_delta(
            start.angle, end.angle, math::TWO_PI_F);
        delta.noise_time = Spacetime::ClockInterval::prepare_delta(
            start.noise_time, end.noise_time, 1.0f);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        delta.phase = Spacetime::ClockInterval::prepare_delta(start.phase,
                                                              end.phase, 1.0f);
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          delta.rotation = Spacetime::ClockInterval::prepare_delta(
              start.rotation, end.rotation, math::TWO_PI_F);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE &&
                           !std::is_same_v<typename Resource::Family,
                                           PeriodicRippleParams> &&
                           !std::is_same_v<
                               typename Resource::Family,
                               SphereAnimation::BurstRippleParams>) {
        delta.phase = Spacetime::ClockInterval::prepare_delta(start.phase,
                                                              end.phase, 1.0f);
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusLensParams>) {
        delta.angle = Spacetime::ClockInterval::prepare_delta(
            start.angle, end.angle, math::TWO_PI_F);
        delta.magnitude = end.magnitude - start.magnitude;
        delta.generated = lens_phase_continuous(
            a.params.template get<Resource::KEY>(),
            b.params.template get<Resource::KEY>(), start, end);
      } else if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        const auto &previous =
            a.params.template get<Resource::KEY>().modulation;
        const auto &current = b.params.template get<Resource::KEY>().modulation;
        phase_compatible =
            phase_compatible && Modulation::phase_compatible(previous, current);
      }
    });
  }

  /** @brief Samples the prepared interval, preserving both endpoints.
   * @param fraction Frame progress in [0, 1].
   * @return Interpolated context borrowing the endpoint resources. */
  HS_FLASH_INLINE Frame at(float fraction) const {
    Frame frame;
    at_into(frame, fraction);
    return frame;
  }

  /** @brief Writes a sampled context into storage separate from the endpoints.
   * @param frame Receives the context borrowing the endpoint resources.
   * @param fraction Frame progress in [0, 1]. */
  HS_FLASH_INLINE void at_into(Frame &frame, float fraction) const {
    if (fraction <= 0.0f) {
      frame = from;
      hold_discrete_lenses(frame, from, to);
      to.params.visit([&]<typename Resource>(const auto &family) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>) {
          auto &resource = frame.resources.template get<Resource::KEY>();
          const auto &end = to.resources.template get<Resource::KEY>();
          const uint32_t DELTA = end.frames - resource.frames;
          resource.radial_tables = DELTA <= 1 ? end.radial_tables : nullptr;
          resource.radial_fraction = -1.0f;
          if (std::memcmp(&family, &from.params.template get<Resource::KEY>(),
                          sizeof(family)) != 0)
            resource.radial_tables = nullptr;
        }
        if constexpr (Resource::KIND == ResourceKind::MODULATION) {
          if (!Modulation::same_topology(
                  from.params.template get<Resource::KEY>().modulation,
                  family.modulation)) {
            frame.params.template get<Resource::KEY>() = family;
            frame.resources.template get<Resource::KEY>() =
                to.resources.template get<Resource::KEY>();
          }
        }
      });
      frame.temporal = {from.palette, from.hue_rotation_lut, 0.0f,
                        same_palette};
      return;
    }
    frame = to;
    if (fraction >= 1.0f) {
      frame.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>)
          frame.resources.template get<Resource::KEY>().radial_fraction = -1.0f;
      });
      return;
    }
    frame.projection_conjugate = projection.sample(fraction);
    frame.outer_conjugate = outer.sample(fraction);
    if (constant_params) {
      frame.params = interior_params;
    } else {
      int family = 0;
      frame.params.visit(
          [&]<typename Resource>(auto &value) __attribute__((always_inline)) {
            value = constant_families & (1u << family)
                        ? interior_params.template get<Resource::KEY>()
                        : interpolate(from.params.template get<Resource::KEY>(),
                                      to.params.template get<Resource::KEY>(),
                                      fraction);
            ++family;
          });
    }
    frame.palette_mapping = Color::PaletteMappingWeights::lerp(
        from.palette_mapping, to.palette_mapping, fraction);
    frame.palette_oscillation_phase = Spacetime::ClockInterval::apply(
        from.palette_oscillation_phase, oscillation_delta, fraction);
    frame.params.visit([&]<typename Resource>(const auto &family) {
      auto &clock = frame.resources.template get<Resource::KEY>();
      const auto &start = from.resources.template get<Resource::KEY>();
      const auto &delta = deltas.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        clock.primary = Spacetime::ClockInterval::apply(
            start.primary, delta.primary, fraction);
        clock.secondary = Spacetime::ClockInterval::apply(
            start.secondary, delta.secondary, fraction);
        clock.angle =
            Spacetime::ClockInterval::apply(start.angle, delta.angle, fraction);
        clock.noise_time = Spacetime::ClockInterval::apply(
            start.noise_time, delta.noise_time, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        clock.phase =
            Spacetime::ClockInterval::apply(start.phase, delta.phase, fraction);
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          clock.rotation = Spacetime::ClockInterval::apply(
              start.rotation, delta.rotation, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>) {
          const auto &end = to.resources.template get<Resource::KEY>();
          const uint32_t DELTA = end.frames - start.frames;
          for (auto &wave : clock.waves)
            wave.age -= DELTA * (1.0f - fraction);
          clock.radial_fraction = DELTA == 1 ? fraction : 1.0f;
          if (DELTA > 1 ||
              std::memcmp(&family, &to.params.template get<Resource::KEY>(),
                          sizeof(family)) != 0 ||
              std::memcmp(&family, &from.params.template get<Resource::KEY>(),
                          sizeof(family)) != 0)
            clock.radial_tables = nullptr;
        } else if constexpr (std::is_same_v<typename Resource::Family,
                                            PeriodicRippleParams>)
          clock.phase = Spacetime::sample_clock(
              start.phase, to.resources.template get<Resource::KEY>().phase,
              fraction, family.period);
        else
          clock.phase = Spacetime::ClockInterval::apply(start.phase,
                                                        delta.phase, fraction);
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusLensParams>) {
        const auto &end = to.params.template get<Resource::KEY>();
        frame.params.template get<Resource::KEY>() =
            delta.generated
                ? sample_lens_phase(
                      end,
                      Spacetime::ClockInterval::apply(start.angle, delta.angle,
                                                      fraction),
                      Spacetime::ClockInterval::apply(
                          start.magnitude, delta.magnitude, fraction))
                : end;
      }
    });
    bool same_sources = true;
    frame.params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        auto &resource = frame.resources.template get<Resource::KEY>();
        resource.previous_modulation =
            Modulation::same_topology(
                from.params.template get<Resource::KEY>().modulation,
                family.modulation)
                ? from.resources.template get<Resource::KEY>().modulation
                : resource.modulation;
        if (family.modulation.enabled &&
            family.modulation.mechanic.amount != 0.0f &&
            resource.modulation != resource.previous_modulation)
          same_sources = same_sources && resource.modulation &&
                         resource.previous_modulation &&
                         Modulation::same_field(*resource.modulation,
                                                *resource.previous_modulation);
      }
    });
    frame.temporal = {from.palette, from.hue_rotation_lut,
                      same_palette && same_sources ? 1.0f : fraction,
                      same_palette, phase_compatible};
  }
};

} // namespace Pullback

namespace Spacetime {

/** @brief Samples spatial state and clocks from a prepared frame interval.
 * @tparam Frame Composed temporal frame context.
 * @param from Previous frame with retained resources.
 * @param to Current frame.
 * @param fraction Frame fraction in [0, 1].
 * @return Context borrowing the endpoint resources. */
template <typename Frame>
HS_FLASH_INLINE Frame sample_frame(const Frame &from, const Frame &to,
                                   float fraction) {
  return Pullback::FrameInterval(from, to).at(fraction);
}

/** @brief Owned palette and hue tables for the previous frame endpoint. */
template <typename Frame, bool Hue>
  requires BorrowsColorResources<Frame>
struct FrameHistory<Frame, true, Hue> {
  Frame previous{};            ///< Retained endpoint.
  BakedPaletteStorage palette; ///< Owned palette entries.
  std::array<Pixel, Hue ? Pullback::Color::HueRotationLutView::SIZE : 0>
      hue_rotation; ///< Retained hue table.
  Pullback::ComposedDetail::ResourceStorage<
      Pullback::ComposedDetail::ResourceModulationHistory,
      typename decltype(Frame::params)::ResourceTypes>
      resources;            ///< Retained layer sources.
  bool initialized = false; ///< Whether a frame was captured.
  /** @brief Retains palette and per-layer sources.
   * @param frame Frame to retain.
   * @param arena Palette storage arena. */
  HS_FLASH_MEMBER void capture(const Frame &frame, Arena &arena) {
    if (!initialized)
      palette.clone_from(*frame.palette, arena);
    else
      palette.rebake_copy(*frame.palette);
    previous = frame;
    previous.palette = &palette.view();
    if constexpr (Hue) {
      std::copy_n(frame.hue_rotation_lut, hue_rotation.size(),
                  hue_rotation.data());
      previous.hue_rotation_lut = hue_rotation.data();
    }
    frame.params.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
        if constexpr (Resource::Family::SOURCE !=
                      Pattern::SourceKind::DISPLACEMENT) {
          auto &owned = resources.template get<Resource::KEY>();
          auto &view = previous.resources.template get<Resource::KEY>();
          view.modulation = Pullback::Modulation::retain(
              frame.resources.template get<Resource::KEY>().modulation,
              owned.modulation, owned.table);
          view.previous_modulation = view.modulation;
        }
      }
    });
    previous.temporal = {previous.palette, previous.hue_rotation_lut, 1.0f};
    initialized = true;
  }
};

} // namespace Spacetime
