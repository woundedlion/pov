/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <optional>
#include <type_traits>
#include <variant>
#include "core/memory.h"
#include "core/platform/cycle_counting.h"
#include "core/render/pullback/composed_policies.h"

/** @file composed_interval.h
 * @brief Retained frame endpoints and column-time sampling for composed effects.
 */

namespace Pullback {

/** @brief Samples independent spin and wander clocks of a keyed rotation.
 * @param clock Receives sampled clocks and conjugate.
 * @param start Previous endpoint.
 * @param end Current endpoint.
 * @param fraction Frame progress. */
HS_HOT_FLASH_MEMBER inline void sample_keyed_rotation(
    ComposedDetail::ClockState<RotationParams, ResourceKind::ROTATION> &clock,
    const ComposedDetail::ClockState<RotationParams, ResourceKind::ROTATION>
        &start,
    const ComposedDetail::ClockState<RotationParams, ResourceKind::ROTATION>
        &end,
    float fraction) {
  clock.wander =
      Spacetime::TemporalRotation(start.wander, end.wander).sample(fraction);
  clock.spin_phase = Spacetime::sample_clock(start.spin_phase, end.spin_phase,
                                             fraction, math::TWO_PI_F);
  clock.conjugate =
      (math::make_rotation(math::Y_AXIS, clock.spin_phase) * clock.wander)
          .conjugate();
}

/** @brief Tests palette entries, hue rotation pixels and gamut masks across endpoints.
 * @param a Outgoing palette, or null.
 * @param b Incoming palette, or null.
 * @param a_hue Outgoing hue pixels, or null.
 * @param b_hue Incoming hue pixels, or null.
 * @param a_gamut Outgoing hue gamut mask, or null.
 * @param b_gamut Incoming hue gamut mask, or null.
 * @return Whether the two endpoints have equivalent palette data. */
HS_SIZE_FLASH_MEMBER inline bool
same_interval_palette(const BakedPalette *a, const BakedPalette *b,
                      const Pixel *a_hue, const Pixel *b_hue,
                      const HueRotationGamutMask *a_gamut,
                      const HueRotationGamutMask *b_gamut) {
  return a && b && a->same_entries(*b) &&
         (a_hue == b_hue || (a_hue && b_hue &&
                             std::memcmp(a_hue, b_hue,
                                         Color::HueRotationLutView::SIZE *
                                             sizeof(Pixel)) == 0)) &&
         (a_gamut == b_gamut ||
          (a_gamut && b_gamut && a_gamut->bits == b_gamut->bits));
}

/** @brief Selects a retained scalar source and checks continuous phase compatibility.
 * @param a Outgoing source configuration.
 * @param b Incoming source configuration.
 * @param previous Outgoing prepared source, or null.
 * @param current Incoming prepared source, or null.
 * @param phase_compatible Accumulated phase compatibility across layers.
 * @return Source to use for the outgoing side of the interval. */
HS_HOT_FLASH_MEMBER inline const Modulation::Prepared *
interval_modulation(const Modulation::Config &a, const Modulation::Config &b,
                    const Modulation::Prepared *previous,
                    const Modulation::Prepared *current,
                    bool &phase_compatible) {
  if (!Modulation::same_topology(a, b)) {
    phase_compatible = false;
    return current;
  }
  phase_compatible = phase_compatible && Modulation::phase_compatible(a, b);
  if (previous && current && previous != current &&
      Modulation::same_endpoint_sample(*previous, *current))
    return current;
  return previous;
}

/** @brief Prepares periodic source clocks independently of the effect's resource list.
 * @param delta Receives the signed clock steps.
 * @param a Outgoing source clocks.
 * @param b Incoming source clocks. */
HS_SIZE_FLASH_MEMBER inline void
prepare_source_clocks(ComposedDetail::SourceClocks &delta,
                      const ComposedDetail::SourceClocks &a,
                      const ComposedDetail::SourceClocks &b) {
  using Spacetime::ClockInterval;
  delta.primary =
      ClockInterval::prepare_delta(a.primary, b.primary, math::TWO_PI_F);
  delta.secondary =
      ClockInterval::prepare_delta(a.secondary, b.secondary, math::TWO_PI_F);
  delta.angle = ClockInterval::prepare_delta(a.angle, b.angle, math::TWO_PI_F);
  delta.noise_time =
      ClockInterval::prepare_delta(a.noise_time, b.noise_time, 1.0f);
}

/** @brief Samples all four source clocks along their prepared periodic arcs.
 * @param out Receives sampled clocks.
 * @param start Outgoing source clocks.
 * @param delta Prepared clock steps.
 * @param fraction Position within the frame interval. */
HS_HOT_FLASH_MEMBER inline void
sample_source_clocks(ComposedDetail::SourceClocks &out,
                     const ComposedDetail::SourceClocks &start,
                     const ComposedDetail::SourceClocks &delta,
                     float fraction) {
  using Spacetime::ClockInterval;
  out.primary = ClockInterval::apply(start.primary, delta.primary, fraction);
  out.secondary =
      ClockInterval::apply(start.secondary, delta.secondary, fraction);
  out.angle = ClockInterval::apply(start.angle, delta.angle, fraction);
  out.noise_time =
      ClockInterval::apply(start.noise_time, delta.noise_time, fraction);
}

/** @brief Per-resource data retained between column samples. */
template <typename Resource> struct IntervalDelta {};

/// Clock storage retained for an interpolated resource.
template <typename Resource>
  requires(Resource::KIND == ResourceKind::SOURCE ||
           Resource::KIND == ResourceKind::WARP ||
           Resource::KIND == ResourceKind::TWIST ||
           std::is_same_v<typename Resource::Family, MobiusGeneratorParams> ||
           (Resource::KIND == ResourceKind::SURFACE &&
            !std::is_same_v<typename Resource::Family, PeriodicRippleParams> &&
            !std::is_same_v<typename Resource::Family,
                            SphereAnimation::BurstRippleParams>))
struct IntervalDelta<Resource> : ComposedDetail::ResourceClock<Resource> {};

/// Prepared scalar source retained for a modulation interval.
template <typename Resource>
  requires(Resource::KIND == ResourceKind::MODULATION)
struct IntervalDelta<Resource> {
  /// Outgoing source, or null when absent.
  const Modulation::Prepared *previous_modulation;
};

/** @brief Shared preparation for column samples of one frame interval. */
template <typename Frame> class FrameInterval {
  using Params = decltype(Frame::params);
  using SourceSelector = typename Params::template Family<"source">;
  static constexpr bool HAS_PLANAR =
      requires { typename SourceSelector::PlanarWeights; };

  const Frame &from;
  const Frame &to;
  Spacetime::TemporalRotation projection;
  Spacetime::TemporalRotation outer;
  ComposedDetail::ResourceStorage<IntervalDelta, typename Params::ResourceTypes>
      deltas;
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceMobiusInterval,
                                  typename Params::ResourceTypes>
      mobius_intervals;
  bool constant_params;
  uint32_t constant_families = 0;
  float oscillation_delta;
  bool same_palette;
  bool phase_compatible = true;

  /** @brief Matches projection controls when that family exists.
   * @param a Previous endpoint.
   * @param b Current endpoint.
   * @return Whether both endpoints share projection controls. */
  static bool same_projection_controls(const Frame &a, const Frame &b) {
    if constexpr (Params::template HAS<"projection">)
      return std::memcmp(&a.params.template get<"projection">(),
                         &b.params.template get<"projection">(),
                         sizeof(b.params.template get<"projection">())) == 0;
    return true;
  }

  /** @brief Borrows shared Mobius exponentials for this sampled context.
   * @param frame Sampled context. */
  HS_FLASH_INLINE void bind_mobius_intervals(Frame &frame) const {
    frame.params.visit([&]<typename Resource>(const auto &) {
      if constexpr (std::is_same_v<typename Resource::Family,
                                   MobiusGeneratorParams>) {
        const auto &cache = mobius_intervals.template get<Resource::KEY>();
        frame.resources.template get<Resource::KEY>().mobius_base =
            cache.active ? &cache.base : nullptr;
      }
    });
  }

public:
  /** @brief Prepares parameter and color-resource identities and clock deltas.
   * @param a Previous frame endpoint.
   * @param b Current frame endpoint. */
  HS_SIZE_FLASH_MEMBER FrameInterval(const Frame &a, const Frame &b)
      : from(a), to(b),
        projection(a.projection_conjugate, b.projection_conjugate),
        outer(a.outer_conjugate, b.outer_conjugate),
        constant_params(std::memcmp(&a.params, &b.params, sizeof(Params)) == 0),
        oscillation_delta(Spacetime::ClockInterval::prepare_delta(
            a.palette_oscillation_phase, b.palette_oscillation_phase, 1.0f)),
        same_palette(same_interval_palette(
            a.palette, b.palette, a.hue_rotation_lut, b.hue_rotation_lut,
            a.hue_rotation_gamut, b.hue_rotation_gamut)) {
    static_assert(std::is_trivially_copyable_v<Params>);
    static_assert(
        []<typename... Resources>(
            ComposedDetail::ResourceList<Resources...> *) {
          return sizeof...(Resources);
        }(static_cast<typename Params::ResourceTypes *>(nullptr)) <= 32,
        "constant_families holds one bit per parameter family");
    int family = 0;
    if (!constant_params)
      b.params.visit(
          [&]<typename Resource>(const auto &) __attribute__((always_inline)) {
            const auto &start = a.params.template get<Resource::KEY>();
            const auto &end = b.params.template get<Resource::KEY>();
            if (std::memcmp(&start, &end, sizeof(end)) == 0)
              constant_families |= 1u << family;
            ++family;
          });
    b.params.visit([&]<typename Resource>(const auto &) HS_HOT_INLINE {
      auto &delta = deltas.template get<Resource::KEY>();
      const auto &start = a.resources.template get<Resource::KEY>();
      const auto &end = b.resources.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        prepare_source_clocks(delta, start, end);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        delta.phase = Spacetime::ClockInterval::prepare_delta(start.phase,
                                                              end.phase, 1.0f);
        if constexpr (std::is_base_of_v<AffineParams,
                                        typename Resource::Family>)
          delta.rotation = Spacetime::ClockInterval::prepare_delta(
              start.rotation, end.rotation, math::TWO_PI_F);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE &&
                           !std::is_same_v<typename Resource::Family,
                                           PeriodicRippleParams> &&
                           !std::is_same_v<
                               typename Resource::Family,
                               SphereAnimation::BurstRippleParams>) {
        constexpr float PERIOD =
            Surface::phase_period<typename Resource::Family>();
        delta.phase = Spacetime::ClockInterval::prepare_delta(
            start.phase, end.phase, PERIOD);
      } else if constexpr (Resource::KIND == ResourceKind::TWIST) {
        delta.phi = Spacetime::ClockInterval::prepare_delta(start.phi, end.phi,
                                                            math::TWO_PI_F);
        delta.lambda = Spacetime::ClockInterval::prepare_delta(
            start.lambda, end.lambda, math::TWO_PI_F);
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusGeneratorParams>) {
        delta.angle = Spacetime::ClockInterval::prepare_delta(
            start.angle, end.angle, math::TWO_PI_F);
        delta.magnitude = end.magnitude - start.magnitude;
        auto &cache = mobius_intervals.template get<Resource::KEY>();
        const auto &a_params = a.params.template get<Resource::KEY>();
        const auto &b_params = b.params.template get<Resource::KEY>();
        cache.active = start.magnitude == end.magnitude &&
                       std::memcmp(&a_params, &b_params, sizeof(b_params)) == 0;
        if (cache.active)
          cache.base = Lens::prepare_mobius(b_params, 0.0f, end.magnitude);
      } else if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        const auto &previous =
            a.params.template get<Resource::KEY>().modulation;
        const auto &current = b.params.template get<Resource::KEY>().modulation;
        delta.previous_modulation =
            interval_modulation(previous, current, start.modulation,
                                end.modulation, phase_compatible);
      }
    });
  }

  /** @brief Samples the prepared interval, preserving both endpoints.
   * @param fraction Frame progress in [0, 1].
   * @return Interpolated context borrowing the endpoints and interval. */
  HS_FLASH_INLINE Frame at(float fraction) const & {
    fraction = sample_fraction(fraction);
    Frame frame = fraction <= 0.0f ? from : to;
    apply_sample(frame, fraction);
    return frame;
  }

  /** @brief Samples a temporary interval without borrowing its cache.
   * @param fraction Frame progress in [0, 1].
   * @return Interpolated context borrowing the endpoint resources. */
  HS_FLASH_INLINE Frame at(float fraction) const && {
    Frame frame = at(fraction);
    frame.params.visit([&]<typename Resource>(const auto &) {
      if constexpr (std::is_same_v<typename Resource::Family,
                                   MobiusGeneratorParams>)
        frame.resources.template get<Resource::KEY>().mobius_base = nullptr;
    });
    return frame;
  }

  /** @brief Writes a sampled context into storage separate from the endpoints.
   * @param frame Receives the context borrowing the endpoint resources.
   * @param fraction Frame progress in [0, 1]. */
  HS_FLASH_INLINE void at_into(Frame &frame, float fraction) const & {
    fraction = sample_fraction(fraction);
    frame = fraction <= 0.0f ? from : to;
    apply_sample(frame, fraction);
  }

  /** @brief Samples a temporary interval without retaining its cache.
   * @param frame Receives the sampled context.
   * @param fraction Frame progress in [0, 1]. */
  HS_FLASH_INLINE void at_into(Frame &frame, float fraction) const && {
    frame = std::move(*this).at(fraction);
  }

private:
  HS_HOT_INLINE float sample_fraction(float fraction) const {
    if constexpr (requires { to.source_morph; })
      if (to.source_morph.snapped)
        return 1.0f;
    if constexpr (HAS_PLANAR)
      if (to.resources.template get<"source">().snapped)
        return 1.0f;
    return fraction;
  }

  HS_FLASH_INLINE void apply_sample(Frame &frame, float fraction) const {
    if (fraction <= 0.0f) {
      to.params.visit([&]<typename Resource>(const auto &family) {
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
      frame.temporal.hue_rotation_gamut = from.hue_rotation_gamut;
      bind_mobius_intervals(frame);
      return;
    }
    if (fraction >= 1.0f) {
      frame.params.visit([&]<typename Resource>(const auto &) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>)
          frame.resources.template get<Resource::KEY>().radial_fraction = -1.0f;
      });
      bind_mobius_intervals(frame);
      return;
    }
    if constexpr (requires { frame.source_morph; }) {
      const auto &a = from.source_morph;
      const auto &b = to.source_morph;
      const auto &edge = b.from ? b : a;
      if (edge.from) {
        auto &morph = frame.source_morph;
        morph = edge;
        if constexpr (requires { frame.hex_tables; })
          frame.hex_tables = b.from ? to.hex_tables : from.hex_tables;
        morph.linear = hs::lerp(a.from ? a.linear : 0.0f,
                                b.from ? b.linear : 1.0f, fraction);
        morph.mix = hs::clamp(morph.ease(morph.linear), 0.0f, 1.0f);
      }
    }
    frame.projection_conjugate = projection.sample(fraction);
    frame.outer_conjugate = outer.sample(fraction);
    if (!constant_params) {
      int family = 0;
      frame.params.visit(
          [&]<typename Resource>(auto &value) __attribute__((always_inline)) {
            if (!(constant_families & (1u << family))) {
              const auto &a = from.params.template get<Resource::KEY>();
              const auto &b = to.params.template get<Resource::KEY>();
              if constexpr (HAS_PLANAR)
                Source::sample_planar_family<Resource, SourceSelector>(
                    value, a, b, fraction,
                    from.resources.template get<"source">().weights,
                    to.resources.template get<"source">().weights);
              else
                sample_parameter_family(value, a, b, fraction);
            }
            ++family;
          });
    }
    if constexpr (HAS_PLANAR)
      frame.resources.template get<"source">().weights =
          SourceSelector::PlanarWeights::lerp_components(
              from.resources.template get<"source">().weights,
              to.resources.template get<"source">().weights, fraction);
    frame.palette_mapping = Color::PaletteMappingWeights::lerp(
        from.palette_mapping, to.palette_mapping, fraction);
    frame.palette_oscillation_phase = Spacetime::ClockInterval::apply(
        from.palette_oscillation_phase, oscillation_delta, fraction);
    frame.params.visit([&]<typename Resource>(const auto &family) {
      auto &clock = frame.resources.template get<Resource::KEY>();
      const auto &start = from.resources.template get<Resource::KEY>();
      const auto &delta = deltas.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::ROTATION) {
        const auto &end = to.resources.template get<Resource::KEY>();
        sample_keyed_rotation(clock, start, end, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        sample_source_clocks(clock, start, delta, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        clock.phase =
            Spacetime::ClockInterval::apply(start.phase, delta.phase, fraction);
        if constexpr (std::is_base_of_v<AffineParams,
                                        typename Resource::Family>)
          clock.rotation = Spacetime::ClockInterval::apply(
              start.rotation, delta.rotation, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>) {
          const auto &end = to.resources.template get<Resource::KEY>();
          const uint32_t delta_value = end.frames - start.frames;
          SphereAnimation::rewind_wave_ages(clock.waves, delta_value, fraction);
          clock.radial_fraction = delta_value == 1 ? fraction : 1.0f;
          if (clock.radial_tables &&
              (delta_value > 1 ||
               std::memcmp(&from.params.template get<Resource::KEY>(),
                           &to.params.template get<Resource::KEY>(),
                           sizeof(family)) != 0))
            clock.radial_tables = nullptr;
        } else if constexpr (std::is_same_v<typename Resource::Family,
                                            PeriodicRippleParams>)
          clock.phase = Spacetime::sample_clock(
              start.phase, to.resources.template get<Resource::KEY>().phase,
              fraction, family.period);
        else
          clock.phase = Spacetime::ClockInterval::apply(start.phase,
                                                        delta.phase, fraction);
      } else if constexpr (Resource::KIND == ResourceKind::TWIST) {
        clock.phi =
            Spacetime::ClockInterval::apply(start.phi, delta.phi, fraction);
        clock.lambda = Spacetime::ClockInterval::apply(start.lambda,
                                                       delta.lambda, fraction);
      } else if constexpr (std::is_same_v<typename Resource::Family,
                                          MobiusGeneratorParams>) {
        clock.angle =
            Spacetime::ClockInterval::apply(start.angle, delta.angle, fraction);
        clock.magnitude = start.magnitude + delta.magnitude * fraction;
      }
    });
    bind_mobius_intervals(frame);
    bool same_sources = true;
    frame.params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        auto &resource = frame.resources.template get<Resource::KEY>();
        resource.previous_modulation =
            deltas.template get<Resource::KEY>().previous_modulation;
        if constexpr (requires { resource.phase_noise_current; }) {
          resource.phase_noise_previous =
              from.resources.template get<Resource::KEY>().phase_noise_current;
          resource.phase_noise_fraction = fraction;
          if (!same_projection_controls(from, to)) {
            resource.phase_noise_previous = {};
            resource.phase_noise_current = {};
          }
        }
        if (ColorModulation::active(family.modulation.enabled,
                                    family.modulation.mechanic) &&
            resource.modulation != resource.previous_modulation)
          same_sources = false;
      }
    });
    frame.temporal = {from.palette,
                      same_palette ? to.hue_rotation_lut
                                   : from.hue_rotation_lut,
                      same_palette && same_sources ? 1.0f : fraction,
                      same_palette, phase_compatible};
    frame.temporal.hue_rotation_gamut =
        same_palette ? to.hue_rotation_gamut : from.hue_rotation_gamut;
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

/** @brief Value-initializes an arena-owned frame context.
 * @return Value-initialized frame context. */
template <typename Frame> HS_SIZE_FLASH_MEMBER Frame initial_frame_context() {
  return {};
}

/** @brief Owned palette and hue tables for the previous frame endpoint. */
template <typename Frame, bool Hue>
  requires BorrowsColorResources<Frame>
struct FrameHistory<Frame, true, Hue> : HueRotationGamutStorage<Hue> {
  Frame previous = initial_frame_context<Frame>(); ///< Retained endpoint.
  BakedPaletteStorage palette;                     ///< Owned palette entries.
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
      if (frame.hue_rotation_lut) {
        std::memcpy(hue_rotation.data(), frame.hue_rotation_lut,
                    sizeof(hue_rotation));
        previous.hue_rotation_lut = hue_rotation.data();
      }
      if (frame.hue_rotation_gamut) {
        this->retained_hue_rotation_gamut = *frame.hue_rotation_gamut;
        previous.hue_rotation_gamut = &this->retained_hue_rotation_gamut;
      }
    }
    frame.params.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
        if constexpr (Resource::Family::DYNAMIC ||
                      Resource::Family::SOURCE !=
                          Pattern::SourceKind::DISPLACEMENT) {
          auto &owned = resources.template get<Resource::KEY>();
          auto &view = previous.resources.template get<Resource::KEY>();
          view.modulation = Pullback::Modulation::retain(
              frame.resources.template get<Resource::KEY>().modulation,
              owned.modulation, owned.table);
          view.previous_modulation = view.modulation;
          if constexpr (requires { view.phase_noise_current; })
            view.phase_noise_previous = view.phase_noise_current;
        }
      }
    });
    previous.temporal = {previous.palette, previous.hue_rotation_lut, 1.0f};
    previous.temporal.hue_rotation_gamut = previous.hue_rotation_gamut;
    initialized = true;
  }
};

} // namespace Spacetime

namespace Pullback {

/**
 * @brief Retained endpoint, column-time interval and interval palette caches
 *        of one composed effect.
 * @details The storage types are bases and members of the runtime State,
 * which every operation receives.
 * @tparam Traits The effect's `ComposedTraits`.
 */
template <typename Traits> class ComposedIntervalCache {
  using FrameState = typename Traits::FrameState;
  using Params = typename Traits::Params;
  using SpecT = typename Traits::Spec;
  static constexpr bool TEMPORAL_SAMPLING = Traits::TEMPORAL_SAMPLING;
  static constexpr bool PALETTE_LAB_CACHE = Traits::PALETTE_LAB_CACHE;
  static constexpr bool PALETTE_RGB_CACHE = Traits::PALETTE_RGB_CACHE;

public:
  /// Retained previous endpoint and its owned color resources.
  using History = Spacetime::FrameHistory<FrameState, TEMPORAL_SAMPLING,
                                          Traits::HUE_ROTATION_LUT>;
  /// Ping-pong interval palette bases.
  using LabLuts =
      OptionalPaletteLabLut<PALETTE_LAB_CACHE, SpecT::PALETTE_LAB_VALUES>;
  /// Ping-pong final hue tables over the interval bases.
  using RgbLuts = OptionalPaletteRgbLut<PALETTE_RGB_CACHE,
                                        Traits::PALETTE_RGB_MIDPOINT_ONLY>;
  /// Interval sampler between the retained and current endpoints.
  using Interval = std::conditional_t<TEMPORAL_SAMPLING,
                                      std::optional<FrameInterval<FrameState>>,
                                      std::monostate>;

  /** @brief Retains the first endpoint when nothing is retained yet.
   * @tparam Storage Runtime State deriving from `History`.
   * @param storage Arena-owned interval storage.
   * @param current This draw's frame. */
  template <typename Storage>
  static HS_O3_FN HS_HOT_INLINE void begin(Storage &storage,
                                           const FrameState &current) {
    History &history = storage;
    if (!history.initialized)
      history.capture(current, persistent_arena);
  }

  /**
   * @brief Bakes the interval palette caches when both endpoints qualify and
   *        opens the interval sampler.
   * @tparam Storage Runtime State deriving from `History`, `LabLuts` and
   *         `RgbLuts`, with an `interval` member.
   * @tparam Eligible Predicate admitting one endpoint's controls to caching.
   * @param storage Arena-owned interval storage.
   * @param current This draw's frame; receives the baked base cache.
   * @param eligible Endpoint control predicate.
   * @return The sampler from the retained to the current endpoint.
   */
  template <typename Storage, typename Eligible>
  static HS_O3_FN HS_HOT_INLINE const FrameInterval<FrameState> &
  bake(Storage &storage, FrameState &current, Eligible &&eligible) {
    const History &history = storage;
    if constexpr (PALETTE_LAB_CACHE) {
      if (eligible(current.params) && eligible(history.previous.params) &&
          Color::same_palette_controls(
              history.previous.params.template get<"color">(),
              current.params.template get<"color">())) {
        HS_PROFILE(fx_palette_base_bake);
        storage.palette_lab_slot ^= 1;
        auto &cache = storage.palette_lab_luts[storage.palette_lab_slot];
        cache.bake(*history.previous.palette, *current.palette);
        if constexpr (PALETTE_RGB_CACHE) {
          bool hue_only = true;
          current.params.visit([&]<typename Resource>(const auto &family) {
            if constexpr (Resource::KIND == ResourceKind::MODULATION)
              if (!ColorModulation::hue_only(family.modulation.enabled,
                                             family.modulation.mechanic))
                hue_only = false;
          });
          if (hue_only) {
            HS_PROFILE(fx_palette_rgb_bake);
            storage.palette_rgb_luts[storage.palette_lab_slot].bake(
                cache, &storage.palette_lab_luts[storage.palette_lab_slot ^ 1],
                &storage.palette_rgb_luts[storage.palette_lab_slot ^ 1]);
          }
        }
        current.palette_lab_lut = &cache;
      }
    }
    storage.interval.emplace(history.previous, current);
    return *storage.interval;
  }

  /** @brief Closes the interval and retains this draw's endpoint.
   * @tparam Storage Runtime State deriving from `History`, with an
   *         `interval` member.
   * @param storage Arena-owned interval storage.
   * @param current This draw's frame. */
  template <typename Storage>
  static HS_O3_FN HS_HOT_INLINE void capture(Storage &storage,
                                             const FrameState &current) {
    storage.interval.reset();
    static_cast<History &>(storage).capture(current, persistent_arena);
  }
};

} // namespace Pullback
