/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <type_traits>
#include "core/math/noise_field.h"
#include "core/render/pullback/composed_descriptors.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/spacetime.h"

/** @file composed_clocks.h
 * @brief `ComposedClocks`, the per-resource clocks of a composed effect.
 */

namespace Pullback {

/**
 * @brief Per-resource phase clocks, rotation walks, burst ripples, noise
 *        generators and chamber noise caches of one composed effect.
 * @details Clocks live in this object; the arena-owned per-resource state is
 * `Storage`, held by the runtime State and bound by seed().
 * @tparam Traits The effect's `ComposedTraits`.
 */
template <typename Traits> class ComposedClocks {
  using Params = typename Traits::Params;
  using FrameState = typename Traits::FrameState;
  using SpecT = typename Traits::Spec;
  static constexpr bool AUTHORED_TOPOLOGY = Traits::AUTHORED_TOPOLOGY;

public:
  /** @brief Owned per-resource state with optional ripple transfer curves.
   * @tparam Resource Keyed resource. */
  template <typename Resource>
  struct Owned
      : ComposedDetail::ResourceNoise<Resource>,
        ComposedDetail::ResourcePhaseNoiseCache<
            Resource, SpecT::PHASE_NOISE_CACHE &&
                          ComposedDetail::RESOURCE_DIRECT_NOISE<Resource>>,
        ComposedDetail::ResourceRippleRadialTable<
            Resource, SpecT::RIPPLE_RADIAL_TABLE &&
                          std::is_same_v<typename Resource::Family,
                                         SphereAnimation::BurstRippleParams>> {
  };
  /// Arena-owned per-resource state.
  using Storage =
      ComposedDetail::ResourceStorage<Owned, typename Params::ResourceTypes>;

  /** @brief Binds the effect's live parameters.
   * @param params Parameter set the clocks read for their lifetime. */
  HS_HOT_INLINE explicit ComposedClocks(const Params &params)
      : params(&params) {}

  /** @brief Clock storage of one resource.
   * @tparam Key Resource key; must exist.
   * @return Mutable clock state. */
  template <ResourceKey Key> auto &get() { return clocks.template get<Key>(); }

  /**
   * @brief Binds the arena-owned per-resource state and seeds rotation
   *        walks, ripple schedulers and noise generators from each
   *        resource's instance identity.
   * @tparam Effect Supplies `WARP_NOISE_SEED`, `SOURCE_NOISE_SEED`,
   *         `SURFACE_NOISE_SEED` and `authored_instance_label`.
   * @param storage Arena-owned per-resource state.
   */
  template <typename Effect>
  HS_O3_FN HS_HOT_INLINE void seed(Storage &storage) {
    owned = &storage;
    params->visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == ResourceKind::ROTATION) {
        constexpr auto INSTANCE =
            AUTHORED_TOPOLOGY
                ? Effect::authored_instance_label(Resource::KEY.view())
                : Resource::KEY.view();
        constexpr uint32_t SEED = instance_hash(INSTANCE, "sphere.rotate.v2");
        SphereAnimation::init_walk(
            clocks.template get<Resource::KEY>().rotation_walk,
            static_cast<int32_t>(SEED));
      }
      if constexpr (std::is_same_v<typename Resource::Family,
                                   SphereAnimation::BurstRippleParams>) {
        constexpr auto INSTANCE =
            AUTHORED_TOPOLOGY
                ? Effect::authored_instance_label(Resource::KEY.view())
                : Resource::KEY.view();
        owned->template get<Resource::KEY>().ripples.seed =
            instance_hash(INSTANCE, "sphere.displace.ripple.v2");
      }
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>) {
        constexpr int32_t SEED = Resource::KIND == ResourceKind::WARP
                                     ? Effect::WARP_NOISE_SEED
                                 : Resource::KIND == ResourceKind::SOURCE
                                     ? Effect::SOURCE_NOISE_SEED
                                     : Effect::SURFACE_NOISE_SEED;
        constexpr int32_t INSTANCE_SEED = [] {
          if constexpr (AUTHORED_TOPOLOGY) {
            constexpr auto OPERATOR =
                ComposedDetail::noise_operator_id<typename Resource::Family>();
            static_assert(
                !OPERATOR.empty(),
                "authored noise resource requires a catalog operator identity");
            return static_cast<int32_t>(instance_hash(
                Effect::authored_instance_label(Resource::KEY.view()),
                OPERATOR));
          } else if constexpr (ComposedDetail::fixed_noise_qualified<Resource>(
                                   typename Params::ResourceTypes{})) {
            constexpr uint32_t HASH = fnv1a(Resource::KEY.view());
            return static_cast<int32_t>(static_cast<uint32_t>(SEED) ^ HASH);
          } else {
            return SEED;
          }
        }();
        Pullback::init_effect_noise(owned->template get<Resource::KEY>().noise,
                                    INSTANCE_SEED);
      }
    });
  }

  /**
   * @brief Steps every phase clock the effect's parameter families define.
   * @details Each clock is compiled in only when its field exists. Affine
   * rotations accumulate independently for each affine warp.
   */
  HS_O3_FN HS_HOT_INLINE void advance() {
    params->visit([&]<typename Resource>(const auto &family) HS_HOT_INLINE {
      auto &clock = clocks.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::ROTATION) {
        SphereAnimation::advance_walk(clock.rotation_walk, family.wander,
                                      family.spin_speed);
        clock.wander = clock.rotation_walk.wander;
        clock.spin_phase = clock.rotation_walk.spin_phase;
        clock.conjugate =
            (math::make_rotation(math::Y_AXIS, clock.spin_phase) * clock.wander)
                .conjugate();
      } else if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        Pattern::advance_clocks(family, clock.primary, clock.secondary,
                                clock.angle);
        if constexpr (requires { family.noise_time_rate; })
          clock.noise_time =
              Pattern::advance_noise_phase(clock.noise_time, family);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     SphereAnimation::BurstRippleParams>)
          SphereAnimation::advance_ripple(
              owned->template get<Resource::KEY>().ripples, family);
        else if constexpr (std::is_same_v<typename Resource::Family,
                                          PeriodicRippleParams>)
          Surface::advance_ripple_phase(clock.phase, family);
        else if constexpr (constexpr float PERIOD = Surface::phase_period<
                               typename Resource::Family>();
                           PERIOD != 1.0f) {
          if constexpr (requires { family.drift; })
            clock.phase =
                Surface::advance_field_phase(clock.phase, family.drift, PERIOD);
          else
            clock.phase =
                Surface::advance_field_phase(clock.phase, family.speed, PERIOD);
        } else
          clock.phase = math::wrap_t(clock.phase + family.speed);
      } else if constexpr (Resource::KIND == ResourceKind::TWIST) {
        Lens::advance_twist(clock, family);
      } else if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        Modulation::advance(clock.modulation_clock, family.modulation);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        if constexpr (std::is_base_of_v<AffineParams,
                                        typename Resource::Family>)
          Warp::advance_affine_rotation(clock.rotation, family);
        clock.phase = math::wrap_t(clock.phase + family.speed);
      }
    });
  }

  /**
   * @brief Writes each resource's frame view: clocks, ripples and radial
   *        tables, noise generators, modulation sources and chamber noise.
   * @param frame Frame whose `params`, `resources` and
   *        `projection_conjugate` are already set.
   */
  HS_O3_FN HS_HOT_INLINE void prepare_views(FrameState &frame) {
    params->visit([&]<typename Resource>(const auto &) {
      auto &resource = frame.resources.template get<Resource::KEY>();
      if constexpr (std::is_same_v<typename Resource::Family,
                                   SphereAnimation::BurstRippleParams>) {
        auto &state = owned->template get<Resource::KEY>();
        static_cast<SphereAnimation::RippleState &>(resource) = state.ripples;
        if constexpr (SpecT::RIPPLE_RADIAL_TABLE) {
          SphereAnimation::prepare_radial_tables(
              frame.params.template get<Resource::KEY>(), state.ripples,
              state.radial_tables, true);
          resource.radial_tables = &state.radial_tables;
        }
      } else
        static_cast<ComposedDetail::ResourceClock<Resource> &>(resource) =
            clocks.template get<Resource::KEY>();
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>)
        resource.noise = &owned->template get<Resource::KEY>().noise;
      if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        const auto &config = params->template get<Resource::KEY>().modulation;
        auto &clock = clocks.template get<Resource::KEY>().modulation_clock;
        Modulation::synchronize(clock, config);
        resource.modulation_clock = clock;
        if constexpr (Resource::Family::DYNAMIC ||
                      Resource::Family::SOURCE !=
                          Pattern::SourceKind::DISPLACEMENT)
          resource.modulation =
              &owned->template get<Resource::KEY>()
                   .modulation_resources.prepare(config, clock);
        resource.previous_modulation = resource.modulation;
        if constexpr (SpecT::PHASE_NOISE_CACHE &&
                      ComposedDetail::RESOURCE_DIRECT_NOISE<Resource>) {
          if (!Spacetime::flat_surface()) {
            auto &state = owned->template get<Resource::KEY>();
            resource.phase_noise_current = state.phase_noise.prepare(
                *resource.modulation, frame.projection_conjugate);
            resource.phase_noise_previous = resource.phase_noise_current;
          }
        }
      }
    });
  }

private:
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceRuntimeClock,
                                  typename Params::ResourceTypes>
      clocks; ///< Member clocks, excluding arena-owned burst schedulers.
  const Params *params;     ///< Live parameters.
  Storage *owned = nullptr; ///< Arena-owned per-resource state.
};

} // namespace Pullback
