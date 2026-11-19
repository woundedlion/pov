/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <bit>
#include <optional>
#include "core/render/pullback/composed_parameters.h"
#include "core/render/pullback/planar_source.h"

/** @file composed_planar.h
 * @brief Compiled planar source selection, field interpolation and sampling.
 */

namespace Pullback::Source {

namespace PlanarDetail {

template <PlanarFamily Family> struct FamilyBinding;
template <> struct FamilyBinding<PlanarFamily::GRID> {
  using Params = GridSourceParams;
  static constexpr ResourceKey KEY = "source_grid";
  template <typename Provider> using Policy = Grid<Provider>;
};
template <> struct FamilyBinding<PlanarFamily::TWIN_WAVE> {
  using Params = TwinWaveSourceParams;
  static constexpr ResourceKey KEY = "source_twin_wave";
  template <typename Provider> using Policy = TwinWave<Provider>;
};
template <> struct FamilyBinding<PlanarFamily::BANDS> {
  using Params = BandsSourceParams;
  static constexpr ResourceKey KEY = "source_bands";
  template <typename Provider> using Policy = Bands<Provider>;
};
template <> struct FamilyBinding<PlanarFamily::SPIRAL> {
  using Params = SpiralSourceParams;
  static constexpr ResourceKey KEY = "source_spiral";
  template <typename Provider> using Policy = Spiral<Provider>;
};
template <> struct FamilyBinding<PlanarFamily::TESSELLATION> {
  using Params = TessellationSourceParams;
  static constexpr ResourceKey KEY = "source_tessellation";
};
template <> struct FamilyBinding<PlanarFamily::PROJECTED_NOISE> {
  using Params = ProjectedNoiseSourceParams;
  static constexpr ResourceKey KEY = "source_projected_noise";
};

template <PlanarKind Kind>
struct KindBinding : FamilyBinding<planar_family(Kind)> {};

template <PlanarKind Kind>
  requires(planar_family(Kind) == PlanarFamily::TESSELLATION)
struct KindBinding<Kind> : FamilyBinding<PlanarFamily::TESSELLATION> {
  static constexpr auto INDEX =
      static_cast<unsigned>(Kind) -
      static_cast<unsigned>(PlanarKind::TESSELLATION_TRIANGULAR_CELLS);
  static constexpr auto SHAPE =
      static_cast<Pattern::TessellationKind>(INDEX / 2);
  static constexpr auto FILL =
      static_cast<Pattern::TessellationFill>(INDEX % 2);

  template <typename Provider>
  using Policy = Tessellation<Provider, SHAPE, FILL>;
};

template <PlanarKind Kind>
  requires(planar_family(Kind) == PlanarFamily::PROJECTED_NOISE)
struct KindBinding<Kind> : FamilyBinding<PlanarFamily::PROJECTED_NOISE> {
  static constexpr auto BASIS = static_cast<math::NoiseBasis>(
      static_cast<unsigned>(Kind) -
      static_cast<unsigned>(PlanarKind::PROJECTED_NOISE_SIMPLEX));
  template <typename Provider> using Policy = ProjectedNoise<Provider, BASIS>;
};

template <PlanarKind Kind> consteval bool has_lattice() {
  using Row = KindBinding<Kind>;
  return requires { Row::SHAPE; } ||
         requires(const typename Row::Params &params) {
           translation_lattice(params);
         };
}

template <PlanarKind Kind> consteval bool rectangular_lattice() {
  if constexpr (requires { KindBinding<Kind>::SHAPE; })
    return KindBinding<Kind>::SHAPE == Pattern::TessellationKind::SQUARE;
  return Kind == PlanarKind::GRID;
}

template <PlanarKind Kind, typename Params>
HS_FLASH_INLINE TranslationLattice lattice_for(const Params &params) {
  using Row = KindBinding<Kind>;
  if constexpr (has_lattice<Kind>()) {
    const auto &family = params.template get<Row::KEY>();
    if constexpr (requires { Row::SHAPE; })
      return tessellation_lattice(Row::SHAPE, family.cell_scale);
    else
      return translation_lattice(family);
  }
  return {};
}

template <typename Resource> consteval uint16_t family_mask() {
  if constexpr (Resource::KEY.view() == "source_grid")
    return live_values(PlanarKind::GRID);
  else if constexpr (Resource::KEY.view() == "source_twin_wave")
    return live_values(PlanarKind::TWIN_WAVE);
  else if constexpr (Resource::KEY.view() == "source_bands")
    return live_values(PlanarKind::BANDS);
  else if constexpr (Resource::KEY.view() == "source_spiral")
    return live_values(PlanarKind::SPIRAL);
  else if constexpr (Resource::KEY.view() == "source_tessellation")
    return PLANAR_TESSELLATION_CELLS | PLANAR_TESSELLATION_EDGES;
  else if constexpr (Resource::KEY.view() == "source_projected_noise")
    return PLANAR_FAMILY_MASKS[static_cast<size_t>(
        PlanarFamily::PROJECTED_NOISE)];
  else
    return 0;
}

} // namespace PlanarDetail

/** @brief Preset selector for one composed set of planar kinds.
 * @tparam Kinds Planar kernels included in the effect.
 */
template <PlanarKind... Kinds> struct PlanarSelectorParams {
  /// Authored source kinds, in weight-vector order.
  static constexpr std::array<PlanarKind, sizeof...(Kinds)> COMPILED_KINDS{
      Kinds...};
  static_assert(sizeof...(Kinds) > 0);
  static_assert(((static_cast<size_t>(Kinds) < PLANAR_KIND_COUNT) && ...));
  static_assert(COMPILED_KINDS.size() <= PLANAR_KIND_COUNT);
  /// Bit mask over the catalog's source enum.
  static constexpr uint16_t KIND_MASK = live_values(Kinds...);
  static_assert(std::popcount(KIND_MASK) == sizeof...(Kinds));
  using PlanarWeights = Color::EnumWeights<sizeof...(Kinds)>; ///< Column mix.
  /// Descriptors for the numeric controls in this family.
  static constexpr std::array<Field<PlanarSelectorParams>, 0> FIELDS{};
  /// Any compiled source reports a translation lattice.
  static constexpr bool HAS_LATTICE =
      (PlanarDetail::has_lattice<Kinds>() || ...);
  /// Any compiled source uses the affine's fallback basis.
  static constexpr bool FALLBACK_LATTICE =
      (!PlanarDetail::has_lattice<Kinds>() || ...);
  /// Path reduction selected from the set of compiled lattices.
  static constexpr auto REPEAT_REDUCTION =
      (PlanarDetail::rectangular_lattice<Kinds>() && ...)
          ? Warp::RepeatReduction::RECTANGULAR
      : HAS_LATTICE ? Warp::RepeatReduction::GENERAL
                    : Warp::RepeatReduction::UNWRAPPED;

  PlanarKind source = COMPILED_KINDS.front(); ///< Preset source kind.

  /** @brief Finds the compact index of a catalog kind.
   * @param kind Authored source kind.
   * @return Weight-vector index or its size when absent. */
  static constexpr size_t index(PlanarKind kind) {
    size_t result = 0;
    while (result < COMPILED_KINDS.size() && COMPILED_KINDS[result] != kind)
      ++result;
    return result;
  }

  /** @brief Constructs the authored kind's exact weights.
   * @return Unit weight on the authored source kind. */
  HS_FLASH_INLINE PlanarWeights weights() const {
    const size_t selected = index(source);
    HS_CHECK(selected < COMPILED_KINDS.size(),
             "planar source kind is not compiled");
    return PlanarWeights::single(selected);
  }

  /** @brief Reports the selected source's lattice.
   * @tparam Params Composed parameter set.
   * @param params Parameters of the full source.
   * @return The selected basis, or a zero basis. */
  template <typename Params>
  HS_FLASH_INLINE TranslationLattice lattice(const Params &params) const {
    TranslationLattice result{};
    ((source == Kinds ? void(result = PlanarDetail::lattice_for<Kinds>(params))
                      : void()),
     ...);
    return result;
  }
};

/** @brief Validates membership in the effect's compiled source set.
 * @tparam Kinds Compiled kernels.
 * @param selector Preset source selector.
 * @return True when the selected kind is compiled. */
template <PlanarKind... Kinds>
constexpr bool valid(const PlanarSelectorParams<Kinds...> &selector) {
  return ((selector.source == Kinds) || ...);
}

namespace PlanarDetail {

template <typename B, PlanarKind Kind>
using ProviderFor = SourceProvider<B, typename KindBinding<Kind>::Params,
                                   KindBinding<Kind>::KEY>;

enum class Endpoint { CURRENT, FROM, TO };

template <typename B, PlanarKind Kind, Endpoint Side>
struct EndpointProvider : ProviderFor<B, Kind> {
  using FrameState = typename B::FrameState;

  static HS_HOT_INLINE const auto &params(const FrameState &frame) {
    const auto *endpoint = Side == Endpoint::FROM ? frame.source_morph.from
                                                  : frame.source_morph.to;
    return endpoint->template get<KindBinding<Kind>::KEY>();
  }

  static HS_HOT_INLINE float noise_scale(const FrameState &frame) {
    return params(frame).noise_scale;
  }

  static HS_HOT_INLINE float noise_contrast(const FrameState &frame) {
    return params(frame).noise_contrast;
  }
};

template <typename B, PlanarKind Kind, Endpoint Side = Endpoint::CURRENT>
using PolicyFor = typename KindBinding<Kind>::template Policy<
    std::conditional_t<Side == Endpoint::CURRENT, ProviderFor<B, Kind>,
                       EndpointProvider<B, Kind, Side>>>;

template <typename Selector>
HS_FLASH_INLINE bool live(const typename Selector::PlanarWeights &weights,
                          uint16_t mask) {
  for (size_t i = 0; i < Selector::COMPILED_KINDS.size(); ++i)
    if ((mask & (1U << static_cast<unsigned>(Selector::COMPILED_KINDS[i]))) &&
        weights.values[i] > 0.0f)
      return true;
  return false;
}

template <typename Family>
HS_FLASH_INLINE uint16_t gated_mask(const Field<Family> &field,
                                    const Family &params, uint16_t mask) {
  if constexpr (std::is_same_v<Family, SpiralSourceParams>) {
    if (field.topology_gate.field &&
        !((1U << params.coordinate_system) & field.topology_gate.values))
      return 0;
  } else if constexpr (std::is_same_v<Family, TessellationSourceParams>) {
    if (field.topology_gate.field) {
      if (field.topology_gate.values == TESSELLATION_CELLS_GATE.values)
        return mask & PLANAR_TESSELLATION_CELLS;
      return mask & PLANAR_TESSELLATION_EDGES;
    }
  }
  return mask;
}

template <typename Family>
constexpr bool different_fields(const Family &a, const Family &b,
                                uint16_t mask) {
  if constexpr (requires { Family::TOPOLOGY; })
    for (const auto &field : Family::TOPOLOGY)
      if (a.*field.member != b.*field.member)
        return true;
  for (const auto &field : Family::FIELDS)
    if ((mask == 0 || gated_mask(field, a, mask)) &&
        a.*field.member != b.*field.member)
      return true;
  return false;
}

} // namespace PlanarDetail

/** @brief Tests whether the source endpoints describe different scalar fields.
 * @tparam Params Composed parameter set.
 * @param from Outgoing parameters.
 * @param to Incoming parameters.
 * @return True when the kind or a live source control differs. */
template <typename Params>
HS_FLASH_MEMBER bool different_endpoint_fields(const Params &from,
                                               const Params &to) {
  const auto &a = from.template get<"source">();
  const auto &b = to.template get<"source">();
  if constexpr (requires {
                  typename std::remove_cvref_t<decltype(a)>::PlanarWeights;
                }) {
    if (a.source != b.source)
      return true;
    const uint16_t kind = 1u << static_cast<unsigned>(a.source);
    bool different = false;
    from.visit([&]<typename Resource>(const auto &family) {
      constexpr uint16_t MASK = PlanarDetail::family_mask<Resource>();
      if constexpr (MASK != 0)
        if (kind & MASK)
          different = PlanarDetail::different_fields(
              family, to.template get<Resource::KEY>(), kind);
    });
    return different;
  } else
    return PlanarDetail::different_fields(a, b, 0);
}

/** @brief Samples a family, holding controls live at only one endpoint.
 * @tparam Resource Keyed parameter resource.
 * @tparam Selector Compiled selector type.
 * @param result Receives the sampled family.
 * @param from Outgoing family controls.
 * @param to Incoming family controls.
 * @param progress Interpolation progress.
 * @param a Outgoing source weights.
 * @param b Incoming source weights. */
template <typename Resource, typename Selector>
HS_FLASH_INLINE void
sample_planar_family(typename Resource::Family &result,
                     const typename Resource::Family &from,
                     const typename Resource::Family &to, float progress,
                     const typename Selector::PlanarWeights &a,
                     const typename Selector::PlanarWeights &b) {
  constexpr uint16_t FAMILY = PlanarDetail::family_mask<Resource>();
  if constexpr (FAMILY != 0) {
    const bool from_live = PlanarDetail::live<Selector>(a, FAMILY);
    const bool to_live = PlanarDetail::live<Selector>(b, FAMILY);
    if (from_live != to_live) {
      result = from_live ? from : to;
      return;
    }
  }
  sample_parameter_family(result, from, to, progress);
  if constexpr (std::is_same_v<typename Resource::Family, SpiralSourceParams> ||
                std::is_same_v<typename Resource::Family,
                               TessellationSourceParams>) {
    if constexpr (FAMILY != 0)
      for (const auto &field : Resource::Family::FIELDS) {
        if (!field.topology_gate.field)
          continue;
        const bool from_live = PlanarDetail::live<Selector>(
            a, PlanarDetail::gated_mask(field, from, FAMILY));
        const bool to_live = PlanarDetail::live<Selector>(
            b, PlanarDetail::gated_mask(field, to, FAMILY));
        if (from_live != to_live)
          result.*field.member = (from_live ? from : to).*field.member;
      }
  }
}

/** @brief Composed signed-field blend over a compiled set of planar kernels.
 * @tparam B Frame binding.
 * @tparam Kinds Planar source kernels compiled into the effect. */
template <typename B, PlanarKind... Kinds>
struct Planar
    : Detail::CombinedApproximation<PlanarDetail::PolicyFor<B, Kinds>...> {
  using FrameState = typename B::FrameState;       ///< Effect frame state.
  using Selector = PlanarSelectorParams<Kinds...>; ///< Keyed source selector.
  /// Kernel inputs indexed in the selector's weight order.
  using Kernels = std::tuple<
      Detail::PolicyPrepared<PlanarDetail::PolicyFor<B, Kinds>, FrameState>...>;
  /** @brief Prepared kernel inputs for the active endpoint or morph pair. */
  struct Prepared {
    Kernels first;                 ///< Current or outgoing inputs.
    std::optional<Kernels> second; ///< Incoming inputs during a morph.
  };

  /** @brief Prepares kernels that have positive weight.
   * @param frame Column frame.
   * @return Prepared inputs to the selected kernels. */
  static HS_HOT_FLASH_MEMBER Prepared prepare(const FrameState &frame) {
    using PlanarDetail::Endpoint;
    const auto &morph = frame.source_morph;
    if (morph.mixed())
      return {prepare_endpoint<Endpoint::FROM>(
                  frame, morph.from->template get<"source">().source),
              prepare_endpoint<Endpoint::TO>(
                  frame, morph.to->template get<"source">().source)};
    const auto &weights = frame.resources.template get<"source">().weights;
    return {Kernels{(
                weights.values[Selector::index(Kinds)] > 0.0f
                    ? Detail::prepare_policy<PlanarDetail::PolicyFor<B, Kinds>>(
                          frame)
                    : Detail::PolicyPrepared<PlanarDetail::PolicyFor<B, Kinds>,
                                             FrameState>{})...},
            std::nullopt};
  }

  /** @brief Samples the weighted signed field.
   * @param input Projected coordinate.
   * @param frame Column frame.
   * @param prepared Per-kind prepared inputs.
   * @return Weighted source sample. */
  static HS_HOT_INLINE float sample(const PlaneSample &input,
                                    const FrameState &frame,
                                    const Prepared &prepared) {
    if (prepared.second)
      return sample_endpoints(input, frame, prepared);
    if (frame.resources.template get<"source">().weights.exact == 0)
      return sample_kind<Selector::COMPILED_KINDS.front()>(input, frame,
                                                           prepared.first);
    return sample_weighted(input, frame, prepared.first);
  }

private:
  template <PlanarDetail::Endpoint Side>
  static HS_FLASH_INLINE Kernels prepare_endpoint(const FrameState &frame,
                                                  PlanarKind kind) {
    return {
        (kind == Kinds
             ? Detail::prepare_policy<PlanarDetail::PolicyFor<B, Kinds, Side>>(
                   frame)
             : Detail::PolicyPrepared<PlanarDetail::PolicyFor<B, Kinds, Side>,
                                      FrameState>{})...};
  }

  template <PlanarKind Kind,
            PlanarDetail::Endpoint Side = PlanarDetail::Endpoint::CURRENT>
  static HS_HOT_INLINE float sample_kind(const PlaneSample &input,
                                         const FrameState &frame,
                                         const Kernels &prepared) {
    using Policy = PlanarDetail::PolicyFor<B, Kind, Side>;
    if constexpr (Detail::PolicyPrepares<Policy, FrameState>)
      return Policy::sample(input, frame,
                            std::get<Selector::index(Kind)>(prepared));
    else
      return Policy::sample(input, frame);
  }

  template <PlanarDetail::Endpoint Side>
  static HS_HOT_INLINE float
  sample_endpoint(const PlaneSample &input, const FrameState &frame,
                  const Kernels &prepared, PlanarKind kind) {
    float value = 0.0f;
    ((kind == Kinds
          ? void(value = sample_kind<Kinds, Side>(input, frame, prepared))
          : void()),
     ...);
    return value;
  }

  static HS_HOT_FLASH_MEMBER float sample_endpoints(const PlaneSample &input,
                                                    const FrameState &frame,
                                                    const Prepared &prepared) {
    using PlanarDetail::Endpoint;
    const auto &morph = frame.source_morph;
    const float from = sample_endpoint<Endpoint::FROM>(
        input, frame, prepared.first,
        morph.from->template get<"source">().source);
    const float to = sample_endpoint<Endpoint::TO>(
        input, frame, *prepared.second,
        morph.to->template get<"source">().source);
    return hs::lerp(from, to, morph.mix);
  }

  static HS_HOT_FLASH_MEMBER float sample_weighted(const PlaneSample &input,
                                                   const FrameState &frame,
                                                   const Kernels &prepared) {
    const auto &weights = frame.resources.template get<"source">().weights;
    float value = 0.0f;
    const auto add = [&]<PlanarKind Kind>() HS_HOT_INLINE {
      constexpr size_t INDEX = Selector::index(Kind);
      const float weight = weights.values[INDEX];
      if (weight > 0.0f) {
        const float field = sample_kind<Kind>(input, frame, prepared);
        if (weights.exact == INDEX)
          value = field;
        else
          value += weight * field;
      }
    };
    (add.template operator()<Kinds>(), ...);
    return value;
  }
};

/** @brief Fixed tessellation whose preset edges blend sampled fields.
 * @tparam State Source provider.
 * @tparam Kind Cell shape.
 * @tparam Fill Cell primitives or edge band. */
template <typename State, Pattern::TessellationKind Kind,
          Pattern::TessellationFill Fill>
struct EndpointTessellation : Tessellation<State, Kind, Fill> {
  using FrameState = typename State::FrameState; ///< Effect frame state.

  /** @brief Samples the live field or its two preset endpoints.
   * @param input Working plane coordinate.
   * @param frame Column frame.
   * @return Signed scalar field value. */
  static HS_HOT_INLINE float sample(const PlaneSample &input,
                                    const FrameState &frame) {
    if (frame.source_morph.mixed())
      return sample_endpoints(input, frame);
    return Tessellation<State, Kind, Fill>::sample(input, frame);
  }

private:
  static HS_HOT_FLASH_MEMBER float sample_endpoints(const PlaneSample &input,
                                                    const FrameState &frame) {
    const auto &morph = frame.source_morph;
    const Pattern::PlanePosition position{input.coords};
    const auto sample = [&](const auto &params) {
      return Pattern::tessellation<Kind, Fill>(
          position, source_kernel_params<TessellationControls<Fill>>(
                        params.template get<"source">()));
    };
    return hs::lerp(sample(*morph.from), sample(*morph.to), morph.mix);
  }
};

} // namespace Pullback::Source

namespace Pullback::ComposedDetail {

/** @brief Keeps the fixed tessellation provider's existing resource key. */
template <typename Provider, Pattern::TessellationKind Kind,
          Pattern::TessellationFill Fill>
struct PolicyResources<Source::EndpointTessellation<Provider, Kind, Fill>>
    : PolicyResources<Provider> {};

/** @brief Collects a selector and one keyed resource per compiled family.
 * @tparam B Frame binding.
 * @tparam Kinds Planar source kernels compiled into the effect. */
template <typename B, Source::PlanarKind... Kinds>
struct PolicyResources<Source::Planar<B, Kinds...>> {
  template <Source::PlanarKind Kind>
  using Row = Source::PlanarDetail::KindBinding<Kind>;
  using Type = typename MergeAll<
      ResourceList<
          ParameterResource<"source", Source::PlanarSelectorParams<Kinds...>,
                            ResourceKind::VALUE>>,
      ResourceList<
          ParameterResource<Row<Kinds>::KEY, typename Row<Kinds>::Params,
                            ResourceKind::SOURCE>>...>::Type;
};

} // namespace Pullback::ComposedDetail
