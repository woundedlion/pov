/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <algorithm>
#include <bit>

#include "core/render/pullback/operators/sample.h"
#include "core/render/pullback/planar_source.h"

/** @file planar.h
 * @brief Interpreter model for the preset-selectable planar source. */

namespace Pullback::Interp::Op {

namespace PlanarDetail {

template <typename FieldT, size_t N>
consteval auto prefixed_ids(std::string_view prefix,
                            const std::array<FieldT, N> &fields) {
  std::array<std::array<char, 48>, N> ids{};
  for (size_t index = 0; index < N; ++index) {
    const std::string_view id = fields[index].id;
    std::copy(prefix.begin(), prefix.end(), ids[index].begin());
    ids[index][prefix.size()] = '-';
    std::copy(id.begin(), id.end(), ids[index].begin() + prefix.size() + 1);
  }
  return ids;
}

template <Source::PlanarFamily Family, typename Params>
inline constexpr auto FIELD_IDS = prefixed_ids(
    Source::PLANAR_FAMILY_IDS[static_cast<size_t>(Family)], Params::FIELDS);

template <size_t... N>
consteval auto join_schema(const std::array<ParamFieldInfo, N> &...families) {
  std::array<ParamFieldInfo, (N + ...)> out{};
  auto next = out.begin();
  ((next = std::copy(families.begin(), families.end(), next)), ...);
  return out;
}

inline constexpr const char SPIRAL_COORDINATE_SYSTEM[] =
    "spiral-coordinate-system";

template <Source::PlanarFamily Family, typename Params>
consteval auto family_schema() {
  constexpr Params DEFAULTS{};
  std::array<ParamFieldInfo, Params::FIELDS.size()> out{};
  for (size_t index = 0; index < out.size(); ++index) {
    const auto &field = Params::FIELDS[index];
    auto gate = field.topology_gate;
    if constexpr (Family == Source::PlanarFamily::SPIRAL) {
      if (gate.field != nullptr)
        gate.field = SPIRAL_COORDINATE_SYSTEM;
    }
    if constexpr (Family == Source::PlanarFamily::TESSELLATION) {
      if (gate.field != nullptr)
        gate = {"source", gate.values == Source::TESSELLATION_CELLS_GATE.values
                              ? Source::PLANAR_TESSELLATION_CELLS
                              : Source::PLANAR_TESSELLATION_EDGES};
    }
    out[index] = {FIELD_IDS<Family, Params>[index].data(),
                  field.name,
                  field.min,
                  field.max,
                  DEFAULTS.*(field.member),
                  field.curve,
                  false,
                  0,
                  0,
                  nullptr,
                  gate.field,
                  gate.values};
  }
  return out;
}

constexpr auto tessellation_index(Source::PlanarKind kind) {
  return static_cast<uint8_t>(kind) -
         static_cast<uint8_t>(
             Source::PlanarKind::TESSELLATION_TRIANGULAR_CELLS);
}

constexpr Pattern::TessellationKind tessellation_kind(Source::PlanarKind kind) {
  return static_cast<Pattern::TessellationKind>(tessellation_index(kind) / 2);
}

constexpr Pattern::TessellationFill tessellation_fill(Source::PlanarKind kind) {
  return static_cast<Pattern::TessellationFill>(tessellation_index(kind) % 2);
}

constexpr math::NoiseBasis noise_basis(Source::PlanarKind kind) {
  return static_cast<math::NoiseBasis>(
      static_cast<uint8_t>(kind) -
      static_cast<uint8_t>(Source::PlanarKind::PROJECTED_NOISE_SIMPLEX));
}

} // namespace PlanarDetail

/** @brief All authorable parameter families of the planar selector. */
struct PlanarSampleParams : SampleCrossingParams {
  using Weights =
      Color::EnumWeights<Source::PLANAR_KIND_COUNT>; ///< Source weights.
  using Family = Source::PlanarFamily; ///< Planar parameter family ids.
  Source::GridSourceParams grid;       ///< Grid controls.
  Source::AuthoredTwinWaveSourceParams twin_wave; ///< Twin-wave controls.
  Source::BandsSourceParams bands;                ///< Parallel band controls.
  Source::SpiralSourceParams spiral;              ///< Spiral controls.
  Source::TessellationSourceParams tessellation;  ///< Tessellation controls.
  Source::NoiseSourceParams noise;                ///< Projected noise controls.
  uint8_t source = static_cast<uint8_t>(Source::PlanarKind::GRID);
  ///< Authored planar source option.
  PreviewEnumWeights<Source::PLANAR_KIND_COUNT> source_blend;
  ///< Transient planar source weights.

  /// Numeric crossing fields.
  static constexpr auto FIELDS = sample_crossing_fields<PlanarSampleParams>();
  /// Crossing selectors and planar source kind.
  static constexpr auto TOPOLOGY = sample_crossing_topology<PlanarSampleParams>(
      TopologyField<PlanarSampleParams>{"source", &PlanarSampleParams::source,
                                        Source::PLANAR_KIND_IDS, 0, true});
  /// Namespaced fields in parameter family order.
  static constexpr auto EXTRA_SCHEMA = PlanarDetail::join_schema(
      PlanarDetail::family_schema<Family::GRID, Source::GridSourceParams>(),
      PlanarDetail::family_schema<Family::TWIN_WAVE,
                                  Source::AuthoredTwinWaveSourceParams>(),
      PlanarDetail::family_schema<Family::BANDS, Source::BandsSourceParams>(),
      PlanarDetail::family_schema<Family::SPIRAL, Source::SpiralSourceParams>(),
      PlanarDetail::family_schema<Family::TESSELLATION,
                                  Source::TessellationSourceParams>(),
      PlanarDetail::family_schema<Family::PROJECTED_NOISE,
                                  Source::NoiseSourceParams>(),
      std::array{ParamFieldInfo{
          PlanarDetail::SPIRAL_COORDINATE_SYSTEM, nullptr, 0, 0, 0,
          FieldCurve::SNAP, true,
          std::size(Source::SpiralSourceParams::COORDINATE_SYSTEM_IDS),
          Source::SpiralSourceParams{}.coordinate_system,
          Source::SpiralSourceParams::COORDINATE_SYSTEM_IDS, nullptr, 0}});
  static_assert(FIELDS.size() + TOPOLOGY.size() + EXTRA_SCHEMA.size() == 28);

  /** @brief Finds a namespaced field's storage.
   * @param index Index within the extra schema.
   * @return The field's address. */
  void *extra_address(size_t index) {
    const auto address = [&]<typename... Params>(Params &...params) -> void * {
      void *found = nullptr;
      const auto visit = [&]<typename P>(P &family) {
        if (index < P::FIELDS.size()) {
          found = &(family.*(P::FIELDS[index].member));
          return true;
        }
        index -= P::FIELDS.size();
        return false;
      };
      if ((visit(params) || ...))
        return found;
      return &spiral.coordinate_system;
    };
    return address(grid, twin_wave, bands, spiral, tessellation, noise);
  }

  /** @brief Gets the family source mask for a field.
   * @param id Field id.
   * @return Live source kinds, or zero for common fields. */
  static constexpr auto field_sources(std::string_view id) {
    return Source::planar_field_sources(id);
  }

  /** @brief Overrides or clears preview source weights.
   * @param blend Option blend, or empty for the authored source.
   * @return True when an override was present. */
  bool set_enum_blend(const std::optional<EnumBlend> &blend) {
    return source_blend.set(blend);
  }

  /** @brief Resolves authored or transient source weights.
   * @return Evaluation weights of the planar kinds. */
  Weights weights() const { return source_blend.get(source); }

  /** @brief Whether any live kind reads this field at this endpoint.
   * @param id Family-qualified field id.
   * @param gate Field's topology gate.
   * @param live_kinds Mask of kinds with nonzero weight.
   * @return True when a live source uses the field. */
  bool field_live(std::string_view id, TopologyGate gate,
                  uint16_t live_kinds) const {
    const uint16_t family = field_sources(id);
    if (family != 0 && (family & live_kinds) == 0)
      return false;
    if (gate.field == nullptr)
      return true;
    if (std::string_view(gate.field) == "source")
      return (live_kinds & gate.values) != 0;
    const uint8_t selector = std::string_view(gate.field) == "coverage-mode"
                                 ? coverage_mode
                                 : spiral.coordinate_system;
    return (gate.values & (uint16_t{1} << selector)) != 0;
  }

  /** @brief Mask of kinds contributing to an endpoint.
   * @return Bit per nonzero source weight. */
  uint16_t live_kinds() const {
    if (!source_blend.active)
      return uint16_t{1} << source;
    uint16_t live = 0;
    for (size_t index = 0; index < source_blend.weights.values.size(); ++index)
      if (source_blend.weights.values[index] != 0.0f)
        live |= static_cast<uint16_t>(uint16_t{1} << index);
    return live;
  }

  /** @brief Samples source weights and retains fields at their live endpoint.
   * @param value Interpolated column parameters.
   * @param from Outgoing frame parameters.
   * @param to Incoming frame parameters.
   * @param fraction Column frame fraction. */
  static void sample_extra(PlanarSampleParams &value,
                           const PlanarSampleParams &from,
                           const PlanarSampleParams &to, float fraction) {
    value.source_blend.sample(from.weights(), to.weights(), fraction);
    const uint16_t from_kinds = from.live_kinds(), to_kinds = to.live_kinds();
    const auto field_endpoints = [&](std::string_view id, TopologyGate gate) {
      return std::pair{from.field_live(id, gate, from_kinds),
                       to.field_live(id, gate, to_kinds)};
    };
    if (const auto [from_live, to_live] =
            field_endpoints(FIELDS[0].id, FIELDS[0].topology_gate);
        from_live != to_live)
      value.edge_width = (from_live ? from : to).edge_width;

    const auto sample_family = [&]<Family F, typename Params>(Params &out,
                                                              const Params &a,
                                                              const Params &b) {
      constexpr auto SCHEMA = PlanarDetail::family_schema<F, Params>();
      out = Fields::interpolate(a, b, fraction);
      for (size_t index = 0; index < SCHEMA.size(); ++index) {
        const auto &info = SCHEMA[index];
        const auto [from_live, to_live] =
            field_endpoints(info.id, {info.gated_by, info.gate_values});
        if (from_live != to_live)
          out.*(Params::FIELDS[index].member) =
              (from_live ? a : b).*(Params::FIELDS[index].member);
      }
    };
    sample_family.template operator()<Family::GRID>(value.grid, from.grid,
                                                    to.grid);
    sample_family.template operator()<Family::TWIN_WAVE>(
        value.twin_wave, from.twin_wave, to.twin_wave);
    sample_family.template operator()<Family::BANDS>(value.bands, from.bands,
                                                     to.bands);
    sample_family.template operator()<Family::SPIRAL>(value.spiral, from.spiral,
                                                      to.spiral);
    sample_family.template operator()<Family::TESSELLATION>(
        value.tessellation, from.tessellation, to.tessellation);
    sample_family.template operator()<Family::PROJECTED_NOISE>(
        value.noise, from.noise, to.noise);

    const uint16_t spiral_mask = live_values(Source::PlanarKind::SPIRAL);
    if ((from_kinds & spiral_mask) == 0)
      value.spiral.coordinate_system = to.spiral.coordinate_system;
    else if ((to_kinds & spiral_mask) == 0)
      value.spiral.coordinate_system = from.spiral.coordinate_system;
  }
};

/** @brief Continuously advanced clocks of the planar source families. */
struct PlanarSourceState {
  SourceClockState grid;      ///< Grid phase clock.
  SourceClockState twin_wave; ///< Twin-wave phase and orientation clocks.
  NoisePhaseState noise;      ///< Projected noise field and phase clock.
};

/** @brief Column-time intervals of the planar source clocks. */
struct PlanarSourceInterval {
  SourceClockInterval grid;       ///< Grid clock interval.
  SourceClockInterval twin_wave;  ///< Twin-wave clock intervals.
  Spacetime::ClockInterval noise; ///< Projected-noise clock interval.
};

/** @brief Prepares planar source clocks between consecutive frames.
 * @param from Outgoing frame state.
 * @param to Incoming frame state.
 * @return Column-time clock intervals. */
inline PlanarSourceInterval temporal_interval(const PlanarSourceState &from,
                                              const PlanarSourceState &to) {
  return {temporal_interval(from.grid, to.grid),
          temporal_interval(from.twin_wave, to.twin_wave),
          temporal_interval(from.noise, to.noise)};
}

/** @brief Samples every planar source clock for a column.
 * @param out Column state.
 * @param interval Endpoint clock intervals.
 * @param params Sampled column controls.
 * @param fraction Column frame fraction. */
inline void temporal_sample(PlanarSourceState &out,
                            const PlanarSourceInterval &interval,
                            const PlanarSampleParams &params, float fraction) {
  temporal_sample(out.grid, interval.grid, params, fraction);
  temporal_sample(out.twin_wave, interval.twin_wave, params, fraction);
  temporal_sample(out.noise, interval.noise, params, fraction);
}

/** @brief Plane to field crossing over an authored blend of planar kinds. */
struct SamplePlanar : ValueStateModel<PlanarSourceState>,
                      SampleCoverageEdgeDistance {
  static constexpr const char *ID = "sample.planar.v1";
  static constexpr const char *NAME = "Planar Source";
  using Input = PlaneSample;
  using Output = FieldSample;
  using Params = PlanarSampleParams;
  using Kind = Source::PlanarKind;
  using Family = Source::PlanarFamily;

  /** @brief One kind's column preparation and field contribution. */
  struct Arm {
    Kind kind = Kind::GRID; ///< Source kind of this arm.
    float weight = 0.0f;    ///< Field weight in [0, 1].
    std::variant<std::monostate, Pattern::PreparedSource,
                 Pattern::PlaneBands::Prepared, Pattern::PlaneSpiral::Prepared,
                 PreparedNoiseSource>
        kernel; ///< Prepared kernel arguments.
  };

  /** @brief Prepared arms of a column's source blend. */
  struct Prepared {
    std::array<Arm, 2> arms; ///< Consecutive nonzero arms.
    uint8_t count = 0;       ///< Number of nonzero arms.
  };

  static void init(State &state, InstanceId id) {
    init_noise_phase(state.noise,
                     {id.instance, SampleProjectedNoise::ID,
                      instance_hash(id.instance, SampleProjectedNoise::ID)});
  }

  static void advance(State &state, const Params &params) {
    SourceClockModel::advance(state.grid, params.grid);
    SourceClockModel::advance(state.twin_wave, params.twin_wave);
    PhaseClockModel<NoisePhaseState>::advance(state.noise, params.noise);
  }

  /** @brief Lattice of the selected source option.
   * @param params Source controls.
   * @return Outgoing lattice until the incoming option is exact. */
  static TranslationLattice translation_lattice(const Params &params) {
    switch (const auto kind = static_cast<Kind>(params.source);
            Source::planar_family(kind)) {
    case Family::GRID:
      return Source::translation_lattice(params.grid);
    case Family::BANDS:
      return Source::translation_lattice(params.bands);
    case Family::TESSELLATION:
      return Source::tessellation_lattice(PlanarDetail::tessellation_kind(kind),
                                          params.tessellation.cell_scale);
    default:
      return {};
    }
  }

  /** @brief Rejects a column interval requiring more than two planar kinds.
   * @param from Outgoing frame parameters.
   * @param to Incoming frame parameters.
   * @return True when the union of live kinds fits the prepared arms. */
  static bool history_continuous(const Params &from, const Params &to) {
    return std::popcount(
               static_cast<uint16_t>(from.live_kinds() | to.live_kinds())) <= 2;
  }

  /** @brief Allows endpoint-only fields and preview kinds to vary by column.
   * @param field Schema field index.
   * @param from Outgoing frame parameters.
   * @param to Incoming frame parameters.
   * @return True when a raw field change preserves source continuity. */
  static bool temporal_hold_field(uint16_t field, const Params &from,
                                  const Params &to) {
    const auto &info = SCHEMA<SamplePlanar>[field];
    if (info.blend)
      return from.source_blend.active || to.source_blend.active;
    const TopologyGate gate{info.gated_by, info.gate_values};
    return !from.field_live(info.id, gate, from.live_kinds()) ||
           !to.field_live(info.id, gate, to.live_kinds());
  }

  static Prepared prepare(const FrameContext &, const Params &params,
                          const State &state) {
    check_sample_topology(params.weight_mode, params.coverage_mode);
    const auto weights = params.weights();
    Prepared result{};
    if (weights.exact != Params::Weights::MIXED) {
      result.arms[0] =
          prepare_arm(static_cast<Kind>(weights.exact), 1.0f, params, state);
      result.count = 1;
      return result;
    }
    for (size_t index = 0; index < weights.values.size(); ++index)
      if (weights.values[index] != 0.0f)
        result.arms[result.count++] = prepare_arm(
            static_cast<Kind>(index), weights.values[index], params, state);
    return result;
  }

  static FieldSample run(const PlaneSample &input, const FrameContext &ctx,
                         const Params &params, const Prepared &prepared) {
    const auto &a = prepared.arms[0];
    float raw = sample(input, params, a);
    if (prepared.count == 2) {
      const auto &b = prepared.arms[1];
      raw = a.weight * raw + b.weight * sample(input, params, b);
    }
    return finish_sample(input, raw, params, ctx);
  }

private:
  static Arm prepare_arm(Kind kind, float weight, const Params &params,
                         const State &state) {
    Arm arm{kind, weight, {}};
    switch (Source::planar_family(kind)) {
    case Family::GRID:
      arm.kernel = Pattern::prepare_phases(
          state.grid.primary, state.grid.secondary, state.grid.angle);
      break;
    case Family::TWIN_WAVE:
      arm.kernel = Pattern::prepare_phases(state.twin_wave.primary,
                                           state.twin_wave.secondary,
                                           state.twin_wave.angle);
      break;
    case Family::BANDS:
      arm.kernel = Pattern::PlaneBands::prepare(params.bands, {});
      break;
    case Family::SPIRAL:
      arm.kernel = Pattern::PlaneSpiral::prepare(params.spiral, {});
      break;
    case Family::TESSELLATION:
      break;
    case Family::PROJECTED_NOISE:
      arm.kernel = PreparedNoiseSource{
          &state.noise.noise,
          math::noise_projected_loop_offset(state.noise.phase)};
      break;
    }
    return arm;
  }

  static float sample(const PlaneSample &input, const Params &params,
                      const Arm &arm) {
    switch (Source::planar_family(arm.kind)) {
    case Family::GRID: {
      const auto &grid = params.grid;
      return Pattern::grid(
          Pattern::PlanePosition{projections::stereo_pattern_args(
              input.coords, grid.pattern_freq)},
          grid, std::get<Pattern::PreparedSource>(arm.kernel));
    }
    case Family::TWIN_WAVE: {
      const auto &twin = params.twin_wave;
      return Pattern::twin_wave(
          Pattern::PlanePosition{projections::stereo_pattern_args(
              input.coords, twin.pattern_freq)},
          std::get<Pattern::PreparedSource>(arm.kernel));
    }
    case Family::BANDS:
      return Pattern::PlaneBands::sample(
          {input.coords}, std::get<Pattern::PlaneBands::Prepared>(arm.kernel));
    case Family::SPIRAL:
      return Pattern::PlaneSpiral::sample(
          {input.coords}, std::get<Pattern::PlaneSpiral::Prepared>(arm.kernel));
    case Family::TESSELLATION:
      return Pattern::tessellation(Pattern::PlanePosition{input.coords},
                                   params.tessellation,
                                   PlanarDetail::tessellation_kind(arm.kind),
                                   PlanarDetail::tessellation_fill(arm.kind));
    case Family::PROJECTED_NOISE:
      if (math::noise_plane_in_domain(input.coords)) {
        const auto &noise = std::get<PreparedNoiseSource>(arm.kernel);
        return Pattern::noise_contour(
            *noise.noise, PlanarDetail::noise_basis(arm.kind),
            math::noise_projected_coordinate(
                input.coords, params.noise.noise_scale, noise.loop_offset),
            params.noise.noise_contrast);
      }
      break;
    }
    return 0.0f;
  }
};

} // namespace Pullback::Interp::Op

#endif // HS_ENABLE_CHAIN_INTERPRETER
