/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/modulation.h"
#include "core/control/fields.h"
#include "core/pattern/noise.h"
#include "core/pattern/sources.h"
#include "core/pattern/sphere_lut.h"
#include "core/render/pullback/runtime_seeds.h"
#include <array>
#include <optional>
#include <variant>

/** @file modulation.h
 * @brief Sphere modulation validation, scheduling and owned resources. */

namespace Pullback::Modulation {
using Config = ColorModulation::Config; ///< Canonical passive configuration.
/** @brief The selected source's versioned clock. */
struct State {
  Pattern::SourceKind source = Pattern::SourceKind::NOISE; ///< Clock recipe.
  uint16_t version = 1;                                    ///< Recipe version.
  Pattern::PhaseState phase; ///< Selected source phase in turns.
  bool initialized = false;  ///< Whether the source was initialized.
  constexpr bool operator==(const State &) const = default;
};
/** @brief Tests a source recipe, including disabled configurations.
 * @param config Canonical recipe.
 * @return Whether parameters and the selected strategy are admitted. */
inline bool valid(const Config &config) {
  if (config.version != 1 || !ColorModulation::valid(config.mechanic) ||
      static_cast<unsigned>(config.evaluation) > 1)
    return false;
  if (config.evaluation == Pattern::Evaluation::BAKED_CUBE &&
      (config.source != Pattern::SourceKind::NOISE ||
       config.noise.scale > 8.0f))
    return false;
  if (config.evaluation == Pattern::Evaluation::DIRECT &&
      config.source == Pattern::SourceKind::NOISE &&
      (config.noise.variant != Pattern::NoiseVariant::CONFIGURED_SINGLE ||
       config.mechanic.kind ==
           ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT))
    return false;
  switch (config.source) {
  case Pattern::SourceKind::NOISE:
    return Pattern::SphereNoise::valid(config.noise);
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    return Pattern::DirectionalSine::valid(config.sine);
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    return Pattern::SphereInterference::valid(config.interference);
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    return Pattern::SphereSmoothGrid::valid(config.grid);
  }
  return false;
}
/** @brief Initializes the selected source clock.
 * @param config Source recipe.
 * @return Initialized typed source state. */
inline State init(const Config &config) {
  HS_CHECK(valid(config), "invalid modulation initialization recipe");
  Pattern::PhaseState phase;
  switch (config.source) {
  case Pattern::SourceKind::NOISE:
    phase = Pattern::SphereNoise::init_state(config.noise);
    break;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    phase = Pattern::DirectionalSine::init_state(config.sine);
    break;
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    phase = Pattern::SphereInterference::init_state(config.interference);
    break;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    phase = Pattern::SphereSmoothGrid::init_state(config.grid);
    break;
  }
  return {config.source, config.version, phase, true};
}
/** @brief Adopts source identity without reinterpreting clocks.
 * @param state Scheduled source state.
 * @param config Current configuration. */
inline void synchronize(State &state, const Config &config) {
  HS_CHECK(valid(config), "invalid modulation synchronization recipe");
  if (!state.initialized || state.source != config.source ||
      state.version != config.version)
    state = init(config);
}
/** @brief Advances the selected clock once, including disabled recipes.
 * @param state Scheduled source state.
 * @param config Current controls. */
inline void advance(State &state, const Config &config) {
  synchronize(state, config);
  switch (config.source) {
  case Pattern::SourceKind::NOISE:
    state.phase = Pattern::SphereNoise::advance(state.phase, config.noise);
    break;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    state.phase = Pattern::DirectionalSine::advance(state.phase, config.sine);
    break;
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    state.phase =
        Pattern::SphereInterference::advance(state.phase, config.interference);
    break;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    state.phase = Pattern::SphereSmoothGrid::advance(state.phase, config.grid);
    break;
  }
}
/** @brief Topology equality required for endpoint color interpolation.
 * @param a Previous recipe.
 * @param b Current recipe.
 * @return Whether both endpoints have identical topology. */
inline bool same_topology(const Config &a, const Config &b) {
  return a.enabled == b.enabled && a.source == b.source &&
         a.version == b.version && a.evaluation == b.evaluation &&
         a.mechanic.kind == b.mechanic.kind;
}
/** @brief Interpolates controls within one recipe topology.
 * @param a Previous controls.
 * @param b Current controls.
 * @param fraction Current weight.
 * @return Current topology, with continuous selected controls interpolated. */
inline Config interpolate(const Config &a, const Config &b, float fraction) {
  if (!same_topology(a, b))
    return b;
  if (fraction <= 0.0f)
    return a;
  if (fraction >= 1.0f)
    return b;
  Config result = b;
  result.mechanic.amount =
      hs::lerp(a.mechanic.amount, b.mechanic.amount, fraction);
  const auto wave = [fraction](const auto &start, const auto &end, auto &out) {
    out.frequency = hs::lerp(start.frequency, end.frequency, fraction);
    out.offset = hs::lerp(start.offset, end.offset, fraction);
    out.speed = hs::lerp(start.speed, end.speed, fraction);
  };
  switch (b.source) {
  case Pattern::SourceKind::NOISE:
    if (a.noise.basis == b.noise.basis && a.noise.variant == b.noise.variant) {
      result.noise.scale =
          Control::apply_curve(Control::FieldCurve::LOG_POSITIVE, a.noise.scale,
                               b.noise.scale, fraction);
      result.noise.speed = hs::lerp(a.noise.speed, b.noise.speed, fraction);
      result.noise.contrast =
          hs::lerp(a.noise.contrast, b.noise.contrast, fraction);
    }
    break;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    wave(a.sine, b.sine, result.sine);
    break;
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    wave(a.interference, b.interference, result.interference);
    result.interference.frequency_b = hs::lerp(
        a.interference.frequency_b, b.interference.frequency_b, fraction);
    result.interference.offset_b =
        hs::lerp(a.interference.offset_b, b.interference.offset_b, fraction);
    break;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    wave(a.grid, b.grid, result.grid);
    result.grid.frequency_b =
        hs::lerp(a.grid.frequency_b, b.grid.frequency_b, fraction);
    result.grid.offset_b = hs::lerp(a.grid.offset_b, b.grid.offset_b, fraction);
    break;
  }
  return result;
}
/** @brief Immutable prepared source variant. */
using Source = std::variant<
    Pattern::SphereNoise::Prepared, Pattern::DirectionalSine::Prepared,
    Pattern::SphereInterference::Prepared, Pattern::SphereSmoothGrid::Prepared>;
/** @brief One endpoint's immutable source inputs and borrowed scalar table. */
struct Prepared {
  Config config; ///< Endpoint recipe and controls.
  State state;   ///< Endpoint source clock.
  Source source; ///< Prepared source, including an owned noise generator.
  Pattern::SphereLutView cube{nullptr, false}; ///< Optional scalar table.
  float (*evaluate)(const Prepared &, const Pattern::SphereDirection &) =
      nullptr; ///< Resolved built-in sampler.
};
/** @brief Samples one admitted direction without reading mutable controls.
 * @param prepared Endpoint source.
 * @param direction Admitted unit coordinate.
 * @return Signed scalar in [-1, 1]. */
inline float sample(const Prepared &prepared,
                    const Pattern::SphereDirection &direction) {
  HS_CHECK(prepared.evaluate, "modulation endpoint is unprepared");
  return Pattern::bounded_scalar(prepared.evaluate(prepared, direction));
}
/** @brief Complete spatial cache identity, independent of mechanic and palette. */
struct FieldIdentity {
  Pattern::SourceKind source;     ///< Source kind.
  uint16_t version;               ///< Source version.
  Pattern::Evaluation evaluation; ///< Evaluation layout selector.
  Pattern::NoiseParams noise;     ///< Noise parameters and generator selection.
  Pattern::SineParams sine;       ///< Directional sine parameters.
  Pattern::InterferenceParams interference; ///< Interference parameters.
  Pattern::GridParams grid;                 ///< Grid parameters.
  int32_t seed;                             ///< Generator seed.
  Pattern::PhaseState phase;                ///< Source state.
  uint16_t bake_layout_version =
      Pattern::SphereLutView::LAYOUT_VERSION; ///< Cube storage version.
  constexpr bool operator==(const FieldIdentity &) const = default;
};
/** @brief Builds the selected source identity without dormant parameters.
 * @param config Endpoint configuration.
 * @param state Endpoint source state.
 * @return Complete source and bake identity. */
inline FieldIdentity field_identity(const Config &config, const State &state) {
  FieldIdentity identity{
      config.source, config.version, config.evaluation, {}, {}, {}, {}, 0,
      state.phase};
  switch (config.source) {
  case Pattern::SourceKind::NOISE:
    identity.noise = config.noise;
    identity.seed = config.seed;
    break;
  case Pattern::SourceKind::DIRECTIONAL_SINE:
    identity.sine = config.sine;
    break;
  case Pattern::SourceKind::SPHERE_INTERFERENCE:
    identity.interference = config.interference;
    break;
  case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
    identity.grid = config.grid;
    break;
  }
  return identity;
}
/** @brief Owned prepared endpoint and optional scalar cube. */
template <bool Baked = true> struct Resources {
  Resources() = default;
  Resources(const Resources &) = delete;
  Resources &operator=(const Resources &) = delete;
  Prepared prepared; ///< Current prepared source.
  std::array<int8_t, Baked ? Pattern::SphereLutView::SIZE : 0>
      table{};                      ///< Scalar cube.
  std::optional<FieldIdentity> key; ///< Resident cube identity.
  uint32_t bake_generation = 0;     ///< Number of completed cube bakes.
  /** @brief Invalidates derived resources after restore or ownership changes. */
  void reset() { key.reset(); }
  /** @brief Prepares immutable inputs and refreshes a changed cube.
   * @param config Current recipe.
   * @param clock Current source state.
   * @return Borrowed prepared endpoint. */
  const Prepared &prepare(const Config &config, const State &clock) {
    HS_CHECK(valid(config), "invalid modulation preparation recipe");
    State selected = clock;
    synchronize(selected, config);
    prepared.config = config;
    prepared.state = selected;
    prepared.cube = {table.data(),
                     config.enabled &&
                         config.evaluation == Pattern::Evaluation::BAKED_CUBE};
    switch (config.source) {
    case Pattern::SourceKind::NOISE: {
      FastNoiseLite noise;
      init_effect_noise(noise, config.seed);
      prepared.source.emplace<0>(
          Pattern::SphereNoise::prepare(config.noise, selected.phase, noise));
      break;
    }
    case Pattern::SourceKind::DIRECTIONAL_SINE:
      prepared.source.emplace<1>(
          Pattern::DirectionalSine::prepare(config.sine, selected.phase));
      break;
    case Pattern::SourceKind::SPHERE_INTERFERENCE:
      prepared.source.emplace<2>(Pattern::SphereInterference::prepare(
          config.interference, selected.phase));
      break;
    case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
      prepared.source.emplace<3>(
          Pattern::SphereSmoothGrid::prepare(config.grid, selected.phase));
      break;
    }
    if (config.evaluation == Pattern::Evaluation::BAKED_CUBE) {
      prepared.evaluate = [](const Prepared &value,
                             const Pattern::SphereDirection &coordinate) {
        return Pattern::sample_sphere_lut(value.cube, coordinate);
      };
    } else {
      switch (config.source) {
      case Pattern::SourceKind::NOISE:
        prepared.evaluate = [](const Prepared &value,
                               const Pattern::SphereDirection &coordinate) {
          return Pattern::SphereNoise::sample(coordinate,
                                              std::get<0>(value.source));
        };
        break;
      case Pattern::SourceKind::DIRECTIONAL_SINE:
        prepared.evaluate = [](const Prepared &value,
                               const Pattern::SphereDirection &coordinate) {
          return Pattern::DirectionalSine::sample(coordinate,
                                                  std::get<1>(value.source));
        };
        break;
      case Pattern::SourceKind::SPHERE_INTERFERENCE:
        prepared.evaluate = [](const Prepared &value,
                               const Pattern::SphereDirection &coordinate) {
          return Pattern::SphereInterference::sample(coordinate,
                                                     std::get<2>(value.source));
        };
        break;
      case Pattern::SourceKind::SPHERE_SMOOTH_GRID:
        prepared.evaluate = [](const Prepared &value,
                               const Pattern::SphereDirection &coordinate) {
          return Pattern::SphereSmoothGrid::sample(coordinate,
                                                   std::get<3>(value.source));
        };
        break;
      }
    }
    if (prepared.cube.active &&
        config.evaluation == Pattern::Evaluation::BAKED_CUBE) {
      HS_CHECK(Baked, "modulation owner has no baked cube capacity");
      if constexpr (Baked) {
        const FieldIdentity IDENTITY = field_identity(config, selected);
        if (!key || *key != IDENTITY) {
          Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
              std::span<int8_t, Pattern::SphereLutView::SIZE>(table),
              std::get<0>(prepared.source));
          key = IDENTITY;
          ++bake_generation;
        }
      }
    }
    return prepared;
  }
};
} // namespace Pullback::Modulation
