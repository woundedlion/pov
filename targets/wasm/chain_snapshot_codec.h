/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file chain_snapshot_codec.h
 * @brief Strict JS wire codec for `ChainSnapshot` runtime/palette snapshots.
 */
#pragma once

#include <emscripten/val.h>
#include <algorithm>
#include <limits>
#include "core/render/pullback/interpreter.h"
#include "workbench/shader/chain_snapshot.h"

namespace hs_wasm {

/**
 * @brief Encodes a `ChainSnapshot` to a plain JS object and strictly decodes
 *        one back.
 */
class ChainSnapshotCodec {
  using Val = emscripten::val;
  using Result = ChainSnapshotRestoreResult;
  using Runtime = Pullback::Interp::RuntimeSnapshot;

  static bool array_type(const Val &value) {
    return Val::global("Array").call<bool>("isArray", value);
  }
  static bool array(const Val &value, size_t max_count) {
    return array_type(value) && value["length"].as<size_t>() <= max_count;
  }
  static bool object(const Val &value) {
    return !value.isNull() && !value.isUndefined() &&
           value.typeOf().as<std::string>() == "object";
  }
  static bool number(const Val &object, const char *key, float &out) {
    const auto value = object[key];
    if (!value.isNumber())
      return false;
    out = value.as<float>();
    return std::isfinite(out);
  }
  template <typename T>
  static bool integer(const Val &object, const char *key, T &out) {
    const auto value = object[key];
    if (!value.isNumber())
      return false;
    const double number = value.as<double>();
    if (!std::isfinite(number) || number != std::floor(number) ||
        number < static_cast<double>(std::numeric_limits<T>::min()) ||
        number > static_cast<double>(std::numeric_limits<T>::max()))
      return false;
    out = static_cast<T>(number);
    return true;
  }
  static bool boolean(const Val &object, const char *key, bool &out) {
    const auto value = object[key];
    if (value.typeOf().as<std::string>() != "boolean")
      return false;
    out = value.as<bool>();
    return true;
  }
  static bool text(const Val &object, const char *key, std::string &out) {
    const auto value = object[key];
    if (!value.isString())
      return false;
    out = value.as<std::string>();
    return out.find('\0') == std::string::npos;
  }
  static Val vector(const math::Vector &value) {
    auto out = Val::array();
    out.set(0, value.x);
    out.set(1, value.y);
    out.set(2, value.z);
    return out;
  }
  static Val quaternion(const math::Quaternion &value) {
    auto out = Val::array();
    out.set(0, value.r);
    out.set(1, value.v.x);
    out.set(2, value.v.y);
    out.set(3, value.v.z);
    return out;
  }
  template <size_t N>
  static bool floats(const Val &value, std::array<float, N> &out) {
    if (!array(value, N) || value["length"].as<size_t>() != N)
      return false;
    for (size_t index = 0; index < N; ++index) {
      if (!value[index].isNumber())
        return false;
      out[index] = value[index].as<float>();
      if (!std::isfinite(out[index]))
        return false;
    }
    return true;
  }
  static bool vector(const Val &value, math::Vector &out) {
    std::array<float, 3> values;
    if (!floats(value, values))
      return false;
    out = math::Vector(values[0], values[1], values[2]);
    return true;
  }
  static bool quaternion(const Val &value, math::Quaternion &out) {
    std::array<float, 4> values;
    if (!floats(value, values))
      return false;
    out = math::Quaternion(values[0], values[1], values[2], values[3]);
    return true;
  }
  static Val spatial(const Pullback::Interp::SpatialWalkSnapshot &state) {
    auto out = Val::object();
    out.set("noiseSeed", state.noise_seed);
    out.set("walkTime", state.walk_time);
    out.set("position", vector(state.position));
    out.set("direction", vector(state.direction));
    out.set("wander", quaternion(state.wander));
    out.set("angularVelocity", state.angular_velocity);
    out.set("spinPhase", state.spin_phase);
    return out;
  }
  static bool spatial(const Val &value,
                      Pullback::Interp::SpatialWalkSnapshot &state) {
    return object(value) && integer(value, "noiseSeed", state.noise_seed) &&
           integer(value, "walkTime", state.walk_time) &&
           vector(value["position"], state.position) &&
           vector(value["direction"], state.direction) &&
           quaternion(value["wander"], state.wander) &&
           number(value, "angularVelocity", state.angular_velocity) &&
           number(value, "spinPhase", state.spin_phase);
  }
  static Val sine(const Pattern::SineParams &params) {
    auto out = Val::object();
    out.set("axis", vector(params.axis));
    out.set("frequency", params.frequency);
    out.set("offset", params.offset);
    out.set("speed", params.speed);
    return out;
  }
  static bool sine(const Val &value, Pattern::SineParams &params) {
    return object(value) && vector(value["axis"], params.axis) &&
           number(value, "frequency", params.frequency) &&
           number(value, "offset", params.offset) &&
           number(value, "speed", params.speed);
  }
  static Val interference(const Pattern::InterferenceParams &params) {
    auto out = sine(params);
    out.set("axisB", vector(params.axis_b));
    out.set("frequencyB", params.frequency_b);
    out.set("offsetB", params.offset_b);
    return out;
  }
  static bool interference(const Val &value,
                           Pattern::InterferenceParams &params) {
    return sine(value, params) && vector(value["axisB"], params.axis_b) &&
           number(value, "frequencyB", params.frequency_b) &&
           number(value, "offsetB", params.offset_b);
  }
  static Val recipe(const ColorModulation::Config &config) {
    auto out = Val::object();
    out.set("enabled", config.enabled);
    out.set("source", static_cast<uint8_t>(config.source));
    out.set("version", config.version);
    out.set("evaluation", static_cast<uint8_t>(config.evaluation));
    out.set("seed", config.seed);
    auto noise = Val::object();
    noise.set("scale", config.noise.scale);
    noise.set("speed", config.noise.speed);
    noise.set("contrast", config.noise.contrast);
    noise.set("basis", static_cast<uint8_t>(config.noise.basis));
    noise.set("variant", static_cast<uint8_t>(config.noise.variant));
    out.set("noise", noise);
    out.set("sine", sine(config.sine));
    out.set("interference", interference(config.interference));
    out.set("grid", interference(config.grid));
    auto mechanic = Val::object();
    mechanic.set("kind", static_cast<uint8_t>(config.mechanic.kind));
    mechanic.set("amount", config.mechanic.amount);
    out.set("mechanic", mechanic);
    return out;
  }
  static bool recipe(const Val &value, ColorModulation::Config &config) {
    uint8_t source, evaluation, basis, variant, mechanic;
    if (!object(value) || !boolean(value, "enabled", config.enabled) ||
        !integer(value, "source", source) ||
        !integer(value, "version", config.version) || config.version != 1 ||
        !integer(value, "evaluation", evaluation) ||
        !integer(value, "seed", config.seed))
      return false;
    const auto noise = value["noise"];
    const auto operation = value["mechanic"];
    if (!object(noise) || !number(noise, "scale", config.noise.scale) ||
        !number(noise, "speed", config.noise.speed) ||
        !number(noise, "contrast", config.noise.contrast) ||
        !integer(noise, "basis", basis) ||
        !integer(noise, "variant", variant) ||
        !sine(value["sine"], config.sine) ||
        !interference(value["interference"], config.interference) ||
        !interference(value["grid"], config.grid) || !object(operation) ||
        !integer(operation, "kind", mechanic) ||
        !number(operation, "amount", config.mechanic.amount))
      return false;
    config.source = static_cast<Pattern::SourceKind>(source);
    config.evaluation = static_cast<Pattern::Evaluation>(evaluation);
    config.noise.basis = static_cast<Pattern::NoiseBasis>(basis);
    config.noise.variant = static_cast<Pattern::NoiseVariant>(variant);
    config.mechanic.kind = static_cast<ColorModulation::Mechanic>(mechanic);
    return Pullback::Modulation::valid(config);
  }
  static Val modulation_clock(const Pullback::Modulation::State &state) {
    auto out = Val::object();
    out.set("source", static_cast<uint8_t>(state.source));
    out.set("version", state.version);
    out.set("phaseTurns", state.phase.phase_turns);
    out.set("initialized", state.initialized);
    return out;
  }

  static bool modulation_clock(const Val &value,
                               Pullback::Modulation::State &state) {
    uint8_t source;
    if (!object(value) || !integer(value, "source", source) || source > 4 ||
        !integer(value, "version", state.version) || state.version != 1 ||
        !number(value, "phaseTurns", state.phase.phase_turns) ||
        state.phase.phase_turns < 0.0f || state.phase.phase_turns >= 1.0f ||
        !boolean(value, "initialized", state.initialized) ||
        (!state.initialized && state.phase.phase_turns != 0.0f))
      return false;
    state.source = static_cast<Pattern::SourceKind>(source);
    return true;
  }

  struct LegacyClock {
    std::string instance;
    Pullback::Modulation::State state;
  };

  static Result migrate_legacy(ChainSnapshot &snapshot,
                               const std::vector<LegacyClock> &clocks) {
    std::vector<std::string> names;
    for (const auto &entry : snapshot.chain)
      names.push_back(entry.instance);
    for (size_t index = 0; index < snapshot.parameters.size(); ++index) {
      const auto &name = snapshot.parameters[index].name;
      bool known = false;
      for (const auto &entry : snapshot.chain)
        known |= name.starts_with(entry.instance + ".");
      if (!known)
        return Result::INVALID_VALUE;
      for (size_t previous = 0; previous < index; ++previous)
        if (snapshot.parameters[previous].name == name)
          return Result::INVALID_VALUE;
    }
    if (snapshot.runtime)
      for (const auto &entry : *snapshot.runtime)
        if (std::find(names.begin(), names.end(), entry.instance) ==
            names.end())
          return Result::INVALID_VALUE;
    const auto allocate_name = [&](const char *base) {
      std::string candidate = base;
      for (size_t suffix = 2;
           std::find(names.begin(), names.end(), candidate) != names.end();
           ++suffix)
        candidate = std::string(base) + "-" + std::to_string(suffix);
      names.push_back(candidate);
      return candidate;
    };
    ChainSnapshot migrated;
    migrated.animations_paused = snapshot.animations_paused;
    migrated.palette_bank = snapshot.palette_bank;
    if (snapshot.runtime)
      migrated.runtime.emplace();
    size_t recipes = 0;
    for (const auto &entry : snapshot.chain) {
      const bool field_layer = entry.operator_id == "modulate.field.v1";
      const bool path_layer = entry.operator_id == "modulate.hue-shift.v1";
      if (field_layer || path_layer) {
        ColorModulation::Config config;
        Pullback::Modulation::State state;
        if (field_layer) {
          const ChainSnapshot::ModulationRecipe *selected = nullptr;
          for (const auto &recipe : snapshot.modulation)
            if (recipe.instance == entry.instance) {
              if (selected != nullptr)
                return Result::INVALID_VALUE;
              selected = &recipe;
            }
          if (selected == nullptr)
            return Result::INVALID_VALUE;
          config = selected->recipe;
          ++recipes;
          for (const auto &value : snapshot.parameters)
            if (value.name.starts_with(entry.instance + "."))
              return Result::INVALID_VALUE;
          if (snapshot.runtime) {
            const LegacyClock *selected_clock = nullptr;
            for (const auto &clock : clocks)
              if (clock.instance == entry.instance) {
                if (selected_clock != nullptr)
                  return Result::INVALID_VALUE;
                selected_clock = &clock;
              }
            size_t runtime_count = 0;
            for (const auto &runtime : *snapshot.runtime)
              runtime_count += runtime.instance == entry.instance ? 1 : 0;
            if (selected_clock == nullptr || runtime_count != 1)
              return Result::INVALID_VALUE;
            state = selected_clock->state;
          }
        } else {
          config.seed = Pullback::HUE_NOISE_SEED;
          config.source = Pattern::SourceKind::DISPLACEMENT;
          config.evaluation = Pattern::Evaluation::DIRECT;
          config.enabled = true;
          config.mechanic.amount = 0.0f;
          for (const auto &value : snapshot.parameters)
            if (value.name.starts_with(entry.instance + ".")) {
              const auto field = value.name.substr(entry.instance.size() + 1);
              if (field == "hue-shift-mode") {
                if (value.value != 0.0f && value.value != 1.0f)
                  return Result::INVALID_VALUE;
                config.enabled = value.value == 1.0f;
              } else if (field == "hue-shift-amount") {
                if (value.value < -4.0f || value.value > 4.0f)
                  return Result::INVALID_VALUE;
                config.mechanic.amount = value.value;
              } else
                return Result::INVALID_VALUE;
            }
          if (snapshot.runtime)
            for (const auto &runtime : *snapshot.runtime)
              if (runtime.instance == entry.instance)
                return Result::INVALID_VALUE;
          state.source = Pattern::SourceKind::DISPLACEMENT;
          state.initialized = true;
        }
        migrated.chain.push_back({entry.instance, "color.modulate.v1"});
        migrated.modulation.push_back({entry.instance, config});
        if (snapshot.runtime)
          migrated.runtime->push_back(
              {entry.instance,
               Pullback::Interp::ModulationClockSnapshot{state}});
        continue;
      }
      const bool color = entry.operator_id == "colorize.generated-palette.v3";
      migrated.chain.push_back(
          {entry.instance,
           color ? "colorize.generated-palette.v4" : entry.operator_id});
      if (!color) {
        for (const auto &value : snapshot.parameters)
          if (value.name.starts_with(entry.instance + "."))
            migrated.parameters.push_back(value);
        if (snapshot.runtime)
          for (const auto &runtime : *snapshot.runtime)
            if (runtime.instance == entry.instance)
              migrated.runtime->push_back(runtime);
        continue;
      }
      const std::string displacement = allocate_name("displacement");
      const std::string spatial = allocate_name("modulation");
      migrated.chain.push_back({displacement, "color.modulate.v1"});
      migrated.chain.push_back({spatial, "color.modulate.v1"});
      ColorModulation::Config path;
      path.seed = Pullback::HUE_NOISE_SEED;
      path.source = Pattern::SourceKind::DISPLACEMENT;
      path.evaluation = Pattern::Evaluation::DIRECT;
      path.enabled = false;
      for (const auto &value : snapshot.parameters)
        if (value.name.starts_with(entry.instance + ".")) {
          const std::string field =
              value.name.substr(entry.instance.size() + 1);
          if (field == "hue-shift-mode") {
            if (value.value != 0.0f && value.value != 1.0f)
              return Result::INVALID_VALUE;
            path.enabled = value.value == 1.0f;
          } else if (field == "hue-shift-amount") {
            if (value.value < -4.0f || value.value > 4.0f)
              return Result::INVALID_VALUE;
            path.mechanic.amount = value.value;
          } else {
            migrated.parameters.push_back(value);
          }
        }
      const ChainSnapshot::ModulationRecipe *recipe = nullptr;
      for (const auto &candidate : snapshot.modulation)
        if (candidate.instance == entry.instance) {
          if (recipe != nullptr)
            return Result::INVALID_VALUE;
          recipe = &candidate;
        }
      if (recipe == nullptr)
        return Result::INVALID_VALUE;
      ++recipes;
      migrated.modulation.push_back({displacement, path});
      migrated.modulation.push_back({spatial, recipe->recipe});
      if (snapshot.runtime) {
        const LegacyClock *clock = nullptr;
        for (const auto &candidate : clocks)
          if (candidate.instance == entry.instance) {
            if (clock != nullptr)
              return Result::INVALID_VALUE;
            clock = &candidate;
          }
        if (clock == nullptr)
          return Result::INVALID_VALUE;
        for (const auto &runtime : *snapshot.runtime)
          if (runtime.instance == entry.instance)
            migrated.runtime->push_back(runtime);
        Pullback::Modulation::State path_clock;
        path_clock.source = Pattern::SourceKind::DISPLACEMENT;
        path_clock.initialized = true;
        migrated.runtime->push_back(
            {displacement,
             Pullback::Interp::ModulationClockSnapshot{path_clock}});
        migrated.runtime->push_back(
            {spatial, Pullback::Interp::ModulationClockSnapshot{clock->state}});
      }
    }
    if (recipes != snapshot.modulation.size())
      return Result::INVALID_VALUE;
    if (migrated.chain.size() > Pullback::Interp::MAX_CHAIN_OPS ||
        migrated.modulation.size() > Pullback::Interp::MAX_CHAIN_OPS)
      return Result::INVALID_LENGTH;
    snapshot = std::move(migrated);
    return Result::APPLIED;
  }

  static Val runtime_state(const Runtime &snapshot) {
    using namespace Pullback::Interp;
    return std::visit(
        [](const auto &state) {
          using State = std::decay_t<decltype(state)>;
          auto out = Val::object();
          if constexpr (std::is_same_v<State, SpatialWalkSnapshot>)
            out = spatial(state);
          else if constexpr (std::is_same_v<State, SourceClockSnapshot>) {
            out.set("primary", state.primary);
            out.set("secondary", state.secondary);
            out.set("angle", state.angle);
          } else if constexpr (std::is_same_v<State, NoiseClockSnapshot>) {
            out.set("phase", state.phase);
            out.set("noiseSeed", state.noise_seed);
          } else if constexpr (std::is_same_v<State, AffineClockSnapshot>) {
            out.set("phase", state.phase);
            out.set("rotation", state.rotation);
          } else if constexpr (std::is_same_v<State, RippleClockSnapshot>) {
            out.set("phase", state.phase);
            out.set("seed", state.seed);
            out.set("burstCount", state.burst_count);
            out.set("frames", state.frames);
            auto waves = Val::array();
            for (const auto &wave : state.waves) {
              auto entry = Val::object();
              entry.set("age", wave.age);
              entry.set("burst", wave.burst);
              waves.call<void>("push", entry);
            }
            out.set("waves", waves);
          } else if constexpr (std::is_same_v<State, ColorClockSnapshot>) {
            out.set("oscillationPhase", state.oscillation_phase);
          } else if constexpr (std::is_same_v<State, ModulationClockSnapshot>) {
            out = modulation_clock(state.modulation);
          } else if constexpr (std::is_same_v<State, SphericalRingsSnapshot>) {
            out.set("walk", spatial(state.walk));
            out.set("phase", state.phase);
          } else if constexpr (!std::is_same_v<State, std::monostate>)
            out.set("phase", state.phase);
          return out;
        },
        snapshot);
  }
  static bool
  runtime_state(const std::string &kind, const Val &value, Runtime &out,
                Pullback::Modulation::State *legacy_modulation = nullptr) {
    using namespace Pullback::Interp;
    if (!object(value))
      return false;
    if (kind == "spatial-walk-v2") {
      SpatialWalkSnapshot state;
      if (!spatial(value, state))
        return false;
      out = state;
    } else if (kind == "source-clock-v1") {
      SourceClockSnapshot state;
      if (!number(value, "primary", state.primary) ||
          !number(value, "secondary", state.secondary) ||
          !number(value, "angle", state.angle))
        return false;
      out = state;
    } else if (kind == "noise-clock-v1") {
      NoiseClockSnapshot state;
      if (!number(value, "phase", state.phase) ||
          !integer(value, "noiseSeed", state.noise_seed))
        return false;
      out = state;
    } else if (kind == "ripple-burst-v1") {
      RippleClockSnapshot state;
      const auto waves = value["waves"];
      if (!number(value, "phase", state.phase) ||
          !integer(value, "seed", state.seed) ||
          !integer(value, "burstCount", state.burst_count) ||
          !integer(value, "frames", state.frames) ||
          !Val::global("Array").call<bool>("isArray", waves) ||
          waves["length"].as<unsigned>() != state.waves.size())
        return false;
      for (unsigned i = 0; i < state.waves.size(); ++i)
        if (!object(waves[i]) || !number(waves[i], "age", state.waves[i].age) ||
            !integer(waves[i], "burst", state.waves[i].burst))
          return false;
      out = state;
    } else if (kind == "affine-clock-v1") {
      AffineClockSnapshot state;
      if (!number(value, "phase", state.phase) ||
          !number(value, "rotation", state.rotation))
        return false;
      out = state;
    } else if (kind == "color-clock-v3" ||
               (kind == "color-clock-v2" && legacy_modulation != nullptr)) {
      ColorClockSnapshot state;
      if (!number(value, "oscillationPhase", state.oscillation_phase) ||
          (kind == "color-clock-v2" &&
           !modulation_clock(value["modulation"], *legacy_modulation)))
        return false;
      out = state;
    } else if (kind == "modulation-clock-v1") {
      ModulationClockSnapshot state;
      if (!modulation_clock(value, state.modulation))
        return false;
      out = state;
    } else if (kind == "spherical-rings-v2") {
      SphericalRingsSnapshot state;
      if (!spatial(value["walk"], state.walk) ||
          !number(value, "phase", state.phase))
        return false;
      out = state;
    } else {
      float phase;
      if (!number(value, "phase", phase))
        return false;
      if (kind == "phase-clock-v1")
        out = PhaseClockSnapshot{phase};
      else
        return false;
    }
    return true;
  }

public:
  /**
   * @brief Encodes @p snapshot as a JS object.
   * @param snapshot Snapshot to encode.
   * @return Object with `schemaVersion`, `animationsPaused`, `chain`,
   *         `parameters` and required `modulation`, plus optional `runtime`
   *         and `paletteBank`.
   */
  static Val encode(const ChainSnapshot &snapshot) {
    auto out = Val::object();
    out.set("schemaVersion", snapshot.schema_version);
    out.set("animationsPaused", snapshot.animations_paused);
    auto chain = Val::array();
    for (size_t index = 0; index < snapshot.chain.size(); ++index) {
      auto entry = Val::object();
      entry.set("instance", snapshot.chain[index].instance);
      entry.set("operator", snapshot.chain[index].operator_id);
      chain.set(index, entry);
    }
    out.set("chain", chain);
    auto parameters = Val::array();
    for (size_t index = 0; index < snapshot.parameters.size(); ++index) {
      auto entry = Val::object();
      entry.set("name", snapshot.parameters[index].name);
      entry.set("value", snapshot.parameters[index].value);
      parameters.set(index, entry);
    }
    out.set("parameters", parameters);
    auto modulation = Val::array();
    for (size_t index = 0; index < snapshot.modulation.size(); ++index) {
      auto entry = Val::object();
      entry.set("instance", snapshot.modulation[index].instance);
      entry.set("recipe", recipe(snapshot.modulation[index].recipe));
      modulation.set(index, entry);
    }
    out.set("modulation", modulation);
    if (snapshot.runtime) {
      auto runtime = Val::array();
      for (size_t index = 0; index < snapshot.runtime->size(); ++index) {
        auto entry = Val::object();
        const auto &source = (*snapshot.runtime)[index];
        entry.set("instance", source.instance);
        entry.set("kind", std::string(Pullback::Interp::runtime_snapshot_kind(
                              source.state)));
        entry.set("state", runtime_state(source.state));
        runtime.set(index, entry);
      }
      out.set("runtime", runtime);
    }
    if (snapshot.palette_bank) {
      const auto &bank = *snapshot.palette_bank;
      auto palette = Val::object();
      palette.set("chroma", bank.chroma);
      auto hues = Val::array();
      auto cycles = Val::array();
      for (size_t index = 0; index < 3; ++index) {
        hues.set(index, bank.hues[index]);
        auto cycle = Val::object();
        cycle.set("frame", bank.cycles[index].frame);
        cycle.set("nextSequence", bank.cycles[index].next_sequence);
        cycle.set("fadeActive", bank.cycles[index].fade_active);
        cycle.set("displayDirty", bank.cycles[index].display_dirty);
        cycles.set(index, cycle);
      }
      palette.set("hues", hues);
      palette.set("cycles", cycles);
      out.set("paletteBank", palette);
    }
    return out;
  }

  /**
   * @brief Decodes an untrusted JS payload into @p out.
   * @details `modulation` is required; `runtime` and `paletteBank` are optional.
   * The first failure is
   * returned; `schemaVersion` is checked before every other field. Within
   * each group (`chain`/`parameters`, `hues`/`cycles`) the arrays' types are
   * checked first, then their lengths, then their entries in wire order.
   * @p out must start default-constructed and is partially filled on failure.
   * @param input JS payload, expected to be a caller-owned clone.
   * @param out Snapshot to fill.
   * @return `APPLIED` when every field decodes (nothing is applied here);
   *         `INVALID_VALUE` for a missing or mistyped field, including a
   *         missing or non-array `modulation`;
   *         `UNSUPPORTED_VERSION` for a `schemaVersion` other than 3 or
   *         `ChainSnapshot::SCHEMA_VERSION`; `INVALID_LENGTH` for an over-cap
   *         `chain`, `parameters`, `modulation` or `runtime` array, or `hues`/`cycles`
   *         arrays not of length 3.
   */
  static Result decode(const Val &input, ChainSnapshot &out) {
    if (!object(input) || !integer(input, "schemaVersion", out.schema_version))
      return Result::INVALID_VALUE;
    const bool legacy = out.schema_version == 3;
    if (!legacy && out.schema_version != ChainSnapshot::SCHEMA_VERSION)
      return Result::UNSUPPORTED_VERSION;
    std::vector<LegacyClock> legacy_clocks;
    if (!boolean(input, "animationsPaused", out.animations_paused))
      return Result::INVALID_VALUE;
    const auto chain = input["chain"];
    const auto parameters = input["parameters"];
    if (!array_type(chain) || !array_type(parameters))
      return Result::INVALID_VALUE;
    if (chain["length"].as<size_t>() > Pullback::Interp::MAX_CHAIN_OPS ||
        parameters["length"].as<size_t>() > Pullback::Interp::MAX_CHAIN_PARAMS)
      return Result::INVALID_LENGTH;
    for (size_t index = 0; index < chain["length"].as<size_t>(); ++index) {
      ChainSnapshot::Entry entry;
      if (!object(chain[index]) ||
          !text(chain[index], "instance", entry.instance) ||
          !text(chain[index], "operator", entry.operator_id))
        return Result::INVALID_VALUE;
      out.chain.push_back(std::move(entry));
    }
    for (size_t index = 0; index < parameters["length"].as<size_t>(); ++index) {
      ChainSnapshot::Parameter entry;
      if (!object(parameters[index]) ||
          !text(parameters[index], "name", entry.name) ||
          !number(parameters[index], "value", entry.value))
        return Result::INVALID_VALUE;
      out.parameters.push_back(std::move(entry));
    }
    const auto modulation = input["modulation"];
    if (!array_type(modulation))
      return Result::INVALID_VALUE;
    if (modulation["length"].as<size_t>() > Pullback::Interp::MAX_CHAIN_OPS)
      return Result::INVALID_LENGTH;
    for (size_t index = 0; index < modulation["length"].as<size_t>(); ++index) {
      ChainSnapshot::ModulationRecipe entry;
      if (!object(modulation[index]) ||
          !text(modulation[index], "instance", entry.instance) ||
          !recipe(modulation[index]["recipe"], entry.recipe))
        return Result::INVALID_VALUE;
      if (legacy && entry.recipe.source == Pattern::SourceKind::DISPLACEMENT)
        return Result::INVALID_VALUE;
      out.modulation.push_back(std::move(entry));
    }
    const auto runtime = input["runtime"];
    if (!runtime.isUndefined()) {
      if (!array_type(runtime))
        return Result::INVALID_VALUE;
      if (runtime["length"].as<size_t>() > Pullback::Interp::MAX_CHAIN_OPS)
        return Result::INVALID_LENGTH;
      out.runtime.emplace();
      for (size_t index = 0; index < runtime["length"].as<size_t>(); ++index) {
        ChainSnapshot::Runtime entry;
        std::string kind;
        Pullback::Modulation::State clock;
        if (!object(runtime[index]) ||
            !text(runtime[index], "instance", entry.instance) ||
            !text(runtime[index], "kind", kind) ||
            !runtime_state(kind, runtime[index]["state"], entry.state,
                           legacy ? &clock : nullptr))
          return Result::INVALID_VALUE;
        if (legacy && kind == "color-clock-v2") {
          if (clock.source == Pattern::SourceKind::DISPLACEMENT)
            return Result::INVALID_VALUE;
          legacy_clocks.push_back({entry.instance, clock});
        }
        out.runtime->push_back(std::move(entry));
      }
    }
    const auto palette = input["paletteBank"];
    if (!palette.isUndefined()) {
      GeneratedPaletteBank::Snapshot bank;
      if (!object(palette) || !number(palette, "chroma", bank.chroma))
        return Result::INVALID_VALUE;
      if (!array_type(palette["hues"]) || !array_type(palette["cycles"]))
        return Result::INVALID_VALUE;
      if (palette["hues"]["length"].as<size_t>() != 3 ||
          palette["cycles"]["length"].as<size_t>() != 3)
        return Result::INVALID_LENGTH;
      for (size_t index = 0; index < 3; ++index) {
        auto holder = Val::object();
        holder.set("value", palette["hues"][index]);
        auto &clock = bank.cycles[index];
        const auto cycle = palette["cycles"][index];
        if (!integer(holder, "value", bank.hues[index]) || !object(cycle) ||
            !integer(cycle, "frame", clock.frame) ||
            !integer(cycle, "nextSequence", clock.next_sequence) ||
            !boolean(cycle, "fadeActive", clock.fade_active) ||
            !boolean(cycle, "displayDirty", clock.display_dirty))
          return Result::INVALID_VALUE;
      }
      out.palette_bank = bank;
    }
    return legacy ? migrate_legacy(out, legacy_clocks) : Result::APPLIED;
  }
};

} // namespace hs_wasm
