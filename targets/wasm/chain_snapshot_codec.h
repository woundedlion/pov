/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <emscripten/val.h>
#include <limits>
#include "core/render/pullback/interpreter.h"
#include "workbench/shader/chain_snapshot.h"

namespace hs_wasm {

class ChainSnapshotCodec {
  using Val = emscripten::val;
  using Result = ChainSnapshotRestoreResult;
  using Runtime = Pullback::Interp::RuntimeSnapshot;

  static bool array(const Val &value, size_t max_count) {
    return Val::global("Array").call<bool>("isArray", value) &&
           value["length"].as<size_t>() <= max_count;
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
    out.set("rawOrientation", quaternion(state.raw_orientation));
    out.set("angularVelocity", state.angular_velocity);
    out.set("spinPhase", state.spin_phase);
    out.set("legacy", state.legacy);
    return out;
  }
  static bool spatial(const Val &value,
                      Pullback::Interp::SpatialWalkSnapshot &state) {
    return object(value) && integer(value, "noiseSeed", state.noise_seed) &&
           integer(value, "walkTime", state.walk_time) &&
           vector(value["position"], state.position) &&
           vector(value["direction"], state.direction) &&
           quaternion(value["wander"], state.wander) &&
           quaternion(value["rawOrientation"], state.raw_orientation) &&
           number(value, "angularVelocity", state.angular_velocity) &&
           number(value, "spinPhase", state.spin_phase) &&
           boolean(value, "legacy", state.legacy);
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
          } else if constexpr (std::is_same_v<State, ColorClockSnapshot>) {
            out.set("oscillationPhase", state.oscillation_phase);
            out.set("hueNoisePhase", state.hue_noise_phase);
            out.set("hueNoiseSeed", state.hue_noise_seed);
          } else if constexpr (std::is_same_v<State, SphericalRingsSnapshot>) {
            out.set("walk", spatial(state.walk));
            out.set("phase", state.phase);
          } else if constexpr (!std::is_same_v<State, std::monostate>)
            out.set("phase", state.phase);
          return out;
        },
        snapshot);
  }
  static bool runtime_state(const std::string &kind, const Val &value,
                            Runtime &out) {
    using namespace Pullback::Interp;
    if (!object(value))
      return false;
    if (kind == "spatial-walk-v1") {
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
    } else if (kind == "affine-clock-v1") {
      AffineClockSnapshot state;
      if (!number(value, "phase", state.phase) ||
          !number(value, "rotation", state.rotation))
        return false;
      out = state;
    } else if (kind == "color-clock-v1") {
      ColorClockSnapshot state;
      if (!number(value, "oscillationPhase", state.oscillation_phase) ||
          !number(value, "hueNoisePhase", state.hue_noise_phase) ||
          !integer(value, "hueNoiseSeed", state.hue_noise_seed))
        return false;
      out = state;
    } else if (kind == "spherical-rings-v1") {
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
      else if (kind == "ripple-clock-v1")
        out = RippleClockSnapshot{phase};
      else
        return false;
    }
    return true;
  }

public:
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

  static Result decode(const Val &input, ChainSnapshot &out) {
    if (!object(input) || !integer(input, "schemaVersion", out.schema_version))
      return Result::INVALID_VALUE;
    if (out.schema_version != 1)
      return Result::UNSUPPORTED_VERSION;
    if (!boolean(input, "animationsPaused", out.animations_paused))
      return Result::INVALID_VALUE;
    const auto chain = input["chain"];
    const auto parameters = input["parameters"];
    if (!array(chain, Pullback::Interp::MAX_CHAIN_OPS) ||
        !array(parameters, Pullback::Interp::MAX_CHAIN_PARAMS))
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
    const auto runtime = input["runtime"];
    if (!runtime.isUndefined()) {
      if (!array(runtime, Pullback::Interp::MAX_CHAIN_OPS))
        return Result::INVALID_LENGTH;
      out.runtime.emplace();
      for (size_t index = 0; index < runtime["length"].as<size_t>(); ++index) {
        ChainSnapshot::Runtime entry;
        std::string kind;
        if (!object(runtime[index]) ||
            !text(runtime[index], "instance", entry.instance) ||
            !text(runtime[index], "kind", kind) ||
            !runtime_state(kind, runtime[index]["state"], entry.state))
          return Result::INVALID_VALUE;
        out.runtime->push_back(std::move(entry));
      }
    }
    const auto palette = input["paletteBank"];
    if (!palette.isUndefined()) {
      GeneratedPaletteBank::Snapshot bank;
      if (!object(palette) || !number(palette, "chroma", bank.chroma) ||
          !array(palette["hues"], 3) ||
          palette["hues"]["length"].as<size_t>() != 3 ||
          !array(palette["cycles"], 3) ||
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
    return Result::APPLIED;
  }
};

} // namespace hs_wasm
