/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file document_parameters.h
 * @brief Typed compiled-document modulation access and admission. */

#include "core/control/params.h"
#include "core/render/pullback/operators.h"
#include <cmath>
#include <string_view>

namespace hs_wasm::document_parameters {
/// Canonical operator schema and typed member addresses.
using Params = Pullback::Interp::Op::ModulationParams;

/** @brief Reads the canonical spelling of a compiled enum value.
 * @param field Canonical enum field metadata.
 * @param index Compiled enum value index.
 * @return Registered enum spelling. */
inline const char *enum_id(const Pullback::Interp::ParamFieldInfo &field,
                           uint8_t index) {
  HS_CHECK(index < field.enum_count, "invalid compiled document enum");
  return field.enum_ids[index];
}

/** @brief Admission reason for a canonical compiled-document field.
 * @param config Current compiled recipe.
 * @param index Field index in the operator schema.
 * @return Active source, inactive source, inactive variant or fixed topology. */
inline const char *reason(const ColorModulation::Config &config, size_t index) {
  const auto &field = Params::EXTRA_SCHEMA[index];
  if (field.enum_count || field.topology)
    return "fixed-topology";
  if (Params::active(config, index))
    return "active-source";
  return config.source == Pattern::SourceKind::NOISE &&
                 std::string_view(field.id) == "pattern-contrast"
             ? "inactive-variant"
             : "inactive-source";
}

/** @brief Whether a canonical numeric field belongs to the selected recipe.
 * @param config Current compiled recipe.
 * @param index Field index in the operator schema.
 * @return Whether the field accepts document writes. */
inline bool writable(const ColorModulation::Config &config, size_t index) {
  const auto &field = Params::EXTRA_SCHEMA[index];
  return !field.enum_count && !field.topology && Params::active(config, index);
}

/** @brief Writes one active numeric member through the canonical schema.
 * @param config Candidate recipe; unchanged on refusal.
 * @param id Exact operator field id.
 * @param value Requested float32 value.
 * @return Applied, unknown, immutable, nonfinite or inadmissible. */
inline ParamSetResult write(ColorModulation::Config &config,
                            std::string_view id, float value) {
  for (size_t i = 0; i < Params::EXTRA_SCHEMA.size(); ++i) {
    const auto &field = Params::EXTRA_SCHEMA[i];
    if (id != field.id)
      continue;
    if (!writable(config, i))
      return ParamSetResult::READONLY;
    if (!std::isfinite(value))
      return ParamSetResult::NON_FINITE;
    if (value < field.min || value > field.max)
      return ParamSetResult::INADMISSIBLE;
    Params candidate;
    candidate.modulation = config;
    *static_cast<float *>(candidate.extra_address(i)) = value;
    if (!Pullback::Modulation::admitted(candidate.modulation))
      return ParamSetResult::INADMISSIBLE;
    config = candidate.modulation;
    return ParamSetResult::APPLIED;
  }
  return ParamSetResult::UNKNOWN_PARAM;
}
} // namespace hs_wasm::document_parameters
