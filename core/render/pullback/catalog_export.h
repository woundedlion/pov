/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "platform/build_features.h"

#if HS_ENABLE_CHAIN_INTERPRETER

#include <charconv>
#include <string>

#include "render/pullback/interpreter.h"

/**
 * @file catalog_export.h
 * @brief Embind-free string builder for the operator catalog JSON, so native
 *        tests golden-pin the exact bytes the tool receives.
 */

namespace Pullback {

namespace Interp {

namespace Detail {

inline void append_json_string(std::string &out, const char *text) {
  out += '"';
  for (const char *cursor = text; *cursor != '\0'; ++cursor) {
    if (*cursor == '"' || *cursor == '\\')
      out += '\\';
    out += *cursor;
  }
  out += '"';
}

inline void append_json_number(std::string &out, float value) {
  char buffer[32];
  const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
  out.append(buffer, result.ptr);
}

inline void append_json_number(std::string &out, size_t value) {
  char buffer[32];
  const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
  out.append(buffer, result.ptr);
}

inline const char *curve_name(FieldCurve curve) {
  switch (curve) {
  case FieldCurve::LERP:
    return "lerp";
  case FieldCurve::LOG_POSITIVE:
    return "log-positive";
  case FieldCurve::SHORTEST_PERIODIC:
    return "shortest-periodic";
  case FieldCurve::SNAP:
    break;
  }
  return "snap";
}

inline void append_block_json(std::string &out, const char *name,
                              const BlockLayout &block) {
  append_json_string(out, name);
  out += ":{\"size\":";
  append_json_number(out, static_cast<size_t>(block.size));
  out += ",\"align\":";
  append_json_number(out, static_cast<size_t>(block.align));
  out += '}';
}

inline void append_param_json(std::string &out, const ParamFieldInfo &field) {
  out += "{\"id\":";
  append_json_string(out, field.id);
  if (field.topology) {
    out += ",\"topology\":true,\"values\":[";
    for (uint8_t value = 0; value < field.enum_count; ++value) {
      if (value > 0)
        out += ',';
      append_json_string(out, field.enum_ids[value]);
    }
    out += "],\"default\":";
    append_json_string(out, field.enum_ids[field.enum_def]);
  } else {
    out += ",\"name\":";
    if (field.display_name != nullptr)
      append_json_string(out, field.display_name);
    else
      out += "null";
    out += ",\"min\":";
    append_json_number(out, field.min);
    out += ",\"max\":";
    append_json_number(out, field.max);
    out += ",\"default\":";
    append_json_number(out, field.def);
    out += ",\"curve\":";
    append_json_string(out, curve_name(field.curve));
  }
  out += '}';
}

} // namespace Detail

/**
 * @brief Appends the catalog JSON — budgets, carriers, and every
 *        operator-table entry — to @p out.
 * @details The authoritative editor-accounting figures ride along as the
 * budgets object. Golden-pinned byte-exact by the native suite.
 */
inline void append_catalog_json(std::string &out) {
  out += "{\"catalog_version\":2,\"budgets\":{\"max_chain_ops\":";
  Detail::append_json_number(out, MAX_CHAIN_OPS);
  out += ",\"arena_bytes\":";
  Detail::append_json_number(out, CHAIN_ARENA_BYTES);
  out += ",\"max_params\":";
  Detail::append_json_number(out, MAX_CHAIN_PARAMS);
  out += ",\"max_instance_id_length\":";
  Detail::append_json_number(out, MAX_INSTANCE_ID);
  out += ",\"per_op_overhead_bytes\":";
  Detail::append_json_number(out, PER_OP_OVERHEAD_BYTES);
  out += ",\"per_param_name_bytes\":";
  Detail::append_json_number(out, PER_PARAM_NAME_BYTES);
  out += "},\"carriers\":[";
  for (size_t index = 0; index < CARRIER_NAMES.size(); ++index) {
    if (index > 0)
      out += ',';
    Detail::append_json_string(out, CARRIER_NAMES[index]);
  }
  out += "],\"operators\":[";
  for (size_t index = 0; index < OPERATOR_TABLE.size(); ++index) {
    const OperatorDescriptor &op = OPERATOR_TABLE[index];
    if (index > 0)
      out += ',';
    out += "{\"id\":";
    Detail::append_json_string(out, op.operator_id);
    out += ",\"name\":";
    Detail::append_json_string(out, op.display_name);
    out += ",\"input\":";
    Detail::append_json_string(out,
                               CARRIER_NAMES[static_cast<size_t>(op.input)]);
    out += ",\"output\":";
    Detail::append_json_string(out,
                               CARRIER_NAMES[static_cast<size_t>(op.output)]);
    out += ",\"blocks\":{";
    Detail::append_block_json(out, "param", op.runtime.param);
    out += ',';
    Detail::append_block_json(out, "prepared", op.runtime.prepared);
    out += ',';
    Detail::append_block_json(out, "state", op.runtime.state);
    out += "},\"params\":[";
    for (uint16_t field = 0; field < op.schema_count; ++field) {
      if (field > 0)
        out += ',';
      Detail::append_param_json(out, op.schema[field]);
    }
    out += "],\"approximate\":";
    out += op.approximate ? "true" : "false";
    out += '}';
  }
  out += "]}";
}

} // namespace Interp

} // namespace Pullback

#endif // HS_ENABLE_CHAIN_INTERPRETER
