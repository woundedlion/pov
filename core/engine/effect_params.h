/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_params.h
 * @brief Runtime parameter registry: ParamDef descriptors and the
 *        fixed-capacity ParamList an Effect owns.
 */

#include "engine/build_features.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>

// External ParamDef storage covers the GUI bridge and firmware builds that opt
// in.
#if HS_ENABLE_PARAM_GUI_BRIDGE || defined(HS_EXTERNAL_PARAM_STORAGE)
#define HS_PARAM_EXTERNAL_STORAGE 1
#else
#define HS_PARAM_EXTERNAL_STORAGE 0
#endif

class Effect;

/**
 * @brief Outcome of a named parameter write (Effect::updateParameter, surfaced
 *        to JS by the WASM bridge's setParameter).
 * @details NO_EFFECT is produced only by the bridge, which may have no effect
 *          installed to forward the write to; Effect::updateParameter itself
 *          never returns it.
 */
enum class ParamSetResult {
  APPLIED,       /**< Value written (floats clamped to [min,max] first). */
  NO_EFFECT,     /**< No effect is installed to receive the write. */
  UNKNOWN_PARAM, /**< No registered parameter has this name. */
  READONLY,      /**< Parameter is engine-written telemetry. */
  NON_FINITE,    /**< Value is NaN or infinite. */
};

/**
 * @brief Defines a runtime-adjustable parameter.
 */
struct ParamDef {
  /** @brief Runtime representation of the parameter target. */
  enum class TargetType : uint8_t {
    FLOAT,
    BOOL,
    INT_I8,
    INT_U8,
    INT_I16,
    INT_U16,
    INT_I32,
    INT_U32,
  };

  const char *name = nullptr; /**< Parameter name. */
  void *target = nullptr;     /**< Writable parameter target. */
#if HS_ENABLE_PARAM_GUI_BRIDGE
  const void *display_target = nullptr; /**< Optional live-value source. */
#endif
  const char *const *options = nullptr; /**< Option labels for an enumerated
                             param (GUI dropdown), or null for a plain param.
                             Must outlive the effect (string literals). */
  const char *const *export_options =
      nullptr;          /**< C++ enum literals indexed like options, or null. */
  float min = 0;        /**< Minimum value (for floats). */
  float max = 1;        /**< Maximum value (for floats). */
  int option_count = 0; /**< Number of labels; > 0 marks an enum target. */
  TargetType target_type = TargetType::FLOAT; /**< Target storage format. */
  bool animated = false; /**< True if an animation drives this member; the GUI
                             surfaces these as auto-pausing sliders. */
  bool readonly = false; /**< True if this is engine-written telemetry; the
                             GUI shows it live but disables editing. */
  bool preset = true;    /**< Whether preset exports include this parameter. */

  template <typename Integer> static float get_integer(const void *source) {
    Integer value;
    std::memcpy(&value, source, sizeof(value));
    return static_cast<float>(value);
  }

  template <typename Integer> void set_integer(float value) {
    const Integer stored = static_cast<Integer>(value);
    std::memcpy(target, &stored, sizeof(stored));
  }

  /** @brief Reads one value source as a float (bool maps to 0/1). */
  float get_from(const void *source) const {
    switch (target_type) {
    case TargetType::FLOAT:
      return *static_cast<const float *>(source);
    case TargetType::BOOL:
      return *static_cast<const bool *>(source) ? 1.0f : 0.0f;
    case TargetType::INT_I8:
      return get_integer<int8_t>(source);
    case TargetType::INT_U8:
      return get_integer<uint8_t>(source);
    case TargetType::INT_I16:
      return get_integer<int16_t>(source);
    case TargetType::INT_U16:
      return get_integer<uint16_t>(source);
    case TargetType::INT_I32:
      return get_integer<int32_t>(source);
    case TargetType::INT_U32:
      return get_integer<uint32_t>(source);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Reads the displayed value as float (bool maps to 0/1).
   * @return The display source's value; a bool yields 1.0 or 0.0.
   */
  float get() const {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    return get_from(display_target != nullptr ? display_target : target);
#else
    return get_from(target);
#endif
  }

  /**
   * @brief Reads the writable target value as float (bool maps to 0/1).
   * @return The value a new renderer must adopt, independent of display lerp.
   */
  float get_requested() const { return get_from(target); }

  /**
   * @brief Write a float value (bool threshold at 0.5).
   * @param v Value to store; a bool target is set true when v > 0.5.
   * @warning Raw write: applies no readonly/finite/[min,max] gate. That
   * contract lives solely in Effect::updateParameter — any write from outside
   * trusted engine code must route through there, not call set() on a handle
   * from ParamList::find().
   */
  void set(float v) {
    switch (target_type) {
    case TargetType::FLOAT:
      *static_cast<float *>(target) = v;
      return;
    case TargetType::BOOL:
      *static_cast<bool *>(target) = v > 0.5f;
      return;
    case TargetType::INT_I8:
      return set_integer<int8_t>(v);
    case TargetType::INT_U8:
      return set_integer<uint8_t>(v);
    case TargetType::INT_I16:
      return set_integer<int16_t>(v);
    case TargetType::INT_U16:
      return set_integer<uint16_t>(v);
    case TargetType::INT_I32:
      return set_integer<int32_t>(v);
    case TargetType::INT_U32:
      return set_integer<uint32_t>(v);
    }
    __builtin_unreachable();
  }

  /**
   * @brief Check if this parameter targets a bool.
   * @return True if the target is a bool pointer, false if a float pointer.
   */
  bool is_bool() const { return target_type == TargetType::BOOL; }

  /**
   * @brief Whether the target stores whole numbers, so the GUI steps by one.
   * @details True for both an enumerated target and a plain integer one; the
   * two are told apart by option_count, which only the former sets.
   */
  bool is_integer() const {
    return target_type != TargetType::FLOAT && target_type != TargetType::BOOL;
  }

  /**
   * @brief Check if this parameter is enumerated (renders as a dropdown).
   * @return True if option labels are attached.
   */
  bool is_enum() const { return option_count > 0; }
};
static_assert(sizeof(void *) != 4 ||
                  sizeof(ParamDef) == (HS_ENABLE_PARAM_GUI_BRIDGE ? 36 : 32),
              "ParamDef must keep its 32-bit device footprint");

/**
 * @brief Fixed-capacity registry of an effect's runtime parameters.
 * @details Stack-allocated array (no heap) to uphold the WASM no-realloc
 * memory-view invariant; capacity is enforced at registration time.
 */
struct ParamList {
  // Effect is the sole trusted mutator; the writable accessors below are
  // private so other callers see only the const overloads and route value
  // writes through updateParameter.
  friend class Effect;

  std::array<ParamDef, 32> elements; /**< Default fixed-capacity storage. */

  const ParamDef *data() const {
#if HS_PARAM_EXTERNAL_STORAGE
    return external_elements != nullptr ? external_elements : elements.data();
#else
    return elements.data();
#endif
  }
  size_t capacity() const {
#if HS_PARAM_EXTERNAL_STORAGE
    return external_elements != nullptr ? external_capacity : elements.size();
#else
    return elements.size();
#endif
  }

  /**
   * @brief Const iterator to the first registered parameter.
   * @return Pointer to the first element.
   */
  const ParamDef *begin() const { return data(); }
  /**
   * @brief Const one-past-the-end iterator over registered parameters.
   * @return Pointer just past the last registered element.
   */
  const ParamDef *end() const { return data() + count; }
  /**
   * @brief Looks up a registered parameter by name (the public, read-only
   * lookup).
   * @param name Parameter name to match (exact string compare).
   * @return Const pointer to the matching parameter, or nullptr if not found.
   */
  const ParamDef *find(const char *name) const {
    for (size_t i = 0; i < count; ++i) {
      if (std::strcmp(data()[i].name, name) == 0)
        return &data()[i];
    }
    return nullptr;
  }
  /**
   * @brief Number of registered parameters.
   * @return The count of registered parameters.
   */
  size_t size() const { return count; }
  /** @brief Monotonic token changed whenever the descriptor schema mutates. */
  uint32_t schema_generation() const {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    return schema_gen;
#else
    return 0;
#endif
  }

private:
  // Storage bookkeeping and the writable accessors, reachable only by the
  // friended Effect (see the note at the top of the struct). Kept private so
  // value writes route through updateParameter.
#if HS_PARAM_EXTERNAL_STORAGE
  ParamDef *external_elements = nullptr;
  size_t external_capacity = 0;
#endif
  size_t count = 0; /**< Number of registered parameters. */

  ParamDef *data() {
#if HS_PARAM_EXTERNAL_STORAGE
    return external_elements != nullptr ? external_elements : elements.data();
#else
    return elements.data();
#endif
  }
  ParamDef *begin() { return data(); }
  ParamDef *end() { return data() + count; }
  ParamDef *find(const char *name) {
    return const_cast<ParamDef *>(std::as_const(*this).find(name));
  }
  void bump_schema_generation() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    ++schema_gen;
#endif
  }
#if HS_ENABLE_PARAM_GUI_BRIDGE
  uint32_t schema_gen = 0;
#endif
};
