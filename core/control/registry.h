/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file registry.h
 * @brief Effect factory records generated from the target effect roster.
 * @details Compiled only when HS_ENABLE_EFFECT_REGISTRY is set.
 */

#include "core/platform/build_features.h"
#include "core/platform/effect_identity.h"

#if HS_ENABLE_EFFECT_REGISTRY

#include "core/platform/constants.h"
#include "core/platform/platform.h"
#include <array>
#include <string_view>
#include <functional>
#include <memory>

class Effect;

/** @brief Storage for a concrete effect type identity. */
template <typename T> struct EffectTypeTag {
  static constexpr char id = 0; ///< Its address is the type's identity.
};

/**
 * @brief RTTI-free identity token for one concrete effect type.
 * @tparam T Effect type at a fixed resolution, e.g. ShaderChain<288, 144>.
 * @return An address unique to T for the module's lifetime.
 * @details Proves an Effect base pointer's concrete type before a downcast.
 */
template <typename T> constexpr const void *effect_type_key() {
  return &EffectTypeTag<T>::id;
}

/**
 * @brief Concrete factory record for one registered effect at a fixed resolution.
 * @details identity comes from EffectRegistration; its fill function
 * sets the resolution-specific fields.
 */
struct FactoryEntry {
  /// Maps a preset index to its stable preset ID.
  using PresetIdFn = std::string_view (*)(size_t);

  hs::EffectIdentity
      identity; /**< Canonical IDs, aliases and immutable RNG key. */
  std::function<std::unique_ptr<Effect>()>
      creator; /**< Allocates a new effect instance. */
  const void *type_key =
      nullptr;     /**< effect_type_key() of the type creator() builds. */
  size_t size = 0; /**< sizeof the effect at this resolution, in bytes. */
  size_t preset_count = 0; /**< Number of authored presets. */
  PresetIdFn preset_id =
      nullptr; /**< Stable ID of a preset by index; null without PRESET_IDS. */
};

// Supported render resolutions; per-resolution registration code expands from
// this X-macro.
/**
 * @brief Applies X to each supported render resolution.
 * @param X Macro invoked as X(W, H).
 */
#define HS_RESOLUTIONS(X)                                                      \
  X(96, 20)                                                                    \
  X(288, 144)

/**
 * @brief Asserts one resolution fits MAX_W x MAX_H.
 * @param W Width in pixels.
 * @param H Height in pixels.
 */
#define HS_REG_RESOLUTION_FITS(W, H)                                           \
  static_assert(W <= MAX_W && H <= MAX_H,                                      \
                "HS_RESOLUTIONS entry exceeds MAX_W/MAX_H");
HS_RESOLUTIONS(HS_REG_RESOLUTION_FITS)
#undef HS_REG_RESOLUTION_FITS

/// Resolution whose instantiation supplies the effect identity, as W, H.
#define HS_REG_IDENTITY_RESOLUTION 96, 20

/**
 * @brief Resolution-specific fill functions for one registered effect.
 * @details One `fill_<W>_<H>` function pointer per HS_RESOLUTIONS entry.
 */
struct EffectRegistration {
  hs::EffectIdentity
      identity; /**< Machine/display identity and compatibility keys. */
  using FillFn = void (*)(
      FactoryEntry &); /**< Populates a FactoryEntry for a given resolution. */
/**
 * @brief Declares the fill pointer for one resolution.
 * @param W Width in pixels.
 * @param H Height in pixels.
 */
#define HS_REG_FILL_FIELD(W, H) FillFn fill_##W##_##H;
  HS_RESOLUTIONS(HS_REG_FILL_FIELD)
#undef HS_REG_FILL_FIELD
};

/**
 * @brief Validates the canonical-ID and compatibility-alias lookup namespace.
 * @param entries Registrations to check.
 * @return False on an empty primary identity or a key shared by different effects.
 */
template <size_t N>
constexpr bool
registration_names_unique(const std::array<EffectRegistration, N> &entries) {
  for (size_t i = 0; i < N; ++i) {
    const auto &identity = entries[i].identity;
    if (identity.display_name.empty() || identity.machine_id.empty())
      return false;
    for (size_t j = 0; j < i; ++j)
      for (std::string_view key : identity.accepted_ids())
        if (!key.empty() && entries[j].identity.matches(key))
          return false;
  }
  return true;
}

/**
 * @brief Traps unless registration_names_unique() holds.
 * @tparam N Registration count.
 * @param entries Registrations to check.
 */
template <size_t N>
void validate_effect_registrations(
    const std::array<EffectRegistration, N> &entries) {
  HS_CHECK(registration_names_unique(entries),
           "duplicate effect registration identity");
}

// Dependent condition delays the diagnostic until an unsupported resolution
// is instantiated.
/// Always false; a dependent static_assert condition.
template <int> constexpr bool unsupported_resolution = false;

/**
 * @brief Selects the fill function pointer matching the given <W,H>.
 * @tparam W Frame width in pixels.
 * @tparam H Frame height in pixels.
 * @param reg Registration holding one fill pointer per supported resolution.
 * @return The fill function pointer for <W,H>.
 * @details An unlisted <W,H> is a compile error.
 */
template <int W, int H>
constexpr auto get_fill_fn(const EffectRegistration &reg) {
#define HS_REG_FILL_BRANCH(w, h)                                               \
  if constexpr (W == (w) && H == (h))                                          \
    return reg.fill_##w##_##h;                                                 \
  else
  HS_RESOLUTIONS(HS_REG_FILL_BRANCH)
#undef HS_REG_FILL_BRANCH
  {
    static_assert(unsupported_resolution<W>,
                  "get_fill_fn: unsupported <W,H> — add it to HS_RESOLUTIONS");
    return EffectRegistration::FillFn{}; // unreachable (static_assert fires)
  }
}

/**
 * @brief Populates resolution-specific fields of one factory entry.
 * @param entry Entry to populate.
 */
template <template <int, int> class ClassName, int W, int H>
void fill_registration(FactoryEntry &entry) {
  entry.creator = []() -> std::unique_ptr<Effect> {
    return std::make_unique<ClassName<W, H>>();
  };
  entry.type_key = effect_type_key<ClassName<W, H>>();
  entry.size = sizeof(ClassName<W, H>);
  if constexpr (requires { ClassName<W, H>::PRESET_IDS; }) {
    entry.preset_count = ClassName<W, H>::PRESET_IDS.size();
    entry.preset_id = [](size_t index) -> std::string_view {
      const auto &ids = ClassName<W, H>::PRESET_IDS;
      return index < ids.size() ? ids[index] : std::string_view{};
    };
  } else if constexpr (requires { ClassName<W, H>::authored_preset_count(); }) {
    entry.preset_count = ClassName<W, H>::authored_preset_count();
  }
}

/**
 * @brief Builds resolution fill pointers for one HS_EFFECT_LIST entry.
 * @tparam ClassName Registered effect family.
 * @return The registration.
 */
template <template <int, int> class ClassName>
constexpr EffectRegistration make_registration() {
  static_assert(
      requires { ClassName<HS_REG_IDENTITY_RESOLUTION>::EFFECT_IDENTITY; },
      "registered effects require a typed EFFECT_IDENTITY");
  constexpr auto identity =
      ClassName<HS_REG_IDENTITY_RESOLUTION>::EFFECT_IDENTITY;
  static_assert(hs::canonical_effect_id(identity.machine_id));
  static_assert(!identity.rng_seed_key.empty());
  static_assert(ClassName<HS_REG_IDENTITY_RESOLUTION>::EFFECT_ID ==
                identity.machine_id);
#define HS_REG_FILL_POINTER(W, H) &fill_registration<ClassName, W, H>,
  return {identity, HS_RESOLUTIONS(HS_REG_FILL_POINTER)};
#undef HS_REG_FILL_POINTER
}

#endif
