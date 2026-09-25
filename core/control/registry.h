/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file registry.h
 * @brief Effect factory records generated from the target effect roster.
 * @details Enabled for the WASM build and registry tests. On firmware it is a
 *          no-op, so effect registration pulls in no std::vector/std::function
 *          overhead.
 */

#include "platform/build_features.h"

#if HS_ENABLE_EFFECT_REGISTRY

#include "platform/constants.h"
#include "platform/platform.h"
#include <array>
#include <string_view>
#include <functional>
#include <memory>

class Effect; // forward decl — defined in canvas.h

/**
 * @brief RTTI-free identity token for one concrete effect type.
 * @tparam T Effect type at a fixed resolution, e.g. Shader<288, 144>.
 * @return An address unique to T for the module's lifetime.
 * @details Lets a holder of an Effect base pointer prove which concrete type the
 *          factory built before downcasting to it, without RTTI and without
 *          trusting a name string.
 */
template <typename T> struct EffectTypeTag {
  static constexpr char id = 0;
};

template <typename T> constexpr const void *effect_type_key() {
  return &EffectTypeTag<T>::id;
}

/**
 * @brief Concrete factory record for one registered effect at a fixed resolution.
 * @details Populated by a registration's fill function; holds the effect's name,
 *          a creator closure that allocates an instance, the concrete type key
 *          the creator produces, and its byte size.
 */
struct FactoryEntry {
  using PresetIdFn = std::string_view (*)(size_t);

  std::string_view name;      /**< Effect class name (string literal). */
  std::string_view stable_id; /**< Persisted effect identity: the class's
                                 EFFECT_ID when it declares one, else its class
                                 name. It feeds stable_effect_seed(), so an
                                 effect without an EFFECT_ID changes its
                                 persisted id and its seed when renamed. */
  std::function<std::unique_ptr<Effect>()>
      creator; /**< Allocates a new effect instance. */
  const void *type_key =
      nullptr;     /**< effect_type_key() of the type creator() builds. */
  size_t size = 0; /**< sizeof the effect at this resolution, in bytes. */
  size_t preset_count = 0; /**< Number of authored presets. */
  PresetIdFn preset_id =
      nullptr; /**< Registry-only stable preset lookup, when declared. */
};

// Single source of truth for the supported render resolutions. Adding a resolution
// is ONE edit here: the EffectRegistration fields, the get_fill_fn dispatch, and
// the factory fill-pointer list below all expand from this X-macro, so they
// cannot drift out of sync.
#define HS_RESOLUTIONS(X)                                                      \
  X(96, 20)                                                                    \
  X(288, 144)

// Every listed resolution must fit the framebuffers the Effect constructor
// bounds, or its factory would only fail once invoked.
#define HS_REG_RESOLUTION_FITS(W, H)                                           \
  static_assert(W <= MAX_W && H <= MAX_H,                                      \
                "HS_RESOLUTIONS entry exceeds MAX_W/MAX_H");
HS_RESOLUTIONS(HS_REG_RESOLUTION_FITS)
#undef HS_REG_RESOLUTION_FITS

// The instantiation stable-id queries are answered from. EFFECT_ID does not
// depend on <W,H>, so every resolution names the same id; spelling one out
// keeps the persisted identity — and the per-effect RNG stream behind it — off
// HS_RESOLUTIONS' ordering.
#define HS_REG_IDENTITY_RESOLUTION 96, 20

/**
 * @brief Resolution-specific fill functions for one registered effect.
 * @details Fill functions are templated per <W,H> but stored as concrete
 *          function pointers, one field per supported resolution — generated
 *          from HS_RESOLUTIONS as `fill_<W>_<H>` (e.g. `fill_96_20`).
 */
struct EffectRegistration {
  std::string_view name; /**< Effect class/header stem. */
  using FillFn = void (*)(
      FactoryEntry &); /**< Populates a FactoryEntry for a given resolution. */
#define HS_REG_FILL_FIELD(W, H) FillFn fill_##W##_##H;
  HS_RESOLUTIONS(HS_REG_FILL_FIELD)
#undef HS_REG_FILL_FIELD
  std::string_view stable_id{}; /**< Persisted identity shared by every
                                        resolution's FactoryEntry. */
};

/** @brief Validates the shared class-name and stable-ID lookup namespace. */
template <size_t N>
constexpr bool
registration_names_unique(const std::array<EffectRegistration, N> &entries) {
  for (size_t i = 0; i < N; ++i) {
    if (entries[i].name.empty() || entries[i].stable_id.empty())
      return false;
    for (size_t j = 0; j < i; ++j)
      if (entries[i].name == entries[j].name ||
          entries[i].stable_id == entries[j].stable_id ||
          entries[i].name == entries[j].stable_id ||
          entries[i].stable_id == entries[j].name)
        return false;
  }
  return true;
}

template <size_t N>
void validate_effect_registrations(
    const std::array<EffectRegistration, N> &entries) {
  HS_CHECK(registration_names_unique(entries),
           "duplicate effect registration identity");
}

// Dependent-false constant so a static_assert in a discarded `if constexpr`
// branch only fires when that branch is actually instantiated. A bare
// `static_assert(false)` would be ill-formed even in the taken branches.
template <int> constexpr bool unsupported_resolution = false;

/**
 * @brief Selects the fill function pointer matching the given <W,H>.
 * @tparam W Frame width in pixels.
 * @tparam H Frame height in pixels.
 * @param reg Registration holding one fill pointer per supported resolution.
 * @return The fill function pointer for <W,H>.
 * @details Resolutions are enumerated from HS_RESOLUTIONS: each generates one
 *          `if constexpr` branch below, and the trailing static_assert turns an
 *          unlisted <W,H> into a COMPILE error instead of silently
 *          mis-instantiating an unrecognised resolution.
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

/** @brief Populates one concrete factory entry from its roster name. */
template <template <int, int> class ClassName, int W, int H>
void fill_registration(FactoryEntry &entry) {
  entry.stable_id = hs::stable_effect_id<ClassName<W, H>>(entry.name);
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

/** @brief Builds resolution fill pointers for one HS_EFFECT_LIST entry. */
template <template <int, int> class ClassName>
constexpr EffectRegistration make_registration(std::string_view name) {
#define HS_REG_FILL_POINTER(W, H) &fill_registration<ClassName, W, H>,
  return {name, HS_RESOLUTIONS(HS_REG_FILL_POINTER)
                    hs::stable_effect_id<ClassName<HS_REG_IDENTITY_RESOLUTION>>(
                        name)};
#undef HS_REG_FILL_POINTER
}

#endif
