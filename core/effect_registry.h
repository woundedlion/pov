/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Self-registering effect factory. Each effect header calls
 * REGISTER_EFFECT(ClassName) which appends to a global registry
 * at static-init time, eliminating the hand-maintained list in
 * wasm_bridge.cpp.
 *
 * On non-WASM targets (Teensy), REGISTER_EFFECT is a no-op to
 * avoid pulling in std::vector/std::function overhead.
 */
#pragma once

#ifdef __EMSCRIPTEN__

#include <vector>
#include <string_view>
#include <functional>
#include <memory>

class Effect; // forward decl — defined in canvas.h

/**
 * @brief Concrete factory record for one registered effect at a fixed resolution.
 * @details Populated by a registration's fill function; holds the effect's name,
 *          a creator closure that allocates an instance, and its byte size.
 */
struct FactoryEntry {
  std::string_view name;                          /**< Effect class name (string literal). */
  std::function<std::unique_ptr<Effect>()> creator; /**< Allocates a new effect instance. */
  size_t size;                                    /**< sizeof the effect at this resolution, in bytes. */
};

// Single source of truth for the supported render resolutions. Adding a
// resolution is now ONE edit here: the EffectRegistration fields, the
// get_fill_fn dispatch, and the REGISTER_EFFECT fill-pointer list below all
// expand from this X-macro, so they can no longer drift out of sync (a missing
// field used to value-init to a null FillFn that only faulted at table build).
#define HS_RESOLUTIONS(X)                                                      \
  X(96, 20)                                                                    \
  X(288, 144)

/**
 * @brief Resolution-specific fill functions for one registered effect.
 * @details Fill functions are templated per <W,H> but stored as concrete
 *          function pointers, one field per supported resolution — generated
 *          from HS_RESOLUTIONS as `fill_<W>_<H>` (e.g. `fill_96_20`).
 */
struct EffectRegistration {
  using FillFn = void(*)(FactoryEntry&); /**< Populates a FactoryEntry for a given resolution. */
#define HS_REG_FILL_FIELD(W, H) FillFn fill_##W##_##H;
  HS_RESOLUTIONS(HS_REG_FILL_FIELD)
#undef HS_REG_FILL_FIELD
};

/**
 * @brief Global table of registered effects, populated at static-init time.
 * @details Entries are appended by REGISTER_EFFECT. The Meyers-singleton vector
 *          avoids static-init-order issues across translation units.
 */
class EffectRegistry {
public:
  /**
   * @brief Accesses the global registration table.
   * @return Reference to the singleton vector of registrations.
   */
  static std::vector<EffectRegistration>& entries() {
    static std::vector<EffectRegistration> s = [] {
      std::vector<EffectRegistration> v;
      constexpr size_t kReserveHint = 64; // ~2x the HS_EFFECT_LIST roster
      v.reserve(kReserveHint);
      return v;
    }();
    return s;
  }
  /**
   * @brief Appends a registration to the global table.
   * @param reg Registration record to append.
   * @return Dummy 0, so this can be used as a static-init expression.
   * @warning The append order here is the static-init order of the
   *          REGISTER_EFFECT objects, which is link/translation-unit-order
   *          dependent and therefore NOT stable across builds, toolchains, or a
   *          reordered link line. The index of an entry in `entries()` carries
   *          no meaning. Today this is safe — dispatch is by name and the
   *          rosters use the `HS_EFFECT_LIST` X-macro as the single source of
   *          truth — but do NOT introduce a positional, persisted, or
   *          transmitted index over this table: drive any stable ordering from
   *          `HS_EFFECT_LIST` instead, or sort `entries()` by name at use.
   */
  static int add(EffectRegistration reg) {
    entries().push_back(reg);
    return 0; // dummy return for static-init expression
  }
};

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
// Dependent-false constant so a static_assert in a discarded `if constexpr`
// branch only fires when that branch is actually instantiated. A bare
// `static_assert(false)` would be ill-formed even in the taken branches.
template <int> constexpr bool unsupported_resolution = false;

template <int W, int H>
constexpr auto get_fill_fn(const EffectRegistration& reg) {
#define HS_REG_FILL_BRANCH(w, h) \
  if constexpr (W == (w) && H == (h)) return reg.fill_##w##_##h; else
  HS_RESOLUTIONS(HS_REG_FILL_BRANCH)
#undef HS_REG_FILL_BRANCH
  {
    static_assert(unsupported_resolution<W>,
                  "get_fill_fn: unsupported <W,H> — add it to HS_RESOLUTIONS");
    return EffectRegistration::FillFn{}; // unreachable (static_assert fires)
  }
}

// Anchor attribute for the self-registration object. `used` (always available
// on this WASM clang) keeps the compiler from eliding the unreferenced static
// and roots it for wasm-ld via llvm.used; `retain` additionally survives an ELF
// linker's --gc-sections. `retain` is newer (Clang 13+ / GCC 11+), so guard it
// behind __has_attribute and fall back to `used` alone — wasm-ld is not an ELF
// linker, so `used` is what actually anchors the registrar here and the fallback
// stays correct on any toolchain that predates `retain`.
#if defined(__has_attribute) && __has_attribute(retain)
#define HS_REGISTRAR_ANCHOR __attribute__((used, retain))
#else
#define HS_REGISTRAR_ANCHOR __attribute__((used))
#endif

/**
 * @brief Self-registers an effect class with the global EffectRegistry.
 * @param ClassName Effect class template (instantiated per supported <W,H>).
 * @details Defines an anonymous-namespace registrar whose static initializer
 *          appends fill functions for every supported resolution. The
 *          used+retain attributes keep the dynamic initializer from being
 *          discarded under LTO / --gc-sections. On non-WASM targets this macro
 *          expands to nothing.
 */
#define REGISTER_EFFECT(ClassName)                                     \
  namespace {                                                          \
  struct ClassName##_Registrar {                                        \
    template <int W, int H>                                            \
    static void fill(FactoryEntry& e) {                                \
      e.name    = #ClassName;                                          \
      e.creator = []() -> std::unique_ptr<Effect> {                    \
        return std::make_unique<ClassName<W, H>>();                     \
      };                                                               \
      e.size = sizeof(ClassName<W, H>);                                \
    }                                                                  \
    /* HS_REGISTRAR_ANCHOR anchors the registrar: nothing references _reg, so  \
     * under LTO / --gc-sections the dynamic initializer (and thus the whole   \
     * registration side effect) could be discarded, silently dropping the     \
     * effect from the registry. See the macro definition above for the        \
     * used / retain rationale and the pre-`retain` fallback. */               \
    HS_REGISTRAR_ANCHOR                                                       \
    static inline int _reg = EffectRegistry::add({                     \
      HS_RESOLUTIONS(HS_DETAIL_REG_FILL_PTR)                          \
    });                                                                \
  };                                                                   \
  }

// Emits one `&fill<W, H>,` per resolution for the REGISTER_EFFECT initializer
// above. Defined outside the macro (preprocessor directives can't live inside a
// macro body) and left defined because REGISTER_EFFECT expands it in every
// effect translation unit. The trailing comma is harmless in a braced-init list.
#define HS_DETAIL_REG_FILL_PTR(W, H) &fill<W, H>,

#else
// Non-WASM targets (Teensy): no-op, so effect registration pulls in no
// std::vector/std::function overhead and effects are selected statically.
#define REGISTER_EFFECT(ClassName)
#endif
