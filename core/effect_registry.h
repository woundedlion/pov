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

/**
 * @brief Resolution-specific fill functions for one registered effect.
 * @details Fill functions are templated per <W,H> but stored as concrete
 *          function pointers, one field per supported resolution.
 */
struct EffectRegistration {
  using FillFn = void(*)(FactoryEntry&); /**< Populates a FactoryEntry for a given resolution. */
  FillFn fill_96_20;   /**< Fill function for the 96x20 resolution. */
  FillFn fill_288_144; /**< Fill function for the 288x144 resolution. */
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
    static std::vector<EffectRegistration> s;
    return s;
  }
  /**
   * @brief Appends a registration to the global table.
   * @param reg Registration record to append.
   * @return Dummy 0, so this can be used as a static-init expression.
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
 * @details Resolutions are enumerated explicitly: a new resolution requires a
 *          new field on EffectRegistration AND a branch here. The static_assert
 *          turns that coupling into a COMPILE error instead of silently
 *          mis-instantiating an unrecognised <W,H>.
 */
template <int W, int H>
constexpr auto get_fill_fn(const EffectRegistration& reg) {
  if constexpr (W == 96 && H == 20) return reg.fill_96_20;
  else if constexpr (W == 288 && H == 144) return reg.fill_288_144;
  else {
    static_assert(W == 96 && H == 20,
                  "get_fill_fn: unsupported <W,H> — add a fill_* field to "
                  "EffectRegistration and a branch here");
    return reg.fill_288_144; // unreachable (static_assert fires)
  }
}

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
    /* used+retain anchor the registrar: nothing references _reg, so under   \
     * LTO / --gc-sections the dynamic initializer (and thus the whole        \
     * registration side effect) could be discarded, silently dropping the    \
     * effect from the registry. `used` keeps the compiler from eliding it    \
     * (and roots it for wasm-ld via llvm.used); `retain` additionally        \
     * survives an ELF linker's --gc-sections. */                             \
    __attribute__((used, retain))                                            \
    static inline int _reg = EffectRegistry::add({                     \
      &fill<96, 20>,                                                   \
      &fill<288, 144>                                                  \
    });                                                                \
  };                                                                   \
  }

#else
// Non-WASM targets (Teensy): no-op, so effect registration pulls in no
// std::vector/std::function overhead and effects are selected statically.
#define REGISTER_EFFECT(ClassName)
#endif
