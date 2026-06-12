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

struct FactoryEntry {
  std::string_view name;
  std::function<std::unique_ptr<Effect>()> creator;
  size_t size;
};

// Each registration stores resolution-specific fill functions.
// Fill functions are templated per <W,H> but stored as concrete
// function pointers for each supported resolution.
struct EffectRegistration {
  using FillFn = void(*)(FactoryEntry&);
  FillFn fill_96_20;
  FillFn fill_288_144;
};

// Global table of registered effects, populated at static-init time by
// REGISTER_EFFECT. Meyers-singleton vector avoids static-init-order issues.
class EffectRegistry {
public:
  static std::vector<EffectRegistration>& entries() {
    static std::vector<EffectRegistration> s;
    return s;
  }
  static int add(EffectRegistration reg) {
    entries().push_back(reg);
    return 0; // dummy return for static-init expression
  }
};

// Helper: select the correct fill function pointer for a given <W,H>.
// The registry stores one fill pointer per supported resolution, so this must
// enumerate them explicitly: a new resolution requires a new field above AND a
// branch here. The static_assert makes that coupling a COMPILE error instead of
// silently mis-instantiating an unrecognised <W,H> at the 288x144 fill.
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
