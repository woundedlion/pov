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

// Helper: select the correct fill function pointer for a given <W,H>
template <int W, int H>
constexpr auto get_fill_fn(const EffectRegistration& reg) {
  if constexpr (W == 96 && H == 20) return reg.fill_96_20;
  else                               return reg.fill_288_144;
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
    static inline int _reg = EffectRegistry::add({                     \
      &fill<96, 20>,                                                   \
      &fill<288, 144>                                                  \
    });                                                                \
  };                                                                   \
  }

#else
// Non-WASM targets: no-op
#define REGISTER_EFFECT(ClassName)
#endif
