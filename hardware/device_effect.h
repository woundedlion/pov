/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file device_effect.h
 * @brief Foreground effect construction for device drivers. */

#include <new>
#include <utility>
#include "core/math/geometry.h"
#include "core/memory.h"

/** @brief Prepares geometry and arenas, then constructs and initializes an effect.
 * @tparam E Effect type.
 * @tparam Args Constructor argument types.
 * @param args Arguments forwarded to the effect constructor.
 * @return Constructed effect owned by the caller. */
template <typename E, typename... Args>
E *construct_device_effect(Args &&...args) {
  math::GeometryResolution<E>::init();
  configure_arenas_default();
  E *effect = new (std::nothrow) E(std::forward<Args>(args)...);
  HS_CHECK(effect != nullptr, "effect allocation failed (OOM)");
  effect->init();
  return effect;
}
