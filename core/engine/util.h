/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file util.h
 * @brief Change-gated application of a live parameter value.
 */

// platform.h first: on device it defines NDEBUG, which must be set before
// <cassert> expands the assert macro, or assert-stripping depends on a prior TU
// having pulled in platform.h.
#include "platform/platform.h"
#include <type_traits>

/**
 * @brief Invokes `apply(current)` only when `current` differs from `last`,
 * then latches `last = current`.
 * @details Live-apply a slider value only on change, so a per-frame push does
 *   no work while the slider sits still. The animation setters it usually feeds
 *   (`set_duration`/`set_period`) guard the no-change case themselves, so this
 *   is a work filter, not a correctness gate.
 * @param current The latest parameter value to test against the cached one.
 * @param last The cached value; updated to `current` when they differ.
 * @param apply Callable receiving the new value; run only on change.
 */
template <typename T, typename F>
inline void apply_if_changed(const T &current, T &last, F &&apply) {
  // T must compare with exact equality: a tolerance comparator (e.g. Vector's
  // operator!=) is non-transitive and could re-fire or latch incorrectly.
  static_assert(std::is_arithmetic_v<T>,
                "apply_if_changed requires an exact-equality T (scalar/int)");
  if (current != last) {
    last = current;
    apply(current);
  }
}
