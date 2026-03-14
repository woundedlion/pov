/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "memory.h"

/**
 * @brief Universal generation wrapper for procedural geometry.
 * @details Resets and scopes both scratch arenas, then invokes:
 *   fn(target, scratch_a, scratch_b, args...)
 *
 * All procedural geometry creation should go through this wrapper to
 * ensure deterministic arena lifecycle. The fn signature must be:
 *   ReturnType fn(Arena& target, Arena& scratch_a, Arena& scratch_b, Args...)
 *
 * @param target The arena for persistent output (e.g., persistent_arena).
 * @param fn     The generation function or callable.
 * @param args   Extra arguments forwarded to fn.
 * @return       Whatever fn returns.
 */
template <typename GenerateFn, typename... Args>
auto generate(Arena &target, GenerateFn &&fn, Args &&...args) {
  scratch_arena_a.reset();
  scratch_arena_b.reset();
  ScopedScratch _a(scratch_arena_a);
  ScopedScratch _b(scratch_arena_b);
  return fn(target, scratch_arena_a, scratch_arena_b,
            std::forward<Args>(args)...);
}
