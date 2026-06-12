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
 * Reentrant: a generator callback may itself call generate(). The full
 * arena reset happens only at the outermost call; a nested call sub-scopes
 * off the caller's live frame (via ScratchScope), so it sees only the
 * scratch headroom above the outer allocations rather than clobbering them.
 *
 * @param target The arena for persistent output (e.g., persistent_arena).
 * @param fn     The generation function or callable.
 * @param args   Extra arguments forwarded to fn.
 * @return       Whatever fn returns.
 */
namespace hs_detail {
// Nesting depth shared across ALL generate() instantiations. It must live
// outside the template: a function-local static inside generate() would give
// each GenerateFn type its own counter, so a nested call (a different lambda
// type) would see depth==0 and reset the arena out from under its caller.
// Single-threaded device, so a plain int needs no synchronization.
inline int &generate_depth() {
  static int depth = 0;
  return depth;
}
} // namespace hs_detail

template <typename GenerateFn, typename... Args>
[[nodiscard]] auto generate(Arena &target, GenerateFn &&fn, Args &&...args) {
  int &depth = hs_detail::generate_depth();
  if (depth == 0) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
  }
  ScratchScope _a(scratch_arena_a);
  ScratchScope _b(scratch_arena_b);
  ++depth;
  struct DepthGuard {
    int &d;
    ~DepthGuard() { --d; }
  } _g{depth};
  return fn(target, scratch_arena_a, scratch_arena_b,
            std::forward<Args>(args)...);
}
