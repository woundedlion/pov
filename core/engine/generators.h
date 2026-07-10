/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "engine/memory.h"

namespace hs_detail {
/**
 * @brief Returns the shared generate() nesting-depth counter.
 * @return Reference to the process-wide recursion depth (starts at 0).
 * @details The counter is shared across ALL generate() instantiations, so it must
 *   live outside the template: a function-local static inside generate() would give
 *   each GenerateFn type its own counter, so a nested call (a different lambda type)
 *   would see depth==0 and reset the arena out from under its caller.
 */
inline int &generate_depth() {
  static int depth = 0;
  return depth;
}
} // namespace hs_detail

constexpr int MAX_GENERATE_DEPTH = 16;

/**
 * @brief Universal generation wrapper for procedural geometry.
 * @tparam GenerateFn The generation function or callable type.
 * @tparam Args       Types of the extra arguments forwarded to fn.
 * @param target The arena for persistent output (e.g., persistent_arena). Must
 *   not be either engine scratch arena (`scratch_arena_a`/`scratch_arena_b`):
 *   the depth-0 reset and `ScratchScope` rewind would destroy the output.
 * @param fn     The generation function or callable.
 * @param args   Extra arguments forwarded to fn.
 * @return Whatever fn returns.
 * @details Resets and scopes both scratch arenas, then invokes
 *   fn(target, scratch_a, scratch_b, args...). All procedural geometry creation
 *   goes through this wrapper for a deterministic arena lifecycle. The fn signature
 *   must be: ReturnType fn(Arena& target, Arena& scratch_a, Arena& scratch_b, Args...)
 *
 *   Reentrant: a generator callback may itself call generate(). The full arena
 *   reset happens only at the outermost call; a nested call sub-scopes off the
 *   caller's live frame (via ScratchScope), so it sees only the scratch headroom
 *   above the outer allocations rather than clobbering them.
 */
template <typename GenerateFn, typename... Args>
auto generate(Arena &target, GenerateFn &&fn, Args &&...args) {
  HS_CHECK(&target != &scratch_arena_a && &target != &scratch_arena_b,
           "generate: target must not alias an engine scratch arena");
  int &depth = hs_detail::generate_depth();
  if (depth == 0) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
  }
  ScratchScope a_guard(scratch_arena_a);
  ScratchScope b_guard(scratch_arena_b);
  ++depth;
  HS_CHECK(depth <= MAX_GENERATE_DEPTH, "generate: recursion too deep");
  /**
   * @brief RAII guard that decrements the nesting depth on scope exit.
   */
  struct DepthGuard {
    int &d; /**< Reference to the shared nesting-depth counter. */
    /**
     * @brief Decrements the referenced depth counter on destruction.
     */
    ~DepthGuard() { --d; }
  } depth_guard{depth};
  return fn(target, scratch_arena_a, scratch_arena_b,
            std::forward<Args>(args)...);
}
