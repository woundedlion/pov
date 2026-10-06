/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Callables death fixtures and guard cases.

/**
 * @brief Death case: calling an empty (default-constructed) Fn must trap.
 * @details Concepts surface — hs::inplace_function routes an empty-state call
 *          through ipf_empty_ops::invoke, which fail-fast traps via check_fail
 *          rather than dereferencing the empty buffer (std::function would throw
 *          bad_function_call; the engine builds without exceptions). The
 *          never-taken opaque(false) assignment keeps the optimizer from proving
 *          the function empty and folding the trap at compile time. The non-trap
 *          value semantics (copy/move/empty operator bool) are covered in-process
 *          by tests/test_concepts.h. Host/WASM only: the device Fn backend
 *          returns a zero-initialized R instead of trapping (row 9 of
 *          docs/ledgers/device_host_divergence_ledger.md).
 */
inline void case_empty_fn_call() {
  Fn<int(int), 16> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: invoking an empty FunctionRef must trap.
 * @details Concepts surface — the empty state's thunk diverges through
 *          function_ref_empty_call rather than calling through a null
 *          context. Unlike the Fn trap, this one ships to the device.
 */
inline void case_empty_function_ref_call() {
  FunctionRef<int(int)> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  if (v == 42)
    std::printf("x");
}
