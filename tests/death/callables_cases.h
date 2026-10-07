/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Callables death cases.

/**
 * @brief Death case: calling an empty (default-constructed) Fn must trap.
 * @details Host/WASM only: the device Fn backend returns a zero-initialized R
 *          instead of trapping.
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
