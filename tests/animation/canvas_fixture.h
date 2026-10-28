/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Stand-in Canvas reference
// ----------------------------------------------------------------------------
// The shared Effect and Canvas stay live for the whole module. No test in
// this module may construct another Effect or Canvas.
// ============================================================================

/**
 * @brief Storage for the module-scoped fake-canvas pointer.
 * @return Reference to the pointer run_animation_tests() sets to the shared
 * fixture for the module's duration and clears (fixture destroyed) on exit.
 */
inline Canvas *&fake_canvas_ptr() {
  static Canvas *p = nullptr;
  return p;
}

/**
 * @brief Returns the module-scoped shared Canvas reused across the tests.
 * @return Reference to the fixture owned by run_animation_tests(); valid for the
 * duration of the module but never dereferenced by the animations under test.
 */
inline Canvas &fake_canvas() { return *fake_canvas_ptr(); }
