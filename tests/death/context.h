/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Death tests for the fail-fast (HS_CHECK / __builtin_trap) seams, plus
 * cross-process cold/warm effect determinism.
 *
 * Each case runs in a re-exec'd child (HS_DEATH_CHILD=harness,
 * HS_DEATH_CASE=<name>). The parent requires the illegal-instruction trap
 * status (SIGILL / STATUS_ILLEGAL_INSTRUCTION) and an HS_CHECK breadcrumb
 * naming the case's exact guard: under
 * -fsanitize-trap=undefined, UB lowers to the same instruction.
 */
#pragma once

/** @file
 * @brief Fail-fast case fixtures and definitions. */

#include "core/animation/orientation.h"
#include <array>
#include <cerrno>
#include <chrono>
#include <thread>
#include <bit>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <optional>
#include <string>
#include <utility>

#include "death_guard_sites.h" // generated HS_CHECK census
#include "tests/test_fixture.h"
#include "tests/test_pullback.h"
#include "tests/effects/fixture.h"
#include "tests/test_effects.h"
#include "tests/test_shapeshifter_oracle.h"
#include "tests/shader_chain_fixture.h"
#include "tests/test_hyper_lattice.h"
#include "tests/prepared_face_test_access.h"
#include "tests/test_harness.h"

#include "core/math/3dmath.h"
#include "core/math/4dmath.h"
#include "core/animation/animation.h"
#include "core/animation/carousel.h"
#include "core/render/canvas.h"
#include "core/color/color.h"
#include "core/color/modulation.h"
#include "core/pattern/noise.h"
#include "core/pattern/sphere_cells.h"
#include "core/pattern/sphere_lut.h"
#include "core/math/geometry.h"
#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/control/choreography.h"
#include "core/control/presets.h"
#include "core/control/registry.h"
#include "core/memory.h"
#include "core/platform/led.h"
#include "core/math/lenses.h"
#include "core/mesh/hankin.h"
#include "core/mesh/mesh.h"
#include "core/mesh/recipe.h"
#include "core/render/plot.h"
#include "core/render/pullback/interpreter.h"
#include "core/render/pullback/catalog_export.h"
#include "core/render/scan.h"
#include "core/render/scan/radial_surface.h"
#include "core/render/scan/radial_surface_candidate_mask.h"
#include "core/render/sdf.h"
#include "core/render/sdf/volume.h"
#include "core/mesh/solids.h"
#include "core/math/spherical_field.h"
#include "core/math/spherical_harmonics.h"
#include "core/spatial/kd_tree.h"
#include "core/spatial/reaction_graph.h"
#include "core/containers/triangular_bitset.h"
#include "core/containers/static_circular_buffer.h"
#include "core/animation/transformer.h"
#include "hardware/pov_sync.h"
#include "targets/wasm/param_marshal.h"
#include "targets/wasm/document_parameters.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <fcntl.h>    // open / O_WRONLY for the capture file
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#include <unistd.h>   // fork / execv / dup2 / close / _exit — shell-free spawn
#else
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <process.h> // _getpid
#endif

namespace hs_test::death_tests {

/**
 * @brief Accessor for the running test binary path, used to re-exec self.
 * @return Reference to the static char-pointer slot; captured in main().
 */
inline const char *&self_exe() {
  static const char *s = nullptr;
  return s;
}

/**
 * @brief Captures the executable path for child processes.
 * @param argv0 Invocation path, used on POSIX.
 */
inline void set_self_exe(const char *argv0) {
#if defined(_WIN32)
  (void)argv0;
  static std::string path;
  path.resize(MAX_PATH);
  for (;;) {
    const DWORD length = GetModuleFileNameA(nullptr, path.data(),
                                            static_cast<DWORD>(path.size()));
    if (length == 0) {
      self_exe() = nullptr;
      return;
    }
    if (length < path.size()) {
      path.resize(length);
      self_exe() = path.c_str();
      return;
    }
    path.resize(path.size() * 2);
  }
#else
  self_exe() = argv0;
#endif
}

/**
 * @brief Launders a value through a volatile to defeat constant-folding.
 * @tparam T Value type to pass through opaquely.
 * @param v Value to make opaque to the optimizer.
 * @return A copy of @p v the optimizer cannot prove constant.
 * @details Keeps the compiler from proving a trap is taken at compile time and
 *          reshaping the case; each case must trap at run time.
 */
template <typename T> inline T opaque(T v) {
  volatile T x = v;
  return x;
}

/** @brief Records successful admission of the death case's valid companion. */
inline void death_admitted() { std::fputs("death: admitted\n", stderr); }
} // namespace hs_test::death_tests
