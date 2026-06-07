/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Death tests — the suite's first coverage of the fail-fast philosophy the
 * project markets. An HS_CHECK violation is a deliberate __builtin_trap() that
 * aborts the whole process; the in-process HS_EXPECT_* harness cannot catch it.
 * So each trap is exercised in a CHILD process: the test binary re-exec's
 * itself with HS_DEATH_CASE=<name> (handled in main() before any module runs),
 * runs exactly one trap-triggering case, and the parent asserts the child died
 * by the *specific* trap status — clang lowers __builtin_trap() to an illegal
 * instruction (x86 ud2), so the child dies by SIGILL (POSIX) /
 * STATUS_ILLEGAL_INSTRUCTION (Windows). Asserting that exact status, not merely
 * "nonzero", means an unrelated child crash (a different signal, an ordinary
 * nonzero exit) can no longer masquerade as a passing death test.
 *
 * Cross-platform by design: the child is selected through an inherited env var
 * and spawned with std::system(), so no fork() (absent on Windows) is needed.
 * A control "spawn check" runs first; if the harness cannot re-exec itself
 * (unknown argv[0], a sandbox that blocks process creation), the death tests
 * are SKIPPED with a notice — EXCEPT under CI (the CI env var is set), where a
 * suite that cannot run is a hard FAILURE rather than a silent green skip, so a
 * regression in the fail-fast layer cannot slip through unobserved.
 */
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "tests/test_harness.h"

#include "core/3dmath.h"
#include "core/memory.h"
#include "core/solids.h"
#include "core/spatial.h"
#include "core/static_circular_buffer.h"

#if !defined(_WIN32)
#include <csignal>    // SIGILL — the expected trap signal
#include <sys/wait.h> // WIFSIGNALED / WTERMSIG / WIFEXITED / WEXITSTATUS
#endif

#if defined(_WIN32)
// Declared locally (no <windows.h>) to suppress the WER "stopped working" box
// for the children we intentionally crash. kernel32 is linked by default for a
// Windows-Clang console build. 0x0001|0x0002 = SEM_FAILCRITICALERRORS |
// SEM_NOGPFAULTERRORBOX.
extern "C" __declspec(dllimport) unsigned int __stdcall
SetErrorMode(unsigned int uMode);
#endif

namespace hs_test {
namespace death {

// argv[0] of the running test binary, captured in main(); used to re-exec self.
inline const char *&self_exe() {
  static const char *s = nullptr;
  return s;
}

// Defeat constant-folding so the optimizer can't prove a trap is taken at
// compile time and reshape the case (each case must trap at *run* time).
template <typename T> inline T opaque(T v) {
  volatile T x = v;
  return x;
}

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Memory surface: arena over-allocation.
inline void case_arena_oom() {
  static uint8_t buf[64];
  Arena a(buf, sizeof(buf));
  void *p = a.allocate(opaque<size_t>(1024)); // > capacity -> HS_CHECK
  if (p == reinterpret_cast<void *>(0x1))
    std::printf("x"); // keep the call live
}

// Arena-container surface: fixed-capacity overflow.
inline void case_arena_vector_overflow() {
  static uint8_t buf[256];
  Arena a(buf, sizeof(buf));
  ArenaVector<int> v(a, 2);
  v.push_back(1);
  v.push_back(2);
  v.push_back(opaque(3)); // exceeds capacity -> HS_CHECK
}

// Math-core surface: normalize a degenerate (zero-length) vector.
inline void case_normalize_zero() {
  Vector z{opaque(0.0f), opaque(0.0f), opaque(0.0f)};
  Vector n = z.normalized(); // length < eps -> HS_CHECK
  if (n.x == 42.0f)
    std::printf("x");
}

// Math-core surface: normalize a NaN vector. A NaN coordinate poisons the
// length to NaN, and `NaN >= epsilon` is false, so the normalize guard fires.
// This is the suite's NaN/Inf fault case: it proves a non-finite producer is
// trapped at the math seam rather than silently propagating NaN into geometry.
inline void case_normalize_nan() {
  const float nan = opaque(std::numeric_limits<float>::quiet_NaN());
  Vector bad{nan, opaque(0.0f), opaque(0.0f)};
  Vector n = bad.normalized(); // length is NaN -> HS_CHECK fails
  if (n.x == 42.0f)
    std::printf("x");
}

// Lookup/registry surface: out-of-range solids index.
inline void case_solids_index_oob() {
  const auto &e = Solids::get_entry(opaque<size_t>(Solids::NUM_ENTRIES));
  if (e.name == nullptr)
    std::printf("x");
}

// Registry surface (by name): an unknown solid name has no valid fallback.
inline void case_solids_unknown_name() {
  configure_arenas_default();
  PolyMesh m = Solids::get_by_name(persistent_arena, scratch_arena_a,
                                   scratch_arena_b, "definitely_not_a_solid");
  if (m.vertices.size() == 0x7fff)
    std::printf("x");
}

// Container surface: StaticCircularBuffer index past the live count.
inline void case_circular_buffer_oob() {
  StaticCircularBuffer<int, 4> cb;
  cb.push_back(10);
  cb.push_back(20);
  int v = cb[opaque<size_t>(5)]; // index >= count -> HS_CHECK
  if (v == 42)
    std::printf("x");
}

// Container surface: SpatialHash insert pool exhaustion (sizing bug).
inline void case_spatial_hash_overflow() {
  SpatialHash sh(1.0f);
  const int n = static_cast<int>(SpatialHash::MAX_ENTRIES) + 1;
  for (int i = 0; i < opaque(n); ++i)
    sh.insert(Vector{0.0f, 0.0f, 1.0f}, i); // exceeds MAX_ENTRIES -> HS_CHECK
}

// Config surface: an over-subscribed arena partition.
inline void case_arena_oversubscribed() {
  // Each request alone fits, but the sum exceeds GLOBAL_ARENA_SIZE -> HS_CHECK.
  configure_arenas(opaque(GLOBAL_ARENA_SIZE), opaque<size_t>(1024),
                   opaque<size_t>(1024));
}

struct Case {
  const char *name;
  void (*fn)();
};

inline const Case *all_cases(int &n) {
  static const Case cases[] = {
      {"arena_oom", case_arena_oom},
      {"arena_vector_overflow", case_arena_vector_overflow},
      {"normalize_zero", case_normalize_zero},
      {"normalize_nan", case_normalize_nan},
      {"solids_index_oob", case_solids_index_oob},
      {"solids_unknown_name", case_solids_unknown_name},
      {"circular_buffer_oob", case_circular_buffer_oob},
      {"spatial_hash_overflow", case_spatial_hash_overflow},
      {"arena_oversubscribed", case_arena_oversubscribed},
  };
  n = static_cast<int>(sizeof(cases) / sizeof(cases[0]));
  return cases;
}

// CHILD entry (called from main when HS_DEATH_CASE is set): run exactly one
// case, then return. The case is expected to trap before returning; returning
// means it did NOT trap, so the child exits 0 and the parent flags it. An
// unknown name (the "__spawn_check__" control) also just returns -> exit 0.
inline void run_child_case(const char *name) {
#if defined(_WIN32)
  SetErrorMode(0x0001u | 0x0002u);
#endif
  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i)
    if (std::strcmp(cs[i].name, name) == 0) {
      cs[i].fn();
      return;
    }
}

// Set HS_DEATH_CASE in this process's env; children inherit it via system().
inline void set_case_env(const char *name) {
#if defined(_WIN32)
  _putenv_s("HS_DEATH_CASE", name);
#else
  setenv("HS_DEATH_CASE", name, 1);
#endif
}

// Spawn the test binary as a child running `name`; return its raw system()
// status. Child stdout/stderr are discarded — we care only about the exit code.
inline int spawn_child(const char *name) {
  set_case_env(name);
  std::string cmd;
#if defined(_WIN32)
  // cmd.exe /c strips one outer quote pair, so wrap the whole command again.
  cmd = "\"\"";
  cmd += self_exe();
  cmd += "\" >NUL 2>&1\"";
#else
  cmd = "\"";
  cmd += self_exe();
  cmd += "\" >/dev/null 2>&1";
#endif
  return std::system(cmd.c_str());
}

// Interpret a std::system() return value as "the child died by the SPECIFIC
// trap status". clang lowers __builtin_trap() to an illegal instruction, so a
// fired HS_CHECK kills the child with SIGILL (POSIX) / STATUS_ILLEGAL_INSTRUCTION
// (Windows). Requiring that exact status — rather than any nonzero exit —
// prevents an unrelated crash or an ordinary nonzero return from being misread
// as a passing death test.
#if defined(_WIN32)
// EXCEPTION_ILLEGAL_INSTRUCTION; an unhandled trap sets it as the process exit
// code, which std::system() returns through cmd.exe.
inline constexpr int kTrapStatus = static_cast<int>(0xC000001D);
#endif

inline bool child_trapped(int rc) {
#if defined(_WIN32)
  return rc == kTrapStatus;
#else
  return rc != -1 && WIFSIGNALED(rc) && WTERMSIG(rc) == SIGILL;
#endif
}

// "Child exited cleanly" (exit 0) — used by the control spawn check.
inline bool child_exited_clean(int rc) {
#if defined(_WIN32)
  return rc == 0;
#else
  return rc != -1 && WIFEXITED(rc) && WEXITSTATUS(rc) == 0;
#endif
}

// Are we running under CI? GitHub Actions (and most CI providers) set CI=true.
// Under CI a death suite that cannot run must FAIL loudly, not skip silently.
inline bool in_ci() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const char *ci = std::getenv("CI");
#pragma clang diagnostic pop
  return ci && ci[0] != '\0';
}

// Report an inability to run the death suite. Loud (counts as a failure) under
// CI, quiet skip otherwise.
inline void report_unrunnable(const char *why, int rc) {
  if (in_ci()) {
    std::printf("  [FAIL] death tests: %s (rc=%d, CI=on)\n", why, rc);
    HS_EXPECT_TRUE(false && "death suite must run under CI");
  } else {
    std::printf("  [skip] death tests: %s (rc=%d)\n", why, rc);
  }
}

inline int run_death_tests() {
  auto scope = hs_test::begin_module("death");

  if (!self_exe() || self_exe()[0] == '\0') {
    report_unrunnable("no argv[0] to re-exec", 0);
    return hs_test::end_module(scope);
  }

  // Control: a child given an unknown case must exit cleanly. If it doesn't,
  // this harness can't reliably spawn itself here (e.g. a sandbox) — skip
  // (or FAIL under CI) rather than emit false results for every case.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    report_unrunnable("cannot re-exec self", control);
    set_case_env("");
    return hs_test::end_module(scope);
  }

  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i) {
    int rc = spawn_child(cs[i].name);
    bool trapped = child_trapped(rc);
    HS_EXPECT_TRUE(trapped);
    std::printf("  [%s] trap fires: %-22s (child rc=%d)\n",
                trapped ? "ok" : "FAIL", cs[i].name, rc);
  }

  set_case_env(""); // leave the env clean for anything that runs after us
  return hs_test::end_module(scope);
}

} // namespace death
} // namespace hs_test
