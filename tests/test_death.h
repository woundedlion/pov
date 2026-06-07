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
 * abnormally.
 *
 * Cross-platform by design: the child is selected through an inherited env var
 * and spawned with std::system(), so no fork() (absent on Windows) is needed.
 * A control "spawn check" runs first; if the harness cannot re-exec itself
 * (unknown argv[0], a sandbox that blocks process creation), the death tests
 * are SKIPPED with a notice rather than reported as failures.
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

#if !defined(_WIN32)
#include <sys/wait.h> // WIFSIGNALED / WIFEXITED / WEXITSTATUS for system() status
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

// Lookup/registry surface: out-of-range solids index.
inline void case_solids_index_oob() {
  const auto &e = Solids::get_entry(opaque<size_t>(Solids::NUM_ENTRIES));
  if (e.name == nullptr)
    std::printf("x");
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
      {"solids_index_oob", case_solids_index_oob},
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

// Interpret a std::system() return value as "the child terminated abnormally".
inline bool child_died(int rc) {
#if defined(_WIN32)
  return rc != 0; // a trap -> nonzero (illegal-instruction) exit code
#else
  if (rc == -1)
    return false; // system() itself failed
  if (WIFSIGNALED(rc))
    return true; // killed by SIGILL/SIGABRT (the trap)
  return WIFEXITED(rc) && WEXITSTATUS(rc) != 0;
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

inline int run_death_tests() {
  auto scope = hs_test::begin_module("death");

  if (!self_exe() || self_exe()[0] == '\0') {
    std::printf("  [skip] death tests: no argv[0] to re-exec\n");
    return hs_test::end_module(scope);
  }

  // Control: a child given an unknown case must exit cleanly. If it doesn't,
  // this harness can't reliably spawn itself here (e.g. a sandbox) — SKIP
  // rather than emit false failures for every case.
  int control = spawn_child("__spawn_check__");
  if (!child_exited_clean(control)) {
    std::printf("  [skip] death tests: cannot re-exec self (control rc=%d)\n",
                control);
    set_case_env("");
    return hs_test::end_module(scope);
  }

  int n;
  const Case *cs = all_cases(n);
  for (int i = 0; i < n; ++i) {
    int rc = spawn_child(cs[i].name);
    bool died = child_died(rc);
    HS_EXPECT_TRUE(died);
    std::printf("  [%s] trap fires: %-22s (child rc=%d)\n",
                died ? "ok" : "FAIL", cs[i].name, rc);
  }

  set_case_env(""); // leave the env clean for anything that runs after us
  return hs_test::end_module(scope);
}

} // namespace death
} // namespace hs_test
