/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Death child processes, coverage and shard records. */

#include "tests/death/registry.h"

namespace hs_test::death_tests {

/**
 * @brief Path the spawned child's stdout/stderr is redirected to.
 * @return Stable NUL-terminated path, built once per process.
 * @details Named from the parent's pid so two test binaries running
 *          concurrently cannot overwrite each other's capture.
 */
inline const char *child_capture_path() {
  static char path[512];
  if (path[0] == '\0') {
#if defined(_WIN32)
    const char *dir = hs_test::env_var("TEMP");
    if (!dir || dir[0] == '\0')
      dir = ".";
    std::snprintf(path, sizeof(path), "%s\\hs_death_%d.out", dir, _getpid());
#else
    const char *dir = hs_test::env_var("TMPDIR");
    if (!dir || dir[0] == '\0')
      dir = "/tmp";
    std::snprintf(path, sizeof(path), "%s/hs_death_%d.out", dir,
                  static_cast<int>(getpid()));
#endif
  }
  return path;
}

/**
 * @brief Loads the tail of the capture file into child_output().
 * @details Keeps the last CHILD_OUTPUT_CAP-1 bytes, where a trapping child's
 *          breadcrumb lands.
 */
inline void load_child_output() {
  char *buf = child_output();
  buf[0] = '\0';
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *f = std::fopen(child_capture_path(), "rb");
#pragma clang diagnostic pop
  if (!f)
    return;
  if (std::fseek(f, 0, SEEK_END) == 0) {
    const long size = std::ftell(f);
    if (size > 0) {
      const long CAP = static_cast<long>(CHILD_OUTPUT_CAP) - 1;
      const long keep = size < CAP ? size : CAP;
      if (std::fseek(f, size - keep, SEEK_SET) == 0)
        buf[std::fread(buf, 1, static_cast<size_t>(keep), f)] = '\0';
    }
  }
  std::fclose(f);
}

/**
 * @brief Tests whether captured child output carries one guard's breadcrumb.
 * @param out Text captured from the child.
 * @param file Repo-relative path of the guard source.
 * @param text Expected prefix of the "(condition) message" tail.
 * @return The matched source line, or zero when no breadcrumb matches.
 */
inline int breadcrumb_names_guard(const char *out, const char *file,
                                  const char *text) {
  std::string normalized(out);
  std::replace(normalized.begin(), normalized.end(), '\\', '/');
  const char *prefix = "HS_CHECK failed: ";
  const size_t file_len = std::strlen(file);
  const size_t text_len = std::strlen(text);
  for (const char *p = std::strstr(normalized.c_str(), prefix); p;
       p = std::strstr(p + 1, prefix)) {
    const char *path = p + std::strlen(prefix);
    const char *q = std::strchr(path, '\n');
    const char *end = q ? q : path + std::strlen(path);
    for (q = path; q < end; ++q) {
      if (*q != ':' || q + 1 == end || q[1] < '0' || q[1] > '9')
        continue;
      const size_t path_len = static_cast<size_t>(q - path);
      if (path_len < file_len ||
          std::strncmp(q - file_len, file, file_len) != 0 ||
          (path_len > file_len &&
           q[-static_cast<ptrdiff_t>(file_len) - 1] != '/'))
        break;
      ++q;
      int line = 0;
      while (q < end && *q >= '0' && *q <= '9') {
        line = line * 10 + (*q - '0');
        ++q;
      }
      if (end - q >= 2 && q[0] == ':' && q[1] == ' ' &&
          static_cast<size_t>(end - q - 2) >= text_len &&
          std::strncmp(q + 2, text, text_len) == 0)
        return line;
      break;
    }
  }
  return 0;
}

/**
 * @brief Sets HS_DEATH_CASE in this process's environment.
 * @param name Case selector to publish; the spawned child inherits it through
 *             the environment.
 */
inline void set_case_env(const char *name) {
#if defined(_WIN32)
  _putenv_s("HS_DEATH_CASE", name);
  _putenv_s("HS_DEATH_CHILD", name[0] ? "harness" : "");
#else
  setenv("HS_DEATH_CASE", name, 1);
  setenv("HS_DEATH_CHILD", name[0] ? "harness" : "", 1);
#endif
}

#if defined(_WIN32)
/** @brief Whether the latest child raised an unhandled illegal instruction. */
inline bool &child_unhandled_illegal_instruction() {
  static bool observed = false;
  return observed;
}
#endif

/**
 * @brief Spawns the test binary as a child running the given death case.
 * @param name Case selector passed to the child via HS_DEATH_CASE.
 * @param timeout_ms Maximum child runtime in milliseconds.
 * @return The child's raw status (debugged process on Windows, fork+execv wait status
 *         on POSIX). -1 on a spawn failure, a wait failure, or a timeout
 *         (the child is killed).
 * @details Child stdout/stderr are redirected to child_capture_path() and the
 *          tail is loaded into child_output(), so the caller can require the
 *          HS_CHECK breadcrumb of the guard the case is supposed to fire.
 */
inline int spawn_child(const char *name, unsigned timeout_ms = 10000) {
  set_case_env(name);
  const char *capture_path = child_capture_path();
  // A spawn that fails before the redirect must leave an EMPTY capture, not
  // the last child's breadcrumb.
  child_output()[0] = '\0';
  std::remove(capture_path);
#if defined(_WIN32)
  child_unhandled_illegal_instruction() = false;
  SECURITY_ATTRIBUTES security{sizeof(SECURITY_ATTRIBUTES), nullptr, TRUE};
  HANDLE capture = CreateFileA(capture_path, GENERIC_WRITE,
                               FILE_SHARE_READ | FILE_SHARE_WRITE, &security,
                               CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (capture == INVALID_HANDLE_VALUE)
    return -1;
  STARTUPINFOA startup{};
  startup.cb = sizeof(startup);
  startup.dwFlags = STARTF_USESTDHANDLES;
  startup.hStdOutput = capture;
  startup.hStdError = capture;
  startup.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
  PROCESS_INFORMATION process{};
  std::string command = std::string("\"") + self_exe() + "\"";
  const BOOL started =
      CreateProcessA(self_exe(), command.data(), nullptr, nullptr, TRUE,
                     DEBUG_ONLY_THIS_PROCESS | CREATE_NO_WINDOW, nullptr,
                     nullptr, &startup, &process);
  CloseHandle(capture);
  if (!started)
    return -1;
  int rc = -1;
  bool timed_out = false;
  bool stopping = false;
  bool initial_breakpoint = true;
  auto deadline =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout_ms);
  for (;;) {
    const auto now = std::chrono::steady_clock::now();
    if (now >= deadline) {
      if (stopping) {
        DebugActiveProcessStop(process.dwProcessId);
        WaitForSingleObject(process.hProcess, 5000);
        break;
      }
      timed_out = true;
      stopping = true;
      TerminateProcess(process.hProcess, 1);
      deadline = now + std::chrono::seconds(5);
    }
    const auto remaining =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            deadline - std::chrono::steady_clock::now())
            .count();
    DEBUG_EVENT event{};
    if (!WaitForDebugEvent(
            &event, static_cast<DWORD>(std::max<int64_t>(1, remaining)))) {
      if (GetLastError() == ERROR_SEM_TIMEOUT)
        continue;
      TerminateProcess(process.hProcess, 1);
      DebugActiveProcessStop(process.dwProcessId);
      WaitForSingleObject(process.hProcess, 5000);
      break;
    }
    DWORD continuation = DBG_CONTINUE;
    bool exited = false;
    switch (event.dwDebugEventCode) {
    case CREATE_PROCESS_DEBUG_EVENT:
      if (event.u.CreateProcessInfo.hFile)
        CloseHandle(event.u.CreateProcessInfo.hFile);
      break;
    case LOAD_DLL_DEBUG_EVENT:
      if (event.u.LoadDll.hFile)
        CloseHandle(event.u.LoadDll.hFile);
      break;
    case EXCEPTION_DEBUG_EVENT: {
      const auto &exception = event.u.Exception;
      const DWORD code = exception.ExceptionRecord.ExceptionCode;
      if (initial_breakpoint && exception.dwFirstChance &&
          code == EXCEPTION_BREAKPOINT) {
        initial_breakpoint = false;
      } else {
        continuation = DBG_EXCEPTION_NOT_HANDLED;
        if (!exception.dwFirstChance && code == EXCEPTION_ILLEGAL_INSTRUCTION)
          child_unhandled_illegal_instruction() = true;
      }
      break;
    }
    case EXIT_PROCESS_DEBUG_EVENT:
      if (!stopping)
        rc = static_cast<int>(event.u.ExitProcess.dwExitCode);
      exited = true;
      break;
    }
    if (!ContinueDebugEvent(event.dwProcessId, event.dwThreadId,
                            continuation)) {
      TerminateProcess(process.hProcess, 1);
      DebugActiveProcessStop(process.dwProcessId);
      WaitForSingleObject(process.hProcess, 5000);
      rc = -1;
      break;
    }
    if (exited) {
      if (WaitForSingleObject(process.hProcess, 5000) != WAIT_OBJECT_0)
        rc = -1;
      break;
    }
  }
  CloseHandle(process.hThread);
  CloseHandle(process.hProcess);
  load_child_output();
  if (timed_out)
    std::fprintf(stderr, "death child timed out after %u ms: %s\n", timeout_ms,
                 name);
  return rc;
#else
  // Shell-free spawn: fork + execv with stdout/stderr redirected to the
  // capture file; returns the raw wait status.
  std::fflush(stdout);
  std::fflush(stderr);
  const char *exe = self_exe();
  pid_t pid = fork();
  if (pid < 0)
    return -1;
  if (pid == 0) {
    int capture =
        open(capture_path, O_WRONLY | O_CREAT | O_EXCL | O_TRUNC, 0600);
    if (capture >= 0) {
      dup2(capture, 1);
      dup2(capture, 2);
      if (capture > 2)
        close(capture);
    }
    const char *argv[] = {exe, nullptr};
    execv(exe, const_cast<char *const *>(argv));
    _exit(127); // exec failed — never returns to the harness
  }
  int status = 0;
  const auto deadline_value =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout_ms);
  for (;;) {
    const pid_t waited = waitpid(pid, &status, WNOHANG);
    if (waited == pid)
      break;
    if (waited < 0 && errno != EINTR)
      return -1;
    if (std::chrono::steady_clock::now() >= deadline_value) {
      kill(pid, SIGKILL);
      while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {
      }
      load_child_output();
      std::fprintf(stderr, "death child timed out after %u ms: %s\n",
                   timeout_ms, name);
      return -1;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  }
  load_child_output();
  return status;
#endif
}

#if defined(_WIN32)
/**
 * @brief The Windows EXCEPTION_ILLEGAL_INSTRUCTION process exit code.
 * @details Accepted only with a second-chance illegal-instruction debug event.
 */
inline constexpr int TRAP_STATUS = static_cast<int>(0xC000001D);
#endif

/**
 * @brief Tests whether a child died from an observable illegal instruction.
 * @param rc The raw spawn_child() return value to interpret.
 * @return True iff the child died by SIGILL or an unhandled Windows trap.
 */
inline bool child_trapped(int rc) {
#if defined(_WIN32)
  return rc == TRAP_STATUS && child_unhandled_illegal_instruction();
#else
  return rc != -1 && WIFSIGNALED(rc) && WTERMSIG(rc) == SIGILL;
#endif
}

/**
 * @brief Tests whether the child exited cleanly (exit code 0).
 * @param rc The raw spawn_child() return value to interpret.
 * @return True iff the child exited normally with status 0.
 */
inline bool child_exited_clean(int rc) {
#if defined(_WIN32)
  return rc == 0;
#else
  return rc != -1 && WIFEXITED(rc) && WEXITSTATUS(rc) == 0;
#endif
}

/**
 * @brief Reports that the death suite could not run.
 * @param why Human-readable reason the suite is unrunnable.
 * @param rc The associated child return code, for diagnostics.
 * @details An unrunnable death tier is a test failure on every host.
 */
inline void report_unrunnable(const char *why, int rc) {
  std::printf("  [FAIL] death tests: %s (rc=%d)\n", why, rc);
  HS_EXPECT_TRUE(false && "death suite must run");
}

/**
 * @brief Counts distinct fired guard lines in @p file.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param file Repo-relative source path.
 * @param lines Matched source line for each successfully trapped case, or zero.
 * @return Number of distinct covered source lines.
 */
inline int pinned_guards_in(const Case *cs, int n, const char *file,
                            const int *lines) {
  int pinned = 0;
  for (int i = 0; i < n; ++i) {
    if (!lines[i] || std::strcmp(cs[i].guard_file, file) != 0)
      continue;
    bool duplicate = false;
    for (int j = 0; j < i && !duplicate; ++j)
      duplicate =
          std::strcmp(cs[j].guard_file, file) == 0 && lines[j] == lines[i];
    if (!duplicate)
      ++pinned;
  }
  return pinned;
}

/**
 * @brief Prints what fraction of the engine's fail-fast surface is pinned.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param lines Matched source line for each successfully trapped case, or zero.
 * @details The denominator is the generated fail-fast site census
 *          (HS_CHECK, HS_AUDIT_CHECK, HS_ASSERT_BLOCK_ALIVE and check_fail).
 *          The numerator is the distinct source lines observed in trapped cases. A case naming a file outside the census
 *          fails the module; the ratio and per-file gaps are not gated.
 */
inline void report_guard_coverage(const Case *cs, int n, const int *lines) {
  int covered = 0;
  int off_census = 0;
  constexpr int GAPS = 5;
  const GuardSiteCount *worst[GAPS] = {};
  int worst_gap[GAPS] = {};
  for (const GuardSiteCount &f : GUARD_SITE_COUNTS) {
    int pinned = pinned_guards_in(cs, n, f.file, lines);
    HS_EXPECT_LE(pinned, f.sites);
    covered += pinned;
    const int gap = f.sites - pinned;
    for (int slot = 0; slot < GAPS; ++slot) {
      if (gap <= worst_gap[slot])
        continue;
      for (int k = GAPS - 1; k > slot; --k) {
        worst[k] = worst[k - 1];
        worst_gap[k] = worst_gap[k - 1];
      }
      worst[slot] = &f;
      worst_gap[slot] = gap;
      break;
    }
  }
  for (int i = 0; i < n; ++i) {
    bool in_census = false;
    for (const GuardSiteCount &f : GUARD_SITE_COUNTS)
      if (std::strcmp(cs[i].guard_file, f.file) == 0) {
        in_census = true;
        break;
      }
    if (!in_census)
      ++off_census;
  }
  std::printf(
      "  guard coverage: %d/%d fail-fast sites pinned by a case (%d%%), "
      "%d case(s) outside the census\n",
      covered, GUARD_SITE_TOTAL,
      GUARD_SITE_TOTAL ? covered * 100 / GUARD_SITE_TOTAL : 0, off_census);
  std::printf("  widest gaps:");
  for (int slot = 0; slot < GAPS && worst[slot]; ++slot)
    std::printf(" %s %d/%d", worst[slot]->file,
                worst[slot]->sites - worst_gap[slot], worst[slot]->sites);
  std::printf("\n");
  HS_EXPECT_EQ(off_census, 0);
}

/** @brief The slice of the death module one process runs. */
struct ShardSpec {
  /** @brief What the process does with the case table. */
  enum class Role {
    WHOLE,  /**< Every case, then the coverage report. */
    SHARD,  /**< Cases with table index % count == index, then a record. */
    CENSUS, /**< No cases; merges the count records, then the coverage report. */
  };
  Role role = Role::WHOLE;   /**< This process's role. */
  int index = 0;             /**< This shard's index, in [0, count). */
  int count = 1;             /**< Number of shards the table is split across. */
  const char *dir = nullptr; /**< Directory holding the shard records. */
};

/**
 * @brief Reads the shard selection from HS_DEATH_SHARD and HS_DEATH_SHARD_DIR.
 * @param out Set to the parsed selection.
 * @return False when HS_DEATH_SHARD is malformed or HS_DEATH_SHARD_DIR is
 *         missing for a sharded role.
 * @details An unset or empty HS_DEATH_SHARD selects Role::WHOLE; "k/N" selects
 *          shard k of N and "census/N" the merge of N shard records.
 */
inline bool shard_spec_from_env(ShardSpec &out) {
  const char *spec = hs_test::env_var("HS_DEATH_SHARD");
  const char *dir = hs_test::env_var("HS_DEATH_SHARD_DIR");
  out = ShardSpec{};
  if (!spec || spec[0] == '\0')
    return true;
  const char *slash = std::strchr(spec, '/');
  if (!slash || !dir || dir[0] == '\0')
    return false;
  out.dir = dir;
  if (std::strncmp(spec, "census/", 7) == 0) {
    char *end = nullptr;
    const long count = std::strtol(slash + 1, &end, 10);
    if (end == slash + 1 || *end != '\0' || count < 1 || count > 1024)
      return false;
    out.count = static_cast<int>(count);
    out.role = ShardSpec::Role::CENSUS;
    return true;
  }
  if (!hs_test::parse_shard_fraction(spec, 1024, out.index, out.count))
    return false;
  out.role = ShardSpec::Role::SHARD;
  return true;
}

/**
 * @brief Path of one shard's record file.
 * @param dir Directory holding the shard records.
 * @param index Shard index.
 * @return The path.
 */
inline std::string shard_record_path(const char *dir, int index) {
  return std::string(dir) + "/death_shard_" + std::to_string(index) + ".txt";
}

/**
 * @brief Writes one shard's results for the census to merge.
 * @param shard The shard that ran.
 * @param n Number of cases in the table.
 * @param ran Whether this shard ran each case.
 * @param lines Matched source line for each trapped case, or zero.
 * @details The header line is "hs-death-shard <index> <count> <n>"; each
 *          following line is "<case index> <line>" for one case this shard ran.
 */
inline void write_shard_record(const ShardSpec &shard, int n,
                               const std::vector<bool> &ran,
                               const std::vector<int> &lines) {
  const std::string path = shard_record_path(shard.dir, shard.index);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *f = std::fopen(path.c_str(), "w");
#pragma clang diagnostic pop
  HS_EXPECT_TRUE(f != nullptr);
  if (!f)
    return;
  std::fprintf(f, "hs-death-shard %d %d %d\n", shard.index, shard.count, n);
  for (int i = 0; i < n; ++i)
    if (ran[static_cast<size_t>(i)])
      std::fprintf(f, "%d %d\n", i, lines[static_cast<size_t>(i)]);
  HS_EXPECT_EQ(std::fclose(f), 0);
}

/**
 * @brief Parses whitespace-separated decimal integers filling a record row.
 * @param row NUL-terminated row text.
 * @param out Receives the parsed values.
 * @param count Number of integers the row must hold.
 * @return True iff the row holds exactly @p count integers.
 */
inline bool parse_record_row(const char *row, int *out, int count) {
  for (int k = 0; k < count; ++k) {
    char *end = nullptr;
    const long value = std::strtol(row, &end, 10);
    if (end == row)
      return false;
    out[k] = static_cast<int>(value);
    row = end;
  }
  while (*row == ' ' || *row == '\r' || *row == '\n')
    ++row;
  return *row == '\0';
}

/**
 * @brief Merges every shard's record into one per-case line table.
 * @param shard The census selection; its count names the records to read.
 * @param cs The case table.
 * @param n Number of cases in it.
 * @param lines Set to the matched source line for each case, or zero.
 * @details Fails unless every shard's record is present and well formed and
 *          the shards together ran each enabled case exactly once. The records
 *          are deleted after reading, so a stale record never feeds a later
 *          census.
 */
inline void merge_shard_records(const ShardSpec &shard, const Case *cs, int n,
                                int *lines) {
  std::vector<int> runs(static_cast<size_t>(n), 0);
  int problems = 0;
  for (int k = 0; k < shard.count; ++k) {
    const std::string path = shard_record_path(shard.dir, k);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    std::FILE *f = std::fopen(path.c_str(), "r");
#pragma clang diagnostic pop
    if (!f) {
      std::printf("  [FAIL] shard %d/%d left no record at %s\n", k, shard.count,
                  path.c_str());
      ++problems;
      continue;
    }
    constexpr char HEADER[] = "hs-death-shard ";
    char row[128] = "";
    int header[3] = {-1, -1, -1};
    if (!std::fgets(row, sizeof(row), f) ||
        std::strncmp(row, HEADER, sizeof(HEADER) - 1) != 0 ||
        !parse_record_row(row + sizeof(HEADER) - 1, header, 3) ||
        header[0] != k || header[1] != shard.count || header[2] != n) {
      std::printf("  [FAIL] shard record %s has a bad header: %s\n",
                  path.c_str(), row);
      ++problems;
    } else {
      while (std::fgets(row, sizeof(row), f)) {
        int fields[2] = {-1, 0};
        const bool parsed = parse_record_row(row, fields, 2);
        const int i = fields[0];
        const int line = fields[1];
        if (!parsed || i < 0 || i >= n || i % shard.count != k) {
          std::printf("  [FAIL] shard record %s has a bad row: %s\n",
                      path.c_str(), row);
          ++problems;
          continue;
        }
        ++runs[static_cast<size_t>(i)];
        lines[i] = line;
      }
    }
    std::fclose(f);
    std::remove(path.c_str());
  }
  for (int i = 0; i < n; ++i) {
    const int expected = case_enabled(cs[i]) ? 1 : 0;
    if (runs[static_cast<size_t>(i)] != expected) {
      std::printf("  [FAIL] case %s ran in %d shard(s), expected %d\n",
                  cs[i].name, runs[static_cast<size_t>(i)], expected);
      ++problems;
    }
  }
  std::printf("  shard census: %d record(s) merged\n", shard.count);
  HS_EXPECT_EQ(problems, 0);
}

/**
 * @brief Compares per-effect frame folds from two fresh child processes.
 * @details Each child folds every effect cold and then warm; both folds must
 *          match within a child and across the two children. Their complete
 *          records must fit in CHILD_OUTPUT_CAP. Each child renders the whole
 *          effect roster, so it gets a longer timeout than a trap case.
 */
inline void compare_determinism_folds() {
  constexpr unsigned TIMEOUT_MS = 60000;
  const int determinism_a = spawn_child(DETERMINISM_PROBE_CASE, TIMEOUT_MS);
  const std::string fold_a = child_output();
  const int determinism_b = spawn_child(DETERMINISM_PROBE_CASE, TIMEOUT_MS);
  const std::string fold_b = child_output();
  if (!child_exited_clean(determinism_a))
    std::printf("  determinism child A logged: %s\n", fold_a.c_str());
  if (!child_exited_clean(determinism_b))
    std::printf("  determinism child B logged: %s\n", fold_b.c_str());
  HS_EXPECT_TRUE(child_exited_clean(determinism_a));
  HS_EXPECT_TRUE(child_exited_clean(determinism_b));
  constexpr size_t EFFECT_COUNT = static_cast<size_t>(HS_EFFECT_COUNT);
  constexpr size_t PREFIX_BYTES = sizeof("capture ") - 1;
  constexpr size_t PAYLOAD_BYTES = 33;
  size_t record_count = 0;
  size_t offset = 0;
  while ((offset = fold_a.find("capture ", offset)) != std::string::npos) {
    HS_CONTEXT("effect", record_count);
    offset += PREFIX_BYTES;
    if (offset + PAYLOAD_BYTES >= fold_a.size()) {
      HS_EXPECT_TRUE(false);
      break;
    }
    const std::string cold = fold_a.substr(offset, 16);
    const std::string warm = fold_a.substr(offset + 17, 16);
    HS_EXPECT_EQ(cold.find_first_not_of("0123456789abcdef"), std::string::npos);
    HS_EXPECT_EQ(warm.find_first_not_of("0123456789abcdef"), std::string::npos);
    HS_EXPECT_EQ(fold_a[offset + 16], ' ');
    HS_EXPECT_TRUE(fold_a[offset + PAYLOAD_BYTES] == '\r' ||
                   fold_a[offset + PAYLOAD_BYTES] == '\n');
    HS_EXPECT_TRUE(cold != "0000000000000000");
    HS_EXPECT_EQ(cold, warm);
    ++record_count;
    offset += PAYLOAD_BYTES;
  }
  HS_EXPECT_EQ(record_count, EFFECT_COUNT);
  HS_EXPECT_EQ(fold_a, fold_b);
}

} // namespace hs_test::death_tests
