/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdio>
#include <optional>
#include <string>
#include <utility>
#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace hs_test {

// POSIX fd primitives with their underscore-prefixed Windows equivalents.
inline int fd_pipe(int fds[2]) {
#if defined(_WIN32)
  return _pipe(fds, 4096, _O_BINARY);
#else
  return pipe(fds);
#endif
}

inline int fd_dup(int fd) {
#if defined(_WIN32)
  return _dup(fd);
#else
  return dup(fd);
#endif
}

inline int fd_dup2(int from, int to) {
#if defined(_WIN32)
  return _dup2(from, to);
#else
  return dup2(from, to);
#endif
}

inline void fd_close(int fd) {
#if defined(_WIN32)
  _close(fd);
#else
  close(fd);
#endif
}

inline long fd_read(int fd, char *buf, size_t n) {
#if defined(_WIN32)
  return _read(fd, buf, static_cast<unsigned int>(n));
#else
  return static_cast<long>(read(fd, buf, n));
#endif
}

inline int fd_fileno(std::FILE *file) {
#if defined(_WIN32)
  return _fileno(file);
#else
  return fileno(file);
#endif
}

/** @brief Runs a callable with C stdout captured in a temporary file.
 * @return Captured bytes, or nullopt if the capture fails.
 */
template <typename Fn> std::optional<std::string> capture_stdout(Fn &&fn) {
  struct Capture {
    std::FILE *file;
    int saved;
    ~Capture() {
      if (saved >= 0) {
        std::fflush(stdout);
        fd_dup2(saved, 1);
        fd_close(saved);
      }
      if (file)
        std::fclose(file);
    }
  } capture{nullptr, -1};
#if defined(_WIN32)
  if (tmpfile_s(&capture.file) != 0)
    return std::nullopt;
#else
  capture.file = std::tmpfile();
#endif
  if (!capture.file || std::fflush(stdout) != 0)
    return std::nullopt;
  capture.saved = fd_dup(1);
  if (capture.saved < 0 || fd_dup2(fd_fileno(capture.file), 1) < 0)
    return std::nullopt;
  std::forward<Fn>(fn)();
  const bool flushed = std::fflush(stdout) == 0;
  const bool restored = fd_dup2(capture.saved, 1) >= 0;
  fd_close(capture.saved);
  capture.saved = -1;
  if (!flushed || !restored || std::fseek(capture.file, 0, SEEK_SET) != 0)
    return std::nullopt;
  std::string text;
  char buffer[4096];
  for (size_t count;
       (count = std::fread(buffer, 1, sizeof(buffer), capture.file));)
    text.append(buffer, count);
  if (std::ferror(capture.file))
    return std::nullopt;
  return text;
}

} // namespace hs_test
