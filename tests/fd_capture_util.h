/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdio>
#if defined(_WIN32)
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace hs_test {

// The fd primitives below differ from POSIX only by the Windows underscore
// prefix; _pipe additionally takes a buffer size and a text/binary mode.
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

inline void fd_dup2(int from, int to) {
#if defined(_WIN32)
  _dup2(from, to);
#else
  dup2(from, to);
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

} // namespace hs_test
