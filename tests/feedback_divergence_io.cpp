/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#define _CRT_SECURE_NO_WARNINGS
#include "tests/test_filter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <vector>

namespace {
int write_count = 0;
int close_count = 0;
int failed_write = 0;
bool failed_close = false;

size_t fault_fwrite(const void *data, size_t size, size_t count, FILE *stream) {
  ++write_count;
  if (write_count == failed_write)
    return std::fwrite(data, size, count - 1, stream);
  return std::fwrite(data, size, count, stream);
}

int fault_fclose(FILE *stream) {
  ++close_count;
  const int status = std::fclose(stream);
  return failed_close ? EOF : status;
}
} // namespace

#define fwrite fault_fwrite
#define fclose fault_fclose
#define main feedback_divergence_main
#include "tests/feedback_divergence.cpp"
#undef main
#undef fclose
#undef fwrite

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  for (int failure = 1; failure <= 6; ++failure) {
    write_count = 0;
    close_count = 0;
    failed_write = failure <= 5 ? failure : 0;
    failed_close = failure == 6;
    const int status = dump(argv[1], 1);
    if (status == 0 || std::filesystem::exists(argv[1]) ||
        write_count != (failed_write ? failed_write : 5) || close_count != 1) {
      std::fprintf(stderr,
                   "failure %d: status=%d, partial dump retained=%d, "
                   "writes=%d, closes=%d\n",
                   failure, status, std::filesystem::exists(argv[1]),
                   write_count, close_count);
      std::remove(argv[1]);
      return 1;
    }
  }
  return 0;
}
