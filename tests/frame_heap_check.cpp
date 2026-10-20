/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
// Replaces the global allocation functions to count every operator new call
// made while a roster effect draws and presents its frames.
#include "targets/effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <new>

namespace {

bool g_counting = false;
int g_allocations = 0;

void *allocate(std::size_t size) {
  if (g_counting)
    ++g_allocations;
  if (void *p = std::malloc(size ? size : 1))
    return p;
  std::fputs("frame_heap_check: out of memory\n", stderr);
  std::abort();
}

void *allocate_aligned(std::size_t size, std::align_val_t alignment) {
  if (g_counting)
    ++g_allocations;
  const auto align = static_cast<std::size_t>(alignment);
#ifdef _WIN32
  void *p = _aligned_malloc(size ? size : 1, align);
#else
  void *p = std::aligned_alloc(align, (size + align - 1) / align * align);
#endif
  if (p)
    return p;
  std::fputs("frame_heap_check: out of memory\n", stderr);
  std::abort();
}

void release_aligned(void *p) {
#ifdef _WIN32
  _aligned_free(p);
#else
  std::free(p);
#endif
}

/** @brief Counts allocations for the lifetime of the scope. */
struct CountAllocations {
  CountAllocations() {
    g_allocations = 0;
    g_counting = true;
  }
  ~CountAllocations() { g_counting = false; }
};

template <template <int, int> class E, int W = hs_test::effects_tests::SMALL_W,
          int H = hs_test::effects_tests::SMALL_H>
void draw_without_heap(const char *name) {
  HS_CONTEXT(name);
  HS_CONTEXT("width", W);
  hs_test::reset_globals();
  hs_test::pin_frame_clock(0);
  auto effect = std::make_unique<E<W, H>>();
  effect->init();
  int allocations = 0;
  for (int frame = 0; frame < hs_test::smoke_frames(); ++frame) {
    hs_test::pin_frame_clock(frame);
    CountAllocations count;
    effect->draw_frame();
    effect->advance_display();
    allocations += g_allocations;
  }
  if (allocations != 0)
    std::printf("  HEAP %-20s %d allocation(s) over %d frames @ %dx%d\n", name,
                allocations, hs_test::smoke_frames(), W, H);
  HS_EXPECT_EQ(allocations, 0);
}

} // namespace

void *operator new(std::size_t size) { return allocate(size); }
void *operator new[](std::size_t size) { return allocate(size); }
void *operator new(std::size_t size, const std::nothrow_t &) noexcept {
  return allocate(size);
}
void *operator new[](std::size_t size, const std::nothrow_t &) noexcept {
  return allocate(size);
}
void *operator new(std::size_t size, std::align_val_t alignment) {
  return allocate_aligned(size, alignment);
}
void *operator new[](std::size_t size, std::align_val_t alignment) {
  return allocate_aligned(size, alignment);
}
void *operator new(std::size_t size, std::align_val_t alignment,
                   const std::nothrow_t &) noexcept {
  return allocate_aligned(size, alignment);
}
void *operator new[](std::size_t size, std::align_val_t alignment,
                     const std::nothrow_t &) noexcept {
  return allocate_aligned(size, alignment);
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }
void operator delete(void *p, const std::nothrow_t &) noexcept { std::free(p); }
void operator delete[](void *p, const std::nothrow_t &) noexcept {
  std::free(p);
}
void operator delete(void *p, std::align_val_t) noexcept { release_aligned(p); }
void operator delete[](void *p, std::align_val_t) noexcept {
  release_aligned(p);
}
void operator delete(void *p, std::size_t, std::align_val_t) noexcept {
  release_aligned(p);
}
void operator delete[](void *p, std::size_t, std::align_val_t) noexcept {
  release_aligned(p);
}
void operator delete(void *p, std::align_val_t,
                     const std::nothrow_t &) noexcept {
  release_aligned(p);
}
void operator delete[](void *p, std::align_val_t,
                       const std::nothrow_t &) noexcept {
  release_aligned(p);
}

int main() {
  hs_test::ModuleFixture fixture("frame_heap (no operator new while drawing)");
  {
    CountAllocations count;
    void *probe = ::operator new(1);
    ::operator delete(probe);
    HS_EXPECT_EQ(g_allocations, 1);
  }
#define HS_DRAW_WITHOUT_HEAP(name) draw_without_heap<name>(#name);
  HS_EFFECT_LIST(HS_DRAW_WITHOUT_HEAP)
#undef HS_DRAW_WITHOUT_HEAP
  draw_without_heap<ShapeShifter, 288, 144>("ShapeShifter");
  return fixture.result() ? 1 : 0;
}
