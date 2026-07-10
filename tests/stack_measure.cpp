/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host stack high-water-mark gate across every effect.
 *
 * Method: stack painting. Descend writing an address-keyed pattern, unwind (the
 * painted region is now free stack below SP), run one effect (which descends into
 * it), then find the deepest point it reached by paint DENSITY (painting leaves
 * small per-frame gaps, so an untouched window is mostly-painted and a written one
 * is not). Keying each byte to its address stops a frame that writes runs of one
 * value from reading as unpainted. Built -Os + the device math flags so codegen
 * tracks the size build.
 *
 * CI gate: fails (non-zero exit) if the worst effect exceeds BUDGET_BYTES.
 */
#include <cstdint>
#include <cstdio>

#include "engine/effects.h"
#include "engine/memory.h"
#include "tests/test_fixture.h"

namespace {

constexpr int W = 288; // Phantasm device canvas
constexpr int H = 144;
constexpr int FRAMES = 8;
constexpr int CHUNK = 2048;

// Below the 16 KB device stack reservation minus the ISR allowance.
constexpr size_t BUDGET_BYTES = 12288;

volatile uint8_t *g_lo;
int g_measured = 0;

// Paint value for a byte, keyed to its address so an incidental workload byte
// matches only ~1/256 of the time.
inline uint8_t paint_byte(const volatile uint8_t *a) {
  uintptr_t x = reinterpret_cast<uintptr_t>(a);
  return static_cast<uint8_t>((x ^ (x >> 7)) * 0x2Bu + 0xA5u);
}

// Descend painting the address-keyed pattern, then unwind. After return the
// painted region starts at g_lo and sits below the caller's SP (free stack).
__attribute__((noinline)) void paint(int chunks) {
  volatile uint8_t buf[CHUNK];
  for (int i = 0; i < CHUNK; ++i)
    buf[i] = paint_byte(&buf[i]);
  uint8_t *lo = const_cast<uint8_t *>(buf);
  if (!g_lo || lo < g_lo)
    g_lo = lo;
  if (chunks > 1)
    paint(chunks - 1);
  asm volatile("" ::"r"(lo) : "memory"); // defeat tail-call / dead-store elision
}

// One effect through the smoke_one(test_effects.h) sequence.
template <typename Effect> __attribute__((noinline)) void run_effect() {
  hs_test::reset_globals();
  Effect effect;
  effect.init();
  for (int f = 0; f < FRAMES; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  volatile auto px = effect.get_pixel(0, 0); // keep the render live
  (void)px;
}

template <typename Effect> size_t measure(const char *name) {
  ++g_measured;
  g_lo = nullptr;
  paint(220); // ~440 KB painted region, then unwind
  volatile uint8_t topmark;
  uint8_t *top = const_cast<uint8_t *>(&topmark);
  run_effect<Effect>();
  // Deepest reach by paint density: an untouched 256-B window is mostly still
  // painted (paint gaps aside), a workload-written one is not.
  constexpr size_t WIN = 256;
  constexpr size_t WRITTEN_MAX = WIN / 4;
  size_t peak = 0;
  for (uint8_t *p = const_cast<uint8_t *>(g_lo); p + WIN < top; p += WIN) {
    size_t s = 0;
    for (size_t i = 0; i < WIN; ++i)
      if (p[i] == paint_byte(p + i))
        ++s;
    if (s < WRITTEN_MAX) {
      peak = static_cast<size_t>(top - p);
      break;
    }
  }
  std::printf("  %-22s peak = %6zu B\n", name, peak);
  return peak;
}

} // namespace

int main() {
  std::printf("=== host stack high-water mark per effect (-Os, x86-64, %dx%d) ===\n",
              W, H);
  size_t worst = 0;
  const char *worst_name = "";
#define HS_MEASURE_ONE(name)                                                   \
  {                                                                            \
    size_t p = measure<name<W, H>>(#name);                                     \
    if (p > worst) {                                                           \
      worst = p;                                                               \
      worst_name = #name;                                                      \
    }                                                                          \
  }
  HS_EFFECT_LIST(HS_MEASURE_ONE)
#undef HS_MEASURE_ONE
  if (g_measured != HS_EFFECT_COUNT) {
    std::printf("measured %d effects but HS_EFFECT_COUNT = %d — roster empty or "
                "measure() calls dropped\n",
                g_measured, HS_EFFECT_COUNT);
    return 1;
  }
  const bool over = worst > BUDGET_BYTES;
  std::printf("\nWORST: %s = %zu B (%.1f KB)   budget %zu B   [%s]\n", worst_name,
              worst, worst / 1024.0, BUDGET_BYTES, over ? "FAIL" : "PASS");
  if (over)
    std::printf("  stack budget exceeded — a deep call chain grew; see "
                "tests/stack_measure.cpp header.\n");
  return over ? 1 : 0;
}
