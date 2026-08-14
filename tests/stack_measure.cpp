/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host stack high-water-mark gate across every effect.
 *
 * Method: stack painting. Descend writing an address-keyed pattern, unwind (the
 * painted region is now free stack below SP), run one effect (which descends into
 * it), then find the deepest point it reached by paint DENSITY. Painting leaves
 * one small gap per paint() frame, so an untouched window is not perfectly
 * painted either; the deepest CONTROL_BYTES of the region, which no effect comes
 * near, calibrates that gap noise floor and the deepest window above it is the
 * reach. Keying each byte to its address stops a frame that writes runs of one
 * value from reading as unpainted. Built -Os + the device math flags so codegen
 * tracks the size build.
 *
 * The per-effect window is HS_SMOKE_FRAMES (default 8); CI drives the 120-frame
 * window the effects sweep uses, so the deepest call chains an effect only
 * reaches late in its lifecycle are inside the measured region. Frames advance
 * under the injected fixed-cadence clock (hs_test::pin_frame_clock), so a
 * time-driven effect recurses the same way whatever the runner's speed.
 *
 * CI gate: fails (non-zero exit) if the worst effect exceeds the device's DTCM
 * stack reservation (see below).
 */
#include <cstdint>
#include <cstdio>
#include <new>

#include "engine/effects.h"
#include "engine/memory.h"
#include "tests/test_fixture.h"

namespace {

constexpr int W = 288; // Phantasm device canvas
constexpr int H = 144;
const int FRAMES = hs_test::smoke_frames();
constexpr int CHUNK = 2048;
constexpr size_t WIN = 256; // classifier granularity

// Deepest slice of the painted region, taken as untouched: it ends ~380 KB
// below the top, far past any plausible frame.
constexpr size_t CONTROL_BYTES = 65536;

#ifndef HS_DEVICE_STACK_FLOOR_BYTES
#error "HS_DEVICE_STACK_FLOOR_BYTES must be defined by tests/CMakeLists.txt"
#endif

// The device's DTCM stack reservation, read from tools/teensy_budgets.json
// (phantasm.regions.ram1.free_min_bytes) by tests/CMakeLists.txt. That is the
// bytes the Teensy size gate holds free for locals, and the ITCM code ceiling
// is derived from the same figure, so this gate and the device budget cannot
// drift apart.
constexpr size_t BUDGET_BYTES = HS_DEVICE_STACK_FLOOR_BYTES;

volatile uint8_t *g_lo;
int g_measured = 0;
int g_unmeasured = 0;

// Paint value for a byte, keyed to its address so an incidental workload byte
// matches only ~1/256 of the time.
inline uint8_t paint_byte(const volatile uint8_t *a) {
  uintptr_t x = reinterpret_cast<uintptr_t>(a);
  return static_cast<uint8_t>((x ^ (x >> 7)) * 0x2Bu + 0xA5u);
}

// Bytes of a WIN-sized window that no longer hold their paint value.
inline size_t mismatch(const uint8_t *p) {
  size_t m = 0;
  for (size_t i = 0; i < WIN; ++i)
    if (p[i] != paint_byte(p + i))
      ++m;
  return m;
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
  asm volatile("" ::"r"(lo)
               : "memory"); // defeat tail-call / dead-store elision
}

template <typename EffectT>
__attribute__((noinline)) Effect *construct_effect() {
  EffectT *concrete = new (std::nothrow) EffectT();
  HS_CHECK(concrete != nullptr, "stack probe effect allocation failed");
  concrete->init();
  return concrete;
}

// One effect through the smoke_one(test_effects.h) sequence.
template <typename EffectT> __attribute__((noinline)) void run_effect() {
  hs_test::reset_globals();
  hs_test::pin_frame_clock(0);
  Effect *effect = construct_effect<EffectT>();
  for (int f = 0; f < FRAMES; ++f) {
    hs_test::pin_frame_clock(f);
    effect->draw_frame();
    effect->advance_display();
  }
  volatile auto px = effect->get_pixel(0, 0); // keep the render live
  (void)px;
  delete effect;
}

template <typename Effect> size_t measure(const char *name) {
  ++g_measured;
  g_lo = nullptr;
  paint(220); // ~440 KB painted region, then unwind
  volatile uint8_t topmark;
  uint8_t *top = const_cast<uint8_t *>(&topmark);
  run_effect<Effect>();
  uint8_t *lo = const_cast<uint8_t *>(g_lo);
  // Worst untouched window in the control band = the paint-gap noise floor.
  size_t floor_mismatch = 0;
  for (uint8_t *p = lo; p + WIN <= lo + CONTROL_BYTES && p + WIN < top;
       p += WIN) {
    size_t m = mismatch(p);
    if (m > floor_mismatch)
      floor_mismatch = m;
  }
  size_t peak = 0;
  if (floor_mismatch >= WIN / 4) {
    std::printf("  %-22s control band dirty (%zu/%zu B) — reach ran past the "
                "painted region\n",
                name, floor_mismatch, WIN);
  } else {
    // Deepest window damaged past the floor. Control-band windows cannot trip
    // it: the floor is their own maximum.
    for (uint8_t *p = lo; p + WIN < top; p += WIN)
      if (mismatch(p) > floor_mismatch) {
        peak = static_cast<size_t>(top - p);
        break;
      }
    std::printf("  %-22s peak = %6zu B\n", name, peak);
  }
  if (peak == 0)
    ++g_unmeasured;
  return peak;
}

} // namespace

int main() {
  std::printf("=== host stack high-water mark per effect (-Os, x86-64, %dx%d, "
              "%d frames) ===\n",
              W, H, FRAMES);
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
    std::printf(
        "measured %d effects but HS_EFFECT_COUNT = %d — roster empty or "
        "measure() calls dropped\n",
        g_measured, HS_EFFECT_COUNT);
    return 1;
  }
  if (g_unmeasured != 0) {
    std::printf("%d effect(s) reported peak = 0 B — the paint-density probe "
                "found no written window, so the budget below scores nothing\n",
                g_unmeasured);
    return 1;
  }
  const bool over = worst > BUDGET_BYTES;
  std::printf("\nWORST: %s = %zu B (%.1f KB)   device floor %zu B   [%s]\n",
              worst_name, worst, worst / 1024.0, BUDGET_BYTES,
              over ? "FAIL" : "PASS");
  if (over)
    std::printf("  stack budget exceeded — a deep call chain grew; see "
                "tests/stack_measure.cpp header.\n");
  return over ? 1 : 0;
}
