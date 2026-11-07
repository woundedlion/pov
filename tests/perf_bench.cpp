/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host render micro-benchmark: per-frame wall time for every effect.
 * Usage: perf_bench [--galaxies-draw on|off full|segment]
 *        perf_bench --temporal-sample EFFECT on|off full|segment
 *        perf_bench --temporal-cost
 *
 * Build the two targets (perf_bench_Os at -Os = the Phantasm size config,
 * perf_bench_O3 at -O3) and diff to see what optimization buys the render. Host
 * x86-64 wall time is a relative proxy for the Teensy device, not an absolute.
 * Not a CTest.
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <cstdint>
#include <cstdio>

#include "targets/effects.h"
#include "core/memory.h"
#include "tests/test_fixture.h"

#if HS_ENABLE_TEST_HOOKS || HS_ENABLE_TEST_ORACLES ||                          \
    HS_ENABLE_STRUCTURAL_AUDITS || HS_ENABLE_PARAM_GUI_BRIDGE ||               \
    HS_ENABLE_EFFECT_REGISTRY || HS_ENABLE_EFFECT_CONTROL_API
#error "perf_bench must use production render capabilities"
#endif

#define HS_TEMPORAL_BENCH_EFFECTS(X)                                           \
  X(PetalFlow)                                                                 \
  X(MeshFeedback)                                                              \
  X(Comets)                                                                    \
  X(RingSpin)                                                                  \
  X(HopfFibration)                                                             \
  X(Galaxies)

namespace {
constexpr int W = 288, H = 144;
constexpr int FRAMES = 60, WARM = 10;
/** @brief Scratch capacity covering uncached feedback under a column clip. */
constexpr size_t BENCH_SCRATCH_BYTES = 64 * 1024;

/** @brief Frames run before timing an effect. */
template <typename Effect> constexpr int WARM_FRAMES = WARM;
/** @brief Galaxies fills its particle pool over its particle lifetime. */
template <> constexpr int WARM_FRAMES<Galaxies<W, H>> = 900;

/** @brief Measures one effect after a fixed warmup.
 * @param name Effect name, or null to suppress the individual measurement.
 * @param temporal Selects column-time sampling.
 * @param clipped Selects the segment clip.
 * @return Microseconds per frame. */
template <typename Effect>
double bench(const char *name, bool temporal = true, bool clipped = false) {
  hs_test::reset_globals();
  configure_arenas(GLOBAL_ARENA_SIZE - 2 * BENCH_SCRATCH_BYTES,
                   BENCH_SCRATCH_BYTES, BENCH_SCRATCH_BYTES);
  Spacetime::set_flat_surface(!temporal);
  hs_test::pin_frame_clock(0);
  Effect effect;
  effect.init();
  if (clipped)
    effect.set_clip(36, 72, 0, 144);
  for (int i = 0; i < WARM_FRAMES<Effect>; ++i) {
    hs_test::pin_frame_clock(i);
    effect.draw_frame();
    effect.advance_display();
  }
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < FRAMES; ++i) {
    hs_test::pin_frame_clock(WARM_FRAMES<Effect> + i);
    effect.draw_frame();
    effect.advance_display();
  }
  auto t1 = std::chrono::steady_clock::now();
  double us =
      std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() /
      1000.0 / FRAMES;
  volatile auto px = effect.get_pixel(0, 0); // keep the render live
  (void)px;
  if (name)
    std::printf("  %-22s %9.1f us/frame\n", name, us);
  return us;
}
/** @brief Reports interleaved temporal measurements for one effect.
 * @param name Effect name. */
template <typename Effect> void temporal_cost(const char *name) {
  constexpr int REPEATS = 7;
  for (bool clipped : {false, true}) {
    std::array<double, REPEATS> on{}, off{};
    for (int repetition = 0; repetition < REPEATS; ++repetition) {
      if ((repetition & 1) == 0) {
        on[repetition] = bench<Effect>(nullptr, true, clipped);
        off[repetition] = bench<Effect>(nullptr, false, clipped);
      } else {
        off[repetition] = bench<Effect>(nullptr, false, clipped);
        on[repetition] = bench<Effect>(nullptr, true, clipped);
      }
    }
    std::sort(on.begin(), on.end());
    std::sort(off.begin(), off.end());
    std::printf("%s,%s,%.3f,%.3f,%.3f,%.3f\n", name,
                clipped ? "segment" : "full", off.front(), off[REPEATS / 2],
                on.front(), on[REPEATS / 2]);
    std::fflush(stdout);
  }
}
} // namespace

namespace hs_test::effects_tests {
/** @brief Times Galaxies' particle and core draw alone, on a filled pool. */
struct GalaxiesWhiteBox {
  /** @brief Microseconds per draw_particles() after the warmup.
   * @param temporal Selects column-time sampling.
   * @param clipped Selects the segment clip.
   * @return Fastest of several timed blocks, per draw. */
  static double draw_us(bool temporal, bool clipped) {
    constexpr int BLOCKS = 9, DRAWS = 20;
    hs_test::reset_globals();
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * BENCH_SCRATCH_BYTES,
                     BENCH_SCRATCH_BYTES, BENCH_SCRATCH_BYTES);
    Spacetime::set_flat_surface(!temporal);
    Galaxies<W, H> effect;
    effect.init();
    if (clipped)
      effect.set_clip(36, 72, 0, 144);
    const int warm = WARM_FRAMES<Galaxies<W, H>>;
    for (int i = 0; i < warm; ++i) {
      hs_test::pin_frame_clock(i);
      effect.draw_frame();
      effect.advance_display();
    }
    double best = 1e30;
    for (int block = 0; block < BLOCKS; ++block) {
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < DRAWS; ++i) {
        {
          Canvas canvas(effect);
          effect.draw_particles(canvas);
        }
        effect.advance_display();
      }
      const auto t1 = std::chrono::steady_clock::now();
      best = std::min(
          best, std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)
                        .count() /
                    1000.0 / DRAWS);
    }
    return best;
  }
};
} // namespace hs_test::effects_tests

int main(int argc, char **argv) {
  const auto usage = [] {
#define HS_TEMPORAL_NAME(Effect) " " #Effect
    std::fprintf(
        stderr,
        "usage: perf_bench\n"
        "       perf_bench --galaxies-draw on|off full|segment\n"
        "       perf_bench --temporal-sample EFFECT on|off full|segment\n"
        "       perf_bench --temporal-cost\n"
        "effects:" HS_TEMPORAL_BENCH_EFFECTS(HS_TEMPORAL_NAME) "\n");
#undef HS_TEMPORAL_NAME
    return 2;
  };
  const auto valid_modes = [](const char *temporal, const char *clip) {
    return (std::strcmp(temporal, "on") == 0 ||
            std::strcmp(temporal, "off") == 0) &&
           (std::strcmp(clip, "full") == 0 ||
            std::strcmp(clip, "segment") == 0);
  };
  if (argc > 1 && std::strcmp(argv[1], "--galaxies-draw") == 0) {
    if (argc != 4 || !valid_modes(argv[2], argv[3]))
      return usage();
    std::printf("%.3f\n", hs_test::effects_tests::GalaxiesWhiteBox::draw_us(
                              std::strcmp(argv[2], "on") == 0,
                              std::strcmp(argv[3], "segment") == 0));
    return 0;
  }
  if (argc > 1 && std::strcmp(argv[1], "--temporal-sample") == 0) {
    if (argc != 5 || !valid_modes(argv[3], argv[4]))
      return usage();
    const bool temporal = std::strcmp(argv[3], "on") == 0;
    const bool clipped = std::strcmp(argv[4], "segment") == 0;
#define HS_TEMPORAL_SAMPLE(Effect)                                             \
  if (std::strcmp(argv[2], #Effect) == 0) {                                    \
    std::printf("%.3f\n", bench<Effect<W, H>>(nullptr, temporal, clipped));    \
    return 0;                                                                  \
  }
    HS_TEMPORAL_BENCH_EFFECTS(HS_TEMPORAL_SAMPLE)
#undef HS_TEMPORAL_SAMPLE
    std::fprintf(stderr, "unknown temporal effect '%s'\n", argv[2]);
    return usage();
  }
  if (argc > 1 && std::strcmp(argv[1], "--temporal-cost") == 0) {
    if (argc != 2)
      return usage();
    std::puts("effect,clip,off_min_us,off_median_us,on_min_us,on_median_us");
#define HS_TEMPORAL_COST(Effect) temporal_cost<Effect<W, H>>(#Effect);
    HS_TEMPORAL_BENCH_EFFECTS(HS_TEMPORAL_COST)
#undef HS_TEMPORAL_COST
    return 0;
  }
  if (argc != 1)
    return usage();
#ifdef __OPTIMIZE_SIZE__
  const char *opt = "-Os";
#else
  const char *opt = "-O3";
#endif
  std::printf("=== perf_bench [%s] us/frame, %dx%d (%d frames) ===\n", opt, W,
              H, FRAMES);
  double total = 0;
#define HS_BENCH_ONE(name) total += bench<name<W, H>>(#name);
  HS_EFFECT_LIST(HS_BENCH_ONE)
#undef HS_BENCH_ONE
  std::printf("  %-22s %9.1f us/frame (sum)\n", "TOTAL", total);
  return 0;
}

#undef HS_TEMPORAL_BENCH_EFFECTS
