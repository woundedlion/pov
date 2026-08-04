/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Emits the segment->canvas mapping hardware/pov_segment_map.h derives, as JSON
 * on stdout: for every (S, N) config the simulator's cross-check sweeps, each
 * segment's arm side, LED-0 row, strip direction and full row traversal, plus
 * the arm-A/arm-B sampled column per canvas width.
 *
 * The committed hardware/pov_segment_map.json is diffed against this program's
 * output by the unit_pov_segment_map_golden CTest, so a convention change in
 * the header fails the suite until the golden is regenerated. The WASM install
 * copies the golden into the daydream checkout, where its
 * tests/segment_crosscheck.test.js reads it as the firmware reference rather
 * than re-porting the header in JS.
 *
 * Regenerate with:
 *   cmake --build --preset tests --target pov_segment_map_gen
 *   build/tests/tests/pov_segment_map_gen > hardware/pov_segment_map.json
 *
 * Every printf here stays on one source line under the column limit so the
 * emitted layout is readable next to the JSON it produces.
 */
#include <cstdio>
#include <iterator>

#include "hardware/pov_segment_map.h"

namespace {

struct Config {
  int leds;     /**< S: total LEDs across both arms. */
  int segments; /**< N: segment count (power of two <= 8). */
};

// Every config tests/segment_crosscheck.test.js sweeps, plus the two the host
// tests in tests/test_pov_segmented.h pin directly.
constexpr Config CONFIGS[] = {
    {288, 2}, {288, 4}, {288, 8}, {8, 4}, {8, 8},
};

// Canvas widths the cross-check pairs with those configs.
constexpr int WIDTHS[] = {8, 96, 288};

constexpr int CONFIG_COUNT = static_cast<int>(std::size(CONFIGS));
constexpr int WIDTH_COUNT = static_cast<int>(std::size(WIDTHS));

/**
 * @brief Emits one segment's mapping as a JSON object on a single line.
 * @param m Segment mapping from pov::segment_map.
 * @param pps LEDs per segment (S / N).
 * @param tail Separator after the object: "," unless this is the last segment.
 */
void emit_segment(const pov::SegmentMap &m, int pps, const char *tail) {
  std::printf("        { \"arm_b\": %s,", m.arm_b ? "true" : "false");
  std::printf(" \"y_base\": %d, \"y_step\": %d,", m.y_base, m.y_step);
  std::printf(" \"y\": [");
  for (int i = 0; i < pps; i++)
    std::printf("%s%d", i ? ", " : "", pov::segment_y(m, i));
  std::printf("] }%s\n", tail);
}

/**
 * @brief Emits one (S, N) config's segment table.
 * @param cfg LED and segment counts.
 * @param tail Separator after the object: "," unless this is the last config.
 */
void emit_config(const Config &cfg, const char *tail) {
  const int S = cfg.leds;
  const int N = cfg.segments;
  // ROWS and PPS as pov_segment_map.h defines them internally.
  const int rows = S / 2;
  const int pps = S / N;
  std::printf("    {\n");
  std::printf("      \"leds\": %d,\n", S);
  std::printf("      \"segment_count\": %d,\n", N);
  std::printf("      \"rows\": %d,\n", rows);
  std::printf("      \"pps\": %d,\n", pps);
  std::printf("      \"segments\": [\n");
  for (int id = 0; id < N; id++)
    emit_segment(pov::segment_map(id, S, N), pps, id == N - 1 ? "" : ",");
  std::printf("      ]\n");
  std::printf("    }%s\n", tail);
}

/**
 * @brief Emits the columns arm A and arm B sample at rotation column 0.
 * @param w Canvas width.
 * @param tail Separator after the object: "," unless this is the last width.
 */
void emit_x_cols(int w, const char *tail) {
  const int arm_a = pov::segment_x_col(false, 0, w);
  const int arm_b = pov::segment_x_col(true, 0, w);
  std::printf("    { \"width\": %d,", w);
  std::printf(" \"arm_a\": %d, \"arm_b\": %d }%s\n", arm_a, arm_b, tail);
}

} // namespace

int main() {
  std::printf("{\n");
  std::printf("  \"source\": \"hardware/pov_segment_map.h\",\n");
  std::printf("  \"generator\": \"tools/pov_segment_map_export.cpp\",\n");
  std::printf("  \"configs\": [\n");
  for (int c = 0; c < CONFIG_COUNT; c++)
    emit_config(CONFIGS[c], c == CONFIG_COUNT - 1 ? "" : ",");
  std::printf("  ],\n");
  std::printf("  \"x_cols\": [\n");
  for (int i = 0; i < WIDTH_COUNT; i++)
    emit_x_cols(WIDTHS[i], i == WIDTH_COUNT - 1 ? "" : ",");
  std::printf("  ]\n");
  std::printf("}\n");
  return 0;
}
