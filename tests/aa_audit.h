/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host-only AA coverage audit for the mesh scan (HS_AA_AUDIT).
 *
 * Brute-forces every column of every scanned row and records the pixels whose
 * distance clears the shade threshold but that the emitted column runs never
 * visit. Zero code when HS_AA_AUDIT is undefined.
 */
#pragma once

#include <cmath>
#include <cstdint>

namespace hs_aa {

inline constexpr int MAX_ROWS = 512;

struct Audit {
  long long probes = 0;      /**< Pixels the emitted runs evaluate. */
  long long painted = 0;     /**< Shaded pixels inside the runs. */
  long long missed = 0;      /**< Shaded pixels outside the runs. */
  long long probe_rows[MAX_ROWS] = {};
  long long missed_rows[MAX_ROWS] = {};
  long long painted_rows[MAX_ROWS] = {};
  double missed_alpha_sum = 0.0;
  double missed_alpha_max = 0.0;
  long long alpha_hist[10] = {}; /**< Missed-pixel alpha, 10 deciles. */
  int max_gap_cols = 0;          /**< Widest column gap to a run edge. */
  int frames = 0;
  bool enabled = false;   /**< Run the brute-force coverage comparison. */
  bool full_scan = false; /**< Force every face to the full-width scan. */
  bool legacy_pad = false; /**< Use the constant AA pad (pre-fix behaviour). */

  void reset() {
    probes = painted = missed = 0;
    missed_alpha_sum = missed_alpha_max = 0.0;
    max_gap_cols = 0;
    frames = 0;
    for (int i = 0; i < MAX_ROWS; ++i)
      missed_rows[i] = painted_rows[i] = probe_rows[i] = 0;
    for (int i = 0; i < 10; ++i)
      alpha_hist[i] = 0;
  }

  void note_missed(int y, float alpha, int gap) {
    ++missed;
    if (y >= 0 && y < MAX_ROWS)
      ++missed_rows[y];
    missed_alpha_sum += alpha;
    if (alpha > missed_alpha_max)
      missed_alpha_max = alpha;
    int b = static_cast<int>(alpha * 10.0f);
    if (b < 0)
      b = 0;
    if (b > 9)
      b = 9;
    ++alpha_hist[b];
    if (gap > max_gap_cols)
      max_gap_cols = gap;
  }

  void note_painted(int y) {
    ++painted;
    if (y >= 0 && y < MAX_ROWS)
      ++painted_rows[y];
  }
};

inline Audit g_audit;

} // namespace hs_aa
