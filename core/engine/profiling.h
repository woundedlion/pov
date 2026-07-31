/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

// Profiling and metrics layer: scan/probe counters and the cycle-counting
// instrumentation. Included from the end of engine/platform.h, which supplies
// HS_OS_CYCLES and hs::log.

#include <cstdint>
#include <cstring>

namespace hs {

/** @brief Scanline profiling counters (platform-agnostic). */
struct ScanMetrics {
  uint32_t plot = 0;          /**< Cycles spent plotting pixels. */
  uint32_t sdf_dist = 0;      /**< Cycles spent evaluating SDF distances. */
  uint32_t frag_shader = 0;   /**< Cycles spent in the fragment shader. */
  uint32_t bounds = 0;        /**< Cycles spent computing bounds. */
  uint32_t face_setup = 0;    /**< Cycles spent on per-face setup. */
  uint32_t scan_loop = 0;     /**< Cycles spent in the scanline loop. */
  uint32_t pixels_tested = 0; /**< Count of pixels tested. */
  uint32_t pixels_culled = 0; /**< Count of pixels culled before shading. */
  uint32_t exact_hits =
      0; /**< Count of full distance evaluations (umbrella: convex, sector, and exact-walk paths). */
  uint32_t convex_hits =
      0; /**< Count of evaluations on the convex half-plane path. */
  uint32_t sector_hits =
      0; /**< Count of evaluations on the concave sector walk (subset of exact_hits). */
  uint32_t lut_hits = 0; /**< Count of class-LUT bilinear serves. */
  uint32_t plot_backstop_hits =
      0; /**< Count of plot() steps_cache capacity-backstop trips. */
  uint32_t shade_candidates =
      0; /**< Count of pixels passing the scan's d < pixel_width test (shading + alpha-rejected). */
  /** @brief Zeroes every counter. */
  void reset() {
    plot = sdf_dist = frag_shader = bounds = face_setup = scan_loop =
        pixels_tested = pixels_culled = exact_hits = convex_hits = sector_hits =
            lut_hits = plot_backstop_hits = shade_candidates = 0;
  }
};
/** @brief Global scanline profiling counters. Compiled in only when
 *  HS_SCAN_METRICS is defined; otherwise HS_SCAN_METRIC(...) expands to nothing
 *  and this would be dead storage. */
#ifdef HS_SCAN_METRICS
inline ScanMetrics g_scan_metrics;
#endif

/**
 * @brief Per-probe cycle buckets splitting one Face::distance probe into its
 *        stages, plus the event counts each bucket divides by.
 * @details Accumulated from raw cycle-counter deltas rather than CycleScope
 * RAII, which at tens of thousands of probes per frame would cost more than the
 * stages being measured. One counter read per stage boundary; `tick` sums a
 * back-to-back read pair per probe so a capture measures its own read cost and
 * the buckets can be discounted by it. Report ratios from such a build, not
 * absolute times.
 */
struct ProbeBreakdown {
  uint32_t point = 0; /**< Cycles: probe entry through the back-face cull. */
  uint32_t project =
      0; /**< Cycles: gnomonic projection through the radius cull. */
  uint32_t edge_lut = 0;    /**< Cycles: class-LUT bilinear serve. */
  uint32_t edge_convex = 0; /**< Cycles: convex half-plane max. */
  uint32_t edge_sector = 0; /**< Cycles: concave sector walk. */
  uint32_t edge_exact = 0;  /**< Cycles: full per-edge walk. */
  uint32_t pack =
      0; /**< Cycles: plane->angle conversion and result packaging. */
  uint32_t alpha = 0;   /**< Cycles: scan-side AA coverage kernel. */
  uint32_t tick = 0;    /**< Cycles: summed back-to-back counter-read pairs. */
  uint32_t n_probe = 0; /**< Probes entered. */
  uint32_t n_cull_cos = 0; /**< Probes leaving at the back-face cull. */
  uint32_t n_cull_r = 0;   /**< Probes leaving at the radius cull. */
  uint32_t n_lut = 0;      /**< Probes served by the class LUT. */
  uint32_t n_convex = 0;   /**< Probes taking the convex path. */
  uint32_t n_sector = 0;   /**< Probes taking the sector walk. */
  uint32_t n_exact = 0;    /**< Probes taking the full edge walk. */
  uint32_t n_alpha = 0;    /**< Probes reaching the AA coverage kernel. */
  /** @brief Zeroes every bucket and count. */
  void reset() {
    point = project = edge_lut = edge_convex = edge_sector = edge_exact = pack =
        alpha = tick = 0;
    n_probe = n_cull_cos = n_cull_r = n_lut = n_convex = n_sector = n_exact =
        n_alpha = 0;
  }
};
/** @brief Global per-probe cycle buckets. Compiled in only when
 *  HS_PROBE_BREAKDOWN is defined; otherwise HS_PROBE_* expand to nothing and
 *  this would be dead storage. */
#ifdef HS_PROBE_BREAKDOWN
inline ProbeBreakdown g_probe_breakdown;
#endif

#ifdef HS_PLOT_COUNTS
/** @brief Count-only Plot rasterizer workload attribution. */
struct PlotCounts {
  uint32_t rings = 0;
  uint32_t edges = 0;
  uint32_t planar = 0;
  uint32_t geodesic = 0;
  uint32_t degenerate = 0;
  uint32_t one_dot = 0;
  uint32_t cull_tests = 0;
  uint32_t culled = 0;
  uint32_t sim_samples = 0;
  uint32_t replay_samples = 0;
  uint32_t planar_unprojects = 0;
  uint32_t planar_arc_samples = 0;
  uint32_t normalizations = 0;
  uint32_t shader_calls = 0;
  uint32_t plotted_samples = 0;
  uint32_t steps_peak = 0;
  uint32_t backstops = 0;

  /** @brief Zeroes every count. */
  void reset() { *this = {}; }
};

inline PlotCounts g_plot_counts;
#endif

} // namespace hs

// Per-pixel scan instrumentation is OFF by default: a g_scan_metrics increment is
// a non-atomic global load-modify-store on a shared cache line for every pixel.
// Define HS_SCAN_METRICS to compile the counters in (the native test build does,
// to assert which Face::distance path each sample took); otherwise
// HS_SCAN_METRIC(...) expands to nothing.
#ifdef HS_SCAN_METRICS
#define HS_SCAN_METRIC(stmt)                                                   \
  do {                                                                         \
    (stmt);                                                                    \
  } while (0)
#else
#define HS_SCAN_METRIC(stmt) ((void)0)
#endif

// Per-probe stage timing, OFF by default: each boundary is a cycle-counter read
// plus a global accumulate, which at the scan's probe rate distorts the very
// stages it splits. Define HS_PROBE_BREAKDOWN to compile the buckets in and read
// RATIOS from the capture, discounted by the self-measured `tick` read cost.
// HS_PROBE_MARK opens a rolling timestamp; HS_PROBE_SPAN closes one stage and
// reopens the next off the same read, so a chain of N stages costs N reads.
#ifdef HS_PROBE_BREAKDOWN
#define HS_PROBE_MARK(var) uint32_t var = HS_OS_CYCLES()
#define HS_PROBE_SPAN(field, var)                                              \
  do {                                                                         \
    uint32_t hs_now = HS_OS_CYCLES();                                          \
    hs::g_probe_breakdown.field += hs_now - (var);                             \
    (var) = hs_now;                                                            \
  } while (0)
#define HS_PROBE_COUNT(field)                                                  \
  do {                                                                         \
    ++hs::g_probe_breakdown.field;                                             \
  } while (0)
#define HS_PROBE_TICK()                                                        \
  do {                                                                         \
    uint32_t hs_a = HS_OS_CYCLES();                                            \
    uint32_t hs_b = HS_OS_CYCLES();                                            \
    hs::g_probe_breakdown.tick += hs_b - hs_a;                                 \
  } while (0)
#else
#define HS_PROBE_MARK(var)
#define HS_PROBE_SPAN(field, var) ((void)0)
#define HS_PROBE_COUNT(field) ((void)0)
#define HS_PROBE_TICK() ((void)0)
#endif

#ifdef HS_PLOT_COUNTS
#define HS_PLOT_COUNT(field) (++hs::g_plot_counts.field)
#define HS_PLOT_ADD(field, value) (hs::g_plot_counts.field += (value))
#define HS_PLOT_MAX(field, value)                                             \
  do {                                                                         \
    uint32_t hs_plot_value = static_cast<uint32_t>(value);                     \
    if (hs_plot_value > hs::g_plot_counts.field)                               \
      hs::g_plot_counts.field = hs_plot_value;                                 \
  } while (0)
#else
#define HS_PLOT_COUNT(field) ((void)0)
#define HS_PLOT_ADD(field, value) ((void)0)
#define HS_PLOT_MAX(field, value) ((void)0)
#endif

// ---------------------------------------------------------------------------
// Cycle-counting instrumentation
//   CycleCounter — named cumulative accumulator (self-registers for bulk log)
//   CycleScope   — RAII guard that accumulates into a CycleCounter
//   HS_PROFILE   — one-liner convenience macro
// ---------------------------------------------------------------------------
namespace hs {

/**
 * @brief Formats v as decimal into buf.
 * @param v Value to format.
 * @param buf Buffer of at least 21 bytes (20 digits + NUL).
 * @return Pointer to the first digit inside buf.
 * @details Manual conversion because newlib-nano's integer printf (the -Os
 *          device build) has no long-long support.
 */
inline const char *u64_dec(uint64_t v, char *buf) {
  char *p = buf + 20;
  *p = '\0';
  do {
    *--p = static_cast<char>('0' + v % 10);
    v /= 10;
  } while (v);
  return p;
}

/**
 * @brief Named cumulative cycle accumulator. Each instance self-registers into
 *        a static intrusive list at construction so log_all()/reset_all() can
 *        walk every counter without a central registry. Counters nest: a
 *        CycleScope sets `parent` to whichever counter was active when this one
 *        started, giving log_all() a call tree with per-parent percentages.
 * @warning REENTRANCY: the registry head and the `active` nesting pointer are
 *        non-atomic statics (like hs::random()'s generator), so construction and
 *        CycleScope enter/exit are main-loop-only — driving a CycleScope from an
 *        ISR would race the list/active pointer and corrupt the call tree.
 */
struct CycleCounter {
  static constexpr uint32_t CYCLES_PER_US =
      600; /**< Core clock: Teensy 4 @ 600 MHz. */

  const char *name;    /**< Counter label used in log output. */
  uint64_t cycles = 0; /**< Accumulated cycle count. 64-bit because a
                                        32-bit accumulator overflows after only
                                        ~7 s of summed time at 600 MHz, which a
                                        multi-frame profiling run easily exceeds. */
  uint32_t count = 0;  /**< Number of timed invocations. */
  CycleCounter *parent = nullptr; /**< Enclosing counter for tree nesting. */
  CycleCounter *next =
      nullptr; /**< Next link in the intrusive registry list. */

  /**
   * @brief Constructs a named counter and self-registers it for bulk logging.
   * @param n Counter label (must outlive the counter; typically a literal).
   */
  explicit CycleCounter(const char *n) : name(n), next(head) { head = this; }

  /** @brief Zeroes this counter's accumulated cycles and call count. */
  void reset() {
    cycles = 0;
    count = 0;
  }

  /** @brief Logs every root counter (no parent) and its subtree as a tree. */
  static void log_all() {
    hs::log("--- Cycle Counters ---");
    for (auto *c = head; c; c = c->next)
      if (!c->parent && c->count)
        log_node(c, 0);
  }

  /** @brief Zeroes every registered counter (between profiling runs). */
  static void reset_all() {
    for (auto *c = head; c; c = c->next)
      c->reset();
  }

  /**
   * @brief Finds the first registered counter whose name ends with @p suffix.
   * @param suffix Name suffix to match (e.g. "_buffer_wait").
   * @return The counter, or nullptr if none is registered yet (counters
   *         self-register on first scope entry).
   */
  static CycleCounter *find_suffix(const char *suffix) {
    const size_t sl = strlen(suffix);
    for (auto *c = head; c; c = c->next) {
      const size_t nl = strlen(c->name);
      if (nl >= sl && memcmp(c->name + nl - sl, suffix, sl) == 0)
        return c;
    }
    return nullptr;
  }

private:
  static inline CycleCounter *head =
      nullptr; /**< Head of the intrusive registry list. */
  static inline CycleCounter *active =
      nullptr; /**< Currently active counter (for nesting). */
  friend struct CycleScope;

  /**
   * @brief Recursively logs one counter node and its children as a tree.
   * @param node Counter node to log.
   * @param depth Tree depth; drives indentation.
   * @details The reported percentage is this node's cycles over its parent's
   *          (or 100% for a root), and cycles are converted to microseconds via
   *          CYCLES_PER_US.
   */
  static void log_node(const CycleCounter *node, int depth) {
    if (!node->count)
      return;
    uint64_t ref = node->parent ? node->parent->cycles : node->cycles;
    uint32_t pct = ref ? (uint32_t)(node->cycles * 100 / ref) : 100;
    char cyc_buf[21], us_buf[21];
    const char *cyc = hs::u64_dec(node->cycles, cyc_buf);
    const char *us = hs::u64_dec(node->cycles / CYCLES_PER_US, us_buf);
    int indent = depth * 2;
    int name_w = 22 - indent;
    if (name_w < 1)
      name_w = 1;
    hs::log("%*s%-*s %s us (%lu%%)  %lu calls  %s cyc", indent, "", name_w,
            node->name, us, (unsigned long)pct, (unsigned long)node->count,
            cyc);
    for (auto *c = head; c; c = c->next)
      if (c->parent == node)
        log_node(c, depth + 1);
  }
};

/**
 * @brief RAII guard that times its enclosing scope and accumulates the elapsed
 *        cycles into a CycleCounter. On construction it makes its counter the
 *        active one (recording the previously-active counter as parent on first
 *        use) and snapshots the cycle counter; the destructor adds the delta and
 *        restores the previous active counter, rebuilding the nesting tree.
 */
struct CycleScope {
  CycleCounter &counter; /**< Counter this scope accumulates into. */
  CycleCounter
      *prev_active; /**< Counter to restore as active on destruction. */
  uint32_t start;   /**< Cycle snapshot taken at construction (32-bit,
                                    matching the hardware DWT CYCCNT register). */

  /**
   * @brief Begins timing the enclosing scope into the given counter.
   * @param c Counter that receives the elapsed cycles.
   * @details Makes c the active counter (recording the previously-active
   *          counter as its parent on first use) and snapshots the cycle
   *          counter.
   */
  explicit CycleScope(CycleCounter &c) : counter(c), start(HS_OS_CYCLES()) {
    prev_active = CycleCounter::active;
    if (!counter.parent && prev_active)
      counter.parent = prev_active;
    CycleCounter::active = &counter;
  }
  /**
   * @brief Adds the elapsed cycles to the counter and restores the previous one.
   * @pre The scope must not span a full CYCCNT wrap (~7 s at 600 MHz). The delta
   *      below is a 32-bit subtraction (matching the hardware register width),
   *      correct modulo 2^32, so a single scope longer than one wrap reads short
   *      by a multiple of 2^32. Accumulation across scopes is wrap-safe: the
   *      32-bit delta widens into the 64-bit `cycles` accumulator.
   */
  ~CycleScope() {
    counter.cycles += (uint32_t)(HS_OS_CYCLES() - start);
    counter.count++;
    CycleCounter::active = prev_active;
  }

  /**
   * @brief Deleted copy constructor; a scope guard must not be copied.
   */
  CycleScope(const CycleScope &) = delete;
  /**
   * @brief Deleted copy assignment; a scope guard must not be copied.
   * @return Never returns; deleted.
   */
  CycleScope &operator=(const CycleScope &) = delete;
};

/**
 * @brief ISR-safe cycle accumulator: plain single-writer fields, no registry.
 * @details CycleCounter/CycleScope are main-loop-only (non-atomic registry and
 *          nesting pointer), so ISR paths accumulate into one of these
 *          instead. Contract: the ISR is the sole writer; a foreground reader
 *          copies and reset()s under a brief IRQ-off window.
 */
struct IsrCycleStats {
  uint64_t cycles = 0;       /**< Accumulated cycles across all scopes. */
  uint32_t count = 0;        /**< Number of timed scopes. */
  uint32_t min = UINT32_MAX; /**< Shortest single scope, in cycles. */
  uint32_t max = 0;          /**< Longest single scope, in cycles. */

  /** @brief Folds one scope's elapsed cycles into the accumulator. */
  void add(uint32_t dt) {
    cycles += dt;
    ++count;
    if (dt < min)
      min = dt;
    if (dt > max)
      max = dt;
  }
  /** @brief Zeroes the accumulator (foreground, IRQs off). */
  void reset() {
    cycles = 0;
    count = 0;
    min = UINT32_MAX;
    max = 0;
  }
};

/**
 * @brief RAII guard timing its enclosing scope into an IsrCycleStats.
 */
struct IsrCycleScope {
  IsrCycleStats &stats; /**< Accumulator receiving the elapsed cycles. */
  uint32_t start;       /**< Cycle snapshot taken at construction. */

  explicit IsrCycleScope(IsrCycleStats &s) : stats(s), start(HS_OS_CYCLES()) {}
  ~IsrCycleScope() { stats.add((uint32_t)(HS_OS_CYCLES() - start)); }
  IsrCycleScope(const IsrCycleScope &) = delete;
  IsrCycleScope &operator=(const IsrCycleScope &) = delete;
};

} // namespace hs

/**
 * @brief Times the enclosing scope into an IsrCycleStats instance.
 * @param stats An hs::IsrCycleStats lvalue expression. One use per block (the
 *        guard has a fixed name; open a nested block for a second scope).
 * @details The ISR-context sibling of HS_PROFILE; compiled in only under
 *          HS_PROFILE_ENABLE.
 */
#ifdef HS_PROFILE_ENABLE
#define HS_ISR_PROFILE(stats) hs::IsrCycleScope hs_isr_scope(stats)
#else
#define HS_ISR_PROFILE(stats) ((void)0)
#endif

/**
 * @brief Times the enclosing scope into a named cycle counter.
 * @param label Counter name (used both as the identifier suffix and log label).
 * @details Compiled in only under HS_PROFILE_ENABLE; off by default so regular
 *          builds pay nothing for the per-scope bookkeeping and CYCCNT read on
 *          every hot-path face/pixel scope.
 */
#ifdef HS_PROFILE_ENABLE
#define HS_PROFILE(label)                                                      \
  static hs::CycleCounter hs_ctr_##label(#label);                              \
  hs::CycleScope hs_scope_##label(hs_ctr_##label)
#else
#define HS_PROFILE(label) ((void)0)
#endif

#if defined(HS_PROFILE_DEEP_ENABLE) && !defined(HS_PROFILE_ENABLE)
#error "HS_PROFILE_DEEP_ENABLE needs HS_PROFILE_ENABLE (the counter registry)"
#endif

/**
 * @brief Times the enclosing scope, but only in a deep-profile build.
 * @param label Counter name (used both as the identifier suffix and log label).
 * @details The form shared render code must use for any scope entered more than
 *          once per draw call — per pixel, per cell, per face. Such a scope sits
 *          on every effect's hot path, so leaving it in the ordinary `profile`
 *          image would tax the whole roster's numbers to instrument one effect
 *          (the feedback composite's per-pixel set costs ~8 cyc per entry,
 *          ~2.5% of the flush). Off unless HS_PROFILE_DEEP_ENABLE is defined on
 *          top of HS_PROFILE_ENABLE, so a deep run is opt-in per capture
 *          (HS_PROFILE_DEEP=1 in profile_one.sh, deep=1 in `just profile`).
 *          Per-frame scopes stay on plain HS_PROFILE — they are what the
 *          standard reports are built from.
 */
#ifdef HS_PROFILE_DEEP_ENABLE
#define HS_PROFILE_DEEP(label) HS_PROFILE(label)
#else
#define HS_PROFILE_DEEP(label) ((void)0)
#endif
