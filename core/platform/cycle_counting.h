/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/platform/profiling.h.

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
 * @param buf Exact-fit buffer: 20 digits plus the NUL.
 * @return Pointer to the first digit inside buf.
 * @details Manual conversion because newlib-nano's integer printf (the -Os
 *          device build) has no long-long support.
 */
inline const char *u64_dec(uint64_t v, char (&buf)[21]) {
  char *p = buf + sizeof(buf) - 1;
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
 *        CycleScope latches `parent` to whichever counter was active the first
 *        time this one started, giving log_all() a call tree with per-parent
 *        percentages. A counter entered under two different callers keeps the
 *        latched parent and sets `mixed_parent`, which log_all() marks in the
 *        report — the tree is a single-parent view, so the marked node's share
 *        of that parent counts cycles the other caller spent.
 * @warning REENTRANCY: the registry head and the `active` nesting pointer are
 *        non-atomic statics (like hs::random()'s generator), so construction and
 *        CycleScope enter/exit are main-loop-only — driving a CycleScope from an
 *        ISR would race the list/active pointer and corrupt the call tree.
 */
struct CycleCounter {
#ifdef CORE_TEENSY
  static constexpr uint32_t CYCLES_PER_US = F_CPU / 1000000;
#else
  static constexpr uint32_t CYCLES_PER_US = 600;
#endif

  const char *name;               /**< Counter label used in log output. */
  uint64_t cycles = 0;            /**< Accumulated cycle count. 64-bit because a
                                         32-bit accumulator overflows after only
                                         ~7 s of summed time at 600 MHz, which a
                                         multi-frame profiling run easily exceeds. */
  uint32_t count = 0;             /**< Number of timed invocations. */
  CycleCounter *parent = nullptr; /**< Enclosing counter for tree nesting. */
  CycleCounter *next =
      nullptr; /**< Next link in the intrusive registry list. */

  /** @brief True once `parent` has been latched, so a null `parent` afterwards
   *         means a genuine root rather than "not nested yet". */
  bool parented = false;
  /** @brief True once entered under a counter other than the latched `parent`.
   *         Every caller's cycles land on `parent`, so its share is overstated
   *         and can exceed 100%. */
  bool mixed_parent = false;

  /**
   * @brief Constructs a named counter and self-registers it for bulk logging.
   * @param n Counter label (must outlive the counter; typically a literal).
   */
  explicit CycleCounter(const char *n) : name(n), next(head) { head = this; }

  /**
   * @brief Unregisters the counter, so no registry walk reaches dead storage.
   * @details Also unlatches every counter that had latched this one as its
   * parent; log_all() and log_node() dereference `parent`, so a surviving edge
   * into this storage is read after the destructor returns. An unlatched
   * counter re-latches on its next entry, unless that latch would close a
   * parent cycle.
   */
  ~CycleCounter() {
    if (active == this)
      active = nullptr;
    for (CycleCounter **p = &head; *p;) {
      if (*p == this) {
        *p = next;
        continue;
      }
      if ((*p)->parent == this) {
        (*p)->parent = nullptr;
        (*p)->parented = false;
      }
      p = &(*p)->next;
    }
  }

  // The registry links every counter by address, so a counter is fixed in place.
  CycleCounter(const CycleCounter &) = delete;
  CycleCounter &operator=(const CycleCounter &) = delete;

  /**
   * @brief Zeroes this counter's accumulated cycles, call count and
   *        mixed-parent flag.
   * @details mixed_parent describes the entries just discarded, so it clears
   * with them; the next run re-sets it if that run mixes callers. The latched
   * parent edge survives (see log_all()).
   */
  void reset() {
    cycles = 0;
    count = 0;
    mixed_parent = false;
  }

  /**
   * @brief Logs every root counter and its subtree as a tree.
   * @details reset() zeroes counts but keeps the latched parent edges, so a
   *          counter whose parent saw no entries this run is logged as a root
   *          rather than dropped with the parent it is no longer reached from.
   */
  static void log_all() {
    hs::log("--- Cycle Counters ---");
    for (auto *c = head; c; c = c->next)
      if (c->count && (!c->parent || !c->parent->count))
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
   *         self-register at construction).
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
   * @brief Whether @p node is reachable by following @p from's parent chain.
   * @param from Counter to start the walk at; may be null.
   * @param node Counter looked for.
   * @return True when the chain reaches @p node.
   * @details ~CycleCounter unlatches its children, so a re-latch is no longer
   *          time-ordered and a later entry can otherwise close a parent cycle.
   */
  static bool in_parent_chain(const CycleCounter *from,
                              const CycleCounter *node) {
    for (const CycleCounter *c = from; c; c = c->parent)
      if (c == node)
        return true;
    return false;
  }

  /**
   * @brief Reports whether another counter with entries shares @p node's name.
   * @param node Counter being logged.
   * @return True when a second registered counter carries the same label.
   * @details A HS_PROFILE scope inside a function template registers one counter
   *          per instantiation, all under the macro's label, so the report would
   *          otherwise show several identical rows each covering a fraction of
   *          the label's cycles.
   */
  static bool duplicate_name(const CycleCounter *node) {
    for (auto *c = head; c; c = c->next)
      if (c != node && c->count && strcmp(c->name, node->name) == 0)
        return true;
    return false;
  }

  /**
   * @brief Recursively logs one counter node and its children as a tree.
   * @param node Counter node to log.
   * @param depth Tree depth; drives indentation.
   * @details The reported percentage is this node's cycles over its parent's
   *          (or 100% for a root), and cycles are converted to microseconds via
   *          CYCLES_PER_US. A mixed_parent node carries a MIXED-PARENT tag: its
   *          cycles include entries made from callers other than the parent it
   *          is printed under. A node sharing its name with another registered
   *          counter carries a DUPLICATE-NAME tag: its row accounts for only one
   *          of them.
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
    hs::log("%*s%-*s %s us (%lu%%)  %lu calls  %s cyc%s%s", indent, "", name_w,
            node->name, us, (unsigned long)pct, (unsigned long)node->count, cyc,
            node->mixed_parent ? "  MIXED-PARENT" : "",
            duplicate_name(node) ? "  DUPLICATE-NAME" : "");
    for (auto *c = head; c; c = c->next)
      if (c->parent == node)
        log_node(c, depth + 1);
  }
};

/**
 * @brief RAII guard that times its enclosing scope and accumulates the elapsed
 *        cycles into a CycleCounter. On construction it makes its counter the
 *        active one (latching the previously-active counter as parent on first
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
   * @details Makes c the active counter (latching the previously-active
   *          counter as its parent on first use) and snapshots the cycle
   *          counter. A later entry under a different counter flags
   *          mixed_parent instead of re-parenting: the tree is a single-parent
   *          view, and log_all() marks the node so the inflated share is
   *          visible rather than silent.
   */
  explicit CycleScope(CycleCounter &c) : counter(c), start(HS_OS_CYCLES()) {
    prev_active = CycleCounter::active;
    // A recursive scope would self-parent, hiding the counter from log_all()'s
    // root walk; a latch that closes a parent cycle would make log_node()'s
    // recursion diverge. The counter stays unlatched in either case and gets
    // another chance on its next entry.
    if (prev_active != &counter) {
      if (!counter.parented) {
        if (!CycleCounter::in_parent_chain(prev_active, &counter)) {
          counter.parented = true;
          counter.parent = prev_active;
        }
      } else if (counter.parent != prev_active) {
        counter.mixed_parent = true;
      }
    }
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
 *          every hot-path face/pixel scope. The enabled expansion is a guard
 *          declaration, so it must open a braced scope: as the unbraced body of
 *          an `if` or a loop it compiles, is destroyed on the same line, and
 *          records ~0 cycles. HS_OS_CYCLES() is 0 on every non-Teensy target,
 *          so host and WASM captures read all-zero.
 */
#ifdef HS_PROFILE_ENABLE
#define HS_PROFILE(label)                                                      \
  hs::CycleScope hs_scope_##label([]() -> hs::CycleCounter & {                 \
    static hs::CycleCounter ctr(#label);                                       \
    return ctr;                                                                \
  }())
#else
#define HS_PROFILE(label) ((void)0)
#endif

#if defined(HS_PROFILE_DEEP_ENABLE) && !defined(HS_PROFILE_ENABLE)
#error "HS_PROFILE_DEEP_ENABLE needs HS_PROFILE_ENABLE (the counter registry)"
#endif

/**
 * @brief Times the enclosing scope, but only in a deep-profile build.
 * @param label Counter name (used both as the identifier suffix and log label).
 * @details Shared per-pixel, per-cell and per-face instrumentation uses this
 * form unless the standard report consumes its counter. Report counters such
 * as filter_blend, scan_face_setup, scan_mesh_raster and plot_ps_* stay on
 * plain HS_PROFILE, including their per-sample instrumentation cost.
 * Deep scopes require HS_PROFILE_DEEP_ENABLE on top of HS_PROFILE_ENABLE:
 * HS_PROFILE_DEEP=1 in profile_one.sh, or the third positional argument of
 * `just profile`, e.g. `just profile MeshFeedback 150 1`.
 * Per-frame counters also stay on plain HS_PROFILE.
 */
#ifdef HS_PROFILE_DEEP_ENABLE
#define HS_PROFILE_DEEP(label) HS_PROFILE(label)
#else
#define HS_PROFILE_DEEP(label) ((void)0)
#endif
