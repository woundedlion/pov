/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/platform/profiling.h.

// ---------------------------------------------------------------------------
// Cycle-counting instrumentation
// ---------------------------------------------------------------------------
namespace hs {

/**
 * @brief Formats v as decimal into buf.
 * @param v Value to format.
 * @param buf Exact-fit buffer: 20 digits plus the NUL.
 * @return Pointer to the first digit inside buf.
 * @details newlib-nano's integer printf has no long-long support.
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
 * @brief Named cumulative cycle accumulator, self-registered in a static
 *        intrusive list at construction.
 * @details `parent` latches to the counter active on first entry; a later
 *          entry under a different counter sets `mixed_parent` instead.
 * @warning The registry head and `active` pointer are non-atomic statics:
 *          construction and CycleScope enter/exit are main-loop-only.
 */
struct CycleCounter {
#ifdef CORE_TEENSY
  static constexpr uint32_t CYCLES_PER_US = F_CPU / 1000000;
#else
  static constexpr uint32_t CYCLES_PER_US = 600;
#endif

  const char *name;               /**< Counter label used in log output. */
  uint64_t cycles = 0;            /**< Accumulated cycle count (32 bits would
                                         wrap after ~7 s at 600 MHz). */
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
   * @details Also unlatches every counter that latched this one as its parent;
   * an unlatched counter re-latches on its next entry.
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
   * @details The latched parent edge survives.
   */
  void reset() {
    cycles = 0;
    count = 0;
    mixed_parent = false;
  }

  /**
   * @brief Logs every root counter and its subtree as a tree.
   * @details A counter whose parent has no entries this run logs as a root.
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
   * @details Guards a re-latch against closing a parent cycle.
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
   *          per instantiation under the same label.
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
   * @details The percentage is relative to the parent's cycles (100% for a
   *          root). MIXED-PARENT marks cycles entered from other callers;
   *          DUPLICATE-NAME marks a row covering one of several same-named
   *          counters.
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
 *        cycles into a CycleCounter, maintaining the nesting tree.
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
   */
  explicit CycleScope(CycleCounter &c) : counter(c), start(HS_OS_CYCLES()) {
    prev_active = CycleCounter::active;
    // No latch on recursion or when it would close a parent cycle; retried on
    // the next entry.
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
   * @pre The scope spans less than one CYCCNT wrap (~7 s at 600 MHz).
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
 * @details The ISR is the sole writer; a foreground reader copies and
 *          reset()s with IRQs off.
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
 * @details Compiled in only under HS_PROFILE_ENABLE. Expands to a guard
 *          declaration: as the unbraced body of an `if` or loop it records ~0
 *          cycles. HS_OS_CYCLES() is 0 off Teensy.
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
 * @details For per-pixel, per-cell and per-face scopes the standard report
 * does not consume. Requires HS_PROFILE_DEEP_ENABLE on top of
 * HS_PROFILE_ENABLE, e.g. `just profile MeshFeedback 150 1`.
 */
#ifdef HS_PROFILE_DEEP_ENABLE
#define HS_PROFILE_DEEP(label) HS_PROFILE(label)
#else
#define HS_PROFILE_DEEP(label) ((void)0)
#endif
