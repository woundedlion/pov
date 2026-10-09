/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ── Multi-board simulator ───────────────────────────────────────────────────

struct SimEffect {
  uint64_t frames = 0;
  void advance_display() { ++frames; }
};

struct SimHandoff : pov::EffectHandoff<SimEffect> {
  Wake last;
  Wake apply_wake(const pov::WakeInputs &inputs) {
    last = pov::EffectHandoff<SimEffect>::apply_wake(inputs);
    if (last.adopted)
      last.live->frames = 0;
    return last;
  }
};

/**
 * @brief One simulated board: its SyncBoard engine plus the host-side state the
 *        simulator models around it.
 * @details The modeled state covers crystal offset/phase, the masked-IRQ latch,
 *          symbol-drop and EMI windows, the foreground build/commit model, and
 *          probes.
 */
struct SimBoard {
  SyncBoard board;
  SimHandoff handoff;
  pov::SyncPulseGate sync_pulse;
  pov::SubmitGate submit_gate;
  SimEffect instance; /**< Address the foreground publishes when built. */
  bool master = false;
  int32_t ppm = 0;      /**< Crystal offset, parts per million. */
  uint64_t phase0 = 0;  /**< Local cycle-counter offset at g = 0. */
  double next_tick = 0; /**< Next flywheel wake, in global cycles. */
  double tick_step = 0; /**< Wake period in global cycles. */
  /** Flywheel phase this board was seeded at, columns ahead of the master. */
  int32_t birth_cols = 0;
  /** Masked edges latch into one delayed delivery per pin. */
  bool edge_latched = false;
  std::vector<std::pair<uint64_t, uint64_t>> masks; /**< [from, to) global. */
  uint64_t drop_from = 0, drop_to = 0;              /**< Symbol drop window. */
  // Foreground model.
  uint32_t seen_gen = 0;
  uint32_t pending_pickup = 0;
  int32_t pending_index = -1;
  uint32_t pending_gen = 0;
  uint64_t pending_ready_g = 0;
  bool have_pending = false;
  uint64_t init_delay = 1000000; /**< ~1.7 ms construction time. */
  bool live = false;
  int32_t live_index = -1;
  uint64_t t = 0; /**< Frames shown: flips since this effect went live. */
  uint64_t swap_g = 0;
  bool trapped = false; /**< Commit deadline missed (device would HS_CHECK). */
  // Probes.
  uint64_t flips = 0;
  bool dark_now = true;
  float envelope = 0.0f;
  int32_t envelope_column = -1;

  /** @brief Reset foreground and probe state at a local reboot timestamp. */
  void reboot(uint32_t local_now) {
    board.seed(local_now, master);
    handoff.adopt(nullptr, 0);
    handoff.clear_pending();
    sync_pulse = {};
    submit_gate = {};
    instance = {};
    seen_gen = 0;
    pending_pickup = 0;
    pending_index = -1;
    pending_gen = 0;
    pending_ready_g = 0;
    have_pending = false;
    live = false;
    live_index = -1;
    t = 0;
    swap_g = 0;
    trapped = false;
    flips = 0;
    dark_now = true;
    envelope = 0.0f;
    envelope_column = -1;
  }

  /**
   * @brief Constructs a board wrapping a SyncBoard engine for config @p c.
   * @param c Sync configuration passed to the embedded SyncBoard.
   */
  explicit SimBoard(const Config &c) : board(c) {}
};

/**
 * @brief Event-driven multi-board simulator.
 * @details Advances the earliest-due flywheel wake across all boards, routes
 *          master pulses and injected EMI to downstream edge ISRs through the
 *          masked-IRQ model, and drives each board's foreground state.
 */
class Sim {
public:
  Config cfg;
  std::deque<SimBoard> boards;
  uint64_t g = 0;                            /**< Global time, cycles. */
  std::vector<std::pair<uint64_t, int>> emi; /**< (g, target), sorted. */
  size_t emi_pos = 0;

  /**
   * @brief Builds @p n boards (index 0 is the master) and seeds each flywheel.
   * @param c Shared sync configuration.
   * @param n Number of boards to construct.
   * @param ppm Per-board crystal offsets, parts per million (length @p n).
   * @param phase0 Common starting clock offset, cycles, at global time 0.
   * @details Each flywheel polls on a ⅛-column grid scaled by its ppm. Board i
   *          is born believing ZERO happened i·W/n columns ago; the master is
   *          born at offset 0.
   */
  Sim(const Config &c, int n, const int32_t *ppm, uint64_t phase0 = 0)
      : cfg(c) {
    const double step0 = double(c.cycles_per_half_rev) / (c.W / 2) / 8.0;
    for (int i = 0; i < n; ++i) {
      boards.emplace_back(c);
      SimBoard &b = boards.back();
      b.master = (i == 0);
      b.ppm = ppm[i];
      b.phase0 = phase0;
      b.birth_cols = i * c.W / n;
      b.tick_step = step0 * 1e6 / (1e6 + ppm[i]);
      b.next_tick = double(7 * (i + 1)) + (0.5 + 0.1 * (i % 5)) * b.tick_step;
      b.board.seed(local_now(b, 0) - static_cast<uint32_t>(b.birth_cols) *
                                         c.cycles_per_column(),
                   b.master);
    }
  }

  /**
   * @brief Computes board @p b's local 32-bit cycle counter at global cycle
   *        @p gg: phase offset plus crystal skew (ppm).
   * @param b The board whose local clock is evaluated.
   * @param gg Global time, cycles.
   * @return The board-local cycle count, truncated to 32 bits to model CYCCNT
   *         wrap.
   */
  static uint32_t local_now(const SimBoard &b, uint64_t gg) {
    const int64_t skew = static_cast<int64_t>(gg) * b.ppm / 1000000;
    return static_cast<uint32_t>(gg + b.phase0 + static_cast<uint64_t>(skew));
  }

  /**
   * @brief Tests whether global cycle @p at falls inside an IRQ-mask window.
   * @param b The board whose mask windows are checked.
   * @param at Global time, cycles.
   * @return The window's end cycle (when delivery resumes) if masked; 0 if
   *         unmasked.
   */
  uint64_t masked_until(const SimBoard &b, uint64_t at) const {
    for (const auto &w : b.masks)
      if (at >= w.first && at < w.second)
        return w.second;
    return 0;
  }

  /**
   * @brief Routes a sync-wire edge at global cycle @p at to downstream board
   *        @p target.
   * @param target Index of the downstream board receiving the edge.
   * @param at Global time of the edge, cycles.
   * @details Dropped if in the board's deafen window, latched/merged if masked,
   *          otherwise fed to its edge ISR at the board-local timestamp.
   */
  void deliver_edge(int target, uint64_t at) {
    SimBoard &b = boards[target];
    if (b.master)
      return; // master's edge ISR is not attached
    if (at >= b.drop_from && at < b.drop_to)
      return;
    if (masked_until(b, at)) {
      b.edge_latched = true; // single latched flag: merged, delayed
      return;
    }
    b.board.on_sync_edge(local_now(b, at));
  }

  /**
   * @brief Computes board @p b's next wake in global cycles.
   * @param b The board whose next wake is evaluated.
   * @return The scheduled wake, pushed to the mask end if its slot lands inside
   *         a mask (coalesced wakes).
   */
  double effective_tick(const SimBoard &b) const {
    const uint64_t m = masked_until(b, static_cast<uint64_t>(b.next_tick));
    return m ? double(m) : b.next_tick;
  }

  /**
   * @brief Advances global time to the earliest-due board wake, delivering any
   *        EMI edges before it, then runs that board's tick.
   */
  void step() {
    int bi = 0;
    double best = effective_tick(boards[0]);
    for (size_t i = 1; i < boards.size(); ++i) {
      const double e = effective_tick(boards[i]);
      if (e < best) {
        best = e;
        bi = static_cast<int>(i);
      }
    }
    const uint64_t tg = static_cast<uint64_t>(best + 0.5);
    // Deliver EMI edges due before this wake (edge ISRs run on their own).
    while (emi_pos < emi.size() && emi[emi_pos].first <= tg) {
      deliver_edge(emi[emi_pos].second, emi[emi_pos].first);
      ++emi_pos;
    }
    run_tick(boards[bi], tg);
    g = tg;
  }

  /**
   * @brief Steps the simulation until @p revs revolutions of global time have
   *        elapsed.
   * @param revs Number of revolutions to advance.
   */
  void run_revs(double revs) {
    const uint64_t until =
        g + static_cast<uint64_t>(revs * 2 * cfg.cycles_per_half_rev);
    while (g < until)
      step();
  }

  /**
   * @brief Steps the simulation until @p pred holds or @p max_revs elapse.
   * @tparam Pred Callable taking `Sim &` and returning bool.
   * @param pred Predicate evaluated after each step.
   * @param max_revs Maximum revolutions to advance before giving up.
   * @return True if @p pred returned true within the budget, false on timeout.
   */
  template <typename Pred>
  [[nodiscard]] bool run_until(Pred pred, double max_revs) {
    const uint64_t until =
        g + static_cast<uint64_t>(max_revs * 2 * cfg.cycles_per_half_rev);
    while (g < until) {
      step();
      if (pred(*this))
        return true;
    }
    return false;
  }

  /**
   * @brief Computes board @p i's current flywheel column position.
   * @param i Index of the board to sample.
   * @return The column position, evaluated at the board's own local clock at
   *         the present global time.
   */
  int32_t board_pos(int i) const {
    return flywheel(boards[i].board).position(local_now(boards[i], g));
  }

  double board_phase(int i) const {
    const SimBoard &board = boards[i];
    const Flywheel &fly = flywheel(board.board);
    const uint32_t now = local_now(board, g);
    const double elapsed =
        double(cfg.cycles_per_half_rev) - fly.cycles_to_next_boundary(now);
    double phase = boundary_column(fly.current_boundary(), cfg.W) +
                   elapsed * (cfg.W / 2) / cfg.cycles_per_half_rev;
    phase = std::fmod(phase, cfg.W);
    return phase < 0.0 ? phase + cfg.W : phase;
  }

  /**
   * @brief Computes the largest circular column distance of any locked board
   *        from the master.
   * @return Worst phase error, or infinity if a required board is unlocked.
   */
  double max_phase_err() const {
    const double master_phase = board_phase(0);
    double worst = 0.0;
    for (size_t i = 1; i < boards.size(); ++i) {
      if (lock(boards[i].board) != LockState::LOCKED)
        return std::numeric_limits<double>::infinity();
      const double direct =
          std::abs(board_phase(static_cast<int>(i)) - master_phase);
      const double d = std::min(direct, cfg.W - direct);
      if (d > worst)
        worst = d;
    }
    return worst;
  }

private:
  /**
   * @brief Runs one flywheel wake for board @p b at global cycle @p tg.
   * @param b The board being woken.
   * @param tg Global time of this wake, cycles.
   * @details Advances the poll grid, delivers any latched edge, drives
   *          board.tick(), fans master pulses to downstream boards, then steps
   *          the foreground build/commit model and updates probes.
   */
  void run_tick(SimBoard &b, uint64_t tg) {
    // Resume the grid at the next slot after this wake; masked slots coalesce.
    do {
      b.next_tick += b.tick_step;
    } while (b.next_tick <= double(tg));
    if (b.edge_latched && !masked_until(b, tg)) {
      b.board.on_sync_edge(local_now(b, tg)); // delayed, merged timestamp
      b.edge_latched = false;
    }

    const uint32_t now = local_now(b, tg);
    BurstSnapshot s;
    const BurstSnapshot *sp = nullptr;
    if (!b.master && b.board.claim_sync_burst(now, &s))
      sp = &s;
    TickActions a;
    pov::run_wake_sequence(
        b.sync_pulse, b.submit_gate, b.handoff,
        [&] { return a = b.board.tick(now, sp); },
        [&] {
          // Foreground model (build pickup at a flip, commit/join swaps).
          const uint32_t bw = b.board.build_word();
          const uint32_t gen = SyncBoard::build_gen_of(bw);
          const bool pickup = gen == b.pending_pickup && a.flip;
          b.pending_pickup = gen;
          if (gen != b.seen_gen && pickup) {
            b.seen_gen = gen;
            // Release + delete the outgoing instance.
            b.handoff.request_release();
            b.handoff.clear_pending();
            b.pending_index = SyncBoard::build_index_of(bw);
            b.pending_gen = gen;
            b.pending_ready_g = tg + b.init_delay;
            b.have_pending = true;
          }
          // Construction completes: the instance is published for the ISR to adopt.
          if (b.have_pending && b.pending_gen == gen && tg >= b.pending_ready_g)
            b.handoff.publish(&b.instance, b.pending_gen);

          return gen;
        },
        [&](bool high) {
          if (high && b.master)
            for (size_t j = 1; j < boards.size(); ++j)
              deliver_edge(static_cast<int>(j), tg);
        },
        [&] { b.trapped = true; },
        [&](SimEffect *, int32_t column) {
          b.envelope = b.board.effect_envelope(column, cfg.W);
          b.envelope_column = column;
        },
        [](pov::SubmitAction, SimEffect *, int32_t) { return true; });
    const auto &w = b.handoff.last;
    if (w.adopted) {
      b.live_index = b.pending_index;
      b.swap_g = tg;
    }
    b.live = w.live != nullptr;
    if (a.flip)
      ++b.flips;
    b.t = b.instance.frames;
    b.dark_now = w.dark;
  }
};

/**
 * @brief Predicate: every board has gone live (boot join complete).
 * @param s The simulation to inspect.
 * @return True iff no board is still pre-live.
 */
inline bool all_boards_live(Sim &s) {
  for (auto &b : s.boards)
    if (!b.live)
      return false;
  return true;
}

/**
 * @brief Runs the sim until every board has gone live (boot join complete).
 * @param sim The simulation to advance.
 * @param cfg The active config (supplies the join budget).
 * @return True if all boards went live within the join budget.
 */
inline bool boot_join(Sim &sim, const Config &cfg) {
  return sim.run_until(all_boards_live, double(cfg.join_grid_revs) + 2);
}

/**
 * @brief Runs the sim to one revolution before the train start.
 * @param sim The simulation to advance.
 * @param cfg The active config (supplies the effect-length budget).
 * @return True if the pre-train point was reached within budget.
 * @details The primary EPOCH copy follows about one revolution later.
 */
inline bool to_pre_train(Sim &sim, const Config &cfg) {
  return sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4);
}

// ── Scenario: clean 4-board run (boot join, phase, flips, wrap) ─────────────

/**
 * @brief Verifies a clean 4-board run: boards born out of phase lock, join live
 *        at the same boundary, and hold phase, frame counters and flip cadence
 *        through a 32-bit clock wrap with no rejections or traps.
 */
inline void test_sim_boot_and_phase() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 40};
  // Local clocks start below the 32-bit wrap so CYCCNT wraps mid-run (§12).
  Sim sim(cfg, 4, ppm, 0xFFFFFFFFull - 10ull * 2 * PERIOD + 12345);

  // Birth phase, before a single symbol: every downstream board is tens of
  // columns from the master, beyond what a LOCKED snap could close.
  for (int i = 1; i < 4; ++i) {
    const int32_t born = circ_dist(sim.board_pos(i), sim.board_pos(0), cfg.W);
    HS_EXPECT_GT(born, cfg.gate_cols);
    HS_EXPECT_GE(born, 40);
  }

  // All boards lock within the first revolution (two boundary symbols).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (lock(b.board) != LockState::LOCKED)
            return false;
        return true;
      },
      1.5));

  // Boot join: every board goes live at the SAME join-grid boundary with
  // identical effect and frame counter (no master head start).
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  for (int i = 0; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);
    // Swaps happen at each board's own crossing of the same boundary —
    // within a couple of columns of global time of each other.
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
  }

  // Run through the cycle-counter wrap and beyond; check phase + flip
  // cadence + frame-counter equality at stable mid-half instants.
  for (int slice = 0; slice < 6; ++slice) {
    uint64_t flips_before[4];
    for (int i = 0; i < 4; ++i)
      flips_before[i] = sim.boards[i].flips;
    sim.run_revs(4.0);
    // Sample at a stable point: master mid-first-half (~x=72), where every
    // board's ZERO flip for this rev has long settled.
    HS_EXPECT_TRUE(
        sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
    HS_EXPECT_LE(sim.max_phase_err(), 0.13);
    for (int i = 0; i < 4; ++i) {
      const uint64_t df = sim.boards[i].flips - flips_before[i];
      HS_EXPECT_GE(df, 8u); // ~2 flips/rev over the ≥4-rev slice
      HS_EXPECT_LE(df, 12u);
      HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    }
  }
  for (int i = 1; i < 4; ++i) {
    const Telemetry &tm = sim.boards[i].board.telemetry_snapshot();
    HS_EXPECT_EQ(tm.symbols_rejected_gate, 0u);
    HS_EXPECT_EQ(tm.symbols_discarded_invalid, 0u);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

/**
 * @brief Verifies eight boards acquire, join, and remain content-coherent.
 */
inline void test_sim_eight_board_boot_and_phase() {
  const Config cfg = test_config();
  const int32_t PPM[8] = {0, 20, -20, 40, -35, 15, -10, 30};
  Sim sim(cfg, 8, PPM);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (lock(b.board) != LockState::LOCKED)
            return false;
        return true;
      },
      1.5));
  HS_EXPECT_TRUE(boot_join(sim, cfg));

  sim.run_revs(8.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 0.13);
  for (size_t i = 1; i < sim.boards.size(); ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

// ── Scenario: epoch commit — lockstep advance, dark window, deadline ───────

/**
 * @brief Verifies epoch commit lockstep: boards hold the outgoing effect at
 *        zero envelope through announce, go dark for the K-rev construction
 *        window, then swap at the same boundary with frame counters re-zeroed
 *        together; the cadence holds across a second epoch.
 */
inline void test_sim_epoch_commit() {
  const Config cfg = test_config(2);
  const int32_t ppm[4] = {0, 30, -25, 10};
  Sim sim(cfg, 4, ppm);

  HS_EXPECT_TRUE(boot_join(sim, cfg));

  auto expect_envelopes = [&](int32_t column) {
    HS_EXPECT_TRUE(sim.run_until(
        [=](Sim &s) {
          for (const auto &board : s.boards)
            if (board.envelope_column != column)
              return false;
          return true;
        },
        1.1));
    for (const auto &board : sim.boards)
      HS_EXPECT_EQ(board.envelope, sim.boards[0].envelope);
  };
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).rev_in_effect ==
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect)));
  expect_envelopes(72);
  HS_EXPECT_GT(sim.boards[0].envelope, 0.0f);
  HS_EXPECT_LT(sim.boards[0].envelope, 1.0f);

  // Run to the train start. Through the announce phase (the first R revs of
  // the B+R+K countdown) every board retains the outgoing effect at zero
  // envelope.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return content(s.boards[0].board).commit_pending; },
      double(cfg.revs_per_effect) + 2));
  sim.run_revs(1.0); // mid-announce
  expect_envelopes(72);
  HS_EXPECT_EQ(sim.boards[0].envelope, 0.0f);
  for (auto &b : sim.boards) {
    HS_EXPECT_TRUE(content(b.board).commit_pending);
    HS_EXPECT_FALSE(b.dark_now);
  }
  // Then all enter the K-revolution construction window together.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).commit_in_revs <= s.cfg.commit_revs;
      },
      double(cfg.epoch_repeats) + 1));
  sim.run_revs(1.0); // mid-construction
  for (auto &b : sim.boards) {
    HS_EXPECT_TRUE(content(b.board).commit_pending);
    HS_EXPECT_TRUE(b.dark_now);
  }

  // Commit: all four swap to effect 1 at the same boundary, frame counters
  // reset together; the EPOCH redundancy repeats were refractory-ignored.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      4.0));
  for (int i = 1; i < 4; ++i) {
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
    HS_EXPECT_EQ(
        sim.boards[i].board.telemetry_snapshot().epochs_refractory_ignored,
        static_cast<uint32_t>(cfg.epoch_repeats));
  }
  expect_envelopes(72);
  HS_EXPECT_GT(sim.boards[0].envelope, 0.0f);
  HS_EXPECT_LT(sim.boards[0].envelope, 1.0f);
  // Post-epoch: full content coherence (index AND t) including the master.
  sim.run_revs(3.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, 1);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  }

  // Second epoch keeps the cadence (roster wraps mod effect_count).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 0)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 6));
}

/** @brief Verifies the master schedules each epoch from its roster duration. */
inline void test_sim_variable_effect_durations() {
  uint32_t effect_revolutions[4] = {24, 52, 32, 64};
  Config cfg = test_config();
  cfg.set_effect_revolutions(effect_revolutions);
  const int32_t ppm[4] = {0, 30, -25, 10};
  Sim sim(cfg, 4, ppm);

  HS_EXPECT_TRUE(boot_join(sim, cfg));
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(effect_revolutions[0]) + 6));

  sim.run_revs(45.0);
  for (auto &b : sim.boards)
    HS_EXPECT_EQ(b.live_index, 1);

  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 2)
            return false;
        return true;
      },
      double(effect_revolutions[1] + cfg.epoch_repeats + cfg.commit_revs) -
          45.0 + 1.0));
}

/**
 * @brief Verifies the widest beacon-to-beacon gap a receiver sees across
 *        commits matches Config::commit_beacon_gap_revs().
 */
inline void test_sim_commit_beacon_gap() {
  const uint32_t effect_revolutions[2] = {49, 50};
  Config cfg = test_config(2);
  cfg.beacon_period_revs = 16;
  cfg.set_effect_revolutions(effect_revolutions);
  cfg.rejoin_budget_revs = cfg.rejoin_bound_revs();
  HS_EXPECT_TRUE(cfg.valid() == nullptr);
  const int32_t ppm[2] = {0, 0};
  Sim sim(cfg, 2, ppm);
  HS_EXPECT_TRUE(boot_join(sim, cfg));

  uint32_t seen = sim.boards[1].board.telemetry_snapshot().beacons_ok;
  uint64_t last_at = 0;
  uint64_t widest = 0;
  (void)sim.run_until(
      [&](Sim &s) {
        const uint32_t ok = s.boards[1].board.telemetry_snapshot().beacons_ok;
        if (ok != seen) {
          if (last_at != 0 && s.g - last_at > widest)
            widest = s.g - last_at;
          last_at = s.g;
          seen = ok;
        }
        return false;
      },
      double(effect_revolutions[0] + effect_revolutions[1]) + 12);
  const uint64_t rev = 2ull * PERIOD;
  const uint64_t expected = cfg.commit_beacon_gap_revs(0);
  HS_EXPECT_EQ(expected, 22u);
  HS_EXPECT_EQ(cfg.commit_beacon_gap_revs(1), 7u);
  HS_EXPECT_EQ((widest + rev / 2) / rev, expected);
}

/**
 * @brief Verifies an effect whose construction outruns the K-revolution window
 *        traps (HS_CHECK on the device) and never silently skews the show
 *        (§6.1).
 */
inline void test_sim_commit_deadline_trap() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 0, 0, 0};
  Sim sim(cfg, 4, ppm);
  sim.boards[2].init_delay =
      static_cast<uint64_t>(3) * 2 * PERIOD; // 3 revs > K
  // Boot joins are not deadline-bound: board 2 goes live at a later join-grid
  // boundary without trapping.
  sim.run_revs(double(cfg.join_grid_revs) * 3);
  HS_EXPECT_TRUE(sim.boards[2].live);
  HS_EXPECT_FALSE(sim.boards[2].trapped);
  // The epoch commit IS deadline-bound.
  HS_EXPECT_TRUE(sim.run_until([](Sim &s) { return s.boards[2].trapped; },
                               double(cfg.revs_per_effect) + 6));
  for (int i : {0, 1, 3})
    HS_EXPECT_FALSE(sim.boards[i].trapped);
}

/** @brief Pickup traps with 0.4 half-revs of commit slack, but accepts 1.5. */
inline void test_sim_commit_pickup_budget() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 0, 0, 0};
  for (bool over_budget : {false, true}) {
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim, cfg));
    sim.boards[2].init_delay = static_cast<uint64_t>(
        (2.0 * cfg.commit_revs - (over_budget ? 0.4 : 1.5)) * PERIOD);
    sim.run_revs(double(cfg.revs_per_effect) + 6);
    HS_EXPECT_EQ(sim.boards[2].trapped, over_budget);
    if (!over_budget)
      HS_EXPECT_EQ(sim.boards[2].live_index, sim.boards[0].live_index);
  }
}

// ── Scenario: masked-IRQ windows (§4.1, §5.2) ───────────────────────────────

/**
 * @brief Verifies masked-IRQ windows (§4.1, §5.2): boundary masks degrade a
 *        symbol to missed, never misclassified, mid-rev masks coalesce wakes,
 *        and a master masked across its own boundary self-censors; phase and
 *        content stay coherent.
 */
inline void test_sim_masked_windows() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 25, -25, 15};
  Sim sim(cfg, 4, ppm);

  // Board 2: recurring masks over the ZERO boundary swallow wakes and merge
  // burst edges, so its decoder sees truncated counts. Masks start after boot
  // join so acquisition is clean.
  const uint64_t rev = 2ull * PERIOD;
  for (int k = 8; k < 28; ++k) {
    const uint64_t b0 = k * rev; // master ZERO crossings ≈ k·rev (ppm 0)
    sim.boards[1].masks.push_back({b0 - COL / 4, b0 + COL / 4});
    sim.boards[2].masks.push_back({b0 - COL / 2, b0 + 2 * COL + COL / 2});
  }
  // Board 3: mid-revolution masks only coalesce wakes.
  for (int k = 8; k < 28; ++k) {
    const uint64_t m0 = k * rev + 40 * COL;
    sim.boards[3].masks.push_back({m0, m0 + 5 * COL});
  }

  sim.run_revs(30.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));

  // Truncated bursts were discarded (count telemetry), never accepted as
  // the wrong boundary: phase stays sub-column-ish and content equal.
  const Telemetry &tm2 = sim.boards[2].board.telemetry_snapshot();
  const Telemetry &tm1 = sim.boards[1].board.telemetry_snapshot();
  HS_EXPECT_GT(tm1.symbols_accepted, 20u);
  HS_EXPECT_EQ(tm1.symbols_discarded_invalid, 0u);
  HS_EXPECT_GT(tm2.symbols_discarded_invalid, 5u);
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(lock(sim.boards[i].board), LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }

  // Master masked across its own boundary: it self-censors (or truncates)
  // rather than emitting late; downstream coasts one half-rev and re-snaps.
  Sim sim2(cfg, 2, ppm);
  sim2.run_revs(8.0);
  const uint64_t b0 = (static_cast<uint64_t>(sim2.g / rev) + 2) * rev;
  sim2.boards[0].masks.push_back({b0 - COL / 4, b0 + 2 * COL});
  sim2.run_revs(6.0);
  const Telemetry &tm0 = sim2.boards[0].board.telemetry_snapshot();
  HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted +
                   tm0.beacons_overrun_dropped + tm0.boundary_bursts_dropped,
               0u);
  HS_EXPECT_LE(sim2.max_phase_err(), 2);
  HS_EXPECT_EQ(sim2.boards[1].board.telemetry_snapshot().symbols_rejected_gate,
               0u);
}

// ── Scenario: EMI on the sync wire (§5.2, §5.3, §9.1) ───────────────────────

/**
 * @brief Verifies EMI on the sync wire (§5.2, §5.3, §9.1): isolated spurious
 *        edges form valid HALF symbols the LOCKED gate rejects, edges injected
 *        inside or near real bursts corrupt the count to invalid (discarded
 *        whole) — none unlock or desync the show.
 */
inline void test_sim_emi() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -30, 5};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;

  // Isolated spurious edges on board 1, away from boundaries: each is a valid
  // HALF symbol the LOCKED gate must reject.
  uint32_t lcg = 12345;
  for (int k = 8; k < 40; ++k) {
    lcg = lcg * 1664525u + 1013904223u;
    const uint64_t off = (20 + lcg % 100) * static_cast<uint64_t>(COL);
    sim.emi.push_back({k * rev + off, 1});
  }
  // Two edges injected INSIDE master ZERO bursts on board 2 (count 3 → 4):
  // even count = invalid, discarded whole; the crossing flip covers it.
  sim.emi.push_back({10 * rev + COL, 2});
  sim.emi.push_back({14 * rev + COL, 2});
  // The nearby forged edge merges with the real HALF into an invalid count.
  sim.emi.push_back({12 * rev + (PERIOD - 2ull * COL), 3});
  std::sort(sim.emi.begin(), sim.emi.end());

  sim.run_revs(42.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));

  const Telemetry &tm1 = sim.boards[1].board.telemetry_snapshot();
  HS_EXPECT_GT(tm1.symbols_rejected_gate, 20u); // isolated EMI all rejected
  const Telemetry &tm2 = sim.boards[2].board.telemetry_snapshot();
  HS_EXPECT_GE(tm2.symbols_discarded_invalid, 2u); // corrupted bursts dropped
  HS_EXPECT_GE(
      sim.boards[3].board.telemetry_snapshot().symbols_discarded_invalid, 1u);
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(lock(sim.boards[i].board), LockState::LOCKED);
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  }
}

// ── Scenario: dropped symbols → coast; dropped epoch → beacon fix (§6.3) ───

/**
 * @brief Verifies dropped-symbol recovery (§6.3): a multi-rev symbol gap coasts
 *        and re-snaps, and a board that misses a whole EPOCH train stays dark
 *        until the next index beacon corrects it and it rejoins on the join
 *        grid.
 */
inline void test_sim_drops_and_missed_epoch() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 35, -35, 20};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;

  HS_EXPECT_TRUE(boot_join(sim, cfg));

  // Board 1 hears nothing for 2 revolutions: it coasts, then re-snaps.
  sim.boards[1].drop_from = sim.g + rev;
  sim.boards[1].drop_to = sim.g + 3 * rev;
  sim.run_revs(5.0);
  HS_EXPECT_GE(sim.boards[1].board.telemetry_snapshot().max_coast_halves, 4u);
  HS_EXPECT_EQ(lock(sim.boards[1].board), LockState::LOCKED);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);

  // Board 3 loses its wire for the entire EPOCH train; the next index beacon
  // corrects it (§6.3.2).
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4));
  sim.boards[3].drop_from = sim.g;
  sim.boards[3].drop_to = sim.g + 7 * rev; // covers B..B+K and the repeats
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return s.boards[0].live_index == 1 && s.boards[1].live_index == 1;
      },
      8.0));
  HS_EXPECT_EQ(sim.boards[3].live_index, 0);
  HS_EXPECT_EQ(sim.boards[3].envelope, 0.0f);
  // Correction: ≤ one beacon period + join grid after the wire returns.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.boards[3].live_index == 1; },
                    double(cfg.beacon_period_revs + cfg.join_grid_revs) + 6));
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_EQ(content(sim.boards[3].board).rev_in_effect & 63,
               content(sim.boards[0].board).rev_in_effect & 63);
  HS_EXPECT_EQ(
      sim.boards[3].board.telemetry_snapshot().beacon_index_corrections, 1u);
  HS_EXPECT_FALSE(sim.boards[3].trapped);
}

// ── Scenario: mid-show reboot — fail-dark, rejoin at correct effect ─────────

/**
 * @brief Verifies mid-show reboot: a board reseeded with fresh engine state
 *        re-acquires phase within ~one boundary symbol, stays dark through
 *        ACQUIRE, then rejoins at the master's current effect from the next
 *        beacon + join-grid boundary — never a wrong frame in between.
 */
inline void test_sim_reboot(const Config &cfg) {
  const int32_t ppm[4] = {0, 15, -15, 30};
  Sim sim(cfg, 4, ppm);
  sim.run_revs(12.0);

  const uint64_t reboot_at = sim.g;
  // Reboot board 2 mid-show: fresh engine state, no identity assumption.
  SimBoard &b2 = sim.boards[2];
  b2.reboot(Sim::local_now(b2, sim.g));

  HS_EXPECT_EQ(lock(b2.board), LockState::ACQUIRE);

  // Phase re-acquires within ~one boundary symbol; the board stays pre-live
  // (hence dark) at every step of ACQUIRE.
  bool lit_in_acquire = false;
  HS_EXPECT_TRUE(sim.run_until(
      [&](Sim &s) {
        lit_in_acquire |= s.boards[2].live;
        return lock(s.boards[2].board) == LockState::LOCKED;
      },
      1.5));
  HS_EXPECT_FALSE(lit_in_acquire);
  // Identity from the next beacon, display from the next join-grid boundary;
  // this reboot lands clear of the commit window.
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.boards[2].live; },
                    double(cfg.beacon_period_revs + cfg.join_grid_revs) + 4));
  HS_EXPECT_LE(sim.g - reboot_at,
               uint64_t(cfg.beacon_period_revs + cfg.join_grid_revs + 1) * 2 *
                   PERIOD);
  HS_EXPECT_EQ(b2.live_index, sim.boards[0].live_index);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_LE(sim.max_phase_err(), 2);
  HS_EXPECT_EQ(content(b2.board).rev_in_effect & 63,
               content(sim.boards[0].board).rev_in_effect & 63);
}

// ── Scenario: forged plausible burst (§8.4 spurious-flip hole, closed) ──────

/**
 * @brief Verifies the forged-plausible-burst defense (§8.4): the strongest
 *        cheap spurious-flip attack is held as a suspect and rejected, never
 *        snapped or flipped; flip cadence and content stay intact.
 */
inline void test_sim_forged_burst() {
  const Config cfg = test_config();
  const int32_t ppm[2] = {0, 10};
  Sim sim(cfg, 2, ppm);
  sim.run_revs(10.0);
  const uint64_t rev = 2ull * PERIOD;

  // An isolated valid-count (HALF) burst 30 columns past a real ZERO boundary,
  // clear of the gap-timeout window and the quiet-before guard (§5.3).
  const uint64_t b0 = (sim.g / rev + 2) * rev; // a future master ZERO
  sim.emi.push_back({b0 + 30 * static_cast<uint64_t>(COL), 1});
  sim.emi_pos = 0;
  std::sort(sim.emi.begin(), sim.emi.end());

  const uint64_t flips_before = sim.boards[1].flips;
  const uint32_t rejected_before =
      sim.boards[1].board.telemetry_snapshot().symbols_rejected_gate;
  sim.run_revs(4.0);
  HS_EXPECT_GT(sim.boards[1].board.telemetry_snapshot().symbols_rejected_gate,
               rejected_before);
  // Flip cadence unbroken: 2 per revolution, no extra content advance.
  const uint64_t df = sim.boards[1].flips - flips_before;
  HS_EXPECT_EQ(df, 8u);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  HS_EXPECT_EQ(sim.boards[1].t, sim.boards[0].t);
}

// ── Scenario: epoch-repeat lockstep (§6.3.1) ────────────────────────────────

/**
 * @brief Verifies a board that misses the EPOCH primary copy at B but accepts
 *        a repeat at B+j commits at the same absolute boundary as its peers,
 *        with equal frame counters (§6.3.1).
 * @details Covers a downstream board deafened for the primary copy and a
 *          master that self-censors it.
 */
inline void test_sim_epoch_repeat_lockstep() {
  const Config cfg = test_config();
  const uint64_t rev = 2ull * PERIOD;

  /**
   * @brief Predicate: all boards have committed to effect 1.
   * @param s The simulation to inspect.
   * @return True iff every board's live_index is 1.
   */
  auto all_on_effect_1 = [](Sim &s) {
    for (auto &b : s.boards)
      if (b.live_index != 1)
        return false;
    return true;
  };
  /**
   * @brief Asserts every downstream board committed at the same boundary as the
   *        master (≤3 columns apart, no trap) and holds an equal frame counter
   *        afterward.
   * @param sim The simulation to inspect (advanced two revs while checking).
   */
  auto expect_lockstep = [&](Sim &sim) {
    for (int i = 1; i < 4; ++i) {
      // Commits land at each board's own crossing of the same boundary —
      // within a couple of columns of global time, never a revolution apart.
      const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                         static_cast<int64_t>(sim.boards[0].swap_g);
      HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
      HS_EXPECT_FALSE(sim.boards[i].trapped);
    }
    sim.run_revs(2.0);
    HS_EXPECT_TRUE(
        sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
    for (int i = 1; i < 4; ++i)
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
  };
  const double commit_revs_max =
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) + 8;

  // Sub-scenario 1: board 2 misses ONLY the primary copy at B.
  {
    const int32_t ppm[4] = {0, 20, -20, 10};
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim, cfg));
    HS_EXPECT_TRUE(to_pre_train(sim, cfg));
    sim.boards[2].drop_from = sim.g + rev - 4 * COL;
    sim.boards[2].drop_to = sim.g + rev + rev / 4; // before the B+1 repeat
    HS_EXPECT_TRUE(sim.run_until(all_on_effect_1, commit_revs_max));
    expect_lockstep(sim);
  }

  // Sub-scenario 2: the master self-censors the primary copy (masked across
  // its own train-start boundary, a designed event — §5.2): it must not
  // schedule a commit its peers cannot match.
  {
    const int32_t ppm[4] = {0, 15, -25, 30};
    Sim sim(cfg, 4, ppm);
    HS_EXPECT_TRUE(boot_join(sim, cfg));
    HS_EXPECT_TRUE(to_pre_train(sim, cfg));
    const uint64_t b0 = sim.g + rev;
    sim.boards[0].masks.push_back({b0 - COL / 4, b0 + 2 * COL});
    HS_EXPECT_TRUE(sim.run_until(all_on_effect_1, commit_revs_max));
    const Telemetry &tm0 = sim.boards[0].board.telemetry_snapshot();
    HS_EXPECT_GT(tm0.emit_censored + tm0.emit_aborted +
                     tm0.beacons_overrun_dropped + tm0.boundary_bursts_dropped,
                 0u);
    expect_lockstep(sim);
  }
}

// ── Scenario: schedule-counter resync from the beacon rev cross-check ───────

/**
 * @brief Verifies the §6.4 beacon rev cross-check resyncs a slipped
 *        schedule counter within one beacon period, restoring exact
 *        j-inference before the next train.
 * @details A board whose rev_in_effect slipped against the master infers j
 *          wrongly at every later epoch train. The cross-check corrects via the
 *          signed mod-64 difference, leaving content t untouched.
 */
inline void test_sim_rev_resync() {
  const Config cfg = test_config();
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);
  const uint64_t rev = 2ull * PERIOD;
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  // Settle past the boot beacons, then slip board 3's counter by +2: left
  // uncorrected it would commit 2 revolutions early.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) { return content(s.boards[0].board).rev_in_effect == 5; },
      12.0));
  content_mut(sim.boards[3].board).rev_in_effect += 2;

  // Detected and corrected at the next beacon (rev 9), within one period.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[3].board).rev_in_effect ==
               content(s.boards[0].board).rev_in_effect;
      },
      double(cfg.beacon_period_revs) + 2));
  HS_EXPECT_GE(sim.boards[3].board.telemetry_snapshot().beacon_rev_mismatches,
               1u);

  // The next epoch commits in lockstep even though board 3 ALSO misses the
  // primary copy, taking the repeat path that depends on j-inference.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        return content(s.boards[0].board).rev_in_effect >=
               s.cfg.revs_per_effect - 1;
      },
      double(cfg.revs_per_effect) + 4));
  sim.boards[3].drop_from = sim.g + rev - 4 * COL;
  sim.boards[3].drop_to = sim.g + rev + rev / 4;
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats)) + 8));
  for (int i = 1; i < 4; ++i) {
    const int64_t dg = static_cast<int64_t>(sim.boards[i].swap_g) -
                       static_cast<int64_t>(sim.boards[0].swap_g);
    HS_EXPECT_LE(dg < 0 ? -dg : dg, int64_t(3) * COL);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
  sim.run_revs(2.0);
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
}

// ── Scenario: 6-bit rev counter wraps WITHIN one effect (§6.4 mod-64) ────────

/**
 * @brief Verifies normal play stays in lockstep — and the epoch still commits
 *        correctly — as rev_in_effect rolls through its 6-bit (mod-64) residue
 *        within a single effect.
 * @details The beacon carries rev mod 64; the cross-check compares
 *          f.rev_count against `content_tracker.rev_in_effect & 63`.
 */
inline void test_sim_rev_wrap_within_effect() {
  Config cfg = test_config();
  cfg.revs_per_effect =
      90; // > 64: rev_in_effect wraps its 6-bit residue mid-effect
  const int32_t ppm[4] = {0, 20, -20, 10};
  Sim sim(cfg, 4, ppm);

  // Boot: all four boards live on effect 0.
  HS_EXPECT_TRUE(boot_join(sim, cfg));
  for (int i = 0; i < 4; ++i)
    HS_EXPECT_EQ(sim.boards[i].live_index, 0);

  // Advance ACROSS the 63→0 seam and hold past it, sampling lockstep at stable
  // mid-half instants (master ~x=72).
  bool crossed_seam = false;
  bool reached_end = false;
  for (int sample = 0; sample < 40; ++sample) {
    HS_EXPECT_TRUE(
        sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
    const uint32_t rev = content(sim.boards[0].board).rev_in_effect;
    HS_EXPECT_LE(sim.max_phase_err(), 2);
    for (int i = 0; i < 4; ++i) {
      HS_EXPECT_EQ(sim.boards[i].live_index, 0);
      HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
      HS_EXPECT_FALSE(sim.boards[i].trapped);
      // The schedule counter tracks the master's exactly through the wrap, and
      // the &63 cross-check raises no spurious rev mismatch.
      HS_EXPECT_EQ(content(sim.boards[i].board).rev_in_effect, rev);
      HS_EXPECT_EQ(
          sim.boards[i].board.telemetry_snapshot().beacon_rev_mismatches, 0u);
    }
    if (rev >= 64)
      crossed_seam = true;
    if (rev >= 80) {
      reached_end = true;
      break;
    }
    sim.run_revs(3.0);
  }
  HS_EXPECT_TRUE(crossed_seam); // the run actually exercised rev_in_effect ≥ 64

  HS_EXPECT_TRUE(reached_end);
  // The epoch still commits in lockstep at the post-seam effect boundary:
  // on_epoch_symbol infers j from a rev whose 6-bit residue has wrapped.
  HS_EXPECT_TRUE(sim.run_until(
      [](Sim &s) {
        for (auto &b : s.boards)
          if (b.live_index != 1)
            return false;
        return true;
      },
      double(cfg.revs_per_effect) + 8));
  HS_EXPECT_TRUE(
      sim.run_until([](Sim &s) { return s.board_pos(0) == 72; }, 1.1));
  for (int i = 1; i < 4; ++i) {
    HS_EXPECT_EQ(sim.boards[i].live_index, sim.boards[0].live_index);
    HS_EXPECT_EQ(sim.boards[i].t, sim.boards[0].t);
    HS_EXPECT_FALSE(sim.boards[i].trapped);
  }
}

// ── Content protocol units: EPOCH burst + boundary fold (§6.3.1) ────────────

/**
 * @brief Verifies §6.3.1 j-inference stays correct when an EPOCH burst is
 *        consumed in the SAME tick() its boundary is folded.
 * @details on_epoch_symbol infers j = rev_in_effect − revs_per_effect, so it
 *          must observe the post-fold rev. Every copy j of the train commits at
 *          the same absolute B+R+K boundary whether the fold was deferred into
 *          the burst's tick or applied a tick earlier.
 */
inline void test_epoch_same_tick_burst_fold() {
  const Config cfg = test_config();
  const uint32_t RPE = cfg.revs_per_effect;
  const uint32_t R = static_cast<uint32_t>(cfg.epoch_repeats);
  const uint32_t K = cfg.commit_revs;

  for (uint32_t j = 1; j <= R; ++j) {
    SyncBoard board(cfg);
    constexpr uint32_t START = 1000u;
    board.seed(START, false);
    flywheel_mut(board).force_lock();
    content_mut(board).identity_known = true;
    content_mut(board).rev_in_effect = RPE + j - 1;
    const uint32_t ZERO = START + 2u * cfg.cycles_per_half_rev;
    const uint32_t COL = cfg.cycles_per_column();
    const BurstSnapshot burst{5, ZERO, ZERO + 4u * COL};
    board.tick(ZERO + 4u * COL + cfg.gap_timeout_cycles(), &burst);
    HS_EXPECT_TRUE(content(board).commit_pending);
    HS_EXPECT_EQ(content(board).commit_in_revs, K + R - j);
  }

  // Returns the absolute rev at which a tracker hearing copy j commits
  // (computed: on_zero_crossing zeroes rev_in_effect on commit). same_tick
  // defers the fold into the burst's tick; otherwise it landed a tick earlier.
  auto commit_rev_for = [&](uint32_t j, bool same_tick) -> uint32_t {
    ContentTracker c;
    c.identity_known = true;
    c.effect_index = 0;
    if (same_tick) {
      c.rev_in_effect = RPE + j - 1;            // B+j fold still pending…
      HS_EXPECT_FALSE(c.on_zero_crossing(cfg)); // …backstop apply_flip folds it
    } else {
      c.rev_in_effect = RPE + j; // already folded a tick earlier
    }
    const uint32_t base_rev = c.rev_in_effect; // == RPE + j either way
    HS_EXPECT_EQ(base_rev, RPE + j);
    HS_EXPECT_TRUE(c.on_epoch_symbol(cfg)); // opens the commit window
    HS_EXPECT_EQ(c.commit_in_revs, K + R - j);
    uint32_t count = 0;
    while (count <= RPE) {
      const bool committed = c.on_zero_crossing(cfg);
      ++count;
      if (committed)
        break;
    }
    return base_rev + count;
  };

  // Every copy commits at the absolute B+R+K boundary, independent of j and of
  // whether the burst shared its tick with the fold.
  for (uint32_t j = 0; j <= R; ++j) {
    HS_EXPECT_EQ(commit_rev_for(j, /*same_tick=*/true), RPE + R + K);
    HS_EXPECT_EQ(commit_rev_for(j, /*same_tick=*/false), RPE + R + K);
  }

  // On the pre-fold rev, on_epoch_symbol infers a repeat copy one short (j−1),
  // scheduling commit_in_revs a revolution too large.
  for (uint32_t j = 1; j <= R; ++j) {
    ContentTracker pre;
    pre.identity_known = true;
    pre.rev_in_effect = RPE + j - 1; // fold NOT applied before the inference
    HS_EXPECT_TRUE(pre.on_epoch_symbol(cfg));
    HS_EXPECT_EQ(pre.commit_in_revs, K + R - (j - 1)); // one too large
  }
}

/** @brief Verifies EPOCH symbols stay rejected for the full refractory window. */
inline void test_epoch_refractory_window() {
  const Config cfg = test_config();
  ContentTracker content;
  content.identity_known = true;
  content.effect_index = 0;
  content.rev_in_effect = cfg.revs_per_effect;

  HS_EXPECT_TRUE(content.on_epoch_symbol(cfg));
  HS_EXPECT_EQ(content.refractory_revs_left, cfg.refractory_revs);
  for (uint32_t rev = 1; rev < cfg.refractory_revs; ++rev) {
    HS_CONTEXT("refractory rev", rev);
    content.on_zero_crossing(cfg);
    HS_EXPECT_FALSE(content.on_epoch_symbol(cfg));
    HS_EXPECT_EQ(content.refractory_revs_left, cfg.refractory_revs - rev);
  }
  content.on_zero_crossing(cfg);
  HS_EXPECT_EQ(content.refractory_revs_left, 0u);
  HS_EXPECT_TRUE(content.on_epoch_symbol(cfg));
}

/**
 * @brief Pins ContentTracker::construction_opens and ::constructing directly:
 *        the window opens exactly once and lasts exactly K revolutions, for
 *        every repeat copy the board may have heard.
 * @details The window is anchored to the absolute commit boundary, so a board
 *          that heard the last repeat (announce phase already spent) opens it
 *          at the accept itself.
 */
inline void test_construction_window_predicates() {
  const Config cfg = test_config();
  const uint32_t RPE = cfg.revs_per_effect;
  const uint32_t R = static_cast<uint32_t>(cfg.epoch_repeats);
  const uint32_t K = cfg.commit_revs;

  for (uint32_t j = 0; j <= R; ++j) {
    HS_CONTEXT("copy", static_cast<long long>(j));
    ContentTracker c;
    c.identity_known = true;
    c.effect_index = 0;
    c.rev_in_effect = RPE + j;

    // Nothing scheduled: neither predicate can fire off commit_pending.
    HS_EXPECT_FALSE(c.construction_opens(cfg));
    HS_EXPECT_FALSE(c.constructing(cfg));

    HS_EXPECT_TRUE(c.on_epoch_symbol(cfg));
    const uint32_t announce = R - j; // crossings before the window opens
    HS_EXPECT_EQ(c.commit_in_revs, K + announce);
    HS_EXPECT_EQ(c.construction_opens(cfg), announce == 0);
    HS_EXPECT_EQ(c.constructing(cfg), announce == 0);

    uint32_t opens = c.construction_opens(cfg) ? 1u : 0u;
    uint32_t dark = c.constructing(cfg) ? 1u : 0u;
    uint32_t steps = 0;
    bool committed = false;
    while (!committed && steps <= K + R + 1) {
      committed = c.on_zero_crossing(cfg);
      ++steps;
      // Opening is a moment inside the window, never outside it.
      HS_EXPECT_TRUE(!c.construction_opens(cfg) || c.constructing(cfg));
      opens += c.construction_opens(cfg) ? 1u : 0u;
      dark += c.constructing(cfg) ? 1u : 0u;
    }
    HS_EXPECT_TRUE(committed);
    HS_EXPECT_EQ(steps, K + announce);
    HS_EXPECT_EQ(opens, 1u);
    HS_EXPECT_EQ(dark, K);

    // The commit closes the window; a board past it is lit again.
    HS_EXPECT_FALSE(c.construction_opens(cfg));
    HS_EXPECT_FALSE(c.constructing(cfg));
  }
}
