/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * @file pov_handoff.h
 * @brief Pure, host-testable effect-handoff state machine for POVSegmented.
 *
 * Split out of pov_segmented.h (Arduino-only) so the concurrency glue between
 * the foreground effect builder and the flywheel ISR — the teardown handshake,
 * the acquire/release publish/adopt of the pending effect, the consumed-
 * generation gate, and the display-window (clip) alternation — is unit-testable
 * on the host, exactly as pov_sync.h is for the sync protocol.
 *
 * Ownership model (spec §8): the foreground constructs and deletes effect
 * instances; the ISR only ever dereferences the instance it has been handed via
 * live(). The commit-time generation check is the guard against a use-after-free
 * of a deleted instance, so the exact memory orderings below are load-bearing:
 *
 *   - publish() stores pending_gen (relaxed) then pending_effect (release);
 *     the release store publishes the (effect, gen) pair and orders every
 *     constructor/draw write before the ISR's acquire load can observe it.
 *   - pending_acquire() loads pending_effect (acquire) then pending_gen
 *     (relaxed); the acquire pairs with publish()'s release.
 *   - request_release()/service_release() are a counter handshake: the
 *     foreground bumps req, the ISR drops the live pointer and copies req to
 *     ack. release_complete() is the foreground's wait predicate.
 *
 * Free of FastLED/Teensy/timer dependencies; templated on the effect pointee so
 * the host tests drive it with a stand-in type.
 */
#pragma once

#include <atomic>
#include <cstdint>

namespace pov {

/**
 * @brief One flywheel wake's sync-layer decisions, as the handoff consumes them.
 */
struct WakeInputs {
  bool commit = false;        /**< B+R+K epoch deadline reached (spec §6.1). */
  bool join_boundary = false; /**< ZERO crossing on the join grid. */
  bool dark = false;          /**< The sync layer wants black this wake. */
  bool flip = false;          /**< A display window opened this wake. */
  bool zero_crossing = false; /**< The opening flip crossed ZERO. */
  uint32_t wire_gen = 0;      /**< Build generation advertised by the wire. */
};

/**
 * @brief Foreground↔ISR effect-handoff state machine.
 * @tparam T Effect pointee type (the device instantiates it with Effect; tests
 *           use a stand-in).
 *
 * Single-writer per field (spec §8): live_effect is ISR-written; pending_effect
 * and pending_gen are foreground-written (published under a brief interrupts-off
 * bracket on the device); the release_req/release_ack pair is the teardown
 * handshake (foreground bumps req, ISR acks); window_sweeps_left is ISR-written,
 * foreground-read.
 */
template <class T> class EffectHandoff {
public:
  /**
   * @brief An acquire-loaded (effect, generation) snapshot of the pending slot.
   */
  struct Pending {
    T *effect;    /**< Constructed instance awaiting adoption, or nullptr. */
    uint32_t gen; /**< Build generation the effect was constructed for. */
  };

  /**
   * @brief What a caller must still do once apply_wake() has run the sequence.
   */
  struct Wake {
    T *live = nullptr; /**< Effect live after this wake's adopt, or nullptr. */
    bool commit_ok = true; /**< False when a commit deadline found nothing
                                committable: the caller must trap. */
    bool adopted = false;  /**< An effect was taken live this wake. */
    bool advance = false;  /**< Call advance_display() on `live`. */
    bool dark = false;     /**< Submit-gate dark input. */
  };

  // ── Foreground side ────────────────────────────────────────────────────

  /**
   * @brief Requests the ISR drop its live pointer before the instance is freed.
   * @details Bumps the request counter; the ISR acknowledges within one wake-up.
   */
  void request_release() {
    release_req.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @brief Foreground wait predicate for the teardown handshake.
   * @return True once the ISR has acknowledged every outstanding release.
   */
  bool release_complete() const {
    return release_ack.load(std::memory_order_relaxed) ==
           release_req.load(std::memory_order_relaxed);
  }

  /**
   * @brief Clears the pending slot before the old instance is destroyed.
   */
  void clear_pending() {
    pending_effect.store(nullptr, std::memory_order_release);
  }

  /**
   * @brief Publishes a constructed (effect, generation) pair to the ISR.
   * @param effect Fully-constructed, first-frame-drawn instance.
   * @param gen Build generation the effect was constructed for.
   * @details The device brackets this call in interrupts-off so the pair
   *          publishes atomically w.r.t. the ISR; the release store still orders
   *          the constructor writes ahead of the ISR's acquire load on its own.
   */
  void publish(T *effect, uint32_t gen) {
    pending_gen.store(gen, std::memory_order_relaxed);
    pending_effect.store(effect, std::memory_order_release);
  }

  /**
   * @brief Whether the ISR has taken a given build generation live.
   * @param gen Build generation to test.
   * @return True if consumed_gen equals @p gen.
   */
  bool consumed(uint32_t gen) const {
    return consumed_gen.load(std::memory_order_relaxed) == gen;
  }

  // ── ISR side ───────────────────────────────────────────────────────────

  /**
   * @brief Services a pending teardown request: drop the live pointer and ack.
   * @details Idempotent — a no-op when no release is outstanding.
   */
  void service_release() {
    if (release_ack.load(std::memory_order_relaxed) !=
        release_req.load(std::memory_order_relaxed)) {
      live_effect = nullptr;
      release_ack.store(release_req.load(std::memory_order_relaxed),
                        std::memory_order_relaxed);
    }
  }

  /**
   * @brief The effect the ISR currently renders (ISR-owned).
   * @return Live effect pointer, or nullptr.
   */
  T *live() const { return live_effect; }

  /**
   * @brief Acquire-loads the pending (effect, generation) pair.
   * @return The pending snapshot.
   * @details The acquire load pairs with publish()'s release store, ordering the
   *          constructed instance's member writes before any dereference.
   */
  Pending pending_acquire() const {
    T *e = pending_effect.load(std::memory_order_acquire);
    return {e, pending_gen.load(std::memory_order_relaxed)};
  }

  /**
   * @brief Whether a pending snapshot may be committed at the epoch deadline.
   * @param p A snapshot from pending_acquire().
   * @param wire_gen Generation advertised by the sync wire this tick.
   * @return True iff the effect is present and its generation matches the wire.
   * @details The caller's HS_CHECK on this is the use-after-free guard: a false
   *          result at a commit boundary means effect init overran its window.
   */
  bool committable(const Pending &p, uint32_t wire_gen) const {
    return p.effect != nullptr && p.gen == wire_gen;
  }

  /**
   * @brief Whether a pending snapshot may be adopted at a join boundary.
   * @param p A snapshot from pending_acquire().
   * @param wire_gen Generation advertised by the sync wire this tick.
   * @return True iff the effect is present, not already consumed, and still
   *         matches the wire's advertised generation.
   * @details A visibility lag that fails the wire match simply joins one grid
   *          step later — join is conditional where commit traps.
   */
  bool joinable(const Pending &p, uint32_t wire_gen) const {
    return p.effect != nullptr &&
           p.gen != consumed_gen.load(std::memory_order_relaxed) &&
           p.gen == wire_gen;
  }

  /**
   * @brief Takes a pending effect live and records its generation as consumed.
   * @param effect Instance to render.
   * @param gen Its build generation.
   */
  void adopt(T *effect, uint32_t gen) {
    live_effect = effect;
    consumed_gen.store(gen, std::memory_order_relaxed);
  }

  /**
   * @brief Publishes which half the now-open display window sweeps.
   * @param zero_crossing True if the flip that opened the window was a ZERO
   *        crossing (opens the arm-A-left [0,W/2) half-rev).
   */
  void set_window_left(bool zero_crossing) {
    window_sweeps_left.store(zero_crossing ? 1u : 0u,
                             std::memory_order_relaxed);
  }

  /**
   * @brief The half the open display window sweeps.
   * @return 1 when the window sweeps arm-A columns [0,W/2), else 0.
   */
  uint8_t window_left() const {
    return window_sweeps_left.load(std::memory_order_relaxed);
  }

  /**
   * @brief Runs the whole ISR-side handoff sequence for one flywheel wake.
   * @param in This wake's sync-layer decisions and the wire's build generation.
   * @return The post-adopt live pointer plus the caller's remaining duties.
   * @details The step order is the contract: the teardown handshake runs before
   * the adopt, so a same-wake adopt cannot have its pointer dropped from under
   * it and the join arm's emptiness test sees the released state; the window
   * half publishes before the caller's advance_display() unblocks a foreground
   * waiting in buffer_free(); and `live` is read after the adopt, so an effect
   * taken live this wake also flips this wake. Commit reports commit_ok=false
   * for the caller to trap on, where join is conditional — a failed join gate
   * waits for the next grid step. always_inline: the flywheel ISR must not pay
   * a call for it.
   */
  __attribute__((always_inline)) Wake apply_wake(const WakeInputs &in) {
    service_release();
    Wake out;
    if (in.commit) {
      const Pending p = pending_acquire();
      out.commit_ok = committable(p, in.wire_gen);
      if (out.commit_ok) {
        adopt(p.effect, p.gen);
        out.adopted = true;
      }
    } else if (in.join_boundary && !in.dark && live_effect == nullptr) {
      const Pending p = pending_acquire();
      if (joinable(p, in.wire_gen)) {
        adopt(p.effect, p.gen);
        out.adopted = true;
      }
    }
    if (in.flip)
      set_window_left(in.zero_crossing);
    out.live = live_effect;
    out.advance = in.flip && out.live != nullptr;
    out.dark = in.dark || out.live == nullptr;
    return out;
  }

private:
  T *live_effect = nullptr; /**< Effect the ISR renders; ISR-owned. */
  std::atomic<T *> pending_effect{
      nullptr}; /**< Next effect; fg-written (release), ISR-read (acquire). */
  std::atomic<uint32_t> pending_gen{
      0}; /**< Build generation of pending_effect; fg-written. */
  std::atomic<uint32_t> consumed_gen{
      0}; /**< Build generation taken live; ISR-written. */
  std::atomic<uint32_t> release_req{
      0}; /**< Teardown request counter; fg-written. */
  std::atomic<uint32_t> release_ack{
      0}; /**< Teardown acknowledge counter; ISR-written. */
  std::atomic<uint8_t> window_sweeps_left{
      1}; /**< ISR-written: 1 when the open window sweeps arm-A [0,W/2). */
};

} // namespace pov
