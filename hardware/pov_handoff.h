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
 *   - publish() stores pending_gen_ (relaxed) then pending_effect_ (release);
 *     the release store publishes the (effect, gen) pair and orders every
 *     constructor/draw write before the ISR's acquire load can observe it.
 *   - pending_acquire() loads pending_effect_ (acquire) then pending_gen_
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
 * @brief Foreground↔ISR effect-handoff state machine.
 * @tparam T Effect pointee type (the device instantiates it with Effect; tests
 *           use a stand-in).
 *
 * Single-writer per field (spec §8): live_ is ISR-written; pending_* are
 * foreground-written (published under a brief interrupts-off bracket on the
 * device); the release_req_/release_ack_ pair is the teardown handshake
 * (foreground bumps req, ISR acks); window_left_ is ISR-written, foreground-read.
 */
template <class T>
class EffectHandoff {
public:
  /**
   * @brief An acquire-loaded (effect, generation) snapshot of the pending slot.
   */
  struct Pending {
    T *effect;    /**< Constructed instance awaiting adoption, or nullptr. */
    uint32_t gen; /**< Build generation the effect was constructed for. */
  };

  // ── Foreground side ────────────────────────────────────────────────────

  /**
   * @brief Requests the ISR drop its live pointer before the instance is freed.
   * @details Bumps the request counter; the ISR acknowledges within one wake-up.
   */
  void request_release() {
    release_req_.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @brief Foreground wait predicate for the teardown handshake.
   * @return True once the ISR has acknowledged every outstanding release.
   */
  bool release_complete() const {
    return release_ack_.load(std::memory_order_relaxed) ==
           release_req_.load(std::memory_order_relaxed);
  }

  /**
   * @brief Clears the pending slot before the old instance is destroyed.
   */
  void clear_pending() {
    pending_effect_.store(nullptr, std::memory_order_release);
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
    pending_gen_.store(gen, std::memory_order_relaxed);
    pending_effect_.store(effect, std::memory_order_release);
  }

  /**
   * @brief Whether the ISR has taken a given build generation live.
   * @param gen Build generation to test.
   * @return True if consumed_gen_ equals @p gen.
   */
  bool consumed(uint32_t gen) const {
    return consumed_gen_.load(std::memory_order_relaxed) == gen;
  }

  // ── ISR side ───────────────────────────────────────────────────────────

  /**
   * @brief Services a pending teardown request: drop the live pointer and ack.
   * @details Idempotent — a no-op when no release is outstanding.
   */
  void service_release() {
    if (release_ack_.load(std::memory_order_relaxed) !=
        release_req_.load(std::memory_order_relaxed)) {
      live_ = nullptr;
      release_ack_.store(release_req_.load(std::memory_order_relaxed),
                         std::memory_order_relaxed);
    }
  }

  /**
   * @brief The effect the ISR currently renders (ISR-owned).
   * @return Live effect pointer, or nullptr.
   */
  T *live() const { return live_; }

  /**
   * @brief Acquire-loads the pending (effect, generation) pair.
   * @return The pending snapshot.
   * @details The acquire load pairs with publish()'s release store, ordering the
   *          constructed instance's member writes before any dereference.
   */
  Pending pending_acquire() const {
    T *e = pending_effect_.load(std::memory_order_acquire);
    return {e, pending_gen_.load(std::memory_order_relaxed)};
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
           p.gen != consumed_gen_.load(std::memory_order_relaxed) &&
           p.gen == wire_gen;
  }

  /**
   * @brief Takes a pending effect live and records its generation as consumed.
   * @param effect Instance to render.
   * @param gen Its build generation.
   */
  void adopt(T *effect, uint32_t gen) {
    live_ = effect;
    consumed_gen_.store(gen, std::memory_order_relaxed);
  }

  /**
   * @brief Publishes which half the now-open display window sweeps.
   * @param zero_crossing True if the flip that opened the window was a ZERO
   *        crossing (opens the arm-A-left [0,W/2) half-rev).
   */
  void set_window_left(bool zero_crossing) {
    window_left_.store(zero_crossing ? 1u : 0u, std::memory_order_relaxed);
  }

  /**
   * @brief The half the open display window sweeps.
   * @return 1 when the window sweeps arm-A columns [0,W/2), else 0.
   */
  uint8_t window_left() const {
    return window_left_.load(std::memory_order_relaxed);
  }

  /**
   * @brief The generation the ISR most recently took live.
   * @return consumed_gen_.
   */
  uint32_t consumed_gen() const {
    return consumed_gen_.load(std::memory_order_relaxed);
  }

private:
  T *live_ = nullptr; /**< Effect the ISR renders; ISR-owned. */
  std::atomic<T *> pending_effect_{nullptr}; /**< Next effect; fg-written (release), ISR-read (acquire). */
  std::atomic<uint32_t> pending_gen_{0};     /**< Build generation of pending_effect_; fg-written. */
  std::atomic<uint32_t> consumed_gen_{0};    /**< Build generation taken live; ISR-written. */
  std::atomic<uint32_t> release_req_{0};     /**< Teardown request counter; fg-written. */
  std::atomic<uint32_t> release_ack_{0};     /**< Teardown acknowledge counter; ISR-written. */
  std::atomic<uint8_t> window_left_{1};      /**< ISR-written: 1 when the open window sweeps arm-A [0,W/2). */
};

} // namespace pov
