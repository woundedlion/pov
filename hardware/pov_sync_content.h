/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file pov_sync_content.h
 * @brief Layer 3 of the Phantasm sync design (spec §6): the index-beacon
 *        codec and the per-board content tracker that the epoch and beacon
 *        paths both drive.
 */
#pragma once

#include "pov_sync_protocol.h"

#include <cstdint>

namespace pov {
namespace sync {

// ── Layer 3: beacon codec (spec §6.4) ───────────────────────────────────────

/**
 * @brief Decoded index beacon: absolute effect index + revolution count
 * (mod 64).
 */
struct BeaconFrame {
  int32_t effect_index = 0; /**< Effect index, 0–63. */
  uint32_t rev_count = 0;   /**< Revolution within the effect, mod 64. */
};

/**
 * @brief Encodes a beacon as five base-8 digits.
 * @param effect_index Effect index (low 6 bits used).
 * @param rev_count Revolution count (low 6 bits used).
 * @param out Out: five digits [index_hi, index_lo, rev_hi, rev_lo, checksum].
 * @details Each digit is transmitted as a burst of (digit+1) pulses.
 */
constexpr void encode_beacon_digits(int32_t effect_index, uint32_t rev_count,
                                    uint8_t out[5]) {
  const uint32_t idx = static_cast<uint32_t>(effect_index) & 63u;
  const uint32_t rev = rev_count & 63u;
  out[0] = static_cast<uint8_t>(idx >> 3);
  out[1] = static_cast<uint8_t>(idx & 7u);
  out[2] = static_cast<uint8_t>(rev >> 3);
  out[3] = static_cast<uint8_t>(rev & 7u);
  // Position-weighted Σ(i+1)·dᵢ mod 8: catches digit transpositions and
  // compensating miscounts a plain sum would miss.
  out[4] = static_cast<uint8_t>(
      (1u * out[0] + 2u * out[1] + 3u * out[2] + 4u * out[3]) & 7u);
}

/**
 * @brief Smallest signed revolution slip that maps @p current's 6-bit residue
 * onto a beacon's @p beacon_rev_count (spec §6.4 rev cross-check).
 * @param beacon_rev_count The beacon's rev_count digit pair (0–63).
 * @param current The board's current rev_in_effect (any magnitude).
 * @return A delta in [-32, 31]. The `+96 %64 -32` fold resolves the 63↔0 mod-64
 * seam, so a small slip across the wrap reads as a small *signed* delta (e.g.
 * current residue 62 vs beacon 0 → +2) rather than a ~64-rev jump. Slips of
 * magnitude ≥ 32 are aliased by the mod-64 (the protocol guarantees agreement
 * within 32 revs), so the caller treats them as the nearest-residue correction.
 */
constexpr int32_t beacon_rev_resync_delta(uint32_t beacon_rev_count,
                                          uint32_t current) {
  return ((static_cast<int32_t>(beacon_rev_count & 63u) -
           static_cast<int32_t>(current & 63u) + 96) %
          64) -
         32;
}

/**
 * @brief Assembles data bursts into beacon frames. Integrity by rejection
 * (spec §6.4): any out-of-range digit, stale partial frame, or checksum
 * mismatch drops the whole frame — the next beacon is one cadence away.
 */
class BeaconParser {
public:
  /**
   * @brief Feed one data burst into the frame assembler.
   * @param s The data burst.
   * @param cfg Protocol configuration.
   * @param out Out: the parsed frame, valid only when the return is true.
   * @param rejected Out: set when a frame (not merely a digit) was discarded.
   * @return True when a complete, checksum-valid frame was parsed into @p out.
   */
  bool feed(const BurstSnapshot &s, const Config &cfg, BeaconFrame *out,
            bool *rejected) {
    *rejected = false;
    if (n > 0 && (s.first_cycles - last_first_cycles) >
                     cfg.interdigit_timeout_cycles()) {
      n = 0; // stale partial frame: a new frame may start with this burst
      *rejected = true;
    }
    if (s.count < 1 || s.count > 8) {
      if (n > 0)
        *rejected = true;
      n = 0;
      return false;
    }
    last_first_cycles = s.first_cycles;
    digits[n++] = static_cast<uint8_t>(s.count - 1);
    if (n < 5)
      return false;
    n = 0;
    // Position-weighted checksum (see encode_beacon_digits).
    if (((1u * digits[0] + 2u * digits[1] + 3u * digits[2] + 4u * digits[3]) &
         7u) != digits[4]) {
      *rejected = true;
      return false;
    }
    out->effect_index = digits[0] * 8 + digits[1];
    out->rev_count = static_cast<uint32_t>(digits[2]) * 8u + digits[3];
    return true;
  }

  /**
   * @brief Discards any partial frame in progress.
   */
  void reset() { n = 0; }
  /**
   * @brief Whether a partial frame is being assembled.
   * @return True if at least one digit has been buffered.
   */
  bool active() const { return n > 0; }

private:
  int32_t n = 0;
  uint8_t digits[5] = {};
  uint32_t last_first_cycles = 0;
};

// ── Layer 3: content tracker (spec §6.1, §6.3) ──────────────────────────────

/**
 * @brief Per-board content state: which effect, which revolution, and the
 * deadline-scheduled epoch commit. Owned by the flywheel ISR via SyncBoard.
 */
struct ContentTracker {
  bool identity_known =
      false;                /**< False until an epoch/beacon names the index. */
  int32_t effect_index = 0; /**< Currently displayed effect index. */
  /**
   * @brief ZERO crossings since effect start.
   * @details For a beacon-joined board this starts from the beacon's mod-64
   * value — it only feeds the master's schedule and cross-check telemetry.
   */
  uint32_t rev_in_effect = 0;
  bool commit_pending = false;       /**< An epoch commit is scheduled. */
  uint32_t commit_in_revs = 0;       /**< ZERO crossings remaining until the
                                        absolute B+R+K boundary, with j already
                                        subtracted for the repeat this board
                                        heard (NOT the announce-phase length). */
  uint32_t refractory_revs_left = 0; /**< EPOCH dedup window (spec §6.1). */

  /**
   * @brief Accepts an EPOCH symbol and may open a commit window.
   * @param cfg Protocol configuration.
   * @return True if it opened a commit window (false inside the refractory
   * window — the §6.3 redundancy repeats land here).
   * @details The commit boundary is ABSOLUTE: B + R + K, where B is the primary
   * copy's boundary (R = epoch_repeats, K = commit_revs). Which copy of the
   * train this is (j, 0 = primary) is inferred from the shared revolution
   * count — the master starts the train exactly when rev_in_effect reaches
   * the active effect's configured duration, and by the time a symbol is
   * consumed the local crossing for its boundary has already incremented
   * rev_in_effect (classification
   * completes ~13 columns after the boundary instant). So every board that
   * hears ANY copy counts down to the same boundary, and hearing a repeat
   * instead of the primary cannot skew the commit (§6.3.1). A board whose
   * revolution count is not absolute (it beacon-joined mid-effect, §6.4)
   * lands outside the train window and falls back to j = 0 — it commits up
   * to j revolutions late, an epoch-bounded degradation confined to that
   * case. The resulting counter slip is then resynced from the next beacon's
   * rev cross-check (handle_beacon_burst), so every subsequent epoch is
   * lockstep.
   */
  bool on_epoch_symbol(const Config &cfg) {
    if (refractory_revs_left > 0)
      return false;
    const uint32_t effect_revolutions =
        cfg.revolutions_for_effect(effect_index);
    uint32_t j = 0;
    if (rev_in_effect >= effect_revolutions &&
        rev_in_effect - effect_revolutions <=
            static_cast<uint32_t>(cfg.epoch_repeats))
      j = rev_in_effect - effect_revolutions;
    commit_pending = true;
    commit_in_revs =
        cfg.commit_revs + static_cast<uint32_t>(cfg.epoch_repeats) - j;
    refractory_revs_left = cfg.refractory_revs;
    return true;
  }

  /**
   * @brief Whether the construction window opens at this crossing/accept.
   * @param cfg Protocol configuration.
   * @return True at the single crossing/accept where the construction window —
   * the last K revolutions before the commit boundary — opens.
   * @details commit_in_revs strictly decreases, so == fires exactly once per
   * window. The announce phase before it is what gives a board that hears only
   * the last repeat the full K-revolution construction budget; the outgoing
   * effect keeps playing through it, so the dark window stays K revolutions.
   */
  bool construction_opens(const Config &cfg) const {
    return commit_pending && commit_in_revs == cfg.commit_revs;
  }

  /**
   * @brief Whether the board is in the construction window.
   * @param cfg Protocol configuration.
   * @return True throughout the construction window (fail-dark, spec §6.1).
   */
  bool constructing(const Config &cfg) const {
    return commit_pending && commit_in_revs <= cfg.commit_revs;
  }

  /**
   * @brief Advances revolution bookkeeping on a ZERO boundary flip.
   * @param cfg Protocol configuration.
   * @return True when the commit deadline is reached — the board must swap to
   * the (already constructed) next effect at exactly this boundary.
   */
  bool on_zero_crossing(const Config &cfg) {
    ++rev_in_effect;
    if (refractory_revs_left > 0)
      --refractory_revs_left;
    if (commit_pending && --commit_in_revs == 0) {
      commit_pending = false;
      effect_index = (effect_index + 1) % cfg.effect_count;
      rev_in_effect = 0;
      return true;
    }
    return false;
  }
};

} // namespace sync
} // namespace pov
