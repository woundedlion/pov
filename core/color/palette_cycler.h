/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

#include "color/generative_palette.h"
#include "color/baked_palette.h"

/**
 * @file palette_cycler.h
 * @brief Time-based palette cycling through a display LUT.
 */

/**
 * @brief Cycles a display LUT through a sequence of palettes over time.
 * @details Adjacent morph-compatible GenerativePalettes fade by key-space
 * morph, other pairs by baked-LUT crossfade. Entries are caller-owned, must
 * outlive the cycler, and must remain unchanged between init() calls. Outside a
 * fade, after step(), the display LUT is a bit-exact bake of the current entry;
 * after advance_without_display() it holds the last displayed bake until the
 * next step().
 */
class PaletteCycler {
public:
  static constexpr int MAX_ENTRIES = 8; ///< Roster capacity.

  PaletteCycler() = default;
  // The display and fade LUTs and the morph slots are arena handles, so a copy
  // would drive the original's bakes.
  PaletteCycler(const PaletteCycler &) = delete;
  PaletteCycler &operator=(const PaletteCycler &) = delete;

  /**
   * @brief Tagged reference to one palette of any kind.
   */
  struct Entry {
    const Palette *source = nullptr; /**< Sampled source; null for prebaked. */
    const GenerativePalette *generative =
        nullptr; /**< Set when the source is generative. */
    const BakedPalette *baked = nullptr; /**< Prebaked LUT; null otherwise. */

    Entry() = default;
    /**
     * @brief References a generative palette, eligible for key morphs.
     * @param palette Borrowed palette; must outlive the cycler.
     */
    Entry(const GenerativePalette &palette)
        : source(&palette), generative(&palette) {}
    /**
     * @brief References a sampled palette source.
     * @param palette Borrowed palette; must outlive the cycler.
     */
    Entry(const Palette &palette) : source(&palette) {}
    /**
     * @brief References a prebaked LUT.
     * @param palette Borrowed table; must outlive the cycler.
     */
    Entry(const BakedPalette &palette) : baked(&palette) {}

    // The referenced palette must outlive the cycler; reject temporaries.
    Entry(const GenerativePalette &&) = delete;
    Entry(const Palette &&) = delete;
    Entry(const BakedPalette &&) = delete;
  };

  /**
   * @brief Conservative arena byte budget for the display LUT allocated by
   * init().
   * @return Byte budget.
   */
  static constexpr size_t display_arena_bytes() {
    return BakedPalette::required_arena_bytes();
  }

  /** @brief Conservative extra arena byte budget when any adjacent pair fades by
   *  LUT crossfade rather than key morph.
   *  @return Byte budget.
   */
  static constexpr size_t crossfade_arena_bytes() {
    return 2 * BakedPalette::required_arena_bytes();
  }

  /**
   * @brief Conservative extra arena byte budget when any adjacent pair
   * key-morphs.
   * @return Byte budget.
   */
  static constexpr size_t morph_arena_bytes() {
    return sizeof(GenerativePalette) + alignof(GenerativePalette);
  }

  /** @brief Conservative arena byte budget for init(). @return Byte budget. */
  static constexpr size_t required_arena_bytes() {
    return display_arena_bytes() + crossfade_arena_bytes() +
           morph_arena_bytes();
  }

  /**
   * @brief Conservative arena byte budget for init_generated().
   * @return Byte budget.
   */
  static constexpr size_t generated_arena_bytes() {
    return display_arena_bytes() + 3 * morph_arena_bytes();
  }

  /**
   * @brief Fills @p out with palette number @p sequence of a generated cycle.
   * @details Successive palettes must be morph-compatible; the cycler
   * fail-fast checks each retarget.
   */
  using NextPaletteFn = void (*)(void *context, uint32_t sequence,
                                 GenerativePalette &out);

  /**
   * @brief Allocates the display LUT and enters the first entry's dwell.
   * @param arena Arena for the display, two endpoint scratch LUTs if any pair
   * needs a LUT crossfade, and a morph scratch palette if any pair key-morphs.
   * @param entry_list Caller-owned entry array; must outlive the cycler.
   * @param count Number of entries in [1, MAX_ENTRIES]; 1 shows a static
   * palette.
   * @param dwell_frames Hold steps before fading, with a minimum of one step;
   * values 0 and 1 both begin the next fade on the next step.
   * @param fade_frames Frames each fade spans, >= 1.
   * @param easing_fn Optional easing over fade progress; null is linear.
   * @param paused_flag Optional pause gate; freezes step() while set and true.
   */
  HS_COLD_MEMBER void init(Arena &arena, const Entry *entry_list, int count,
                           int dwell_frames, int fade_frames,
                           float (*easing_fn)(float) = nullptr,
                           const bool *paused_flag = nullptr) {
    HS_CHECK(entry_list != nullptr && count >= 1 && count <= MAX_ENTRIES,
             "PaletteCycler entry count must be in [1, MAX_ENTRIES]");
    HS_CHECK(dwell_frames >= 0 && fade_frames >= 1,
             "PaletteCycler dwell must be >= 0 and fade >= 1 frames");
    for (int i = 0; i < count; ++i)
      HS_CHECK(entry_list[i].source != nullptr ||
                   entry_list[i].baked != nullptr,
               "PaletteCycler entry references no palette");
    entries = entry_list;
    entry_count = count;
    provider = nullptr;
    provider_context = nullptr;
    from_slot = nullptr;
    to_slot = nullptr;
    morph = nullptr;
    next_sequence = 0;
    dwell = dwell_frames;
    fade = fade_frames;
    easing = easing_fn;
    paused = paused_flag;
    current = 0;
    frame = 0;
    fade_active = false;
    display_dirty = false;

    key_morph_mask = 0;
    bool crossfades = false;
    if (entry_count > 1) {
      for (int i = 0; i < entry_count; ++i) {
        const Entry &a = entries[i];
        const Entry &b = entries[next_of(i)];
        if (a.generative != nullptr && b.generative != nullptr &&
            a.generative->morph_compatible(*b.generative))
          key_morph_mask |= 1u << i;
        else
          crossfades = true;
      }
    }

    bake_entry(display, arena, entries[0]);
    ++generation;
    if (crossfades) {
      bake_entry(fade_from, arena, entries[0]);
      bake_entry(fade_to, arena, entries[0]);
    }
    if (key_morph_mask != 0)
      morph = allocate_palette(arena);
  }

  /**
   * @brief Allocates the display LUT and enters a provider-generated cycle.
   * @param arena Arena the display, morph scratch, and two palette slots are
   * allocated from.
   * @param next_fn Fills palette number `sequence`; called with 0 and 1 here,
   * then once per completed fade. Successors must be morph-compatible.
   * @param context Opaque state forwarded to @p next_fn; caller-owned.
   * @param dwell_frames Hold steps before fading, with a minimum of one step;
   * values 0 and 1 both begin the next fade on the next step.
   * @param fade_frames Frames each fade spans, >= 1.
   * @param easing_fn Optional easing over fade progress; null is linear.
   * @param paused_flag Optional pause gate; freezes step() while set and true.
   */
  HS_COLD_MEMBER void init_generated(Arena &arena, NextPaletteFn next_fn,
                                     void *context, int dwell_frames,
                                     int fade_frames,
                                     float (*easing_fn)(float) = nullptr,
                                     const bool *paused_flag = nullptr) {
    HS_CHECK(next_fn != nullptr,
             "PaletteCycler generated cycle needs a provider");
    HS_CHECK(dwell_frames >= 0 && fade_frames >= 1,
             "PaletteCycler dwell must be >= 0 and fade >= 1 frames");
    entries = nullptr;
    entry_count = 0;
    key_morph_mask = 0;
    provider = next_fn;
    provider_context = context;
    dwell = dwell_frames;
    fade = fade_frames;
    easing = easing_fn;
    paused = paused_flag;
    current = 0;
    frame = 0;
    fade_active = false;
    display_dirty = false;
    next_sequence = 2;

    from_slot = allocate_palette(arena);
    to_slot = allocate_palette(arena);
    morph = allocate_palette(arena);
    provider(provider_context, 0, *from_slot);
    provider(provider_context, 1, *to_slot);
    HS_CHECK(from_slot->morph_compatible(*to_slot),
             "PaletteCycler generated palettes must be morph-compatible");
    display.bake(arena, *from_slot);
    ++generation;
  }

  /**
   * @brief Advances the cycle by one frame.
   * @details The cycle clock freezes while a wired pause flag is set. A dirty
   * display still rebuilds at its held phase. The step completing a
   * fade rebakes the display directly from the target entry, so the landing
   * is bit-exact regardless of easing endpoint behavior.
   */
  HS_COLD_MEMBER void step() {
    if (!advance_clock(true))
      return;
    const float w = fade_weight();
    if (provider != nullptr) {
      morph->morph_palettes(*from_slot, *to_slot, w);
      rebake_display(*morph);
    } else if ((key_morph_mask & (1u << current)) != 0) {
      HS_AUDIT_CHECK(
          entries[current].generative->morph_compatible(
              *entries[next_of(current)].generative),
          "PaletteCycler borrowed palette policy changed after init");
      morph->morph_palettes(*entries[current].generative,
                            *entries[next_of(current)].generative, w);
      rebake_display(*morph);
    } else {
      display.rebake_crossfade(fade_from, fade_to, w);
      ++generation;
    }
    display_dirty = false;
  }

  /** @brief Advances the cycle without rebuilding the display LUT.
   *  @details Timeline and provider state still land exactly. A later step()
   *  immediately rebuilds the display at the current phase. */
  HS_COLD_MEMBER void advance_without_display() {
    if (advance_clock(false))
      display_dirty = true;
  }

  /**
   * @brief The display LUT effects shade from.
   * @return The display LUT.
   */
  const BakedPalette &palette() const { return display; }

  /** @brief Counter incremented on every display-LUT bake; a cache key for
   *  tables derived from palette().
   *  @return Bake count.
   */
  uint32_t bake_generation() const { return generation; }

  /** @brief Sets generated endpoint chroma and rebakes the current display.
   *  @param chroma Gamut-relative chroma in [0, 1].
   *  @details Valid after init_generated(). Skips the rebake while a display
   *  rebuild is already pending, which the next step() serves at the new
   *  chroma. */
  HS_COLD_MEMBER void set_generated_chroma(float chroma) {
    HS_CHECK(from_slot != nullptr && to_slot != nullptr && morph != nullptr,
             "PaletteCycler chroma needs a generated cycle");
    from_slot->set_constant_chroma(chroma);
    to_slot->set_constant_chroma(chroma);
    if (display_dirty)
      return;
    if (!fade_active) {
      rebake_display(*from_slot);
      return;
    }
    const float weight = fade_weight();
    morph->morph_palettes(*from_slot, *to_slot, weight);
    rebake_display(*morph);
  }

  /** @brief Index of the roster entry currently dwelt on or faded away from.
   *  @details Always 0 for a generated cycle (init_generated()).
   *  @return Roster index.
   */
  int current_index() const { return current; }

  /**
   * @brief True while a fade toward the next entry is in flight.
   * @return Whether a fade is active.
   */
  bool fading() const { return fade_active; }

  /** @brief Serializable timeline state of a generated cycle. */
  struct GeneratedClock {
    uint32_t frame = 0;         ///< Frame within the current dwell or fade.
    uint32_t next_sequence = 2; ///< Sequence of the next provider call; >= 2.
    bool fade_active = false;   ///< Whether a fade is in flight.
    bool display_dirty = false; ///< Whether the display LUT is stale.
  };

  /**
   * @brief Captures the generated-cycle timeline.
   * @return The current clock.
   */
  GeneratedClock generated_clock() const {
    return {static_cast<uint32_t>(frame), next_sequence, fade_active,
            display_dirty};
  }

  /**
   * @brief Whether @p clock is a reachable state for the given timing.
   * @param clock Clock to check.
   * @param fade_frames Frames each fade spans.
   * @param dwell_frames Hold frames between fades; minimum 1 is applied.
   * @return True when `next_sequence` >= 2 and `frame` is inside its phase.
   */
  static bool valid_generated_clock(const GeneratedClock &clock,
                                    uint32_t fade_frames,
                                    uint32_t dwell_frames) {
    return clock.next_sequence >= 2 &&
           clock.frame < (clock.fade_active
                              ? fade_frames
                              : std::max(uint32_t{1}, dwell_frames));
  }

  /**
   * @brief Restores a generated cycle's timeline and endpoints, and rebakes
   * the display LUT.
   * @param clock Timeline to restore; check with valid_generated_clock().
   * @param from Fade w = 0 endpoint.
   * @param to Fade w = 1 endpoint.
   * @pre init_generated() has run.
   */
  HS_COLD_MEMBER void restore_generated(const GeneratedClock &clock,
                                        const GenerativePalette &from,
                                        const GenerativePalette &to) {
    HS_CHECK(from_slot != nullptr && to_slot != nullptr && morph != nullptr,
             "PaletteCycler restore needs a generated cycle");
    *from_slot = from;
    *to_slot = to;
    frame = static_cast<int>(clock.frame);
    next_sequence = clock.next_sequence;
    fade_active = clock.fade_active;
    display_dirty = clock.display_dirty;
    if (fade_active) {
      const float weight = fade_weight();
      morph->morph_palettes(*from_slot, *to_slot, weight);
      rebake_display(*morph);
    } else {
      rebake_display(*from_slot);
    }
  }

private:
  __attribute__((always_inline)) inline bool
  advance_clock(bool update_display) {
    if (update_display && display_dirty && !fade_active) {
      if (provider != nullptr)
        rebake_display(*from_slot);
      else
        rebake_display_entry(entries[current]);
      display_dirty = false;
    }
    if ((paused != nullptr && *paused) ||
        (provider == nullptr && entry_count < 2))
      return update_display && display_dirty && fade_active;
    ++frame;
    if (!fade_active) {
      if (frame >= dwell) {
        if (provider == nullptr)
          begin_fade();
        fade_active = true;
        frame = 0;
      }
      return false;
    }
    if (frame >= fade) {
      finish_fade(update_display);
      return false;
    }
    return true;
  }

  __attribute__((always_inline)) inline float fade_weight() const {
    const float progress = static_cast<float>(frame) / static_cast<float>(fade);
    return easing != nullptr ? easing(progress) : progress;
  }

  int next_of(int index) const {
    return index + 1 == entry_count ? 0 : index + 1;
  }

  HS_COLD_MEMBER void rebake_display(const GenerativePalette &source) {
    display.rebake(source);
    ++generation;
  }

  HS_COLD_MEMBER void rebake_display_entry(const Entry &entry) {
    rebake_entry(display, entry);
    ++generation;
  }

  HS_COLD_MEMBER static GenerativePalette *allocate_palette(Arena &arena) {
    return arena.make<GenerativePalette>();
  }

  HS_COLD_MEMBER void finish_fade(bool update_display = true) {
    if (provider != nullptr) {
      if (update_display)
        rebake_display(*to_slot);
      GenerativePalette *retired = from_slot;
      from_slot = to_slot;
      to_slot = retired;
      provider(provider_context, next_sequence++, *to_slot);
      HS_CHECK(from_slot->morph_compatible(*to_slot),
               "PaletteCycler generated palette breaks morph compatibility");
    } else {
      if (update_display)
        rebake_display_entry(entries[next_of(current)]);
      current = next_of(current);
    }
    display_dirty = !update_display;
    fade_active = false;
    frame = 0;
  }

  HS_COLD_MEMBER static void bake_entry(BakedPaletteStorage &lut, Arena &arena,
                                        const Entry &entry) {
    if (entry.baked != nullptr)
      lut.clone_from(*entry.baked, arena);
    else if (entry.generative != nullptr)
      lut.bake(arena, *entry.generative);
    else
      lut.bake(arena, *entry.source);
  }

  // Generative entries bake through their concrete type: the mirror and loop
  // domain shortcuts, and a loop's exact seam entry, need it.
  HS_COLD_MEMBER static void rebake_entry(BakedPaletteStorage &lut,
                                          const Entry &entry) {
    if (entry.baked != nullptr)
      lut.rebake_copy(*entry.baked);
    else if (entry.generative != nullptr)
      lut.rebake(*entry.generative);
    else
      lut.rebake(*entry.source);
  }

  HS_COLD_MEMBER void begin_fade() {
    if ((key_morph_mask & (1u << current)) != 0)
      return;
    rebake_entry(fade_from, entries[current]);
    rebake_entry(fade_to, entries[next_of(current)]);
  }

  const Entry *entries = nullptr; /**< Caller-owned entry array. */
  BakedPaletteStorage display; /**< The LUT handed to shaders via palette(). */
  BakedPaletteStorage fade_from; /**< Crossfade w = 0 endpoint scratch. */
  BakedPaletteStorage fade_to;   /**< Crossfade w = 1 endpoint scratch. */
  /** @brief Arena-owned key-morph scratch rebaked into the display; allocated
   *  by init() only when some pair key-morphs. */
  GenerativePalette *morph = nullptr;
  NextPaletteFn provider = nullptr; /**< Generated-cycle palette source. */
  void *provider_context = nullptr; /**< Caller state handed to the provider. */
  GenerativePalette *from_slot = nullptr; /**< Generated fade w = 0 endpoint. */
  GenerativePalette *to_slot = nullptr;   /**< Generated fade w = 1 endpoint. */
  uint32_t next_sequence = 0; /**< Sequence number of the next provider call. */
  uint32_t generation = 0;    /**< Display-LUT bakes performed so far. */
  float (*easing)(float) = nullptr; /**< Fade easing; null = linear. */
  const bool *paused = nullptr; /**< Optional pause gate; null = always runs. */
  int entry_count = 0;
  int dwell = 0;   /**< Hold steps between fades; effective minimum 1. */
  int fade = 0;    /**< Frames each fade spans. */
  int frame = 0;   /**< Frame counter within the current phase. */
  int current = 0; /**< Entry dwelt on or faded away from. */
  uint8_t key_morph_mask = 0; /**< Bit i: entry i fades to its successor by
                                 key morph rather than LUT crossfade. */
  static_assert(MAX_ENTRIES <= 8 * static_cast<int>(sizeof(key_morph_mask)),
                "key_morph_mask needs one bit per entry; widen it alongside "
                "MAX_ENTRIES");
  bool fade_active = false;
  bool display_dirty = false;
};

/** @brief Shared triadic, complementary, and analogous palette cyclers. */
class GeneratedPaletteBank {
public:
  /// Base-hue advance per generated palette, in 1/256-turn wheel steps.
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr int DWELL_FRAMES = 0;  ///< Hold frames between fades.
  static constexpr int FADE_FRAMES = 600; ///< Frames each fade spans.

  /** @brief Serializable state of all three cycles. */
  struct Snapshot {
    float chroma = 0.0f; ///< Shared chroma control in [0, 1].
    /// Current target hue per cycle (triadic, complementary, analogous),
    /// wheel steps.
    std::array<uint32_t, 3> hues{};
    /// Timeline per cycle, same order as `hues`.
    std::array<PaletteCycler::GeneratedClock, 3> cycles{};
  };

  /**
   * @brief Captures the bank's state.
   * @return Chroma, hues and clocks of all three cycles.
   */
  Snapshot snapshot() const {
    return {chroma,
            {triadic_hue, complementary_hue, analogous_hue},
            {triadic.generated_clock(), complementary.generated_clock(),
             analogous.generated_clock()}};
  }

  /**
   * @brief Whether @p snapshot is a state this bank can reach.
   * @param snapshot Snapshot to check.
   * @return True when chroma is in [0, 1] and each clock is valid and agrees
   * with its hue.
   */
  static bool valid_snapshot(const Snapshot &snapshot) {
    if (!std::isfinite(snapshot.chroma) || snapshot.chroma < 0.0f ||
        snapshot.chroma > 1.0f)
      return false;
    for (size_t index = 0; index < snapshot.cycles.size(); ++index) {
      const auto &clock = snapshot.cycles[index];
      if (!PaletteCycler::valid_generated_clock(clock, FADE_FRAMES,
                                                DWELL_FRAMES) ||
          snapshot.hues[index] != (clock.next_sequence - 1) * HUE_STEP)
        return false;
    }
    return true;
  }

  /**
   * @brief Restores all three cycles and rebakes their displays.
   * @param snapshot State to restore; check with valid_snapshot().
   * @pre init() has run.
   */
  HS_COLD_MEMBER void restore_snapshot(const Snapshot &snapshot) {
    chroma = snapshot.chroma;
    triadic_hue = snapshot.hues[0];
    complementary_hue = snapshot.hues[1];
    analogous_hue = snapshot.hues[2];
    const auto restore = [&](PaletteCycler &cycler, size_t index,
                             PaletteHarmony harmony) {
      uint32_t from_hue = snapshot.hues[index] - HUE_STEP;
      uint32_t to_hue = snapshot.hues[index];
      GenerativePalette from, to;
      next_palette(from_hue, 0, harmony, chroma, from);
      next_palette(to_hue, 0, harmony, chroma, to);
      cycler.restore_generated(snapshot.cycles[index], from, to);
    };
    restore(triadic, 0, PaletteHarmony::TRIADIC);
    restore(complementary, 1, PaletteHarmony::COMPLEMENTARY);
    restore(analogous, 2, PaletteHarmony::ANALOGOUS);
  }

  GeneratedPaletteBank() = default;
  // init() hands this to every cycler as its provider context.
  GeneratedPaletteBank(const GeneratedPaletteBank &) = delete;
  GeneratedPaletteBank &operator=(const GeneratedPaletteBank &) = delete;

  /**
   * @brief Conservative arena byte budget for init(), one generated cycler per
   * harmony.
   * @return Byte budget.
   */
  static constexpr size_t required_arena_bytes() {
    return 3 * PaletteCycler::generated_arena_bytes();
  }

  /**
   * @brief Starts all three cycles from hue 0.
   * @param arena Arena for the cyclers' LUTs and slots; see
   *        required_arena_bytes().
   * @param chroma Chroma control in [0, 1].
   * @param easing Fade easing; null = linear.
   */
  HS_COLD_MEMBER void init(Arena &arena, float chroma, float (*easing)(float)) {
    triadic_hue = 0;
    complementary_hue = 0;
    analogous_hue = 0;
    this->chroma = chroma;
    triadic.init_generated(arena, next_triadic, this, DWELL_FRAMES, FADE_FRAMES,
                           easing);
    complementary.init_generated(arena, next_complementary, this, DWELL_FRAMES,
                                 FADE_FRAMES, easing);
    analogous.init_generated(arena, next_analogous, this, DWELL_FRAMES,
                             FADE_FRAMES, easing);
  }

  /**
   * @brief Advances all three cycles one frame; only @p visible rebakes.
   * @tparam PaletteMode Enum with TRIADIC, COMPLEMENTARY and ANALOGOUS.
   * @param visible Cycle whose display LUT is kept current.
   */
  template <typename PaletteMode> void step(PaletteMode visible) {
    step_one(triadic, visible == PaletteMode::TRIADIC);
    step_one(complementary, visible == PaletteMode::COMPLEMENTARY);
    step_one(analogous, visible == PaletteMode::ANALOGOUS);
  }

  /**
   * @brief Display LUT of one cycle.
   * @tparam PaletteMode Enum with TRIADIC, COMPLEMENTARY and ANALOGOUS.
   * @param mode Cycle to read.
   * @return That cycler's display table.
   */
  template <typename PaletteMode>
  HS_COLD_MEMBER const BakedPalette &palette(PaletteMode mode) const {
    switch (mode) {
    case PaletteMode::TRIADIC:
      return triadic.palette();
    case PaletteMode::COMPLEMENTARY:
      return complementary.palette();
    case PaletteMode::ANALOGOUS:
      return analogous.palette();
    }
    HS_CHECK(false, "GeneratedPaletteBank::palette: unknown palette mode");
    return triadic.palette();
  }

  /**
   * @brief Sets the shared chroma on all three cycles.
   * @param chroma Chroma control in [0, 1].
   */
  HS_COLD_MEMBER void set_chroma(float chroma) {
    if (chroma == this->chroma)
      return;
    this->chroma = chroma;
    triadic.set_generated_chroma(chroma);
    complementary.set_generated_chroma(chroma);
    analogous.set_generated_chroma(chroma);
  }

  /**
   * @brief Builds the palette for one step of a cycle.
   * @param hue In/out base hue, wheel steps; advanced by `HUE_STEP` unless
   *        @p sequence is 0.
   * @param sequence Provider sequence number.
   * @param harmony Hue relationship.
   * @param chroma Chroma control in [0, 1].
   * @param out Receives the palette.
   */
  static void next_palette(uint32_t &hue, uint32_t sequence,
                           PaletteHarmony harmony, float chroma,
                           GenerativePalette &out) {
    if (sequence > 0)
      hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, harmony, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(hue), chroma)};
  }

private:
  static void step_one(PaletteCycler &cycler, bool visible) {
    if (visible)
      cycler.step();
    else
      cycler.advance_without_display();
  }

  static void next_triadic(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    auto &bank = *static_cast<GeneratedPaletteBank *>(context);
    next_palette(bank.triadic_hue, sequence, PaletteHarmony::TRIADIC,
                 bank.chroma, out);
  }

  static void next_complementary(void *context, uint32_t sequence,
                                 GenerativePalette &out) {
    auto &bank = *static_cast<GeneratedPaletteBank *>(context);
    next_palette(bank.complementary_hue, sequence,
                 PaletteHarmony::COMPLEMENTARY, bank.chroma, out);
  }

  static void next_analogous(void *context, uint32_t sequence,
                             GenerativePalette &out) {
    auto &bank = *static_cast<GeneratedPaletteBank *>(context);
    next_palette(bank.analogous_hue, sequence, PaletteHarmony::ANALOGOUS,
                 bank.chroma, out);
  }

  PaletteCycler triadic;
  PaletteCycler complementary;
  PaletteCycler analogous;
  uint32_t triadic_hue = 0;
  uint32_t complementary_hue = 0;
  uint32_t analogous_hue = 0;
  float chroma = 0.62f;
};
