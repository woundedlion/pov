/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "color/color.h"

/**
 * @file palette_cycler.h
 * @brief Time-based palette orchestration: a display LUT cycling through an
 *        arbitrary sequence of palettes.
 */

/**
 * @brief Cycles a display LUT through a sequence of palettes over time.
 * @details Dwells on each entry, then fades into the next over a fixed frame
 * count. Adjacent morph-compatible GenerativePalettes fade by key-space morph
 * (perceptually coherent hue travel); every other pair fades by baked-LUT
 * crossfade, so any mix of generative, composed, and prebaked palettes cycles
 * correctly. Entries are caller-owned and must outlive the cycler. Effects
 * call step() once per frame and shade from palette(); outside a fade the
 * display LUT is a bit-exact bake of the current entry.
 */
class PaletteCycler {
public:
  static constexpr int MAX_ENTRIES = 8;

  /**
   * @brief Tagged reference to one palette of any kind.
   */
  struct Entry {
    const Palette *source = nullptr; /**< Sampled source; null for prebaked. */
    const GenerativePalette *generative =
        nullptr; /**< Set when the source is generative. */
    const BakedPalette *baked = nullptr; /**< Prebaked LUT; null otherwise. */

    Entry() = default;
    Entry(const GenerativePalette &palette)
        : source(&palette), generative(&palette) {}
    Entry(const Palette &palette) : source(&palette) {}
    Entry(const BakedPalette &palette) : baked(&palette) {}

    // Borrow contract: the referenced palette is sampled every fade, so it must
    // outlive the cycler; these deleted overloads reject a temporary.
    Entry(const GenerativePalette &&) = delete;
    Entry(const Palette &&) = delete;
    Entry(const BakedPalette &&) = delete;
  };

  /** @brief Arena bytes for the display LUT; always consumed by init(). */
  static constexpr size_t display_arena_bytes() {
    return BakedPalette::required_arena_bytes();
  }

  /** @brief Extra arena bytes consumed when any adjacent pair fades by
   *  LUT crossfade rather than key morph. */
  static constexpr size_t crossfade_arena_bytes() {
    return 2 * BakedPalette::required_arena_bytes();
  }

  /** @brief Extra arena bytes consumed when any adjacent pair key-morphs. */
  static constexpr size_t morph_arena_bytes() {
    return sizeof(GenerativePalette) + alignof(GenerativePalette);
  }

  /** @brief Worst-case arena bytes init() can consume. */
  static constexpr size_t required_arena_bytes() {
    return display_arena_bytes() + crossfade_arena_bytes() +
           morph_arena_bytes();
  }

  /** @brief Arena bytes init_generated() consumes. */
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
   * @param arena Arena the display (and, if any pair needs a LUT crossfade,
   * two endpoint scratch LUTs) is allocated from.
   * @param entry_list Caller-owned entry array; must outlive the cycler.
   * @param count Number of entries in [1, MAX_ENTRIES]; 1 shows a static
   * palette.
   * @param dwell_frames Frames to hold each entry before fading; 0 chains
   * fades back to back for a continuously shifting palette.
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
    next_sequence = 0;
    dwell = dwell_frames;
    fade = fade_frames;
    easing = easing_fn;
    paused = paused_flag;
    current = 0;
    frame = 0;
    fade_active = false;

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
   * @param dwell_frames Frames to hold each palette before fading; 0 chains
   * fades back to back.
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
    next_sequence = 2;

    from_slot = allocate_palette(arena);
    to_slot = allocate_palette(arena);
    morph = allocate_palette(arena);
    provider(provider_context, 0, *from_slot);
    provider(provider_context, 1, *to_slot);
    HS_CHECK(from_slot->morph_compatible(*to_slot),
             "PaletteCycler generated palettes must be morph-compatible");
    display.bake(arena, *from_slot);
  }

  /**
   * @brief Advances the cycle by one frame.
   * @details Freezes while a wired pause flag is set. The step completing a
   * fade rebakes the display directly from the target entry, so the landing
   * is bit-exact regardless of easing endpoint behavior.
   */
  HS_COLD_MEMBER void step() {
    if ((paused != nullptr && *paused) ||
        (provider == nullptr && entry_count < 2))
      return;
    ++frame;
    if (!fade_active) {
      if (frame >= dwell) {
        if (provider == nullptr)
          begin_fade();
        fade_active = true;
        frame = 0;
      }
      return;
    }
    if (frame >= fade) {
      finish_fade();
      return;
    }
    const float progress = static_cast<float>(frame) / static_cast<float>(fade);
    const float w = easing != nullptr ? easing(progress) : progress;
    if (provider != nullptr) {
      morph->lerp(*from_slot, *to_slot, w);
      display.rebake(*morph);
    } else if ((key_morph_mask & (1u << current)) != 0) {
      morph->lerp(*entries[current].generative,
                  *entries[next_of(current)].generative, w);
      display.rebake(*morph);
    } else {
      display.rebake_crossfade(fade_from, fade_to, w);
    }
  }

  /** @brief The display LUT effects shade from. */
  const BakedPalette &palette() const { return display; }

  /** @brief Index of the entry currently dwelt on or faded away from. */
  int current_index() const { return current; }

  /** @brief True while a fade toward the next entry is in flight. */
  bool fading() const { return fade_active; }

private:
  int next_of(int index) const {
    return index + 1 == entry_count ? 0 : index + 1;
  }

  HS_COLD_MEMBER static GenerativePalette *allocate_palette(Arena &arena) {
    return new (arena.allocate(sizeof(GenerativePalette),
                               alignof(GenerativePalette))) GenerativePalette();
  }

  HS_COLD_MEMBER void finish_fade() {
    if (provider != nullptr) {
      display.rebake(*to_slot);
      GenerativePalette *retired = from_slot;
      from_slot = to_slot;
      to_slot = retired;
      provider(provider_context, next_sequence++, *to_slot);
      HS_CHECK(from_slot->morph_compatible(*to_slot),
               "PaletteCycler generated palette breaks morph compatibility");
    } else {
      bake_entry_into_display(entries[next_of(current)]);
      current = next_of(current);
    }
    fade_active = false;
    frame = 0;
  }

  HS_COLD_MEMBER static void bake_entry(BakedPalette &lut, Arena &arena,
                                        const Entry &entry) {
    // clone_from copies the prebaked table verbatim; sampling it through
    // get() would round-trip every entry through the float interpolator.
    if (entry.baked != nullptr)
      lut.clone_from(*entry.baked, arena);
    else
      lut.bake(arena, *entry.source);
  }

  HS_COLD_MEMBER void bake_entry_into_display(const Entry &entry) {
    if (entry.baked != nullptr)
      display.rebake_copy(*entry.baked);
    else
      display.rebake(*entry.source);
  }

  HS_COLD_MEMBER void begin_fade() {
    if ((key_morph_mask & (1u << current)) != 0)
      return;
    if (entries[current].baked != nullptr)
      fade_from.rebake_copy(*entries[current].baked);
    else
      fade_from.rebake(*entries[current].source);
    const Entry &to = entries[next_of(current)];
    if (to.baked != nullptr)
      fade_to.rebake_copy(*to.baked);
    else
      fade_to.rebake(*to.source);
  }

  const Entry *entries = nullptr; /**< Caller-owned entry array. */
  BakedPalette display;   /**< The LUT handed to shaders via palette(). */
  BakedPalette fade_from; /**< Crossfade w = 0 endpoint scratch. */
  BakedPalette fade_to;   /**< Crossfade w = 1 endpoint scratch. */
  /** @brief Arena-owned key-morph scratch rebaked into the display; allocated
   *  by init() only when some pair key-morphs. */
  GenerativePalette *morph = nullptr;
  NextPaletteFn provider = nullptr; /**< Generated-cycle palette source. */
  void *provider_context = nullptr; /**< Caller state handed to the provider. */
  GenerativePalette *from_slot = nullptr; /**< Generated fade w = 0 endpoint. */
  GenerativePalette *to_slot = nullptr;   /**< Generated fade w = 1 endpoint. */
  uint32_t next_sequence = 0; /**< Sequence number of the next provider call. */
  float (*easing)(float) = nullptr; /**< Fade easing; null = linear. */
  const bool *paused = nullptr; /**< Optional pause gate; null = always runs. */
  int entry_count = 0;
  int dwell = 0;              /**< Frames held on each entry between fades. */
  int fade = 0;               /**< Frames each fade spans. */
  int frame = 0;              /**< Frame counter within the current phase. */
  int current = 0;            /**< Entry dwelt on or faded away from. */
  uint8_t key_morph_mask = 0; /**< Bit i: entry i fades to its successor by
                                 key morph rather than LUT crossfade. */
  bool fade_active = false;
};
