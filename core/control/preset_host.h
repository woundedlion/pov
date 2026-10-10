/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file preset_host.h
 * @brief PresetHost: tracks an effect's current preset and handles preset
 *        selection.
 */

#include "control/param_host.h"
#include "platform/platform.h"
#include <cstddef>
#include <cstdint>

/**
 * @brief Current preset index and preset selection for an effect.
 * @details configure_presets() sets the preset count. Every change goes
 * through apply_preset(), which may refuse it. Manual selection pauses the
 * parameter animations; advance_preset() does not.
 */
class PresetHost : public ParamHost {
public:
  /** @brief Number of presets exposed for manual navigation. */
  size_t getPresetCount() const { return preset_count; }
  /** @brief Index of the preset currently displayed. */
  size_t getPresetIndex() const { return displayed_preset_index(); }
  /** @brief Selects a preset and pauses animations; false if refused. */
  bool selectPreset(size_t index) {
    if (!change_preset(index, PresetChangeOrigin::MANUAL))
      return false;
    setAnimationsPaused(true);
    return true;
  }
  /** @brief Selects one preset without changing the animation pause state. */
  bool synchronizePreset(size_t index) {
    return preset_count > 0 &&
           ((index == preset_index && index == displayed_preset_index()) ||
            change_preset(index, PresetChangeOrigin::SYNCHRONIZED));
  }
  /** @brief Selects and pauses the next preset. */
  bool nextPreset() {
    return preset_count > 0 &&
           selectPreset((getPresetIndex() + 1) % preset_count);
  }
  /** @brief Selects and pauses the previous preset. */
  bool previousPreset() {
    return preset_count > 0 &&
           selectPreset((getPresetIndex() + preset_count - 1) % preset_count);
  }

protected:
  ~PresetHost() = default;

  /** @brief What caused a preset change; RESTORED is a cancelled fade
      returning to the preset it was leaving. */
  enum class PresetChangeOrigin : uint8_t {
    AUTOMATIC,
    MANUAL,
    SYNCHRONIZED,
    RESTORED
  };

  /** @brief A preset change passed to apply_preset() and preset_changed(). */
  struct PresetChange {
    size_t from;               /**< Previously committed preset index. */
    size_t to;                 /**< Candidate preset index. */
    PresetChangeOrigin origin; /**< Source of the preset change. */
  };

  /** @brief Preset whose parameters currently drive the display. */
  virtual size_t displayed_preset_index() const { return preset_index; }

  /** @brief Sets the preset count; call once with a positive count. */
  HS_FLASH_MEMBER void configure_presets(size_t count) {
    HS_CHECK(count > 0, "preset count must be positive: count=%lu",
             static_cast<unsigned long>(count));
    HS_CHECK(preset_count == 0,
             "presets already configured: count=%lu previous=%lu",
             static_cast<unsigned long>(count),
             static_cast<unsigned long>(preset_count));
    preset_count = count;
  }

  /** @brief Moves to the next preset, wrapping, as an AUTOMATIC change. */
  HS_FLASH_MEMBER bool advance_preset() {
    return preset_count > 0 && change_preset((preset_index + 1) % preset_count,
                                             PresetChangeOrigin::AUTOMATIC);
  }

  /**
   * @brief Applies a preset change before the index is updated.
   * @return True to accept the change, false to refuse it.
   * @details A refusal leaves the index unchanged but does not undo anything
   *          this function already did.
   */
  virtual bool apply_preset(const PresetChange &) { return false; }
  /** @brief Runs after a successful preset change has been committed. */
  virtual void preset_changed(const PresetChange &) {}

  bool change_preset(size_t index, PresetChangeOrigin origin) {
    if (index >= preset_count)
      return false;
    const PresetChange change{preset_index, index, origin};
    if (!apply_preset(change))
      return false;
    preset_index = index;
    preset_changed(change);
    return true;
  }

  size_t preset_count = 0;
  size_t preset_index = 0;
};
