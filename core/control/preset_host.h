/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file preset_host.h
 * @brief PresetHost: the shared preset controller — index bookkeeping, the
 *        vetoable change hook, and the manual navigation surface.
 */

#include "control/param_host.h"
#include "platform/platform.h"
#include <cstddef>
#include <cstdint>

/**
 * @brief Preset index controller for effects that expose selectable presets.
 * @details configure_presets() enables it; apply_preset() vetoes or accepts a
 * candidate before its index is committed. Manual navigation pauses the
 * parameter animations, choreography (advancePreset) leaves them running.
 */
class PresetHost : public ParamHost {
public:
  /** @brief Number of presets exposed for manual navigation. */
  size_t getPresetCount() const { return preset_count; }
  /** @brief Index of the currently selected preset. */
  size_t getPresetIndex() const { return preset_index; }
  /** @brief Selects and pauses one preset, or reports that it is unavailable. */
  bool selectPreset(size_t index) {
    if (!change_preset(index, PresetChangeOrigin::MANUAL))
      return false;
    setAnimationsPaused(true);
    return true;
  }
  /** @brief Selects one preset without changing the animation pause state. */
  bool synchronizePreset(size_t index) {
    return preset_count > 0 &&
           (index == preset_index ||
            change_preset(index, PresetChangeOrigin::SYNCHRONIZED));
  }
  /** @brief Selects and pauses the next preset. */
  bool nextPreset() {
    return preset_count > 0 && selectPreset((preset_index + 1) % preset_count);
  }
  /** @brief Selects and pauses the previous preset. */
  bool previousPreset() {
    return preset_count > 0 &&
           selectPreset((preset_index + preset_count - 1) % preset_count);
  }

protected:
  /** @brief Whether a preset move came from choreography or a user control. */
  enum class PresetChangeOrigin : uint8_t { AUTOMATIC, MANUAL, SYNCHRONIZED };

  /** @brief One validated preset transition handed to the effect hook. */
  struct PresetChange {
    size_t from;               /**< Previously committed preset index. */
    size_t to;                 /**< Candidate preset index. */
    PresetChangeOrigin origin; /**< Source of the preset change. */
  };

  /** @brief Enables the shared preset controller for this effect. */
  HS_FLASH_MEMBER void configure_presets(size_t count) {
    HS_CHECK(count > 0, "preset count must be positive");
    HS_CHECK(preset_count == 0, "presets already configured");
    preset_count = count;
  }

  /** @brief Advances choreography through the shared preset controller. */
  HS_FLASH_MEMBER bool advancePreset() {
    return preset_count > 0 && change_preset((preset_index + 1) % preset_count,
                                             PresetChangeOrigin::AUTOMATIC);
  }

  /**
   * @brief Applies a candidate preset before its index is committed.
   * @return True to commit the candidate index, false to veto the change.
   * @details A veto leaves the committed index untouched and is not rolled
   *          back, so reject before mutating any effect state.
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
