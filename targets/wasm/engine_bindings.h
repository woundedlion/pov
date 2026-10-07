/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file engine_bindings.h
 * @brief JS-facing render bridge: the HolosphereEngine class and its embind
 *        registration.
 *
 * Owns the stack high-water instrumentation and the JS readback buffers.
 */
#pragma once

#include <emscripten/bind.h>
#include <emscripten/stack.h>
#include "targets/effects.h"
#include "core/control/registry.h"
#include "core/platform/platform.h"
#include "hardware/pov_segment_map.h"
#include "targets/wasm/arena_metrics.h"
#include "targets/wasm/workbench_bindings.h"
#include "targets/wasm/effect_factory.h"
#include "targets/wasm/param_marshal.h"
#include "targets/wasm/wasm_predicates.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstring>
#include <climits>
#include <memory>
#include <string>

// ---- Stack canary painting for high water mark tracking ----
inline constexpr uint8_t STACK_CANARY = 0xCD;

/**
 * @brief Lowest stack address the canary may occupy.
 * @return The stack end, raised past the runtime's stack cookie.
 * @details Clears both placements of the two -sASSERTIONS cookie words, so
 *          checkStackCookie() never reports a false overflow.
 */
static uintptr_t stack_canary_floor() {
  return emscripten_stack_get_end() + 4 * sizeof(uint32_t);
}

/**
 * @brief Paints the unused portion of the stack with a canary byte for high
 *        water mark tracking.
 * @details Writes STACK_CANARY from the current stack pointer down to the canary
 *          floor. Safe to call at any time — only touches memory below the
 *          current stack pointer.
 */
static void stack_paint_canary() {
  uintptr_t sp = emscripten_stack_get_current();
  uintptr_t floor = stack_canary_floor();
  if (sp > floor) {
    std::memset(reinterpret_cast<void *>(floor), STACK_CANARY, sp - floor);
  }
}

/**
 * @brief Computes the stack high water mark by scanning the canary region.
 * @return Number of stack bytes that have been touched, in bytes (the high
 *         water mark), found by scanning from the canary floor upward to the
 *         first overwritten canary byte.
 * @details A lower bound: a frame that wrote 0xCD, or reserved space it never
 *          stored to, reads back as still-canary.
 */
static size_t stack_high_water_mark() {
  uintptr_t base = emscripten_stack_get_base();
  const uint8_t *p = reinterpret_cast<const uint8_t *>(stack_canary_floor());
  const uint8_t *top = reinterpret_cast<const uint8_t *>(base);
  constexpr size_t WORD = sizeof(uint64_t);
  constexpr uint64_t CANARY_WORD = 0x0101010101010101ull * STACK_CANARY;
  while (p < top && (reinterpret_cast<uintptr_t>(p) % WORD) != 0 &&
         *p == STACK_CANARY)
    p++;
  while (p + WORD <= top) {
    uint64_t w;
    std::memcpy(&w, p, WORD);
    if (w != CANARY_WORD)
      break;
    p += WORD;
  }
  while (p < top && *p == STACK_CANARY)
    p++;
  return static_cast<size_t>(top - p);
}

// Running max of effect construction + init() stack depth; survives repaints.
static size_t init_stack_peak = 0;

// Bound on one effect's exposed parameters.
inline constexpr size_t MAX_PARAMS = hs_wasm::ParamStreams::CAPACITY;

#if HS_ENABLE_CHAIN_INTERPRETER
static_assert(MAX_PARAMS >= Pullback::Interp::MAX_CHAIN_PARAMS,
              "MAX_PARAMS must cover the full chain parameter schema");
#endif

static_assert(MAX_PARAMS >= Effect::ParamList::FIXED_CAPACITY,
              "MAX_PARAMS must cover ParamList's default array size");

/**
 * @brief Outcome of a HolosphereEngine::setClip() call.
 * @details NO_EFFECT is the ordinary state after a RESIZED setResolution();
 *          INVALID_BOUNDS is a caller bug. APPLIED and FULL_FRAME_KEPT are both
 *          successes, but only APPLIED means the band is in force. Exposed to JS
 *          as the Module.ClipSetResult embind enum. Compare this and every other
 *          embind result enum by value, never by truthiness: every enum value
 *          is a truthy object.
 */
enum class ClipSetResult {
  APPLIED,         /**< Band installed; rendering is narrowed to it. */
  NO_EFFECT,       /**< No effect is installed to receive the clip. */
  INVALID_BOUNDS,  /**< Bounds malformed or out of range for the resolution. */
  FULL_FRAME_KEPT, /**< Bounds accepted but ignored: the effect reports
                        needs_full_frame() or persists_pixels(), so the clip stays at the
                        full canvas. */
};

/**
 * @brief Outcome of a HolosphereEngine::setResolution() call.
 * @details Exposed to JS as the Module.ResolutionSetResult embind enum.
 */
enum class ResolutionSetResult {
  RESIZED,        /**< Resolution switched; the effect was torn down, so
                       setEffect() and any clip must be re-applied. */
  ALREADY_ACTIVE, /**< Request matched the active resolution; pure no-op —
                       the current effect and clip stay live. */
  UNSUPPORTED,    /**< Size is not an HS_RESOLUTIONS row; ignored, prior
                       valid state kept. */
};

/**
 * @brief Outcome of a HolosphereEngine::setEffect() call.
 * @details Both rejections keep the prior effect. UNSUPPORTED_RESOLUTION means
 *          no name can succeed until a supported setResolution(). Exposed to JS
 *          as the Module.EffectSetResult embind enum.
 */
enum class EffectSetResult {
  INSTALLED,              /**< Fresh effect instantiated at default parameters
                               and the full-canvas clip. */
  UNKNOWN_EFFECT,         /**< Name not registered at the active resolution;
                               prior effect kept. */
  UNSUPPORTED_RESOLUTION, /**< Active resolution has no factory; prior state
                               kept. */
};

// True while a HolosphereEngine is constructed-but-not-deleted. The engine is a
// singleton: its Effect aliases shared static buffers and its arenas are
// module-global, so a second instance would corrupt the first's frames.
static bool engine_alive = false;

/**
 * @brief JS-facing render engine driving one resolution/effect at a time.
 * @details Owns the current effect and the stable readback buffers.
 *          At most one instance may be live: delete() the current engine before
 *          constructing another, or the constructor traps. Test isLive() before construction.
 */
class HolosphereEngine {
public:
  /**
   * @brief Constructs the engine with a valid default resolution and effect.
   * @details Pre-sizes the JS-facing readback buffers to their maximum extent so
   *          their backing storage never moves (the WASM memory-view contract).
   */
  HolosphereEngine() {
    HS_CHECK(!engine_alive,
             "HolosphereEngine is a singleton: delete() the live instance "
             "before constructing another (its Effect and arenas are shared "
             "module-global storage)");
    // Module global: reset so a fresh engine never inherits a predecessor's.
    Render::pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
    apply_display_geometry(0.0f, math::PI_F);
    stack_paint_canary();

    // Pre-size once: getPixels/getParamValues views alias this storage, so it
    // must never reallocate.
    pixel_buffer.assign(MAX_W * MAX_H * CHANNELS, 0);

    // Bootstrap on the first row of each roster.
    const bool bootstrap_resolution_set =
        setResolution(hs_wasm::WASM_RESOLUTIONS[0].w,
                      hs_wasm::WASM_RESOLUTIONS[0].h) ==
        ResolutionSetResult::RESIZED;
    HS_CHECK(bootstrap_resolution_set,
             "the first HS_RESOLUTIONS row must be dispatchable here");
    const bool bootstrap_effect_set =
        setEffect(hs_wasm::EFFECT_REGISTRATIONS[0].name.data()) ==
        EffectSetResult::INSTALLED;
    HS_CHECK(bootstrap_effect_set,
             "the first HS_EFFECT_LIST entry must be registered and buildable "
             "at the first HS_RESOLUTIONS row");

    engine_alive = true;
  }

  /**
   * @brief Destroys the engine and admits the next construction.
   * @details Reached from JS via delete(). Re-partitions the module-global
   *          engine arenas, so their metrics report the released state.
   */
  ~HolosphereEngine() {
#if HS_ENABLE_CHAIN_INTERPRETER
    HS_CHECK(!snapshot_decode_active,
             "delete() from a caller accessor during engine decode");
#endif
    teardown_effect();
    binding_state->alive = false;
    engine_alive = false;
  }

  /**
   * @brief Whether an engine instance is currently constructed.
   * @return True from the end of a successful construction until that
   *         instance's delete().
   * @details Exposed to JS as the static Module.HolosphereEngine.isLive().
   *          Construction over a live instance traps.
   */
  static bool isLive() { return engine_alive; }

  /**
   * @brief Sets the missing arc at each pole as a percentage in [0, 25].
   * @return False for non-finite or out-of-range inputs; otherwise true.
   * @details Changes rebuild geometry caches, retaining parameters, preset,
   *          pause and clip; ShaderChain restores its full executable snapshot,
   *          while other effects restart animation and trail history.
   *          Identical geometry leaves the effect and parameter generation intact.
   */
  bool setDisplayCaps(double top_percent, double bottom_percent) {
    if (!std::isfinite(top_percent) || !std::isfinite(bottom_percent) ||
        top_percent < 0.0 || top_percent > 25.0 || bottom_percent < 0.0 ||
        bottom_percent > 25.0)
      return false;
    const float NORTH = static_cast<float>(top_percent / 100.0) * math::PI_F;
    const float SOUTH =
        (1.0f - static_cast<float>(bottom_percent / 100.0)) * math::PI_F;
    if (NORTH == getDisplayNorthPhi() && SOUTH == getDisplaySouthPhi())
      return true;
    if (!current_effect) {
      apply_display_geometry(NORTH, SOUTH);
      return true;
    }
    hs_wasm::dispatch_resolution(
        pixel_width, pixel_height,
        [&]<int W, int H>() { rebuild_display_geometry<W, H>(NORTH, SOUTH); });
    return true;
  }

  /** @brief Polar angle in radians of the first displayed LED row. */
  float getDisplayNorthPhi() const { return math::DISPLAY_NORTH_PHI; }

  /** @brief Polar angle in radians of the last displayed LED row. */
  float getDisplaySouthPhi() const { return math::DISPLAY_SOUTH_PHI; }

  /**
   * @brief Switches the active canvas resolution.
   * @param w Requested canvas width in pixels.
   * @param h Requested canvas height in pixels.
   * @return RESIZED if the resolution switched, ALREADY_ACTIVE if the request
   *         matched the active resolution (a pure no-op — nothing is torn
   *         down), or UNSUPPORTED if the request was rejected and the previous
   *         valid state was kept.
   * @details RESIZED tears down the current effect: call setEffect() before the
   *          next drawFrame() (else it renders blank) and re-apply any setClip(),
   *          whose bounds are not rescaled. The teardown re-partitions the
   *          engine arenas. An outstanding getPixels() view keeps aliasing live
   *          memory at the previous length; re-fetch it.
   */
  ResolutionSetResult setResolution(double w, double h) {
    if (!hs_wasm::wasm_resolution_supported(w, h)) {
      hs::log("WASM: Unsupported resolution %gx%g — ignored", w, h);
      return ResolutionSetResult::UNSUPPORTED;
    }

    if (w == pixel_width && h == pixel_height)
      return ResolutionSetResult::ALREADY_ACTIVE;

    pixel_width = static_cast<int>(w);
    pixel_height = static_cast<int>(h);
    const int count = pixel_width * pixel_height * CHANNELS;
    std::fill_n(pixel_buffer.data(), count, uint16_t{0});

    if (current_effect)
      teardown_effect();
    return ResolutionSetResult::RESIZED;
  }

  /**
   * @brief Tears down the current effect and instantiates the named one at the
   *        active resolution.
   * @param name Effect class name or stable EFFECT_ID to instantiate.
   * @return INSTALLED iff an effect was actually instantiated, else the
   *         rejection reason — UNKNOWN_EFFECT for an unknown/stale effect name
   *         or UNSUPPORTED_RESOLUTION.
   * @details A rejected name keeps the prior effect. INSTALLED resets the clip
   *          to the full canvas and parameters to defaults; re-apply setClip()
   *          and setParameter(). The animation pause state is retained.
   */
  EffectSetResult setEffect(const std::string &name) {
    const FactoryEntry *entry = nullptr;
    const bool dispatched = hs_wasm::dispatch_resolution(
        pixel_width, pixel_height, [&]<int W, int H>() {
          entry = hs_wasm::find_factory_entry<W, H>(name);
        });
    if (!dispatched) {
      hs::log("WASM: setEffect at unsupported resolution %dx%d — keeping "
              "current effect",
              pixel_width, pixel_height);
      return EffectSetResult::UNSUPPORTED_RESOLUTION;
    }

    if (!entry) {
      hs::log("WASM: setEffect unknown effect '%s' — keeping current effect",
              name.c_str());
      return EffectSetResult::UNKNOWN_EFFECT;
    }

    teardown_effect();

    hs_wasm::dispatch_resolution(pixel_width, pixel_height, []<int W, int H>() {
      math::init_geometry_luts<W,
                               H>(); // eager-fill LUTs before the first frame
    });
    // Per-load RNG stream keyed by the effect's stable id; seeded after
    // teardown so a teardown draw cannot advance it.
    hs::random().seed(hs::stable_effect_seed(entry->stable_id));
    current_effect = entry->creator();
    current_factory_entry = entry;
    binding_state->effect = current_effect.get();
    binding_state->entry = entry;
    binding_state->width = pixel_width;
    binding_state->height = pixel_height;
    current_effect->setAnimationsPaused(binding_state->paused);
    current_effect->init();
    hs_wasm::check_param_capacity(*current_effect);
    param_generation.replace(current_effect->getParameterSchemaGeneration());
    const size_t init_hwm = stack_high_water_mark();
    if (init_hwm > init_stack_peak)
      init_stack_peak = init_hwm;
    stack_paint_canary();
    return EffectSetResult::INSTALLED;
  }

  /**
   * @brief Restricts rendering to a clip band for the current effect.
   * @param x0 Inclusive left column of the clip band in [0, pixel_width].
   * @param x1 Exclusive right column, with x0 <= x1 <= pixel_width.
   * @param y0 Inclusive top row of the clip band in [0, pixel_height].
   * @param y1 Exclusive bottom row of the clip band, with y0 <= y1 <= pixel_height.
   *           All four must be integral; a fractional or NaN number is
   *           INVALID_BOUNDS.
   * @return APPLIED if the band was installed, FULL_FRAME_KEPT if the bounds
   *         were accepted but the effect keeps the full-canvas clip, otherwise
   *         NO_EFFECT or INVALID_BOUNDS.
   * @details Malformed input is rejected without trapping.
   *          See docs/specs/segmented_stateful_effects_spec.md.
   */
  ClipSetResult setClip(double x0, double x1, double y0, double y1) {
    if (!current_effect)
      return ClipSetResult::NO_EFFECT;
    if (!hs_wasm::clip_bounds_valid(x0, x1, y0, y1, pixel_width,
                                    pixel_height)) {
      hs::log("WASM: setClip bounds out of range (x0=%g,x1=%g,y0=%g,y1=%g) — "
              "ignored",
              x0, x1, y0, y1);
      return ClipSetResult::INVALID_BOUNDS;
    }
    // A band-clipped stateful effect has stale cv.prev outside its band.
    if (!pov::segment_clip_applies(current_effect->needs_full_frame(),
                                   current_effect->persists_pixels()))
      return ClipSetResult::FULL_FRAME_KEPT;
    current_effect->set_clip(static_cast<int>(y0), static_cast<int>(y1),
                             static_cast<int>(x0), static_cast<int>(x1));
    return ClipSetResult::APPLIED;
  }

  /**
   * @brief Renders one frame of the current effect into the JS-facing buffer.
   * @details Copies the effect's canvas into pixel_buffer as 16-bit linear RGB
   *          triples; clears the active readback if no effect is set. Only the
   *          display clip is copied; pixels outside it keep their previous
   *          readback values.
   */
  void drawFrame() {
    if (!current_effect) {
      // No active effect: clear the active prefix so getPixels() hands JS a
      // blank frame at the current resolution, not stale content.
      const int count = pixel_width * pixel_height * CHANNELS;
      std::fill_n(pixel_buffer.data(), count, uint16_t{0});
      return;
    }

    current_effect->draw_frame();
    current_effect->advance_display();

    static_assert(static_cast<long long>(MAX_W) * MAX_H * CHANNELS <= INT_MAX,
                  "drawFrame pixel-index accumulators are int");
    // Display bounds, not render bounds: the margin expansion feeds stateful
    // filters and is never shown.
    const ClipRegion &band = current_effect->clip();
    if (!current_effect->overrides_get_pixel() &&
        current_effect->output_envelope_u16() == 65535u) {
      // Fast path: display_buffer()[i] == get_pixel(x, y), so copy directly.
      const Pixel *buf = current_effect->display_buffer();
      static_assert(sizeof(Pixel) == 3 * sizeof(uint16_t),
                    "fast-path memcpy assumes packed RGB16 Pixel layout");
      if (band.is_full()) {
        const int count = pixel_width * pixel_height;
        std::memcpy(pixel_buffer.data(), buf,
                    static_cast<size_t>(count) * sizeof(Pixel));
      } else {
        const size_t row_bytes =
            static_cast<size_t>(band.x_end - band.x_start) * sizeof(Pixel);
        for (int y = band.y_start; y < band.y_end; y++) {
          const int first = y * pixel_width + band.x_start;
          std::memcpy(pixel_buffer.data() +
                          static_cast<size_t>(first) * CHANNELS,
                      buf + first, row_bytes);
        }
      }
    } else {
      for (int y = band.y_start; y < band.y_end; y++) {
        int idx = (y * pixel_width + band.x_start) * CHANNELS;
        for (int x = band.x_start; x < band.x_end; x++) {
          const Pixel source = current_effect->get_pixel(x, y);
          const Pixel p = current_effect->apply_output_envelope(source);
          pixel_buffer[idx++] = p.r;
          pixel_buffer[idx++] = p.g;
          pixel_buffer[idx++] = p.b;
        }
      }
    }
  }

  /**
   * @brief Reports the active effect's POV column-strobe mode for the simulator.
   * @return true if the effect strobes each column to black after it is shown
   *         (discrete columns with dark gaps); false if columns persist and
   *         smear horizontally into the next (a continuous, gap-free band).
   *         false when no effect is set.
   * @details See Effect::strobe_columns (core/render/canvas.h).
   */
  bool strobeColumns() const {
    return current_effect ? current_effect->strobe_columns() : false;
  }

  /**
   * @brief Exposes the raw pixel buffer to JS as a zero-copy Uint16Array view.
   * @return Typed memory view over the active resolution's R,G,B pixels within
   *         the stable MAX_W*MAX_H*3 backing buffer.
   * @details The view aliases WASM memory; it is not a copy. Heap growth
   *          detaches it (buffer.byteLength === 0); a successful setResolution()
   *          leaves it attached at the wrong length. A caller caching the view
   *          across frames must test both buffer.byteLength !== 0 and
   *          length === getBufferLength().
   */
  emscripten::val getPixels() {
    return emscripten::val(emscripten::typed_memory_view(
        pixel_width * pixel_height * CHANNELS, pixel_buffer.data()));
  }

  /**
   * @brief Returns the length of the active pixel buffer view.
   * @return Number of uint16 elements in the active view (pixel_width *
   *         pixel_height * 3, three channels per pixel).
   * @details A cached getPixels() view whose length differs is stale.
   */
  int getBufferLength() const { return pixel_width * pixel_height * CHANNELS; }

  /**
   * @brief Updates one named effect parameter.
   * @param name Parameter name to update.
   * @param value New parameter value, in the parameter's native units.
   * @return APPLIED if the write was accepted, otherwise the rejection reason:
   *         NO_EFFECT (no effect is set), UNKNOWN_PARAM (the name is unknown to
   *         the effect), READONLY (engine-written telemetry the GUI must not
   *         poke), NON_FINITE (rejected before it can poison render math), or
   *         INADMISSIBLE (an unlisted option ID or a cross-parameter constraint
   *         refuses the value). Exposed to JS as the Module.ParamSetResult
   *         embind enum. An APPLIED float is silently clamped to the
   *         parameter's registered [min,max]; read the effective value back via
   *         getParamValues().
   * @details An APPLIED write to an *animated* param engages the animation
   *          pause, as setAnimationsPaused(true) would; read it back through
   *          getAnimationsPaused().
   */
  ParamSetResult setParameter(const std::string &name, float value) {
    if (!current_effect)
      return ParamSetResult::NO_EFFECT;
    const ParamSetResult result =
        current_effect->updateParameter(name.c_str(), value);
    binding_state->paused = current_effect->animations_paused();
    return result;
  }

  /**
   * @brief Pauses or resumes the current effect's parameter animations.
   * @param paused true to pause animations, false to resume.
   * @details Retained across effect and resolution changes and applied to the
   *          next effect when no effect is currently loaded.
   */
  void setAnimationsPaused(bool paused) {
    binding_state->paused = paused;
    if (current_effect)
      current_effect->setAnimationsPaused(paused);
  }

  /**
   * @brief Reports whether the engine's parameter animations are paused.
   * @return true while animations are frozen.
   * @details Also engaged by an APPLIED write to an animated param.
   */
  bool getAnimationsPaused() const { return binding_state->paused; }

  /**
   * @brief Number of presets the current effect exposes for manual navigation.
   * @return The preset count, 0 when no effect is set or the effect authored
   *         none.
   */
  uint32_t getPresetCount() const {
    return current_effect
               ? static_cast<uint32_t>(current_effect->getPresetCount())
               : 0;
  }

  /**
   * @brief Index of the effect's currently selected preset.
   * @return The index, 0 when no effect is set — indistinguishable from a real
   *         selection of preset 0, so a caller tells the two apart by
   *         getPresetCount() != 0.
   * @details Also moves with engine-driven preset advancement, without any JS
   *          call.
   */
  uint32_t getPresetIndex() const {
    return current_effect
               ? static_cast<uint32_t>(current_effect->getPresetIndex())
               : 0;
  }

  /**
   * @brief Stable preset IDs in numeric navigation order.
   * @return The preset ID array, empty when no effect or ID table is set.
   */
  emscripten::val getPresetIds() const {
    emscripten::val ids = emscripten::val::array();
    if (!current_effect || !current_factory_entry ||
        !current_factory_entry->preset_id)
      return ids;
    for (size_t index = 0; index < current_effect->getPresetCount(); ++index)
      ids.set(index, std::string(current_factory_entry->preset_id(index)));
    return ids;
  }

  /**
   * @brief Selects a preset through its persisted identity.
   * @param preset_id Persisted preset identity.
   * @return true when the preset was found and applied; false when no effect
   * or ID table is set, the ID is empty or unknown, or selectPreset refuses it.
   * @details Engages the animation pause as selectPreset does. A rejected call
   * leaves the preset and pause untouched.
   */
  bool selectPresetById(const std::string &preset_id) {
    if (!current_effect || !current_factory_entry ||
        !current_factory_entry->preset_id || preset_id.empty())
      return false;
    for (size_t index = 0; index < current_effect->getPresetCount(); ++index)
      if (current_factory_entry->preset_id(index) == preset_id)
        return selectPreset(static_cast<double>(index));
    return false;
  }

  /**
   * @brief Selects one preset and freezes the effect's parameter animations.
   * @param index Preset to select, in [0, getPresetCount()). Must be integral;
   *        a fractional, NaN or out-of-range number is rejected.
   * @return true when the preset was applied; false when no effect is set, the
   *         index is malformed or out of range, or the effect refused the
   *         preset.
   * @details Engages the animation pause as setAnimationsPaused(true) would.
   *          Parameter values move with the preset; re-read them via
   *          getParamValues(). The index is a double so wrapped and NaN inputs
   *          are rejected rather than coerced.
   */
  bool selectPreset(double index) {
    if (!current_effect || !preset_index_accepted(index, "selectPreset") ||
        !current_effect->selectPreset(static_cast<size_t>(index)))
      return false;
    binding_state->paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Selects one preset without touching the animation pause state.
   * @param index Preset to select, in [0, getPresetCount()). Must be integral;
   *        a fractional, NaN or out-of-range number is rejected.
   * @return true when the preset is the active one already or was applied;
   *         false when no effect is set, the index is malformed or out of
   *         range, or the effect refused the preset.
   * @details Leaves choreography running. A request for the active index is a
   *          success no-op.
   */
  bool synchronizePreset(double index) {
    if (!current_effect || !preset_index_accepted(index, "synchronizePreset") ||
        !current_effect->synchronizePreset(static_cast<size_t>(index)))
      return false;
    binding_state->paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Selects the next preset, wrapping past the last, and freezes
   *        animations.
   * @return true when the preset was applied; false when no effect is set, the
   *         effect has no presets, or it refused the preset.
   * @details Same animation pause as selectPreset().
   */
  bool nextPreset() {
    if (!current_effect || !current_effect->nextPreset())
      return false;
    binding_state->paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Selects the previous preset, wrapping past the first, and freezes
   *        animations.
   * @return true when the preset was applied; false when no effect is set, the
   *         effect has no presets, or it refused the preset.
   * @details Same animation pause as selectPreset().
   */
  bool previousPreset() {
    if (!current_effect || !current_effect->previousPreset())
      return false;
    binding_state->paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Sets near-pole azimuthal shading decimation.
   * @param aggressiveness Columns per shade are this over sin(colatitude);
   *        0 disables. NaN and negative inputs clamp to 0; positive infinity to 8.
   * @details Engine-scoped: the constructor restores HS_POLE_LOD_DEFAULT, and
   *          each WASM instance carries its own.
   */
  void setPoleLod(float aggressiveness) {
    Render::pole_lod_aggressiveness =
        hs_wasm::clamp_pole_lod_aggressiveness(aggressiveness);
  }

  /**
   * @brief Current near-pole azimuthal decimation aggressiveness.
   * @return The clamped value of the last setPoleLod() on this engine, else
   *         HS_POLE_LOD_DEFAULT.
   */
  float getPoleLod() const { return Render::pole_lod_aggressiveness; }

  /**
   * @brief Builds the GUI's parameter descriptor list.
   * @return JS array with one {name, value, requestedValue, acceptedValue,
   *         animated, readonly, preset} object per param in declaration order,
   *         with {warning} when present, {min, max} on every non-boolean param,
   *         {step} on every whole-number param, and {options} (plus
   *         {exportOptions} for C++ enum literals) on every enum param; empty
   *         when no effect is set.
   * @details `value` is the rendered state for display; `requestedValue` is the
   *          writable target used to seed another renderer; `acceptedValue` is
   *          the accepted target and `warning` describes any adjustment or
   *          rejection. Boolean values are JS booleans; every other value is a
   *          number. An enum's optionValues maps its labels to numeric IDs;
   *          when absent, values index options directly. preset marks the
   *          params a preset export carries. The order matches
   *          getParamValues(); pin getParamGeneration() beside a snapshot to
   *          detect a rebind.
   */
  emscripten::val getParameterDefinitions() {
    if (!current_effect)
      return emscripten::val::array();

    current_effect->refresh_parameter_display();
    emscripten::val result = emscripten::val::array();
    // Walks the registered ParamList in declaration order.
    hs_wasm::collect_param_views(*current_effect, param_streams.views);

    int i = 0;
    for (const auto &v : param_streams.views) {
      emscripten::val entry = emscripten::val::object();
      entry.set("name", emscripten::val(v.name));

      if (v.is_bool) {
        // Booleans carry no range.
        entry.set("value", emscripten::val(v.value > 0.5f));
        entry.set("requestedValue", emscripten::val(v.requested_value > 0.5f));
        entry.set("acceptedValue", emscripten::val(v.accepted_value > 0.5f));
      } else {
        entry.set("value", v.value);
        entry.set("requestedValue", v.requested_value);
        entry.set("acceptedValue", v.accepted_value);
        entry.set("min", v.min);
        entry.set("max", v.max);
        if (v.is_integer)
          entry.set("step", 1);
        if (v.option_count > 0) {
          emscripten::val opts = emscripten::val::array();
          for (int k = 0; k < v.option_count; ++k)
            opts.set(k, emscripten::val(v.options[k]));
          entry.set("options", opts);
          if (v.option_values != nullptr) {
            emscripten::val ids = emscripten::val::array();
            for (int k = 0; k < v.option_count; ++k)
              ids.set(k, static_cast<double>(v.option_values[k]));
            entry.set("optionValues", ids);
          }
          if (v.export_options != nullptr) {
            emscripten::val export_opts = emscripten::val::array();
            for (int k = 0; k < v.option_count; ++k)
              export_opts.set(k, emscripten::val(v.export_options[k]));
            entry.set("exportOptions", export_opts);
          }
        }
      }
      entry.set("animated", emscripten::val(v.animated));
      entry.set("readonly", emscripten::val(v.readonly));
      entry.set("preset", emscripten::val(v.preset));
      if (const char *warning = current_effect->parameter_warning(v.name))
        entry.set("warning", emscripten::val(warning));
      result.set(i++, entry);
    }
    return result;
  }

  /**
   * @brief Streams the current param values to the GUI per frame.
   * @return Zero-copy Float32Array view over the current param values, in the
   *         same order as getParameterDefinitions(); empty array if no effect is
   *         set.
   * @details Same memory-view contract as getPixels(): the view aliases WASM
   *          memory and must be consumed before the next allocation. Never
   *          reallocates (size <= MAX_PARAMS), so it detaches no other view.
   */
  emscripten::val getParamValues() {
    if (!current_effect) {
      // clear() keeps the reserved capacity, so the zero-length view's backing
      // pointer stays valid.
      param_streams.values.clear();
      return emscripten::val(emscripten::typed_memory_view(
          param_streams.values.size(), param_streams.values.data()));
    }

    current_effect->refresh_parameter_display();

    // Same order as getParameterDefinitions().
    hs_wasm::fill_param_values(*current_effect, param_streams.values);
    return emscripten::val(emscripten::typed_memory_view(
        param_streams.values.size(), param_streams.values.data()));
  }

  /**
   * @brief Identity token joining a getParameterDefinitions() snapshot to a
   *        later getParamValues() read.
   * @return A fresh value after every effect replacement or descriptor rebind.
   * @details Parameter counts can repeat even when names and order change. Pin
   *          this at definition-snapshot time and compare it beside each value
   *          read; a change means the definitions must be rebuilt.
   */
  uint32_t getParamGeneration() {
    if (current_effect)
      param_generation.observe(current_effect->getParameterSchemaGeneration());
    return param_generation.generation();
  }

  /**
   * @brief Reports engine arena and stack metrics for the JS memory HUD.
   * @return JS object of the engine arenas' metrics ({usage,
   *         high_water_mark, lifetime_high_water_mark, capacity}) plus a
   *         "stack" entry ({high_water_mark, init_high_water_mark, capacity}),
   *         all in bytes.
   * @details Excludes the tooling arenas. An arena's `high_water_mark` covers
   *          only the window since its last peak reset or rebind, and an
   *          overrun check compares it with `capacity`.
   *          `lifetime_high_water_mark` is the sizing peak and may exceed the
   *          current `capacity` after configure_arenas() re-splits the arenas.
   *          On the stack entry, `high_water_mark` is the canary's live
   *          reading, which a repaint resets; `init_high_water_mark` is the
   *          latched deepest effect construction + init().
   */
  emscripten::val getArenaMetrics() const {
    emscripten::val metrics = collect_engine_arena_metrics();

    // Stack region. No running usage: the live depth at this call is outside any
    // render, so only the high-water marks are meaningful.
    {
      uintptr_t base = emscripten_stack_get_base();
      uintptr_t end = emscripten_stack_get_end();
      emscripten::val m = emscripten::val::object();
      m.set("high_water_mark", static_cast<size_t>(stack_high_water_mark()));
      m.set("init_high_water_mark", init_stack_peak);
      m.set("capacity", static_cast<size_t>(base >= end ? base - end : 0));
      metrics.set("stack", m);
    }

    return metrics;
  }

  /**
   * @brief Returns the effect-size map for the active resolution.
   * @return JS object mapping each effect name to its sizeof, in bytes, at the current
   *         resolution; empty map if unsupported/uninitialized.
   */
  emscripten::val getEffectSizes() const {
    emscripten::val sizes = emscripten::val::object();
    hs_wasm::dispatch_resolution(
        pixel_width, pixel_height,
        [&]<int W, int H>() { sizes = get_effect_sizes_helper<W, H>(); });
    return sizes;
  }

  /**
   * @brief Returns the authored preset counts for the active resolution.
   * @return JS object mapping every effect name to its preset count; empty map
   *         if unsupported or uninitialized.
   */
  emscripten::val getEffectPresetCounts() const {
    emscripten::val counts = emscripten::val::object();
    hs_wasm::dispatch_resolution(
        pixel_width, pixel_height, [&]<int W, int H>() {
          for (const auto &entry : hs_wasm::get_factory<W, H>())
            counts.set(std::string(entry.name), entry.preset_count);
        });
    return counts;
  }

#if HS_ENABLE_CHAIN_INTERPRETER
  /** @brief Acquires a chain authoring handle, or null for other effects. */
  std::shared_ptr<ShaderChainBindings> getShaderChainBindings() const {
    return acquire_shader_chain_bindings(binding_state);
  }
#endif

  /**
   * @brief Enumerates the resolutions the factory can build.
   * @return JS array of [W, H] pairs from WASM_RESOLUTIONS (generated from
   *         HS_RESOLUTIONS).
   */
  static emscripten::val getSupportedResolutions() {
    emscripten::val out = emscripten::val::array();
    int i = 0;
    for (const hs_wasm::WasmResolution &row : hs_wasm::WASM_RESOLUTIONS) {
      emscripten::val pair = emscripten::val::array();
      pair.set(0, row.w);
      pair.set(1, row.h);
      out.set(i++, pair);
    }
    return out;
  }

private:
  void teardown_effect() {
    ++binding_state->generation;
    binding_state->effect = nullptr;
    binding_state->entry = nullptr;
    param_streams.views.clear();
    current_effect.reset();
    current_factory_entry = nullptr;
    configure_arenas_default();
    param_generation.replace(0);
    stack_paint_canary();
  }
  static void apply_display_geometry(float north, float south) {
    HS_CHECK(math::set_display_geometry(north, south),
             "Validated display geometry must be accepted");
#define HS_REFRESH_DISPLAY_GEOMETRY(W, H) math::init_geometry_luts<W, H>();
    HS_RESOLUTIONS(HS_REFRESH_DISPLAY_GEOMETRY)
#undef HS_REFRESH_DISPLAY_GEOMETRY
  }

  template <int W, int H>
  void rebuild_display_geometry(float north, float south) {
    const std::string NAME(current_factory_entry->name);
    const size_t PRESET = current_effect->getPresetIndex();
    const bool HAS_PRESET = current_effect->getPresetCount() != 0;
    const bool PAUSED = binding_state->paused;
    const ClipRegion CLIP = current_effect->clip();
    std::vector<std::pair<std::string, float>> parameters;
    current_effect->refresh_parameter_display();
    for (const auto &parameter : current_effect->getParameters())
      if (!parameter.readonly)
        parameters.emplace_back(parameter.name, parameter.get_requested());

    WorkbenchBindings::RebuildRestore restore_workbench;
#if HS_ENABLE_CHAIN_INTERPRETER
    if (auto bindings = getShaderChainBindings())
      restore_workbench = bindings->capture_rebuild_state();
#endif

    teardown_effect();
    apply_display_geometry(north, south);
    HS_CHECK(setEffect(NAME) == EffectSetResult::INSTALLED,
             "Geometry rebuild must reinstall the current effect");
    if (HAS_PRESET)
      HS_CHECK(current_effect->synchronizePreset(PRESET),
               "Geometry rebuild must restore the selected preset");

    const bool parameters_restored =
        restore_workbench && restore_workbench(binding_state);
    HS_CHECK(!restore_workbench || parameters_restored,
             "Geometry rebuild must restore the workbench configuration");
    if (!parameters_restored)
      current_effect->replay_parameter_writes(parameters);
    setAnimationsPaused(PAUSED);
    setClip(CLIP.x_start, CLIP.x_end, CLIP.y_start, CLIP.y_end);
    param_generation.observe(current_effect->getParameterSchemaGeneration());
  }

  /**
   * @brief Builds the {effect name -> sizeof, in bytes} map for the (W,H) factory.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @return JS object mapping each effect name to its sizeof, in bytes, for the GUI.
   */
  template <int W, int H> static emscripten::val get_effect_sizes_helper() {
    emscripten::val s = emscripten::val::object();
    const auto &factory = hs_wasm::get_factory<W, H>();
    for (const auto &entry : factory)
      s.set(std::string(entry.name), static_cast<int>(entry.size));
    return s;
  }

  /**
   * @brief Range-checks a preset index arriving from the untyped JS boundary.
   * @param index Requested index, as the double the binding takes.
   * @param call Entry-point name for the rejection log.
   * @return true iff the index is integral and inside the effect's roster.
   * @details Logs, never traps: a trap aborts the whole WASM module. Requires a
   *          live effect.
   */
  bool preset_index_accepted(double index, const char *call) const {
    if (hs_wasm::preset_index_valid(index, current_effect->getPresetCount()))
      return true;
    hs::log("WASM: %s index out of range (%g) — ignored", call, index);
    return false;
  }

  /** Channels per pixel in the readback buffer (linear RGB triples). */
  static constexpr int CHANNELS = 3;

  std::unique_ptr<Effect>
      current_effect; /**< Currently active effect, or null. */
  std::shared_ptr<WorkbenchBindingState> binding_state =
      std::make_shared<WorkbenchBindingState>();
  const FactoryEntry *current_factory_entry =
      nullptr; /**< Registry metadata for current_effect; WASM lifetime-static. */
  std::vector<uint16_t> pixel_buffer; /**< 16-bit linear RGB readback buffer. */
  hs_wasm::ParamStreams param_streams;
  int pixel_width = 0;  /**< Active canvas width in pixels. */
  int pixel_height = 0; /**< Active canvas height in pixels. */
  hs_wasm::ParamGenerationTracker
      param_generation; /**< Effect and descriptor-schema identity token. */
};

/** @brief Registers the render bridge's enums and class with Embind. */
static void bind_engine() {
  bind_workbench_adapters();
  // Boot-time defaults; the geometry getters report the active display caps.
  emscripten::constant("DISPLAY_PROFILE", HS_DISPLAY_PROFILE);
  emscripten::constant("DISPLAY_NORTH_PHI",
                       math::DisplayGeometry<144>::NORTH_PHI);
  emscripten::constant("DISPLAY_SOUTH_PHI",
                       math::DisplayGeometry<144>::SOUTH_PHI);
  emscripten::enum_<ParamSetResult>("ParamSetResult")
      .value("APPLIED", ParamSetResult::APPLIED)
      .value("NO_EFFECT", ParamSetResult::NO_EFFECT)
      .value("UNKNOWN_PARAM", ParamSetResult::UNKNOWN_PARAM)
      .value("READONLY", ParamSetResult::READONLY)
      .value("NON_FINITE", ParamSetResult::NON_FINITE)
      .value("INADMISSIBLE", ParamSetResult::INADMISSIBLE)
      .value("MALFORMED_PAYLOAD", ParamSetResult::MALFORMED_PAYLOAD)
      .value("TOO_LONG", ParamSetResult::TOO_LONG);

  emscripten::enum_<ClipSetResult>("ClipSetResult")
      .value("APPLIED", ClipSetResult::APPLIED)
      .value("NO_EFFECT", ClipSetResult::NO_EFFECT)
      .value("INVALID_BOUNDS", ClipSetResult::INVALID_BOUNDS)
      .value("FULL_FRAME_KEPT", ClipSetResult::FULL_FRAME_KEPT);

  emscripten::enum_<ResolutionSetResult>("ResolutionSetResult")
      .value("RESIZED", ResolutionSetResult::RESIZED)
      .value("ALREADY_ACTIVE", ResolutionSetResult::ALREADY_ACTIVE)
      .value("UNSUPPORTED", ResolutionSetResult::UNSUPPORTED);

  emscripten::enum_<EffectSetResult>("EffectSetResult")
      .value("INSTALLED", EffectSetResult::INSTALLED)
      .value("UNKNOWN_EFFECT", EffectSetResult::UNKNOWN_EFFECT)
      .value("UNSUPPORTED_RESOLUTION", EffectSetResult::UNSUPPORTED_RESOLUTION);

#if HS_ENABLE_CHAIN_INTERPRETER
  emscripten::enum_<Pullback::Interp::ChainStatus>("ChainStatus")
      .value("OK", Pullback::Interp::ChainStatus::OK)
      .value("NOT_CHAIN_EFFECT",
             Pullback::Interp::ChainStatus::NOT_CHAIN_EFFECT)
      .value("MALFORMED_PAYLOAD",
             Pullback::Interp::ChainStatus::MALFORMED_PAYLOAD)
      .value("EMPTY", Pullback::Interp::ChainStatus::EMPTY)
      .value("TOO_LONG", Pullback::Interp::ChainStatus::TOO_LONG)
      .value("UNKNOWN_OPERATOR",
             Pullback::Interp::ChainStatus::UNKNOWN_OPERATOR)
      .value("DUPLICATE_INSTANCE",
             Pullback::Interp::ChainStatus::DUPLICATE_INSTANCE)
      .value("MALFORMED_INSTANCE",
             Pullback::Interp::ChainStatus::MALFORMED_INSTANCE)
      .value("ENTRY_FAMILY", Pullback::Interp::ChainStatus::ENTRY_FAMILY)
      .value("EXIT_FAMILY", Pullback::Interp::ChainStatus::EXIT_FAMILY)
      .value("CARRIER_MISMATCH",
             Pullback::Interp::ChainStatus::CARRIER_MISMATCH)
      .value("ARENA_OVERFLOW", Pullback::Interp::ChainStatus::ARENA_OVERFLOW)
      .value("PARAM_OVERFLOW", Pullback::Interp::ChainStatus::PARAM_OVERFLOW)
      .value("MIGRATE_FAILED", Pullback::Interp::ChainStatus::MIGRATE_FAILED);
#endif

  emscripten::class_<HolosphereEngine>("HolosphereEngine")
      .constructor<>()
      .function("setResolution", &HolosphereEngine::setResolution)
      .function("setDisplayCaps", &HolosphereEngine::setDisplayCaps)
      .function("getDisplayNorthPhi", &HolosphereEngine::getDisplayNorthPhi)
      .function("getDisplaySouthPhi", &HolosphereEngine::getDisplaySouthPhi)
      .function("setEffect", &HolosphereEngine::setEffect)
      .function("drawFrame", &HolosphereEngine::drawFrame)
      .function("getPixels", &HolosphereEngine::getPixels)
      .function("getBufferLength", &HolosphereEngine::getBufferLength)
      .function("setParameter", &HolosphereEngine::setParameter)
      .function("setAnimationsPaused", &HolosphereEngine::setAnimationsPaused)
      .function("getAnimationsPaused", &HolosphereEngine::getAnimationsPaused)
      .function("getPresetCount", &HolosphereEngine::getPresetCount)
      .function("getPresetIndex", &HolosphereEngine::getPresetIndex)
      .function("getPresetIds", &HolosphereEngine::getPresetIds)
      .function("selectPreset", &HolosphereEngine::selectPreset)
      .function("selectPresetById", &HolosphereEngine::selectPresetById)
      .function("synchronizePreset", &HolosphereEngine::synchronizePreset)
      .function("nextPreset", &HolosphereEngine::nextPreset)
      .function("previousPreset", &HolosphereEngine::previousPreset)
      .function("setPoleLod", &HolosphereEngine::setPoleLod)
      .function("getPoleLod", &HolosphereEngine::getPoleLod)
      .function("getParameterDefinitions",
                &HolosphereEngine::getParameterDefinitions)
      .function("getParamValues", &HolosphereEngine::getParamValues)
      .function("getParamGeneration", &HolosphereEngine::getParamGeneration)
      .function("getArenaMetrics", &HolosphereEngine::getArenaMetrics)
      .function("getEffectSizes", &HolosphereEngine::getEffectSizes)
      .function("getEffectPresetCounts",
                &HolosphereEngine::getEffectPresetCounts)
#if HS_ENABLE_CHAIN_INTERPRETER
      .function("getShaderChainBindings",
                &HolosphereEngine::getShaderChainBindings)
#endif
      .class_function("getSupportedResolutions",
                      &HolosphereEngine::getSupportedResolutions)
      .class_function("isLive", &HolosphereEngine::isLive)
      .function("setClip", &HolosphereEngine::setClip)
      .function("strobeColumns", &HolosphereEngine::strobeColumns);
}
