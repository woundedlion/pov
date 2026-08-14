/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file engine_bindings.h
 * @brief JS-facing render bridge: the HolosphereEngine class and its embind
 *        registration.
 *
 * Owns the stack high-water instrumentation and the readback buffers JS renders
 * from, over the factory and HS_RESOLUTIONS dispatch in effect_factory.h.
 * Included only by targets/wasm/wasm.cpp.
 */
#pragma once

#include <emscripten/bind.h>
#include <emscripten/stack.h>
#include "core/engine/effects.h" // Includes all effect headers (triggers REGISTER_EFFECT)
#include "core/engine/effect_registry.h"
#include "core/engine/platform.h"
#include "targets/wasm/arena_metrics.h"
#include "targets/wasm/effect_factory.h" // pure, host-tested factory + dispatch
#include "targets/wasm/param_marshal.h"  // pure, host-tested param marshaling
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <algorithm> // std::fill_n — blank-frame clear in drawFrame
#include <cmath>
#include <string_view>
#include <cstring>
#include <climits> // INT_MAX — drawFrame pixel-index accumulator bound
#include <memory>
#include <string>

using namespace emscripten;

// ---- Stack canary painting for high water mark tracking ----
static constexpr uint8_t STACK_CANARY = 0xCD;

/**
 * @brief Paints the unused portion of the stack with a canary byte for high
 *        water mark tracking.
 * @details Writes STACK_CANARY from the current stack pointer down to the stack
 *          end. Safe to call at any time — only touches memory below the current
 *          stack pointer.
 */
static void stack_paint_canary() {
  uintptr_t sp = emscripten_stack_get_current();
  uintptr_t end = emscripten_stack_get_end();
  if (sp > end) {
    std::memset(reinterpret_cast<void *>(end), STACK_CANARY, sp - end);
  }
}

/**
 * @brief Computes the stack high water mark by scanning the canary region.
 * @return Number of stack bytes that have been touched, in bytes (the high
 *         water mark), found by scanning from the stack end upward to the first
 *         overwritten canary byte.
 * @details This is a LOWER bound, not an exact figure: a frame that wrote the
 *          coincidental byte 0xCD, or that reserved space it never actually
 *          stored to, reads back as still-canary and under-reports. Use it as a
 *          conservative "at least this deep" signal, not a precise measurement.
 *          The scan walks the whole still-canary region, so it runs word-wise
 *          and falls back to bytes only at the alignment head and the boundary
 *          word.
 */
static size_t stack_high_water_mark() {
  uintptr_t base = emscripten_stack_get_base();
  uintptr_t end = emscripten_stack_get_end();
  const uint8_t *p = reinterpret_cast<const uint8_t *>(end);
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

// Deepest stack an effect's construction + init() has reached, as a running max
// over every load. Latched because the very next repaint erases the canary
// evidence: without it the only readable mark is the render path's, and a
// stack-hungry constructor — what this instrumentation exists to catch — leaves
// nothing a gate can fail on.
static size_t init_stack_peak = 0;

// Upper bound on a single effect's exposed parameters, including effect-local
// fixed storage larger than ParamList's default inline array. Used
// to pre-reserve the getParamValues() backing store so it never reallocates.
static constexpr size_t MAX_PARAMS = 128;

static_assert(MAX_PARAMS >=
                  std::tuple_size<decltype(Effect::ParamList::elements)>::value,
              "MAX_PARAMS must cover ParamList's default array size");

/**
 * @brief Outcome of a HolosphereEngine::setClip() call.
 * @details Each outcome wants a different caller response, so they are separate
 *          enumerators rather than one bool: NO_EFFECT is the ordinary state
 *          after an init or setResolution that carried no effect name, while
 *          INVALID_BOUNDS is a caller bug. APPLIED and FULL_FRAME_KEPT are both
 *          successes, but only APPLIED means the band is in force — a segment
 *          pool needs the two apart to tell parallel speedup from N workers each
 *          computing the same full frame. Exposed to JS as the
 *          Module.ClipSetResult embind enum; compare against its values, never
 *          by truthiness (every enum value is a truthy object).
 */
enum class ClipSetResult {
  APPLIED,         /**< Band installed; rendering is narrowed to it. */
  NO_EFFECT,       /**< No effect is installed to receive the clip. */
  INVALID_BOUNDS,  /**< Bounds malformed or out of range for the resolution. */
  FULL_FRAME_KEPT, /**< Bounds accepted but ignored: the effect reports
                        Effect::needs_full_frame(), so the clip stays at the
                        full canvas. */
};

/**
 * @brief Outcome of a HolosphereEngine::setResolution() call.
 * @details Each outcome wants a different caller response, so they are separate
 *          enumerators rather than one bool: RESIZED tears the current effect
 *          down, so the caller must re-apply setEffect() — and any setClip() —
 *          before the next drawFrame(); ALREADY_ACTIVE tears nothing down, so
 *          the current effect, parameters, and clip all stay live; UNSUPPORTED
 *          keeps the prior valid state alive. Exposed to JS as the
 *          Module.ResolutionSetResult embind enum; compare against its values,
 *          never by truthiness (every enum value is a truthy object).
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
 * @details The two rejections leave the prior effect alive but demand different
 *          caller responses: UNKNOWN_EFFECT is a stale/typo'd name to drop from
 *          the UI, while UNSUPPORTED_RESOLUTION means no name can succeed until
 *          a supported setResolution() lands. Exposed to JS as the
 *          Module.EffectSetResult embind enum; compare against its values,
 *          never by truthiness (every enum value is a truthy object).
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
 * @brief Outcome of restoreFullConfigSnapshot().
 * @details Exposed to JS as the Module.FullConfigRestoreResult embind enum;
 *          compare against its values, never by truthiness. Every value but
 *          APPLIED leaves the effect exactly as it was.
 */
enum class FullConfigRestoreResult : uint8_t {
  APPLIED,             /**< Snapshot installed. */
  NOT_SHADERBALL,      /**< The loaded effect has no full configuration. */
  UNSUPPORTED_VERSION, /**< schemaVersion is not the current schema. */
  INVALID_LENGTH,      /**< Snapshot missing, or an array whose length is not
                            the field count. */
  INVALID_VALUE,       /**< A field or runtime value outside what its slot
                            admits. */
  INVALID_ACCEPTED,    /**< Fields each in range but a combination the effect
                            will not render. */
  INVALID_PENDING,     /**< Pending list is not a set of in-range field indices,
                            or does not match where accepted and requested
                            differ. */
};

/**
 * @brief JS-facing render engine driving one resolution/effect at a time.
 * @details Owns the current effect and the stable readback buffers. Every public
 *          method is exported to JavaScript via EMSCRIPTEN_BINDINGS below.
 *          At most one instance may be live: delete() the current engine before
 *          constructing another, or the constructor traps.
 */
class HolosphereEngine {
public:
  /**
   * @brief Constructs the engine with a valid default resolution and effect.
   * @details Pre-sizes the JS-facing readback buffers to their maximum extent so
   *          their backing storage never moves (the WASM memory-view contract),
   *          verifies the self-registering effect count against the static
   *          roster, and installs a default effect that JS overrides almost
   *          immediately.
   */
  HolosphereEngine() {
    HS_CHECK(!engine_alive,
             "HolosphereEngine is a singleton: delete() the live instance "
             "before constructing another (its Effect and arenas are shared "
             "module-global storage)");
    // The scan reads pole_lod_aggressiveness as a module global; claim it for
    // this instance so a fresh engine never inherits a predecessor's setting.
    pole_lod_aggressiveness = HS_POLE_LOD_DEFAULT;
    stack_paint_canary();

    // SSOT guard: the self-registering effect count must match the static roster
    // (HS_EFFECT_LIST / HS_EFFECT_COUNT) or the live set and the native smoke
    // suite silently diverge.
    HS_CHECK(EffectRegistry::entries().size() ==
                 static_cast<size_t>(HS_EFFECT_COUNT),
             "each effect header must be included by exactly one translation "
             "unit");

    // Pre-size the view-backed readback buffers ONCE: under ALLOW_MEMORY_GROWTH
    // a reallocation detaches the ArrayBuffer behind a typed_memory_view, so the
    // buffers returned as views (getPixels/getParamValues) must never move.
    pixel_buffer.assign(MAX_W * MAX_H * CHANNELS, 0);
    param_values.reserve(MAX_PARAMS);
    param_views.reserve(
        MAX_PARAMS); // not view-backed; reserve is amortization only

    // Bootstrap on the first row of each roster; daydream overrides both almost
    // immediately.
    const bool bootstrap_resolution_set =
        setResolution(WASM_RESOLUTIONS[0].w, WASM_RESOLUTIONS[0].h) ==
        ResolutionSetResult::RESIZED;
    HS_CHECK(bootstrap_resolution_set,
             "the first HS_RESOLUTIONS row must be dispatchable here");
    const bool bootstrap_effect_set =
        setEffect(WASM_EFFECT_NAMES[0]) == EffectSetResult::INSTALLED;
    HS_CHECK(bootstrap_effect_set,
             "the first HS_EFFECT_LIST entry must be registered and buildable "
             "at the first HS_RESOLUTIONS row");

    // Last: a bootstrap trap the JS caller catches must not leave the singleton
    // flag latched against every later construction.
    engine_alive = true;
  }

  /**
   * @brief Destroys the engine and admits the next construction.
   * @details JS reaches this through the embind-generated delete().
   */
  ~HolosphereEngine() { engine_alive = false; }

  /**
   * @brief Switches the active canvas resolution.
   * @param w Requested canvas width in pixels.
   * @param h Requested canvas height in pixels.
   * @return RESIZED if the resolution switched, ALREADY_ACTIVE if the request
   *         matched the active resolution (a pure no-op — nothing is torn
   *         down), or UNSUPPORTED if the request was rejected and the previous
   *         valid state was kept. Exposed to JS as the
   *         Module.ResolutionSetResult embind enum; compare against its
   *         values, never by truthiness.
   * @details A RESIZED result tears down the current effect (a new one cannot
   *          be carried across pixel dimensions), so the caller must call
   *          setEffect() again before the next drawFrame() or it renders a
   *          blank frame. Callers that use a sub-canvas clip (segmented
   *          rendering) must likewise re-apply setClip() after a RESIZED
   *          setResolution(): a prior clip was expressed in the old
   *          resolution's pixel bounds and is not rescaled here, so it must be
   *          recomputed for the new dimensions. The teardown re-partitions the
   *          engine arenas, so getArenaMetrics() reports the released state
   *          rather than the destroyed effect's usage. An outstanding
   *          getPixels() view is left aliasing live memory at the previous
   *          resolution's length rather than detached, so it must be re-fetched
   *          and its length compared against getBufferLength().
   */
  ResolutionSetResult setResolution(int w, int h) {
    if (w == pixel_width && h == pixel_height)
      return ResolutionSetResult::ALREADY_ACTIVE;

    // Reject unsupported sizes and keep the prior valid state alive rather than
    // switching to a null effect that renders blank with no signal to JS.
    if (!wasm_resolution_supported(w, h)) {
      hs::log("WASM: Unsupported resolution %dx%d — ignored", w, h);
      return ResolutionSetResult::UNSUPPORTED;
    }

    pixel_width = w;
    pixel_height = h;
    const int count = pixel_width * pixel_height * CHANNELS;
    std::fill_n(pixel_buffer.data(), count, uint16_t{0});

    if (current_effect) {
      current_effect = nullptr;
      current_effect_type_key = nullptr;
      // Same teardown setEffect() performs: without the re-partition the
      // destroyed effect's arena usage keeps reading as live from
      // getArenaMetrics().
      configure_arenas_default();
      param_generation.replace(0);
      stack_paint_canary(); // repaint to reset stack HWM after teardown
    }
    return ResolutionSetResult::RESIZED;
  }

  /**
   * @brief Tears down the current effect and instantiates the named one at the
   *        active resolution.
   * @param name Effect name to instantiate.
   * @return INSTALLED iff an effect was actually instantiated, else the
   *         rejection reason — UNKNOWN_EFFECT for an unknown/stale effect name
   *         or UNSUPPORTED_RESOLUTION — so the frontend can detect a no-op
   *         instead of believing the switch took. Exposed to JS as the
   *         Module.EffectSetResult embind enum; compare against its values,
   *         never by truthiness.
   * @details Validates the name against the factory for the current resolution
   *          BEFORE tearing anything down, so a typo'd/stale UI string keeps the
   *          prior valid state alive rather than blanking the engine. An
   *          INSTALLED switch instantiates a fresh effect at the full-canvas
   *          clip and with default parameter values, so callers that use a
   *          sub-canvas clip (segmented rendering) or tuned parameters must
   *          re-apply setClip() and setParameter() after every INSTALLED
   *          setEffect(). The engine-owned animation pause state is retained.
   */
  EffectSetResult setEffect(const std::string &name) {
    // Validate against the current resolution's factory BEFORE tearing anything
    // down, so a typo'd name keeps the prior valid state alive.
    const FactoryEntry *entry = nullptr;
    const bool dispatched =
        dispatch_resolution(pixel_width, pixel_height, [&]<int W, int H>() {
          entry = find_factory_entry<W, H>(name);
        });
    if (!dispatched) {
      hs::log("WASM: setEffect at unsupported resolution %dx%d — keeping "
              "current effect",
              pixel_width, pixel_height);
      return EffectSetResult::UNSUPPORTED_RESOLUTION;
    }

    hs::log("WASM: setEffect called with %s", name.c_str());
    if (!entry) {
      hs::log("WASM: setEffect unknown effect '%s' — keeping current effect",
              name.c_str());
      return EffectSetResult::UNKNOWN_EFFECT;
    }

    // Per-load RNG stream, mirroring the device's per-epoch reseed: load 0
    // (the constructor's bootstrap effect) keeps the identity 1337 stream;
    // each later load derives a fresh seed. Replica-safe with no protocol
    // change — workers process identical message streams, so their counters
    // agree by construction.
    hs::random().seed(hs::epoch_seed(effect_loads++));

    current_effect.reset();
    current_effect_type_key = nullptr;
    configure_arenas_default(); // Reset before init so effects can override

    stack_paint_canary(); // reset stack HWM by repainting unused region

    dispatch_resolution(pixel_width, pixel_height, []<int W, int H>() {
      init_geometry_luts<W, H>(); // eager-fill LUTs before the first frame
    });
    current_effect = entry->creator();
    current_effect_type_key = entry->type_key;
    current_effect->setAnimationsPaused(animations_paused);
    current_effect->init();
    param_generation.replace(current_effect->getParameterSchemaGeneration());
    const size_t init_hwm = stack_high_water_mark();
    if (init_hwm > init_stack_peak)
      init_stack_peak = init_hwm;
    hs::log("WASM: init stack HWM = %u bytes", (unsigned)init_hwm);
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
   *         were accepted but the effect reports Effect::needs_full_frame() and
   *         so keeps the full-canvas clip, otherwise the rejection reason:
   *         NO_EFFECT (no effect is set to receive the clip) or INVALID_BOUNDS
   *         (malformed/out of range, then ignored). Exposed to JS as the
   *         Module.ClipSetResult embind enum; compare against its values, never
   *         by truthiness. Both APPLIED and FULL_FRAME_KEPT are successes.
   *         NO_EFFECT is the ordinary answer between a resolution change and the
   *         setEffect that follows it, so a caller that faults on a rejection
   *         must fault on INVALID_BOUNDS only.
   * @details Args are x-pair-first to match the (x, y) convention: embind binds
   *          positionally, so a y-first order would let a transposed JS call pass
   *          the range check and silently clip the wrong axis. Rejects malformed
   *          input at the untyped JS boundary rather than trapping, since a trap
   *          there aborts the whole WASM module. Bounds arrive as doubles: an
   *          i32 embind parameter coerces NaN and out-of-range JS numbers to 0
   *          with no range check in a release build, installing an empty band
   *          under an APPLIED result. Segment workers always pass
   *          valid, ordered, in-range bounds. A cross-segment stateful effect
   *          (Effect::needs_full_frame()) keeps the full-canvas clip instead of
   *          narrowing to the band — see docs/segmented_stateful_effects_spec.md.
   */
  ClipSetResult setClip(double x0, double x1, double y0, double y1) {
    if (!current_effect)
      return ClipSetResult::NO_EFFECT;
    // Reject malformed bounds from the untyped JS boundary (negatives would feed
    // ClipRegion's modulo arithmetic); reject-and-return, never trap.
    if (!hs_wasm::clip_bounds_valid(x0, x1, y0, y1, pixel_width,
                                    pixel_height)) {
      hs::log("WASM: setClip bounds out of range (x0=%g,x1=%g,y0=%g,y1=%g) — "
              "ignored",
              x0, x1, y0, y1);
      return ClipSetResult::INVALID_BOUNDS;
    }
    // Cross-segment stateful effects must render the FULL canvas in every worker
    // (a band-clipped worker has stale cv.prev outside its band, so trails seam);
    // keep the full clip. See docs/segmented_stateful_effects_spec.md.
    if (current_effect->needs_full_frame())
      return ClipSetResult::FULL_FRAME_KEPT;
    current_effect->set_clip(static_cast<int>(y0), static_cast<int>(y1),
                             static_cast<int>(x0), static_cast<int>(x1));
    return ClipSetResult::APPLIED;
  }

  /**
   * @brief Renders one frame of the current effect into the JS-facing buffer.
   * @details Copies the effect's full canvas into pixel_buffer as 16-bit linear
   *          RGB triples; no-op if no effect is set. The readback copies the
   *          FULL canvas regardless of any active clip region — a clip restricts
   *          rendering, not this readback.
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

    // Readback copies the FULL canvas regardless of any clip; segment_worker.js
    // extracts its quadrant JS-side (README §10.7).
    static_assert(static_cast<long long>(MAX_W) * MAX_H * CHANNELS <= INT_MAX,
                  "drawFrame pixel-index accumulators are int");
    const int count = pixel_width * pixel_height;
    if (!current_effect->overrides_get_pixel()) {
      // Fast path: display_buffer()[i] == get_pixel(x, y), so copy directly.
      const Pixel *buf = current_effect->display_buffer();
      static_assert(sizeof(Pixel) == 3 * sizeof(uint16_t),
                    "fast-path memcpy assumes packed RGB16 Pixel layout");
      std::memcpy(pixel_buffer.data(), buf,
                  static_cast<size_t>(count) * sizeof(Pixel));
    } else {
      int idx = 0;
      for (int y = 0; y < pixel_height; y++) {
        for (int x = 0; x < pixel_width; x++) {
          const Pixel &p = current_effect->get_pixel(x, y);
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
   * @details The simulator's dot mesh inherently renders discrete columns with
   *          gaps — already the strobe == true look — so daydream reads this to
   *          decide whether to fill the inter-column gap (false) or leave it
   *          dark (true). See docs/strobe_columns_audit.md.
   */
  bool strobeColumns() const {
    return current_effect ? current_effect->strobe_columns() : false;
  }

  /**
   * @brief Exposes the raw pixel buffer to JS as a zero-copy Uint16Array view.
   * @return Typed memory view over the active resolution's R,G,B pixels within
   *         the stable MAX_W*MAX_H*3 backing buffer.
   * @details WASM memory-view contract: the returned view aliases WASM linear
   *          memory, it is NOT a copy, and two independent events invalidate it.
   *          With ALLOW_MEMORY_GROWTH=1, any subsequent heap growth detaches the
   *          underlying ArrayBuffer and leaves this view zero-length
   *          (buffer.byteLength === 0). A successful setResolution() detaches
   *          nothing — the backing vector is pre-sized once and never
   *          reallocated — but moves the active prefix this view spans, so an
   *          outstanding view keeps aliasing live memory at the wrong length. A
   *          caller caching the view across frames MUST therefore test both:
   *          buffer.byteLength !== 0 and length === getBufferLength().
   */
  val getPixels() {
    return val(typed_memory_view(pixel_width * pixel_height * CHANNELS,
                                 pixel_buffer.data()));
  }

  /**
   * @brief Returns the length of the active pixel buffer view.
   * @return Number of uint16 elements in the active view (pixel_width *
   *         pixel_height * 3, three channels per pixel).
   * @details The staleness test for a cached getPixels() view: a resolution
   *          change moves this length without detaching the view, so a held view
   *          whose length differs describes a prior resolution.
   */
  int getBufferLength() const { return pixel_width * pixel_height * CHANNELS; }

  /**
   * @brief Updates one named effect parameter.
   * @param name Parameter name to update.
   * @param value New parameter value, in the parameter's native units.
   * @return APPLIED if the write was accepted, otherwise the rejection reason:
   *         NO_EFFECT (no effect is set), UNKNOWN_PARAM (the name is unknown to
   *         the effect), READONLY (engine-written telemetry the GUI must not
   *         poke), or NON_FINITE (rejected before it can poison render math).
   *         Exposed to JS as the Module.ParamSetResult embind enum; compare
   *         against its values, never by truthiness (every enum value is a
   *         truthy object). An APPLIED float is silently clamped to the
   *         parameter's registered [min,max] — APPLIED does NOT imply the
   *         stored value equals the requested one. A consumer that needs the
   *         effective value should read it back via getParamValues().
   * @details An APPLIED write to an *animated* param engages the engine's
   *          animation pause, exactly as setAnimationsPaused(true) would: the
   *          value the caller just set would otherwise be overwritten by the
   *          animation on the next frame. The pause is retained across
   *          setEffect(), so a caller mirroring it must read it back through
   *          getAnimationsPaused() rather than track it separately.
   */
  ParamSetResult setParameter(const std::string &name, float value) {
    if (!current_effect)
      return ParamSetResult::NO_EFFECT;
    // The rejection classification is single-sourced in Effect::updateParameter,
    // not re-checked here.
    const ParamSetResult result =
        current_effect->updateParameter(name.c_str(), value);
    animations_paused = current_effect->animations_paused();
    return result;
  }

  /**
   * @brief Pauses or resumes the current effect's parameter animations.
   * @param paused true to pause animations, false to resume.
   * @details Retained across effect and resolution changes and applied to the
   *          next effect when no effect is currently loaded. setParameter()
   *          engages the same pause on an animated write, so this is not the
   *          only writer.
   */
  void setAnimationsPaused(bool paused) {
    animations_paused = paused;
    if (current_effect)
      current_effect->setAnimationsPaused(paused);
  }

  /**
   * @brief Reports whether the engine's parameter animations are paused.
   * @return true while animations are frozen.
   * @details The state setAnimationsPaused() writes and an animated
   *          setParameter() write engages, so a caller reads the rule here
   *          instead of re-implementing it.
   */
  bool getAnimationsPaused() const { return animations_paused; }

  /**
   * @brief Number of presets the current effect exposes for manual navigation.
   * @return The preset count, 0 when no effect is set or the effect authored
   *         none. A GUI reads this to decide whether to offer preset controls
   *         at all.
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
   * @details Tracks engine-driven advancement as well as the calls below: an
   *          effect whose choreography advances its own presets moves this
   *          without any JS call, so a mirroring GUI polls it rather than
   *          tracking the index it last wrote.
   */
  uint32_t getPresetIndex() const {
    return current_effect
               ? static_cast<uint32_t>(current_effect->getPresetIndex())
               : 0;
  }

  /**
   * @brief Selects one preset and freezes the effect's parameter animations.
   * @param index Preset to select, in [0, getPresetCount()).
   * @return true when the preset was applied; false when no effect is set, the
   *         index is out of range, or the effect refused the preset.
   * @details Engages the animation pause exactly as setAnimationsPaused(true)
   *          would, so the selected preset's values are not immediately
   *          overwritten by the animation — the manual-navigation counterpart
   *          to synchronizePreset(), which leaves the pause alone. The pause is
   *          retained across setEffect(), so a caller mirroring it reads it
   *          back through getAnimationsPaused(). Parameter values move with the
   *          preset; re-read them via getParamValues().
   */
  bool selectPreset(uint32_t index) {
    if (!current_effect || !current_effect->selectPreset(index))
      return false;
    animations_paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Selects one preset without touching the animation pause state.
   * @param index Preset to select, in [0, getPresetCount()).
   * @return true when the preset is the active one already or was applied;
   *         false when no effect is set, the index is out of range, or the
   *         effect refused the preset.
   * @details The call a GUI uses to follow engine-driven preset advancement:
   *          selectPreset() would freeze the choreography it is trying to
   *          mirror, this one leaves it running. A request for the active index
   *          is a success no-op.
   */
  bool synchronizePreset(uint32_t index) {
    if (!current_effect || !current_effect->synchronizePreset(index))
      return false;
    animations_paused = current_effect->animations_paused();
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
    animations_paused = current_effect->animations_paused();
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
    animations_paused = current_effect->animations_paused();
    return true;
  }

  /**
   * @brief Sets near-pole azimuthal shading decimation.
   * @param aggressiveness Columns per shade are this over sin(colatitude);
   *        0 disables. Non-finite and negative inputs clamp to 0.
   * @details The physically-neutral setting is 1.0: at that value one shade
   *          covers the columns sharing a physical LED footprint. Exposed so
   *          the value can be tuned against real hardware. The setting is
   *          engine-scoped: the constructor restores HS_POLE_LOD_DEFAULT, and
   *          each WASM instance carries its own, so a segmented pool must
   *          re-send it to every worker (README §10.7).
   */
  void setPoleLod(float aggressiveness) {
    pole_lod_aggressiveness =
        hs_wasm::clamp_pole_lod_aggressiveness(aggressiveness);
  }

  /**
   * @brief Current near-pole azimuthal decimation aggressiveness.
   * @return The clamped value of the last setPoleLod() on this engine, else
   *         HS_POLE_LOD_DEFAULT.
   */
  float getPoleLod() const { return pole_lod_aggressiveness; }

  /**
   * @brief Builds the GUI's parameter descriptor list.
   * @return JS array with one {name, value, animated, readonly, preset} object
   *         per param in the effect's declaration order, plus {min, max} on
   *         every non-boolean param, {step} on every whole-number param, and
   *         {options} — with {exportOptions} alongside it when the param
   *         declares C++ enum literals — on every enum param; empty array when
   *         no effect is set.
   * @details `value` is the current rendered state for GUI display, while
   *          `requestedValue` is the writable target used to seed another
   *          renderer. A boolean param's values are JS booleans and it carries
   *          no range; every other value is a number. step is 1 on an enum or
   *          integer target and absent on a float one, so the GUI knows which
   *          controls admit only whole values. An enum's value indexes its
   *          options array; an integer param carries a range instead of labels
   *          and exports as a plain numeric literal. preset marks the params a
   *          preset export carries.
   *          The order matches getParamValues(); pin getParamGeneration() beside
   *          a snapshot to detect a rebind.
   */
  val getParameterDefinitions() {
    if (!current_effect)
      return val::array();

    current_effect->refresh_parameter_display();
    check_param_capacity();
    val result = val::array();
    // Both streams walk the effect's registered ParamList in order. The
    // generation token covers effect replacement and dynamic schema rebinds.
    hs_wasm::collect_param_views(*current_effect, param_views);

    int i = 0;
    for (const auto &v : param_views) {
      val entry = val::object();
      entry.set("name", val(v.name));

      if (v.is_bool) {
        // Emit a JS boolean so the frontend renders a checkbox; toggles omit
        // min/max (no range).
        entry.set("value", val(v.value > 0.5f));
        entry.set("requestedValue", val(v.requested_value > 0.5f));
        entry.set("acceptedValue", val(v.accepted_value > 0.5f));
      } else {
        entry.set("value", v.value);
        entry.set("requestedValue", v.requested_value);
        entry.set("acceptedValue", v.accepted_value);
        entry.set("min", v.min);
        entry.set("max", v.max);
        if (v.is_integer)
          entry.set("step", 1);
        if (v.option_count > 0) {
          // Enum: label array indexed by value; the frontend renders a dropdown.
          val opts = val::array();
          for (int k = 0; k < v.option_count; ++k)
            opts.set(k, val(v.options[k]));
          entry.set("options", opts);
          if (v.export_options != nullptr) {
            val export_opts = val::array();
            for (int k = 0; k < v.option_count; ++k)
              export_opts.set(k, val(v.export_options[k]));
            entry.set("exportOptions", export_opts);
          }
        }
      }
      entry.set("animated", val(v.animated));
      entry.set("readonly", val(v.readonly));
      entry.set("preset", val(v.preset));
      if (const char *warning = current_effect->parameter_warning(v.name))
        entry.set("warning", val(warning));
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
   *          memory and must be consumed before the next allocation. param_values
   *          never reallocates here (size <= MAX_PARAMS), so emitting it triggers
   *          no heap growth that could detach other outstanding views.
   */
  val getParamValues() {
    if (!current_effect) {
      // Empty Float32Array (not a JS Array) so callers get a consistent typed
      // view whether or not an effect is set. clear() retains the ctor's
      // reserve, which keeps the zero-length view's backing pointer valid.
      param_values.clear();
      return val(typed_memory_view(param_values.size(), param_values.data()));
    }

    current_effect->refresh_parameter_display();
    // The no-reallocation contract the view rides on: the ctor reserved
    // MAX_PARAMS and clear() retains that capacity, so a fill within that bound
    // never moves the backing store.
    check_param_capacity();

    // Same order as getParameterDefinitions().
    hs_wasm::fill_param_values(*current_effect, param_values);
    return val(typed_memory_view(param_values.size(), param_values.data()));
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
   * @return JS object of the three engine arenas' metrics ({usage,
   *         high_water_mark, capacity}) plus a "stack" entry
   *         ({high_water_mark, init_high_water_mark, capacity}), all in bytes.
   * @details Read once per frame by the HUD, on the main thread and in every
   *          segment worker. The tooling arenas are not included: an engine
   *          instance never moves them, and MeshOps.getArenaMetrics() reports
   *          all six on demand. `high_water_mark` is the canary's live reading,
   *          which a repaint resets and the render path then dominates;
   *          `init_high_water_mark` is the latched deepest effect construction +
   *          init(), which no repaint erases, so the two peaks can be gated
   *          separately.
   */
  val getArenaMetrics() {
    val metrics = collect_engine_arena_metrics();

    // Stack region. No running usage: the live depth at this call is outside any
    // render, so only the high-water marks are meaningful.
    {
      uintptr_t base = emscripten_stack_get_base();
      uintptr_t end = emscripten_stack_get_end();
      val m = val::object();
      m.set("high_water_mark", static_cast<size_t>(stack_high_water_mark()));
      m.set("init_high_water_mark", init_stack_peak);
      m.set("capacity", static_cast<size_t>(base >= end ? base - end : 0));
      metrics.set("stack", m);
    }

    return metrics;
  }

  /**
   * @brief Builds the {effect name -> hint size} map for the (W,H) factory.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @return JS object mapping each effect name to its hint size, for the GUI.
   */
  template <int W, int H> val get_effect_sizes_helper() {
    val s = val::object();
    const auto &factory = get_factory<W, H>();
    for (const auto &entry : factory)
      s.set(std::string(entry.name), static_cast<int>(entry.size));
    return s;
  }

  /**
   * @brief Returns the effect-size map for the active resolution.
   * @return JS object mapping each effect name to its hint size at the current
   *         resolution; empty map if unsupported/uninitialized.
   */
  val getEffectSizes() {
    val sizes = val::object(); // unsupported/uninitialized — empty map
    dispatch_resolution(pixel_width, pixel_height, [&]<int W, int H>() {
      sizes = get_effect_sizes_helper<W, H>();
    });
    return sizes;
  }

  /**
   * @brief Returns the current ShaderBall's versioned full-state snapshot.
   * @return JS object {schemaVersion, accepted, requested, pendingFieldIds,
   *         hasRuntime, runtime}, or null when the loaded effect is not
   *         ShaderBall.
   * @details `accepted` and `requested` are CONFIG_FIELD_COUNT-long arrays of
   *          uint32-encoded field values in ConfigFieldId order;
   *          getFullConfigFieldDefinitions() names the indices.
   *          `pendingFieldIds` lists the fields carrying an unresolved edit,
   *          which at the current schema are exactly the fields where the two
   *          arrays differ. `runtime` is the animation clock state and is
   *          meaningful only when `hasRuntime`. The whole configuration crosses
   *          as one object because ShaderBall vets slots and params together:
   *          replaying the parameter stream entry by entry walks through
   *          combinations it refuses.
   */
  val getFullConfigSnapshot() {
    val output = val::null();
    with_shaderball([&]<typename SB>(SB &shaderball) {
      const typename SB::FullConfigSnapshot snapshot =
          shaderball.capture_full_config_snapshot();
      output = val::object();
      output.set("schemaVersion", snapshot.schema_version);
      output.set("accepted", uint32_array(snapshot.accepted));
      output.set("requested", uint32_array(snapshot.requested));
      val pending_ids = val::array();
      size_t pending_count = 0;
      for (size_t index = 0; index < snapshot.pending.size(); ++index)
        if (snapshot.pending[index] != 0)
          pending_ids.set(pending_count++, index);
      output.set("pendingFieldIds", pending_ids);
      output.set("hasRuntime", snapshot.has_runtime);
      output.set("runtime", float_array(snapshot.runtime));
    });
    return output;
  }

  /**
   * @brief Atomically restores a current ShaderBall full-state snapshot.
   * @param input Object in getFullConfigSnapshot()'s shape.
   * @return APPLIED, or the reason the snapshot was refused.
   * @details Rejections leave the effect exactly as it was, so a failed restore
   *          needs no rollback. NOT_SHADERBALL covers the loaded effect;
   *          UNSUPPORTED_VERSION a schemaVersion other than the current one;
   *          INVALID_LENGTH a missing snapshot or an array whose length is not
   *          the field count; INVALID_VALUE a field or runtime value outside
   *          what its slot admits; INVALID_ACCEPTED fields each in range but a
   *          combination the effect will not render; INVALID_PENDING a pending
   *          list that is not a set of in-range field indices, or one that does
   *          not match where accepted and requested differ.
   */
  FullConfigRestoreResult restoreFullConfigSnapshot(const val &input) {
    FullConfigRestoreResult result = FullConfigRestoreResult::NOT_SHADERBALL;
    with_shaderball([&]<typename SB>(SB &shaderball) {
      typename SB::FullConfigSnapshot snapshot;
      if (input.isUndefined() || input.isNull()) {
        result = FullConfigRestoreResult::INVALID_LENGTH;
        return;
      }
      const val schema_version = input["schemaVersion"];
      if (!schema_version.isNumber()) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      const double schema_number = schema_version.as<double>();
      if (!whole_uint32(schema_number)) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      snapshot.schema_version = static_cast<uint32_t>(schema_number);
      if (!SB::config_version_supported(snapshot.schema_version)) {
        result = FullConfigRestoreResult::UNSUPPORTED_VERSION;
        return;
      }
      ArrayDecode decoded =
          decode_uint32_array(input["accepted"], snapshot.accepted);
      if (decoded == ArrayDecode::OK)
        decoded = decode_uint32_array(input["requested"], snapshot.requested);
      if (decoded == ArrayDecode::OK)
        decoded = decode_runtime(input["runtime"], snapshot.runtime);
      if (decoded != ArrayDecode::OK) {
        result = decoded == ArrayDecode::BAD_LENGTH
                     ? FullConfigRestoreResult::INVALID_LENGTH
                     : FullConfigRestoreResult::INVALID_VALUE;
        return;
      }
      const val has_runtime = input["hasRuntime"];
      if (!has_runtime.isTrue() && !has_runtime.isFalse()) {
        result = FullConfigRestoreResult::INVALID_VALUE;
        return;
      }
      snapshot.has_runtime = has_runtime.as<bool>();
      const val pending_ids = input["pendingFieldIds"];
      if (!is_array(pending_ids)) {
        result = FullConfigRestoreResult::INVALID_LENGTH;
        return;
      }
      const size_t pending_count = pending_ids["length"].as<size_t>();
      if (pending_count > snapshot.pending.size()) {
        result = FullConfigRestoreResult::INVALID_PENDING;
        return;
      }
      for (size_t index = 0; index < pending_count; ++index) {
        const val field_id = pending_ids[index];
        if (!field_id.isNumber()) {
          result = FullConfigRestoreResult::INVALID_PENDING;
          return;
        }
        const double field_number = field_id.as<double>();
        if (!whole_uint32(field_number) ||
            field_number >= snapshot.pending.size()) {
          result = FullConfigRestoreResult::INVALID_PENDING;
          return;
        }
        uint8_t &pending = snapshot.pending[static_cast<size_t>(field_number)];
        if (pending != 0) {
          result = FullConfigRestoreResult::INVALID_PENDING;
          return;
        }
        pending = 1;
      }
      result = map_restore_result<SB>(
          shaderball.restore_full_config_snapshot(snapshot));
    });
    return result;
  }

  /**
   * @brief Returns stable ShaderBall field IDs and diagnostic names.
   * @return JS array of {id, name} in ConfigFieldId order, or null when the
   *         loaded effect is not ShaderBall.
   * @details `id` is the index into a snapshot's accepted/requested arrays and
   *          the value pendingFieldIds carries; `name` is the field's dotted
   *          config path. A caller labels a field through this rather than a
   *          hardcoded index, which moves when the schema gains a field.
   */
  val getFullConfigFieldDefinitions() {
    val output = val::null();
    with_shaderball([&]<typename SB>(SB &) {
      output = val::array();
      for (size_t index = 0; index < SB::CONFIG_FIELD_COUNT; ++index) {
        val field = val::object();
        field.set("id", index);
        field.set("name", SB::config_field_name(
                              static_cast<typename SB::ConfigFieldId>(index)));
        output.set(index, field);
      }
    });
    return output;
  }

  /**
   * @brief Reserved compatibility accessor.
   * @return Empty string.
   */
  std::string getConfigImportNotice() {
    std::string notice;
    with_shaderball([&]<typename SB>(SB &shaderball) {
      notice = shaderball.config_import_notice();
    });
    return notice;
  }

  /**
   * @brief Reserved compatibility no-op.
   */
  void clearConfigImportNotice() {
    with_shaderball([&]<typename SB>(SB &shaderball) {
      shaderball.clear_config_import_notice();
    });
  }

  /**
   * @brief Enumerates the resolutions the factory can build.
   * @return JS array of [W, H] pairs, generated from the same
   *         HS_RESOLUTIONS list that setResolution()/getEffectSizes()
   *         dispatch through.
   * @details Callers (e.g. the CI smoke test) can enumerate this instead of
   *          hand-mirroring the list, so the supported set can never silently
   *          drift.
   */
  static val getSupportedResolutions() {
    val out = val::array();
    int i = 0;
#define X(W, H)                                                                \
  {                                                                            \
    val pair = val::array();                                                   \
    pair.set(0, (W));                                                          \
    pair.set(1, (H));                                                          \
    out.set(i++, pair);                                                        \
  }
    HS_RESOLUTIONS(X)
#undef X
    return out;
  }

private:
  /**
   * @brief Runs a callback on the live effect iff it is a ShaderBall.
   * @param callback Templated callable invoked as callback(ShaderBall<W,H>&).
   * @return true iff the callback ran.
   * @details The downcast is admitted by the factory's own type key, so it is
   *          legal by construction: a registry name that no longer maps to
   *          ShaderBall<W,H> reports "not a ShaderBall" instead of casting to
   *          the wrong type.
   */
  template <typename Callback> bool with_shaderball(Callback &&callback) {
    if (!current_effect)
      return false;
    bool invoked = false;
    dispatch_resolution(pixel_width, pixel_height, [&]<int W, int H>() {
      if (current_effect_type_key != effect_type_key<ShaderBall<W, H>>())
        return;
      callback(*static_cast<ShaderBall<W, H> *>(current_effect.get()));
      invoked = true;
    });
    return invoked;
  }

  /** @brief Why a snapshot array failed to decode. */
  enum class ArrayDecode : uint8_t {
    OK,
    BAD_LENGTH, /**< Not an array, or not the field count long. */
    BAD_VALUE,  /**< An element outside what its slot admits. */
  };

  static bool is_array(const val &value) {
    return val::global("Array").call<bool>("isArray", value);
  }

  static bool whole_uint32(double value) {
    return std::isfinite(value) && value >= 0.0 &&
           value <= static_cast<double>(UINT32_MAX) &&
           value == std::floor(value);
  }

  template <size_t N>
  static val uint32_array(const std::array<uint32_t, N> &values) {
    val output = val::array();
    for (size_t index = 0; index < N; ++index)
      output.set(index, values[index]);
    return output;
  }

  template <size_t N>
  static val float_array(const std::array<float, N> &values) {
    val output = val::array();
    for (size_t index = 0; index < N; ++index)
      output.set(index, values[index]);
    return output;
  }

  template <size_t N>
  static ArrayDecode decode_uint32_array(const val &input,
                                         std::array<uint32_t, N> &output) {
    if (!is_array(input) || input["length"].as<size_t>() != N)
      return ArrayDecode::BAD_LENGTH;
    for (size_t index = 0; index < N; ++index) {
      const val element = input[index];
      if (!element.isNumber())
        return ArrayDecode::BAD_VALUE;
      const double number = element.as<double>();
      if (!whole_uint32(number))
        return ArrayDecode::BAD_VALUE;
      output[index] = static_cast<uint32_t>(number);
    }
    return ArrayDecode::OK;
  }

  template <size_t N>
  static ArrayDecode decode_runtime(const val &input,
                                    std::array<float, N> &output) {
    if (!is_array(input) || input["length"].as<size_t>() != N)
      return ArrayDecode::BAD_LENGTH;
    for (size_t index = 0; index < N; ++index) {
      const val element = input[index];
      if (!element.isNumber())
        return ArrayDecode::BAD_VALUE;
      output[index] = element.as<float>();
      if (!std::isfinite(output[index]))
        return ArrayDecode::BAD_VALUE;
    }
    return ArrayDecode::OK;
  }

  template <typename SB>
  static FullConfigRestoreResult
  map_restore_result(typename SB::ConfigRestoreResult result) {
    switch (result) {
    case SB::ConfigRestoreResult::APPLIED:
      return FullConfigRestoreResult::APPLIED;
    case SB::ConfigRestoreResult::UNSUPPORTED_VERSION:
      return FullConfigRestoreResult::UNSUPPORTED_VERSION;
    case SB::ConfigRestoreResult::INVALID_LENGTH:
      return FullConfigRestoreResult::INVALID_LENGTH;
    case SB::ConfigRestoreResult::INVALID_VALUE:
      return FullConfigRestoreResult::INVALID_VALUE;
    case SB::ConfigRestoreResult::INVALID_ACCEPTED:
      return FullConfigRestoreResult::INVALID_ACCEPTED;
    case SB::ConfigRestoreResult::INVALID_PENDING:
      return FullConfigRestoreResult::INVALID_PENDING;
    }
    __builtin_unreachable();
  }

  /**
   * @brief Traps when the live effect exposes more params than were reserved.
   * @details Both param streams ride the constructor's MAX_PARAMS reserve, so
   *          both check it and neither builds a partial stream past it.
   *          use_parameter_storage() lets an effect exceed ParamList's default
   *          inline array, so the live count is the bound.
   */
  void check_param_capacity() const {
    HS_CHECK(current_effect->getParameters().size() <= MAX_PARAMS,
             "effect exposes %zu params, past the %zu reserved",
             current_effect->getParameters().size(), MAX_PARAMS);
  }

  /** Channels per pixel in the readback buffer (linear RGB triples). */
  static constexpr int CHANNELS = 3;

  std::unique_ptr<Effect>
      current_effect; /**< Currently active effect, or null. */
  const void *current_effect_type_key =
      nullptr; /**< Concrete type the factory built current_effect as; null when
                    there is no effect. Gates every downcast off the base. */
  std::vector<uint16_t> pixel_buffer; /**< 16-bit linear RGB readback buffer. */
  std::vector<float> param_values;    /**< Backing store for getParamValues. */
  std::vector<hs_wasm::ParamView>
      param_views;           /**< Scratch for getParameterDefinitions. */
  int pixel_width = 0;       /**< Active canvas width in pixels. */
  int pixel_height = 0;      /**< Active canvas height in pixels. */
  uint32_t effect_loads = 0; /**< Effect loads so far; epoch for the per-load
                                 RNG seed (0 = the constructor's bootstrap). */
  hs_wasm::ParamGenerationTracker
      param_generation; /**< Effect and descriptor-schema identity token. */
  bool animations_paused = false; /**< Pause state applied to every effect. */
};

/** @brief Registers the render bridge's enums and class with Embind. */
static void bind_engine() {
  enum_<ParamSetResult>("ParamSetResult")
      .value("APPLIED", ParamSetResult::APPLIED)
      .value("NO_EFFECT", ParamSetResult::NO_EFFECT)
      .value("UNKNOWN_PARAM", ParamSetResult::UNKNOWN_PARAM)
      .value("READONLY", ParamSetResult::READONLY)
      .value("NON_FINITE", ParamSetResult::NON_FINITE);

  enum_<ClipSetResult>("ClipSetResult")
      .value("APPLIED", ClipSetResult::APPLIED)
      .value("NO_EFFECT", ClipSetResult::NO_EFFECT)
      .value("INVALID_BOUNDS", ClipSetResult::INVALID_BOUNDS)
      .value("FULL_FRAME_KEPT", ClipSetResult::FULL_FRAME_KEPT);

  enum_<ResolutionSetResult>("ResolutionSetResult")
      .value("RESIZED", ResolutionSetResult::RESIZED)
      .value("ALREADY_ACTIVE", ResolutionSetResult::ALREADY_ACTIVE)
      .value("UNSUPPORTED", ResolutionSetResult::UNSUPPORTED);

  enum_<EffectSetResult>("EffectSetResult")
      .value("INSTALLED", EffectSetResult::INSTALLED)
      .value("UNKNOWN_EFFECT", EffectSetResult::UNKNOWN_EFFECT)
      .value("UNSUPPORTED_RESOLUTION", EffectSetResult::UNSUPPORTED_RESOLUTION);

  enum_<FullConfigRestoreResult>("FullConfigRestoreResult")
      .value("APPLIED", FullConfigRestoreResult::APPLIED)
      .value("NOT_SHADERBALL", FullConfigRestoreResult::NOT_SHADERBALL)
      .value("UNSUPPORTED_VERSION",
             FullConfigRestoreResult::UNSUPPORTED_VERSION)
      .value("INVALID_LENGTH", FullConfigRestoreResult::INVALID_LENGTH)
      .value("INVALID_VALUE", FullConfigRestoreResult::INVALID_VALUE)
      .value("INVALID_ACCEPTED", FullConfigRestoreResult::INVALID_ACCEPTED)
      .value("INVALID_PENDING", FullConfigRestoreResult::INVALID_PENDING);

  class_<HolosphereEngine>("HolosphereEngine")
      .constructor<>()
      .function("setResolution", &HolosphereEngine::setResolution)
      .function("setEffect", &HolosphereEngine::setEffect)
      .function("drawFrame", &HolosphereEngine::drawFrame)
      .function("getPixels", &HolosphereEngine::getPixels)
      .function("getBufferLength", &HolosphereEngine::getBufferLength)
      .function("setParameter", &HolosphereEngine::setParameter)
      .function("setAnimationsPaused", &HolosphereEngine::setAnimationsPaused)
      .function("getAnimationsPaused", &HolosphereEngine::getAnimationsPaused)
      .function("getPresetCount", &HolosphereEngine::getPresetCount)
      .function("getPresetIndex", &HolosphereEngine::getPresetIndex)
      .function("selectPreset", &HolosphereEngine::selectPreset)
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
      .function("getFullConfigSnapshot",
                &HolosphereEngine::getFullConfigSnapshot)
      .function("restoreFullConfigSnapshot",
                &HolosphereEngine::restoreFullConfigSnapshot)
      .function("getFullConfigFieldDefinitions",
                &HolosphereEngine::getFullConfigFieldDefinitions)
      .function("getConfigImportNotice",
                &HolosphereEngine::getConfigImportNotice)
      .function("clearConfigImportNotice",
                &HolosphereEngine::clearConfigImportNotice)
      .class_function("getSupportedResolutions",
                      &HolosphereEngine::getSupportedResolutions)
      .function("setClip", &HolosphereEngine::setClip)
      .function("strobeColumns", &HolosphereEngine::strobeColumns);
}
