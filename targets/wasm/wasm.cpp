/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/stack.h>
#include "engine/effects.h" // Includes all effect headers (triggers REGISTER_EFFECT)
#include "core/engine/effect_registry.h"
#include "engine/platform.h"
#include "targets/wasm/param_marshal.h"    // pure, host-tested param marshaling
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <algorithm> // std::fill_n — blank-frame clear in drawFrame
#include <string_view>
#include <cstring>
#include <cstdlib> // std::malloc for the lazily-allocated tooling arenas
#include <cmath>   // std::isfinite — validate MeshOps args at the JS boundary
#include <climits>  // INT_MAX — drawFrame pixel-index accumulator bound
#include <initializer_list> // all_finite() variadic-arg gate for free exports

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
 */
static size_t stack_high_water_mark() {
  uintptr_t base = emscripten_stack_get_base();
  uintptr_t end = emscripten_stack_get_end();
  const uint8_t *p = reinterpret_cast<const uint8_t *>(end);
  const uint8_t *top = reinterpret_cast<const uint8_t *>(base);
  while (p < top && *p == STACK_CANARY)
    p++;
  return static_cast<size_t>(top - p);
}

using namespace emscripten;

// Upper bound on a single effect's exposed parameters. Mirrors the fixed
// `std::array<ParamDef, 32>` backing `Effect::ParamList` (core/render/canvas.h). Used
// to pre-reserve the getParamValues() backing store so it never reallocates.
static constexpr size_t MAX_PARAMS = 32;

static_assert(MAX_PARAMS ==
                  std::tuple_size<decltype(Effect::ParamList::elements)>::value,
              "MAX_PARAMS must match Effect::ParamList's fixed array size");

// Arenas for the JS mesh-editor tools (8 MB build + two 4 MB scratch), used
// only by MeshOpsWrapper. malloc'd lazily on first MeshOps use (start at
// capacity 0) so engine/worker instances that never touch MeshOps don't reserve
// 16 MB; the block lives for the module's lifetime (reset via clearToolingMemory).
static constexpr size_t TOOLING_ARENA_BYTES = 8 * 1024 * 1024;
static constexpr size_t TOOLING_SCRATCH_BYTES = 4 * 1024 * 1024;
static Arena tooling_arena(nullptr, 0);
// Transient single-op scratch, shared module-globally. Every MeshOps entry
// point reset()s both at its head; valid only within one synchronous call. A
// ToolingOpGuard at each entry traps any re-entrant or interleaved use (e.g. a
// future worker or async refactor) before it can alias this scratch.
static Arena tooling_scratch_a(nullptr, 0);
static Arena tooling_scratch_b(nullptr, 0);

// Bumped on every clearToolingMemory(). Each wrapper records the generation it
// was built under and traps via check_live() if a wipe reclaimed its storage.
static uint32_t tooling_generation = 0;

// Set for the duration of one MeshOps entry point. The tooling scratch arenas
// are module-global and reset() at each entry's head, so a second op entered
// before the first returns would alias the first's scratch and corrupt its
// geometry. Synchronous single-threaded calls never overlap today; the guard
// makes that contract enforced rather than implicit, so a future worker/async
// refactor traps here instead of silently aliasing.
static bool tooling_op_active = false;
struct ToolingOpGuard {
  ToolingOpGuard() {
    HS_CHECK(!tooling_op_active,
             "re-entrant MeshOps call aliases module-global tooling scratch");
    tooling_op_active = true;
  }
  ~ToolingOpGuard() { tooling_op_active = false; }
};

/**
 * @brief Allocates and binds the tooling arenas on first MeshOps use.
 * @details A no-op once bound, so it is cheap to call at the head of every
 *          MeshOps entry point. Reading an unbound arena's metrics
 *          (collect_arena_metrics) is safe and reports 0/0/0, so engine
 *          instances that never call MeshOps never trigger this allocation.
 */
static void ensure_tooling_arenas() {
  if (tooling_arena.get_capacity() != 0)
    return;
  const size_t total = TOOLING_ARENA_BYTES + 2 * TOOLING_SCRATCH_BYTES;
  uint8_t *block = static_cast<uint8_t *>(std::malloc(total));
  HS_CHECK(block != nullptr); // 16 MB tooling block — fail-fast on OOM
  tooling_arena.rebind(block, TOOLING_ARENA_BYTES);
  tooling_scratch_a.rebind(block + TOOLING_ARENA_BYTES, TOOLING_SCRATCH_BYTES);
  tooling_scratch_b.rebind(block + TOOLING_ARENA_BYTES + TOOLING_SCRATCH_BYTES,
                           TOOLING_SCRATCH_BYTES);
}

/**
 * @brief Builds a {usage, high_water_mark, capacity} report for the four engine
 *        arenas.
 * @return JS object mapping each arena name to its {usage, high_water_mark,
 *         capacity} metrics, in bytes.
 * @details Shared by HolosphereEngine and MeshOpsWrapper. Callers that also want
 *          stack metrics append them to the returned object.
 */
static val collect_arena_metrics() {
  val metrics = val::object();
  auto add_metrics = [&](const char *name, Arena &arena) {
    val m = val::object();
    m.set("usage", arena.get_offset());
    m.set("high_water_mark", arena.get_high_water_mark());
    m.set("capacity", arena.get_capacity());
    metrics.set(name, m);
  };
  add_metrics("scratch_arena_a", scratch_arena_a);
  add_metrics("scratch_arena_b", scratch_arena_b);
  add_metrics("persistent_arena", persistent_arena);
  add_metrics("tooling_arena", tooling_arena);
  return metrics;
}

/**
 * @brief Builds a concrete factory table from the self-registering entries.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @return Reference to the lazily-built, static per-(W,H) factory table.
 */
template <int W, int H> const std::vector<FactoryEntry> &get_factory() {
  static std::vector<FactoryEntry> table = []() {
    const auto &regs = EffectRegistry::entries();
    std::vector<FactoryEntry> t(regs.size());
    for (size_t i = 0; i < regs.size(); ++i)
      get_fill_fn<W, H>(regs[i])(t[i]);
    // Duplicate names silently shadow (lookups return the first match); the
    // names aren't known until the fill functions run, so trap here. The O(n²)
    // pairwise scan is deliberate: the table is built once per (W,H) and the
    // roster is tiny (~27), so a sort-and-adjacent-compare or a set buys no
    // measurable time and adds code.
    for (size_t i = 0; i < t.size(); ++i)
      for (size_t j = i + 1; j < t.size(); ++j)
        HS_CHECK(t[i].name != t[j].name,
                 "get_factory: duplicate effect name in registry");
    return t;
  }();
  return table;
}

/**
 * @brief Reports whether the (W,H) factory contains an effect with the given
 *        name.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param name Effect name to look up.
 * @return true if an effect with this name is registered for (W,H).
 * @details Cheap linear scan used by setEffect() to validate a stale/typo'd UI
 *          string BEFORE it tears down the running effect, so an unknown name is
 *          a transactional no-op rather than a blanked engine.
 */
template <int W, int H> bool factory_has_effect(std::string_view name) {
  for (const auto &entry : get_factory<W, H>())
    if (name == entry.name)
      return true;
  return false;
}

/**
 * @brief Instantiates the named effect from the (W,H) factory.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param name Effect name to instantiate.
 * @return Owning pointer to the new effect, or null if the name is unknown (a
 *         typo'd/stale UI string).
 */
template <int W, int H>
std::unique_ptr<Effect> create_effect(std::string_view name) {
  const auto &factory = get_factory<W, H>();
  for (const auto &entry : factory) {
    if (name == entry.name)
      return entry.creator();
  }
  // Unknown name: return null (setEffect's guard makes it a safe no-op) rather
  // than silently substituting a different effect.
  hs::log("WASM: create_effect: unknown effect name (no effect created)");
  return nullptr;
}

// ---------------------------------------------------------------------------
// The (W,H) resolutions the WASM factory can build. Aliased to the registry's
// HS_RESOLUTIONS (core/engine/effect_registry.h) so the runtime dispatch here and the
// per-resolution fill functions the registry generates share one list and
// cannot drift: a resolution the registry can build is dispatchable here with no
// second edit. setResolution()/setEffect()/getEffectSizes() all expand this via
// the X-macro below. A new resolution is one edit, in HS_RESOLUTIONS (the effect
// templates must also be instantiable at that <W,H>).
// ---------------------------------------------------------------------------
#define HS_WASM_RESOLUTIONS(X) HS_RESOLUTIONS(X)

// Pin every resolution row to the MAX_W×MAX_H pixel-buffer bound.
#define X(W, H)                                                                 \
  static_assert((W) <= MAX_W && (H) <= MAX_H,                                   \
                "HS_WASM_RESOLUTIONS row exceeds the MAX_W×MAX_H pixel buffer");
HS_WASM_RESOLUTIONS(X)
#undef X

/**
 * @brief Invokes f.operator()<W,H>() for the single HS_WASM_RESOLUTIONS row
 *        matching the runtime (w,h).
 * @param f A C++20 templated callable (e.g. `[]<int W, int H>(){...}`) run with
 *          the matching row's dimensions as compile-time template arguments.
 * @return true iff a row matched and f was invoked; false otherwise.
 * @details One shared expansion of the resolution list for every runtime
 *          dispatch site (setEffect validate/create, resolution checks), so they
 *          can never drift and a new per-resolution step lands in exactly one
 *          place. Fully inlined at -O2 — identical codegen to an open-coded
 *          X-macro, no runtime cost.
 */
template <typename F> static bool dispatch_resolution(int w, int h, F &&f) {
#define X(W, H)                                                                 \
  if (w == (W) && h == (H)) {                                                   \
    f.template operator()<(W), (H)>();                                          \
    return true;                                                                \
  }
  HS_WASM_RESOLUTIONS(X)
#undef X
  return false;
}

/**
 * @brief Reports whether (w,h) is a resolution the factory can build.
 * @param w Candidate canvas width in pixels.
 * @param h Candidate canvas height in pixels.
 * @return true iff (w,h) is one of the HS_WASM_RESOLUTIONS rows.
 */
static bool wasm_resolution_supported(int w, int h) {
  return dispatch_resolution(w, h, []<int W, int H>() {});
}

/**
 * @brief JS-facing render engine driving one resolution/effect at a time.
 * @details Owns the current effect and the stable readback buffers. Every public
 *          method is exported to JavaScript via EMSCRIPTEN_BINDINGS below.
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
    stack_paint_canary();

    // SSOT guard: the self-registering effect count must match the static roster
    // (HS_EFFECT_LIST / HS_EFFECT_COUNT) or the live set and the native smoke
    // suite silently diverge.
    HS_CHECK(EffectRegistry::entries().size() ==
             static_cast<size_t>(HS_EFFECT_COUNT));

    // Pre-size the view-backed readback buffers ONCE: under ALLOW_MEMORY_GROWTH
    // a reallocation detaches the ArrayBuffer behind a typed_memory_view, so the
    // buffers returned as views (getPixels/getParamValues) must never move.
    pixelBuffer.assign(MAX_W * MAX_H * CHANNELS, 0);
    paramValues.reserve(MAX_PARAMS);
    paramViews.reserve(MAX_PARAMS); // not view-backed; reserve is amortization only

    // Bootstrap default; daydream overrides it almost immediately.
    setResolution(96, 20);
    HS_CHECK(setEffect("DistortedRing"));
  }

  /**
   * @brief Switches the active canvas resolution.
   * @param w Requested canvas width in pixels.
   * @param h Requested canvas height in pixels.
   * @return true if the resolution is now active; false if the request was
   *         rejected (unsupported size) and the previous valid state was kept.
   * @details A successful resolution change tears down the current effect (a new
   *          one cannot be carried across pixel dimensions), so the caller must
   *          call setEffect() again before the next drawFrame() or it renders a
   *          blank frame. Callers that use a sub-canvas clip (segmented rendering)
   *          must likewise re-apply setClip() after a successful setResolution():
   *          a prior clip was expressed in the old resolution's pixel bounds and
   *          is not rescaled here, so it must be recomputed for the new
   *          dimensions.
   */
  bool setResolution(int w, int h) {
    if (w == pixel_width && h == pixel_height)
      return true;

    // Reject unsupported sizes and keep the prior valid state alive rather than
    // switching to a null effect that renders blank with no signal to JS.
    if (!wasm_resolution_supported(w, h)) {
      hs::log("WASM: Unsupported resolution %dx%d — ignored", w, h);
      return false;
    }

    pixel_width = w;
    pixel_height = h;

    if (currentEffect) {
      currentEffect = nullptr;
      stack_paint_canary(); // repaint to reset stack HWM after teardown
    }
    return true;
  }

  /**
   * @brief Tears down the current effect and instantiates the named one at the
   *        active resolution.
   * @param name Effect name to instantiate.
   * @return true iff an effect was actually instantiated; false for an
   *         unknown/stale effect name or an unsupported resolution, so the
   *         frontend can detect a no-op instead of believing the switch took.
   * @details Validates the name against the factory for the current resolution
   *          BEFORE tearing anything down, so a typo'd/stale UI string keeps the
   *          prior valid state alive rather than blanking the engine.
   */
  bool setEffect(std::string name) {
    // Pass name as the %s argument, never as the format string itself: a '%' in
    // an effect name would otherwise consume uninitialized varargs.
    hs::log("WASM: setEffect called with %s", name.c_str());

    // Validate against the current resolution's factory BEFORE tearing anything
    // down, so a typo'd name keeps the prior valid state alive.
    bool name_valid = false;
    dispatch_resolution(pixel_width, pixel_height, [&]<int W, int H>() {
      name_valid = factory_has_effect<W, H>(name);
    });
    if (!name_valid) {
      hs::log("WASM: setEffect unknown effect '%s' — keeping current effect",
              name.c_str());
      return false;
    }

    currentEffect.reset();
    // Configured before construction; effect constructors must not allocate from
    // the engine arenas (Teensy configures after construction — see Phantasm.ino).
    configure_arenas_default();

    stack_paint_canary(); // reset stack HWM by repainting unused region

    bool created = dispatch_resolution(
        pixel_width, pixel_height, [&]<int W, int H>() {
          init_geometry_luts<W, H>(); // eager-fill LUTs before the first frame
          currentEffect = create_effect<W, H>(name);
        });
    if (!created) {
      hs::log("WASM: Unsupported resolution for factory!"); // unreachable guard
      return false;
    }
    if (!currentEffect) {
      return false; // unreachable: name was validated above
    }
    currentEffect->init();
    hs::log("WASM: init stack HWM = %u bytes", (unsigned)stack_high_water_mark());
    stack_paint_canary();
    return true;
  }

  /**
   * @brief Restricts rendering to a clip band for the current effect.
   * @param x0 Inclusive left column of the clip band in [0, pixel_width].
   * @param x1 Exclusive right column, with x0 <= x1 <= pixel_width.
   * @param y0 Inclusive top row of the clip band in [0, pixel_height].
   * @param y1 Exclusive bottom row of the clip band, with y0 <= y1 <= pixel_height.
   * @return true if the bounds were accepted (band applied, OR intentionally
   *         ignored for a full-frame stateful effect); false if no effect is set
   *         or the bounds are malformed/out of range (then ignored).
   * @details Args are x-pair-first to match the (x, y) convention: embind binds
   *          positionally, so a y-first order would let a transposed JS call pass
   *          the range check and silently clip the wrong axis. Rejects malformed
   *          input at the untyped JS boundary rather than trapping, since a trap
   *          there aborts the whole WASM module. Segment workers always pass
   *          valid, ordered, in-range bounds. A cross-segment stateful effect
   *          (Effect::needs_full_frame()) keeps the full-canvas clip instead of
   *          narrowing to the band — see docs/segmented_stateful_effects_spec.md.
   */
  bool setClip(int x0, int x1, int y0, int y1) {
    if (!currentEffect)
      return false;
    // Reject malformed bounds from the untyped JS boundary (negatives would feed
    // ClipRegion's modulo arithmetic); reject-and-return, never trap.
    if (!hs_wasm::clip_bounds_valid(x0, x1, y0, y1, pixel_width,
                                    pixel_height)) {
      hs::log("WASM: setClip bounds out of range (x0=%d,x1=%d,y0=%d,y1=%d) — "
              "ignored",
              x0, x1, y0, y1);
      return false;
    }
    // Cross-segment stateful effects must render the FULL canvas in every worker
    // (a band-clipped worker has stale cv.prev outside its band, so trails seam);
    // keep the full clip. See docs/segmented_stateful_effects_spec.md.
    if (currentEffect->needs_full_frame())
      return true;
    currentEffect->set_clip(y0, y1, x0, x1);
    return true;
  }

  /**
   * @brief Renders one frame of the current effect into the JS-facing buffer.
   * @details Copies the effect's full canvas into pixelBuffer as 16-bit linear
   *          RGB triples; no-op if no effect is set. The readback copies the
   *          FULL canvas regardless of any active clip region — a clip restricts
   *          rendering, not this readback.
   */
  void drawFrame() {
    if (!currentEffect) {
      // No active effect: clear the active prefix so getPixels() hands JS a
      // blank frame at the current resolution, not stale content.
      const int count = pixel_width * pixel_height * CHANNELS;
      std::fill_n(pixelBuffer.data(), count, uint16_t{0});
      return;
    }

    currentEffect->render_us = 0.0;
    currentEffect->draw_frame();
    currentEffect->advance_display();

    // Readback copies the FULL canvas regardless of any clip; segment_worker.js
    // extracts its quadrant JS-side (README §10.7).
    static_assert(static_cast<long long>(MAX_W) * MAX_H * CHANNELS <= INT_MAX,
                  "drawFrame pixel-index accumulators are int");
    int idx = 0;
    const int count = pixel_width * pixel_height;
    if (!currentEffect->overrides_get_pixel()) {
      // Fast path: display_buffer()[i] == get_pixel(x, y), so copy directly.
      const Pixel *buf = currentEffect->display_buffer();
      static_assert(sizeof(Pixel) == 3 * sizeof(uint16_t),
                    "fast-path memcpy assumes packed RGB16 Pixel layout");
      std::memcpy(pixelBuffer.data(), buf, static_cast<size_t>(count) * sizeof(Pixel));
    } else {
      for (int y = 0; y < pixel_height; y++) {
        for (int x = 0; x < pixel_width; x++) {
          const Pixel &p = currentEffect->get_pixel(x, y);
          pixelBuffer[idx++] = p.r;
          pixelBuffer[idx++] = p.g;
          pixelBuffer[idx++] = p.b;
        }
      }
    }
  }

  /**
   * @brief Returns the last frame's render time.
   * @return Render duration of the most recent drawFrame() in microseconds, or
   *         0.0 if no effect is set.
   * @note render_us is accumulated by the scan/plot ScopedRenderTimer, so
   *       full-screen shader effects that bypass that path (e.g. Raymarch,
   *       Liquid2D, Flyby — they shade every pixel directly rather than via
   *       Scan/Plot) report 0.0 even while actively rendering. A 0 here means
   *       "not instrumented," not "free"; it is not a reliable cross-effect cost
   *       metric for those effects.
   */
  double getRenderUs() {
    return currentEffect ? currentEffect->render_us : 0.0;
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
    return currentEffect ? currentEffect->strobe_columns() : false;
  }

  /**
   * @brief Exposes the raw pixel buffer to JS as a zero-copy Uint16Array view.
   * @return Typed memory view over the active resolution's R,G,B pixels within
   *         the stable MAX_W*MAX_H*3 backing buffer.
   * @details WASM memory-view contract: the returned view aliases WASM linear
   *          memory, it is NOT a copy. With ALLOW_MEMORY_GROWTH=1, any subsequent
   *          heap growth detaches the underlying ArrayBuffer and leaves this view
   *          zero-length (buffer.byteLength === 0). Callers MUST re-fetch the
   *          view after anything that may allocate (resolution/effect change) and
   *          may only cache it across frames while guarding for detachment. The
   *          backing vector is pre-sized once and never reallocated, so the only
   *          possible detachment source is heap growth elsewhere.
   *          daydream.js::refreshPixelView mirrors this expectation.
   */
  val getPixels() {
    return val(typed_memory_view(pixel_width * pixel_height * CHANNELS,
                                 pixelBuffer.data()));
  }

  /**
   * @brief Returns the length of the active pixel buffer view.
   * @return Number of uint16 elements in the active view (pixel_width *
   *         pixel_height * 3, three channels per pixel).
   */
  int getBufferLength() { return pixel_width * pixel_height * CHANNELS; }

  /**
   * @brief Updates one named effect parameter.
   * @param name Parameter name to update.
   * @param value New parameter value, in the parameter's native units.
   * @return true if the write was accepted; false otherwise.
   * @details Like setClip()/setResolution(), the boolean is intentionally
   *          coarse — it does NOT let JS distinguish the rejection reasons. A
   *          false collapses four distinct cases: no effect is set, the name is
   *          unknown to the effect, the parameter is readonly (engine-written
   *          telemetry the GUI must not poke), or the value is non-finite
   *          (rejected before it can poison render math; see
   *          Canvas::updateParameter). A true means the value was applied, but
   *          note an accepted float is silently clamped to the parameter's
   *          registered [min,max] — true does NOT imply the stored value equals
   *          the requested one. A consumer that needs the effective value should
   *          read it back via getParamValues() rather than trust this flag.
   */
  bool setParameter(std::string name, float value) {
    if (!currentEffect)
      return false;
    // Finiteness is single-sourced in Canvas::updateParameter, not re-checked
    // here.
    return currentEffect->updateParameter(name.c_str(), value);
  }

  /**
   * @brief Pauses or resumes the current effect's parameter animations.
   * @param paused true to pause animations, false to resume. No-op if no effect
   *        is set.
   */
  void setAnimationsPaused(bool paused) {
    if (currentEffect) {
      currentEffect->setAnimationsPaused(paused);
    }
  }

  /**
   * @brief Builds the GUI's parameter descriptor list.
   * @return JS array with one {name, value, animated, readonly, (+min/max for
   *         floats)} object per param in the effect's declaration order; empty
   *         array when no effect is set.
   */
  val getParameterDefinitions() {
    if (!currentEffect)
      return val::array();

    val result = val::array();
    // Single source of order shared with getParamValues() (param_marshal.h), so
    // the value stream cannot index-drift from these definitions.
    hs_wasm::collect_param_views(*currentEffect, paramViews);

    int i = 0;
    for (const auto &v : paramViews) {
      val entry = val::object();
      entry.set("name", val(v.name));

      if (v.is_bool) {
        // Emit a JS boolean so the frontend renders a checkbox; toggles omit
        // min/max (no range).
        entry.set("value", val(v.value > 0.5f));
      } else {
        entry.set("value", v.value);
        entry.set("min", v.min);
        entry.set("max", v.max);
      }
      entry.set("animated", val(v.animated));
      entry.set("readonly", val(v.readonly));
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
   *          memory and must be consumed before the next allocation. paramValues
   *          never reallocates here (size <= MAX_PARAMS), so emitting it triggers
   *          no heap growth that could detach other outstanding views.
   */
  val getParamValues() {
    if (!currentEffect) {
      // Empty Float32Array (not a JS Array) so callers get a consistent typed
      // view whether or not an effect is set. The zero-length view still needs a
      // valid backing pointer, which holds only because the ctor reserved
      // MAX_PARAMS and clear() retains that capacity.
      HS_CHECK(paramValues.capacity() >= MAX_PARAMS);
      paramValues.clear();
      return val(typed_memory_view(paramValues.size(), paramValues.data()));
    }

    // Same order as getParameterDefinitions(); clear retains MAX_PARAMS capacity.
    hs_wasm::fill_param_values(*currentEffect, paramValues);
    return val(typed_memory_view(paramValues.size(), paramValues.data()));
  }

  /**
   * @brief Reports engine arena and stack metrics for the JS memory HUD.
   * @return JS object of arena metrics plus a "stack" entry, each in the
   *         {usage, high_water_mark, capacity} shape, all in bytes.
   */
  val getArenaMetrics() {
    val metrics = collect_arena_metrics();

    // Stack metrics (same format as arenas)
    {
      uintptr_t base = emscripten_stack_get_base();
      uintptr_t end = emscripten_stack_get_end();
      uintptr_t sp = emscripten_stack_get_current();
      val m = val::object();
      m.set("usage", static_cast<size_t>(base >= sp ? base - sp : 0));
      m.set("high_water_mark", static_cast<size_t>(stack_high_water_mark()));
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
   * @brief Enumerates the resolutions the factory can build.
   * @return JS array of [W, H] pairs, generated from the same
   *         HS_WASM_RESOLUTIONS list that setResolution()/getEffectSizes()
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
    HS_WASM_RESOLUTIONS(X)
#undef X
    return out;
  }

private:
  /** Channels per pixel in the readback buffer (linear RGB triples). */
  static constexpr int CHANNELS = 3;

  std::unique_ptr<Effect> currentEffect; /**< Currently active effect, or null. */
  std::vector<uint16_t> pixelBuffer; /**< 16-bit linear RGB readback buffer. */
  std::vector<float> paramValues;    /**< Backing store for getParamValues. */
  std::vector<hs_wasm::ParamView> paramViews; /**< Scratch for getParameterDefinitions. */
  int pixel_width = 0;  /**< Active canvas width in pixels. */
  int pixel_height = 0; /**< Active canvas height in pixels. */
};

// ==========================================================================================
// MESH OPS BINDINGS
// ==========================================================================================

#include "mesh/solids.h"

/**
 * @brief JS-facing wrapper around a PolyMesh and the Conway/Goldberg operators.
 * @details Named to avoid collision with the MeshOps namespace. Each wrapper's
 *          mesh is built into the tooling arena and records the generation it was
 *          built under, so a wipe via clearToolingMemory() is detected by
 *          check_live().
 */
struct MeshOpsWrapper {
private:
  PolyMesh mesh; /**< The wrapped mesh, stored in the tooling arena. */
  /**
   * Generation of the tooling arena this mesh was built into; compared against
   * the live counter on every use (see check_live()).
   */
  uint32_t generation_ = tooling_generation;

public:
  /**
   * @brief Constructs a wrapper taking ownership of an existing mesh.
   * @param m Mesh to move into this wrapper.
   */
  MeshOpsWrapper(PolyMesh &&m) : mesh(std::move(m)) {}

  /**
   * @brief Traps if this wrapper outlived a clearToolingMemory() wipe.
   * @details Its mesh would alias reclaimed arena storage, which release builds
   *          would otherwise read back as silently wrong geometry. Called at
   *          every entry point that touches `mesh`. This is JS-editor tooling,
   *          never the render loop, so the compare is free.
   */
  void check_live() const {
    HS_CHECK(generation_ == tooling_generation,
             "MeshOps wrapper used after clearToolingMemory()");
  }

  /**
   * @brief Resets all tooling arenas to empty and invalidates live wrappers.
   * @details Reclaims the storage behind every live wrapper and bumps the
   *          generation so any wrapper built before this wipe traps on next use
   *          (check_live). JS-callable.
   *
   *          Despite the name, this does NOT shrink the module's linear memory:
   *          the 16 MB tooling block is retained for the module's lifetime and
   *          only its arena bump-pointers are reset. A JS caller will not see
   *          memory usage drop; this frees the arenas for reuse, not the OS.
   */
  static void clearToolingMemory() {
    tooling_arena.reset();
    tooling_arena.reset_high_water_mark();
    tooling_scratch_a.reset();
    tooling_scratch_a.reset_high_water_mark();
    tooling_scratch_b.reset();
    tooling_scratch_b.reset_high_water_mark();
    ++tooling_generation;
  }

  /**
   * @brief Builds a wrapper for the named base solid.
   * @param name Solid name to look up in the Solids registry.
   * @return Owning pointer to the new wrapper, or null for an unknown name.
   * @details Rejects an unknown name at the untrusted JS boundary rather than
   *          tripping get_by_name()'s fail-fast HS_CHECK and aborting the module.
   */
  static std::unique_ptr<MeshOpsWrapper> fromSolidName(std::string name) {
    if (!Solids::has_name(name)) {
      hs::log("WASM: fromSolidName unknown solid '%s' — ignored", name.c_str());
      return nullptr;
    }
    ToolingOpGuard guard;
    ensure_tooling_arenas();
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::get_by_name(
        tooling_arena, tooling_scratch_a, tooling_scratch_b, name));
  }

  /**
   * @brief Returns the mesh vertices as a JS Float32Array.
   * @return Float32Array of flattened [x,y,z] triples, copied out of the mesh.
   */
  val getVertices() const {
    check_live();
    std::vector<float> data;
    data.reserve(mesh.vertices.size() * 3);
    for (const auto &v : mesh.vertices) {
      data.push_back(v.x);
      data.push_back(v.y);
      data.push_back(v.z);
    }
    return val::global("Float32Array")
        .new_(val(typed_memory_view(data.size(), data.data())));
  }

  /**
   * @brief Returns the mesh faces as flat index + per-face side-count buffers.
   * @return JS object `{ indices: Uint16Array, counts: Uint8Array }`; JS
   *         unflattens the parallel arrays into per-face index lists. Both are
   *         copied out of WASM memory (same tooling-arena lifetime contract as
   *         getVertices()), so they are safe to hold across later mesh ops.
   */
  val getFaces() const {
    check_live();
    size_t total = 0;
    for (size_t i = 0; i < mesh.get_face_counts_size(); ++i)
      total += mesh.get_face_counts_data()[i];
    HS_CHECK(total == mesh.get_faces_size(),
             "getFaces: face_counts sum disagrees with the flat index count");
    val out = val::object();
    out.set("indices", val::global("Uint16Array").new_(val(typed_memory_view(
                           mesh.get_faces_size(), mesh.get_faces_data()))));
    out.set("counts",
            val::global("Uint8Array").new_(val(typed_memory_view(
                mesh.get_face_counts_size(), mesh.get_face_counts_data()))));
    return out;
  }

  /**
   * @brief Classifies faces by topology and returns the per-face codes.
   * @return JS Int32Array of one topology code per face, copied out of the
   *         mesh's now-populated topology buffer.
   * @details Same tooling-arena lifetime contract as getVertices(): the
   *          `topology` buffer lives in tooling_arena and is invalidated by the
   *          next mesh op / arena reset. The `.new_(Int32Array)(view)` form
   *          *copies* the typed_memory_view into a fresh JS array (it does not
   *          alias WASM memory), so the result is safe to hold across later
   *          calls — but if this is ever changed to return the view directly
   *          (as getPixels/bakeLut do), it MUST then be read before the next
   *          allocation, per that memory-view contract.
   */
  val classifyFaces() {
    check_live();
    ToolingOpGuard guard;
    ensure_tooling_arenas();
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    MeshOps::classify_faces_by_topology(mesh, tooling_scratch_a,
                                        tooling_scratch_b, tooling_arena);
    // Copies (see the contract note above); does not alias WASM memory.
    return val::global("Int32Array")
        .new_(
            val(typed_memory_view(mesh.topology.size(), mesh.topology.data())));
  }

  // --- Conway/Goldberg operators -------------------------------------------

  /**
   * @brief Runs a mesh operator and wraps the result.
   * @tparam Op Callable of signature (const PolyMesh&, Arena&, Arena&) ->
   *         PolyMesh.
   * @param op Operator to run against this wrapper's mesh.
   * @return Owning pointer to a new wrapper holding the finalized result mesh.
   * @details Captures the shared operator boilerplate: reset both tooling scratch
   *          arenas, run the op into a fresh PolyMesh, finalize it into
   *          tooling_arena, and hand back a new wrapper. This is tooling, never
   *          the render loop, and every lambda inlines, so there is no added cost.
   */
  template <typename Op>
  std::unique_ptr<MeshOpsWrapper> apply(Op &&op) const {
    check_live();
    ToolingOpGuard guard;
    ensure_tooling_arenas();
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        op(mesh, tooling_scratch_a, tooling_scratch_b), tooling_arena));
  }

  /**
   * @brief Validates that an operator argument is finite.
   * @param arg Operator argument crossing the untrusted JS boundary.
   * @param op Operator name, for the rejection log message.
   * @return true if arg is finite; false (after logging) otherwise.
   * @details A non-finite fraction would flow straight into the geometry math and
   *          silently corrupt the mesh, so it is rejected (log + null) rather than
   *          producing NaN geometry.
   */
  bool finite_arg(float arg, const char *op) const {
    if (std::isfinite(arg))
      return true;
    hs::log("WASM: MeshOps::%s got a non-finite argument — ignored", op);
    return false;
  }

/**
 * @brief Defines a zero-argument Conway/Goldberg operator method.
 * @param name MeshOps operator name; becomes the generated method name.
 * @details The generated method runs MeshOps::name(mesh) via apply() and returns
 *          a new wrapper holding the result.
 */
#define MESHOP_0(name)                                                         \
  std::unique_ptr<MeshOpsWrapper> name() const {                              \
    return apply(                                                              \
        [](const PolyMesh &m, Arena &a, Arena &b) { return MeshOps::name(m, a, b); });           \
  }
/**
 * @brief Defines a one-float-argument Conway/Goldberg operator method.
 * @param name MeshOps operator name; becomes the generated method name.
 * @details The generated method rejects a non-finite arg (finite_arg) before
 *          running MeshOps::name(mesh, arg) via apply(), returning a new wrapper
 *          or null.
 */
#define MESHOP_1F(name)                                                        \
  std::unique_ptr<MeshOpsWrapper> name(float arg) const {                     \
    if (!finite_arg(arg, #name))                                              \
      return nullptr;                                                          \
    return apply([arg](const PolyMesh &m, Arena &a, Arena &b) {               \
      return MeshOps::name(m, a, b, arg);                                     \
    });                                                                        \
  }

/**
 * @brief Defines a one-float-argument operator whose argument is a [0,1]
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @details Like MESHOP_1F but clamps the fraction to [0,1] at the JS boundary
 *          (logging when it changed) so a direct/API caller passing a finite
 *          out-of-range value stays within the operator's documented domain —
 *          and, for operators whose fraction reaches an always-on HS_CHECK
 *          (truncate/bevel), cannot trip that trap and abort the whole module.
 */
#define MESHOP_1U(name)                                                        \
  std::unique_ptr<MeshOpsWrapper> name(float arg) const {                     \
    if (!finite_arg(arg, #name))                                              \
      return nullptr;                                                          \
    if (arg < 0.0f || arg > 1.0f)                                            \
      hs::log("WASM: MeshOps::%s clamped %g to [0,1]", #name, arg);           \
    float t = arg < 0.0f ? 0.0f : (arg > 1.0f ? 1.0f : arg);                  \
    return apply([t](const PolyMesh &m, Arena &a, Arena &b) {                 \
      return MeshOps::name(m, a, b, t);                                       \
    });                                                                        \
  }

/**
 * @brief Single source of truth for the Conway/Goldberg operator roster.
 * @param _OP0  Macro applied to each zero-argument operator name.
 * @param _OP1F Macro applied to each one-float-argument operator name.
 * @param _OP1U Macro applied to each [0,1]-fraction operator name.
 * @details Expanded twice: with MESHOP_0/MESHOP_1F/MESHOP_1U to generate the
 *          wrapper methods (below), and with MESHOP_BIND to generate the embind
 *          .function() bindings (in EMSCRIPTEN_BINDINGS), so an operator cannot
 *          be added to one site and silently left unreachable from the other.
 *          truncate, bevel, and chamfer use _OP1U: each has a documented [0,1]
 *          domain (truncate/bevel additionally reach an always-on engine trap).
 *          expand takes an unbounded factor, so it stays _OP1F. relax, hankin,
 *          and snub have bespoke signatures/validation (snub takes two float
 *          controls), so their wrapper methods are hand-written; their names
 *          live in MESHOP_IRREGULAR_LIST below and their bindings expand from
 *          it, so a new irregular op is bound the moment it joins the list.
 */
#define MESHOP_LIST(_OP0, _OP1F, _OP1U)                                         \
  _OP0(kis) _OP0(ambo) _OP0(gyro) _OP0(dual) _OP0(meta)                         \
  _OP0(needle) _OP0(zip)                                                        \
  _OP1F(expand)                                                                 \
  _OP1U(truncate) _OP1U(bevel) _OP1U(chamfer)

// Irregular ops: hand-written wrapper methods (custom signatures/validation),
// enumerated here so their embind bindings expand from one list.
#define MESHOP_IRREGULAR_LIST(_) _(relax) _(hankin) _(snub)

  MESHOP_LIST(MESHOP_0, MESHOP_1F, MESHOP_1U)

#undef MESHOP_0
#undef MESHOP_1F
#undef MESHOP_1U

  /**
   * Engine-enforced ceiling on relax smoothing passes. The editor
   * (solids.html) caps its slider at 500; this 1000 is deliberate headroom
   * above that cap, an independent defense-in-depth limit for direct/API
   * callers that bypass the editor. The two numbers are intentionally
   * different, not a mismatch.
   */
  static constexpr int MAX_RELAX_ITERATIONS = 1000;

  /**
   * @brief Applies relax smoothing passes to the mesh.
   * @param iterations Number of smoothing passes; floored at 0 and clamped to
   *        MAX_RELAX_ITERATIONS.
   * @return Owning pointer to a new wrapper holding the relaxed mesh.
   * @details Explicit (not a MESHOP_* macro) because its int iteration count
   *          crosses the JS boundary unbounded: relax(1e9) would freeze the main
   *          thread for billions of passes, so the count is clamped rather than
   *          trusted.
   */
  std::unique_ptr<MeshOpsWrapper> relax(int iterations) const {
    int clamped = hs_wasm::clamp_relax_iterations(iterations, MAX_RELAX_ITERATIONS);
    if (clamped != iterations)
      hs::log("WASM: MeshOps::relax clamped %d iterations to %d", iterations,
              clamped);
    iterations = clamped;
    return apply([iterations](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::relax(m, a, b, iterations);
    });
  }

  /**
   * @brief Applies the Hankin interlace operator to the mesh.
   * @param radians Interlace angle in radians (the unit MeshOps::hankin expects).
   * @return Owning pointer to a new wrapper holding the result, or null if the
   *         angle is non-finite.
   * @details Explicit (not a MESHOP_* macro) so the radians unit contract the JS
   *          caller relies on is carried here.
   */
  std::unique_ptr<MeshOpsWrapper> hankin(float radians) const {
    if (!finite_arg(radians, "hankin"))
      return nullptr;
    return apply([radians](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::hankin(m, a, b, radians);
    });
  }

  /**
   * @brief Applies the chiral snub operator with explicit inset and twist.
   * @param t Inset factor of each face toward its centroid, clamped to [0, 1]
   *          (its documented domain) at the JS boundary.
   * @param twist Per-face rotation about the face normal, in radians (0 = none);
   *          unbounded, so only finiteness is checked.
   * @return Owning pointer to a new wrapper, or null if either arg is non-finite.
   * @details Explicit (not a MESHOP_* macro) because MeshOps::snub takes TWO
   *          float controls, which neither the zero-arg nor the one-float
   *          generator can express. Binding it via MESHOP_0 (as it was) hardcodes
   *          the (0.5, 0.0) defaults and leaves both controls unreachable from JS;
   *          this 2-arg form exposes them to the solids editor.
   */
  std::unique_ptr<MeshOpsWrapper> snub(float t, float twist) const {
    if (!finite_arg(t, "snub") || !finite_arg(twist, "snub"))
      return nullptr;
    if (t < 0.0f || t > 1.0f)
      hs::log("WASM: MeshOps::snub clamped t=%g to [0,1]", t);
    float ct = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
    return apply([ct, twist](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::snub(m, a, b, ct, twist);
    });
  }
  /**
   * @brief Lists all available solids for the editor's solid picker.
   * @return JS array of {name, category} objects, one per registered solid.
   */
  static val getRegistry() {
    val registry = val::array();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      const auto &entry = Solids::get_entry(i);
      val item = val::object();
      item.set("name", val(entry.name));
      item.set("category", entry.category == Solids::Category::Simple
                               ? "Simple"
                               : "Complex");
      registry.set(i, item);
    }
    return registry;
  }

#ifdef HS_WASM_DEV_BINDINGS
  /**
   * @brief Measures the maximum vertex/face/index counts across all solids.
   * @return JS object with {max_v, v_name, max_f, f_name, max_i, i_name} giving
   *         the largest counts and the solids that produce them.
   * @details Dev-only roster measurement for sizing MAX_VERTS-style constants;
   *          no UI consumer. Off by default; enable the HS_WASM_DEV_BINDINGS
   *          CMake option to compile + re-export it
   *          (`cmake --preset wasm-release -DHS_WASM_DEV_BINDINGS=ON`; see
   *          CMakeLists.txt). Measures each solid in the scratch arenas only —
   *          never tooling_arena, which backs live wrappers the JS side holds.
   */
  static val getMaxBounds() {
    int max_v = 0;
    int max_f = 0;
    int max_i = 0;
    const char *mv_name = "";
    const char *mf_name = "";
    const char *mi_name = "";

    ToolingOpGuard guard;
    ensure_tooling_arenas();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      // Measure in the scratch arenas only — never tooling_arena, which backs
      // live wrappers the JS side holds.
      tooling_scratch_a.reset();
      tooling_scratch_b.reset();
      PolyMesh temp =
          Solids::get_entry(i).generate(tooling_scratch_a, tooling_scratch_b);

      int v = static_cast<int>(temp.vertices.size());
      int f = static_cast<int>(temp.get_face_counts_size());
      int idxs = static_cast<int>(temp.get_faces_size());

      if (v > max_v) {
        max_v = v;
        mv_name = Solids::get_entry(i).name;
      }
      if (f > max_f) {
        max_f = f;
        mf_name = Solids::get_entry(i).name;
      }
      if (idxs > max_i) {
        max_i = idxs;
        mi_name = Solids::get_entry(i).name;
      }
    }

    tooling_scratch_a.reset();
    tooling_scratch_b.reset();

    val stats = val::object();
    stats.set("max_v", max_v);
    stats.set("v_name", val(mv_name));
    stats.set("max_f", max_f);
    stats.set("f_name", val(mf_name));
    stats.set("max_i", max_i);
    stats.set("i_name", val(mi_name));
    return stats;
  }
#endif

  /**
   * @brief Reports the engine arena metrics for the mesh tooling HUD.
   * @return JS object of {usage, high_water_mark, capacity} metrics per arena,
   *         in bytes.
   */
  static val getArenaMetrics() { return collect_arena_metrics(); }
};

/**
 * @brief WASM bridge that bakes a GenerativePalette LUT for the daydream palette
 *        tool, so the tool previews the engine's exact perceptual color math
 *        rather than a hand-ported JS reimplementation that can silently drift.
 * @details The tool owns its own deterministic profile randomization (its stable
 *          PRNG keeps the hue slider from reshuffling structure, and the engine's
 *          global-RNG draws cannot be reproduced anyway), so it resolves the
 *          three (h,s,v) key triples itself and asks here only for the
 *          deterministic OKLCH authoring + gradient evaluation. No global RNG is
 *          touched, so calling this never perturbs a live engine's render stream.
 */
struct PaletteOps {
private:
  // 256 sRGB entries (R,G,B) backing the typed_memory_view bakeLut returns.
  // Sized once at construction so the view's ArrayBuffer never reallocates
  // between calls (same contract as HolosphereEngine::getPixels): JS must read
  // the result before the next bakeLut call.
  std::vector<uint8_t> lut;

public:
  PaletteOps() : lut(256 * 3, 0) {}

  /**
   * @brief Bakes a 256-entry sRGB LUT for a generative palette.
   * @param gradientShape GradientShape as an int (STRAIGHT=0, CIRCULAR=1,
   *        VIGNETTE=2, FALLOFF=3).
   * @param h1 First key hue in [0,255]; s1/v1 its saturation/value.
   * @param s1 First key saturation in [0,255].
   * @param v1 First key value in [0,255].
   * @param h2 Second key hue; s2/v2 its saturation/value.
   * @param s2 Second key saturation.
   * @param v2 Second key value.
   * @param h3 Third key hue; s3/v3 its saturation/value.
   * @param s3 Third key saturation.
   * @param v3 Third key value.
   * @return JS Uint8Array view over 256*3 sRGB bytes; entry i is the palette
   *         sampled at t = i/255. Aliases the shared `lut` buffer (same
   *         memory-view contract as getPixels): read it before the next bakeLut
   *         call, which overwrites the buffer in place.
   */
  val bakeLut(int gradientShape, int h1, int s1, int v1, int h2, int s2, int v2,
              int h3, int s3, int v3) {
    // Out-of-range gradientShape is UB when cast into the enum; clamp and log
    // rather than trap at the JS boundary.
    const int straight = static_cast<int>(GradientShape::STRAIGHT);
    const int falloff = static_cast<int>(GradientShape::FALLOFF);
    if (hs_wasm::gradient_shape_out_of_range(gradientShape, straight, falloff)) {
      hs::log("WASM: bakeLut gradientShape %d out of range — using STRAIGHT",
              gradientShape);
    }
    gradientShape =
        hs_wasm::clamp_gradient_shape(gradientShape, straight, falloff);
    // h/s/v keys arrive as untyped JS ints; clamp to the documented [0,255]
    // rather than letting the uint8_t cast wrap mod 256.
    bool clamped = false;
    auto u8 = [&clamped](int v) -> uint8_t {
      if (hs_wasm::hsv_key_out_of_range(v))
        clamped = true;
      return hs_wasm::clamp_hsv_key(v);
    };
    GenerativePalette pal = GenerativePalette::from_hsv_keys(
        static_cast<GradientShape>(gradientShape), u8(h1), u8(s1), u8(v1),
        u8(h2), u8(s2), u8(v2), u8(h3), u8(s3), u8(v3));
    if (clamped)
      hs::log("WASM: bakeLut hsv key out of [0,255] — clamped");
    for (int i = 0; i < 256; ++i) {
      CRGB c = static_cast<CRGB>(pal.get(i / 255.0f));
      lut[3 * i + 0] = c.r;
      lut[3 * i + 1] = c.g;
      lut[3 * i + 2] = c.b;
    }
    return val(typed_memory_view(lut.size(), lut.data()));
  }
};

/**
 * @brief Packs a Vector into a JS {x,y,z} object for the spline exports.
 */
static val vector_to_xyz(const Vector &r) {
  val v = val::object();
  v.set("x", r.x);
  v.set("y", r.y);
  v.set("z", r.z);
  return v;
}

/**
 * @brief True iff every argument is finite (no NaN/Inf).
 * @details The exported free functions below take raw JS floats straight into
 *          engine math, unlike the class methods (which gate via finite_arg). A
 *          non-finite value can trip an HS_CHECK deep inside the engine (e.g.
 *          Vector::normalized() on the slerp spline path) and abort the *entire*
 *          WASM module — the exact JS-boundary trap the rest of the bridge
 *          avoids. Gate the free functions the same way: reject non-finite input
 *          and return a benign zero result instead of trapping.
 */
static bool all_finite(std::initializer_list<float> args) {
  for (float a : args)
    if (!std::isfinite(a))
      return false;
  return true;
}

/**
 * @brief Evaluates a four-control-point cubic spline (cubic_fast/cubic_slerp)
 *        from the 13 flat floats Embind passes, returning the point as {x,y,z}.
 *
 * Centralizes the {p0..p3} pack and the {x,y,z} marshaling so each binding is
 * a single call and a third interpolator never has to copy it again.
 */
static val eval_cubic_spline(Vector (*fn)(const Vector &, const Vector &,
                                          const Vector &, const Vector &, float),
                             float p0x, float p0y, float p0z, float p1x,
                             float p1y, float p1z, float p2x, float p2y,
                             float p2z, float p3x, float p3y, float p3z,
                             float t, bool reject_degenerate = false) {
  // Reject non-finite input at the boundary: it would abort the slerp normalize.
  if (!all_finite({p0x, p0y, p0z, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z,
                   t})) {
    hs::log("WASM: cubic spline got a non-finite argument — returning zero");
    return vector_to_xyz(Vector(0.0f, 0.0f, 0.0f));
  }
  // slerp normalizes each control point; a finite zero-length one would trip the
  // strict normalized() and abort the module, so reject it at the boundary too.
  if (reject_degenerate) {
    const Vector pts[4] = {{p0x, p0y, p0z}, {p1x, p1y, p1z},
                           {p2x, p2y, p2z}, {p3x, p3y, p3z}};
    for (const Vector &p : pts) {
      if (p.x * p.x + p.y * p.y + p.z * p.z < math::EPS_NORMALIZE_SQ) {
        hs::log("WASM: cubic spline got a degenerate control point — returning zero");
        return vector_to_xyz(Vector(0.0f, 0.0f, 0.0f));
      }
    }
  }
  return vector_to_xyz(fn({p0x, p0y, p0z}, {p1x, p1y, p1z}, {p2x, p2y, p2z},
                          {p3x, p3y, p3z}, t));
}

/**
 * @brief Registers the HolosphereEngine, MeshOps, and spline bindings with
 *        Embind so JavaScript can construct and call them.
 */
EMSCRIPTEN_BINDINGS(holosphere_engine) {
  class_<HolosphereEngine>("HolosphereEngine")
      .constructor<>()
      .function("setResolution", &HolosphereEngine::setResolution)
      .function("setEffect", &HolosphereEngine::setEffect)
      .function("drawFrame", &HolosphereEngine::drawFrame)
      .function("getPixels", &HolosphereEngine::getPixels)
      .function("getBufferLength", &HolosphereEngine::getBufferLength)
      .function("setParameter", &HolosphereEngine::setParameter)
      .function("setAnimationsPaused", &HolosphereEngine::setAnimationsPaused)
      .function("getParameterDefinitions",
                &HolosphereEngine::getParameterDefinitions)
      .function("getParamValues", &HolosphereEngine::getParamValues)
      .function("getArenaMetrics", &HolosphereEngine::getArenaMetrics)
      .function("getEffectSizes", &HolosphereEngine::getEffectSizes)
      .class_function("getSupportedResolutions",
                      &HolosphereEngine::getSupportedResolutions)
      .function("setClip", &HolosphereEngine::setClip)
      .function("getRenderUs", &HolosphereEngine::getRenderUs)
      .function("strobeColumns", &HolosphereEngine::strobeColumns);

  // No public .constructor<>(): all construction goes through fromSolidName so
  // JS cannot wrap an empty mesh past the operator boundary's check_live().
  class_<MeshOpsWrapper>("MeshOps")
      .class_function("clearToolingMemory", &MeshOpsWrapper::clearToolingMemory)
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName)
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
#ifdef HS_WASM_DEV_BINDINGS
      .class_function("getMaxBounds", &MeshOpsWrapper::getMaxBounds)
#endif
      .class_function("getArenaMetrics", &MeshOpsWrapper::getArenaMetrics)
      .function("getVertices", &MeshOpsWrapper::getVertices)
      .function("getFaces", &MeshOpsWrapper::getFaces)
      .function("classifyFaces", &MeshOpsWrapper::classifyFaces)
      // Bound from the same MESHOP_LIST that generates the wrapper methods, plus
      // MESHOP_IRREGULAR_LIST for the hand-written ops.
#define MESHOP_BIND(name) .function(#name, &MeshOpsWrapper::name)
      MESHOP_LIST(MESHOP_BIND, MESHOP_BIND, MESHOP_BIND)
      MESHOP_IRREGULAR_LIST(MESHOP_BIND);
#undef MESHOP_BIND
#undef MESHOP_IRREGULAR_LIST
#undef MESHOP_LIST

  class_<PaletteOps>("PaletteOps")
      .constructor<>()
      .function("bakeLut", &PaletteOps::bakeLut);

  // Spline evaluation — thin wrappers returning val {x,y,z}.
  /**
   * @brief Registers spline_cubic_fast: cubic Bézier point at parameter t over
   *        control points p0..p3.
   * @return JS {x,y,z} object for the evaluated point.
   */
  function("spline_cubic_fast",
           optional_override([](float p0x, float p0y, float p0z, float p1x,
                                float p1y, float p1z, float p2x, float p2y,
                                float p2z, float p3x, float p3y, float p3z,
                                float t) -> val {
             return eval_cubic_spline(&Spline::cubic_fast, p0x, p0y, p0z, p1x,
                                      p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, t);
           }));
  /**
   * @brief Registers spline_cubic_slerp: spherically-interpolated cubic at
   *        parameter t over control points p0..p3.
   * @return JS {x,y,z} object for the evaluated point.
   */
  function("spline_cubic_slerp",
           optional_override([](float p0x, float p0y, float p0z, float p1x,
                                float p1y, float p1z, float p2x, float p2y,
                                float p2z, float p3x, float p3y, float p3z,
                                float t) -> val {
             return eval_cubic_spline(&Spline::cubic_slerp, p0x, p0y, p0z, p1x,
                                      p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, t,
                                      true);
           }));
  /**
   * @brief Registers spline_catmull_rom_tangents: Catmull-Rom tangent estimation.
   * @return JS {cp1:{x,y,z}, cp2:{x,y,z}} giving the two Bézier control points
   *         for the segment start->end.
   */
  function("spline_catmull_rom_tangents",
           optional_override([](float prevx, float prevy, float prevz,
                                float startx, float starty, float startz,
                                float endx, float endy, float endz, float nextx,
                                float nexty, float nextz, float tension) -> val {
             if (!all_finite({prevx, prevy, prevz, startx, starty, startz, endx,
                              endy, endz, nextx, nexty, nextz, tension})) {
               hs::log("WASM: catmull_rom_tangents got a non-finite argument — "
                       "returning zero");
               val v = val::object();
               v.set("cp1", vector_to_xyz(Vector(0.0f, 0.0f, 0.0f)));
               v.set("cp2", vector_to_xyz(Vector(0.0f, 0.0f, 0.0f)));
               return v;
             }
             Vector cp1, cp2;
             Spline::catmull_rom_tangents({prevx, prevy, prevz},
                                          {startx, starty, startz},
                                          {endx, endy, endz},
                                          {nextx, nexty, nextz}, tension, cp1,
                                          cp2);
             val v = val::object();
             v.set("cp1", vector_to_xyz(cp1));
             v.set("cp2", vector_to_xyz(cp2));
             return v;
           }));

  // ── Color / palette / lissajous exports ────────────────────────────────────
  // The real engine math, exported so the JS tool ports can cross-check it.

  // sRGB transfer function (color.js srgbToLinearFloat / linearToSrgbFloat).
  function("srgb_to_linear_float",
           optional_override([](float s) -> float {
             if (all_finite({s}))
               return srgb_to_linear_float(s);
             hs::log("WASM: srgb_to_linear_float got a non-finite argument — "
                     "returning zero");
             return 0.0f;
           }));
  function("linear_to_srgb_float",
           optional_override([](float l) -> float {
             if (all_finite({l}))
               return linear_to_srgb_float(l);
             hs::log("WASM: linear_to_srgb_float got a non-finite argument — "
                     "returning zero");
             return 0.0f;
           }));

  // The interpolated sRGB->16-bit-linear LUT the cosine palette path uses.
  function("srgb_to_linear_interp",
           optional_override([](float s) -> int {
             if (!all_finite({s})) {
               hs::log("WASM: srgb_to_linear_interp got a non-finite argument — "
                       "returning zero");
               return 0;
             }
             return static_cast<int>(srgb_to_linear_interp(s));
           }));

  // OKLab matrices (color.js linearRgbToOklab / oklabToLinearRgb).
  function("linear_rgb_to_oklab",
           optional_override([](float r, float g, float b) -> val {
             if (!all_finite({r, g, b})) {
               hs::log("WASM: linear_rgb_to_oklab got a non-finite argument — "
                       "returning zero");
               val v = val::object();
               v.set("L", 0.0f);
               v.set("a", 0.0f);
               v.set("b", 0.0f);
               return v;
             }
             OKLab lab = linear_rgb_to_oklab(r, g, b);
             val v = val::object();
             v.set("L", lab.L);
             v.set("a", lab.a);
             v.set("b", lab.b);
             return v;
           }));
  function("oklab_to_linear_rgb",
           optional_override([](float L, float a, float b) -> val {
             if (!all_finite({L, a, b})) {
               hs::log("WASM: oklab_to_linear_rgb got a non-finite argument — "
                       "returning zero");
               val v = val::object();
               v.set("r", 0.0f);
               v.set("g", 0.0f);
               v.set("b", 0.0f);
               return v;
             }
             float r, g, bb;
             oklab_to_linear_rgb({L, a, b}, r, g, bb);
             val v = val::object();
             v.set("r", r);
             v.set("g", g);
             v.set("b", bb);
             return v;
           }));

  // HSV -> sRGB integer sextant path (palette_math.js hsvToRgb), via the engine's
  // CRGB(CHSV) constructor. Returns sRGB bytes. The uint8_t cast wraps mod 256,
  // mirroring palette_math.js's `h &= 0xff` byte masking (device CHSV semantics);
  // bakeLut above instead clamps its HSV keys, so the two paths handle
  // out-of-range inputs differently by design.
  function("hsv_to_rgb",
           optional_override([](int h, int s, int v) -> val {
             CRGB c = CRGB(CHSV(static_cast<uint8_t>(h), static_cast<uint8_t>(s),
                                static_cast<uint8_t>(v)));
             val o = val::object();
             o.set("r", static_cast<int>(c.r));
             o.set("g", static_cast<int>(c.g));
             o.set("b", static_cast<int>(c.b));
             return o;
           }));

  // ProceduralPalette cosine formula (palette_math.js ProceduralPalette). Returns
  // the engine's 16-bit linear color so the JS test can pin both the cosine
  // formula and the sRGB->linear interp (paired with srgb_to_linear_interp).
  function("procedural_palette_linear",
           optional_override([](float a0, float a1, float a2, float b0, float b1,
                                float b2, float c0, float c1, float c2, float d0,
                                float d1, float d2, float t) -> val {
             if (!all_finite({a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2,
                              t})) {
               hs::log("WASM: procedural_palette_linear got a non-finite "
                       "argument — returning zero");
               val o = val::object();
               o.set("r", 0);
               o.set("g", 0);
               o.set("b", 0);
               return o;
             }
             ProceduralPalette pal({a0, a1, a2}, {b0, b1, b2}, {c0, c1, c2},
                                   {d0, d1, d2});
             Color4 col = pal.get(t);
             val o = val::object();
             o.set("r", static_cast<int>(col.color.r));
             o.set("g", static_cast<int>(col.color.g));
             o.set("b", static_cast<int>(col.color.b));
             return o;
           }));

  // Lissajous curve (lissajous_math.js lissajous), via geometry.h.
  function("lissajous",
           optional_override([](float m1, float m2, float a, float t) -> val {
             if (!all_finite({m1, m2, a, t})) {
               hs::log("WASM: lissajous got a non-finite argument — returning "
                       "zero");
               return vector_to_xyz(Vector(0.0f, 0.0f, 0.0f));
             }
             return vector_to_xyz(lissajous(m1, m2, a, t));
           }));
}

#endif
