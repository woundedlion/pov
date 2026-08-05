/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/stack.h>
#include "core/engine/effects.h" // Includes all effect headers (triggers REGISTER_EFFECT)
#include "core/engine/effect_registry.h"
#include "core/color/palettes.h" // HS_PROCEDURAL_PALETTE_LIST — named-palette export
#include "core/engine/platform.h"
#include "core/mesh/solids.h"
#include "targets/wasm/param_marshal.h"   // pure, host-tested param marshaling
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <algorithm> // std::fill_n — blank-frame clear in drawFrame
#include <string_view>
#include <cstring>
#include <cstdlib> // std::malloc for the lazily-allocated tooling arenas
#include <cmath>   // std::isfinite — validate MeshOps args at the JS boundary
#include <climits> // INT_MAX — drawFrame pixel-index accumulator bound
#include <initializer_list> // all_finite() variadic-arg gate for free exports
#include <iterator>         // std::size — X-macro roster tables

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

// Worst-case scratch bytes one operator allocates per element of its largest
// stage, counting vertices, faces and flat face indices alike. The heaviest
// operator (meta = kis of dual of ambo) keeps an intermediate mesh, its output and
// its index scratch live in one arena, at 12 bytes per vertex and 6 per index;
// this bound covers that with margin.
static constexpr size_t TOOLING_BYTES_PER_MESH_ELEMENT = 64;

// Largest element count any operator stage may reach. Half-edge construction and
// the topology classifier both narrow face/index counts to uint16_t behind an
// always-on HS_CHECK, so a mesh past this must be rejected at the JS boundary
// rather than allowed to reach one. The scratch arenas hold exactly this many
// elements, so the same ceiling also keeps Arena::allocate's trap out of reach.
static constexpr size_t MAX_MESH_CONNECTIVITY_ELEMENTS = UINT16_MAX;
static_assert(MAX_MESH_CONNECTIVITY_ELEMENTS <=
                  TOOLING_SCRATCH_BYTES / TOOLING_BYTES_PER_MESH_ELEMENT,
              "a stage at the 16-bit ceiling must still fit a scratch arena");

// Tooling-arena bytes a finalized mesh retains per element: one Vector of
// vertex, one uint8_t of side count, one uint16_t of index and the uint16_t
// topology code classifyFaces() later binds into the same arena, with slack for
// the four blocks' alignment padding.
static constexpr size_t TOOLING_ARENA_BYTES_PER_MESH_ELEMENT = 20;
static_assert(sizeof(Vector) + sizeof(uint8_t) + 2 * sizeof(uint16_t) <
                  TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
              "finalized mesh element must fit its predicted arena bytes");

// Widest face a mesh can hold: PolyMesh stores per-face side counts as uint8_t
// and narrow_face_count traps past this, so an operator that would emit a wider
// face must be rejected at the JS boundary.
static constexpr size_t MAX_MESH_FACE_DEGREE = UINT8_MAX;

/**
 * @brief How far one mesh operator grows its input, for the boundary guards.
 * @details Every field is a multiple of an input measurement that some stage of
 *          the operator reaches; see MESHOP_LIST for the per-operator values and
 *          where they come from. `elements` must be nonzero. A zero degree or
 *          valence means the operator emits nothing of that kind; a zero
 *          `valence` additionally skips the valence scan, the one measurement
 *          that costs a pass over the flat index list.
 */
struct MeshOpBounds {
  size_t elements;    /**< Multiple of the largest input element count. */
  size_t face_degree; /**< Multiple of the widest input face's side count. */
  size_t valence;     /**< Multiple of the highest input vertex valence. */
};

// Bumped on every clearToolingMemory(). Each wrapper records the generation it
// was built under and rejects via wrapper_live() if a wipe reclaimed its storage.
static uint32_t tooling_generation = 0;

// Set for the duration of one MeshOps entry point. The tooling scratch arenas
// are module-global and reset() at each entry's head, so a second op entered
// before the first returns would alias the first's scratch and corrupt its
// geometry. Synchronous single-threaded calls never overlap today; the guard
// makes that contract enforced rather than implicit, so a future worker/async
// refactor traps here instead of silently aliasing.
//
// A trap inside an op compiles to wasm `unreachable`, which unwinds nothing, so
// ~ToolingOpGuard() never runs and this stays set. clearToolingMemory() clears
// it, so a JS caller that catches the RuntimeError can recover through there.
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
  HS_CHECK(block != nullptr,
           "tooling arena block allocation of %zu bytes failed", total);
  tooling_arena.rebind(block, TOOLING_ARENA_BYTES);
  tooling_scratch_a.rebind(block + TOOLING_ARENA_BYTES, TOOLING_SCRATCH_BYTES);
  tooling_scratch_b.rebind(block + TOOLING_ARENA_BYTES + TOOLING_SCRATCH_BYTES,
                           TOOLING_SCRATCH_BYTES);
}

/**
 * @brief Adds one arena's {usage, high_water_mark, capacity} entry to a report.
 * @param metrics Report object to extend.
 * @param name Key the entry is stored under.
 * @param arena Arena to measure.
 */
static void add_arena_metrics(val &metrics, const char *name, Arena &arena) {
  val m = val::object();
  m.set("usage", arena.get_offset());
  m.set("high_water_mark", arena.get_high_water_mark());
  m.set("capacity", arena.get_capacity());
  metrics.set(name, m);
}

/**
 * @brief Builds a {usage, high_water_mark, capacity} report for the three engine
 *        arenas.
 * @return JS object mapping each engine arena name to its metrics, in bytes.
 * @details The per-frame HUD path. Every entry costs an embind round-trip, so
 *          this covers only the arenas an engine instance can move; the tooling
 *          arenas are reported by collect_arena_metrics().
 */
static val collect_engine_arena_metrics() {
  val metrics = val::object();
  add_arena_metrics(metrics, "scratch_arena_a", scratch_arena_a);
  add_arena_metrics(metrics, "scratch_arena_b", scratch_arena_b);
  add_arena_metrics(metrics, "persistent_arena", persistent_arena);
  return metrics;
}

/**
 * @brief Builds a {usage, high_water_mark, capacity} report for the three engine
 *        arenas and the three tooling arenas.
 * @return JS object mapping each arena name to its {usage, high_water_mark,
 *         capacity} metrics, in bytes.
 * @details Read on demand through MeshOps.getArenaMetrics(). The tooling scratch
 *          arenas are the regions TOOLING_BYTES_PER_MESH_ELEMENT is sized
 *          against, and an operator that overruns one takes the module down, so
 *          they are reported alongside the rest rather than left unobservable.
 */
static val collect_arena_metrics() {
  val metrics = collect_engine_arena_metrics();
  add_arena_metrics(metrics, "tooling_arena", tooling_arena);
  add_arena_metrics(metrics, "tooling_scratch_a", tooling_scratch_a);
  add_arena_metrics(metrics, "tooling_scratch_b", tooling_scratch_b);
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
    return t;
  }();
  return table;
}

/**
 * @brief Looks up an effect by name in the (W,H) factory.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param name Effect name to look up.
 * @return Pointer to the matching entry, or null if the name is unknown (a
 *         typo'd/stale UI string). The entry belongs to the static per-(W,H)
 *         table and stays valid for the module's lifetime.
 * @details Cheap linear scan used by setEffect() to validate a stale/typo'd UI
 *          string BEFORE it tears down the running effect, so an unknown name is
 *          a transactional no-op rather than a blanked engine; the returned
 *          entry then creates the effect without a second lookup.
 */
template <int W, int H>
const FactoryEntry *find_factory_entry(std::string_view name) {
  for (const auto &entry : get_factory<W, H>())
    if (name == entry.name)
      return &entry;
  return nullptr;
}

// ---------------------------------------------------------------------------
// The (W,H) resolutions the WASM factory can build are the registry's
// HS_RESOLUTIONS rows (core/engine/effect_registry.h), so the runtime dispatch
// here and the per-resolution fill functions the registry generates share one
// list and cannot drift: a resolution the registry can build is dispatchable
// here with no second edit. setResolution()/setEffect()/getEffectSizes() all
// expand it via the X-macro below. A new resolution is one edit, in
// HS_RESOLUTIONS (the effect templates must also be instantiable at that
// <W,H>).
// ---------------------------------------------------------------------------

// Pin every resolution row to the MAX_W×MAX_H pixel-buffer bound.
#define X(W, H)                                                                \
  static_assert((W) <= MAX_W && (H) <= MAX_H,                                  \
                "HS_RESOLUTIONS row exceeds the MAX_W×MAX_H pixel buffer");
HS_RESOLUTIONS(X)
#undef X

/** @brief One HS_RESOLUTIONS row as runtime values. */
struct WasmResolution {
  int w; /**< Canvas width in pixels. */
  int h; /**< Canvas height in pixels. */
};

// The rows as data, so the constructor can bootstrap on the first one instead of
// naming a preset that HS_RESOLUTIONS may later drop.
static constexpr WasmResolution WASM_RESOLUTIONS[] = {
#define X(W, H) {(W), (H)},
    HS_RESOLUTIONS(X)
#undef X
};
static_assert(std::size(WASM_RESOLUTIONS) > 0,
              "HS_RESOLUTIONS must carry a row to bootstrap on");

// The registered effect names, in HS_EFFECT_LIST order. Only the first is used
// (the constructor's bootstrap), for the same reason as WASM_RESOLUTIONS.
static constexpr const char *WASM_EFFECT_NAMES[] = {
#define X(name) #name,
    HS_EFFECT_LIST(X)
#undef X
};
static_assert(std::size(WASM_EFFECT_NAMES) == HS_EFFECT_COUNT,
              "HS_EFFECT_LIST and HS_EFFECT_COUNT must agree");

/**
 * @brief Invokes f.operator()<W,H>() for the single HS_RESOLUTIONS row
 *        matching the runtime (w,h).
 * @tparam F Type of the templated callable.
 * @param w Runtime canvas width in pixels to match against the row list.
 * @param h Runtime canvas height in pixels to match against the row list.
 * @param f A C++20 templated callable (e.g. `[]<int W, int H>(){...}`) run with
 *          the matching row's dimensions as compile-time template arguments.
 * @return true iff a row matched and f was invoked; false otherwise.
 * @details One shared expansion of the resolution list for every runtime
 *          dispatch site (setEffect validate/create, resolution checks), so they
 *          can never drift and a new per-resolution step lands in exactly one
 *          place.
 */
template <typename F> static bool dispatch_resolution(int w, int h, F &&f) {
#define X(W, H)                                                                \
  if (w == (W) && h == (H)) {                                                  \
    f.template operator()<(W), (H)>();                                         \
    return true;                                                               \
  }
  HS_RESOLUTIONS(X)
#undef X
  return false;
}

/**
 * @brief Reports whether (w,h) is a resolution the factory can build.
 * @param w Candidate canvas width in pixels.
 * @param h Candidate canvas height in pixels.
 * @return true iff (w,h) is one of the HS_RESOLUTIONS rows.
 */
static bool wasm_resolution_supported(int w, int h) {
  return dispatch_resolution(w, h, []<int W, int H>() {});
}

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

// True while a HolosphereEngine is constructed-but-not-deleted. The engine is a
// singleton: its Effect aliases shared static buffers and its arenas are
// module-global, so a second instance would corrupt the first's frames.
static bool engine_alive = false;

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
    engine_alive = true;
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
        setResolution(WASM_RESOLUTIONS[0].w, WASM_RESOLUTIONS[0].h);
    HS_CHECK(bootstrap_resolution_set,
             "the first HS_RESOLUTIONS row must be dispatchable here");
    const bool bootstrap_effect_set = setEffect(WASM_EFFECT_NAMES[0]);
    HS_CHECK(bootstrap_effect_set,
             "the first HS_EFFECT_LIST entry must be registered and buildable "
             "at the first HS_RESOLUTIONS row");
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
   * @return true if the resolution is now active; false if the request was
   *         rejected (unsupported size) and the previous valid state was kept.
   *         A request matching the active resolution also returns true but is a
   *         pure no-op — nothing is torn down — so true alone does not imply the
   *         teardown described below happened.
   * @details A successful resolution change tears down the current effect (a new
   *          one cannot be carried across pixel dimensions), so the caller must
   *          call setEffect() again before the next drawFrame() or it renders a
   *          blank frame. Callers that use a sub-canvas clip (segmented rendering)
   *          must likewise re-apply setClip() after a successful setResolution():
   *          a prior clip was expressed in the old resolution's pixel bounds and
   *          is not rescaled here, so it must be recomputed for the new
   *          dimensions. An outstanding getPixels() view is left aliasing live
   *          memory at the previous resolution's length rather than detached, so
   *          it must be re-fetched and its length compared against
   *          getBufferLength().
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

    if (current_effect) {
      current_effect = nullptr;
      ++param_generation;
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
   *          prior valid state alive rather than blanking the engine. A
   *          successful switch instantiates a fresh effect at the full-canvas
   *          clip and with default parameter values, so callers that use a
   *          sub-canvas clip (segmented rendering) or tuned parameters must
   *          re-apply setClip() and setParameter() after every successful
   *          setEffect(). The engine-owned animation pause state is retained.
   */
  bool setEffect(const std::string &name) {
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
      return false;
    }

    hs::log("WASM: setEffect called with %s", name.c_str());
    if (!entry) {
      hs::log("WASM: setEffect unknown effect '%s' — keeping current effect",
              name.c_str());
      return false;
    }

    // Per-load RNG stream, mirroring the device's per-epoch reseed: load 0
    // (the constructor's bootstrap effect) keeps the identity 1337 stream;
    // each later load derives a fresh seed. Replica-safe with no protocol
    // change — workers process identical message streams, so their counters
    // agree by construction.
    hs::random().seed(hs::epoch_seed(effect_loads++));

    current_effect.reset();
    configure_arenas_default(); // Reset before init so effects can override

    stack_paint_canary(); // reset stack HWM by repainting unused region

    dispatch_resolution(pixel_width, pixel_height, []<int W, int H>() {
      init_geometry_luts<W, H>(); // eager-fill LUTs before the first frame
    });
    current_effect = entry->creator();
    current_effect->setAnimationsPaused(animations_paused);
    current_effect->init();
    ++param_generation;
    const size_t init_hwm = stack_high_water_mark();
    if (init_hwm > init_stack_peak)
      init_stack_peak = init_hwm;
    hs::log("WASM: init stack HWM = %u bytes", (unsigned)init_hwm);
    stack_paint_canary();
    return true;
  }

  /**
   * @brief Restricts rendering to a clip band for the current effect.
   * @param x0 Inclusive left column of the clip band in [0, pixel_width].
   * @param x1 Exclusive right column, with x0 <= x1 <= pixel_width.
   * @param y0 Inclusive top row of the clip band in [0, pixel_height].
   * @param y1 Exclusive bottom row of the clip band, with y0 <= y1 <= pixel_height.
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
   *          there aborts the whole WASM module. Segment workers always pass
   *          valid, ordered, in-range bounds. A cross-segment stateful effect
   *          (Effect::needs_full_frame()) keeps the full-canvas clip instead of
   *          narrowing to the band — see docs/segmented_stateful_effects_spec.md.
   */
  ClipSetResult setClip(int x0, int x1, int y0, int y1) {
    if (!current_effect)
      return ClipSetResult::NO_EFFECT;
    // Reject malformed bounds from the untyped JS boundary (negatives would feed
    // ClipRegion's modulo arithmetic); reject-and-return, never trap.
    if (!hs_wasm::clip_bounds_valid(x0, x1, y0, y1, pixel_width,
                                    pixel_height)) {
      hs::log("WASM: setClip bounds out of range (x0=%d,x1=%d,y0=%d,y1=%d) — "
              "ignored",
              x0, x1, y0, y1);
      return ClipSetResult::INVALID_BOUNDS;
    }
    // Cross-segment stateful effects must render the FULL canvas in every worker
    // (a band-clipped worker has stale cv.prev outside its band, so trails seam);
    // keep the full clip. See docs/segmented_stateful_effects_spec.md.
    if (current_effect->needs_full_frame())
      return ClipSetResult::FULL_FRAME_KEPT;
    current_effect->set_clip(y0, y1, x0, x1);
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
   *          next effect when no effect is currently loaded.
   */
  void setAnimationsPaused(bool paused) {
    animations_paused = paused;
    if (current_effect)
      current_effect->setAnimationsPaused(paused);
  }

  /**
   * @brief Sets near-pole azimuthal shading decimation.
   * @param aggressiveness Columns per shade are this over sin(colatitude);
   *        0 disables. Non-finite and negative inputs clamp to 0.
   * @details The physically-neutral setting is 1.0: at that value one shade
   *          covers the columns sharing a physical LED footprint. Exposed so
   *          the value can be tuned against real hardware.
   */
  void setPoleLod(float aggressiveness) {
    pole_lod_aggressiveness =
        hs_wasm::clamp_pole_lod_aggressiveness(aggressiveness);
  }

  /** @brief Current near-pole azimuthal decimation aggressiveness. */
  float getPoleLod() const { return pole_lod_aggressiveness; }

  /**
   * @brief Builds the GUI's parameter descriptor list.
   * @return JS array with one {name, value, animated, readonly, (+min/max for
   *         floats)} object per param in the effect's declaration order; empty
   *         array when no effect is set.
   */
  val getParameterDefinitions() {
    if (!current_effect)
      return val::array();

    val result = val::array();
    // Both streams walk the effect's registered ParamList in order
    // (param_marshal.h); the definitions and getParamValues() agree only
    // because that list is fixed after init, with getParamGeneration() covering
    // replacement of the effect itself.
    hs_wasm::collect_param_views(*current_effect, param_views);

    int i = 0;
    for (const auto &v : param_views) {
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
    // The no-reallocation contract both paths ride on: the ctor reserved
    // MAX_PARAMS and clear() retains that capacity, which also keeps the
    // zero-length view's backing pointer valid.
    HS_CHECK(param_values.capacity() >= MAX_PARAMS,
             "param_values capacity %zu below the reserved %zu",
             param_values.capacity(), MAX_PARAMS);
    if (!current_effect) {
      // Empty Float32Array (not a JS Array) so callers get a consistent typed
      // view whether or not an effect is set.
      param_values.clear();
      return val(typed_memory_view(param_values.size(), param_values.data()));
    }

    // Same order as getParameterDefinitions().
    hs_wasm::fill_param_values(*current_effect, param_values);
    return val(typed_memory_view(param_values.size(), param_values.data()));
  }

  /**
   * @brief Identity token joining a getParameterDefinitions() snapshot to a
   *        later getParamValues() read.
   * @return A fresh value after every effect load or teardown.
   * @details Neither stream carries a per-effect identity, so nothing but this
   *          counter distinguishes a value read that describes the snapshotted
   *          effect from one taken after a switch. Parameter counts repeat
   *          across the roster, so a length comparison is not a substitute. Pin
   *          this at snapshot time and compare it against a re-read beside each
   *          value read; a change means the snapshot is stale and the
   *          definitions must be rebuilt before the values are applied.
   */
  uint32_t getParamGeneration() const { return param_generation; }

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
  /** Channels per pixel in the readback buffer (linear RGB triples). */
  static constexpr int CHANNELS = 3;

  std::unique_ptr<Effect>
      current_effect;                 /**< Currently active effect, or null. */
  std::vector<uint16_t> pixel_buffer; /**< 16-bit linear RGB readback buffer. */
  std::vector<float> param_values;    /**< Backing store for getParamValues. */
  std::vector<hs_wasm::ParamView>
      param_views;           /**< Scratch for getParameterDefinitions. */
  int pixel_width = 0;       /**< Active canvas width in pixels. */
  int pixel_height = 0;      /**< Active canvas height in pixels. */
  uint32_t effect_loads = 0; /**< Effect loads so far; epoch for the per-load
                                 RNG seed (0 = the constructor's bootstrap). */
  uint32_t param_generation =
      0; /**< Identity token for the current effect or no-effect state. */
  bool animations_paused = false; /**< Pause state applied to every effect. */
};

// ==========================================================================================
// MESH OPS BINDINGS
// ==========================================================================================

/**
 * @brief Why the most recent MeshOps call answered null.
 * @details A bare null collapses reasons that demand opposite caller actions —
 *          shrinking the op chain versus calling clearToolingMemory() — so the
 *          reason is recorded here and read back via MeshOps.getLastResult().
 *          Exposed to JS as the Module.MeshOpResult embind enum; compare
 *          against its values, never by truthiness.
 */
enum class MeshOpResult {
  OK,                    /**< The call produced a result. */
  UNKNOWN_NAME,          /**< No registry entry carries that name. */
  CONNECTIVITY_OVERFLOW, /**< A stage would pass the 16-bit element ceiling. */
  FACE_DEGREE_OVERFLOW,  /**< A stage would emit a face past the 8-bit side
                              count. */
  ARENA_EXHAUSTED,       /**< The result would not fit tooling_arena's
                              remaining bytes. */
  NON_FINITE_ARG,        /**< An operator argument was NaN or infinite. */
  ANGLE_OUT_OF_DOMAIN,   /**< An angle argument sat outside its operator's
                              domain. */
  STALE_WRAPPER,         /**< The wrapper's storage was reclaimed by a
                              clearToolingMemory(). */
};

// Outcome of the most recent MeshOps call that could answer null
// (getLastResult()).
static MeshOpResult last_mesh_op_result = MeshOpResult::OK;

/**
 * @brief JS-facing wrapper around a PolyMesh and the Conway/Goldberg operators.
 * @details Named to avoid collision with the MeshOps namespace. Each wrapper's
 *          mesh is built into the tooling arena and records the generation it was
 *          built under, so a wipe via clearToolingMemory() is detected by
 *          wrapper_live().
 */
struct MeshOpsWrapper {
private:
  PolyMesh mesh; /**< The wrapped mesh, stored in the tooling arena. */
  /**
   * Generation of the tooling arena this mesh was built into; compared against
   * the live counter on every use (see wrapper_live()).
   */
  uint32_t generation = tooling_generation;

public:
  /**
   * @brief Constructs a wrapper taking ownership of an existing mesh.
   * @param m Mesh to move into this wrapper.
   */
  MeshOpsWrapper(PolyMesh &&m) : mesh(std::move(m)) {}

private:
  /**
   * @brief Reports whether this wrapper outlived a clearToolingMemory() wipe.
   * @return true while its mesh still owns live arena storage; false — after
   *         logging and recording STALE_WRAPPER — once a wipe reclaimed it.
   * @details A stale wrapper's mesh aliases reclaimed arena storage, which
   *          release builds would otherwise read back as silently wrong
   *          geometry. Called at every entry point that touches `mesh`, which
   *          rejects rather than traps: a JS caller holding a wrapper across an
   *          interleaved clearToolingMemory() is an ordering slip the page can
   *          recover from, not an engine invariant violation.
   */
  bool wrapper_live() const {
    if (generation == tooling_generation)
      return true;
    hs::log("WASM: MeshOps wrapper used after clearToolingMemory() — ignored");
    last_mesh_op_result = MeshOpResult::STALE_WRAPPER;
    return false;
  }

public:
  /**
   * @brief Resets all tooling arenas to empty and invalidates live wrappers.
   * @details Reclaims the storage behind every live wrapper and bumps the
   *          generation so any wrapper built before this wipe is rejected on
   *          next use (wrapper_live). JS-callable.
   *
   *          Despite the name, this does NOT shrink the module's linear memory:
   *          the 16 MB tooling block is retained for the module's lifetime and
   *          only its arena bump-pointers are reset. A JS caller will not see
   *          memory usage drop; this frees the arenas for reuse, not the OS.
   *
   *          Also clears the re-entrancy flag, which a trap inside an op leaves
   *          set (see tooling_op_active). No op is in flight once JS regains
   *          control, so this is the recovery entry point after a trap.
   */
  static void clearToolingMemory() {
    tooling_op_active = false;
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
   * @return Owning pointer to the new wrapper, or null for an unknown name or a
   *         solid that would not fit what is left of tooling_arena;
   *         getLastResult() names which.
   * @details Rejects an unknown name at the untrusted JS boundary rather than
   *          tripping get_by_name()'s fail-fast HS_CHECK and aborting the module.
   *          Generates into the scratch arenas and prices the finalized copy
   *          against tooling_arena's remaining bytes before committing it, since
   *          that arena accumulates one finalized mesh per live wrapper until
   *          clearToolingMemory() and Arena::allocate traps when it runs out.
   */
  static std::unique_ptr<MeshOpsWrapper>
  fromSolidName(const std::string &name) {
    last_mesh_op_result = MeshOpResult::OK;
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: fromSolidName unknown solid '%s' — ignored", name.c_str());
      last_mesh_op_result = MeshOpResult::UNKNOWN_NAME;
      return nullptr;
    }
    ToolingOpGuard guard;
    ensure_tooling_arenas();
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    const PolyMesh generated =
        entry->generate(tooling_scratch_a, tooling_scratch_b);
    if (hs_wasm::mesh_op_output_over_arena(
            generated.vertices.size(), generated.get_face_counts_size(),
            generated.get_faces_size(), 1, TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
            tooling_arena.get_offset(), tooling_arena.get_capacity())) {
      hs::log("WASM: fromSolidName '%s' does not fit the tooling arena (%zu of "
              "%zu bytes used) — ignored; call clearToolingMemory() to reclaim "
              "it, which invalidates every live mesh",
              name.c_str(), tooling_arena.get_offset(),
              tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return nullptr;
    }
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(generated, tooling_arena));
  }

  /**
   * @brief Returns the mesh vertices as a JS Float32Array.
   * @return Float32Array of flattened [x,y,z] triples, copied out of the mesh,
   *         or null if a clearToolingMemory() reclaimed this wrapper's storage;
   *         getLastResult() then reports STALE_WRAPPER.
   */
  val getVertices() const {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
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
   *         Null if a clearToolingMemory() reclaimed this wrapper's storage;
   *         getLastResult() then reports STALE_WRAPPER.
   */
  val getFaces() const {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
    size_t total = 0;
    for (size_t i = 0; i < mesh.get_face_counts_size(); ++i)
      total += mesh.get_face_counts_data()[i];
    HS_CHECK(total == mesh.get_faces_size(),
             "getFaces: face_counts sum disagrees with the flat index count");
    val out = val::object();
    out.set("indices", val::global("Uint16Array")
                           .new_(val(typed_memory_view(
                               mesh.get_faces_size(), mesh.get_faces_data()))));
    out.set("counts",
            val::global("Uint8Array")
                .new_(val(typed_memory_view(mesh.get_face_counts_size(),
                                            mesh.get_face_counts_data()))));
    return out;
  }

  /**
   * @brief Classifies faces by topology and returns the per-face codes.
   * @return JS Int32Array of one topology code per face, copied out of the
   *         mesh's now-populated topology buffer, or null when this wrapper's
   *         storage was reclaimed, the mesh is past
   *         MAX_MESH_CONNECTIVITY_ELEMENTS, or its topology block would not fit
   *         what is left of tooling_arena; getLastResult() names which. Null
   *         rather than an empty array so a caller can tell "no classification"
   *         from "no faces" with a plain truthiness test.
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
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return val::null();
    if (hs_wasm::tooling_mesh_over_ceiling(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), MAX_MESH_CONNECTIVITY_ELEMENTS)) {
      hs::log("WASM: classifyFaces mesh of %zu verts / %zu faces / %zu indices "
              "is past the %zu-element 16-bit connectivity range — ignored",
              mesh.vertices.size(), mesh.get_face_counts_size(),
              mesh.get_faces_size(), MAX_MESH_CONNECTIVITY_ELEMENTS);
      last_mesh_op_result = MeshOpResult::CONNECTIVITY_OVERFLOW;
      return val::null();
    }
    ensure_tooling_arenas();
    if (hs_wasm::mesh_op_output_over_arena(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), 1, TOOLING_ARENA_BYTES_PER_MESH_ELEMENT,
            tooling_arena.get_offset(), tooling_arena.get_capacity())) {
      hs::log("WASM: classifyFaces topology block does not fit the tooling "
              "arena (%zu of %zu bytes used) — ignored; call "
              "clearToolingMemory() to reclaim it, which invalidates every "
              "live mesh",
              tooling_arena.get_offset(), tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return val::null();
    }
    ToolingOpGuard guard;
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

private:
  /**
   * @brief Highest vertex valence in this wrapper's mesh.
   * @return Faces meeting at the most-incident vertex, or 0 for an empty mesh.
   * @details Counts incidences in tooling_scratch_a (4 bytes per vertex, released
   *          before the operator runs), so the whole measurement is one pass over
   *          the flat index list on top of clearing the counters.
   */
  size_t max_vertex_valence() {
    const size_t verts = mesh.vertices.size();
    if (verts == 0)
      return 0;
    tooling_scratch_a.reset();
    uint32_t *incidence = tooling_scratch_a.allocate_n<uint32_t>(verts);
    const size_t valence = hs_wasm::mesh_max_vertex_valence(
        mesh.get_faces_data(), mesh.get_faces_size(), incidence, verts);
    tooling_scratch_a.reset();
    return valence;
  }

  /**
   * @brief Runs a mesh operator and wraps the result.
   * @tparam Op Callable of signature (const PolyMesh&, Arena&, Arena&) ->
   *         PolyMesh.
   * @param bounds Operator's growth factors (see MESHOP_LIST).
   * @param op Operator to run against this wrapper's mesh.
   * @return Owning pointer to a new wrapper holding the finalized result mesh, or
   *         null if this wrapper's storage was reclaimed, some stage of this
   *         operator would pass MAX_MESH_CONNECTIVITY_ELEMENTS or
   *         MAX_MESH_FACE_DEGREE, or its output would not fit what is left of
   *         tooling_arena; getLastResult() names which.
   * @details Captures the shared operator boilerplate: reset both tooling scratch
   *          arenas, run the op into a fresh PolyMesh, finalize it into
   *          tooling_arena, and hand back a new wrapper. An input the operator
   *          would grow past the engine's 16-bit connectivity range is rejected
   *          here, at the untrusted JS boundary, rather than reaching
   *          build_half_edge_mesh's fail-fast trap partway through a composition;
   *          likewise an input holding a face or a vertex valence the operator
   *          would widen past narrow_face_count's 8-bit side count, and a result
   *          that would overrun tooling_arena, which accumulates one finalized
   *          mesh per chained wrapper until clearToolingMemory().
   */
  template <typename Op>
  std::unique_ptr<MeshOpsWrapper> apply(MeshOpBounds bounds, Op &&op) {
    last_mesh_op_result = MeshOpResult::OK;
    if (!wrapper_live())
      return nullptr;
    if (hs_wasm::mesh_op_expansion_over_ceiling(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), bounds.elements,
            MAX_MESH_CONNECTIVITY_ELEMENTS)) {
      hs::log("WASM: MeshOps input mesh of %zu verts / %zu faces / %zu indices "
              "expands %zux, past the %zu-element 16-bit connectivity range — "
              "ignored",
              mesh.vertices.size(), mesh.get_face_counts_size(),
              mesh.get_faces_size(), bounds.elements,
              MAX_MESH_CONNECTIVITY_ELEMENTS);
      last_mesh_op_result = MeshOpResult::CONNECTIVITY_OVERFLOW;
      return nullptr;
    }
    ensure_tooling_arenas();
    // Before max_vertex_valence(), the first tooling-scratch use on this path.
    ToolingOpGuard guard;
    if (hs_wasm::mesh_op_output_over_arena(
            mesh.vertices.size(), mesh.get_face_counts_size(),
            mesh.get_faces_size(), bounds.elements,
            TOOLING_ARENA_BYTES_PER_MESH_ELEMENT, tooling_arena.get_offset(),
            tooling_arena.get_capacity())) {
      hs::log("WASM: MeshOps result does not fit the tooling arena (%zu of %zu "
              "bytes used) — ignored; call clearToolingMemory() to reclaim it, "
              "which invalidates every live mesh",
              tooling_arena.get_offset(), tooling_arena.get_capacity());
      last_mesh_op_result = MeshOpResult::ARENA_EXHAUSTED;
      return nullptr;
    }
    const size_t face_degree = hs_wasm::mesh_max_face_degree(
        mesh.get_face_counts_data(), mesh.get_face_counts_size());
    const size_t valence = bounds.valence == 0 ? 0 : max_vertex_valence();
    if (hs_wasm::mesh_op_face_degree_overflows(
            face_degree, valence, bounds.face_degree, bounds.valence,
            MAX_MESH_FACE_DEGREE)) {
      hs::log("WASM: MeshOps input mesh (widest face %zu sides x%zu, highest "
              "valence %zu x%zu) would emit a face past the %zu-side limit — "
              "ignored",
              face_degree, bounds.face_degree, valence, bounds.valence,
              MAX_MESH_FACE_DEGREE);
      last_mesh_op_result = MeshOpResult::FACE_DEGREE_OVERFLOW;
      return nullptr;
    }
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
    last_mesh_op_result = MeshOpResult::NON_FINITE_ARG;
    return false;
  }

public:
/**
 * @brief Defines a zero-argument Conway/Goldberg operator method.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details The generated method runs MeshOps::name(mesh) via apply() and returns
 *          a new wrapper holding the result.
 */
#define MESHOP_0(name, elements, degree, valence)                              \
  std::unique_ptr<MeshOpsWrapper> name() {                                     \
    return apply({elements, degree, valence},                                  \
                 [](const PolyMesh &m, Arena &a, Arena &b) {                   \
                   return MeshOps::name(m, a, b);                              \
                 });                                                           \
  }
/**
 * @brief Defines a one-float-argument operator whose argument is a [0,1]
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Rejects a non-finite arg (finite_arg), then clamps the fraction to
 *          [0,1] (logging when it changed) so a direct/API caller passing a
 *          finite out-of-range value stays within the operator's documented
 *          domain and cannot trip truncate's or bevel's always-on HS_CHECK and
 *          abort the whole module.
 */
#define MESHOP_1U(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(float arg) {                            \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    if (hs_wasm::unit_fraction_out_of_range(arg))                              \
      hs::log("WASM: MeshOps::%s clamped %g to [0,1]", #name, arg);            \
    float t = hs_wasm::clamp_unit_fraction(arg);                               \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

/**
 * @brief Defines a one-float-argument operator whose argument is a [0,1)
 *        fraction, clamped at the JS boundary.
 * @param name MeshOps operator name; becomes the generated method name.
 * @param elements Multiple of the largest input element count (see MESHOP_LIST).
 * @param degree Multiple of the widest input face's side count.
 * @param valence Multiple of the highest input vertex valence.
 * @details Like MESHOP_1U, but these operators assert `t < 1.0f`, so 1 is
 *          outside the domain and clamping to it would land on the trap instead
 *          of avoiding it.
 */
#define MESHOP_1H(name, elements, degree, valence)                             \
  std::unique_ptr<MeshOpsWrapper> name(float arg) {                            \
    if (!finite_arg(arg, #name))                                               \
      return nullptr;                                                          \
    if (hs_wasm::half_open_fraction_out_of_range(arg))                         \
      hs::log("WASM: MeshOps::%s clamped %g to [0,1)", #name, arg);            \
    float t = hs_wasm::clamp_half_open_fraction(arg);                          \
    return apply({elements, degree, valence},                                  \
                 [t](const PolyMesh &m, Arena &a, Arena &b) {                  \
                   return MeshOps::name(m, a, b, t);                           \
                 });                                                           \
  }

  /**
 * @brief Single source of truth for the Conway/Goldberg operator roster and each
 *        operator's growth factors.
 * @param _OP0  Macro applied to each zero-argument operator.
 * @param _OP1U Macro applied to each [0,1]-fraction operator.
 * @param _OP1H Macro applied to each [0,1)-fraction operator.
 * @details Expanded twice: with MESHOP_0/MESHOP_1U/MESHOP_1H to generate the
 *          wrapper methods (below), and with MESHOP_BIND to generate the embind
 *          .function() bindings (in EMSCRIPTEN_BINDINGS), so an operator cannot
 *          be added to one site and silently left unreachable from the other.
 *          Each fraction operator's macro matches the domain its always-on
 *          engine trap asserts: truncate and bevel accept 1 and use _OP1U;
 *          chamfer and expand assert t < 1 and use _OP1H. relax, hankin,
 *          and snub have bespoke signatures/validation (snub takes two float
 *          controls, and clamps its inset to [0,1) for the same trap), so their
 *          wrapper methods are hand-written; their names
 *          live in MESHOP_IRREGULAR_LIST below and their bindings expand from
 *          it, so a new irregular op is bound the moment it joins the list.
 *
 *          The trailing arguments are the operator's MeshOpBounds, in order.
 *
 *          `elements` is the largest multiple of the input's biggest element count
 *          that any of its intermediate or output stages reaches. Read off the
 *          output bindings in core/mesh/conway.h with I as the input flat index
 *          count — dual I, ambo 2I, kis/truncate 3I, expand/chamfer 4I, snub 5I —
 *          and multiplied through for the compositions: gyro = d(snub) 5I,
 *          needle = k(d) 3I, zip = d(k) 3I, meta = k(d(a)) 6I, bevel = t(a) 6I.
 *          Every stage's vertex and face count stays under its operator's index
 *          expansion, so one factor bounds all three counts.
 *
 *          `degree` and `valence` are the multiples that reach narrow_face_count,
 *          from the three emitters: emit_shrunk_face and emit_primary_faces widen
 *          a face to verts_per_side x its side count (1 for ambo/expand/chamfer/
 *          snub, 2 for truncate), emit_vertex_orbit_faces emits one face per
 *          vertex at its valence, and Hankin doubles both (star faces 2x a face's
 *          sides, rosettes 2x a valence). kis emits only triangles, so it needs
 *          neither measurement. Compositions inherit their widest stage: needle's
 *          dual stage is valence-wide, zip's is 2x valence (its kis stage doubles
 *          valence first), meta's and gyro's ambo/snub stages are 1x both, and
 *          bevel's truncate stage doubles what its ambo stage already widened.
 *
 *          Irregular ops pass theirs at their own apply() call: relax {1, 1, 0}
 *          (topology preserving), hankin {4, 2, 2}, snub {5, 1, 1}.
 */
  // clang-format off
#define MESHOP_LIST(_OP0, _OP1U, _OP1H)                                         \
  _OP0(kis, 3, 0, 0) _OP0(ambo, 2, 1, 1) _OP0(gyro, 5, 1, 1)                    \
  _OP0(dual, 1, 0, 1) _OP0(meta, 6, 1, 1) _OP0(needle, 3, 0, 1)                 \
  _OP0(zip, 3, 1, 2)                                                            \
  _OP1U(truncate, 3, 2, 1) _OP1U(bevel, 6, 2, 2)                                \
  _OP1H(chamfer, 4, 1, 0) _OP1H(expand, 4, 1, 1)
// clang-format on

// Irregular ops: hand-written wrapper methods (custom signatures/validation),
// enumerated here so their embind bindings expand from one list.
#define MESHOP_IRREGULAR_LIST(_) _(relax) _(hankin) _(snub)

  MESHOP_LIST(MESHOP_0, MESHOP_1U, MESHOP_1H)

#undef MESHOP_0
#undef MESHOP_1U
#undef MESHOP_1H

  /**
   * Engine-enforced ceiling on relax smoothing passes: an independent
   * defense-in-depth limit for direct/API callers that bypass the editor's
   * 500-pass slider cap.
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
  std::unique_ptr<MeshOpsWrapper> relax(int iterations) {
    int clamped =
        hs_wasm::clamp_relax_iterations(iterations, MAX_RELAX_ITERATIONS);
    if (clamped != iterations)
      hs::log("WASM: MeshOps::relax clamped %d iterations to %d", iterations,
              clamped);
    iterations = clamped;
    return apply({1, 1, 0},
                 [iterations](const PolyMesh &m, Arena &a, Arena &b) {
                   return MeshOps::relax(m, a, b, iterations);
                 });
  }

  /**
   * Inclusive upper bound of the Hankin contact-angle domain: at pi/2 the
   * contact rays leave the edge perpendicular to it, and past that they tilt
   * back into the neighbouring face, mirroring an angle already in domain.
   */
  static constexpr float MAX_HANKIN_ANGLE = PI_F / 2.0f;

  /**
   * @brief Applies the Hankin interlace operator to the mesh.
   * @param radians Interlace angle in radians (the unit MeshOps::hankin
   *        expects), in the operator's [0, MAX_HANKIN_ANGLE] domain.
   * @return Owning pointer to a new wrapper holding the result, or null if the
   *         angle is non-finite or out of domain; getLastResult() names which.
   * @details Explicit (not a MESHOP_* macro) so the radians unit contract the JS
   *          caller relies on is carried here. The angle is rejected rather than
   *          clamped: MeshOps::hankin aliases an out-of-domain angle onto an
   *          in-domain pattern, so clamping would hand back geometry the caller
   *          did not ask for.
   */
  std::unique_ptr<MeshOpsWrapper> hankin(float radians) {
    if (!finite_arg(radians, "hankin"))
      return nullptr;
    if (hs_wasm::hankin_angle_out_of_range(radians, MAX_HANKIN_ANGLE)) {
      hs::log("WASM: MeshOps::hankin angle %g outside [0, %g] — ignored",
              radians, MAX_HANKIN_ANGLE);
      last_mesh_op_result = MeshOpResult::ANGLE_OUT_OF_DOMAIN;
      return nullptr;
    }
    return apply({4, 2, 2}, [radians](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::hankin(m, a, b, radians);
    });
  }

  /**
   * @brief Applies the chiral snub operator with explicit inset and twist.
   * @param t Inset factor of each face toward its centroid, clamped to [0, 1)
   *          (its documented domain) at the JS boundary.
   * @param twist Per-face rotation about the face normal, in radians (0 = none);
   *          unbounded, so only finiteness is checked.
   * @return Owning pointer to a new wrapper, or null if either arg is non-finite.
   * @details Explicit (not a MESHOP_* macro) because MeshOps::snub takes TWO
   *          float controls, which neither the zero-arg nor the one-float
   *          generator can express. Binding it via MESHOP_0 hardcodes the
   *          (0.5, 0.0) defaults and leaves both controls unreachable from JS;
   *          this 2-arg form exposes them to the solids editor.
   */
  std::unique_ptr<MeshOpsWrapper> snub(float t, float twist) {
    if (!finite_arg(t, "snub") || !finite_arg(twist, "snub"))
      return nullptr;
    if (hs_wasm::half_open_fraction_out_of_range(t))
      hs::log("WASM: MeshOps::snub clamped t=%g to [0,1)", t);
    float ct = hs_wasm::clamp_half_open_fraction(t);
    return apply({5, 1, 1}, [ct, twist](const PolyMesh &m, Arena &a, Arena &b) {
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

  /**
   * @brief Maps an authored Solids::Op to the editor's lowercase op string.
   * @param op Authored operator from a recipe step.
   * @return The op name matching the solids editor's vocabulary.
   */
  static const char *op_name(Solids::Op op) {
    switch (op) {
    case Solids::Op::TRUNCATE:
      return "truncate";
    case Solids::Op::EXPAND:
      return "expand";
    case Solids::Op::SNUB:
      return "snub";
    case Solids::Op::CHAMFER:
      return "chamfer";
    case Solids::Op::HANKIN:
      return "hankin";
    case Solids::Op::RELAX:
      return "relax";
    case Solids::Op::KIS:
      return "kis";
    case Solids::Op::DUAL:
      return "dual";
    case Solids::Op::AMBO:
      return "ambo";
    case Solids::Op::BEVEL:
      return "bevel";
    case Solids::Op::GYRO:
      return "gyro";
    case Solids::Op::META:
      return "meta";
    case Solids::Op::NEEDLE:
      return "needle";
    case Solids::Op::ZIP:
      return "zip";
    }
    HS_CHECK(false, "unhandled Solids::Op %d in op_name", static_cast<int>(op));
  }

  /**
   * @brief Returns a solid's authored recipe chain for the editor.
   * @param name Registry name to look up.
   * @return JS object {seed: string, ops: [{op: string, param, twist}]}, or
   *         null for an unknown name or an entry without a recipe.
   * @details Pure table read: no arenas, no wrapper, no clearToolingMemory()
   *          pairing. Params cross in engine-native units (radians for hankin,
   *          raw t, relax iteration count), matching the MeshOps op bindings.
   *          An unknown name is rejected at the untrusted JS boundary
   *          (fromSolidName's precedent) rather than trapping; a recipe-less
   *          known entry is the normal not-morphable case and returns null
   *          without logging.
   */
  static val getRecipe(const std::string &name) {
    const Solids::Entry *entry = Solids::find_entry(name);
    if (!entry) {
      hs::log("WASM: getRecipe unknown solid '%s' — ignored", name.c_str());
      return val::null();
    }
    if (!entry->recipe)
      return val::null();
    const Solids::Recipe &recipe = *entry->recipe;
    // recipe.seed indexes simple_registry specifically (get_entry spans all
    // three registries and would misresolve it).
    val out = val::object();
    out.set("seed", val(Solids::simple_registry[recipe.seed].name));
    val ops = val::array();
    for (size_t i = 0; i < recipe.count; ++i) {
      const Solids::OpStep &step = recipe.steps[i];
      val item = val::object();
      item.set("op", val(op_name(step.op)));
      item.set("param", step.param);
      item.set("twist", step.twist);
      ops.set(i, item);
    }
    out.set("ops", ops);
    return out;
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

  /**
   * @brief Reports why the most recent MeshOps call answered null.
   * @return OK when that call produced a result, otherwise its rejection reason.
   * @details Covers fromSolidName, getVertices, getFaces, classifyFaces and the
   *          operator methods. Read it immediately after the null; the next such
   *          call overwrites it.
   */
  static MeshOpResult getLastResult() { return last_mesh_op_result; }
};

// 256 sRGB entries (R,G,B) backing the typed_memory_view PaletteOps::bakeLut
// returns. File-scope and fixed-size so the view neither outlives its storage
// (a PaletteOps delete() would free a member buffer while the JS view stayed
// non-empty and readable over reusable heap) nor reallocates between calls.
static constexpr size_t PALETTE_LUT_BYTES = 256 * 3;
static uint8_t palette_lut[PALETTE_LUT_BYTES];

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
   *         sampled at t = i/255. Aliases the module-lifetime palette_lut buffer
   *         (same memory-view contract as getPixels): read it before the next
   *         bakeLut call on any PaletteOps, which overwrites the buffer in place.
   */
  val bakeLut(int gradientShape, int h1, int s1, int v1, int h2, int s2, int v2,
              int h3, int s3, int v3) {
    // The JS caller passes GradientShape by its integer value; pin the mapping so
    // a reorder in color.h can't silently remap shapes.
    static_assert(static_cast<int>(GradientShape::STRAIGHT) == 0 &&
                      static_cast<int>(GradientShape::CIRCULAR) == 1 &&
                      static_cast<int>(GradientShape::VIGNETTE) == 2 &&
                      static_cast<int>(GradientShape::FALLOFF) == 3,
                  "GradientShape integer values are part of the JS ABI");
    // Out-of-range gradientShape is UB when cast into the enum; clamp and log
    // rather than trap at the JS boundary.
    const int straight = static_cast<int>(GradientShape::STRAIGHT);
    const int falloff = static_cast<int>(GradientShape::FALLOFF);
    if (hs_wasm::gradient_shape_out_of_range(gradientShape, straight,
                                             falloff)) {
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
      palette_lut[3 * i + 0] = c.r;
      palette_lut[3 * i + 1] = c.g;
      palette_lut[3 * i + 2] = c.b;
    }
    return val(typed_memory_view(PALETTE_LUT_BYTES, palette_lut));
  }
};

/**
 * @brief Packs a Vector into a JS {x,y,z} object for the free-function exports.
 * @param r Vector to pack.
 * @return JS object with the x, y and z components as numbers.
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
 * @param args Floats to test, as a brace-init list at the call site.
 * @return true if every argument is finite; false if any is NaN or Inf.
 * @details The exported free functions below take raw JS floats straight into
 *          engine math, unlike the class methods (which gate via finite_arg). A
 *          non-finite value can trip an HS_CHECK deep inside the engine and
 *          abort the *entire* WASM module — the exact JS-boundary trap the rest
 *          of the bridge avoids. Gate the free functions the same way: reject
 *          non-finite input and return a benign zero result instead of trapping.
 */
static bool all_finite(std::initializer_list<float> args) {
  for (float a : args)
    if (!std::isfinite(a))
      return false;
  return true;
}

/**
 * @brief Registers the HolosphereEngine and MeshOps bindings with Embind so
 *        JavaScript can construct and call them.
 */
EMSCRIPTEN_BINDINGS(holosphere_engine) {
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

  class_<HolosphereEngine>("HolosphereEngine")
      .constructor<>()
      .function("setResolution", &HolosphereEngine::setResolution)
      .function("setEffect", &HolosphereEngine::setEffect)
      .function("drawFrame", &HolosphereEngine::drawFrame)
      .function("getPixels", &HolosphereEngine::getPixels)
      .function("getBufferLength", &HolosphereEngine::getBufferLength)
      .function("setParameter", &HolosphereEngine::setParameter)
      .function("setAnimationsPaused", &HolosphereEngine::setAnimationsPaused)
      .function("setPoleLod", &HolosphereEngine::setPoleLod)
      .function("getPoleLod", &HolosphereEngine::getPoleLod)
      .function("getParameterDefinitions",
                &HolosphereEngine::getParameterDefinitions)
      .function("getParamValues", &HolosphereEngine::getParamValues)
      .function("getParamGeneration", &HolosphereEngine::getParamGeneration)
      .function("getArenaMetrics", &HolosphereEngine::getArenaMetrics)
      .function("getEffectSizes", &HolosphereEngine::getEffectSizes)
      .class_function("getSupportedResolutions",
                      &HolosphereEngine::getSupportedResolutions)
      .function("setClip", &HolosphereEngine::setClip)
      .function("strobeColumns", &HolosphereEngine::strobeColumns);

  enum_<MeshOpResult>("MeshOpResult")
      .value("OK", MeshOpResult::OK)
      .value("UNKNOWN_NAME", MeshOpResult::UNKNOWN_NAME)
      .value("CONNECTIVITY_OVERFLOW", MeshOpResult::CONNECTIVITY_OVERFLOW)
      .value("FACE_DEGREE_OVERFLOW", MeshOpResult::FACE_DEGREE_OVERFLOW)
      .value("ARENA_EXHAUSTED", MeshOpResult::ARENA_EXHAUSTED)
      .value("NON_FINITE_ARG", MeshOpResult::NON_FINITE_ARG)
      .value("ANGLE_OUT_OF_DOMAIN", MeshOpResult::ANGLE_OUT_OF_DOMAIN)
      .value("STALE_WRAPPER", MeshOpResult::STALE_WRAPPER);

  // No public .constructor<>(): all construction goes through fromSolidName so
  // JS cannot wrap an empty mesh past the operator boundary's wrapper_live().
  class_<MeshOpsWrapper>("MeshOps")
      .class_function("clearToolingMemory", &MeshOpsWrapper::clearToolingMemory)
      .class_function("getLastResult", &MeshOpsWrapper::getLastResult)
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName)
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
      .class_function("getRecipe", &MeshOpsWrapper::getRecipe)
#ifdef HS_WASM_DEV_BINDINGS
      .class_function("getMaxBounds", &MeshOpsWrapper::getMaxBounds)
#endif
      .class_function("getArenaMetrics", &MeshOpsWrapper::getArenaMetrics)
      .function("getVertices", &MeshOpsWrapper::getVertices)
      .function("getFaces", &MeshOpsWrapper::getFaces)
      .function("classifyFaces", &MeshOpsWrapper::classifyFaces)
  // Bound from the same MESHOP_LIST that generates the wrapper methods, plus
  // MESHOP_IRREGULAR_LIST for the hand-written ops.
// Variadic so it can take both MESHOP_LIST's (name, expansion) entries and
// MESHOP_IRREGULAR_LIST's bare names.
#define MESHOP_BIND(name, ...) .function(#name, &MeshOpsWrapper::name)
          MESHOP_LIST(MESHOP_BIND, MESHOP_BIND, MESHOP_BIND)
              MESHOP_IRREGULAR_LIST(MESHOP_BIND);
#undef MESHOP_BIND
#undef MESHOP_IRREGULAR_LIST
#undef MESHOP_LIST

  class_<PaletteOps>("PaletteOps")
      .constructor<>()
      .function("bakeLut", &PaletteOps::bakeLut);

  // ── Color / palette / geometry exports ─────────────────────────────────────
  // The real engine math, exported so the JS tool ports can cross-check it.

  // sRGB transfer function (color.js srgbToLinearFloat / linearToSrgbFloat).
  function("srgb_to_linear_float", optional_override([](float s) -> float {
             if (all_finite({s}))
               return srgb_to_linear_float(s);
             hs::log("WASM: srgb_to_linear_float got a non-finite argument — "
                     "returning zero");
             return 0.0f;
           }));
  function("linear_to_srgb_float", optional_override([](float l) -> float {
             if (all_finite({l}))
               return linear_to_srgb_float(l);
             hs::log("WASM: linear_to_srgb_float got a non-finite argument — "
                     "returning zero");
             return 0.0f;
           }));

  // The interpolated sRGB->16-bit-linear LUT the cosine palette path uses.
  function("srgb_to_linear_interp", optional_override([](float s) -> int {
             if (!all_finite({s})) {
               hs::log(
                   "WASM: srgb_to_linear_interp got a non-finite argument — "
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
  function("hsv_to_rgb", optional_override([](int h, int s, int v) -> val {
             CRGB c =
                 CRGB(CHSV(static_cast<uint8_t>(h), static_cast<uint8_t>(s),
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
  function(
      "procedural_palette_linear",
      optional_override([](float a0, float a1, float a2, float b0, float b1,
                           float b2, float c0, float c1, float c2, float d0,
                           float d1, float d2, float t) -> val {
        if (!all_finite({a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2, t})) {
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

  // The named procedural palettes (palette_math.js NAMED_PROCEDURAL_PALETTES),
  // in core/color/palettes.h declaration order. Enumerated from the same X-macro
  // the Palettes:: instances are declared from, so the browser tool's mirror is
  // compared against the literals the engine compiles, not a second hand-copy.
  function("named_procedural_palettes", optional_override([]() -> val {
             val out = val::array();
             int index = 0;
             auto push = [&](const char *name, std::array<float, 3> a,
                             std::array<float, 3> b, std::array<float, 3> c,
                             std::array<float, 3> d) {
               val entry = val::object();
               entry.set("name", std::string(name));
               const std::array<float, 3> *coeff[] = {&a, &b, &c, &d};
               const char *keys[] = {"a", "b", "c", "d"};
               for (int k = 0; k < 4; ++k) {
                 val vec = val::array();
                 for (int ch = 0; ch < 3; ++ch)
                   vec.set(ch, (*coeff[k])[ch]);
                 entry.set(keys[k], vec);
               }
               out.set(index++, entry);
             };
#define HS_EXPORT_PALETTE(name, A, B, C, D)                                    \
  push(#name, {HS_PALETTE_VEC3 A}, {HS_PALETTE_VEC3 B}, {HS_PALETTE_VEC3 C},   \
       {HS_PALETTE_VEC3 D});
             HS_PROCEDURAL_PALETTE_LIST(HS_EXPORT_PALETTE)
#undef HS_EXPORT_PALETTE
             return out;
           }));

  // GenerativePalette's seeded profile resolution (palette_math.js
  // GenerativePalette's constructor, calcHues and seededRandInt). Returns the
  // nine HSV key components the engine's profile constructor authors its keys
  // from, so the JS mirror is checked against the engine rather than a golden.
  // Feed them straight back into PaletteOps.bakeLut for the matching gradient.
  //
  // The seed is clamped into [0,255] rather than passed through: a negative
  // seed would send every draw to the global RNG and advance the shared hue
  // cursor, so this export can never perturb a live render stream.
  function(
      "generative_palette_hsv_keys",
      optional_override([](int harmony, int brightness, int sat,
                           int seed) -> val {
        // The JS caller passes each enum by its integer value; pin the mappings
        // so a reorder in color.h can't silently remap profiles.
        static_assert(static_cast<int>(HarmonyType::TRIADIC) == 0 &&
                          static_cast<int>(HarmonyType::SPLIT_COMPLEMENTARY) ==
                              1 &&
                          static_cast<int>(HarmonyType::COMPLEMENTARY) == 2 &&
                          static_cast<int>(HarmonyType::ANALOGOUS) == 3,
                      "HarmonyType integer values are part of the JS ABI");
        static_assert(
            static_cast<int>(BrightnessProfile::ASCENDING) == 0 &&
                static_cast<int>(BrightnessProfile::DESCENDING) == 1 &&
                static_cast<int>(BrightnessProfile::FLAT) == 2 &&
                static_cast<int>(BrightnessProfile::BELL) == 3 &&
                static_cast<int>(BrightnessProfile::CUP) == 4,
            "BrightnessProfile integer values are part of the JS ABI");
        static_assert(
            static_cast<int>(SaturationProfile::PASTEL) == 0 &&
                static_cast<int>(SaturationProfile::MID) == 1 &&
                static_cast<int>(SaturationProfile::VIBRANT) == 2,
            "SaturationProfile integer values are part of the JS ABI");
        // Out-of-range indices are UB when cast into the enum; fold to the first
        // enumerator and log rather than trap at the JS boundary.
        auto pick = [](int index, int count, const char *name) -> int {
          if (index < 0 || index >= count) {
            hs::log("WASM: generative_palette_hsv_keys %s %d out of range — "
                    "using 0",
                    name, index);
            return 0;
          }
          return index;
        };
        harmony = pick(harmony, 4, "harmony");
        brightness = pick(brightness, 5, "brightness");
        sat = pick(sat, 3, "sat");
        if (hs_wasm::hsv_key_out_of_range(seed))
          hs::log("WASM: generative_palette_hsv_keys seed out of [0,255] — "
                  "clamped");
        GenerativePalette::HsvKeys k = GenerativePalette::resolve_hsv_keys(
            static_cast<HarmonyType>(harmony),
            static_cast<BrightnessProfile>(brightness),
            static_cast<SaturationProfile>(sat), hs_wasm::clamp_hsv_key(seed));
        val o = val::object();
        o.set("h1", static_cast<int>(k.h1));
        o.set("s1", static_cast<int>(k.s1));
        o.set("v1", static_cast<int>(k.v1));
        o.set("h2", static_cast<int>(k.h2));
        o.set("s2", static_cast<int>(k.s2));
        o.set("v2", static_cast<int>(k.v2));
        o.set("h3", static_cast<int>(k.h3));
        o.set("s3", static_cast<int>(k.s3));
        o.set("v3", static_cast<int>(k.v3));
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

  // Mobius sphere map (mobius_transforms.js coefficients), via transformers.h.
  // The eight coefficient floats are taken in the order mobiusCodeString emits
  // them, so the tool's MobiusParams initializer ordering is pinned too.
  function(
      "mobius_transform",
      optional_override([](float x, float y, float z, float ar, float ai,
                           float br, float bi, float cr, float ci, float dr,
                           float di) -> val {
        if (!all_finite({x, y, z, ar, ai, br, bi, cr, ci, dr, di})) {
          hs::log("WASM: mobius_transform got a non-finite argument — "
                  "returning zero");
          return vector_to_xyz(Vector(0.0f, 0.0f, 0.0f));
        }
        return vector_to_xyz(mobius_transform(
            Vector(x, y, z), MobiusParams(ar, ai, br, bi, cr, ci, dr, di)));
      }));
}

#endif
