/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/stack.h>
#include "effects.h" // Includes all effect headers (triggers REGISTER_EFFECT)
#include "core/effect_registry.h"
#include "platform.h"
#include <string_view>
#include <cstring>

// ---- Stack canary painting for high water mark tracking ----
static constexpr uint8_t STACK_CANARY = 0xCD;

/// Paint the unused portion of the stack (current SP down to stack end) with a
/// canary byte. Safe to call at any time — only touches memory below the
/// current stack pointer.
static void stack_paint_canary() {
  uintptr_t sp = emscripten_stack_get_current();
  uintptr_t end = emscripten_stack_get_end();
  if (sp > end) {
    std::memset(reinterpret_cast<void *>(end), STACK_CANARY, sp - end);
  }
}

/// Scan from the stack end upward to find the first overwritten canary byte.
/// Returns the number of bytes that have been touched (high water mark).
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
// `std::array<ParamDef, 32>` backing `Effect::ParamList` (core/canvas.h). Used
// to pre-reserve the getParamValues() backing store so it never reallocates.
static constexpr size_t MAX_PARAMS = 32;

// Dedicated arena for JavaScript Tools
static uint8_t tooling_buf[8 * 1024 * 1024];
static uint8_t tooling_scratch_buf_a[4 * 1024 * 1024];
static uint8_t tooling_scratch_buf_b[4 * 1024 * 1024];
Arena tooling_arena(tooling_buf, sizeof(tooling_buf));
Arena tooling_scratch_a(tooling_scratch_buf_a, sizeof(tooling_scratch_buf_a));
Arena tooling_scratch_b(tooling_scratch_buf_b, sizeof(tooling_scratch_buf_b));

// Build a concrete factory table from the self-registering entries
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

template <int W, int H>
std::unique_ptr<Effect> create_effect(std::string_view name) {
  const auto &factory = get_factory<W, H>();
  for (const auto &entry : factory) {
    if (name == entry.name)
      return entry.creator();
  }
  // Unknown name = typo'd/stale UI string. Substituting factory[0] would render
  // a different effect with no signal; instead surface it and return null —
  // setEffect()'s `if (currentEffect)` guard makes a null a safe no-op.
  hs::log("WASM: create_effect: unknown effect name (no effect created)");
  return nullptr;
}

// ---------------------------------------------------------------------------
// Single source of truth for the (W,H) resolutions the WASM factory can build.
// setResolution()/setEffect()/getEffectSizes() all dispatch through this list
// via the X-macro below, so the supported set can never drift between them.
// To support a new resolution, add one row here (the effect templates must also
// be instantiable at that <W,H>).
// ---------------------------------------------------------------------------
#define HS_WASM_RESOLUTIONS(X)                                                  \
  X(96, 20)                                                                     \
  X(288, 144)

static bool wasm_resolution_supported(int w, int h) {
#define X(W, H)                                                                 \
  if (w == (W) && h == (H))                                                     \
    return true;
  HS_WASM_RESOLUTIONS(X)
#undef X
  return false;
}

class HolosphereEngine {
public:
  HolosphereEngine() {
    // Paint stack canary for HWM tracking
    stack_paint_canary();

    // SSOT guard: the self-registering effect count must match the static roster
    // in core/effects.h (HS_EFFECT_LIST / HS_EFFECT_COUNT). A mismatch means an
    // effect was added to the roster without REGISTER_EFFECT (or registered
    // without a roster entry) — the live set here and the native smoke suite
    // (which generates one case per roster entry) would silently diverge. Trap
    // at startup; this is a cold one-time check with no per-frame cost.
    HS_CHECK(EffectRegistry::entries().size() ==
             static_cast<size_t>(HS_EFFECT_COUNT));

    // Pre-size the JS-facing readback buffers ONCE to their maximum extent so
    // their backing storage never moves and the steady-state render/sync path
    // performs no allocation. This is the core of the WASM memory-view contract
    // (see getPixels()): under ALLOW_MEMORY_GROWTH=1 any heap reallocation
    // detaches the ArrayBuffer behind a typed_memory_view, so the buffers
    // exposed to JS must be stable for the lifetime of the engine.
    pixelBuffer.assign(MAX_W * MAX_H * 3, 0); // 16-bit linear RGB; never resized
    paramValues.reserve(MAX_PARAMS);          // never reallocated past this

    // Initialize with a valid default effect. JS overrides this almost
    // immediately (daydream defaults to IslamicStars / the URL ?effect=), but a
    // real registered name keeps the engine valid for that first instant and for
    // any headless/tool use.
    setResolution(96, 20);
    setEffect("DistortedRing");
  }

  // Returns true if the resolution is now active, false if the request was
  // rejected (unsupported size) and the previous valid state was kept.
  bool setResolution(int w, int h) {
    if (w == pixel_width && h == pixel_height)
      return true; // already at this (valid) resolution

    // Reject sizes the factory can't build rather than switching to them and
    // nulling currentEffect — that previously left the engine rendering blank
    // with no signal to JS. Keep the prior valid resolution/effect alive and
    // report the failure so the caller can surface it.
    if (!wasm_resolution_supported(w, h)) {
      hs::log("WASM: Unsupported resolution %dx%d — ignored", w, h);
      return false;
    }

    pixel_width = w;
    pixel_height = h;

    // NOTE: pixelBuffer is pre-sized to MAX_W*MAX_H*3 in the constructor and is
    // deliberately NEVER resized here. Resizing could move its backing store
    // (and/or grow the WASM heap), detaching every outstanding getPixels()
    // view. getPixels() instead returns a view over just the active
    // pixel_width*pixel_height*3 prefix of this stable buffer.

    // Re-create current effect if exists
    if (currentEffect) {
      currentEffect = nullptr;
    }
    return true;
  }

  void setEffect(std::string name) {
    // hs::log is printf-style: pass name as an arg, never as the format string
    // (an effect name containing '%' would otherwise read nonexistent varargs).
    hs::log("WASM: setEffect called with %s", name.c_str());

    currentEffect.reset();
    configure_arenas_default();

    // Reset stack HWM by repainting unused region
    stack_paint_canary();

    bool created = false;
#define X(W, H)                                                                \
  if (pixel_width == (W) && pixel_height == (H)) {                             \
    currentEffect = create_effect<W, H>(name);                                 \
    created = true;                                                            \
  }
    HS_WASM_RESOLUTIONS(X)
#undef X
    if (!created) {
      // Unreachable in practice: setResolution() only admits supported sizes.
      hs::log("WASM: Unsupported resolution for factory!");
      return;
    }
    if (currentEffect) {
      currentEffect->init();
      // Log init stack HWM, then repaint to isolate render HWM
      char hwm_buf[80];
      snprintf(hwm_buf, sizeof(hwm_buf), "WASM: init stack HWM = %u bytes",
               (unsigned)stack_high_water_mark());
      hs::log(hwm_buf);
      stack_paint_canary();
    }
  }

  bool setClip(int y0, int y1, int x0, int x1) {
    if (!currentEffect)
      return false;
    // Clip bounds cross the untyped JS boundary. Reject malformed input that
    // would otherwise feed negatives into ClipRegion's modulo arithmetic
    // (constants.h render_x_*), yielding a wrong clip band instead of a clean
    // segment. Mirror setResolution: reject-and-return rather than trap, since a
    // trap at the JS boundary aborts the whole WASM module. Segment workers
    // always pass valid, in-range, ordered bounds, so this rejects only
    // malformed external calls.
    if (!(x0 >= 0 && y0 >= 0 && x0 <= x1 && x1 <= pixel_width && y0 <= y1 &&
          y1 <= pixel_height)) {
      hs::log("WASM: setClip bounds out of range (%d,%d,%d,%d) — ignored", y0,
              y1, x0, x1);
      return false;
    }
    currentEffect->set_clip(y0, y1, x0, x1);
    return true;
  }

  void drawFrame() {
    if (!currentEffect)
      return;

    currentEffect->render_us = 0.0;
    currentEffect->draw_frame();
    currentEffect->advance_display();

    // Output 16-bit Linear values directly
    int idx = 0;
    for (int y = 0; y < pixel_height; y++) {
      for (int x = 0; x < pixel_width; x++) {
        const Pixel &p = currentEffect->get_pixel(x, y);
        pixelBuffer[idx++] = p.r;
        pixelBuffer[idx++] = p.g;
        pixelBuffer[idx++] = p.b;
      }
    }
  }

  double getRenderUs() {
    return currentEffect ? currentEffect->render_us : 0.0;
  }

  // Expose the raw pixel buffer to JS as a zero-copy Uint16Array view.
  //
  // WASM memory-view contract: the returned view aliases WASM linear memory, it
  // is NOT a copy. With ALLOW_MEMORY_GROWTH=1, any subsequent heap growth
  // detaches the underlying ArrayBuffer and leaves this view zero-length
  // (buffer.byteLength === 0). Callers MUST therefore re-fetch the view after
  // anything that may allocate (resolution/effect change) and may only cache it
  // across frames while guarding for detachment. The backing vector is pre-sized
  // once and never reallocated, so the only possible detachment source is heap
  // growth elsewhere. daydream.js::refreshPixelView mirrors this expectation.
  val getPixels() {
    // View spans only the active resolution's pixels (R,G,B per pixel) within
    // the stable MAX_W*MAX_H*3 backing buffer.
    return val(typed_memory_view(pixel_width * pixel_height * 3,
                                 pixelBuffer.data()));
  }

  int getBufferLength() { return pixel_width * pixel_height * 3; }

  void setParameter(std::string name, float value) {
    if (currentEffect) {
      currentEffect->updateParameter(name.c_str(), value);
    }
  }

  void setAnimationsPaused(bool paused) {
    if (currentEffect) {
      currentEffect->setAnimationsPaused(paused);
    }
  }

  val getParameterDefinitions() {
    if (!currentEffect)
      return val::array();

    val result = val::array();
    const auto &params = currentEffect->getParameters();

    int i = 0;
    for (const auto &def : params) {
      val entry = val::object();
      entry.set("name", val(def.name));

      if (def.is_bool()) {
        // Emit JS boolean so frontend typeof check works
        entry.set("value", val(def.get() > 0.5f));
      } else {
        entry.set("value", def.get());
        entry.set("min", def.min);
        entry.set("max", def.max);
      }
      // Animation-driven params surface as auto-pausing sliders in the GUI;
      // read-only params are shown live but disabled for editing.
      entry.set("animated", val(def.animated));
      entry.set("readonly", val(def.readonly));
      result.set(i++, entry);
    }
    return result;
  }

  val getParamValues() {
    if (!currentEffect)
      return val::array();

    const auto &params = currentEffect->getParameters();
    paramValues.clear(); // retains capacity reserved (MAX_PARAMS) in ctor

    for (const auto &def : params) {
      paramValues.push_back(def.get());
    }
    // Same memory-view contract as getPixels(): this view aliases WASM memory
    // and must be consumed before the next allocation. paramValues never
    // reallocates here (params.size() <= MAX_PARAMS), so emitting it triggers no
    // heap growth that could detach other outstanding views.
    return val(typed_memory_view(paramValues.size(), paramValues.data()));
  }

  val getArenaMetrics() {
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

    // Stack metrics (same format as arenas)
    {
      uintptr_t base = emscripten_stack_get_base();
      uintptr_t end = emscripten_stack_get_end();
      uintptr_t sp = emscripten_stack_get_current();
      val m = val::object();
      m.set("usage", static_cast<unsigned>(base - sp));
      m.set("high_water_mark", static_cast<unsigned>(stack_high_water_mark()));
      m.set("capacity", static_cast<unsigned>(base - end));
      metrics.set("stack", m);
    }

    return metrics;
  }

  template <int W, int H> val get_effect_sizes_helper() {
    val s = val::object();
    const auto &factory = get_factory<W, H>();
    for (const auto &entry : factory)
      s.set(std::string(entry.name), static_cast<int>(entry.size));
    return s;
  }

  val getEffectSizes() {
#define X(W, H)                                                                \
  if (pixel_width == (W) && pixel_height == (H))                               \
    return get_effect_sizes_helper<W, H>();
    HS_WASM_RESOLUTIONS(X)
#undef X
    return val::object(); // unsupported/uninitialized — empty map
  }

private:
  std::unique_ptr<Effect> currentEffect;
  std::vector<uint16_t> pixelBuffer; // 16-bit
  std::vector<float> paramValues;    // Backing store for getParamValues
  int pixel_width = 0;
  int pixel_height = 0;
};

// ==========================================================================================
// MESH OPS BINDINGS
// ==========================================================================================

#include "solids.h"

// Wrapper to avoid collision with MeshOps namespace
struct MeshOpsWrapper {
  PolyMesh mesh;

  MeshOpsWrapper() {}
  MeshOpsWrapper(PolyMesh &&m) : mesh(std::move(m)) {}

  // Call this from JS whenever you want to wipe the UI memory clean!
  static void clearToolingMemory() {
    tooling_arena.reset();
    tooling_arena.reset_high_water_mark();
    tooling_scratch_a.reset();
    tooling_scratch_a.reset_high_water_mark();
    tooling_scratch_b.reset();
    tooling_scratch_b.reset_high_water_mark();
  }

  // Factory
  static std::unique_ptr<MeshOpsWrapper> fromSolid(int index) {
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::get(
        tooling_arena, tooling_scratch_a, tooling_scratch_b, index));
  }

  static std::unique_ptr<MeshOpsWrapper> fromSolidName(std::string name) {
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::get_by_name(
        tooling_arena, tooling_scratch_a, tooling_scratch_b, name));
  }

  static std::unique_ptr<MeshOpsWrapper> fromData(val vertices, val faces) {
    PolyMesh m;

    // Vertices: Float32Array [x, y, z, ...]
    std::vector<float> vData = convertJSArrayToNumberVector<float>(vertices);
    for (size_t i = 0; i < vData.size(); i += 3) {
      m.vertices.emplace_back(vData[i], vData[i + 1], vData[i + 2]);
    }

    // Faces: Flat Int32Array with -1 delimiter
    std::vector<int> fData = convertJSArrayToNumberVector<int>(faces);

    // First pass: count faces to allocate face_counts
    int num_faces = 0;
    for (int idx : fData) {
      if (idx == -1)
        num_faces++;
    }
    // If last face doesn't end with -1 but has data
    if (!fData.empty() && fData.back() != -1)
      num_faces++;

    m.faces.bind(tooling_arena,
                 fData.size() - num_faces); // Exact size without delimiters
    m.face_counts.bind(tooling_arena, num_faces);

    int current_count = 0;
    for (int idx : fData) {
      if (idx == -1) {
        if (current_count > 0) {
          m.face_counts.push_back((uint8_t)current_count);
          current_count = 0;
        }
      } else {
        m.faces.push_back(idx);
        current_count++;
      }
    }
    if (current_count > 0) {
      m.face_counts.push_back((uint8_t)current_count);
    }

    return std::make_unique<MeshOpsWrapper>(std::move(m));
  }

  // Accessors for JS
  val getVertices() const {
    std::vector<float> data;
    data.reserve(mesh.vertices.size() * 3);
    for (const auto &v : mesh.vertices) {
      data.push_back(v.x);
      data.push_back(v.y);
      data.push_back(v.z);
    }
    // Create JS Float32Array from memory view (copying data)
    return val::global("Float32Array")
        .new_(val(typed_memory_view(data.size(), data.data())));
  }

  val getFaces() const {
    val faces_arr = val::array();
    int flat_idx = 0;
    for (size_t i = 0; i < mesh.face_counts.size(); ++i) {
      val face = val::array();
      int count = mesh.face_counts[i];
      for (int c = 0; c < count; ++c) {
        face.call<void>("push", mesh.faces[flat_idx++]);
      }
      faces_arr.set(i, face);
    }
    return faces_arr;
  }

  val getFaceCounts() const {
    // Array backing is flat so we can pass memory view directly
    return val::global("Uint8Array")
        .new_(val(typed_memory_view(mesh.face_counts.size(),
                                    mesh.face_counts.data())));
  }

  val getFlatIndices() const {
    // Array backing is flat so we can pass memory view directly
    return val::global("Int32Array")
        .new_(val(typed_memory_view(mesh.faces.size(), mesh.faces.data())));
  }

  val classifyFaces() {
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    MeshOps::classify_faces_by_topology(mesh, tooling_scratch_a,
                                        tooling_scratch_b, tooling_arena);
    return val::global("Int32Array")
        .new_(
            val(typed_memory_view(mesh.topology.size(), mesh.topology.data())));
  }

  // --- Conway/Goldberg operators -------------------------------------------
  // Every operator has the same shape: reset both tooling scratch arenas, run
  // the op into a fresh PolyMesh, finalize it into tooling_arena, and hand back
  // a new wrapper. apply() captures that boilerplate so each operator is just
  // its MeshOps call; the MESHOP_* macros remove the remaining stutter. This is
  // tooling — invoked from the JS mesh editor, never the render loop — and
  // every lambda inlines, so there is no added cost.
  template <typename Op>
  std::unique_ptr<MeshOpsWrapper> apply(Op &&op) const {
    tooling_scratch_a.reset();
    tooling_scratch_b.reset();
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        op(mesh, tooling_scratch_a, tooling_scratch_b), tooling_arena));
  }

#define MESHOP_0(name)                                                         \
  std::unique_ptr<MeshOpsWrapper> name() const {                              \
    return apply(                                                              \
        [](const PolyMesh &m, Arena &a, Arena &b) { return MeshOps::name(m, a, b); });           \
  }
#define MESHOP_1(name, T)                                                      \
  std::unique_ptr<MeshOpsWrapper> name(T arg) const {                        \
    return apply([arg](const PolyMesh &m, Arena &a, Arena &b) {               \
      return MeshOps::name(m, a, b, arg);                                     \
    });                                                                        \
  }

  MESHOP_0(kis)
  MESHOP_0(ambo)
  MESHOP_0(gyro)
  MESHOP_0(snub)
  MESHOP_0(dual)
  MESHOP_0(meta)
  MESHOP_0(needle)
  MESHOP_0(zip)
  MESHOP_1(truncate, float)
  MESHOP_1(expand, float)
  MESHOP_1(chamfer, float)
  MESHOP_1(bitruncate, float)
  MESHOP_1(bevel, float)
  MESHOP_1(relax, int)

#undef MESHOP_0
#undef MESHOP_1

  // Hankin shares one MeshOps entry point but with two unit conventions, so it
  // stays explicit rather than fitting the MESHOP_* shape.
  std::unique_ptr<MeshOpsWrapper> hankin(float angle) const {
    return apply([angle](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::hankin(m, a, b, angle * (PI_F / 180.0f));
    });
  }
  std::unique_ptr<MeshOpsWrapper> hankin_rad(float radians) const {
    return apply([radians](const PolyMesh &m, Arena &a, Arena &b) {
      return MeshOps::hankin(m, a, b, radians);
    });
  }
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

  static val getMaxBounds() {
    int max_v = 0;
    int max_f = 0;
    int max_i = 0;
    const char *mv_name = "";
    const char *mf_name = "";
    const char *mi_name = "";

    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      tooling_scratch_a.reset();
      tooling_scratch_b.reset();
      PolyMesh temp =
          Solids::get(tooling_arena, tooling_scratch_a, tooling_scratch_b, i);

      int v = temp.vertices.size();
      int f = temp.face_counts.size();
      int idxs = temp.faces.size();

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

    val stats = val::object();
    stats.set("max_v", max_v);
    stats.set("v_name", val(mv_name));
    stats.set("max_f", max_f);
    stats.set("f_name", val(mf_name));
    stats.set("max_i", max_i);
    stats.set("i_name", val(mi_name));
    return stats;
  }

  static val getArenaMetrics() {
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
};

// Expose to JavaScript
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
      .function("setClip", &HolosphereEngine::setClip)
      .function("getRenderUs", &HolosphereEngine::getRenderUs);

  class_<MeshOpsWrapper>("MeshOps")
      .constructor<>()
      .class_function("clearToolingMemory", &MeshOpsWrapper::clearToolingMemory)
      .class_function("fromSolid", &MeshOpsWrapper::fromSolid)
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName)
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
      .class_function("getMaxBounds", &MeshOpsWrapper::getMaxBounds)
      .class_function("getArenaMetrics", &MeshOpsWrapper::getArenaMetrics)
      .class_function("fromData", &MeshOpsWrapper::fromData)
      .function("getVertices", &MeshOpsWrapper::getVertices)
      .function("getFaces", &MeshOpsWrapper::getFaces)
      .function("getFaceCounts", &MeshOpsWrapper::getFaceCounts)
      .function("getFlatIndices", &MeshOpsWrapper::getFlatIndices)
      .function("classifyFaces", &MeshOpsWrapper::classifyFaces)
      .function("kis", &MeshOpsWrapper::kis)
      .function("ambo", &MeshOpsWrapper::ambo)
      .function("gyro", &MeshOpsWrapper::gyro)
      .function("snub", &MeshOpsWrapper::snub)
      .function("dual", &MeshOpsWrapper::dual)
      .function("truncate", &MeshOpsWrapper::truncate)
      .function("chamfer", &MeshOpsWrapper::chamfer)
      .function("bitruncate", &MeshOpsWrapper::bitruncate)
      .function("expand", &MeshOpsWrapper::expand)
      .function("hankin", &MeshOpsWrapper::hankin_rad)
      .function("meta", &MeshOpsWrapper::meta)
      .function("needle", &MeshOpsWrapper::needle)
      .function("zip", &MeshOpsWrapper::zip)
      .function("bevel", &MeshOpsWrapper::bevel)
      .function("relax", &MeshOpsWrapper::relax);

  register_vector<float>("VectorFloat");
  register_vector<int>("VectorInt");
  register_vector<uint8_t>("VectorUInt8");

  // Spline evaluation — thin wrappers returning val {x,y,z}
  function("spline_cubic_fast",
           optional_override([](float p0x, float p0y, float p0z, float p1x,
                                float p1y, float p1z, float p2x, float p2y,
                                float p2z, float p3x, float p3y, float p3z,
                                float t) -> val {
             Vector r = Spline::cubic_fast({p0x, p0y, p0z}, {p1x, p1y, p1z},
                                           {p2x, p2y, p2z}, {p3x, p3y, p3z}, t);
             val v = val::object();
             v.set("x", r.x);
             v.set("y", r.y);
             v.set("z", r.z);
             return v;
           }));
  function("spline_cubic_slerp",
           optional_override([](float p0x, float p0y, float p0z, float p1x,
                                float p1y, float p1z, float p2x, float p2y,
                                float p2z, float p3x, float p3y, float p3z,
                                float t) -> val {
             Vector r =
                 Spline::cubic_slerp({p0x, p0y, p0z}, {p1x, p1y, p1z},
                                     {p2x, p2y, p2z}, {p3x, p3y, p3z}, t);
             val v = val::object();
             v.set("x", r.x);
             v.set("y", r.y);
             v.set("z", r.z);
             return v;
           }));
}

#endif
