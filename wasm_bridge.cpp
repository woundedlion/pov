/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "effects.h"  // Includes all effect headers
#include "platform.h" //
#include <string_view>

using namespace emscripten;

// Dedicated arena for JavaScript Tools
Arena tooling_arena(8 * 1024 * 1024);
Arena tooling_scratch_a(4 * 1024 * 1024);
Arena tooling_scratch_b(4 * 1024 * 1024);

template <int W, int H>
std::unique_ptr<Effect> create_effect(std::string_view name) {
  struct FactoryEntry {
    std::string_view name;
    std::function<std::unique_ptr<Effect>()> creator;
  };
  static const FactoryEntry factory[] = {
      {"Test", []() { return std::make_unique<Test<W, H>>(); }},
      {"Comets", []() { return std::make_unique<Comets<W, H>>(); }},
      {"RingSpin", []() { return std::make_unique<RingSpin<W, H>>(); }},
      {"MobiusGrid", []() { return std::make_unique<MobiusGrid<W, H>>(); }},
      {"IslamicStars", []() { return std::make_unique<IslamicStars<W, H>>(); }},
      {"MindSplatter", []() { return std::make_unique<MindSplatter<W, H>>(); }},
      {"BZReactionDiffusion",
       []() { return std::make_unique<BZReactionDiffusion<W, H>>(); }},
      {"DreamBalls", []() { return std::make_unique<DreamBalls<W, H>>(); }},
      {"Dynamo", []() { return std::make_unique<Dynamo<W, H>>(); }},
      {"FlowField", []() { return std::make_unique<FlowField<W, H>>(); }},
      {"GSReactionDiffusion",
       []() { return std::make_unique<GSReactionDiffusion<W, H>>(); }},
      {"GnomonicStars",
       []() { return std::make_unique<GnomonicStars<W, H>>(); }},
      {"HankinSolids", []() { return std::make_unique<HankinSolids<W, H>>(); }},
      {"HopfFibration",
       []() { return std::make_unique<HopfFibration<W, H>>(); }},
      {"LSystem", []() { return std::make_unique<LSystem<W, H>>(); }},
      {"Metaballs", []() { return std::make_unique<Metaballs<W, H>>(); }},
      {"Moire", []() { return std::make_unique<Moire<W, H>>(); }},
      {"PetalFlow", []() { return std::make_unique<PetalFlow<W, H>>(); }},
      {"RingShower", []() { return std::make_unique<RingShower<W, H>>(); }},
      {"SphericalHarmonics",
       []() { return std::make_unique<SphericalHarmonics<W, H>>(); }},
      {"SpinShapes", []() { return std::make_unique<SpinShapes<W, H>>(); }},
      {"TestShapes", []() { return std::make_unique<TestShapes<W, H>>(); }},
      {"TestSlewRate", []() { return std::make_unique<TestSlewRate<W, H>>(); }},
      {"TestTemporal", []() { return std::make_unique<TestTemporal<W, H>>(); }},
      {"Thrusters", []() { return std::make_unique<Thrusters<W, H>>(); }},
      {"Voronoi", []() { return std::make_unique<Voronoi<W, H>>(); }},
      {"FlamingMesh", []() { return std::make_unique<FlamingMesh<W, H>>(); }},
      {"Liquid2D", []() { return std::make_unique<Liquid2D<W, H>>(); }},
      {"ChaoticStrings",
       []() { return std::make_unique<ChaoticStrings<W, H>>(); }}};

  for (const auto &entry : factory) {
    if (name == entry.name) {
      return entry.creator();
    }
  }
  return std::make_unique<Test<W, H>>(); // Fallback
}

class HolosphereEngine {
public:
  HolosphereEngine() {
    // Seed randomness
    srand(static_cast<unsigned int>(time(NULL)));

    // Initialize with default
    setResolution(96, 20);
    setEffect("Test");
  }

  void setResolution(int w, int h) {
    if (w == pixel_width && h == pixel_height)
      return;

    // basic validation
    if (w > MAX_W || h > MAX_H) {
      hs::log("WASM: Resolution too large!");
      return;
    }

    pixel_width = w;
    pixel_height = h;

    // Resize buffer (16-bit linear: uint16_t per channel)
    pixelBuffer.resize(pixel_width * pixel_height * 3);

    // Re-create current effect if exists
    if (currentEffect) {
      currentEffect = nullptr;
    }
  }

  void setEffect(std::string name) {
    char buf[128];
    snprintf(buf, sizeof(buf), "WASM: setEffect called with %s", name.c_str());
    hs::log(buf);

    currentEffect.reset();
    scratch_arena_a.reset();
    scratch_arena_a.reset_high_water_mark();
    scratch_arena_b.reset();
    scratch_arena_b.reset_high_water_mark();
    persistent_arena.reset();
    persistent_arena.reset_high_water_mark();

    PersistentTracker::clear_registry();

    if (pixel_width == 96 && pixel_height == 20)
      currentEffect = create_effect<96, 20>(name);
    else if (pixel_width == 288 && pixel_height == 144)
      currentEffect = create_effect<288, 144>(name);
    else {
      hs::log("WASM: Unsupported resolution for factory!");
      return;
    }
  }

  void drawFrame() {
    if (!currentEffect)
      return;

    currentEffect->draw_frame();
    currentEffect->advance_display();

    // Output 16-bit Linear values directly
    int idx = 0;
    for (int y = 0; y < pixel_height; y++) {
      for (int x = 0; x < pixel_width; x++) {
        // Get the pixel (Pixel16)
        const Pixel &p = currentEffect->get_pixel(x, y);

        // Direct 16-bit copy (Linear)
        pixelBuffer[idx++] = p.r;
        pixelBuffer[idx++] = p.g;
        pixelBuffer[idx++] = p.b;
      }
    }
  }

  // Expose the raw memory view to JS to avoid copying overhead
  val getPixels() {
    // Return Uint16Array view
    return val(typed_memory_view(pixelBuffer.size(), pixelBuffer.data()));
  }

  int getBufferLength() { return pixelBuffer.size(); }

  void setParameter(std::string name, float value) {
    if (currentEffect) {
      currentEffect->updateParameter(name.c_str(), value);
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
      entry.set("value", def.get());

      if (!def.is_bool()) {
        entry.set("min", def.min);
        entry.set("max", def.max);
      }
      result.set(i++, entry);
    }
    return result;
  }

  val getParamValues() {
    if (!currentEffect)
      return val::array();

    const auto &params = currentEffect->getParameters();
    paramValues.clear();
    paramValues.reserve(params.size());

    for (const auto &def : params) {
      paramValues.push_back(def.get());
    }
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

    return metrics;
  }

private:
  std::unique_ptr<Effect> currentEffect;
  std::vector<uint16_t> pixelBuffer; // 16-bit
  std::vector<float> paramValues;    // Backing store for getParamValues
  int pixel_width = 0;
  int pixel_height = 0;
};

// ==========================================================================================
// 3. MESH OPS BINDINGS
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
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::get(tooling_arena, ctx, index));
  }

  static std::unique_ptr<MeshOpsWrapper> fromSolidName(std::string name) {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::get_by_name(tooling_arena, ctx, name));
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

    m.faces.initialize(tooling_arena,
                       fData.size() -
                           num_faces); // Exact size without delimiters
    m.face_counts.initialize(tooling_arena, num_faces);

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
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    MeshOps::classify_faces_by_topology(mesh, ctx);
    return val::global("Int32Array")
        .new_(
            val(typed_memory_view(mesh.topology.size(), mesh.topology.data())));
  }

  // Operations
  std::unique_ptr<MeshOpsWrapper> kis() const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::kis(mesh, ctx), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> ambo() const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::ambo(mesh, ctx), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> gyro() const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::gyro(mesh, ctx), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> snub() const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::snub(mesh, ctx), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> dual() const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::dual(mesh, ctx), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> truncate(float t) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::truncate(mesh, ctx, t), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> expand(float t) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(
        Solids::finalize_solid(MeshOps::expand(mesh, ctx, t), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> hankin(float angle) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        MeshOps::hankin(mesh, ctx, angle * (PI_F / 180.0f)), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> hankin_rad(float radians) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        MeshOps::hankin(mesh, ctx, radians), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> bitruncate(float t) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        MeshOps::bitruncate(mesh, ctx, t), tooling_arena));
  }
  std::unique_ptr<MeshOpsWrapper> canonicalize(int iterations) const {
    MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
    return std::make_unique<MeshOpsWrapper>(Solids::finalize_solid(
        MeshOps::canonicalize(mesh, ctx, iterations), tooling_arena));
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
      MemoryCtx ctx(tooling_scratch_a, tooling_scratch_b);
      PolyMesh temp = Solids::get(tooling_arena, ctx, i);

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
      .function("getParameterDefinitions",
                &HolosphereEngine::getParameterDefinitions)
      .function("getParamValues", &HolosphereEngine::getParamValues)
      .function("getArenaMetrics", &HolosphereEngine::getArenaMetrics);

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
      .function("bitruncate", &MeshOpsWrapper::bitruncate)
      .function("expand", &MeshOpsWrapper::expand)
      .function("hankin", &MeshOpsWrapper::hankin_rad)
      .function("canonicalize", &MeshOpsWrapper::canonicalize);

  register_vector<float>("VectorFloat");
  register_vector<int>("VectorInt");
  register_vector<uint8_t>("VectorUInt8");
}

#endif
