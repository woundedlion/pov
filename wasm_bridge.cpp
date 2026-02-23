/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */

/*
 Manual Build:
 cmd /c "call c:\work\emsdk\emsdk_env.bat && emcc -std=c++20 -O3 wasm_bridge.cpp
 -I. -lembind -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s
 EXPORT_NAME=createHolosphereModule -o holosphere_wasm.js"

 CMake Build:
 mkdir build; cd build
 emcmake cmake ..
 cmake --build .
 cmake --install .
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "effects.h"  // Includes all effect headers
#include "platform.h" //
#include <string>

using namespace emscripten;

// Define Global Arenas for Webassembly
#ifdef __EMSCRIPTEN__
// 64MB for Geometry, 128MB for Scratch (WASM limit)
Arena geometry_arena(64 * 1024 * 1024);
Arena scratch_arena(64 * 1024 * 1024);
Arena scratch_arena_b(64 * 1024 * 1024);
#else
// Native Simulator gets 256MB each
Arena geometry_arena(256 * 1024 * 1024);
Arena scratch_arena(128 * 1024 * 1024);
Arena scratch_arena_b(128 * 1024 * 1024);
#endif

// Define the resolution used for the Wasm engine (defaults)
// Define the resolution used for the Wasm engine (defaults)
// Moved to class members

template <int W, int H>
std::unique_ptr<Effect> create_effect(const std::string &name) {
  static const std::map<std::string, std::function<std::unique_ptr<Effect>()>>
      factory = {
          {"Test", []() { return std::make_unique<Test<W, H>>(); }},
          {"Comets", []() { return std::make_unique<Comets<W, H>>(); }},
          {"RingSpin", []() { return std::make_unique<RingSpin<W, H>>(); }},
          {"MobiusGrid", []() { return std::make_unique<MobiusGrid<W, H>>(); }},
          {"IslamicStars",
           []() { return std::make_unique<IslamicStars<W, H>>(); }},
          {"MindSplatter",
           []() { return std::make_unique<MindSplatter<W, H>>(); }},
          {"BZReactionDiffusion",
           []() { return std::make_unique<BZReactionDiffusion<W, H>>(); }},
          {"DreamBalls", []() { return std::make_unique<DreamBalls<W, H>>(); }},
          {"Dynamo", []() { return std::make_unique<Dynamo<W, H>>(); }},
          {"FlowField", []() { return std::make_unique<FlowField<W, H>>(); }},
          {"GSReactionDiffusion",
           []() { return std::make_unique<GSReactionDiffusion<W, H>>(); }},
          {"GnomonicStars",
           []() { return std::make_unique<GnomonicStars<W, H>>(); }},
          {"HankinSolids",
           []() { return std::make_unique<HankinSolids<W, H>>(); }},
          {"HopfFibration",
           []() { return std::make_unique<HopfFibration<W, H>>(); }},
          {"LSystem", []() { return std::make_unique<LSystem<W, H>>(); }},
          {"MetaballEffect",
           []() { return std::make_unique<MetaballEffect<W, H>>(); }},
          {"Moire", []() { return std::make_unique<Moire<W, H>>(); }},
          {"PetalFlow", []() { return std::make_unique<PetalFlow<W, H>>(); }},
          {"RingShower", []() { return std::make_unique<RingShower<W, H>>(); }},
          {"SphericalHarmonics",
           []() { return std::make_unique<SphericalHarmonics<W, H>>(); }},
          {"SpinShapes", []() { return std::make_unique<SpinShapes<W, H>>(); }},
          {"TestShapes", []() { return std::make_unique<TestShapes<W, H>>(); }},
          {"TestSlewRate",
           []() { return std::make_unique<TestSlewRate<W, H>>(); }},
          {"TestTemporal",
           []() { return std::make_unique<TestTemporal<W, H>>(); }},
          {"Thrusters", []() { return std::make_unique<Thrusters<W, H>>(); }},
          {"Voronoi", []() { return std::make_unique<Voronoi<W, H>>(); }},
          {"FlamingMesh",
           []() { return std::make_unique<FlamingMesh<W, H>>(); }},
          {"ChaoticStrings",
           []() { return std::make_unique<ChaoticStrings<W, H>>(); }}};

  auto it = factory.find(name);
  if (it != factory.end()) {
    return it->second();
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
    hs::log(("WASM: setResolution called with " + std::to_string(w) + "x" +
             std::to_string(h))
                .c_str());
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
    hs::log(("WASM: setEffect called with " + name).c_str());

    // Ensure strictly sequential destruction and creation to minimize peak
    // memory
    currentEffect.reset();

    // Use factory
    if (pixel_width == 96 && pixel_height == 20)
      currentEffect = create_effect<96, 20>(name);
    else if (pixel_width == 288 && pixel_height == 144)
      currentEffect = create_effect<288, 144>(name);
    else if (pixel_width == 576 && pixel_height == 288)
      currentEffect = create_effect<576, 288>(name);
    else {
      hs::log("WASM: Unsupported resolution for factory!");
      return;
    }
  }

  void drawFrame() {
    if (!currentEffect)
      return;

    // Draw
    currentEffect->draw_frame();
    currentEffect->advance_display();

    // Copy to buffer
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
      currentEffect->updateParameter(name, value);
    }
  }

  val getParameterDefinitions() {
    if (!currentEffect)
      return val::array();

    val result = val::array();
    const auto &params = currentEffect->getParameters();

    int i = 0;
    for (const auto &[key, def] : params) {
      val entry = val::object();
      entry.set("name", def.name);

      if (def.type == Effect::ParamType::BOOL) {
        entry.set("value", *static_cast<bool *>(def.target));
        // For boolean, we don't necessarily need min/max, but keeping them
        // doesn't hurt if 0/1 logic is used downstream But specifically for
        // dat.gui checkbox, value determines type.
      } else {
        entry.set("value", *static_cast<float *>(def.target));
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

    for (const auto &[key, def] : params) {
      if (def.type == Effect::ParamType::BOOL) {
        paramValues.push_back(*static_cast<bool *>(def.target) ? 1.0f : 0.0f);
      } else {
        paramValues.push_back(*static_cast<float *>(def.target));
      }
    }
    return val(typed_memory_view(paramValues.size(), paramValues.data()));
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

  // Factory
  static MeshOpsWrapper *fromSolid(int index) {
    return new MeshOpsWrapper(
        Solids::get(geometry_arena, scratch_arena, index));
  }

  static MeshOpsWrapper *fromSolidName(std::string name) {
    return new MeshOpsWrapper(
        Solids::get_by_name(geometry_arena, scratch_arena, name));
  }

  static MeshOpsWrapper *fromData(val vertices, val faces) {
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

    m.faces.initialize(geometry_arena,
                       fData.size() -
                           num_faces); // Exact size without delimiters
    m.face_counts.initialize(geometry_arena, num_faces);

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

    return new MeshOpsWrapper(std::move(m));
  }

  // Accessors for JS
  val getVertices() const {
    std::vector<float> data;
    data.reserve(mesh.vertices.size() * 3);
    for (const auto &v : mesh.vertices) {
      data.push_back(v.i);
      data.push_back(v.j);
      data.push_back(v.k);
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

  val classifyFaces() const {
    std::vector<int> colors =
        MeshOps::classify_faces_by_topology(mesh, scratch_arena);
    return val::global("Int32Array")
        .new_(val(typed_memory_view(colors.size(), colors.data())));
  }

  // Operations
  MeshOpsWrapper *kis() const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::kis(mesh, ctx));
  }
  MeshOpsWrapper *ambo() const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::ambo(mesh, ctx));
  }
  MeshOpsWrapper *gyro() const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::gyro(mesh, ctx));
  }
  MeshOpsWrapper *snub() const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::snub(mesh, ctx));
  }
  MeshOpsWrapper *dual() const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::dual(mesh, ctx));
  }
  MeshOpsWrapper *truncate(float t) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::truncate(mesh, ctx, t));
  }
  MeshOpsWrapper *expand(float t) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::expand(mesh, ctx, t));
  }
  MeshOpsWrapper *hankin(float angle) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(
        MeshOps::hankin(mesh, ctx, angle * (PI_F / 180.0f)));
  }
  MeshOpsWrapper *hankin_rad(float radians) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::hankin(mesh, ctx, radians));
  }
  MeshOpsWrapper *bitruncate(float t) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::bitruncate(mesh, ctx, t));
  }
  MeshOpsWrapper *canonicalize(int iterations) const {
    ScratchContext ctx(scratch_arena, scratch_arena_b);
    return new MeshOpsWrapper(MeshOps::canonicalize(mesh, ctx, iterations));
  }
  static val getRegistry() {
    val registry = val::array();
    for (int i = 0; i < Solids::NUM_ENTRIES; ++i) {
      const auto &entry = Solids::get_entry(i);
      val item = val::object();
      item.set("name", std::string(entry.name));
      item.set("category", entry.category == Solids::Category::Simple
                               ? "Simple"
                               : "Complex");
      registry.set(i, item);
    }
    return registry;
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
      .function("getParamValues", &HolosphereEngine::getParamValues);

  class_<MeshOpsWrapper>("MeshOps")
      .constructor<>()
      .class_function("fromSolid", &MeshOpsWrapper::fromSolid,
                      allow_raw_pointers())
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName,
                      allow_raw_pointers())
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
      .class_function("fromData", &MeshOpsWrapper::fromData,
                      allow_raw_pointers())
      .function("getVertices", &MeshOpsWrapper::getVertices)
      .function("getFaces", &MeshOpsWrapper::getFaces)
      .function("getFaceCounts", &MeshOpsWrapper::getFaceCounts)
      .function("getFlatIndices", &MeshOpsWrapper::getFlatIndices)
      .function("classifyFaces", &MeshOpsWrapper::classifyFaces)
      .function("kis", &MeshOpsWrapper::kis, allow_raw_pointers())
      .function("ambo", &MeshOpsWrapper::ambo, allow_raw_pointers())
      .function("gyro", &MeshOpsWrapper::gyro, allow_raw_pointers())
      .function("snub", &MeshOpsWrapper::snub, allow_raw_pointers())
      .function("dual", &MeshOpsWrapper::dual, allow_raw_pointers())
      .function("truncate", &MeshOpsWrapper::truncate, allow_raw_pointers())
      .function("bitruncate", &MeshOpsWrapper::bitruncate, allow_raw_pointers())
      .function("expand", &MeshOpsWrapper::expand, allow_raw_pointers())
      .function("hankin", &MeshOpsWrapper::hankin_rad, allow_raw_pointers())
      .function("canonicalize", &MeshOpsWrapper::canonicalize,
                allow_raw_pointers());

  register_vector<float>("VectorFloat");
  register_vector<int>("VectorInt");
  register_vector<uint8_t>("VectorUInt8");
}

#endif
