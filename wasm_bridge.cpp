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
  MeshOpsWrapper(const PolyMesh &m) : mesh(m) {}

  // Factory
  static MeshOpsWrapper fromSolid(int index) {
    return MeshOpsWrapper(Solids::get(index));
  }

  static MeshOpsWrapper fromSolidName(std::string name) {
    return MeshOpsWrapper(Solids::get_by_name(name));
  }

  static MeshOpsWrapper fromData(val vertices, val faces) {
    PolyMesh m;

    // Vertices: Float32Array [x, y, z, ...]
    std::vector<float> vData = convertJSArrayToNumberVector<float>(vertices);
    for (size_t i = 0; i < vData.size(); i += 3) {
      m.vertices.emplace_back(vData[i], vData[i + 1], vData[i + 2]);
    }

    // Faces: Flat Int32Array with -1 delimiter
    std::vector<int> fData = convertJSArrayToNumberVector<int>(faces);
    std::vector<int> currentFace;
    for (int idx : fData) {
      if (idx == -1) {
        if (!currentFace.empty())
          m.faces.push_back(currentFace);
        currentFace.clear();
      } else {
        currentFace.push_back(idx);
      }
    }
    if (!currentFace.empty())
      m.faces.push_back(currentFace);

    return MeshOpsWrapper(m);
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
    val faces = val::array();
    for (size_t i = 0; i < mesh.faces.size(); ++i) {
      val face = val::array();
      for (int idx : mesh.faces[i]) {
        face.call<void>("push", idx);
      }
      faces.set(i, face);
    }
    return faces;
  }

  val getFaceCounts() const {
    std::vector<uint8_t> counts;
    counts.reserve(mesh.faces.size());
    for (const auto &f : mesh.faces)
      counts.push_back((uint8_t)f.size());
    return val::global("Uint8Array")
        .new_(val(typed_memory_view(counts.size(), counts.data())));
  }

  val getFlatIndices() const {
    std::vector<int32_t> indices;
    for (const auto &f : mesh.faces) {
      indices.insert(indices.end(), f.begin(), f.end());
    }
    return val::global("Int32Array")
        .new_(val(typed_memory_view(indices.size(), indices.data())));
  }

  val classifyFaces() const {
    std::vector<int> colors = MeshOps::classify_faces_by_topology(mesh);
    return val::global("Int32Array")
        .new_(val(typed_memory_view(colors.size(), colors.data())));
  }

  // Operations
  MeshOpsWrapper kis() const { return MeshOpsWrapper(MeshOps::kis(mesh)); }
  MeshOpsWrapper ambo() const { return MeshOpsWrapper(MeshOps::ambo(mesh)); }
  MeshOpsWrapper gyro() const { return MeshOpsWrapper(MeshOps::gyro(mesh)); }
  MeshOpsWrapper snub() const { return MeshOpsWrapper(MeshOps::snub(mesh)); }
  MeshOpsWrapper dual() const { return MeshOpsWrapper(MeshOps::dual(mesh)); }
  MeshOpsWrapper truncate(float t) const {
    return MeshOpsWrapper(MeshOps::truncate(mesh, t));
  }
  MeshOpsWrapper expand(float t) const {
    return MeshOpsWrapper(MeshOps::expand(mesh, t));
  }
  MeshOpsWrapper hankin(float angle) const {
    return MeshOpsWrapper(MeshOps::hankin(mesh, angle * (PI_F / 180.0f)));
  }
  MeshOpsWrapper hankin_rad(float radians) const {
    return MeshOpsWrapper(MeshOps::hankin(mesh, radians));
  }
  MeshOpsWrapper bitruncate(float t) const {
    return MeshOpsWrapper(MeshOps::bitruncate(mesh, t));
  }
  MeshOpsWrapper canonicalize(int iterations) const {
    return MeshOpsWrapper(MeshOps::canonicalize(mesh, iterations));
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
      .class_function("fromSolid", &MeshOpsWrapper::fromSolid)
      .class_function("fromSolidName", &MeshOpsWrapper::fromSolidName)
      .class_function("getRegistry", &MeshOpsWrapper::getRegistry)
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
