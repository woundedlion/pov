/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
 

 /*
  Manual Build:
  cmd /c "call c:\work\emsdk\emsdk_env.bat && emcc -std=c++20 -O3 wasm_bridge.cpp -I. -lembind -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s EXPORT_NAME=createHolosphereModule -o holosphere_wasm.js"

  CMake Build:
  mkdir build; cd build
  emcmake cmake ..
  cmake --build .
  cmake --install .
  */


#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "effects.h"      // Includes all effect headers
#include "platform.h"     //
#include <string>


using namespace emscripten;

// Define the resolution used for the Wasm engine (defaults)
// Define the resolution used for the Wasm engine (defaults)
// Moved to class members


template <int W, int H>
std::unique_ptr<Effect> create_effect(const std::string& name) {
    if (name == "Test") return std::make_unique<Test<W, H>>();
    if (name == "Comets") return std::make_unique<Comets<W, H>>();
    if (name == "RingSpin") return std::make_unique<RingSpin<W, H>>();
    if (name == "MobiusGrid") return std::make_unique<MobiusGrid<W, H>>();
    if (name == "IslamicStars") return std::make_unique<IslamicStars<W, H>>();
    if (name == "MindSplatter") return std::make_unique<MindSplatter<W, H>>();
    if (name == "BZReactionDiffusion") return std::make_unique<BZReactionDiffusion<W, H>>();
    if (name == "DreamBalls") return std::make_unique<DreamBalls<W, H>>();
    if (name == "Dynamo") return std::make_unique<Dynamo<W, H>>();
    if (name == "FlowField") return std::make_unique<FlowField<W, H>>();
    if (name == "GSReactionDiffusion") return std::make_unique<GSReactionDiffusion<W, H>>();
    if (name == "GnomonicStars") return std::make_unique<GnomonicStars<W, H>>();
    if (name == "HankinSolids") return std::make_unique<HankinSolids<W, H>>();
    if (name == "HopfFibration") return std::make_unique<HopfFibration<W, H>>();
    if (name == "LSystem") return std::make_unique<LSystem<W, H>>();
    if (name == "MetaballEffect") return std::make_unique<MetaballEffect<W, H>>();
    if (name == "Moire") return std::make_unique<Moire<W, H>>();
    if (name == "PetalFlow") return std::make_unique<PetalFlow<W, H>>();
    if (name == "RingShower") return std::make_unique<RingShower<W, H>>();
    if (name == "SphericalHarmonics") return std::make_unique<SphericalHarmonics<W, H>>();
    if (name == "SpinShapes") return std::make_unique<SpinShapes<W, H>>();
    if (name == "TestShapes") return std::make_unique<TestShapes<W, H>>();
    if (name == "TestSlewRate") return std::make_unique<TestSlewRate<W, H>>();
    if (name == "TestTemporal") return std::make_unique<TestTemporal<W, H>>();
    if (name == "Thrusters") return std::make_unique<Thrusters<W, H>>();
    if (name == "Voronoi") return std::make_unique<Voronoi<W, H>>();
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
        hs::log(("WASM: setResolution called with " + std::to_string(w) + "x" + std::to_string(h)).c_str());
        if (w == pixel_width && h == pixel_height) return;

        // basic validation
        if (w > MAX_W || h > MAX_H) {
             hs::log("WASM: Resolution too large!");
             return;
        }

        pixel_width = w;
        pixel_height = h;
        
        // Resize buffer (3 bytes per pixel: R, G, B)
        pixelBuffer.resize(pixel_width * pixel_height * 3);
        
        // Re-create current effect if exists
        if (currentEffect) {
            currentEffect = nullptr;
        }
    }

    void setEffect(std::string name) {
        hs::log(("WASM: setEffect called with " + name).c_str());

        // Dispatch based on resolution
        if (pixel_width == 96 && pixel_height == 20) {
            hs::log("WASM: Creating effect <96, 20>");
            currentEffect = create_effect<96, 20>(name);
        }
        else if (pixel_width == 288 && pixel_height == 144) {
             hs::log("WASM: Creating effect <288, 144>");
            currentEffect = create_effect<288, 144>(name);
        }
        else {
             hs::log(("Unsupported resolution: " + std::to_string(pixel_width) + "x" + std::to_string(pixel_height)).c_str());
        }
        if (currentEffect) {
             hs::log("WASM: Effect created successfully");
        } else {
             hs::log("WASM: Effect creation failed (or resolution mismatch)");
        }
    }

    void setDebug(bool enabled) {
        hs::log(("WASM: setDebug called with " + std::to_string(enabled)).c_str());
        if (currentEffect) {
            currentEffect->debug_visuals = enabled;
        }
    }

    void drawFrame() {
        currentEffect->draw_frame();
        currentEffect->advance_display();

        int idx = 0;
        for (int y = 0; y < pixel_height; y++) {
            for (int x = 0; x < pixel_width; x++) {
                // Safeguard boundaries
                const Pixel& p = currentEffect->get_pixel(x, y); 

                // Convert Linear (Engine) -> sRGB (Display)
                // Using operator CRGB() which uses the LUT
                CRGB srgb = (CRGB)p;
                pixelBuffer[idx++] = srgb.r;
                pixelBuffer[idx++] = srgb.g;
                pixelBuffer[idx++] = srgb.b;
            }
        }
    }

    // Returns the memory address of the pixel buffer
    uintptr_t getBufferPointer() const {
        return (uintptr_t)pixelBuffer.data();
    }

    val getPixels() {
        return val(typed_memory_view(pixelBuffer.size(), pixelBuffer.data()));
    }

    size_t getBufferLength() const {
        return pixelBuffer.size();
    }

private:
    std::unique_ptr<Effect> currentEffect;
    std::vector<uint8_t> pixelBuffer;
    int pixel_width = 0;
    int pixel_height = 0;
};

// Expose to JavaScript
EMSCRIPTEN_BINDINGS(holosphere_engine) {
    class_<HolosphereEngine>("HolosphereEngine")
        .constructor<>()
        .function("setResolution", &HolosphereEngine::setResolution)
        .function("setEffect", &HolosphereEngine::setEffect)
        .function("setDebug", &HolosphereEngine::setDebug)
        .function("drawFrame", &HolosphereEngine::drawFrame)
        .function("getPixels", &HolosphereEngine::getPixels)
        .function("getBufferLength", &HolosphereEngine::getBufferLength);
}

#endif
