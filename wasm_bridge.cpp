/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
 

 /*
 cmd /c "call c:\work\emsdk\emsdk_env.bat && emcc -std=c++20 -O3 wasm_bridge.cpp -I. -lembind -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s EXPORT_NAME=createHolosphereModule -o holosphere_wasm.js"
 */


#include <emscripten/bind.h>
#include "effects.h"      // Includes all effect headers
#include "platform.h"     //


using namespace emscripten;

// Define the resolution used for the Wasm engine (must match JS Daydream.W)
constexpr int WASM_WIDTH = 96; 

class HolosphereEngine {
public:
    HolosphereEngine() {
        // Initialize with a default effect
        currentEffect = std::make_unique<Test<WASM_WIDTH>>();
        
        // Allocate a buffer for pixels (3 bytes per pixel: R, G, B)
        // We use a flat vector to be friendly to JS TypedArrays
        pixelBuffer.resize(WASM_WIDTH * H * 3);
    }

    void setEffect(std::string name) {
        // Factory pattern to switch effects based on string name
        // mirroring the logic in Holosphere.ino
        if (name == "Test") currentEffect = std::make_unique<Test<WASM_WIDTH>>();
        else if (name == "Comets") currentEffect = std::make_unique<Comets<WASM_WIDTH>>();
        else if (name == "RingSpin") currentEffect = std::make_unique<RingSpin<WASM_WIDTH>>();
        // ... Add other effects from effects.h ...
        else currentEffect = std::make_unique<Test<WASM_WIDTH>>(); // Fallback
    }

    void drawFrame() {
        if (!currentEffect) return;

        // 1. Run the C++ physics/math step
        currentEffect->draw_frame();

        // 2. Extract pixels from the Effect/Canvas and put them in a flat buffer
        // Note: Effect::get_pixel is defined in effects_engine.h
        // We iterate virtual coordinates to flatten them for JS
        int idx = 0;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < WASM_WIDTH; x++) {
                // get_pixel returns a CRGB (uint8_t r, g, b)
                CRGB p = currentEffect->get_pixel(x, y); 
                pixelBuffer[idx++] = p.r;
                pixelBuffer[idx++] = p.g;
                pixelBuffer[idx++] = p.b;
            }
        }
    }

    // Returns the memory address of the pixel buffer
    // JS will use this to create a Uint8Array view directly into Wasm memory
    uintptr_t getBufferPointer() const {
        return (uintptr_t)pixelBuffer.data();
    }

    int getBufferLength() const {
        return pixelBuffer.size();
    }

private:
    std::unique_ptr<Effect> currentEffect;
    std::vector<uint8_t> pixelBuffer;
};

// Expose to JavaScript
EMSCRIPTEN_BINDINGS(my_module) {
    class_<HolosphereEngine>("HolosphereEngine")
        .constructor<>()
        .function("setEffect", &HolosphereEngine::setEffect)
        .function("drawFrame", &HolosphereEngine::drawFrame)
        .function("getBufferPointer", &HolosphereEngine::getBufferPointer)
        .function("getBufferLength", &HolosphereEngine::getBufferLength);
}