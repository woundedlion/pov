/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file palette_bindings.h
 * @brief JS-facing palette bridge: the PaletteOps class and its embind
 *        registration.
 *
 * Bakes a GenerativePalette into the sRGB LUT the daydream palette tool previews,
 * so the tool reads the engine's own perceptual color math. Included only by
 * targets/wasm/wasm.cpp.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/color/color.h"
#include "core/engine/platform.h"
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <cstddef>
#include <cstdint>

using namespace emscripten;

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

/** @brief Registers the palette bridge's class with Embind. */
static void bind_palette_ops() {
  class_<PaletteOps>("PaletteOps")
      .constructor<>()
      .function("bakeLut", &PaletteOps::bakeLut);
}
