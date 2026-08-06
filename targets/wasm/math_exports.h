/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file math_exports.h
 * @brief Free-function embind exports of the engine's color, palette and
 *        geometry math.
 *
 * The real engine math, exported so the JS tool ports (color.js,
 * palette_math.js, lissajous_math.js, mobius_transforms.js) can cross-check
 * their mirrors against it. Included only by targets/wasm/wasm.cpp.
 */
#pragma once

#include <emscripten/bind.h>
#include "core/color/color.h"
#include "core/color/palettes.h" // HS_PROCEDURAL_PALETTE_LIST — named-palette export
#include "core/engine/platform.h"
#include "core/engine/transformers.h"
#include "core/math/3dmath.h"
#include "core/math/geometry.h"
#include "targets/wasm/wasm_predicates.h" // pure, host-tested boundary predicates
#include <cmath> // std::isfinite — all_finite() gate for the free exports
#include <initializer_list> // all_finite() variadic-arg gate for free exports
#include <array>
#include <string>

using namespace emscripten;

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

/** @brief Registers the free color/palette/geometry exports with Embind. */
static void bind_math_exports() {
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

  // HSV -> sRGB integer sextant path via the engine's CRGB(CHSV) constructor.
  // Returns sRGB bytes; there is no JS mirror, so color_parity_wasm.test.js pins
  // them to golden values. The uint8_t casts wrap h/s/v mod 256 (device CHSV
  // semantics) and the pin covers out-of-range rows; recipe compilation validates
  // its HSV keys, so the two paths handle out-of-range inputs differently by
  // design.
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
