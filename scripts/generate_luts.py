#!/usr/bin/env python3
"""Regenerate core/color_luts.h — the sRGB transfer-function lookup tables.

The device and simulator both convert between 8-bit sRGB and the engine's
16-bit linear working space through two checked-in tables (so the hot path is a
single load, never a powf). This script is their generator of record; it mirrors
the reference implementations in core/color.h:

  srgb_to_linear(s) = s/12.92                       if s <= 0.04045
                      ((s+0.055)/1.055) ** 2.4       otherwise
  linear_to_srgb(l) = l*12.92                        if l <= 0.0031308
                      1.055 * l**(1/2.4) - 0.055      otherwise

quantized by round(clamp(v,0,1) * max), i.e. the +0.5 rule the C++
float_to_pixel16 / CRGB byte conversions use:

  srgb_to_linear_lut[256]    : sRGB byte i  -> round(srgb_to_linear(i/255)   * 65535)
  linear_to_srgb_lut[65536]  : linear u16 i -> round(linear_to_srgb(i/65535) * 255)

Precision note: the checked-in tables were generated in IEEE double precision and
reproduce bit-for-bit only in doubles — Python's native float (and numpy float64)
is correct; computing in float32 flips one srgb_to_linear entry at an LSB
rounding boundary.

Usage:
  python scripts/generate_luts.py > core/color_luts.h
  clang-format -i core/color_luts.h   # normalize array layout to repo style
"""

import sys


def srgb_to_linear(s):
    return s / 12.92 if s <= 0.04045 else ((s + 0.055) / 1.055) ** 2.4


def linear_to_srgb(l):
    return l * 12.92 if l <= 0.0031308 else 1.055 * l ** (1.0 / 2.4) - 0.055


def quantize(v, max_val):
    return int(min(max(v, 0.0), 1.0) * max_val + 0.5)


def srgb_to_linear_lut():
    return [quantize(srgb_to_linear(i / 255.0), 65535) for i in range(256)]


def linear_to_srgb_lut():
    return [quantize(linear_to_srgb(i / 65535.0), 255) for i in range(65536)]


def emit_array(out, decl, values, per_row):
    out.write(decl + " = {\n")
    for start in range(0, len(values), per_row):
        row = values[start:start + per_row]
        out.write("    " + ", ".join(str(v) for v in row) + ",\n")
    out.write("};\n")


def main():
    out = sys.stdout
    out.write("#pragma once\n")
    out.write('#include "platform.h"\n')
    out.write("// Generated LUTs for color conversion\n")
    out.write("// Source: scripts/generate_luts.py (Precise sRGB Transfer Function)\n")
    out.write("\n")
    out.write("// sRGB (0-255) -> Linear (0-65535)\n")
    emit_array(out, "inline const uint16_t srgb_to_linear_lut[256] PROGMEM",
               srgb_to_linear_lut(), 11)
    out.write("\n")
    out.write("// Linear (0-65535) -> sRGB (0-255)\n")
    emit_array(out, "inline const uint8_t linear_to_srgb_lut[65536] PROGMEM",
               linear_to_srgb_lut(), 15)


if __name__ == "__main__":
    main()
