#!/usr/bin/env python3
"""Regenerate core/color/color_luts.h — the sRGB transfer-function lookup tables.

The device and simulator both convert between 8-bit sRGB and the engine's
16-bit linear working space through two checked-in tables (so the hot path is a
single load, never a powf). This script is their generator of record; it mirrors
the reference implementations in core/color/color.h:

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
  python scripts/generate_luts.py > core/color/color_luts.h

The generator self-formats: it pipes its output through clang-format (using the
repo .clang-format) so the result is already in committed style — no separate
manual format step to forget. If clang-format is not on PATH (set CLANG_FORMAT
to override), it emits UNFORMATTED output and prints a loud warning to stderr,
so the missing step can never pass silently. The numeric LUT data is identical
either way (clang-format only reflows whitespace), which is all the CI token
diff checks.
"""

import os
import shutil
import subprocess
import sys
from io import StringIO
from pathlib import Path


def srgb_to_linear(s):
    return s / 12.92 if s <= 0.04045 else ((s + 0.055) / 1.055) ** 2.4


def linear_to_srgb(l):
    return l * 12.92 if l <= 0.0031308 else 1.055 * l ** (1.0 / 2.4) - 0.055


# Axis sizes / quantization maxima. The 8-bit sRGB byte axis has 256 levels
# (max code 255); the 16-bit linear axis has 65536 levels (max code 65535).
# Naming them keeps the range() count, the i/<max> input normalization, the
# quantize() output ceiling, and the emitted C array bound from drifting out
# of agreement — each size otherwise appears literally in four places.
SRGB_LEVELS = 256
LINEAR_LEVELS = 65536
SRGB_MAX = SRGB_LEVELS - 1
LINEAR_MAX = LINEAR_LEVELS - 1


def quantize(v, max_val):
    return int(min(max(v, 0.0), 1.0) * max_val + 0.5)


def srgb_to_linear_lut():
    return [quantize(srgb_to_linear(i / SRGB_MAX), LINEAR_MAX) for i in range(SRGB_LEVELS)]


def linear_to_srgb_lut():
    return [quantize(linear_to_srgb(i / LINEAR_MAX), SRGB_MAX) for i in range(LINEAR_LEVELS)]


def emit_array(out, decl, values, per_row):
    out.write(decl + " = {\n")
    for start in range(0, len(values), per_row):
        row = values[start:start + per_row]
        out.write("    " + ", ".join(str(v) for v in row) + ",\n")
    out.write("};\n")


def render(out, fwd, rev):
    out.write("#pragma once\n")
    out.write('#include "platform.h"\n')
    out.write("// Generated LUTs for color conversion\n")
    out.write("// Source: scripts/generate_luts.py (Precise sRGB Transfer Function)\n")
    out.write("\n")
    out.write("// sRGB (0-255) -> Linear (0-65535)\n")
    # The per_row counts (11 u16 / 15 u8) only shape the *unformatted* fallback
    # layout: when clang-format is present it repacks each row to the column
    # limit, so these values never reach committed output. They are chosen to
    # keep the no-clang-format fallback readable (~one terminal line per row).
    emit_array(out, f"inline const uint16_t srgb_to_linear_lut[{SRGB_LEVELS}] PROGMEM",
               fwd, 11)
    out.write("\n")
    out.write("// Linear (0-65535) -> sRGB (0-255)\n")
    emit_array(out, f"inline const uint8_t linear_to_srgb_lut[{LINEAR_LEVELS}] PROGMEM",
               rev, 15)


def clang_format(text):
    """Format the generated header through clang-format using the repo style.

    Returns the formatted text, or None if clang-format is unavailable. Exits
    non-zero if clang-format is present but fails. --assume-filename points at
    the real header path so clang-format locates the repo .clang-format and
    selects the C++ language.
    """
    cf = os.environ.get("CLANG_FORMAT") or shutil.which("clang-format")
    if not cf:
        return None
    header = Path(__file__).resolve().parent.parent / "core" / "color" / "color_luts.h"
    result = subprocess.run(
        [cf, "--assume-filename=" + str(header)],
        input=text, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write("generate_luts: clang-format failed:\n" + result.stderr)
        sys.exit(1)
    return result.stdout


def check(fwd, rev):
    """Self-validate the tables: monotonicity and round-trip fidelity.

    Both transfer functions are monotonic non-decreasing, so each table must be
    too; a libm change that shifted a single entry would break this. The sRGB
    round trip (byte -> linear -> byte) must return within +/-1 code. Returns the
    number of failures (0 == pass).
    """
    fails = 0
    for name, table in (("srgb_to_linear_lut", fwd), ("linear_to_srgb_lut", rev)):
        for i in range(1, len(table)):
            if table[i] < table[i - 1]:
                sys.stderr.write(
                    f"generate_luts: {name} not monotonic at {i}: "
                    f"{table[i - 1]} -> {table[i]}\n")
                fails += 1
    for s in range(SRGB_LEVELS):
        rt = rev[fwd[s]]
        if abs(rt - s) > 1:
            sys.stderr.write(
                f"generate_luts: round trip {s} -> {fwd[s]} -> {rt} "
                f"exceeds +/-1 code\n")
            fails += 1
    return fails


def main():
    fwd = srgb_to_linear_lut()
    rev = linear_to_srgb_lut()
    if "--check" in sys.argv[1:]:
        sys.exit(1 if check(fwd, rev) else 0)
    if check(fwd, rev):
        sys.stderr.write("generate_luts: self-test failed; refusing to emit\n")
        sys.exit(1)
    buf = StringIO()
    render(buf, fwd, rev)
    text = buf.getvalue()
    formatted = clang_format(text)
    if formatted is None:
        # No clang-format: emit raw and warn loudly so the layout-normalize step
        # is never skipped silently. The numeric data is correct regardless.
        sys.stderr.write(
            "generate_luts: WARNING - clang-format not found; emitting "
            "UNFORMATTED output. Run `clang-format -i core/color/color_luts.h` before "
            "committing, or set CLANG_FORMAT to the binary path.\n")
        sys.stdout.write(text)
    else:
        sys.stdout.write(formatted)


if __name__ == "__main__":
    main()
