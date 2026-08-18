#!/usr/bin/env python3
"""Generates core/color/gamut_lut.h: the sRGB gamut boundary bracket table used
by gamut_clip_preserve_chroma.

Preserve-chroma clipping holds L and hue fixed and scales (a, b) uniformly, so
the clipped chroma is exactly min(C, C_max(hue, L)). C_max is a static property
of the sRGB gamut, so it is bracketed here instead of solved per pixel.

The table is indexed by the diamond angle of (b, a) -- the trig-free angle from
diamond_angle() in core/math/3dmath.h -- and by L. A cell stores the MINIMUM and
the MAXIMUM C_max over the region it covers, so the true C_max of any ray in the
cell lies inside the stored bracket. The per-pixel path bisects that bracket
against the channel cubics, so residual error is the bracket width halved once
per step rather than the whole width.

A cell minimum alone is not enough: the gamut cusp jumps to a different RGB-cube
edge as hue crosses a cube vertex, a feature narrower than any affordable cell,
so the minimum under-saturates by up to 0.041 chroma at every resolution tried.
The bracket keeps that feature inside [min, max] and the bisection resolves it.

BOUNDARY DEFINITION. C_max is the FIRST EXIT: the smallest C > 0 that leaves the
gamut. That is not the same as "the largest in-gamut C" -- linear_rgb_in_gamut's
+-1e-4 slack lets a channel graze zero, leave the tolerance and re-enter, so the
in-gamut set along a ray is occasionally disconnected (~0.1% of sampled rays;
the widest gap observed spans C 0.251 to 0.286 at L 0.419). A plain bisection
lands on either side of such a gap depending on where its midpoints fall, which
makes it discontinuous in L by up to 0.038. First exit is the boundary of the
connected component containing C = 0, is continuous, and agrees with the
bisection everywhere the set is connected.

Usage: python tools/gen_gamut_lut.py [output_path]
       python tools/gen_gamut_lut.py --check

--check regenerates the table in memory and diffs it against the committed
header in full, pins the constants mirrored below against core/color/color.h,
pins the mirrored diamond angle against core/math/3dmath.h, and round-trips the
angle parameterization. Wired as ctest unit_gamut_lut and CI
gamut-lut-provenance.
"""

import argparse
import difflib
import itertools
import os
import re
import sys

import numpy as np

# Flash master resolution. init_gamut_lut() downsamples by integer factors, so
# these bound the finest grid any effect can request.
#
# Resolution only sets how wide the bracket the per-pixel path refines starts,
# so it buys accuracy, not reach. Worst first-exit deficit over the color
# suite's sweep, against that suite's 5e-3 bound:
#   512x256 0.00132 | 256x128 0.00139 | 128x64 0.00235
#   128x32  0.00291 | 64x64 0.00347 | 64x32 0.00360 | 32x16 0.00530
# The deficit is not what binds. From 128x32 down the bracket is wide enough
# that the four-step walk strides over a disconnected in-gamut interval and
# lands past the first exit, oversaturating by up to 0.03 chroma. This grid
# keeps two steps of margin on that, at a quarter of 512x256's flash.
ANGLE_STEPS = 256
L_STEPS = 128
# 65535 / 0.5: OKLab chroma inside sRGB stays below 0.5, so this spends the full
# uint16 range on the live domain at ~7.6e-6 resolution.
SCALE = 131070.0
# Sub-samples per cell per axis, closed at both ends. Set so the sub-sample
# spacing, not the cell size, bounds how far a cell extreme can be missed.
SUBSAMPLES = 16
# Absolute chroma slack widening the bracket before quantization, covering the
# residual between the sub-sampled cell extremes and the true continuous ones.
GUARD = 1e-4
# First-exit solver: coarse scan to bracket, then bisection inside the bracket.
C_HI = 0.45
COARSE = 192
BISECT_ITERS = 28
# Points checked below a candidate to prove the in-gamut set is connected up to
# it. A ray that fails re-solves with the full coarse scan.
CONNECT_CHECKS = 24

# Matches core/color/color.h linear_rgb_in_gamut(). --check pins these against
# that header.
GAMUT_EPS = 1e-4
GAMUT_LO = -GAMUT_EPS
GAMUT_HI = 1.0 + GAMUT_EPS

# Mirrors color.h oklab_to_lms_cbrt(): (a, b) coefficients per cone row.
OKLAB_TO_LMS = ((0.3963377774, 0.2158037573),
                (-0.1055613458, -0.0638541728),
                (-0.0894841775, -1.2914855480))
# Mirrors color.h lms_cbrt_to_linear_rgb(): (l, m, s) coefficients per channel.
LMS_TO_RGB = ((+4.0767416621, -3.3077115913, +0.2309699292),
              (-1.2684380046, +2.6097574011, -0.3413193965),
              (-0.0041960863, -0.7034186147, +1.7076147010))


def oklab_to_lms_cbrt(L, a, b):
    """Inverse OKLab matrix; mirrors color.h oklab_to_lms_cbrt()."""
    l_ = L + OKLAB_TO_LMS[0][0] * a + OKLAB_TO_LMS[0][1] * b
    m_ = L + OKLAB_TO_LMS[1][0] * a + OKLAB_TO_LMS[1][1] * b
    s_ = L + OKLAB_TO_LMS[2][0] * a + OKLAB_TO_LMS[2][1] * b
    return l_, m_, s_


def lms_cbrt_to_linear_rgb(l_, m_, s_):
    """Cube plus RGB matrix; mirrors color.h lms_cbrt_to_linear_rgb()."""
    l = l_ * l_ * l_
    m = m_ * m_ * m_
    s = s_ * s_ * s_
    r = LMS_TO_RGB[0][0] * l + LMS_TO_RGB[0][1] * m + LMS_TO_RGB[0][2] * s
    g = LMS_TO_RGB[1][0] * l + LMS_TO_RGB[1][1] * m + LMS_TO_RGB[1][2] * s
    b = LMS_TO_RGB[2][0] * l + LMS_TO_RGB[2][1] * m + LMS_TO_RGB[2][2] * s
    return r, g, b


def in_gamut(L, a, b):
    """Elementwise linear_rgb_in_gamut() of the OKLab triple (L, a, b)."""
    r, g, bl = lms_cbrt_to_linear_rgb(*oklab_to_lms_cbrt(L, a, b))
    inside = (r >= GAMUT_LO) & (r <= GAMUT_HI)
    inside &= (g >= GAMUT_LO) & (g <= GAMUT_HI)
    inside &= (bl >= GAMUT_LO) & (bl <= GAMUT_HI)
    return inside


def diamond_direction(t):
    """Unit (a, b) direction whose diamond angle is t.

    Inverts diamond_angle(y=b, x=a): the diamond |x| + |y| = 1 is walked
    counter-clockwise from +x, one unit of t per quadrant.
    """
    q = np.floor(t).astype(np.int64) % 4
    f = t - np.floor(t)
    xs = np.array([1.0, 0.0, -1.0, 0.0])
    ys = np.array([0.0, 1.0, 0.0, -1.0])
    x = xs[q] * (1.0 - f) + xs[(q + 1) % 4] * f
    y = ys[q] * (1.0 - f) + ys[(q + 1) % 4] * f
    n = np.hypot(x, y)
    return x / n, y / n


def diamond_angle(y, x):
    """Diamond pseudo-angle of (x, y) in [0, 4); mirrors 3dmath.h
    diamond_angle(), the parameterization the runtime cell lookup indexes by."""
    d = np.abs(x) + np.abs(y)
    r = y / np.where(d < 1e-20, 1.0, d)
    a = np.where(y >= 0.0,
                 np.where(x >= 0.0, r, 2.0 - r),
                 np.where(x < 0.0, 2.0 - r, 4.0 + r))
    return np.where((d < 1e-20) | (a >= 4.0), 0.0, a)


def _bisect(L, a_dir, b_dir, lo, hi, iters=BISECT_ITERS):
    """Refines a bracket [lo, hi] straddling a gamut crossing."""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        inside = in_gamut(L, a_dir * mid, b_dir * mid)
        lo = np.where(inside, mid, lo)
        hi = np.where(inside, hi, mid)
    return lo


def _scan_first_exit(L, a_dir, b_dir):
    """First exit by a full coarse scan; correct but pays COARSE evaluations."""
    Cs = np.linspace(0.0, C_HI, COARSE + 1)
    out = np.zeros(L.shape + (COARSE + 1,), dtype=bool)
    for k in range(COARSE + 1):
        out[..., k] = ~in_gamut(L, a_dir * Cs[k], b_dir * Cs[k])
    any_out = out.any(axis=-1)
    k = np.where(any_out, np.argmax(out, axis=-1), COARSE)
    lo = Cs[np.maximum(k - 1, 0)]
    hi = Cs[k]
    return np.where(any_out, _bisect(L, a_dir, b_dir, lo, hi), C_HI)


def c_max(L, a_dir, b_dir):
    """First-exit chroma along each ray.

    Fast path: plain bisection, then verify the in-gamut set is connected below
    the answer. The rare ray that fails re-solves with the full coarse scan.
    """
    cand = _bisect(L, a_dir, b_dir, np.zeros_like(L), np.full_like(L, C_HI))

    connected = np.ones(L.shape, dtype=bool)
    for i in range(1, CONNECT_CHECKS + 1):
        t = cand * (i / (CONNECT_CHECKS + 1.0))
        connected &= in_gamut(L, a_dir * t, b_dir * t)

    if not connected.all():
        idx = ~connected
        cand[idx] = _scan_first_exit(L[idx], a_dir[idx], b_dir[idx])
    return cand


def window_min(arr, axis, window, stride):
    """Minimum over overlapping closed cells: cell i spans [i*stride,
    i*stride+window) along `axis`."""
    n = (arr.shape[axis] - window) // stride + 1
    idx = np.arange(n)[:, None] * stride + np.arange(window)[None, :]
    return np.take(arr, idx, axis=axis).min(axis=axis + 1)


def build_table():
    """Returns the (L_STEPS, ANGLE_STEPS, 2) uint16 bracket table plus the worst
    bracket width, which bounds the pre-refinement error."""
    n_l = L_STEPS * SUBSAMPLES + 1
    l_samples = np.arange(n_l) / (L_STEPS * SUBSAMPLES)

    table = np.zeros((L_STEPS, ANGLE_STEPS, 2), dtype=np.uint16)
    worst_width = 0.0

    # Angle cells run in groups so one chunk stays a sane array size.
    group = 32
    for a0 in range(0, ANGLE_STEPS, group):
        n_cells = min(group, ANGLE_STEPS - a0)
        n_ang = n_cells * SUBSAMPLES + 1
        base = a0 * SUBSAMPLES
        angle_samples = (base + np.arange(n_ang)) * (
            4.0 / (ANGLE_STEPS * SUBSAMPLES))
        a_dir, b_dir = diamond_direction(angle_samples)

        LL = np.ascontiguousarray(
            np.broadcast_to(l_samples[None, :], (n_ang, n_l)))
        AD = np.ascontiguousarray(np.broadcast_to(a_dir[:, None], LL.shape))
        BD = np.ascontiguousarray(np.broadcast_to(b_dir[:, None], LL.shape))
        fe = c_max(LL, AD, BD)

        cell_min = window_min(fe, 1, SUBSAMPLES + 1, SUBSAMPLES)
        cell_min = window_min(cell_min, 0, SUBSAMPLES + 1, SUBSAMPLES)
        cell_max = -window_min(-fe, 1, SUBSAMPLES + 1, SUBSAMPLES)
        cell_max = -window_min(-cell_max, 0, SUBSAMPLES + 1, SUBSAMPLES)
        worst_width = max(worst_width, float((cell_max - cell_min).max()))

        q_lo = np.floor(np.maximum(cell_min - GUARD, 0.0) * SCALE)
        q_hi = np.ceil((cell_max + GUARD) * SCALE)
        table[:, a0:a0 + n_cells, 0] = np.clip(q_lo, 0, 65535).astype(
            np.uint16).T
        table[:, a0:a0 + n_cells, 1] = np.clip(q_hi, 0, 65535).astype(
            np.uint16).T

        sys.stderr.write("  angle cells %d/%d\n" % (a0 + n_cells, ANGLE_STEPS))

    return table, worst_width


def render(table):
    """Returns the generated header text."""
    flat = table.ravel()
    lines = []
    # 11 five-digit values plus indent stays inside the 80-column convention.
    per_line = 11
    for i in range(0, flat.size, per_line):
        lines.append("    " + ", ".join("%d" % v for v in flat[i:i + per_line])
                     + ",")
    body = "\n".join(lines)

    header = f"""/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdint>

#include "platform/platform.h"

// Generated by tools/gen_gamut_lut.py. Do not edit.
//
// sRGB gamut boundary chroma C_max, indexed by the diamond angle of (b, a) and
// by L. Each cell holds the minimum and the maximum C_max over the region it
// covers, so the true C_max of any ray in the cell lies inside the stored
// bracket; the per-pixel path bisects the bracket against the channel cubics.
// C_max is the first exit from the gamut, not the largest in-gamut chroma; the
// generator explains why the two differ.
//
// This is the FLASH master at full resolution. init_gamut_lut() downsamples it
// into the arena by integer factors, taking the minimum of the merged minima
// and the maximum of the merged maxima.

/** @brief Diamond-angle buckets spanning [0, 4) in the flash master. */
inline constexpr int GAMUT_LUT_ANGLE_STEPS = {ANGLE_STEPS};
/** @brief Lightness buckets spanning [0, 1] in the flash master. */
inline constexpr int GAMUT_LUT_L_STEPS = {L_STEPS};
/** @brief Stored value per chroma unit; 65535 / 0.5 spends the full uint16
 *  range on the reachable chroma domain. */
inline constexpr float GAMUT_LUT_SCALE = {SCALE:.1f}f;
/** @brief Reciprocal of GAMUT_LUT_SCALE, for the per-pixel decode. */
inline constexpr float GAMUT_LUT_INV_SCALE = 1.0f / GAMUT_LUT_SCALE;
/** @brief Flash master entry count; two entries (min, max) per cell, L-major,
 *  so one L row is contiguous. */
inline constexpr int GAMUT_LUT_ENTRIES =
    GAMUT_LUT_ANGLE_STEPS * GAMUT_LUT_L_STEPS * 2;

/**
 * @brief Arena bytes an (angle_steps x l_steps) bracket copy occupies.
 * @param angle_steps Diamond-angle buckets requested.
 * @param l_steps Lightness buckets requested.
 * @return Byte size of that grid, for arena budget static_asserts.
 */
inline constexpr size_t gamut_lut_bytes(int angle_steps, int l_steps) {{
  return static_cast<size_t>(angle_steps) * l_steps * 2 * sizeof(uint16_t);
}}

/** @brief Flash-resident master; index [(l * ANGLE_STEPS + angle) * 2], min
 *  first then max. */
inline const uint16_t GAMUT_LUT[GAMUT_LUT_ENTRIES] HS_PROGMEM_UNIQUE(GAMUT_LUT) = {{
{body}
}};
"""
    return header


def _function_body(text, signature):
    """Source between the opening brace of `signature` and its column-0 close."""
    start = text.index(signature)
    open_brace = text.index("{", start)
    return text[open_brace:text.index("\n}", open_brace)]


NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def _split_paren(text):
    """Splits `text`, which opens with '(', into the parenthesized part and the
    remainder."""
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return text[1:i], text[i + 1:]
    raise ValueError("unbalanced parentheses in %r" % text)


def _cpp_expr_to_python(expr, names):
    """Rewrites one C++ float expression into the equivalent Python source.

    Covers the grammar the mirrored helpers are written in: arithmetic over the
    named locals, float literal suffixes, std::fabs and nested ternaries. Any
    other call or identifier raises ValueError rather than being guessed at.
    """
    expr = expr.strip()
    depth = 0
    q = -1
    for i, ch in enumerate(expr):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "?" and depth == 0:
            q = i
            break
    if q >= 0:
        depth = 0
        level = 0
        colon = -1
        for i in range(q + 1, len(expr)):
            ch = expr[i]
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif depth == 0 and ch == "?":
                level += 1
            elif depth == 0 and ch == ":":
                if level == 0:
                    colon = i
                    break
                level -= 1
        if colon < 0:
            raise ValueError("unpaired ternary in %r" % expr)
        return "(%s) if (%s) else (%s)" % (
            _cpp_expr_to_python(expr[q + 1:colon], names),
            _cpp_expr_to_python(expr[:q], names),
            _cpp_expr_to_python(expr[colon + 1:], names))

    expr = expr.replace("std::fabs", "abs")
    expr = re.sub(r"(?<=[0-9.])f\b", "", expr)
    tokens = re.findall(r"\d[\w.]*|[A-Za-z_][\w:]*", expr)
    unknown = [n for n in tokens if not n[0].isdigit() and n not in names]
    if unknown:
        raise ValueError("unsupported identifiers %s in %r" % (unknown, expr))
    return expr


def _cpp_float_fn(body, params):
    """Compiles a straight-line C++ float function body into a callable.

    Accepts float locals, brace-less single-statement ifs and returns -- the
    shape the mirrored helpers are written in. Anything else raises ValueError,
    so a rewrite on the C++ side surfaces as a parse failure rather than as a
    silently stale mirror.
    """
    text = re.sub(r"//[^\n]*", " ", body).strip().lstrip("{")
    text = " ".join(text.split())
    names = set(params) | {"abs"}
    program = []
    while text:
        conds = []
        while re.match(r"if\s*\(", text):
            cond, text = _split_paren(text[text.index("("):])
            conds.append(compile(_cpp_expr_to_python(cond, names),
                                 "<cpp>", "eval"))
            text = text.lstrip()
        stmt, _, text = text.partition(";")
        stmt = stmt.strip()
        text = text.strip()
        ret = re.match(r"^return\s+(.+)$", stmt)
        decl = re.match(r"^(?:const\s+)?float\s+(\w+)\s*=\s*(.+)$", stmt)
        if ret:
            target, source = None, ret.group(1)
        elif decl:
            target, source = decl.group(1), decl.group(2)
        else:
            raise ValueError("unsupported statement %r" % stmt)
        code = compile(_cpp_expr_to_python(source, names), "<cpp>", "eval")
        program.append((conds, target, code))
        if target is not None:
            names.add(target)

    def call(*args):
        env = {"abs": abs}
        env.update(zip(params, args))
        for conds, target, code in program:
            if not all(eval(c, env) for c in conds):
                continue
            if target is None:
                return eval(code, env)
            env[target] = eval(code, env)
        raise ValueError("control fell off the end of the parsed body")

    return call


def check_angle_roundtrip():
    """Pins diamond_direction() as the inverse of diamond_angle() above.

    Every cell of the table is filled along a direction from
    diamond_direction(), and every runtime lookup indexes it by the diamond
    angle. Both sides live here, so this catches only a half-applied edit to
    this file; the C++ the runtime actually calls is pinned separately by
    check_mirrors().
    """
    n = ANGLE_STEPS * SUBSAMPLES
    t = np.arange(n) * (4.0 / n)
    a_dir, b_dir = diamond_direction(t)
    err = float(np.abs(diamond_angle(b_dir, a_dir) - t).max())
    if err <= 1e-12:
        return True
    sys.stderr.write("diamond_direction() no longer inverts diamond_angle()"
                     " (worst error %g)\n  re-derive both from"
                     " core/math/3dmath.h\n" % err)
    return False


def check_diamond_angle_mirror(math_h_path):
    """Diffs the mirrored diamond angle against 3dmath.h's definition.

    Parses the C++ into a callable and sweeps it against the mirror, so a
    re-parameterization on the side every runtime cell lookup indexes by fails
    here instead of silently re-indexing the whole table.
    """
    with open(math_h_path, "r", encoding="utf-8") as f:
        text = f.read()
    try:
        sig = re.search(r"inline float diamond_angle\(([^)]*)\)", text)
        if sig is None:
            raise ValueError("no float diamond_angle(...) declaration")
        params = tuple(p.split()[-1] for p in sig.group(1).split(","))
        if params != ("y", "x"):
            raise ValueError("parameters are %s, not ('y', 'x')" % (params,))
        cpp = _cpp_float_fn(_function_body(text, "inline float diamond_angle("),
                            params)
    except ValueError as exc:
        sys.stderr.write("could not parse 3dmath.h diamond_angle(): %s\n"
                         "  re-derive the mirror above from that definition\n"
                         % exc)
        return False

    theta = np.arange(721) * (2.0 * np.pi / 721.0)
    # Three radii to cover the scale invariance, plus the degenerate origin and
    # a point inside the 1e-20 guard.
    x = np.concatenate([r * np.cos(theta) for r in (1e-3, 1.0, 7.5)]
                       + [np.array([0.0, 1e-30, -1e-30])])
    y = np.concatenate([r * np.sin(theta) for r in (1e-3, 1.0, 7.5)]
                       + [np.array([0.0, -1e-30, 1e-30])])
    mirror = diamond_angle(y, x)
    err = max(abs(cpp(float(yi), float(xi)) - float(mi))
              for xi, yi, mi in zip(x, y, mirror))
    if err <= 1e-12:
        return True
    sys.stderr.write("diamond_angle() mirror drift against 3dmath.h"
                     " (worst error %g)\n" % err)
    return False


def check_mirrors(color_h_path, math_h_path):
    """Diffs the mirrored OKLab matrices and gamut epsilon against color.h, and
    the mirrored diamond angle against 3dmath.h."""
    with open(color_h_path, "r", encoding="utf-8") as f:
        text = f.read()

    ok = True
    pairs = (("inline void oklab_to_lms_cbrt(", OKLAB_TO_LMS),
             ("inline void lms_cbrt_to_linear_rgb(", LMS_TO_RGB))
    for signature, expected in pairs:
        body = _function_body(text, signature)
        found = [float(s + v)
                 for s, v in re.findall(r"([-+])\s*(\d+\.\d+)f", body)]
        want = [c for row in expected for c in row]
        if found != want:
            sys.stderr.write("mirror drift in %s\n  color.h: %s\n  here:    %s\n"
                             % (signature, found, want))
            ok = False

    gamut = _function_body(text, "inline bool linear_rgb_in_gamut(")
    m = re.search(r"lo\s*=\s*(%s)f\s*,\s*hi\s*=\s*(%s)f\s*\+\s*(%s)f"
                  % (NUM, NUM, NUM), gamut)
    if not m:
        sys.stderr.write("could not parse linear_rgb_in_gamut() slack\n")
        return False
    lo, one, eps = (float(m.group(i)) for i in (1, 2, 3))
    if (lo, one + eps) != (GAMUT_LO, GAMUT_HI):
        sys.stderr.write("gamut slack drift\n  color.h: %r %r\n  here:    %r %r\n"
                         % (lo, one + eps, GAMUT_LO, GAMUT_HI))
        ok = False
    return check_diamond_angle_mirror(math_h_path) and ok


def check_provenance(committed_path):
    """Diffs a fresh table against the committed header in full.

    Whole text, not numeric tokens: the array names, the flash-section marker,
    the constants, the includes and the doc comments are as much a divergence
    as a shifted value, and the header is excluded from the clang-format job
    (it is emitted, not hand-formatted), so this gate is the only thing that
    reads them.
    """
    if not os.path.exists(committed_path):
        sys.stderr.write("missing %s\n" % committed_path)
        return False
    # newline="": universal newlines would read a CRLF-committed header as LF
    # and compare equal to the LF-pinned render().
    with open(committed_path, "r", encoding="utf-8", newline="") as f:
        committed = f.read()

    generated = render(build_table()[0])
    if generated == committed:
        return True
    diff = difflib.unified_diff(committed.splitlines(True),
                                generated.splitlines(True),
                                committed_path, "generated")
    sys.stderr.writelines(itertools.islice(diff, 40))
    sys.stderr.write(
        "%s is out of sync with tools/gen_gamut_lut.py\n"
        "Regenerate with: python tools/gen_gamut_lut.py\n" % committed_path)
    return False


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    parser = argparse.ArgumentParser(
        description="Generate the sRGB gamut boundary bracket table.")
    parser.add_argument(
        "output_path", nargs="?",
        default=os.path.join(root, "core", "color", "gamut_lut.h"),
        help="header to write, or to diff against under --check")
    parser.add_argument(
        "--check", action="store_true",
        help="diff a fresh table against the header and pin the mirrored"
             " constants against color.h and 3dmath.h instead of writing")
    args = parser.parse_args()

    if args.check:
        ok = check_angle_roundtrip()
        ok = check_mirrors(os.path.join(root, "core", "color", "color.h"),
                           os.path.join(root, "core", "math", "3dmath.h")) and ok
        ok = check_provenance(args.output_path) and ok
        sys.exit(0 if ok else 1)

    if not check_angle_roundtrip():
        sys.exit(1)
    table, worst_width = build_table()
    with open(args.output_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(render(table))
    sys.stderr.write(
        "wrote %s (%d x %d cells, %d entries)\n"
        "worst bracket width at full resolution: %.6f\n"
        % (args.output_path, ANGLE_STEPS, L_STEPS, table.size, worst_width))


if __name__ == "__main__":
    main()
