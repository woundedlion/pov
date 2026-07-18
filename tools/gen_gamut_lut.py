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
"""

import os
import sys

import numpy as np

# Flash master resolution. init_gamut_lut() downsamples by integer factors, so
# these bound the finest grid any effect can request.
ANGLE_STEPS = 512
L_STEPS = 256
# 65535 / 0.5: OKLab chroma inside sRGB stays below 0.5, so this spends the full
# uint16 range on the live domain at ~7.6e-6 resolution.
SCALE = 131070.0
# Sub-samples per cell per axis, closed at both ends.
SUBSAMPLES = 8
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

# Matches core/color/color.h linear_rgb_in_gamut().
GAMUT_LO = -1e-4
GAMUT_HI = 1.0 + 1e-4


def oklab_to_lms_cbrt(L, a, b):
    """Inverse OKLab matrix; mirrors color.h oklab_to_lms_cbrt()."""
    l_ = L + 0.3963377774 * a + 0.2158037573 * b
    m_ = L - 0.1055613458 * a - 0.0638541728 * b
    s_ = L - 0.0894841775 * a - 1.2914855480 * b
    return l_, m_, s_


def lms_cbrt_to_linear_rgb(l_, m_, s_):
    """Cube plus RGB matrix; mirrors color.h lms_cbrt_to_linear_rgb()."""
    l = l_ * l_ * l_
    m = m_ * m_ * m_
    s = s_ * s_ * s_
    r = +4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s
    g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s
    b = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
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


def emit(table, out_path):
    """Writes the generated header."""
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
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdint>

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
inline const uint16_t GAMUT_LUT[GAMUT_LUT_ENTRIES] = {{
{body}
}};
"""
    with open(out_path, "w", newline="\n") as f:
        f.write(header)


def main():
    if len(sys.argv) > 1:
        out_path = sys.argv[1]
    else:
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        out_path = os.path.join(root, "core", "color", "gamut_lut.h")

    table, worst_width = build_table()
    emit(table, out_path)
    sys.stderr.write(
        "wrote %s (%d x %d cells, %d entries)\n"
        "worst bracket width at full resolution: %.6f\n"
        % (out_path, ANGLE_STEPS, L_STEPS, table.size, worst_width))


if __name__ == "__main__":
    main()
