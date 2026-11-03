#!/usr/bin/env python3
"""Generate peak-speed corrections and covariant derivative bounds for fields."""

import argparse
import difflib
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "core/math/tangent_field_constants.h"
LATITUDE_STEPS = 1024
LONGITUDE_STEPS = 512
PEAK_MARGIN = 1.01
DERIVATIVE_MARGIN = 1.001


def tetrahedral_axes():
    """Mirror the constexpr half-turn extraction in tangent_field.h."""
    mirrors = np.array([[1, 0, 0], [-0.5, 0.8660254038, 0],
                        [0, -0.5773502692, 0.8164965809]], dtype=np.float32)
    group = [np.eye(3, dtype=np.float32)]
    for matrix in group:
        for normal in mirrors:
            candidate = matrix - np.float32(2) * np.outer(normal, normal @ matrix)
            if not any(np.max(np.abs(candidate - other)) < 2e-5 for other in group):
                group.append(candidate)
    axes = []
    for matrix in group:
        if abs(np.trace(matrix) + 1) > 2e-5 or np.max(np.abs(matrix - matrix.T)) > 2e-5:
            continue
        columns = matrix + np.eye(3, dtype=np.float32)
        axis = columns[:, np.argmax(np.sum(columns * columns, axis=0))]
        axes.append(axis.astype(float))
    if len(axes) != 3:
        raise ValueError("tetrahedral half-turn extraction failed")
    return np.array(axes)


def mode_polynomial(l, m):
    """Reduced Legendre polynomial, from the defining recurrence."""
    previous = Polynomial([0])
    current = Polynomial([1])
    for degree in range(m, l):
        previous, current = current, ((2 * degree + 1) * Polynomial([0, 1]) * current
                                      - (degree + m) * previous) / (degree + 1 - m)
    return current


def mode_derivatives(l, m, p):
    """Ambient gradient and Hessian of Q_lm(y) Re((x+iz)^m)."""
    x, y, z = p.T
    w = x + 1j * z
    q = mode_polynomial(l, m)
    power = w ** m
    first = m * w ** max(0, m - 1)
    second = m * (m - 1) * w ** max(0, m - 2)
    g = np.array([q(y) * first.real, q.deriv()(y) * power.real,
                  -q(y) * first.imag]).T
    h = np.zeros((len(p), 3, 3))
    h[:, 0, 0] = q(y) * second.real
    h[:, 2, 2] = -h[:, 0, 0]
    h[:, 0, 2] = h[:, 2, 0] = -q(y) * second.imag
    h[:, 0, 1] = h[:, 1, 0] = q.deriv()(y) * first.real
    h[:, 2, 1] = h[:, 1, 2] = -q.deriv()(y) * first.imag
    h[:, 1, 1] = q.deriv(2)(y) * power.real
    return g, h


def product_derivatives(values, gradients, hessians):
    """Gradient and Hessian of the product of three scalar factors."""
    g = sum(gradients[i] * (values[(i + 1) % 3] * values[(i + 2) % 3])[:, None]
            for i in range(3))
    h = sum(hessians[i] * (values[(i + 1) % 3] * values[(i + 2) % 3])[:, None, None]
            for i in range(3))
    for i in range(3):
        for j in range(3):
            if i != j:
                h += np.einsum("ni,nj,n->nij", gradients[i], gradients[j], values[3 - i - j])
    return g, h


def invariant_derivatives(kind, p):
    """Ambient derivatives of the fold-aligned invariant polynomials."""
    if kind == "tetrahedral":
        axes = tetrahedral_axes()
        return product_derivatives([p @ a for a in axes],
                                   [np.broadcast_to(a, p.shape) for a in axes],
                                   [np.zeros((len(p), 3, 3)) for _ in axes])
    if kind == "octahedral":
        h = np.zeros((len(p), 3, 3))
        h[:, range(3), range(3)] = 12 * p * p
        return 4 * p ** 3, h
    phi2 = ((1 + math.sqrt(5)) / 2) ** 2
    diagonals = np.array([[phi2, -1, 0], [0, phi2, -1], [-1, 0, phi2]])
    return product_derivatives([np.sum(p * p * a, axis=1) for a in diagonals],
                               [2 * p * a for a in diagonals],
                               [np.broadcast_to(np.diag(2 * a), (len(p), 3, 3))
                                for a in diagonals])


def sampling_grid():
    """Uniform angular grid including both poles and every coordinate axis."""
    latitude = np.linspace(-math.pi / 2, math.pi / 2, LATITUDE_STEPS + 1)
    longitude = np.arange(LONGITUDE_STEPS) * (2 * math.pi / LONGITUDE_STEPS)
    lat, lon = np.meshgrid(latitude, longitude, indexing="ij")
    lat, lon = lat.ravel(), lon.ravel()
    radius, y = np.cos(lat), np.sin(lat)
    cosine, sine = np.cos(lon), np.sin(lon)
    p = np.array([radius * cosine, y, radius * sine]).T
    e = np.array([-sine, np.zeros_like(sine), cosine]).T
    f = np.array([-y * cosine, radius, -y * sine]).T
    return p, e, f


def measure(derivatives, grid):
    """Peak gradient and operator norm of its covariant derivative."""
    p, e, f = grid
    g, h = derivatives
    radial = np.sum(p * g, axis=1)
    speed = np.sqrt(np.sum((g - p * radial[:, None]) ** 2, axis=1)).max()
    a = np.einsum("ni,nij,nj->n", e, h, e) - radial
    b = np.einsum("ni,nij,nj->n", e, h, f)
    d = np.einsum("ni,nij,nj->n", f, h, f) - radial
    derivative = (0.5 * np.abs(a + d) + np.sqrt(0.25 * (a - d) ** 2 + b * b)).max()
    return speed, derivative


def build_constants():
    """Measure all terms and retain a guard on the sampled maxima."""
    grid = sampling_grid()
    modes = []
    for l in range(9):
        for m in range(l + 1):
            if l == 0:
                modes.append((0.0, 0.0))
                continue
            peak, derivative = measure(mode_derivatives(l, m, grid[0]), grid)
            seed = math.prod(range(1, 2 * m, 2))
            scale = math.sqrt((2 * l + 1) / (4 * math.pi)
                              * math.factorial(l - m) / math.factorial(l + m))
            if m:
                scale *= math.sqrt(2)
            correction = 1 / (PEAK_MARGIN * peak * scale * seed)
            modes.append((correction, DERIVATIVE_MARGIN * derivative / (PEAK_MARGIN * peak)))
    invariants = []
    for kind in ("tetrahedral", "octahedral", "dodecahedral"):
        peak, derivative = measure(invariant_derivatives(kind, grid[0]), grid)
        invariants.append((1 / (PEAK_MARGIN * peak),
                           DERIVATIVE_MARGIN * derivative / (PEAK_MARGIN * peak)))
    return modes, invariants


def render(constants):
    """Render the deterministic generated C++ header."""
    modes, invariants = constants
    lines = ["/*", " * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.",
             " * Licensed under the PolyForm Noncommercial License 1.0.0", " */",
             "#pragma once", "", "/** @file tangent_field_constants.h",
             " * @brief Generated by tools/gen_tangent_field_constants.py. Do not edit.",
             " */", "", "#include <array>", "", "namespace math::tangent_field {",
             "/** @brief Peak correction and normalized covariant derivative bound. */",
             "struct TermConstants {", "  float scale;     ///< Peak-speed scale correction.",
             "  float lipschitz; ///< Covariant derivative operator-norm bound.", "};",
             "/** @brief Mode constants indexed by l*(l+1)/2+m. */",
             "inline constexpr std::array<TermConstants, 45> MODE_CONSTANTS = {{"]
    entries = [f"{{{correction:.9e}f, {derivative:.9e}f}}," for correction, derivative in modes]
    for index in range(0, len(entries), 2):
        lines.append("    " + " ".join(entries[index:index + 2]))
    lines += ["}};", "/** @brief Tetrahedral, octahedral and dodecahedral constants. */",
              "inline constexpr std::array<TermConstants, 3> INVARIANT_CONSTANTS = {{"]
    for correction, derivative in invariants:
        lines.append(f"    {{{correction:.9e}f, {derivative:.9e}f}},")
    lines += ["}};", "} // namespace math::tangent_field", ""]
    return "\n".join(lines)


def main():
    """Write or verify the generated header."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=HEADER)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = render(build_constants())
    if args.check:
        committed = args.output.read_text(encoding="utf-8")
        if committed != content:
            print("".join(difflib.unified_diff(committed.splitlines(True), content.splitlines(True))))
            raise SystemExit(1)
    else:
        args.output.write_text(content, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
