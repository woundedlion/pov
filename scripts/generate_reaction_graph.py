#!/usr/bin/env python3
"""Regenerate core/engine/reaction_graph.cpp — the Fibonacci-lattice K-NN adjacency table.

The reaction-diffusion effects (BZ / Gray-Scott) diffuse over a 7680-point
Fibonacci lattice on the unit sphere. core/engine/reaction_graph.cpp is the checked-in
K-nearest-neighbor table for that lattice: neighbors[i][k] is the index of the
k-th nearest node to node i, ordered by increasing distance. This script is its
generator of record. The runtime node() in core/engine/reaction_graph.h MUST reproduce
the same lattice positions this script uses, or find_nearest_node / CubemapLUT
would descend a table built from points that disagree with the ones it samples.

Provenance — DOUBLE PRECISION. node(i) is evaluated in IEEE double here, and the
committed table reproduces bit-for-bit ONLY in double: computing y / radius in
float32 (as a naive runtime port would) flips the sort order of a handful of rows
at near-tie distances. core/engine/reaction_graph.h::node() therefore folds y, radius,
and theta in double and narrows once into the float Vector, symmetric with this
generator. RD_N, RD_K, golden_angle and two_pi are parsed out of that header
rather than restated here, so only the node() formula itself is mirrored by hand;
any edit to it requires regenerating the table.

Lattice (mirrors core/engine/reaction_graph.h::node):

  y(i)      = 1 - i/(RD_N-1) * 2                  # +1 at the north pole down to -1
  radius(i) = sqrt(1 - y*y)
  theta(i)  = fmod(golden_angle * i, 2π)
  node(i)   = (cos(theta)*radius, y, sin(theta)*radius)

neighbors[i] = the RD_K indices j != i minimizing |node(i) - node(j)|^2, sorted
by (distance, index) so ties break deterministically toward the lower index.

Usage:
  python scripts/generate_reaction_graph.py > core/engine/reaction_graph.cpp

The CI provenance gate (.github/workflows/ci.yml :: reaction-graph-provenance)
re-runs this and diffs the numeric tokens against the committed file, so a drift
in either the script or the table fails loudly.
"""

import math
import re
import sys
from pathlib import Path

HEADER = (Path(__file__).resolve().parent.parent / "core" / "engine" /
          "reaction_graph.h")


def _header_constant(text, pattern, name):
    match = re.search(pattern, text)
    if not match:
        raise RuntimeError(
            f"{name} not found in {HEADER}; this generator reads the runtime "
            "header as the single source of truth for the lattice constants")
    return match.group(1)


def _load_header_constants():
    """Lattice constants read out of core/engine/reaction_graph.h.

    Restating them here would let a header-only edit produce a table the CI
    provenance diff still accepts (the C++ array bound comes from the header, so
    a raised RD_N zero-pads the missing rows instead of failing to compile).
    """
    text = HEADER.read_text(encoding="utf-8")
    return (
        int(_header_constant(text, r"constexpr\s+int\s+RD_N\s*=\s*(\d+)\s*;",
                             "RD_N")),
        int(_header_constant(text, r"constexpr\s+int\s+RD_K\s*=\s*(\d+)\s*;",
                             "RD_K")),
        float(_header_constant(
            text, r"constexpr\s+double\s+golden_angle\s*=\s*([0-9.eE+-]+)\s*;",
            "golden_angle")),
        float(_header_constant(
            text, r"constexpr\s+double\s+two_pi\s*=\s*([0-9.eE+-]+)\s*;",
            "two_pi")),
    )


RD_N, RD_K, GOLDEN_ANGLE, TWO_PI = _load_header_constants()


def node(i):
    """Fibonacci-lattice node i as a unit vector, in double precision."""
    y = 1.0 - (i / (RD_N - 1)) * 2.0
    radius = math.sqrt(1.0 - y * y)
    theta = math.fmod(GOLDEN_ANGLE * i, TWO_PI)
    return (math.cos(theta) * radius, y, math.sin(theta) * radius)


def build_neighbors():
    pos = [node(i) for i in range(RD_N)]

    # The lattice is ordered north-to-south, so a node's index tracks its y
    # coordinate monotonically (dy per index step = 2/(RD_N-1)). Every one of the
    # RD_K nearest neighbors sits within ~11 deg (chord <~ 0.193), hence within
    # |dy| <~ 0.193, so a fixed index window around i must contain them all. W is
    # sized so the window's y half-width (W * dy_step) comfortably exceeds that,
    # and the assertion below proves sufficiency per row — if a future RD_N change
    # outgrows the window the generator aborts instead of emitting a wrong table.
    dy_step = 2.0 / (RD_N - 1)
    W = 1000
    window_half_y = W * dy_step

    table = []
    for i in range(RD_N):
        xi, yi, zi = pos[i]
        lo = max(0, i - W)
        hi = min(RD_N, i + W + 1)
        cand = []
        for j in range(lo, hi):
            if j == i:
                continue
            xj, yj, zj = pos[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            cand.append((dx * dx + dy * dy + dz * dz, j))
        cand.sort()  # by (squared distance, index): deterministic tie-break
        chosen = cand[:RD_K]
        # Completeness guard: any node strictly closer than the farthest chosen
        # neighbor has |dy| <= its chord < sqrt(farthest_d2) <= window_half_y, so
        # it cannot lie outside the index window. If this ever fails, W is too
        # small for the current RD_N — fail loudly rather than ship a biased table.
        farthest_d2 = chosen[-1][0]
        if math.sqrt(farthest_d2) >= window_half_y:
            raise RuntimeError(
                f"node {i}: neighbor window too small (chord "
                f"{math.sqrt(farthest_d2):.4f} >= half-width {window_half_y:.4f}); "
                "raise W in generate_reaction_graph.py")
        table.append([j for _, j in chosen])
    return table


def main():
    if RD_N > 32767:
        raise RuntimeError(
            f"RD_N ({RD_N}) exceeds INT16_MAX (32767): node index must fit "
            "int16_t; widen the neighbors element type in "
            "engine/reaction_graph.h before raising RD_N")
    table = build_neighbors()
    # LF on every host: the provenance gate diffs numeric tokens only, so a
    # Windows CRLF flip of the whole table would ride green.
    sys.stdout.reconfigure(newline="\n")
    out = sys.stdout
    out.write("// Auto-generated by scripts/generate_reaction_graph.py -- do not edit\n")
    out.write('#include "engine/reaction_graph.h"\n')
    out.write("\n")
    out.write("PROGMEM const int16_t ReactionGraph::neighbors[RD_N][RD_K] = {\n")
    for idx, row in enumerate(table):
        comma = "," if idx + 1 < len(table) else ""
        out.write("  {" + ", ".join(str(v) for v in row) + "}" + comma + "\n")
    out.write("};\n")


if __name__ == "__main__":
    main()
