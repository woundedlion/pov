/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/math/3dmath.h"

/**
 * @file translation_lattice.h
 * @brief Translation symmetry a planar source reports to the warp before it.
 */

namespace Pullback {

/** @brief Basis of a planar source's translation symmetry group, in plane
    units; a zero basis reports no lattice. */
struct TranslationLattice {
  math::Complex a; ///< First basis vector.
  math::Complex b; ///< Second basis vector.

  /** @brief Whether a lattice is reported.
      @return False for the zero basis. */
  __attribute__((always_inline)) constexpr bool present() const {
    return a != math::Complex{} || b != math::Complex{};
  }
};

/** @brief Axis-aligned square lattice.
    @param period Side of one cell, in plane units.
    @return Basis (period, 0), (0, period). */
constexpr TranslationLattice square_lattice(float period) {
  return {{period, 0.0f}, {0.0f, period}};
}

/** @brief Hexagonal lattice with one basis vector on the x axis.
    @param period Basis vector length, in plane units.
    @return Basis (period, 0), (period / 2, period * sqrt(3) / 2). */
constexpr TranslationLattice hexagonal_lattice(float period) {
  return {{period, 0.0f}, {0.5f * period, 0.8660254037844386f * period}};
}

} // namespace Pullback
