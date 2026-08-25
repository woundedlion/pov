/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file recipe_types.h
 * @brief The recipe model: the authored Conway operator set, one authored step,
 *        and the chain of steps a registry generator mirrors.
 */

#include <cstdint>

namespace MeshOps {
struct RelaxBake;
}

namespace Solids {

/**
 * @brief Authored Conway operator, including the composite ops.
 * @details expand_to_primitives() (mesh/recipe.h) lowers composites to the
 * primitives plus AMBO; authored recipes keep the composite.
 */
enum class Op : uint8_t {
  TRUNCATE,
  EXPAND,
  SNUB,
  CHAMFER,
  HANKIN,
  RELAX,
  KIS,
  DUAL,
  AMBO,
  BEVEL,
  GYRO,
  META,
  NEEDLE,
  ZIP
};

/**
 * @brief One authored step in a recipe's op chain.
 */
struct OpStep {
  Op op; /**< Operator applied at this step. */
  /**
   * @brief t / contact angle (radians) / RELAX iterations.
   * @details Unread on a RELAX step carrying a `bake`; such steps leave it at
   * zero rather than naming a count the replay never runs. A bake-less RELAX
   * step must name at least one iteration: the zero default would replay as a
   * normalize-only pass-through, so apply_step traps on it.
   */
  float param = 0.0f;
  float twist = 0.0f; /**< SNUB face rotation, radians. */
  /**
   * @brief RELAX bake this step lands on, mirroring the generator's
   * relax_baked() call; null replays `param` live iterations.
   * @details When set, every replay of the step (bitwise gate, on-screen build
   * leg, and eager clean endpoint) resolves through the baked converged mesh
   * instead of `param` smoothing steps, so the recipe lands on the shipped
   * geometry rather than a mid-convergence freeze. Left null for a mid-chain
   * relax on a mesh with no bake, which keeps iterating.
   */
  const MeshOps::RelaxBake *bake = nullptr;
};

/**
 * @brief Declarative op chain mirroring a registry generator.
 * @details Parallel declaration of a generator's chain, proven bitwise-equal
 * to it in tests/test_solids.h. The generators stay the source of truth for
 * shipping geometry.
 */
struct Recipe {
  uint8_t seed;        /**< simple_registry index of the base solid. */
  const OpStep *steps; /**< Authored op chain, applied left to right. */
  uint8_t count;       /**< Number of steps. */
};

} // namespace Solids
