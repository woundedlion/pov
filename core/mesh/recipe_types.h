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

#include <cstddef>
#include <cstdint>

namespace MeshOps {
struct RelaxBake;
}

namespace Solids {

/**
 * @brief Authored Conway operator, including the composite ops.
 * @details Authored recipes keep the composite; expand_to_primitives() lowers
 * it.
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
   * @details Unread, and left at zero, on a RELAX step carrying a `bake`. A
   * bake-less RELAX step must name at least one iteration.
   */
  float param = 0.0f;
  float twist = 0.0f; /**< SNUB face rotation, radians. */
  /**
   * @brief RELAX bake this step lands on, mirroring the generator's
   * relax_baked() call; null replays `param` live iterations.
   * @details When set, every replay of the step resolves through the baked
   * mesh. Left null for a relax on a mesh with no bake.
   */
  const MeshOps::RelaxBake *bake = nullptr;
};

/**
 * @brief Declarative op chain mirroring a registry generator.
 * @details The generators are the source of truth for shipping geometry.
 */
struct Recipe {
  uint8_t seed;        /**< simple_registry index of the base solid. */
  const OpStep *steps; /**< Authored op chain, applied left to right. */
  uint8_t count;       /**< Number of steps. */
};

/** @brief Builds a recipe whose count is deduced from its step array. */
template <size_t N>
constexpr Recipe make_recipe(uint8_t seed, const OpStep (&steps)[N]) {
  static_assert(N <= UINT8_MAX);
  return {seed, steps, static_cast<uint8_t>(N)};
}

} // namespace Solids
