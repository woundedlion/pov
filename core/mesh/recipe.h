/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstddef>
#include <cstdint>

#include "mesh/conway_graph.h"
#include "mesh/solids.h"

namespace Solids {

/**
 * @brief Whether a lowered primitive step has a morph leg kind covering it.
 * @param step Lowered primitive step.
 * @return True when a leg can sweep the step; false leaves the whole recipe to
 *   the caller's whole-generate fallback.
 * @details TRUNCATE below ConwayGraph::T_EPS clamps both leg endpoints to the
 * same value (a still image), and above 0.5 crosses the ambo short-circuit,
 * where the leg would clean-swap to a self-intersecting form. CHAMFER has no
 * swept characterization. KIS and DUAL are partition ops with no sweep.
 * EXPAND has a leg kind but no recipe and no sweep coverage on a hankin seed.
 */
inline bool is_morphable_step(const OpStep &step) {
  switch (step.op) {
  case Op::TRUNCATE:
    return step.param >= ConwayGraph::T_EPS && step.param <= 0.5f;
  case Op::SNUB:
  case Op::HANKIN:
  case Op::RELAX:
  case Op::AMBO:
    return true;
  default:
    return false;
  }
}

/**
 * @brief Lowers a recipe's authored steps to the eight primitives plus AMBO.
 * @param recipe Recipe whose authored steps are lowered.
 * @param out Output buffer receiving the lowered steps.
 * @param cap Capacity of `out`; traps on overflow.
 * @return Number of steps written.
 * @details Composites lower to the exact chains their MeshOps counterparts run
 * (composition is right-to-left, so the step order is application order):
 * gyro = snub(0.5, 0), dual; meta = ambo, dual, kis; needle = dual, kis;
 * zip = kis, dual; bevel(t) = ambo, truncate(t). bevel at exactly t == 0.5
 * emits ambo, ambo instead, matching truncate's ambo short-circuit. Primitive
 * steps pass through unchanged.
 */
inline size_t expand_to_primitives(const Recipe &recipe, OpStep *out,
                                   size_t cap) {
  size_t n = 0;
  auto emit = [&](OpStep step) {
    HS_CHECK(n < cap, "expand_to_primitives: output capacity exceeded");
    out[n++] = step;
  };
  for (size_t i = 0; i < recipe.count; ++i) {
    const OpStep &step = recipe.steps[i];
    switch (step.op) {
    case Op::GYRO:
      emit({Op::SNUB, 0.5f, 0.0f});
      emit({Op::DUAL});
      break;
    case Op::META:
      emit({Op::AMBO});
      emit({Op::DUAL});
      emit({Op::KIS});
      break;
    case Op::NEEDLE:
      emit({Op::DUAL});
      emit({Op::KIS});
      break;
    case Op::ZIP:
      emit({Op::KIS});
      emit({Op::DUAL});
      break;
    case Op::BEVEL:
      emit({Op::AMBO});
      if (step.param == 0.5f)
        emit({Op::AMBO});
      else
        emit({Op::TRUNCATE, step.param});
      break;
    default:
      emit(step);
      break;
    }
  }
  return n;
}

/**
 * @brief Replays a recipe's authored chain from its simple_registry seed.
 * @param recipe Recipe to replay; composites run through the SolidBuilder
 *   composite methods, so the replay matches the generator functions exactly.
 * @param a Output arena for even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rebuilt PolyMesh.
 */
FLASHMEM static PolyMesh build_recipe(const Recipe &recipe, Arena &a,
                                      Arena &b) {
  SolidBuilder builder(simple_registry[recipe.seed].generate(a, b), a, b);
  for (size_t i = 0; i < recipe.count; ++i) {
    const OpStep &step = recipe.steps[i];
    switch (step.op) {
    case Op::TRUNCATE:
      builder.truncate(step.param);
      break;
    case Op::EXPAND:
      builder.expand(step.param);
      break;
    case Op::SNUB:
      builder.snub(step.param, step.twist);
      break;
    case Op::CHAMFER:
      builder.chamfer(step.param);
      break;
    case Op::HANKIN:
      builder.hankin(step.param);
      break;
    case Op::RELAX:
      builder.relax(static_cast<int>(step.param));
      break;
    case Op::KIS:
      builder.kis();
      break;
    case Op::DUAL:
      builder.dual();
      break;
    case Op::AMBO:
      builder.ambo();
      break;
    case Op::BEVEL:
      builder.bevel(step.param);
      break;
    case Op::GYRO:
      builder.gyro();
      break;
    case Op::META:
      builder.meta();
      break;
    case Op::NEEDLE:
      builder.needle();
      break;
    case Op::ZIP:
      builder.zip();
      break;
    }
  }
  return builder.build();
}

/**
 * @brief Replays a primitive step list from a simple_registry seed.
 * @param seed simple_registry index of the base solid.
 * @param steps Primitive steps (the eight primitives plus AMBO); a composite
 *   step traps.
 * @param count Number of steps.
 * @param a Output arena for even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rebuilt PolyMesh.
 */
FLASHMEM static PolyMesh build_steps(uint8_t seed, const OpStep *steps,
                                     size_t count, Arena &a, Arena &b) {
  SolidBuilder builder(simple_registry[seed].generate(a, b), a, b);
  for (size_t i = 0; i < count; ++i) {
    const OpStep &step = steps[i];
    switch (step.op) {
    case Op::TRUNCATE:
      builder.truncate(step.param);
      break;
    case Op::EXPAND:
      builder.expand(step.param);
      break;
    case Op::SNUB:
      builder.snub(step.param, step.twist);
      break;
    case Op::CHAMFER:
      builder.chamfer(step.param);
      break;
    case Op::HANKIN:
      builder.hankin(step.param);
      break;
    case Op::RELAX:
      builder.relax(static_cast<int>(step.param));
      break;
    case Op::KIS:
      builder.kis();
      break;
    case Op::DUAL:
      builder.dual();
      break;
    case Op::AMBO:
      builder.ambo();
      break;
    default:
      HS_CHECK(false, "build_steps: non-primitive op");
      break;
    }
  }
  return builder.build();
}

} // namespace Solids
