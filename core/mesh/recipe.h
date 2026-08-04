/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file recipe.h
 * @brief Lowering and replay of Solids::Recipe op chains, plus the leg-kind
 *        coverage rules the morph animation depends on.
 */

#include <cstddef>
#include <cstdint>

#include "mesh/conway_graph.h"
#include "mesh/solids.h"

namespace Solids {

/** Largest chamfer thickness a leg is characterized at: the sweep, birth-limit
 * and birth-epsilon probes (tests/test_opchain_probe.h) run T_EPS -> 0.63 on
 * every chamfer seed. */
inline constexpr float CHAMFER_T_MAX = 0.63f;

/**
 * @brief Whether a lowered primitive step has a morph leg kind covering it.
 * @param step Lowered primitive step.
 * @return True when a leg can sweep or gate the step; false leaves the whole
 *   recipe to the caller's whole-generate fallback.
 * @details TRUNCATE below ConwayGraph::T_EPS_TRUNCATE_MIN sweeps too few pixels
 * to read as motion. A sub-T_EPS arrival (0.01) still sweeps: the leg births at
 * the derived per-arrival floor (min(T_EPS, arrival * T_EPS_TRUNCATE_FRAC))
 * rather than clamping both endpoints to T_EPS. Above 0.5 is a far-side leg: it
 * sweeps through the ambo pinch on the constant-topology truncate branch (the
 * two truncate50d recipes arrive at 0.873), up to
 * T_EPS_TRUNCATE_FAR_MAX; behaviour in [T_EPS_TRUNCATE_MIN, 0.5] is unchanged.
 * CHAMFER is characterized up to CHAMFER_T_MAX. KIS and DUAL run as gated swaps
 * (docs/opchain_morph_spec.md, section 3.3). EXPAND has a leg kind but no
 * recipe and no sweep coverage on a hankin seed.
 */
inline bool is_morphable_step(const OpStep &step) {
  switch (step.op) {
  case Op::TRUNCATE:
    return step.param >= ConwayGraph::T_EPS_TRUNCATE_MIN &&
           step.param <= ConwayGraph::T_EPS_TRUNCATE_FAR_MAX;
  case Op::CHAMFER:
    return step.param >= ConwayGraph::T_EPS && step.param <= CHAMFER_T_MAX;
  case Op::SNUB:
  case Op::HANKIN:
  case Op::RELAX:
  case Op::AMBO:
  case Op::KIS:
  case Op::DUAL:
    return true;
  case Op::EXPAND:
  case Op::BEVEL:
  case Op::GYRO:
  case Op::META:
  case Op::NEEDLE:
  case Op::ZIP:
    return false;
  }
  return false;
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
HS_COLD_MEMBER inline size_t expand_to_primitives(const Recipe &recipe,
                                                  OpStep *out, size_t cap) {
  size_t n = 0;
  auto emit = [&](OpStep step) HS_COLD_MEMBER {
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
 * @brief Applies one step to a builder.
 * @param builder Builder the step runs on.
 * @param step Step to apply.
 * @param allow_composite Whether composite steps (bevel, gyro, meta, needle,
 *   zip) are legal; a composite step traps when false.
 * @details Single dispatch shared by both replay paths, so an authored replay
 * and a lowered primitive replay cannot diverge on a step.
 */
FLASHMEM static void apply_step(SolidBuilder &builder, const OpStep &step,
                                bool allow_composite) {
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
    if (step.bake)
      builder.relax_baked(*step.bake);
    else
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
    HS_CHECK(allow_composite, "apply_step: non-primitive op");
    builder.bevel(step.param);
    break;
  case Op::GYRO:
    HS_CHECK(allow_composite, "apply_step: non-primitive op");
    builder.gyro();
    break;
  case Op::META:
    HS_CHECK(allow_composite, "apply_step: non-primitive op");
    builder.meta();
    break;
  case Op::NEEDLE:
    HS_CHECK(allow_composite, "apply_step: non-primitive op");
    builder.needle();
    break;
  case Op::ZIP:
    HS_CHECK(allow_composite, "apply_step: non-primitive op");
    builder.zip();
    break;
  }
}

/**
 * @brief Replays a recipe's authored chain from its simple_registry seed.
 * @param recipe Recipe to replay; composites run through the SolidBuilder
 *   composite methods, so the replay matches the generator functions exactly.
 * @param a Output arena for even pipeline stages.
 * @param b Scratch arena for odd pipeline stages.
 * @return The rebuilt PolyMesh.
 */
[[maybe_unused]] FLASHMEM static PolyMesh build_recipe(const Recipe &recipe,
                                                       Arena &a, Arena &b) {
  HS_CHECK(recipe.seed < std::size(simple_registry),
           "build_recipe: seed outside simple_registry");
  SolidBuilder builder(simple_registry[recipe.seed].generate(a, b), a, b);
  for (size_t i = 0; i < recipe.count; ++i)
    apply_step(builder, recipe.steps[i], true);
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
  HS_CHECK(seed < std::size(simple_registry),
           "build_steps: seed outside simple_registry");
  SolidBuilder builder(simple_registry[seed].generate(a, b), a, b);
  for (size_t i = 0; i < count; ++i)
    apply_step(builder, steps[i], false);
  return builder.build();
}

} // namespace Solids
