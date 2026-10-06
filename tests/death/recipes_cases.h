/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Recipes death fixtures and guard cases.

/** @brief Death case: an out-of-range solids index must trap. */
inline void case_solids_index_oob() {
  const auto &e = Solids::get_entry(opaque<size_t>(Solids::NUM_ENTRIES));
  if (e.name == nullptr)
    std::printf("x");
}

/** @brief Death case: looking up an unknown solid name must trap. */
inline void case_solids_unknown_name() {
  PolyMesh m = Solids::get_by_name(persistent_arena, scratch_arena_a,
                                   scratch_arena_b, "definitely_not_a_solid");
  if (m.vertices.size() == 0x7fff)
    std::printf("x");
}

inline void case_recipe_bake_wrong_op() {
  const MeshOps::RelaxBake bake{};
  run_invalid_recipe_step({Solids::Op::AMBO, 0.0f, 0.0f, &bake});
}

inline void case_recipe_twist_wrong_op() {
  run_invalid_recipe_step({Solids::Op::AMBO, 0.0f, opaque(0.1f)});
}

inline void case_recipe_bake_live_iterations() {
  const MeshOps::RelaxBake bake{};
  run_invalid_recipe_step({Solids::Op::RELAX, opaque(1.0f), 0.0f, &bake});
}

/**
 * @brief Death case: a HANKIN step with no contact angle must trap.
 * @details A zero angle collapses every star point onto its corner.
 */
inline void case_apply_step_hankin_no_angle() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::HANKIN, opaque(0.0f)}};
  PolyMesh mesh = Solids::build_steps(opaque<uint8_t>(1), steps, 1, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}

/**
 * @brief Death case: a BEVEL step with no depth must trap.
 * @details The composite lowers to ambo, truncate(t); a zero depth is a
 *          depthless truncate.
 */
inline void case_apply_step_bevel_no_depth() {
  static uint8_t a_buf[64 * 1024];
  static uint8_t b_buf[64 * 1024];
  Arena a(a_buf, sizeof(a_buf));
  Arena b(b_buf, sizeof(b_buf));
  const Solids::OpStep steps[] = {{Solids::Op::BEVEL, opaque(0.0f)}};
  const Solids::Recipe recipe = {opaque<uint8_t>(1), steps, 1};
  PolyMesh mesh = Solids::build_recipe(recipe, a, b);
  if (mesh.vertices.size() == opaque<size_t>(0x7fff))
    std::printf("x");
}
