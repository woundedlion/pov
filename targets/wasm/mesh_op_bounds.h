/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/**
 * @file mesh_op_bounds.h
 * @brief Pure (no-Emscripten) mesh-operator roster, growth factors and
 *        byte-per-element budgets.
 *
 * The MeshOps boundary guards price operator chains against these; a value
 * below an operator's real expansion or footprint lets a chain reach an engine
 * trap.
 */
#pragma once

#include <cstddef>
#include <cstring>
#include <iterator>

namespace hs_wasm {

/**
 * @brief How far one mesh operator grows its input, for the boundary guards.
 * @details Each field multiplies an input measurement. `elements` must be
 *          nonzero. A zero degree or valence means those faces have a fixed
 *          side count or are absent; a zero `valence` also skips the valence
 *          scan.
 */
struct MeshOpBounds {
  size_t elements;    /**< Multiple of the largest input element count. */
  size_t face_degree; /**< Multiple of the widest input face's side count. */
  size_t valence;     /**< Multiple of the highest input vertex valence. */
};

/**
 * @brief Single source of truth for the Conway/Goldberg operator roster and each
 *        operator's growth factors.
 * @param OP0  Macro applied to each zero-argument operator.
 * @param OP1U Macro applied to each [0,1]-fraction operator.
 * @param OP1H Macro applied to each [0,1)-fraction operator.
 * @details Each fraction operator's macro matches the domain its engine trap
 *          asserts: truncate and bevel accept 1 (OP1U); chamfer and expand
 *          assert t < 1 (OP1H). Operators with bespoke signatures are in
 *          MESHOP_IRREGULAR_LIST. The trailing arguments are the operator's
 *          MeshOpBounds, in order.
 *
 *          `elements` is the largest multiple of the input's largest count
 *          (vertices, faces or flat indices) that any stage reaches;
 *          compositions multiply through. `face_degree`
 *          and `valence` are the multiples that reach narrow_face_count.
 */
// clang-format off
#define MESHOP_LIST(OP0, OP1U, OP1H)                                         \
  OP0(kis, 3, 0, 0) OP0(ambo, 2, 1, 1) OP0(gyro, 5, 1, 1)                    \
  OP0(dual, 1, 0, 1) OP0(meta, 6, 1, 1) OP0(needle, 3, 0, 1)                 \
  OP0(zip, 3, 1, 2)                                                            \
  OP1U(truncate, 3, 2, 1) OP1U(bevel, 6, 2, 2)                                \
  OP1H(chamfer, 4, 1, 0) OP1H(expand, 4, 1, 1)
// clang-format on

// Irregular ops: hand-written wrapper methods (custom signatures/validation),
// enumerated here so their embind bindings expand from one list.
#define MESHOP_IRREGULAR_LIST(_) _(relax) _(hankin) _(snub)

/** @brief relax growth factors: it preserves topology. */
inline constexpr MeshOpBounds RELAX_BOUNDS{1, 1, 0};
/** @brief hankin growth factors: it doubles both face degree and valence. */
inline constexpr MeshOpBounds HANKIN_BOUNDS{4, 2, 2};
/** @brief snub growth factors. */
inline constexpr MeshOpBounds SNUB_BOUNDS{5, 1, 1};

/**
 * @brief Worst-case scratch bytes one operator allocates per element of its
 *        largest stage, counting vertices, faces and flat face indices alike.
 * @details The heaviest operator keeps an intermediate mesh, its output and its
 *          index scratch live in one arena.
 */
inline constexpr size_t TOOLING_SCRATCH_BYTES_PER_MESH_ELEMENT = 64;

/**
 * @brief Tooling-arena bytes a finalized mesh retains per element.
 * @details Covers the vertex, side count, index and the topology code
 *          classifyFaces() binds into the same arena, plus alignment slack.
 */
inline constexpr size_t TOOLING_RETAINED_BYTES_PER_MESH_ELEMENT = 20;

/** @brief One operator's declared growth factors, by name. */
struct MeshOpBoundsEntry {
  const char *name;    /**< MeshOps operator name. */
  MeshOpBounds bounds; /**< Factors the boundary guard prices it against. */
};

#define MESHOP_BOUNDS_ENTRY(name, elements, degree, valence)                   \
  {#name, {(elements), (degree), (valence)}},

/**
 * @brief Every operator's declared growth factors as data.
 * @details Regular operators expand from MESHOP_LIST; irregular rows carry the
 *          constants their hand-written wrappers pass.
 */
// clang-format off
inline constexpr MeshOpBoundsEntry MESHOP_BOUNDS[] = {
    MESHOP_LIST(MESHOP_BOUNDS_ENTRY, MESHOP_BOUNDS_ENTRY, MESHOP_BOUNDS_ENTRY)
    {"relax", RELAX_BOUNDS},
    {"hankin", HANKIN_BOUNDS},
    {"snub", SNUB_BOUNDS}};
// clang-format on

#undef MESHOP_BOUNDS_ENTRY

/** @brief Number of rows in MESHOP_BOUNDS. */
inline constexpr size_t MESHOP_BOUNDS_COUNT = std::size(MESHOP_BOUNDS);

// The irregular rows are hand-written: an operator joining
// MESHOP_IRREGULAR_LIST needs a row in MESHOP_BOUNDS.
#define MESHOP_COUNT_ONE(...) +1
static_assert(
    MESHOP_BOUNDS_COUNT ==
        static_cast<size_t>(0 MESHOP_LIST(MESHOP_COUNT_ONE, MESHOP_COUNT_ONE,
                                          MESHOP_COUNT_ONE)
                                MESHOP_IRREGULAR_LIST(MESHOP_COUNT_ONE)),
    "MESHOP_BOUNDS needs one row per roster operator");
#undef MESHOP_COUNT_ONE

/** @brief True iff every roster row declares the nonzero element factor the
 *         guards divide by. */
constexpr bool mesh_op_elements_all_nonzero() {
  for (const MeshOpBoundsEntry &entry : MESHOP_BOUNDS)
    if (entry.bounds.elements == 0)
      return false;
  return true;
}
static_assert(mesh_op_elements_all_nonzero(),
              "a zero elements factor rejects every mesh outright");

/**
 * @brief Looks up an operator's declared growth factors by name.
 * @param name MeshOps operator name.
 * @return Pointer to the operator's row, or null when the name is not on the
 *         roster.
 */
inline const MeshOpBoundsEntry *find_mesh_op_bounds(const char *name) {
  for (const MeshOpBoundsEntry &entry : MESHOP_BOUNDS)
    if (std::strcmp(entry.name, name) == 0)
      return &entry;
  return nullptr;
}

} // namespace hs_wasm
