/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * @file mesh_op_bounds.h
 * @brief Pure (no-Emscripten) mesh-operator roster, growth factors and
 *        byte-per-element budgets.
 *
 * Every MeshOps boundary guard in wasm.cpp prices a JS-driven operator chain
 * against these numbers, so one smaller than an operator's real expansion or
 * footprint lets a chain reach an engine trap the guard exists to keep out of
 * reach. They carry no Emscripten dependency, so they live here and are measured
 * against the real operators without the toolchain — see
 * tests/test_wasm_predicates.h. wasm.cpp keeps the logging/embind shell on top.
 */
#pragma once

#include <cstddef>
#include <cstring>  // std::strcmp — roster lookup by operator name
#include <iterator> // std::size — roster table length

namespace hs_wasm {

/**
 * @brief How far one mesh operator grows its input, for the boundary guards.
 * @details Every field is a multiple of an input measurement that some stage of
 *          the operator reaches; see MESHOP_LIST for the per-operator values and
 *          where they come from. `elements` must be nonzero. A zero degree or
 *          valence means the operator emits nothing of that kind; a zero
 *          `valence` additionally skips the valence scan, the one measurement
 *          that costs a pass over the flat index list.
 */
struct MeshOpBounds {
  size_t elements;    /**< Multiple of the largest input element count. */
  size_t face_degree; /**< Multiple of the widest input face's side count. */
  size_t valence;     /**< Multiple of the highest input vertex valence. */
};

/**
 * @brief Single source of truth for the Conway/Goldberg operator roster and each
 *        operator's growth factors.
 * @param _OP0  Macro applied to each zero-argument operator.
 * @param _OP1U Macro applied to each [0,1]-fraction operator.
 * @param _OP1H Macro applied to each [0,1)-fraction operator.
 * @details Expanded three times: with MESHOP_0/MESHOP_1U/MESHOP_1H to generate
 *          the wrapper methods in wasm.cpp, with MESHOP_BIND to generate the
 *          embind .function() bindings there, and with MESHOP_BOUNDS_ENTRY below
 *          to publish the factors as data, so an operator cannot be added to one
 *          site and silently left unreachable or unmeasured at another.
 *          Each fraction operator's macro matches the domain its always-on
 *          engine trap asserts: truncate and bevel accept 1 and use _OP1U;
 *          chamfer and expand assert t < 1 and use _OP1H. relax, hankin,
 *          and snub have bespoke signatures/validation (snub takes two float
 *          controls, and clamps its inset to [0,1) for the same trap), so their
 *          wrapper methods are hand-written; their names live in
 *          MESHOP_IRREGULAR_LIST below and their bindings expand from it, so a
 *          new irregular op is bound the moment it joins the list, and their
 *          factors are the RELAX_BOUNDS/HANKIN_BOUNDS/SNUB_BOUNDS constants.
 *
 *          The trailing arguments are the operator's MeshOpBounds, in order.
 *
 *          `elements` is the largest multiple of the input's biggest element
 *          count that any of its intermediate or output stages reaches. Read off
 *          the output bindings in core/mesh/conway.h with I as the input flat
 *          index count — dual I, ambo 2I, kis/truncate 3I, expand/chamfer 4I,
 *          snub 5I — and multiplied through for the compositions: gyro = d(snub)
 *          5I, needle = k(d) 3I, zip = d(k) 3I, meta = k(d(a)) 6I, bevel = t(a)
 *          6I. Every stage's vertex and face count stays under its operator's
 *          index expansion, so one factor bounds all three counts.
 *
 *          `degree` and `valence` are the multiples that reach
 *          narrow_face_count, from the three emitters: emit_shrunk_face and
 *          emit_primary_faces widen a face to verts_per_side x its side count
 *          (1 for ambo/expand/chamfer/snub, 2 for truncate),
 *          emit_vertex_orbit_faces emits one face per vertex at its valence, and
 *          Hankin doubles both (star faces 2x a face's sides, rosettes 2x a
 *          valence). kis emits only triangles, so it needs neither measurement.
 *          Compositions inherit their widest stage: needle's dual stage is
 *          valence-wide, zip's is 2x valence (its kis stage doubles valence
 *          first), meta's and gyro's ambo/snub stages are 1x both, and bevel's
 *          truncate stage doubles what its ambo stage already widened.
 *
 *          Faces an operator emits at a fixed side count regardless of its input
 *          (kis triangles, gyro pentagons, chamfer hexagons) cannot scale into
 *          the 8-bit side count, so no factor covers them.
 */
// clang-format off
#define MESHOP_LIST(_OP0, _OP1U, _OP1H)                                         \
  _OP0(kis, 3, 0, 0) _OP0(ambo, 2, 1, 1) _OP0(gyro, 5, 1, 1)                    \
  _OP0(dual, 1, 0, 1) _OP0(meta, 6, 1, 1) _OP0(needle, 3, 0, 1)                 \
  _OP0(zip, 3, 1, 2)                                                            \
  _OP1U(truncate, 3, 2, 1) _OP1U(bevel, 6, 2, 2)                                \
  _OP1H(chamfer, 4, 1, 0) _OP1H(expand, 4, 1, 1)
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
 *          index scratch live in one arena. Measured against the real operators
 *          by tests/test_wasm_predicates.h; wasm.cpp ties it to the tooling
 *          scratch arena size.
 */
inline constexpr size_t TOOLING_BYTES_PER_MESH_ELEMENT = 64;

/**
 * @brief Tooling-arena bytes a finalized mesh retains per element.
 * @details One Vector of vertex, one uint8_t of side count, one uint16_t of
 *          index and the uint16_t topology code classifyFaces() later binds
 *          into the same arena, with slack for the four blocks' alignment
 *          padding. Measured against the real operators by
 *          tests/test_wasm_predicates.h.
 */
inline constexpr size_t TOOLING_ARENA_BYTES_PER_MESH_ELEMENT = 20;

/** @brief One operator's declared growth factors, by name. */
struct MeshOpBoundsEntry {
  const char *name;    /**< MeshOps operator name. */
  MeshOpBounds bounds; /**< Factors the boundary guard prices it against. */
};

#define MESHOP_BOUNDS_ENTRY(name, elements, degree, valence)                   \
  {#name, {(elements), (degree), (valence)}},

/**
 * @brief Every operator's declared growth factors as data.
 * @details The roster's regular operators expand from MESHOP_LIST; the irregular
 *          three carry the constants their hand-written wrappers pass. Exists so
 *          a host test can measure the real operators against the same factors
 *          the Emscripten-only guards use.
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
