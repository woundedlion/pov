/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Host unit tests for the WASM boundary predicates
 * (targets/wasm/wasm_predicates.h) and the mesh-operator growth factors those
 * predicates are fed (targets/wasm/mesh_op_bounds.h). These checks gate untyped
 * integers crossing the embind boundary before they reach engine code that would
 * trap (clip bounds) or run unbounded (relax iterations). They compile only
 * under Emscripten inside wasm.cpp, so the pure predicates are extracted and
 * exercised here without the toolchain; the growth factors are measured against
 * the real operators.
 */
#pragma once

#include <algorithm>
#include <cstdint>
#include <iterator>
#include "core/mesh/solids.h" // pulls in the MeshOps operators the probes call
#include "targets/wasm/mesh_op_bounds.h"
#include "targets/wasm/wasm_predicates.h"
#include "tests/mesh_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace wasm_predicates_tests {

/**
 * @brief Exercises the near-pole LOD aggressiveness clamp, NaN included.
 */
inline void test_pole_lod_clamp() {
  HS_EXPECT_EQ(hs_wasm::clamp_pole_lod_aggressiveness(-1.0f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_pole_lod_aggressiveness(1.5f), 1.5f);
  HS_EXPECT_EQ(hs_wasm::clamp_pole_lod_aggressiveness(1e9f),
               hs_wasm::MAX_POLE_LOD_AGGRESSIVENESS);
  HS_EXPECT_EQ(hs_wasm::clamp_pole_lod_aggressiveness(NAN), 0.0f);
}

/**
 * @brief Exercises clip_bounds_valid across in-range, edge, and malformed bands.
 */
inline void test_clip_bounds() {
  constexpr int W = 96, H = 20;

  // A typical sub-canvas band is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(0, 48, 0, 10, W, H));
  // Full-canvas band (exclusive ends equal the extent) is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(0, W, 0, H, W, H));
  // Empty but ordered band (x0 == x1) is accepted.
  HS_EXPECT_TRUE(hs_wasm::clip_bounds_valid(5, 5, 3, 3, W, H));

  // Negative origin rejected (would feed ClipRegion modulo arithmetic).
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(-1, 10, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, -1, 10, W, H));
  // Inverted order rejected on each axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(10, 5, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, 10, 5, W, H));
  // Out of canvas rejected on each axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, W + 1, 0, 10, W, H));
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, 10, 0, H + 1, W, H));
  // A transposed (y-first) call that swaps the extents must fail the check
  // rather than clip the wrong axis.
  HS_EXPECT_TRUE(!hs_wasm::clip_bounds_valid(0, H, 0, W, W, H));
}

/**
 * @brief Exercises clamp_relax_iterations across negative, in-range, and over.
 */
inline void test_relax_clamp() {
  constexpr int MAX = 1000;

  // In-range counts pass through unchanged.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(0, MAX), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1, MAX), 1);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(MAX, MAX), MAX);
  // Negative floors at 0.
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1, MAX), 0);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(-1000000, MAX), 0);
  // Over-large counts clamp to the cap (relax(1e9) would freeze the thread).
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(MAX + 1, MAX), MAX);
  HS_EXPECT_EQ(hs_wasm::clamp_relax_iterations(1000000000, MAX), MAX);
}

/**
 * @brief Exercises the [0,1] operator-fraction range check and clamp.
 */
inline void test_unit_fraction_clamp() {
  // In-range fractions, including the boundaries, pass through unchanged.
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(0.0f));
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(0.5f));
  HS_EXPECT_TRUE(!hs_wasm::unit_fraction_out_of_range(1.0f));
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(0.0f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(0.5f), 0.5f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1.0f), 1.0f);

  // Out-of-range fractions are flagged and saturate, so truncate/bevel cannot
  // trip their always-on HS_CHECK and abort the module.
  HS_EXPECT_TRUE(hs_wasm::unit_fraction_out_of_range(-0.001f));
  HS_EXPECT_TRUE(hs_wasm::unit_fraction_out_of_range(1.001f));
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(-0.001f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1.001f), 1.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(-1e30f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_unit_fraction(1e30f), 1.0f);
}

/**
 * @brief Exercises the [0,1) operator-fraction range check and clamp.
 */
inline void test_half_open_fraction_clamp() {
  // 1 is outside a half-open domain, unlike the inclusive clamp's.
  HS_EXPECT_TRUE(!hs_wasm::half_open_fraction_out_of_range(0.0f));
  HS_EXPECT_TRUE(!hs_wasm::half_open_fraction_out_of_range(0.5f));
  HS_EXPECT_TRUE(hs_wasm::half_open_fraction_out_of_range(1.0f));
  HS_EXPECT_TRUE(hs_wasm::half_open_fraction_out_of_range(-0.001f));
  HS_EXPECT_TRUE(hs_wasm::half_open_fraction_out_of_range(1e30f));

  HS_EXPECT_EQ(hs_wasm::clamp_half_open_fraction(0.0f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_half_open_fraction(0.5f), 0.5f);
  HS_EXPECT_EQ(hs_wasm::clamp_half_open_fraction(-0.001f), 0.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_half_open_fraction(-1e30f), 0.0f);

  // The saturated value is what chamfer/snub/expand assert against: strictly
  // below 1, and the largest float that is, so the clamp loses no usable range.
  HS_EXPECT_TRUE(hs_wasm::clamp_half_open_fraction(1.0f) < 1.0f);
  HS_EXPECT_TRUE(hs_wasm::clamp_half_open_fraction(1e30f) < 1.0f);
  HS_EXPECT_EQ(hs_wasm::clamp_half_open_fraction(1.0f),
               hs_wasm::LARGEST_FRACTION_BELOW_ONE);
  HS_EXPECT_TRUE(hs_wasm::LARGEST_FRACTION_BELOW_ONE > 0.999999f);
  HS_EXPECT_EQ(
      hs_wasm::clamp_half_open_fraction(hs_wasm::LARGEST_FRACTION_BELOW_ONE),
      hs_wasm::LARGEST_FRACTION_BELOW_ONE);
}

/**
 * @brief Exercises the per-operator tooling mesh-size ceiling.
 */
inline void test_tooling_mesh_ceiling() {
  constexpr size_t MAX = 65536;

  // A mesh at the ceiling on every count is still accepted.
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(0, 0, 0, MAX));
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(12, 20, 60, MAX));
  HS_EXPECT_TRUE(!hs_wasm::tooling_mesh_over_ceiling(MAX, MAX, MAX, MAX));

  // Each count is checked independently.
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(MAX + 1, 20, 60, MAX));
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(12, MAX + 1, 60, MAX));
  HS_EXPECT_TRUE(hs_wasm::tooling_mesh_over_ceiling(12, 20, MAX + 1, MAX));
}

/**
 * @brief Exercises the per-operator expansion ceiling.
 */
inline void test_mesh_op_expansion_ceiling() {
  constexpr size_t MAX = 65535; // the engine's 16-bit connectivity range

  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(8, 6, 24), 24u);
  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(70000, 6, 24), 70000u);
  HS_EXPECT_EQ(hs_wasm::mesh_largest_element_count(8, 70000, 24), 70000u);

  // A cube is far inside every operator's bound.
  for (size_t e :
       {size_t(1), size_t(2), size_t(3), size_t(4), size_t(5), size_t(6)}) {
    HS_EXPECT_TRUE(!hs_wasm::mesh_op_expansion_over_ceiling(8, 6, 24, e, MAX));
  }

  // A non-expanding operator (dual, relax) admits a mesh right at the ceiling.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX, 1, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX + 1, 1, MAX));

  // The admissible input scales down by the expansion: the same mesh passes at
  // 3x (kis) and fails at 6x (meta, bevel).
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 3, 3, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 3, 6, MAX));
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 6, 6, MAX));

  // The bound is per-stage, so it bites one element past MAX / expansion.
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 5, 5, MAX));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_expansion_over_ceiling(1, 1, MAX / 5 + 1, 5, MAX));

  // Any of the three counts can be the binding one.
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(30000, 6, 24, 4, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(8, 30000, 24, 4, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(8, 6, 30000, 4, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_expansion_over_ceiling(8, 6, 24, 0, MAX));
}

/**
 * @brief Exercises the tooling-arena room check for an operator's output.
 */
inline void test_mesh_op_arena_room() {
  constexpr size_t CAP = 8 * 1024 * 1024;
  constexpr size_t PER = 16;

  // A cube's output is nothing next to an empty 8 MB arena.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_output_over_arena(8, 6, 24, 6, PER, 0, CAP));
  // Priced at expansion * bytes_per_element, a whole-arena mesh still fits.
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_output_over_arena(1, 1, CAP / PER, 1, PER, 0, CAP));
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_output_over_arena(1, 1, CAP / PER + 1, 1, PER, 0, CAP));
  // The same mesh no longer fits once the arena is mostly spent.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_output_over_arena(1, 1, 1024, 1, PER,
                                                     CAP - 1024 * PER, CAP));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_output_over_arena(1, 1, 1024, 1, PER,
                                                    CAP - 1024 * PER + 1, CAP));
  // Expansion scales the price, so a 6x operator is rejected where 1x fits.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_output_over_arena(1, 1, 1024, 1, PER,
                                                     CAP - 2048 * PER, CAP));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_output_over_arena(1, 1, 1024, 6, PER,
                                                    CAP - 2048 * PER, CAP));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_output_over_arena(8, 6, 24, 0, PER, 0, CAP));
  // A full or unbound (capacity 0) arena rejects everything.
  HS_EXPECT_TRUE(
      hs_wasm::mesh_op_output_over_arena(8, 6, 24, 1, PER, CAP, CAP));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_output_over_arena(8, 6, 24, 1, PER, 0, 0));
}

/**
 * @brief Exercises the mesh measurements the face-degree guard reads.
 */
inline void test_mesh_degree_measurements() {
  // A cube: six quads, eight vertices of valence 3.
  const uint8_t cube_counts[6] = {4, 4, 4, 4, 4, 4};
  const uint16_t cube_faces[24] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 5, 4,
                                   2, 3, 7, 6, 0, 3, 7, 4, 1, 2, 6, 5};
  uint32_t incidence[8] = {};

  HS_EXPECT_EQ(hs_wasm::mesh_max_face_degree(cube_counts, 6), 4u);
  HS_EXPECT_EQ(hs_wasm::mesh_max_vertex_valence(cube_faces, 24, incidence, 8),
               3u);

  // An empty mesh measures zero on both axes.
  HS_EXPECT_EQ(hs_wasm::mesh_max_face_degree(cube_counts, 0), 0u);
  HS_EXPECT_EQ(hs_wasm::mesh_max_vertex_valence(cube_faces, 0, incidence, 8),
               0u);

  // A mixed-degree list reports the widest face.
  const uint8_t mixed[4] = {3, 200, 6, 3};
  HS_EXPECT_EQ(hs_wasm::mesh_max_face_degree(mixed, 4), 200u);

  // The scan reports the single hottest vertex, not an average, and re-clears
  // its counters so a reused buffer cannot inflate the next measurement.
  const uint16_t fan[7] = {0, 0, 0, 0, 0, 1, 2};
  HS_EXPECT_EQ(hs_wasm::mesh_max_vertex_valence(fan, 7, incidence, 3), 5u);
  HS_EXPECT_EQ(hs_wasm::mesh_max_vertex_valence(fan, 7, incidence, 3), 5u);

  // An index past the vertex count is skipped rather than read out of bounds.
  const uint16_t stray[3] = {0, 1, 9};
  HS_EXPECT_EQ(hs_wasm::mesh_max_vertex_valence(stray, 3, incidence, 2), 1u);
}

/**
 * @brief Exercises the per-operator face-degree overflow check.
 */
inline void test_mesh_op_face_degree() {
  constexpr size_t MAX = 255; // PolyMesh's uint8_t per-face side count

  // A cube is inside every operator's degree bound.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(4, 3, 2, 2, MAX));
  // kis measures nothing and can never overflow: it emits only triangles.
  HS_EXPECT_TRUE(
      !hs_wasm::mesh_op_face_degree_overflows(255, 100000, 0, 0, MAX));

  // A vertex-orbit operator (dual, needle) turns valence into a side count.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(0, MAX, 0, 1, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(0, MAX + 1, 0, 1, MAX));
  // The kis x7 cube that trips narrow_face_count: valence 384 through dual.
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(3, 384, 0, 1, MAX));

  // A doubling operator (zip's valence, truncate's face degree) halves the room.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(0, 127, 0, 2, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(0, 128, 0, 2, MAX));
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(127, 3, 2, 1, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(128, 3, 2, 1, MAX));

  // The two axes are independent: a wide face passes an operator that only
  // widens valence, and vice versa.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(200, 3, 0, 2, MAX));
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(3, 200, 1, 0, MAX));
  // Hankin doubles both, so either axis alone can reject it.
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(200, 3, 2, 2, MAX));
  HS_EXPECT_TRUE(hs_wasm::mesh_op_face_degree_overflows(3, 200, 2, 2, MAX));

  // A factor of 1 on the face degree cannot fire: PolyMesh stores side counts
  // in a uint8_t, so the input's widest face is already within the ceiling.
  HS_EXPECT_TRUE(!hs_wasm::mesh_op_face_degree_overflows(MAX, 3, 1, 1, MAX));
}

// ---------------------------------------------------------------------------
// The declared growth factors against the real operators
// ---------------------------------------------------------------------------

/** @brief Arena backing a probe's input solid; outlives the operator run. */
inline uint8_t mesh_op_input_buf[64 * 1024];
/** @brief Arena receiving a probed operator's output mesh. */
inline uint8_t mesh_op_target_buf[256 * 1024];
/** @brief Arena holding a probed operator's intermediate stages. */
inline uint8_t mesh_op_temp_buf[256 * 1024];

/** @brief One roster operator bound to a concrete call. */
struct MeshOpProbe {
  const char *name; /**< Roster name, matched against MESHOP_BOUNDS. */
  PolyMesh (*run)(const PolyMesh &, Arena &, Arena &); /**< Operator call. */
};

// Fraction arguments avoid the topology-degenerate ends of each domain:
// truncate at 0.5 aliases to ambo, and the [0,1) operators trap at 1.
inline const MeshOpProbe MESH_OP_PROBES[] = {
    {"kis", [](const PolyMesh &m, Arena &a,
               Arena &b) { return MeshOps::kis(m, a, b); }},
    {"ambo", [](const PolyMesh &m, Arena &a,
                Arena &b) { return MeshOps::ambo(m, a, b); }},
    {"gyro", [](const PolyMesh &m, Arena &a,
                Arena &b) { return MeshOps::gyro(m, a, b); }},
    {"dual", [](const PolyMesh &m, Arena &a,
                Arena &b) { return MeshOps::dual(m, a, b); }},
    {"meta", [](const PolyMesh &m, Arena &a,
                Arena &b) { return MeshOps::meta(m, a, b); }},
    {"needle", [](const PolyMesh &m, Arena &a,
                  Arena &b) { return MeshOps::needle(m, a, b); }},
    {"zip", [](const PolyMesh &m, Arena &a,
               Arena &b) { return MeshOps::zip(m, a, b); }},
    {"truncate", [](const PolyMesh &m, Arena &a,
                    Arena &b) { return MeshOps::truncate(m, a, b, 0.25f); }},
    {"bevel", [](const PolyMesh &m, Arena &a,
                 Arena &b) { return MeshOps::bevel(m, a, b, 0.25f); }},
    {"chamfer", [](const PolyMesh &m, Arena &a,
                   Arena &b) { return MeshOps::chamfer(m, a, b, 0.5f); }},
    {"expand", [](const PolyMesh &m, Arena &a,
                  Arena &b) { return MeshOps::expand(m, a, b, 0.5f); }},
    {"relax", [](const PolyMesh &m, Arena &a,
                 Arena &b) { return MeshOps::relax(m, a, b, 2); }},
    {"hankin", [](const PolyMesh &m, Arena &a,
                  Arena &b) { return MeshOps::hankin(m, a, b, 0.5f); }},
    {"snub", [](const PolyMesh &m, Arena &a, Arena &b) {
       return MeshOps::snub(m, a, b, 0.5f, 0.0f);
     }}};

// Widest face an operator emits at a side count fixed by its construction
// rather than by the input (kis triangles, gyro pentagons, chamfer hexagons).
// Those cannot scale into the 8-bit side count, so no growth factor covers them.
inline constexpr size_t FIXED_EMITTER_DEGREE = 6;

/**
 * @brief Runs one operator and checks its real growth against its declared
 *        factors.
 * @param probe Operator to run.
 * @param in Input mesh, allocated from an arena neither @p target nor @p temp.
 * @param target Arena receiving the output mesh; reset here.
 * @param temp Arena for the operator's intermediates; reset here.
 * @details A declared factor below the real expansion is the dangerous
 *          direction: the boundary guard then admits a mesh whose operator
 *          reaches build_half_edge_mesh's or narrow_face_count's always-on
 *          HS_CHECK and takes the module down.
 */
inline void check_mesh_op_growth(const MeshOpProbe &probe, const PolyMesh &in,
                                 Arena &target, Arena &temp) {
  const hs_wasm::MeshOpBoundsEntry *row =
      hs_wasm::find_mesh_op_bounds(probe.name);
  HS_EXPECT_TRUE(row != nullptr);
  if (row == nullptr)
    return;

  uint32_t incidence[64] = {};
  HS_EXPECT_LE(in.vertices.size(), std::size(incidence));
  const size_t in_elements = hs_wasm::mesh_largest_element_count(
      in.vertices.size(), in.get_face_counts_size(), in.get_faces_size());
  const size_t in_degree = hs_wasm::mesh_max_face_degree(
      in.get_face_counts_data(), in.get_face_counts_size());
  const size_t in_valence = hs_wasm::mesh_max_vertex_valence(
      in.get_faces_data(), in.get_faces_size(), incidence, in.vertices.size());

  target.reset();
  temp.reset();
  const PolyMesh out = probe.run(in, target, temp);

  const size_t out_elements = hs_wasm::mesh_largest_element_count(
      out.vertices.size(), out.get_face_counts_size(), out.get_faces_size());
  const size_t out_degree = hs_wasm::mesh_max_face_degree(
      out.get_face_counts_data(), out.get_face_counts_size());

  HS_EXPECT_GT(out_elements, 0u);
  HS_EXPECT_LE(out_elements, row->bounds.elements * in_elements);
  HS_EXPECT_LE(out_degree, std::max({FIXED_EMITTER_DEGREE,
                                     row->bounds.face_degree * in_degree,
                                     row->bounds.valence * in_valence}));
}

/**
 * @brief Runs every roster operator against one solid.
 * @tparam Solid Solids::* descriptor to build the input from.
 * @param input Arena the input solid is rebuilt in.
 * @param target Arena receiving each operator's output.
 * @param temp Arena for each operator's intermediates.
 */
template <typename Solid>
inline void probe_solid_growth(Arena &input, Arena &target, Arena &temp) {
  input.reset();
  PolyMesh in;
  build_solid<Solid>(in, input);
  for (const MeshOpProbe &probe : MESH_OP_PROBES)
    check_mesh_op_growth(probe, in, target, temp);
}

/**
 * @brief Measures every roster operator's growth on solids of differing face
 *        degree and vertex valence.
 */
inline void test_mesh_op_growth_factors() {
  Arena input(mesh_op_input_buf, sizeof(mesh_op_input_buf));
  Arena target(mesh_op_target_buf, sizeof(mesh_op_target_buf));
  Arena temp(mesh_op_temp_buf, sizeof(mesh_op_temp_buf));

  // An operator on the roster with no probe here would ship an unmeasured
  // factor.
  HS_EXPECT_EQ(std::size(MESH_OP_PROBES), hs_wasm::MESHOP_BOUNDS_COUNT);

  // Face degree and vertex valence differ across these five (3/3, 4/3, 3/4,
  // 3/5, 5/3), so a factor keyed to the wrong measurement cannot pass on all.
  probe_solid_growth<Solids::Tetrahedron>(input, target, temp);
  probe_solid_growth<Solids::Cube>(input, target, temp);
  probe_solid_growth<Solids::Octahedron>(input, target, temp);
  probe_solid_growth<Solids::Icosahedron>(input, target, temp);
  probe_solid_growth<Solids::Dodecahedron>(input, target, temp);
}

/**
 * @brief Exercises the Hankin contact-angle domain check.
 */
inline void test_hankin_angle_domain() {
  constexpr float MAX = 1.5707963f; // pi/2
  constexpr float FIFTY_FOUR_DEG = 0.9424778f;

  // In-domain angles, including both endpoints, are accepted.
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(0.0f, MAX));
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(FIFTY_FOUR_DEG, MAX));
  HS_EXPECT_TRUE(!hs_wasm::hankin_angle_out_of_range(MAX, MAX));

  // Negative and past-pi/2 angles are rejected rather than mirroring an
  // in-domain pattern.
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(-0.001f, MAX));
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(MAX + 0.01f, MAX));
  HS_EXPECT_TRUE(hs_wasm::hankin_angle_out_of_range(12.0f, MAX));
}

/**
 * @brief Module entry point: runs the boundary-predicate checks.
 * @return The module's failure count.
 */
inline int run_wasm_predicates_tests() {
  hs_test::ModuleFixture fixture("wasm_predicates");
  test_pole_lod_clamp();
  test_clip_bounds();
  test_relax_clamp();
  test_unit_fraction_clamp();
  test_half_open_fraction_clamp();
  test_tooling_mesh_ceiling();
  test_mesh_op_expansion_ceiling();
  test_mesh_op_arena_room();
  test_mesh_degree_measurements();
  test_mesh_op_face_degree();
  test_mesh_op_growth_factors();
  test_hankin_angle_domain();
  return fixture.result();
}

} // namespace wasm_predicates_tests
} // namespace hs_test
