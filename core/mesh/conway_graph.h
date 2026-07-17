/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cstdint>
#include <string_view>

#include "mesh/solids.h"

/**
 * @brief Edge graph for animated Conway-operator transitions between the
 * simple-registry solids (docs/conway_morph_spec.md).
 * @details A node is a simple-registry solid; an edge is one animated operator
 * sweep between two parameter values on one seed. The table plus the pure walk
 * and seed-reconciliation helpers here are consumed by Animation::ConwayMorph
 * and HankinSolids; everything is constexpr and unit-testable with no effect
 * or animation dependency.
 */
namespace ConwayGraph {

// Node ids are simple-registry indices (Solids::simple_registry order).
inline constexpr uint8_t TETRAHEDRON = 0;
inline constexpr uint8_t CUBE = 1;
inline constexpr uint8_t OCTAHEDRON = 2;
inline constexpr uint8_t DODECAHEDRON = 3;
inline constexpr uint8_t ICOSAHEDRON = 4;
inline constexpr uint8_t TRUNCATED_TETRAHEDRON = 5;
inline constexpr uint8_t CUBOCTAHEDRON = 6;
inline constexpr uint8_t TRUNCATED_CUBE = 7;
inline constexpr uint8_t TRUNCATED_OCTAHEDRON = 8;
inline constexpr uint8_t RHOMBICUBOCTAHEDRON = 9;
inline constexpr uint8_t TRUNCATED_CUBOCTAHEDRON = 10;
inline constexpr uint8_t SNUB_CUBE = 11;
inline constexpr uint8_t ICOSIDODECAHEDRON = 12;
inline constexpr uint8_t TRUNCATED_DODECAHEDRON = 13;
inline constexpr uint8_t TRUNCATED_ICOSAHEDRON = 14;
inline constexpr uint8_t RHOMBICOSIDODECAHEDRON = 15;
inline constexpr uint8_t TRUNCATED_ICOSIDODECAHEDRON = 16;
inline constexpr uint8_t SNUB_DODECAHEDRON = 17;

inline constexpr int NUM_NODES = 18;

// Pin the node ids to the registry so a registry reorder fails to compile.
static_assert(std::string_view(Solids::simple_registry[TETRAHEDRON].name) ==
              "tetrahedron");
static_assert(std::string_view(Solids::simple_registry[CUBE].name) == "cube");
static_assert(std::string_view(Solids::simple_registry[OCTAHEDRON].name) ==
              "octahedron");
static_assert(std::string_view(Solids::simple_registry[DODECAHEDRON].name) ==
              "dodecahedron");
static_assert(std::string_view(Solids::simple_registry[ICOSAHEDRON].name) ==
              "icosahedron");
static_assert(
    std::string_view(Solids::simple_registry[TRUNCATED_TETRAHEDRON].name) ==
    "truncatedTetrahedron");
static_assert(std::string_view(Solids::simple_registry[CUBOCTAHEDRON].name) ==
              "cuboctahedron");
static_assert(std::string_view(Solids::simple_registry[TRUNCATED_CUBE].name) ==
              "truncatedCube");
static_assert(
    std::string_view(Solids::simple_registry[TRUNCATED_OCTAHEDRON].name) ==
    "truncatedOctahedron");
static_assert(
    std::string_view(Solids::simple_registry[RHOMBICUBOCTAHEDRON].name) ==
    "rhombicuboctahedron");
static_assert(
    std::string_view(Solids::simple_registry[TRUNCATED_CUBOCTAHEDRON].name) ==
    "truncatedCuboctahedron");
static_assert(std::string_view(Solids::simple_registry[SNUB_CUBE].name) ==
              "snubCube");
static_assert(
    std::string_view(Solids::simple_registry[ICOSIDODECAHEDRON].name) ==
    "icosidodecahedron");
static_assert(
    std::string_view(Solids::simple_registry[TRUNCATED_DODECAHEDRON].name) ==
    "truncatedDodecahedron");
static_assert(
    std::string_view(Solids::simple_registry[TRUNCATED_ICOSAHEDRON].name) ==
    "truncatedIcosahedron");
static_assert(
    std::string_view(Solids::simple_registry[RHOMBICOSIDODECAHEDRON].name) ==
    "rhombicosidodecahedron");
static_assert(std::string_view(
                  Solids::simple_registry[TRUNCATED_ICOSIDODECAHEDRON].name) ==
              "truncatedIcosidodecahedron");
static_assert(
    std::string_view(Solids::simple_registry[SNUB_DODECAHEDRON].name) ==
    "snubDodecahedron");

/**
 * @brief Parameterized Conway operator a graph edge sweeps.
 * @details Chamfer is deliberately unused: no simple-registry endpoint is a
 * chamfered form.
 */
enum class MorphOp : uint8_t { TRUNCATE, EXPAND, SNUB };

/**
 * @brief Reseed primitive tabled on an edge (spec section 2.1).
 * @details ADOPT on a bridge row replaces the family seed with the arrived
 * solid at forward completion; ADOPT on an ambo-chain row (seed_solid is
 * cuboctahedron/icosidodecahedron) means the leg seed is derived as
 * ambo(held platonic seed) at construction (see SeedFix::DERIVE_AMBO).
 * DUAL_SWAP marks the ambo-crossover rows; the actual swap is decided at leg
 * start by seed_fix_at_start against the held seed identity.
 */
enum class Reseed : uint8_t { NONE, ADOPT, DUAL_SWAP };

/** Sweep parameter clamp: legs run inside [T_EPS, t_end], and truncate legs
 * additionally stay below 0.5 - T_EPS_AMBO (the ambo short-circuit changes
 * emission order and face count). */
inline constexpr float T_EPS = 0.02f;
/** Truncate clamp at the ambo (t = 0.5) end. Tighter than T_EPS: near 0.5 no
 * face degenerates (only the residual seed-edge segments shrink), and the
 * clean-swap gap this leaves is sub-pixel at H = 144 where T_EPS's would be a
 * visible ~2 px boundary jump plus a gradient rescale at every leg boundary
 * touching the ambo form. */
inline constexpr float T_EPS_AMBO = 0.005f;
/** Snub clamp at the jitterbug bridge's octahedron end (t = 0.5, where the 12
 * vertices merge pairwise onto the octahedron's 6): the leg stops where the
 * collapsing edge measures 0.02 chord — T_EPS-sized, ~1-2 px at H = 144 —
 * then clean-swaps to the held octahedron. */
inline constexpr float T_EPS_JITTERBUG = 0.5104592f;

/** Operator-sweep frames per leg. */
inline constexpr int SWEEP_FRAMES = 48;
/** Settle (relax-slerp) frames appended to a settling leg; 0 for the rest. */
inline constexpr int SETTLE_FRAMES = 12;

/** Snub-cube twist (matches the registry chain). */
inline constexpr float SNUB_CUBE_TWIST = 0.28f;
/** Snub-dodecahedron cosmetic sweep twist; registry uses 0 (tuned from
 * renders). */
inline constexpr float SNUB_DODECAHEDRON_TWIST = 0.0f;
/** Tetra -> icosa bridge snub twist; relax canonicalizes any value (tuned from
 * renders: -0.40 cuts the settle rotation from 23.4 to 17.3 degrees). */
inline constexpr float SNUB_BRIDGE_TWIST = -0.40f;

/** Jitterbug icosa point: snub(tetrahedron, t, twist) at these values is the
 * exact regular icosahedron — all 30 edges equal with no relax (double-refined
 * t = 0.7299092432622736, twist = -0.3881395153701886). */
inline constexpr float T_JITTERBUG_ICOSA = 0.72990924f;
inline constexpr float TWIST_JITTERBUG_ICOSA = -0.38813952f;
/** Jitterbug octa twist: snub(tetrahedron, 0.5, -pi/3) merges its 12 vertices
 * pairwise onto the octahedron's 6 (the jitterbug closure). */
inline constexpr float TWIST_JITTERBUG_OCTA = -PI_F / 3.0f;

/** Truncation parameter of the truncated cube. */
inline constexpr float T_TRUNC_CUBE = 1.0f / (2.0f + SQRT2);
/** Truncation parameter of the truncated tetra/octa/icosahedron. */
inline constexpr float T_TRUNC_THIRD = 1.0f / 3.0f;

/**
 * @brief One bidirectional graph edge: a single operator sweep on one seed.
 * @details Nodes and seed_solid are simple-registry indices. `settle` means
 * the to_node end is the relax-canonical form (relax(50) at that endpoint).
 * `bridge` marks a symmetry-family crossing, used by the walk weighting.
 */
struct EdgeSpec {
  uint8_t from_node;  /**< Node at the t_from end. */
  uint8_t to_node;    /**< Node at the t_to end. */
  uint8_t seed_solid; /**< Solid the op sweeps on (leg seed identity). */
  MorphOp op;         /**< Operator whose parameter the leg sweeps. */
  float t_from; /**< Parameter at from_node (untabled clamp applies at runtime).
                 */
  float t_to;   /**< Parameter at to_node. */
  float twist_from; /**< Snub twist at from_node (snub legs only). */
  float twist_to;   /**< Snub twist at to_node (snub legs only). */
  bool settle;      /**< to_node end is the relax(50) canonical form. */
  Reseed reseed;    /**< Reseed primitive tabled for this row. */
  bool bridge;      /**< Crosses symmetry families (walk weighting). */
};

/** The 23 edges (spec section 3); every simple-registry solid is a node. */
inline constexpr EdgeSpec EDGES[] = {
    // Octahedral family (seed: cube, octahedron)
    {CUBE, TRUNCATED_CUBE, CUBE, MorphOp::TRUNCATE, 0.0f, T_TRUNC_CUBE, 0.0f,
     0.0f, false, Reseed::NONE, false},
    {CUBE, CUBOCTAHEDRON, CUBE, MorphOp::TRUNCATE, 0.0f, 0.5f, 0.0f, 0.0f,
     false, Reseed::NONE, false},
    {TRUNCATED_CUBE, CUBOCTAHEDRON, CUBE, MorphOp::TRUNCATE, T_TRUNC_CUBE, 0.5f,
     0.0f, 0.0f, false, Reseed::NONE, false},
    {CUBE, RHOMBICUBOCTAHEDRON, CUBE, MorphOp::EXPAND, 0.0f,
     MeshOps::EXPAND_DEFAULT_T, 0.0f, 0.0f, false, Reseed::NONE, false},
    {CUBE, SNUB_CUBE, CUBE, MorphOp::SNUB, 0.0f, T_SNUB_CUBE, 0.0f,
     SNUB_CUBE_TWIST, true, Reseed::NONE, false},
    {CUBOCTAHEDRON, TRUNCATED_CUBOCTAHEDRON, CUBOCTAHEDRON, MorphOp::TRUNCATE,
     0.0f, T_TRUNC_CUBE, 0.0f, 0.0f, true, Reseed::ADOPT, false},
    {OCTAHEDRON, TRUNCATED_OCTAHEDRON, OCTAHEDRON, MorphOp::TRUNCATE, 0.0f,
     T_TRUNC_THIRD, 0.0f, 0.0f, false, Reseed::NONE, false},
    {TRUNCATED_OCTAHEDRON, CUBOCTAHEDRON, OCTAHEDRON, MorphOp::TRUNCATE,
     T_TRUNC_THIRD, 0.5f, 0.0f, 0.0f, false, Reseed::NONE, false},
    {OCTAHEDRON, CUBOCTAHEDRON, OCTAHEDRON, MorphOp::TRUNCATE, 0.0f, 0.5f, 0.0f,
     0.0f, false, Reseed::DUAL_SWAP, false},

    // Icosahedral family (seed: dodecahedron, icosahedron)
    {DODECAHEDRON, TRUNCATED_DODECAHEDRON, DODECAHEDRON, MorphOp::TRUNCATE,
     0.0f, T_TRUNC_ICOS, 0.0f, 0.0f, false, Reseed::NONE, false},
    {DODECAHEDRON, ICOSIDODECAHEDRON, DODECAHEDRON, MorphOp::TRUNCATE, 0.0f,
     0.5f, 0.0f, 0.0f, false, Reseed::NONE, false},
    {TRUNCATED_DODECAHEDRON, ICOSIDODECAHEDRON, DODECAHEDRON, MorphOp::TRUNCATE,
     T_TRUNC_ICOS, 0.5f, 0.0f, 0.0f, false, Reseed::NONE, false},
    {ICOSAHEDRON, TRUNCATED_ICOSAHEDRON, ICOSAHEDRON, MorphOp::TRUNCATE, 0.0f,
     T_TRUNC_THIRD, 0.0f, 0.0f, false, Reseed::NONE, false},
    {TRUNCATED_ICOSAHEDRON, ICOSIDODECAHEDRON, ICOSAHEDRON, MorphOp::TRUNCATE,
     T_TRUNC_THIRD, 0.5f, 0.0f, 0.0f, false, Reseed::NONE, false},
    {ICOSAHEDRON, ICOSIDODECAHEDRON, ICOSAHEDRON, MorphOp::TRUNCATE, 0.0f, 0.5f,
     0.0f, 0.0f, false, Reseed::DUAL_SWAP, false},
    {DODECAHEDRON, RHOMBICOSIDODECAHEDRON, DODECAHEDRON, MorphOp::EXPAND, 0.0f,
     MeshOps::EXPAND_DEFAULT_T, 0.0f, 0.0f, true, Reseed::NONE, false},
    {DODECAHEDRON, SNUB_DODECAHEDRON, DODECAHEDRON, MorphOp::SNUB, 0.0f, 0.5f,
     0.0f, SNUB_DODECAHEDRON_TWIST, true, Reseed::NONE, false},
    {ICOSIDODECAHEDRON, TRUNCATED_ICOSIDODECAHEDRON, ICOSIDODECAHEDRON,
     MorphOp::TRUNCATE, 0.0f, T_TRUNC_ICOS, 0.0f, 0.0f, true, Reseed::ADOPT,
     false},

    // Tetrahedral family and bridges
    {TETRAHEDRON, TRUNCATED_TETRAHEDRON, TETRAHEDRON, MorphOp::TRUNCATE, 0.0f,
     T_TRUNC_THIRD, 0.0f, 0.0f, false, Reseed::NONE, false},
    {TETRAHEDRON, OCTAHEDRON, TETRAHEDRON, MorphOp::TRUNCATE, 0.0f, 0.5f, 0.0f,
     0.0f, false, Reseed::ADOPT, true},
    {TRUNCATED_TETRAHEDRON, OCTAHEDRON, TETRAHEDRON, MorphOp::TRUNCATE,
     T_TRUNC_THIRD, 0.5f, 0.0f, 0.0f, false, Reseed::ADOPT, true},
    {TETRAHEDRON, ICOSAHEDRON, TETRAHEDRON, MorphOp::SNUB, 0.0f, 0.5f, 0.0f,
     SNUB_BRIDGE_TWIST, true, Reseed::ADOPT, true},
    // Jitterbug bridge: a snub sweep between two interior parameter points,
    // the exact regular icosahedron and the pairwise-merged octahedron.
    {ICOSAHEDRON, OCTAHEDRON, TETRAHEDRON, MorphOp::SNUB, T_JITTERBUG_ICOSA,
     0.5f, TWIST_JITTERBUG_ICOSA, TWIST_JITTERBUG_OCTA, false, Reseed::ADOPT,
     true},
};

inline constexpr int NUM_EDGES = static_cast<int>(std::size(EDGES));
static_assert(NUM_EDGES == 23);

/**
 * @brief Whether an edge is the icosahedron <-> octahedron jitterbug bridge:
 * a snub sweep whose t = 0.5 end collapses onto the octahedron by pairwise
 * vertex merge. The T_EPS_JITTERBUG clamp and the ambo endpoint swap key on
 * it.
 */
constexpr bool is_jitterbug_edge(const EdgeSpec &e) {
  return e.op == MorphOp::SNUB && e.to_node == OCTAHEDRON;
}

/** Largest node degree in the table (cuboctahedron, icosidodecahedron,
 * octahedron). */
inline constexpr int MAX_DEGREE = 5;

/**
 * @brief Whether an edge is incident to a node.
 * @param edge Index into EDGES.
 * @param node Simple-registry node id.
 * @return True iff the node is either endpoint.
 */
constexpr bool edge_touches(int edge, int node) {
  return EDGES[edge].from_node == node || EDGES[edge].to_node == node;
}

/**
 * @brief The opposite endpoint of an edge.
 * @param edge Index into EDGES; must be incident to @p node.
 * @param node Simple-registry node id at one endpoint.
 * @return The node id at the other endpoint.
 */
constexpr int edge_other_end(int edge, int node) {
  return EDGES[edge].from_node == node ? EDGES[edge].to_node
                                       : EDGES[edge].from_node;
}

/**
 * @brief Number of edges incident to a node.
 * @param node Simple-registry node id.
 * @return The node's degree in the graph.
 */
constexpr int node_degree(int node) {
  int n = 0;
  for (int e = 0; e < NUM_EDGES; ++e)
    if (edge_touches(e, node))
      ++n;
  return n;
}

/**
 * @brief Collects the edges incident to a node, in table order.
 * @param node Simple-registry node id.
 * @param out Receives edge indices; must hold at least MAX_DEGREE entries.
 * @return Number of edges written.
 */
constexpr int edges_from(int node, uint8_t *out) {
  int n = 0;
  for (int e = 0; e < NUM_EDGES; ++e)
    if (edge_touches(e, node))
      out[n++] = static_cast<uint8_t>(e);
  return n;
}

// Degree pins: every node is covered, the six out-and-back leaves are leaves.
static_assert(node_degree(TETRAHEDRON) == 3);
static_assert(node_degree(CUBE) == 4 && node_degree(OCTAHEDRON) == 5);
static_assert(node_degree(DODECAHEDRON) == 4 && node_degree(ICOSAHEDRON) == 4);
static_assert(node_degree(CUBOCTAHEDRON) == 5 &&
              node_degree(ICOSIDODECAHEDRON) == 5);
static_assert(node_degree(RHOMBICUBOCTAHEDRON) == 1 &&
              node_degree(TRUNCATED_CUBOCTAHEDRON) == 1 &&
              node_degree(SNUB_CUBE) == 1 &&
              node_degree(RHOMBICOSIDODECAHEDRON) == 1 &&
              node_degree(TRUNCATED_ICOSIDODECAHEDRON) == 1 &&
              node_degree(SNUB_DODECAHEDRON) == 1);

/** Symmetry families a node belongs to. */
enum class Family : uint8_t { TETRAHEDRAL, OCTAHEDRAL, ICOSAHEDRAL };

/**
 * @brief Symmetry family of a node.
 * @param node Simple-registry node id.
 * @return The rotation-group family the node's chain lives in.
 */
constexpr Family family(int node) {
  switch (node) {
  case TETRAHEDRON:
  case TRUNCATED_TETRAHEDRON:
    return Family::TETRAHEDRAL;
  case CUBE:
  case OCTAHEDRON:
  case CUBOCTAHEDRON:
  case TRUNCATED_CUBE:
  case TRUNCATED_OCTAHEDRON:
  case RHOMBICUBOCTAHEDRON:
  case TRUNCATED_CUBOCTAHEDRON:
  case SNUB_CUBE:
    return Family::OCTAHEDRAL;
  default:
    return Family::ICOSAHEDRAL;
  }
}

// Every bridge row crosses families and every non-bridge row does not.
constexpr bool bridge_flags_consistent() {
  for (int e = 0; e < NUM_EDGES; ++e) {
    bool crosses = family(EDGES[e].from_node) != family(EDGES[e].to_node);
    if (crosses != EDGES[e].bridge)
      return false;
  }
  return true;
}
static_assert(bridge_flags_consistent());

// ---------------------------------------------------------------------------
// Walk policy
// ---------------------------------------------------------------------------

/** Relative pick weight of a non-bridge edge. */
inline constexpr uint32_t WALK_BASE_WEIGHT = 3;
/** Bridge weight before the family is ripe for a change. */
inline constexpr uint32_t WALK_BRIDGE_EARLY_WEIGHT = 1;
/** Bridge weight once BRIDGE_RIPE_LEGS legs were spent in the family. */
inline constexpr uint32_t WALK_BRIDGE_RIPE_WEIGHT = 12;
/** Legs in one family after which bridges become the preferred move. */
inline constexpr int BRIDGE_RIPE_LEGS = 4;

/** Recency scale: non-bridge candidate weights are base * SCALE^3 divided by
 * (1 + target visit count)^3, so seldom-visited targets strongly outweigh the
 * hubs a degree-proportional walk over-samples (10k-leg simulation: max/min
 * node share drops from ~11-15x to ~4x, near the ~3x structural floor set by
 * the cut-vertex hubs' forced pendant transits). Bridges stay recency-exempt
 * at base * SCALE^2 — between fresh and stale sibling weights — which keeps
 * the family-change cadence at the unweighted walk's ~6-7 legs. */
inline constexpr uint32_t WALK_RECENCY_SCALE = 12;
/** Visit-count ceiling: reaching it halves every node's count, so the counts
 * track a sliding window instead of converging and flattening the weighting
 * out on a long walk. */
inline constexpr uint8_t WALK_VISIT_CAP = 8;

/**
 * @brief Pick weight of one candidate edge.
 * @param edge Index into EDGES.
 * @param legs_in_family Completed legs since the last family change.
 * @return Relative weight; bridges outweigh siblings once the family is ripe,
 * so the walk changes family every ~4-6 legs.
 */
constexpr uint32_t edge_weight(int edge, int legs_in_family) {
  if (!EDGES[edge].bridge)
    return WALK_BASE_WEIGHT;
  return legs_in_family >= BRIDGE_RIPE_LEGS ? WALK_BRIDGE_RIPE_WEIGHT
                                            : WALK_BRIDGE_EARLY_WEIGHT;
}

/**
 * @brief Records a node visit for the walk's recency weighting.
 * @param visits Per-node visit counts (NUM_NODES entries), caller-owned.
 * @param node Node just arrived at.
 */
constexpr void record_visit(uint8_t *visits, int node) {
  ++visits[node];
  if (visits[node] >= WALK_VISIT_CAP)
    for (int i = 0; i < NUM_NODES; ++i)
      visits[i] /= 2;
}

/**
 * @brief Random-walk edge choice: weighted random incident edge biased toward
 * less-visited targets, never the immediate backtrack except at a degree-1
 * node.
 * @param node Current node id.
 * @param prev_edge Edge the walk arrived on, or -1 for the first leg.
 * @param legs_in_family Completed legs since the last family change.
 * @param visits Per-node visit counts maintained via record_visit().
 * @param rnd Uniform random 32-bit value (e.g. hs::random()()).
 * @return Index into EDGES of the next leg.
 */
constexpr int pick_next_edge(int node, int prev_edge, int legs_in_family,
                             const uint8_t *visits, uint32_t rnd) {
  uint8_t cand[MAX_DEGREE];
  int n = edges_from(node, cand);
  if (n == 1)
    return cand[0]; // degree-1: out-and-back is the only legal move

  uint32_t total = 0;
  uint32_t weights[MAX_DEGREE] = {};
  constexpr uint32_t S = WALK_RECENCY_SCALE;
  for (int i = 0; i < n; ++i) {
    if (cand[i] == prev_edge)
      continue;
    uint32_t w = 0;
    if (EDGES[cand[i]].bridge) {
      w = edge_weight(cand[i], legs_in_family) * S * S;
    } else {
      const uint32_t rec = 1u + visits[edge_other_end(cand[i], node)];
      w = edge_weight(cand[i], legs_in_family) * S * S * S / (rec * rec * rec);
      if (w == 0)
        w = 1;
    }
    weights[i] = w;
    total += w;
  }

  uint32_t pick = rnd % total;
  for (int i = 0; i < n; ++i) {
    if (weights[i] == 0)
      continue;
    if (pick < weights[i])
      return cand[i];
    pick -= weights[i];
  }
  return cand[n - 1];
}

/** Deterministic profile tour: a fixed edge cycle from the tetrahedron that
 * visits all 18 nodes, traverses every settle edge and every family bridge,
 * and returns to the tetrahedron with the registry seed state, so the cycle
 * wraps seamlessly. Replaces per-node mod arithmetic whose phases never lined
 * up with dodecahedron's rhombicosidodecahedron edge (node 15 stayed
 * uncovered indefinitely). */
inline constexpr uint8_t ORDERED_TOUR[] = {
    18, 20, 6,  7,  5,  5,  1,  3,  3,  4,  4, 0,  2,  8,  19,
    21, 12, 13, 17, 17, 10, 15, 15, 16, 16, 9, 11, 14, 22, 19};
inline constexpr int ORDERED_TOUR_LEN =
    static_cast<int>(std::size(ORDERED_TOUR));

/**
 * @brief Whether ORDERED_TOUR is a closed walk from TETRAHEDRON covering
 * every node.
 */
constexpr bool ordered_tour_valid() {
  bool seen[NUM_NODES] = {};
  int node = TETRAHEDRON;
  seen[node] = true;
  for (int i = 0; i < ORDERED_TOUR_LEN; ++i) {
    const int e = ORDERED_TOUR[i];
    if (e >= NUM_EDGES || !edge_touches(e, node))
      return false;
    node = edge_other_end(e, node);
    seen[node] = true;
  }
  if (node != TETRAHEDRON)
    return false;
  for (int i = 0; i < NUM_NODES; ++i)
    if (!seen[i])
      return false;
  return true;
}
static_assert(ordered_tour_valid());

/**
 * @brief Whether ORDERED_TOUR traverses every settle edge and every family
 * bridge — all three crossings: tetra <-> octa, tetra <-> icosa, and the
 * icosa <-> octa jitterbug (the section-6 heavy legs the profile regime must
 * exercise).
 */
constexpr bool ordered_tour_covers_heavy_legs() {
  bool has[NUM_EDGES] = {};
  for (int i = 0; i < ORDERED_TOUR_LEN; ++i)
    has[ORDERED_TOUR[i]] = true;
  for (int e = 0; e < NUM_EDGES; ++e)
    if ((EDGES[e].settle || EDGES[e].bridge) && !has[e])
      return false;
  return true;
}
static_assert(ordered_tour_covers_heavy_legs());

/**
 * @brief Deterministic edge choice for HS_PROFILE_ORDERED_CYCLE builds: one
 * pass of ORDERED_TOUR covers all 18 nodes, then the cycle repeats.
 * @param node Current node id (positional: the tour expects the walk at its
 * leg_index-th node; the caller's edge_touches check traps a desync).
 * @param prev_edge Edge the walk arrived on; unused.
 * @param leg_index Monotonic leg counter.
 * @return Index into EDGES of the tour's next leg.
 */
constexpr int pick_next_edge_ordered(int node, int prev_edge,
                                     uint32_t leg_index) {
  (void)node;
  (void)prev_edge;
  return ORDERED_TOUR[leg_index % static_cast<uint32_t>(ORDERED_TOUR_LEN)];
}

// ---------------------------------------------------------------------------
// Seed reconciliation
// ---------------------------------------------------------------------------

/**
 * @brief Whether a solid is one of the five Platonic seeds.
 * @param solid Simple-registry index.
 * @return True for tetra/cube/octa/dodeca/icosa.
 */
constexpr bool is_platonic(int solid) { return solid <= ICOSAHEDRON; }

/**
 * @brief Platonic dual partner (tetrahedron is self-dual).
 * @param solid Platonic simple-registry index.
 * @return The dual's registry index, or -1 for a non-Platonic input.
 */
constexpr int dual_platonic(int solid) {
  switch (solid) {
  case TETRAHEDRON:
    return TETRAHEDRON;
  case CUBE:
    return OCTAHEDRON;
  case OCTAHEDRON:
    return CUBE;
  case DODECAHEDRON:
    return ICOSAHEDRON;
  case ICOSAHEDRON:
    return DODECAHEDRON;
  default:
    return -1;
  }
}

/**
 * @brief Action required on the held seed before a leg may start.
 * @details The persistent seed mesh is always one of the five Platonic solids;
 * DERIVE_AMBO legs build their working seed as ambo(held seed) at construction
 * without replacing it, which is what lets the out-and-back through
 * truncatedCuboctahedron / truncatedIcosidodecahedron return to a node whose
 * departures all need a Platonic seed.
 */
enum class SeedFix : uint8_t {
  KEEP,        /**< Held seed already matches the edge's seed_solid. */
  DUAL_SWAP,   /**< seed := dual(seed) at the ambo crossover (t = 0.5). */
  DERIVE_AMBO, /**< Leg seed := ambo(held seed); held seed unchanged. */
  REGEN_TETRA, /**< seed := registry tetrahedron (reverse family bridge). */
  INVALID,     /**< No legal reconciliation; the caller must trap. */
};

/**
 * @brief Resolves the seed reconciliation for departing on an edge.
 * @param edge Index into EDGES of the leg about to start.
 * @param held_identity Registry index of the solid the persistent seed mesh
 * represents (always Platonic).
 * @return The required SeedFix; INVALID when the walk state is corrupt.
 */
constexpr SeedFix seed_fix_at_start(int edge, int held_identity) {
  const EdgeSpec &e = EDGES[edge];
  if (e.seed_solid == held_identity)
    return SeedFix::KEEP;
  if (!is_platonic(e.seed_solid)) {
    if (e.seed_solid == CUBOCTAHEDRON &&
        (held_identity == CUBE || held_identity == OCTAHEDRON))
      return SeedFix::DERIVE_AMBO;
    if (e.seed_solid == ICOSIDODECAHEDRON &&
        (held_identity == DODECAHEDRON || held_identity == ICOSAHEDRON))
      return SeedFix::DERIVE_AMBO;
    return SeedFix::INVALID;
  }
  if (e.seed_solid == dual_platonic(held_identity))
    return SeedFix::DUAL_SWAP;
  if (e.seed_solid == TETRAHEDRON &&
      (held_identity == OCTAHEDRON || held_identity == ICOSAHEDRON))
    return SeedFix::REGEN_TETRA;
  return SeedFix::INVALID;
}

} // namespace ConwayGraph
