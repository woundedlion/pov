/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_conway_morph.h.

// ---------------------------------------------------------------------------
// Morph-frame scratch high-water gate at HankinSolids' shipping split
// A morph frame runs one op plus MeshOps::compile in the scratch pair under
// LIFO scopes;
// host high-water marks are a conservative upper bound on the device figure.
// ---------------------------------------------------------------------------

/** The arena split is canvas-independent; this instantiation names it. */
using HankinFx = HankinSolids<96, 20>;

constexpr size_t MORPH_SCRATCH_A_BUDGET =
    HankinFx::SCRATCH_A_BYTES; /**< HankinSolids scratch_a split. */
constexpr size_t MORPH_SCRATCH_B_BUDGET =
    HankinFx::SCRATCH_B_BYTES; /**< HankinSolids scratch_b split. */

/**
 * @brief Verifies every edge's op-plus-compile scratch peak fits
 *        HankinSolids' 24 KB / 32 KB scratch split.
 * @details The seed is persistent and topology is t-constant; one mid-sweep
 * sample measures the op-plus-compile peak. OpLeg's constructor checks the
 * blended-LUT term; effect smoke tests exercise the additional draw stack.
 * Reports the worst arena pair across the edge table.
 */
inline void test_edge_morph_frames_fit_scratch_budget() {
  constexpr size_t HALF = sizeof(morph_aux_buf) / 2;
  size_t worst_a = 0, worst_b = 0;
  int worst_a_edge = 0, worst_b_edge = 0;

  for (int ei = 0; ei < ConwayGraph::NUM_EDGES; ++ei) {
    const ConwayGraph::EdgeSpec &e = ConwayGraph::EDGES[ei];

    Arena persist(morph_persist_buf, sizeof(morph_persist_buf));
    PolyMesh seed;
    {
      Arena ga(morph_aux_buf, HALF);
      Arena gb(morph_aux_buf + HALF, HALF);
      seed = Solids::finalize_solid(
          Solids::simple_registry[e.seed_solid].generate(ga, gb), persist);
    }

    float t_lo, t_hi;
    edge_sweep_interval(e, t_lo, t_hi);
    constexpr float U = 0.5f;
    const float t = t_lo + (t_hi - t_lo) * U;
    const float twist = e.twist_from + (e.twist_to - e.twist_from) * U;

    Arena a(morph_target_buf, MORPH_SCRATCH_A_BUDGET);
    Arena b(morph_temp_buf, MORPH_SCRATCH_B_BUDGET);
    {
      ScratchScope frame_a(a);
      ScratchScope frame_b(b);
      PolyMesh swept = run_edge_op(e, seed, a, b, t, twist);
      MeshState frame;
      MeshOps::compile(swept, frame, a, b);
    }

    const size_t a_peak = a.get_high_water_mark();
    const size_t b_peak = b.get_high_water_mark();
    if (a_peak > worst_a) {
      worst_a = a_peak;
      worst_a_edge = ei;
    }
    if (b_peak > worst_b) {
      worst_b = b_peak;
      worst_b_edge = ei;
    }
    HS_EXPECT_LE(a_peak, MORPH_SCRATCH_A_BUDGET);
    HS_EXPECT_LE(b_peak, MORPH_SCRATCH_B_BUDGET);
  }

  const ConwayGraph::EdgeSpec &wa = ConwayGraph::EDGES[worst_a_edge];
  const ConwayGraph::EdgeSpec &wb = ConwayGraph::EDGES[worst_b_edge];
  std::printf(
      "  [morph scratch] worst a=%zu B (%s -> %s) / budget=%zu B, "
      "worst b=%zu B (%s -> %s) / budget=%zu B\n",
      worst_a, Solids::simple_registry[wa.from_node].name,
      Solids::simple_registry[wa.to_node].name, (size_t)MORPH_SCRATCH_A_BUDGET,
      worst_b, Solids::simple_registry[wb.from_node].name,
      Solids::simple_registry[wb.to_node].name, (size_t)MORPH_SCRATCH_B_BUDGET);
}
