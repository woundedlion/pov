/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Smooth kis/needle reconcile: the identity meshes (dt / dtd) and the authored
// kis/needle meshes are topology-exact (identical V/E/F/I); the residual gap is
// closed by a per-vertex great-circle slerp along the nearest-vertex bijection.
// ---------------------------------------------------------------------------

/** One smooth-kis/needle macro site: how X is reached, and which identity vs
 * authored meshes the reconcile leg spans. */
struct ReconcileSite {
  const char *name;
  const Solids::Recipe *recipe;
  bool dtd; /**< true = standalone kis (dtd); false = trailing dual,kis (dt). */
};

inline const ReconcileSite RECONCILE_SITES[] = {
    {"needle", &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE, false},
    {"gyro_kis", &Solids::TRUNCATED_OCTAHEDRON_GYRO_KIS_HK17_RECIPE, false},
    {"truncatedCube_hk90_meta", &Solids::TRUNCATED_CUBE_HK90_META_RECIPE,
     false},
    {"icosahedron_kis_gyro", &Solids::ICOSAHEDRON_KIS_GYRO_RECIPE, true},
};

/**
 * @brief Verifies the Conway-identity reconcile is well-posed on every smooth
 *        kis/needle seed: the identity mesh (dt/dtd) and the authored kis/needle
 *        mesh share V/E/F, the nearest-vertex map is a bijection, and the
 *        residual chord the reconcile slerp closes is bounded.
 * @details The brute-force argmax-dot map here is an independent oracle;
 * MeshOps::reconcile_vertices (the shipping z-banded search) is run on the same
 * endpoints and must reproduce it vertex for vertex.
 */
inline void test_reconcile_bijection_wellposed() {
  // Bound the gap closed by the reconcile slerp.
  constexpr float MAX_RESIDUAL_CHORD = 0.12f;
  for (const ReconcileSite &site : RECONCILE_SITES) {
    const int failed_before = hs_test::stats().failed;
    Arena a(morph_target_buf, sizeof(morph_target_buf));
    Arena b(morph_temp_buf, sizeof(morph_temp_buf));
    Arena aux(morph_aux_buf, sizeof(morph_aux_buf));
    Arena keep(morph_persist_buf, sizeof(morph_persist_buf));

    Solids::OpStep lowered[8];
    const size_t count = Solids::expand_to_primitives(*site.recipe, lowered, 8);
    // Locate the kis and its preceding dual (dt pair) or its standalone form.
    size_t kis = count;
    for (size_t k = 0; k < count; ++k)
      if (lowered[k].op == Solids::Op::KIS) {
        kis = k;
        break;
      }
    HS_EXPECT_LT(kis, count);
    HS_EXPECT_TRUE(site.dtd || kis > 0);
    if (kis >= count || (!site.dtd && kis == 0))
      continue;
    const size_t x_prefix = site.dtd ? kis : kis - 1; // steps to reach X

    PolyMesh X =
        Solids::build_steps(site.recipe->seed, lowered, x_prefix, a, b);
    PolyMesh Xc; // stable copy: identity/authored reuse a,b as ping-pong
    MeshOps::clone(X, Xc, aux);

    // identity: dt = dual(truncate(X)); dtd = dual(truncate(dual(X))). Clone the
    // identity out of a/b before building authored, which reuses those arenas.
    PolyMesh identity;
    PolyMesh authored;
    if (site.dtd) {
      PolyMesh dX = MeshOps::dual(Xc, a, b);
      PolyMesh dXc;
      MeshOps::clone(dX, dXc, aux);
      PolyMesh t = MeshOps::truncate(dXc, a, b, RECONCILE_TRUNCATE_T);
      MeshOps::clone(MeshOps::dual(t, b, a), identity, keep);
      authored = MeshOps::kis(Xc, a, b); // kis(X)
    } else {
      PolyMesh t = MeshOps::truncate(Xc, a, b, RECONCILE_TRUNCATE_T);
      MeshOps::clone(MeshOps::dual(t, b, a), identity, keep);
      authored = MeshOps::needle(Xc, a, b); // kis(dual(X))
    }

    const size_t V = identity.vertices.size();
    HS_EXPECT_EQ(V, authored.vertices.size());
    HS_EXPECT_EQ(identity.face_counts.size(), authored.face_counts.size());
    HS_EXPECT_EQ(identity.faces.size(), authored.faces.size());
    if (V != authored.vertices.size() ||
        identity.face_counts.size() != authored.face_counts.size() ||
        identity.faces.size() != authored.faces.size())
      continue;

    // Per-vertex nearest (argmax-dot) map identity -> authored;
    // must be injective.
    std::vector<bool> used(V, false);
    std::vector<int> match(V, -1);
    float worst_chord = 0.0f;
    bool injective = true;
    for (size_t i = 0; i < V; ++i) {
      int best = -1;
      float best_dot = -2.0f;
      for (size_t j = 0; j < V; ++j) {
        const float d = math::dot(identity.vertices[i], authored.vertices[j]);
        if (d > best_dot) {
          best_dot = d;
          best = static_cast<int>(j);
        }
      }
      if (best < 0 || used[best])
        injective = false;
      else
        used[best] = true;
      match[i] = best;
      if (best >= 0)
        worst_chord = fold_worst(
            worst_chord, math::distance_between(identity.vertices[i],
                                                authored.vertices[best]));
    }
    HS_EXPECT_TRUE(injective);
    HS_EXPECT_LT(worst_chord, MAX_RESIDUAL_CHORD);

    if (injective) {
      PolyMesh reconciled;
      MeshOps::reconcile_vertices(identity, authored, reconciled, b, aux);
      HS_EXPECT_EQ(reconciled.vertices.size(), V);
      HS_EXPECT_EQ(reconciled.face_counts.size(), identity.face_counts.size());
      HS_EXPECT_EQ(reconciled.faces.size(), identity.faces.size());
      if (reconciled.vertices.size() != V ||
          reconciled.face_counts.size() != identity.face_counts.size() ||
          reconciled.faces.size() != identity.faces.size())
        continue;
      bool connectivity_kept = true;
      for (size_t f = 0; f < identity.face_counts.size(); ++f)
        connectivity_kept &=
            reconciled.face_counts[f] == identity.face_counts[f];
      for (size_t i = 0; i < identity.faces.size(); ++i)
        connectivity_kept &= reconciled.faces[i] == identity.faces[i];
      HS_EXPECT_TRUE(connectivity_kept);
      bool matches_oracle = true;
      for (size_t i = 0; i < V; ++i)
        matches_oracle &= reconciled.vertices[i] ==
                          authored.vertices[static_cast<size_t>(match[i])];
      HS_EXPECT_TRUE(matches_oracle);
    }

    if (hs_test::stats().failed != failed_before)
      std::printf("    [reconcile] %s FAILED (V=%zu inj=%d worst_chord=%.4f)\n",
                  site.name, V, injective ? 1 : 0, (double)worst_chord);
    else
      std::printf("  [reconcile] %s: V=%zu F=%zu bijective worst_chord=%.4f\n",
                  site.name, V, identity.face_counts.size(),
                  (double)worst_chord);
  }
}
