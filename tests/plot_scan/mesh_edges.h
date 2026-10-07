/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file mesh_edges.h
 * @brief Plot mesh edge extraction and dissolve partition tests.
 */

/**
 * @brief The complementary masks of a Segue::Dissolve partition a wireframe's
 *        edges exactly: every edge is drawn by one sprite and skipped by the
 *        other, at every phase.
 */
inline void test_mesh_dissolve_masks_partition_edges() {
  constexpr int W = 96, H = 48;
  configure_arenas_default();

  alignas(32) static uint8_t seed_a[24 * 1024];
  alignas(32) static uint8_t seed_b[24 * 1024];
  alignas(32) static uint8_t geom[16 * 1024];
  Arena sa(seed_a, sizeof(seed_a));
  Arena sb(seed_b, sizeof(seed_b));
  Arena ga(geom, sizeof(geom));

  MeshState mesh;
  build_icosahedron_meshstate(sa, sb, ga, mesh);

  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(ga, mesh.faces.size());
  Plot::Mesh::extract_edges(mesh, edges);
  const size_t num_edges = edges.size();

  Segue::Dissolve dissolve;
  hs::random().seed(0xD155);
  dissolve.retarget(math::Y_AXIS);

  auto drawn_set = [&](const DissolveMask &mask) {
    std::vector<bool> seen(num_edges, false);
    auto shade = [&](const math::Vector &, Fragment &f) {
      const int ei = static_cast<int>(f.v2);
      HS_EXPECT_TRUE(ei >= 0 && static_cast<size_t>(ei) < num_edges);
      if (ei < 0 || static_cast<size_t>(ei) >= num_edges)
        return;
      seen[static_cast<size_t>(ei)] = true;
      f.color = Color4(Pixel(65535, 65535, 65535), 0.9f);
    };
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters{
        Filter::Screen::AntiAlias<W, H>()};
    Canvas c(fx);
    Plot::Mesh::draw<W, H>(filters, c, mesh, edges, shade, {}, &mask);
    return seen;
  };

  const float phases[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (uint32_t frame = 0; frame < 3; ++frame) {
    for (float p : phases) {
      const auto masks = dissolve.mask_pair(p, frame);
      auto in_set = drawn_set(masks.incoming);
      auto out_set = drawn_set(masks.outgoing);
      size_t in_count = 0;
      for (size_t e = 0; e < num_edges; ++e) {
        HS_EXPECT_TRUE(in_set[e] != out_set[e]);
        in_count += in_set[e] ? 1 : 0;
      }
      if (p == 0.0f)
        HS_EXPECT_EQ(in_count, size_t{0});
      if (p == 1.0f)
        HS_EXPECT_EQ(in_count, num_edges);
      if (p == 0.5f) {
        HS_EXPECT_GT(in_count, size_t{0});
        HS_EXPECT_GT(num_edges - in_count, size_t{0});
      }
    }
  }
}

/** @brief Four-regular extraction covers every edge once and medial indices name original edges. */
inline void test_four_regular_and_medial_edge_extraction() {
  configure_arenas_default();
  alignas(std::max_align_t) static uint8_t storage[16 * 1024];
  Arena arena(storage, sizeof(storage));
  MeshState mesh;
  build_meshstate_solid<Solids::Octahedron>(mesh, arena);
  ArenaVector<Plot::Mesh::Edge> unique(arena, mesh.faces.size());
  ArenaVector<Plot::Mesh::Edge> woven(arena, mesh.faces.size());
  ArenaVector<Plot::Mesh::Edge> medial(arena, mesh.faces.size());
  Plot::Mesh::extract_edges(mesh, unique);
  Plot::Mesh::extract_four_regular_edges(mesh, woven, scratch_arena_a);
  HS_EXPECT_SIZE_OR_RETURN(woven, unique.size());
  for (const auto &edge : unique) {
    int matches = 0;
    for (const auto &candidate : woven)
      matches += ((edge.u == candidate.u && edge.v == candidate.v) ||
                  (edge.u == candidate.v && edge.v == candidate.u));
    HS_EXPECT_EQ(matches, 1);
    const auto index = Plot::Mesh::find_edge_index(woven, edge.v, edge.u);
    HS_EXPECT_LT(index, woven.size());
    if (index >= woven.size())
      continue;
    const auto &found = woven[index];
    HS_EXPECT_TRUE((found.u == edge.u && found.v == edge.v) ||
                   (found.u == edge.v && found.v == edge.u));
  }
  Plot::Mesh::extract_medial_edges(mesh, unique, medial);
  HS_EXPECT_SIZE_OR_RETURN(medial, mesh.faces.size());
  std::vector<int> incidence(unique.size(), 0);
  for (const auto &edge : medial) {
    HS_EXPECT_LT(edge.u, unique.size());
    HS_EXPECT_LT(edge.v, unique.size());
    HS_EXPECT_NE(edge.u, edge.v);
    if (edge.u >= unique.size() || edge.v >= unique.size())
      continue;
    const auto &a = unique[edge.u];
    const auto &b = unique[edge.v];
    HS_EXPECT_EQ((a.u == b.u) + (a.u == b.v) + (a.v == b.u) + (a.v == b.v), 1);
    ++incidence[edge.u];
    ++incidence[edge.v];
  }
  for (const int count : incidence)
    HS_EXPECT_EQ(count, 4);
}
