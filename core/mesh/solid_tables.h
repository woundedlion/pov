/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/mesh/solid_generators.h.

// ==========================================================================================
// DATA DEFINITIONS (Hardcoded Platonic Solids)
// ==========================================================================================

/**
 * @brief Tetrahedron geometry data.
 */
struct Tetrahedron {
  static constexpr int NUM_VERTS = 4;
  static constexpr std::array<math::Vector, NUM_VERTS> vertices = {
      math::Vector(0.5773502691896258f, 0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(0.5773502691896258f, -0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, 0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, -0.5773502691896258f,
                   0.5773502691896258f)};
  static constexpr int NUM_FACES = 4;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {3, 3, 3, 3};
  static constexpr std::array<int, 12> faces = {0, 3, 1, 0, 2, 3,
                                                0, 1, 2, 1, 3, 2};
};

/**
 * @brief Cube geometry data.
 */
struct Cube {
  static constexpr int NUM_VERTS = 8;
  static constexpr std::array<math::Vector, NUM_VERTS> vertices = {
      math::Vector(-0.5773502691896258f, -0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(0.5773502691896258f, -0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(0.5773502691896258f, 0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, 0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, -0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(0.5773502691896258f, -0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(0.5773502691896258f, 0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(-0.5773502691896258f, 0.5773502691896258f,
                   0.5773502691896258f)};
  static constexpr int NUM_FACES = 6;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {4, 4, 4,
                                                                 4, 4, 4};
  static constexpr std::array<int, 24> faces = {
      0, 3, 2, 1, 0, 1, 5, 4, 0, 4, 7, 3, 6, 5, 1, 2, 6, 2, 3, 7, 6, 7, 4, 5};
};

/**
 * @brief Octahedron geometry data.
 */
struct Octahedron {
  static constexpr int NUM_VERTS = 6;
  static constexpr std::array<math::Vector, NUM_VERTS> vertices = {
      math::Vector(1.0000000000000000f, 0.0000000000000000f,
                   0.0000000000000000f),
      math::Vector(-1.0000000000000000f, 0.0000000000000000f,
                   0.0000000000000000f),
      math::Vector(0.0000000000000000f, 1.0000000000000000f,
                   0.0000000000000000f),
      math::Vector(0.0000000000000000f, -1.0000000000000000f,
                   0.0000000000000000f),
      math::Vector(0.0000000000000000f, 0.0000000000000000f,
                   1.0000000000000000f),
      math::Vector(0.0000000000000000f, 0.0000000000000000f,
                   -1.0000000000000000f)};
  static constexpr int NUM_FACES = 8;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {3, 3, 3, 3,
                                                                 3, 3, 3, 3};
  static constexpr std::array<int, 24> faces = {
      4, 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 5, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3};
};

/**
 * @brief Icosahedron geometry data.
 */
struct Icosahedron {
  static constexpr int NUM_VERTS = 12;
  static constexpr std::array<math::Vector, NUM_VERTS> vertices = {
      math::Vector(-0.5257311121191336f, 0.0000000000000000f,
                   0.8506508083520400f),
      math::Vector(0.5257311121191336f, 0.0000000000000000f,
                   0.8506508083520400f),
      math::Vector(-0.5257311121191336f, 0.0000000000000000f,
                   -0.8506508083520400f),
      math::Vector(0.5257311121191336f, 0.0000000000000000f,
                   -0.8506508083520400f),
      math::Vector(0.0000000000000000f, 0.8506508083520400f,
                   0.5257311121191336f),
      math::Vector(0.0000000000000000f, 0.8506508083520400f,
                   -0.5257311121191336f),
      math::Vector(0.0000000000000000f, -0.8506508083520400f,
                   0.5257311121191336f),
      math::Vector(0.0000000000000000f, -0.8506508083520400f,
                   -0.5257311121191336f),
      math::Vector(0.8506508083520400f, 0.5257311121191336f,
                   0.0000000000000000f),
      math::Vector(-0.8506508083520400f, 0.5257311121191336f,
                   0.0000000000000000f),
      math::Vector(0.8506508083520400f, -0.5257311121191336f,
                   0.0000000000000000f),
      math::Vector(-0.8506508083520400f, -0.5257311121191336f,
                   0.0000000000000000f)};
  static constexpr int NUM_FACES = 20;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  static constexpr std::array<int, 60> faces = {
      0, 1, 4, 0, 4, 9, 9,  4, 5, 4,  8, 5, 4,  1,  8,  8, 1, 10, 8,  10,
      3, 5, 8, 3, 5, 3, 2,  2, 3, 7,  7, 3, 10, 7,  10, 6, 7, 6,  11, 11,
      6, 0, 0, 6, 1, 6, 10, 1, 9, 11, 0, 9, 2,  11, 9,  5, 2, 7,  11, 2};
};

/**
 * @brief Dodecahedron geometry data.
 */
struct Dodecahedron {
  static constexpr int NUM_VERTS = 20;
  static constexpr std::array<math::Vector, NUM_VERTS> vertices = {
      math::Vector(0.5773502691896258f, 0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(0.5773502691896258f, 0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(0.5773502691896258f, -0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(0.5773502691896258f, -0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, 0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(-0.5773502691896258f, 0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(-0.5773502691896258f, -0.5773502691896258f,
                   0.5773502691896258f),
      math::Vector(-0.5773502691896258f, -0.5773502691896258f,
                   -0.5773502691896258f),
      math::Vector(0.3568220897730897f, 0.9341723589627157f,
                   0.0000000000000000f),
      math::Vector(-0.3568220897730897f, 0.9341723589627157f,
                   0.0000000000000000f),
      math::Vector(0.3568220897730897f, -0.9341723589627157f,
                   0.0000000000000000f),
      math::Vector(-0.3568220897730897f, -0.9341723589627157f,
                   0.0000000000000000f),
      math::Vector(0.9341723589627157f, 0.0000000000000000f,
                   0.3568220897730897f),
      math::Vector(0.9341723589627157f, 0.0000000000000000f,
                   -0.3568220897730897f),
      math::Vector(-0.9341723589627157f, 0.0000000000000000f,
                   0.3568220897730897f),
      math::Vector(-0.9341723589627157f, 0.0000000000000000f,
                   -0.3568220897730897f),
      math::Vector(0.0000000000000000f, 0.3568220897730897f,
                   0.9341723589627157f),
      math::Vector(0.0000000000000000f, -0.3568220897730897f,
                   0.9341723589627157f),
      math::Vector(0.0000000000000000f, 0.3568220897730897f,
                   -0.9341723589627157f),
      math::Vector(0.0000000000000000f, -0.3568220897730897f,
                   -0.9341723589627157f)};
  static constexpr int NUM_FACES = 12;
  static constexpr std::array<uint8_t, NUM_FACES> face_counts = {
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  static constexpr std::array<int, 60> faces = {
      0, 8,  9,  4,  16, 0,  12, 13, 1,  8,  0,  16, 17, 2,  12,
      8, 1,  18, 5,  9,  12, 2,  10, 3,  13, 16, 4,  14, 6,  17,
      9, 5,  15, 14, 4,  6,  11, 10, 2,  17, 3,  19, 18, 1,  13,
      7, 15, 5,  18, 19, 7,  11, 6,  14, 15, 7,  19, 3,  10, 11};
};

/**
 * @brief Compile-time consistency check for a hardcoded solid's tables.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces arrays.
 * @return True when index, edge-incidence, Euler and unit-length checks pass.
 * @details Checks that face_counts spans the flat face list exactly, every face
 * index addresses a listed vertex, no directed edge repeats and every directed
 * edge has its reverse (so each undirected edge joins exactly two faces in
 * opposite orientation, making E = sum/2), Euler's formula holds, and every
 * squared vertex length passes the 1 ± 1e-4 comparisons. Finiteness, vertex fans,
 * connectivity, face planarity, convexity and non-self-intersection are not checked.
 * Compares squared lengths so the check is
 * constant-evaluable.
 */
template <typename StaticMeshT> constexpr bool solid_tables_consistent() {
  size_t total_indices = 0;
  for (uint8_t count : StaticMeshT::face_counts)
    total_indices += count;
  if (total_indices != StaticMeshT::faces.size() || total_indices % 2 != 0)
    return false;

  for (int index : StaticMeshT::faces)
    if (index < 0 || static_cast<size_t>(index) >= StaticMeshT::vertices.size())
      return false;

  // Directed edge incidence; the V*V scratch keeps the pass O(indices + V^2).
  constexpr size_t V = StaticMeshT::vertices.size();
  std::array<bool, V * V> used = {};
  size_t base = 0;
  for (uint8_t count : StaticMeshT::face_counts) {
    for (uint8_t i = 0; i < count; ++i) {
      const size_t a = static_cast<size_t>(StaticMeshT::faces[base + i]);
      const size_t b = static_cast<size_t>(
          StaticMeshT::faces[base + (i + 1 == count ? 0 : i + 1)]);
      if (a == b || used[a * V + b])
        return false;
      used[a * V + b] = true;
    }
    base += count;
  }
  for (size_t a = 0; a < V; ++a)
    for (size_t b = 0; b < a; ++b)
      if (used[a * V + b] != used[b * V + a])
        return false;

  // V - E + F == 2, in a form free of unsigned wrap.
  if (StaticMeshT::vertices.size() + StaticMeshT::face_counts.size() !=
      total_indices / 2 + 2)
    return false;

  constexpr float LEN_SQ_TOL = 1e-4f;
  for (const math::Vector &v : StaticMeshT::vertices) {
    const float len_sq = v.x * v.x + v.y * v.y + v.z * v.z;
    if (len_sq < 1.0f - LEN_SQ_TOL || len_sq > 1.0f + LEN_SQ_TOL)
      return false;
  }
  return true;
}

static_assert(solid_tables_consistent<Tetrahedron>(),
              "Tetrahedron tables are inconsistent");
static_assert(solid_tables_consistent<Cube>(), "Cube tables are inconsistent");
static_assert(solid_tables_consistent<Octahedron>(),
              "Octahedron tables are inconsistent");
static_assert(solid_tables_consistent<Icosahedron>(),
              "Icosahedron tables are inconsistent");
static_assert(solid_tables_consistent<Dodecahedron>(),
              "Dodecahedron tables are inconsistent");

/**
 * @brief Materializes a compile-time static mesh into a runtime PolyMesh.
 * @tparam StaticMeshT Type exposing constexpr vertices/face_counts/faces
 * arrays.
 * @param target Arena that backs the returned mesh's storage.
 * @return A PolyMesh holding copies of the static mesh's data in target.
 * @details Face indices narrow to uint16_t, trapping past MeshLimits::MAX_VERTEX_INDEX.
 */
template <typename StaticMeshT> PolyMesh to_polymesh(Arena &target) {
  PolyMesh mesh;
  mesh.vertices.bind(target, StaticMeshT::vertices.size());
  mesh.vertices.append_bulk(StaticMeshT::vertices.data(),
                            StaticMeshT::vertices.size());
  mesh.face_counts.bind(target, StaticMeshT::face_counts.size());
  mesh.face_counts.append_bulk(StaticMeshT::face_counts.data(),
                               StaticMeshT::face_counts.size());
  mesh.faces.bind(target, StaticMeshT::faces.size());
  for (const auto &f : StaticMeshT::faces)
    mesh.faces.push_back(MeshOps::narrow_index(static_cast<size_t>(f)));
  return mesh;
}
