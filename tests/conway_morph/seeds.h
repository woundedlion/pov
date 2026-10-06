/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Seeds: the solids the edge table sweeps from.
// ---------------------------------------------------------------------------

/** @brief Sweep seeds of the OpLeg edge table. */
enum class MorphSeed {
  TETRAHEDRON,
  CUBE,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  CUBOCTAHEDRON,
  ICOSIDODECAHEDRON,
};

inline constexpr MorphSeed MORPH_SEEDS[] = {
    MorphSeed::TETRAHEDRON,       MorphSeed::CUBE,
    MorphSeed::OCTAHEDRON,        MorphSeed::DODECAHEDRON,
    MorphSeed::ICOSAHEDRON,       MorphSeed::CUBOCTAHEDRON,
    MorphSeed::ICOSIDODECAHEDRON,
};

/**
 * @brief Seed name for failure diagnostics.
 * @param s Seed identifier.
 * @return Static name string.
 */
inline const char *seed_name(MorphSeed s) {
  switch (s) {
  case MorphSeed::TETRAHEDRON:
    return "tetrahedron";
  case MorphSeed::CUBE:
    return "cube";
  case MorphSeed::OCTAHEDRON:
    return "octahedron";
  case MorphSeed::DODECAHEDRON:
    return "dodecahedron";
  case MorphSeed::ICOSAHEDRON:
    return "icosahedron";
  case MorphSeed::CUBOCTAHEDRON:
    return "cuboctahedron";
  case MorphSeed::ICOSIDODECAHEDRON:
    return "icosidodecahedron";
  }
  return "?";
}

/**
 * @brief Builds a sweep seed mesh.
 * @param s Seed identifier.
 * @param target Arena receiving the seed mesh.
 * @param temp Scratch arena (holds the platonic base for the ambo seeds).
 * @return The seed PolyMesh in `target`.
 */
inline PolyMesh build_morph_seed(MorphSeed s, Arena &target, Arena &temp) {
  PolyMesh m;
  switch (s) {
  case MorphSeed::TETRAHEDRON:
    build_solid<Solids::Tetrahedron>(m, target);
    return m;
  case MorphSeed::CUBE:
    build_solid<Solids::Cube>(m, target);
    return m;
  case MorphSeed::OCTAHEDRON:
    build_solid<Solids::Octahedron>(m, target);
    return m;
  case MorphSeed::DODECAHEDRON:
    build_solid<Solids::Dodecahedron>(m, target);
    return m;
  case MorphSeed::ICOSAHEDRON:
    build_solid<Solids::Icosahedron>(m, target);
    return m;
  case MorphSeed::CUBOCTAHEDRON:
    build_solid<Solids::Cube>(m, temp);
    return MeshOps::ambo(m, target, temp);
  case MorphSeed::ICOSIDODECAHEDRON:
    build_solid<Solids::Dodecahedron>(m, temp);
    return MeshOps::ambo(m, target, temp);
  }
  return m;
}
