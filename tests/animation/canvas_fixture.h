/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_animation.h.

// ============================================================================
// Stand-in Canvas reference
// ----------------------------------------------------------------------------
// One genuine Canvas over a tiny Effect, shared across tests. The animations
// under test take a Canvas& but never dereference it. No test may construct a
// second Canvas on this shared effect: that would queue_frame() and spin the
// next ctor on a display ISR the host never runs.
// ============================================================================

/**
 * @brief Storage for the module-scoped fake-canvas pointer.
 * @return Reference to the pointer run_animation_tests() sets to the shared
 * fixture for the module's duration and clears (fixture destroyed) on exit.
 */
inline Canvas *&fake_canvas_ptr() {
  static Canvas *p = nullptr;
  return p;
}

/**
 * @brief Returns the module-scoped shared Canvas reused across the tests.
 * @return Reference to the fixture owned by run_animation_tests(); valid for the
 * duration of the module but never dereferenced by the animations under test.
 */
inline Canvas &fake_canvas() { return *fake_canvas_ptr(); }

/**
 * @brief Builds a hand-rolled octahedron (6 axis vertices, 8 triangular faces).
 * @param mesh Output mesh to populate with the octahedron geometry.
 * @param arena Arena backing the mesh's vertex/face/face-count buffers.
 */
inline void build_octahedron(PolyMesh &mesh, Arena &arena) {
  static const math::Vector verts[6] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                        {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  static const uint16_t tris[8][3] = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4},
                                      {3, 0, 4}, {2, 0, 5}, {1, 2, 5},
                                      {3, 1, 5}, {0, 3, 5}};
  mesh.vertices.bind(arena, 6);
  mesh.face_counts.bind(arena, 8);
  mesh.faces.bind(arena, 24);
  for (const auto &v : verts)
    mesh.vertices.push_back(v);
  for (const auto &t : tris) {
    mesh.face_counts.push_back(3);
    mesh.faces.push_back(t[0]);
    mesh.faces.push_back(t[1]);
    mesh.faces.push_back(t[2]);
  }
}
