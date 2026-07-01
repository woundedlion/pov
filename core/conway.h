/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "mesh.h"

/**
 * @brief Operations on meshes (Conway operators and utilities).
 */
namespace MeshOps {

// `narrow_index` (the trapping size_t -> uint16_t topology-index cast) lives in
// mesh.h so the Conway and Hankin operators share a single guarded narrowing.

/**
 * @brief Compute the centroid of a face by walking its half-edge loop.
 * @tparam MeshT Mesh type exposing a `vertices` indexable by topology index.
 * @param he_mesh Half-edge connectivity describing the face loops.
 * @param mesh Source mesh supplying the vertex positions.
 * @param face_index Index of the face whose centroid is computed.
 * @param out_count Out-param set to the number of vertices walked (sides).
 * @return Average of the face's vertex positions, or origin for an empty face.
 */
template <typename MeshT>
inline Vector face_centroid(const HalfEdgeMesh &he_mesh,
                            const MeshT &mesh, size_t face_index,
                            int &out_count) {
  const HEFace &face = he_mesh.faces[face_index];
  Vector c(0, 0, 0);
  out_count = 0;
  uint16_t he_idx = face.half_edge;
  uint16_t start = he_idx;
  if (he_idx != HE_NONE) {
    // Anti-hang guard: a corrupt .next chain would otherwise spin forever.
    const int max_sides = static_cast<int>(he_mesh.half_edges.size());
    do {
      HS_CHECK(out_count < max_sides, "face_centroid: corrupt face loop");
      c = c + mesh.vertices[he_mesh.half_edges[he_idx].vertex];
      out_count++;
      he_idx = he_mesh.half_edges[he_idx].next;
    } while (he_idx != HE_NONE && he_idx != start);
  }
  if (out_count > 0)
    c = c / static_cast<float>(out_count);
  return c;
}

/**
 * @brief Newell's method face normal, robust to non-planar faces and collinear
 *   vertex triplets.
 * @tparam MeshT Mesh type exposing a `vertices` indexable by topology index.
 * @param he_mesh Half-edge connectivity describing the face loops.
 * @param mesh Source mesh supplying the vertex positions.
 * @param face_index Index of the face whose normal is computed.
 * @return Unnormalized Newell normal (the caller normalizes if needed); the
 *   origin vector for an empty face.
 * @details Walks the half-edge loop, summing the Newell contribution of each
 *   edge (consecutive vertex pair).
 */
template <typename MeshT>
inline Vector face_normal(const HalfEdgeMesh &he_mesh, const MeshT &mesh,
                          size_t face_index) {
  const HEFace &face = he_mesh.faces[face_index];
  Vector n(0, 0, 0);
  uint16_t he_idx = face.half_edge;
  if (he_idx == HE_NONE) return n;
  uint16_t start = he_idx;
  // Anti-hang guard plus an HE_NONE trap: this loop dereferences
  // half_edges[he.next] in-body, so a corrupt .next chain could spin or index HE_NONE.
  const int max_sides = static_cast<int>(he_mesh.half_edges.size());
  int sides = 0;
  do {
    HS_CHECK(sides++ < max_sides, "face_normal: corrupt face loop");
    const HalfEdge &he = he_mesh.half_edges[he_idx];
    HS_CHECK(he.next != HE_NONE, "face_normal: HE_NONE in face loop");
    const Vector &curr = mesh.vertices[he.vertex];
    const Vector &next = mesh.vertices[he_mesh.half_edges[he.next].vertex];
    n.x += (curr.y - next.y) * (curr.z + next.z);
    n.y += (curr.z - next.z) * (curr.x + next.x);
    n.z += (curr.x - next.x) * (curr.y + next.y);
    he_idx = he.next;
  } while (he_idx != HE_NONE && he_idx != start);
  return n;
}

/**
 * @brief Walk all half-edges orbiting a vertex, invoking visitor(curr_idx) for
 *   each.
 * @tparam OrbitMode Traversal direction: 'P' = prev->pair (dual, ambo,
 *   truncate); 'N' = pair->next (expand, snub).
 * @tparam VisitorFn Callable accepting the current half-edge index (uint16_t).
 * @param he_mesh Half-edge connectivity to walk.
 * @param start_idx Half-edge index at which the orbit begins.
 * @param visitor Invoked once per visited half-edge with its index.
 * @return Number of half-edges visited in the orbit.
 */
template <char OrbitMode, typename VisitorFn>
inline int vertex_orbit(const HalfEdgeMesh &he_mesh, uint16_t start_idx,
                        VisitorFn &&visitor) {
  uint16_t curr_idx = start_idx;
  int count = 0;
  // Anti-hang guard: a corrupt/non-manifold half-edge graph would otherwise spin forever.
  const int max_orbit = static_cast<int>(he_mesh.half_edges.size());
  do {
    HS_CHECK(count < max_orbit, "vertex_orbit: corrupt orbit");
    const HalfEdge &curr_he = he_mesh.half_edges[curr_idx];
    if (curr_he.face == HE_NONE)
      break;

    visitor(curr_idx);
    count++;

    if constexpr (OrbitMode == 'P') {
      // prev->pair orbit
      if (he_mesh.half_edges[curr_he.prev].pair == HE_NONE)
        break;
      curr_idx = he_mesh.half_edges[curr_he.prev].pair;
    } else {
      // pair->next orbit
      if (curr_he.pair == HE_NONE)
        break;
      curr_idx = he_mesh.half_edges[curr_he.pair].next;
    }
  } while (curr_idx != HE_NONE && curr_idx != start_idx);
  return count;
}

/**
 * @brief Emit one output face per source-vertex orbit (the shared dual/ambo/
 *   expand/snub scaffold).
 * @tparam DIR Orbit direction passed through to vertex_orbit ('P' or 'N').
 * @tparam ValueFn Callable mapping a half-edge index to an output vertex index.
 * @param he_mesh Half-edge connectivity to walk.
 * @param out_mesh Destination mesh; one face is appended per source vertex.
 * @param visited_verts Caller-allocated bool[V] scratch (filled here).
 * @param orbit_buf Caller-allocated uint16_t[I] scratch for the current orbit.
 * @param V Vertex count (size of visited_verts).
 * @param I Half-edge/index count (size of orbit_buf).
 * @param reverse Flip the winding for pair->next ('N') orbits whose natural
 *   order is opposite the desired front.
 * @param value_of Maps each visited half-edge index to its output vertex index.
 * @details Walks every half-edge once, and for each not-yet-visited origin
 *   vertex builds its orbit via vertex_orbit<DIR>, mapping each visited
 *   half-edge through value_of, then emits the collected indices as one face
 *   (skipping degenerate <3-gons).
 */
template <char DIR, typename ValueFn>
inline void emit_vertex_orbit_faces(const HalfEdgeMesh &he_mesh,
                                    PolyMesh &out_mesh, bool *visited_verts,
                                    uint16_t *orbit_buf, size_t V, size_t I,
                                    bool reverse, ValueFn &&value_of) {
  std::fill_n(visited_verts, V, false);
  for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
    uint16_t he_start_idx = static_cast<uint16_t>(i);
    const HalfEdge &he_start = he_mesh.half_edges[he_start_idx];

    uint16_t origin_idx = he_mesh.half_edges[he_start.prev].vertex;
    if (visited_verts[origin_idx])
      continue;
    visited_verts[origin_idx] = true;

    int orbit_count = 0;
    vertex_orbit<DIR>(he_mesh, he_start_idx, [&](uint16_t idx) {
      HS_CHECK(orbit_count < (int)I, "vertex orbit overflow");
      orbit_buf[orbit_count++] = narrow_index(value_of(idx));
    });

    if (orbit_count >= 3) {
      out_mesh.face_counts.push_back(narrow_face_count(orbit_count));
      if (reverse) {
        for (int k = orbit_count - 1; k >= 0; --k)
          out_mesh.faces.push_back(orbit_buf[k]);
      } else {
        for (int k = 0; k < orbit_count; ++k)
          out_mesh.faces.push_back(orbit_buf[k]);
      }
    }
  }
}

/**
 * @brief Count the sides of a face by walking its half-edge loop.
 * @param he_mesh Half-edge connectivity describing the face loop.
 * @param start_he First half-edge of the face, or HE_NONE for an empty face.
 * @return Number of sides; 0 for an empty face (start_he == HE_NONE).
 * @details The shared "how many sides does this face have" walk used by ambo/
 *   truncate before they re-walk to emit; matches the counting half of
 *   face_centroid without the position accumulation.
 */
inline int face_side_count(const HalfEdgeMesh &he_mesh, uint16_t start_he) {
  int count = 0;
  if (start_he != HE_NONE) {
    // Anti-hang guard: corrupt topology would otherwise spin forever.
    const int max_sides = static_cast<int>(he_mesh.half_edges.size());
    uint16_t he_idx = start_he;
    do {
      HS_CHECK(count < max_sides, "face_side_count: corrupt face loop");
      count++;
      he_idx = he_mesh.half_edges[he_idx].next;
    } while (he_idx != HE_NONE && he_idx != start_he);
  }
  return count;
}

/**
 * @brief Walk every undirected interior edge once, invoking visitor(he_idx, he)
 *   for each edge that has a pair.
 * @tparam VisitorFn Callable accepting (uint16_t he_idx, const HalfEdge &he).
 * @param he_mesh Half-edge connectivity to walk.
 * @param visited_edges Caller-allocated bool[I] scratch (filled here).
 * @param I Half-edge count (size of visited_edges).
 * @param visitor Invoked once per undirected interior edge.
 * @details The shared edge-dedup scaffold behind expand/chamfer/snub's "new
 *   face per edge" pass: each half-edge and its pair are marked visited so an
 *   undirected edge is emitted exactly once, and boundary edges (no pair) are
 *   skipped, matching the callers that never produced an edge face for them.
 */
template <typename VisitorFn>
inline void for_each_edge(const HalfEdgeMesh &he_mesh, bool *visited_edges,
                          size_t I, VisitorFn &&visitor) {
  std::fill_n(visited_edges, I, false);
  for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
    uint16_t he_idx = static_cast<uint16_t>(i);
    const HalfEdge &he = he_mesh.half_edges[he_idx];

    if (visited_edges[he_idx])
      continue;

    visited_edges[he_idx] = true;
    if (he.pair == HE_NONE)
      continue;
    visited_edges[he.pair] = true;

    visitor(he_idx, he);
  }
}

/**
 * @brief Copies a mesh's topology (as views) and vertices into a target mesh.
 * @param local_state The source mesh state (centered around origin).
 * @param world_state The destination mesh state to populate.
 * @param arena The memory arena from which to allocate vertices.
 * @details Base case of the variadic overload below: with no transformers,
 *   vertices are copied verbatim.
 * @note Mesh utility, not a Conway operator.
 * @note Input must be owned-mode (its topology is exposed as borrowed views in
 *   the output). The output is borrowed-mode, so it cannot be re-fed as input:
 *   transform(transform(x)) is unsupported and traps at the owned-mode check.
 */
inline void transform(const MeshState &local_state, MeshState &world_state,
                      Arena& arena) {
  HS_CHECK(local_state.face_counts.is_bound(),
           "MeshOps::transform: source mesh must be owned-mode");
  // Unbind any stale owned topology so the borrowed-mode views set below are not
  // shadowed by it on a reused output.
  world_state.face_counts = {};
  world_state.faces = {};
  world_state.face_offsets = {};

  world_state.face_counts_view = ArenaSpan(local_state.face_counts);
  world_state.faces_view = ArenaSpan(local_state.faces);
  world_state.face_offsets_view = ArenaSpan(local_state.face_offsets);

  world_state.vertices.bind(arena, local_state.vertices.size());

  world_state.vertices.append_bulk(local_state.vertices.data(),
                                   local_state.vertices.size());
}

/**
 * @brief Copies a mesh's topology (as views) and vertices into a target mesh,
 *   applying each transformer to every vertex in order.
 * @tparam T1 Type of the first vertex transformer.
 * @tparam Transformers Types of the remaining vertex transformers.
 * @param mesh The source mesh state.
 * @param transformed The destination mesh state to populate.
 * @param arena The memory arena from which to allocate vertices.
 * @param first_transformer First vertex transformer, applied first.
 * @param transformers Remaining vertex transformers, applied left to right and
 *   unrolled at compile time via a fold.
 */
template <typename T1, typename... Transformers>
inline void transform(const MeshState &mesh, MeshState &transformed, Arena& arena,
                      const T1 &first_transformer,
                      const Transformers &...transformers) {
  HS_CHECK(mesh.face_counts.is_bound(),
           "MeshOps::transform: source mesh must be owned-mode");
  // Unbind any stale owned topology (see the base overload).
  transformed.face_counts = {};
  transformed.faces = {};
  transformed.face_offsets = {};

  transformed.face_counts_view = ArenaSpan(mesh.face_counts);
  transformed.faces_view = ArenaSpan(mesh.faces);
  transformed.face_offsets_view = ArenaSpan(mesh.face_offsets);

  transformed.vertices.bind(arena, mesh.vertices.size());

  for (size_t i = 0; i < mesh.vertices.size(); ++i) {
    Vector v = mesh.vertices[i];

    v = first_transformer(v);

    (..., (v = transformers(v)));

    transformed.vertices.push_back(v);
  }
}

// ---------------------------------------------------------------------------
// Conway operators
//
// All operators take a const PolyMesh& input. PRIMITIVE operators (dual,
// kis, ambo, truncate, expand, chamfer, snub, relax) return a fresh PolyMesh in
// `target`. Both arenas are checkpointed via ScratchScope, so all scratch (the
// HalfEdgeMesh build plus, for most primitives, the per-orbit index/flag buffers
// — kis and relax allocate none; see the SCRATCH ARENA CONTRACT below) is
// reclaimed when the operator returns — only the output mesh persists. COMPOSED
// operators
// return in `temp`, not `target`; see COMPOSITION POLARITY below.
//
// SCRATCH ARENA CONTRACT (load-bearing — do not "standardize" blindly):
// The HalfEdgeMesh always builds in `temp`. The per-orbit index/flag buffers,
// however, are deliberately split:
//   - dual / ambo / truncate / expand  -> index buffers in `target`
//   - chamfer / snub                   -> index buffers in `temp`
//   - kis / relax                      -> no extra index buffers
// This is NOT drift. SolidBuilder ping-pongs the two arenas (target/temp swap
// each op) WITHOUT resetting between ops, so every intermediate mesh accumulates
// until the chain ends — and the two arenas can be small and ASYMMETRIC
// (HankinSolids runs the whole chain through a 16 KB / 32 KB pair). At the peak
// op `temp` already carries the live input mesh + the HalfEdgeMesh, so placing
// the index buffers in `target` SPLITS the transient load across both arenas
// instead of piling it onto the already-loaded side. Moving a buffer between
// arenas shifts that arena's high-water mark and can overflow the tight budget,
// so any change here must be measured against the configured arena split, not
// applied for uniformity.
//
// MANIFOLD PRECONDITION (per operator):
//   - dual / ambo / truncate / expand / chamfer / snub -> require a closed
//     manifold; require_closed_manifold traps on a boundary mesh.
//   - kis                                              -> per-face, no manifold
//     requirement.
//   - relax                                            -> tolerates a boundary
//     mesh (skips unpaired-twin vertices, partial relaxation).
//
// COMPOSITION POLARITY (load-bearing — do not "fix" for contract uniformity):
// Composed operators (gyro, meta, needle, zip, bevel) are written
// as op2(op1(mesh, target, temp), temp, target) so they reuse the SAME ping-pong
// internally and need NO extra arena: op1 writes to `target`, then op2 reads
// from `target` and writes to the swapped side. The consequence is that a
// two-op composition lands its output in `temp`, the OPPOSITE arena from a
// primitive — this is intended, not a contract slip. SolidBuilder swaps once per
// op (solids.h), which assumes the primitive polarity, so the single op
// FOLLOWING a composed op runs with its input and output on the same arena
// (no asymmetric split) for that one step before alternation self-restores. In
// the shipping recipes this happens once — the relax() after bevel() in
// IslamicStarPatterns::cube_relax_bevel33_relax_hk675_expand5 — and is
// covered by a high-water regression test
// (test_islamic_recipes_fit_islamicstars_budget in tests/test_solids.h) that
// asserts every Islamic recipe fits the 120 KB / 120 KB scratch split
// IslamicStars actually ships them through. Making SolidBuilder skip the swap
// after a composed op would restore healthy alternation but RELOCATES every
// downstream allocation between the two arenas, so re-run that test before
// adopting, not applied blindly.
//
// Per-vertex orbit construction uses an arena scratch buffer sized to the
// maximum possible valence (= total half-edges), so high-valence vertices
// never overflow.
// ---------------------------------------------------------------------------

/**
 * @brief Normalizes all vertices in the mesh to the unit sphere.
 * @tparam MeshT Mesh type exposing an iterable `vertices` of Vector.
 * @param mesh Mesh whose vertices are normalized in place.
 * @note Mesh utility, not a Conway operator.
 * @note Traps on a zero-length vertex; guard centroid-derived inputs with
 *   normalized_or (e.g. a centrally-symmetric face centroid; see kis/dual).
 */
template <typename MeshT> static void normalize(MeshT &mesh) {
  for (auto &v : mesh.vertices) {
    v.normalize();
  }
}

/**
 * @brief Computes the dual of a mesh (each face becomes a vertex and vice
 *   versa).
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh and its index scratch.
 * @param temp Arena holding the transient HalfEdgeMesh.
 * @return Fresh dual PolyMesh allocated in `target`.
 */
HS_COLD static PolyMesh dual(const PolyMesh &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();

  out_mesh.vertices.bind(target, F);
  out_mesh.face_counts.bind(target, V);
  out_mesh.faces.bind(target, I);

  {
    ScratchScope temp_guard(temp);
    ScratchScope target_guard(target);

    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "dual");

    for (size_t i = 0; i < he_mesh.faces.size(); ++i) {
      int count;
      Vector c = face_centroid(he_mesh, mesh, i, count);
      // Fall back to the face's first vertex on a zero-length (centrally-symmetric)
      // centroid, where strict normalized() would trap.
      Vector first_v =
          mesh.vertices[he_mesh.half_edges[he_mesh.faces[i].half_edge].vertex];
      out_mesh.vertices.push_back(normalized_or(c, first_v));
    }

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    // Per-orbit scratch buffer sized to the absolute upper bound on valence.
    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    // One face per source vertex: its orbit's incident faces become the dual
    // face's vertices.
    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) { return he_mesh.half_edges[idx].face; });
  }
  return out_mesh;
}

/**
 * @brief Kis operator: raises a pyramid on each face (apex at the face
 *   centroid).
 * @param mesh Source mesh.
 * @param target Arena receiving the output mesh.
 * @param temp Arena holding the transient scratch.
 * @return Fresh kis PolyMesh allocated in `target`.
 */
HS_COLD static PolyMesh kis(const PolyMesh &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  const auto *face_counts = mesh.get_face_counts_data();
  const auto *faces = mesh.get_faces_data();

  out_mesh.vertices.bind(target, V + F);
  out_mesh.face_counts.bind(target, I);
  out_mesh.faces.bind(target, 3 * I);

  {
    ScratchScope temp_guard(temp);

    for (size_t i = 0; i < V; ++i)
      out_mesh.vertices.push_back(mesh.vertices[i]);

    size_t offset = 0;
    for (size_t fi = 0; fi < F; ++fi) {
      int count = face_counts[fi];
      HS_CHECK(count > 0, "kis: zero-vertex face");
      Vector centroid(0, 0, 0);
      for (int k = 0; k < count; ++k) {
        centroid = centroid + mesh.vertices[faces[offset + k]];
      }
      centroid = centroid / static_cast<float>(count);

      out_mesh.vertices.push_back(normalized_or(centroid, mesh.vertices[faces[offset]]));
      int center_idx = narrow_index(out_mesh.vertices.size() - 1);

      for (int k = 0; k < count; ++k) {
        uint16_t vi = faces[offset + k];
        uint16_t vj = faces[offset + (k + 1) % count];
        out_mesh.face_counts.push_back(3);
        out_mesh.faces.push_back(vi);
        out_mesh.faces.push_back(vj);
        out_mesh.faces.push_back(static_cast<uint16_t>(center_idx));
      }
      offset += count;
    }
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Ambo operator: truncates vertices to edge midpoints (rectification).
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh and its index scratch.
 * @param temp Arena holding the transient HalfEdgeMesh.
 * @return Fresh ambo PolyMesh allocated in `target`.
 */
HS_COLD static PolyMesh ambo(const PolyMesh &mesh, Arena &target, Arena &temp) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, E);
  out_mesh.face_counts.bind(target, F + V);
  out_mesh.faces.bind(target, 2 * I);

  {
    ScratchScope temp_guard(temp);
    ScratchScope target_guard(target);

    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "ambo");

    uint16_t *edge_to_vert =
        static_cast<uint16_t *>(target.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(edge_to_vert, I, HE_NONE);

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    // 3. Populate Vertices
    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      if (edge_to_vert[i] == HE_NONE) {
        const HalfEdge &he = he_mesh.half_edges[i];
        uint16_t v1 = he.vertex;
        uint16_t v2 = he_mesh.half_edges[he.prev].vertex;

        Vector mid = (mesh.vertices[v1] + mesh.vertices[v2]) * 0.5f;
        out_mesh.vertices.push_back(normalized_or(mid, mesh.vertices[v1]));
        uint16_t new_idx = narrow_index(out_mesh.vertices.size() - 1);

        edge_to_vert[i] = new_idx;
        if (he.pair != HE_NONE)
          edge_to_vert[he.pair] = new_idx;
      }
    }

    // 4. Reconstruct Original Faces (Shrunk)
    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      uint16_t start = he_mesh.faces[fi].half_edge;
      int count = face_side_count(he_mesh, start);
      if (count >= 3) {
        out_mesh.face_counts.push_back(narrow_face_count(count));
        uint16_t he_idx = start;
        int emitted = 0; // anti-hang guard: re-walk emits exactly `count` sides
        do {
          HS_CHECK(emitted++ < count, "face re-walk overran side count");
          out_mesh.faces.push_back(edge_to_vert[he_idx]);
          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    // 5. Build Vertex Orbits (New Faces)
    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) { return edge_to_vert[idx]; });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Recover a truncate edge's two cut vertices in half-edge walk order.
 * @param he_mesh Half-edge connectivity describing the edge.
 * @param edge_to_vert Per-half-edge table of the edge's two cut-vertex indices,
 *   stored canonically as {near-k1, near-k2} with k1 = min(endpoint indices).
 * @param he_idx Half-edge (tail->head) whose oriented cut pair is recovered.
 * @return {tail-side cut, head-side cut}; the stored order is reversed when the
 *   tail is not the canonical k1.
 * @details Centralizes winding so every caller emits a consistent order instead
 *   of re-deriving the `vi==k1` test inline.
 */
inline std::pair<uint16_t, uint16_t> truncate_oriented_cut(
    const HalfEdgeMesh &he_mesh,
    const std::pair<uint16_t, uint16_t> *edge_to_vert, uint16_t he_idx) {
  const HalfEdge &he = he_mesh.half_edges[he_idx];
  uint16_t tail = he_mesh.half_edges[he.prev].vertex;
  uint16_t head = he.vertex;
  const std::pair<uint16_t, uint16_t> &cut = edge_to_vert[he_idx];
  return (tail <= head) ? std::pair<uint16_t, uint16_t>{cut.first, cut.second}
                        : std::pair<uint16_t, uint16_t>{cut.second, cut.first};
}

/**
 * @brief Truncate operator: cuts corners off the polyhedron.
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh and its index scratch.
 * @param temp Arena holding the transient HalfEdgeMesh.
 * @param t Truncation depth, the fraction along each edge at which the two cut
 *   points sit, in [0..1]. Each edge `(k1,k2)` yields `k1+(k2-k1)*t` and
 *   `k2+(k1-k2)*t`. For `t<0.5` the cut points stay on their own half; at
 *   exactly `0.5` both reach the midpoint and this short-circuits to `ambo`;
 *   for `t>0.5` the two points cross past each other, producing intentional
 *   self-intersecting cut faces (used by the `*_truncate50d_*` solids). `t`
 *   outside `[0..1]` would place a cut point beyond the edge endpoints, so it
 *   traps per the fail-fast doctrine.
 * @return Fresh truncated PolyMesh allocated in `target` (or the ambo result
 *   when t == 0.5).
 */
HS_COLD static PolyMesh truncate(const PolyMesh &mesh, Arena &target, Arena &temp,
                                  float t = 0.25f) {
  HS_CHECK(t >= 0.0f && t <= 1.0f, "truncate: t out of [0,1]");
  if (std::abs(t - 0.5f) < math::TOLERANCE) {
    return ambo(mesh, target, temp);
  }

  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, 2 * E);
  out_mesh.face_counts.bind(target, F + V);
  out_mesh.faces.bind(target, 3 * I);

  {
    ScratchScope temp_guard(temp);
    ScratchScope target_guard(target);

    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "truncate");

    // Per-edge cut-vertex pair, unset = {HE_NONE, HE_NONE} (matches the other
    // operators' he->new-vertex maps so the unset test reads identically).
    std::pair<uint16_t, uint16_t> *edge_to_vert =
        static_cast<std::pair<uint16_t, uint16_t> *>(
            target.allocate(
                I * sizeof(std::pair<uint16_t, uint16_t>),
                alignof(std::pair<uint16_t, uint16_t>)));
    std::fill_n(edge_to_vert, I, std::pair<uint16_t, uint16_t>(HE_NONE, HE_NONE));

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
      if (edge_to_vert[i].first == HE_NONE) {
        const HalfEdge &he = he_mesh.half_edges[i];
        uint16_t vi = he_mesh.half_edges[he.prev].vertex;
        uint16_t vj = he.vertex;
        uint16_t k1 = std::min(vi, vj);
        uint16_t k2 = std::max(vi, vj);

        Vector new_u =
            mesh.vertices[k1] + (mesh.vertices[k2] - mesh.vertices[k1]) * t;
        Vector new_v =
            mesh.vertices[k2] + (mesh.vertices[k1] - mesh.vertices[k2]) * t;

        out_mesh.vertices.push_back(new_u);
        uint16_t idx_u = narrow_index(out_mesh.vertices.size() - 1);

        out_mesh.vertices.push_back(new_v);
        uint16_t idx_v = narrow_index(out_mesh.vertices.size() - 1);

        edge_to_vert[i] = {idx_u, idx_v};
        if (he.pair != HE_NONE)
          edge_to_vert[he.pair] = {idx_u, idx_v};
      }
    }

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      uint16_t start = he_mesh.faces[fi].half_edge;
      int count = face_side_count(he_mesh, start);
      if (count >= 3) {
        out_mesh.face_counts.push_back(narrow_face_count(count * 2));
        uint16_t he_idx = start;
        int emitted = 0; // anti-hang guard: re-walk emits exactly `count` sides
        do {
          HS_CHECK(emitted++ < count, "face re-walk overran side count");
          auto [tail_cut, head_cut] =
              truncate_oriented_cut(he_mesh, edge_to_vert, he_idx);
          out_mesh.faces.push_back(tail_cut);
          out_mesh.faces.push_back(head_cut);

          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    emit_vertex_orbit_faces<'P'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/false,
        [&](uint16_t idx) {
          // Orbit walks the half-edges leaving this vertex; each contributes
          // its tail-side cut (the cut point nearest the shared vertex).
          return truncate_oriented_cut(he_mesh, edge_to_vert, idx).first;
        });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Expand operator: separates faces (e = aa).
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh and its index scratch.
 * @param temp Arena holding the transient HalfEdgeMesh.
 * @param t Expansion factor. Default 2-sqrt(2) ~= 0.5857.
 * @return Fresh expanded PolyMesh allocated in `target`.
 */
HS_COLD static PolyMesh expand(const PolyMesh &mesh, Arena &target, Arena &temp,
                                float t = 2.0f - sqrtf(2.0f)) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, I);
  out_mesh.face_counts.bind(target, F + V + E);
  out_mesh.faces.bind(target, 4 * I);

  {
    ScratchScope temp_guard(temp);
    ScratchScope target_guard(target);

    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "expand");
    // he->new-vertex map, unset = HE_NONE (matches ambo/snub/chamfer so the
    // unset test reads identically across every operator).
    uint16_t *he_to_vert_idx = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_vert_idx, I, HE_NONE);

    bool *visited_verts = static_cast<bool *>(
        target.allocate(V * sizeof(bool), alignof(bool)));

    bool *visited_edges = static_cast<bool *>(
        target.allocate(I * sizeof(bool), alignof(bool)));

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        target.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      // Emit the shrunk primary face only when well-formed (>=3 sides); new
      // vertices and the he->vertex mapping are still created unconditionally.
      const bool well_formed = count >= 3;
      if (well_formed)
        out_mesh.face_counts.push_back(narrow_face_count(count));
      he_idx = start;
      if (he_idx != HE_NONE) {
        int walked = 0; // anti-hang guard: face has `count` half-edges
        do {
          HS_CHECK(walked++ < count, "face re-walk overran side count");
          Vector v = mesh.vertices[he_mesh.half_edges[he_idx].vertex];
          Vector new_v = v + (centroid - v) * t;
          out_mesh.vertices.push_back(new_v);
          uint16_t idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_vert_idx[he_idx] = idx;

          if (well_formed)
            out_mesh.faces.push_back(idx);
          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    emit_vertex_orbit_faces<'N'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/true,
        [&](uint16_t idx) {
          return he_to_vert_idx[he_mesh.half_edges[idx].prev];
        });

    for_each_edge(he_mesh, visited_edges, I,
                  [&](uint16_t he_idx, const HalfEdge &he) {
                    out_mesh.face_counts.push_back(4);
                    out_mesh.faces.push_back(he_to_vert_idx[he.prev]);
                    out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
                    out_mesh.faces.push_back(
                        he_to_vert_idx[he_mesh.half_edges[he.pair].prev]);
                    out_mesh.faces.push_back(he_to_vert_idx[he_idx]);
                  });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Chamfer operator: replaces edges with hexagonal faces.
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh.
 * @param temp Arena holding the transient HalfEdgeMesh and index scratch.
 * @param t Thickness factor for the new hexagons [0..1].
 * @return Fresh chamfered PolyMesh allocated in `target`.
 */
HS_COLD static PolyMesh chamfer(const PolyMesh &mesh, Arena &target, Arena &temp,
                                 float t = 0.5f) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, V + 2 * E);
  out_mesh.face_counts.bind(target, F + E);
  out_mesh.faces.bind(target, I + 6 * E);

  {
    ScratchScope temp_guard(temp);
    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "chamfer");

    uint16_t *he_to_new_v =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_new_v, I, HE_NONE);

    // 1. Copy original vertices
    for (size_t i = 0; i < V; ++i) {
      out_mesh.vertices.push_back(mesh.vertices[i]);
    }

    // 2. Generate new vertices and shrunk faces
    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      // Emit the shrunk primary face only when well-formed (>=3 sides); new
      // vertices and the he->vertex mapping are created unconditionally.
      const bool well_formed = count >= 3;
      if (well_formed)
        out_mesh.face_counts.push_back(narrow_face_count(count));
      he_idx = start;

      if (he_idx != HE_NONE) {
        int walked = 0; // anti-hang guard: face has `count` half-edges
        do {
          HS_CHECK(walked++ < count, "face re-walk overran side count");
          uint16_t vi =
              he_mesh.half_edges[he_mesh.half_edges[he_idx].prev].vertex;
          Vector v = mesh.vertices[vi];
          Vector new_v = v + (centroid - v) * t;

          out_mesh.vertices.push_back(new_v);
          uint16_t idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_new_v[he_idx] = idx;

          if (well_formed)
            out_mesh.faces.push_back(idx);

          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    // 3. Generate Hexagon faces for edges
    bool *visited_edges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));

    for_each_edge(he_mesh, visited_edges, I,
                  [&](uint16_t he_idx, const HalfEdge &he) {
                    out_mesh.face_counts.push_back(6);

                    uint16_t A = he_mesh.half_edges[he.prev].vertex;
                    uint16_t B = he.vertex;

                    out_mesh.faces.push_back(A);
                    out_mesh.faces.push_back(
                        he_to_new_v[he_mesh.half_edges[he.pair].next]);
                    out_mesh.faces.push_back(he_to_new_v[he.pair]);
                    out_mesh.faces.push_back(B);
                    out_mesh.faces.push_back(he_to_new_v[he.next]);
                    out_mesh.faces.push_back(he_to_new_v[he_idx]);
                  });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Edge-length relaxation by spring forces on the unit sphere.
 * @param mesh Source mesh, copied into the output before relaxing.
 * @param target Arena receiving the output mesh.
 * @param temp Arena holding the transient movement buffer and HalfEdgeMesh.
 * @param iterations Maximum spring-relaxation passes; stops early on
 *   convergence.
 * @return Fresh relaxed PolyMesh allocated in `target`.
 * @note Unlike its sibling operators, relax tolerates a boundary mesh: a vertex
 *   with no outgoing twin is skipped, yielding a partial relaxation instead of a
 *   closed-manifold trap.
 */
HS_COLD static PolyMesh relax(const PolyMesh &mesh, Arena &target, Arena &temp,
                               int iterations = 8) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  const auto *face_counts = mesh.get_face_counts_data();
  const auto *faces = mesh.get_faces_data();

  out_mesh.vertices.bind(target, V);
  out_mesh.face_counts.bind(target, F);
  out_mesh.faces.bind(target, I);

  for (size_t i = 0; i < V; ++i)
    out_mesh.vertices.push_back(mesh.vertices[i]);
  for (size_t i = 0; i < F; ++i) {
    HS_CHECK(face_counts[i] >= 3, "relax: degenerate face (< 3 sides)");
    out_mesh.face_counts.push_back(face_counts[i]);
  }
  for (size_t i = 0; i < I; ++i)
    out_mesh.faces.push_back(faces[i]);

  {
    ScratchScope temp_guard(temp);

    ArenaVector<Vector> movements;
    movements.bind(temp, V);
    for (size_t i = 0; i < V; ++i)
      movements.push_back(Vector(0, 0, 0));

    HalfEdgeMesh he_mesh(temp, out_mesh);

    // Per-vertex orbit start: an outgoing half-edge (pair of an interior
    // incoming edge). vertices[i].half_edge holds only the last incoming edge
    // written; if that one is a boundary its 1-ring would be dropped even when a
    // paired incoming edge exists, so scan for one.
    ArenaVector<uint16_t> orbit_start;
    orbit_start.bind(temp, V);
    for (size_t i = 0; i < V; ++i)
      orbit_start.push_back(HE_NONE);
    for (size_t e = 0; e < he_mesh.half_edges.size(); ++e) {
      const HalfEdge &he = he_mesh.half_edges[e];
      if (he.pair != HE_NONE && orbit_start[he.vertex] == HE_NONE)
        orbit_start[he.vertex] = he.pair;
    }

    for (int iter = 0; iter < iterations; ++iter) {
      float total_len = 0;
      int edge_count = 0;

      for (size_t i = 0; i < he_mesh.half_edges.size(); ++i) {
        const HalfEdge &he = he_mesh.half_edges[i];
        int u = he_mesh.half_edges[he.prev].vertex;
        int v = he.vertex;
        if (u < v) {
          total_len +=
              distance_between(out_mesh.vertices[u], out_mesh.vertices[v]);
          edge_count++;
        }
      }

      if (edge_count == 0)
        break;
      float target_len = total_len / edge_count;

      for (size_t i = 0; i < V; ++i) {
        Vector force(0, 0, 0);

        uint16_t start_out = orbit_start[i]; // outgoing edge; HE_NONE if boundary
        if (start_out != HE_NONE) {
          vertex_orbit<'N'>(he_mesh, start_out, [&](uint16_t idx) {
            int ni = he_mesh.half_edges[idx].vertex; // head == 1-ring neighbor
            Vector vec = out_mesh.vertices[ni] - out_mesh.vertices[i];
            // Skip near-zero edges before 1/dist can spike the force.
            float len_sq = dot(vec, vec);
            if (len_sq > math::EPS_LEN_SQ) {
              float dist = sqrtf(len_sq);
              float diff = dist - target_len;
              force = force + (vec * (1.0f / dist)) * (diff * 0.1f);
            }
          });
        }
        movements[i] = force;
      }

      float max_move_sq = 0.0f;
      for (size_t i = 0; i < V; ++i) {
        max_move_sq = std::max(max_move_sq, dot(movements[i], movements[i]));
        out_mesh.vertices[i] =
            (out_mesh.vertices[i] + movements[i]).normalized();
      }

      // Stop early once the largest per-vertex spring force settles below ~1e-4
      // (a sub-pixel angular step on the unit sphere).
      constexpr float RELAX_CONVERGE_EPS_SQ = 1e-8f;
      if (max_move_sq < RELAX_CONVERGE_EPS_SQ)
        break;
    }
  }
  // iterations == 0 skips the in-loop sphere projection; normalize so the
  // pass-through still lands on the unit sphere (idempotent otherwise).
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Snub operator: creates a chiral semi-regular polyhedron.
 * @param mesh Source mesh; must be a closed manifold.
 * @param target Arena receiving the output mesh.
 * @param temp Arena holding the transient HalfEdgeMesh and index scratch.
 * @param t Inset factor of each face toward its centroid [0..1].
 * @param twist Per-face rotation about the face normal, in radians; 0 disables
 *   the twist pass.
 * @return Fresh snub PolyMesh allocated in `target`.
 * @details Uses Newell's method for face normals, robust to non-planar faces on
 *   the unit sphere and to collinear vertex triplets.
 */
HS_COLD static PolyMesh snub(const PolyMesh &mesh, Arena &target, Arena &temp,
                              float t = 0.5f, float twist = 0.0f) {
  PolyMesh out_mesh;
  size_t V = mesh.vertices.size();
  size_t F = mesh.get_face_counts_size();
  size_t I = mesh.get_faces_size();
  size_t E = I / 2;

  out_mesh.vertices.bind(target, I);
  out_mesh.face_counts.bind(target, F + V + 2 * E);
  out_mesh.faces.bind(target, 5 * I);

  {
    ScratchScope temp_guard(temp);

    HalfEdgeMesh he_mesh(temp, mesh);
    require_closed_manifold(he_mesh, "snub");
    uint16_t *he_to_vert_idx =
        static_cast<uint16_t *>(temp.allocate(
            I * sizeof(uint16_t), alignof(uint16_t)));
    std::fill_n(he_to_vert_idx, I, HE_NONE);

    uint16_t *orbit_buf = static_cast<uint16_t *>(
        temp.allocate(I * sizeof(uint16_t), alignof(uint16_t)));

    for (size_t fi = 0; fi < he_mesh.faces.size(); ++fi) {
      int count;
      Vector centroid = face_centroid(he_mesh, mesh, fi, count);
      uint16_t start = he_mesh.faces[fi].half_edge;
      uint16_t he_idx = start;

      // Newell's method face normal — robust for sphere-projected faces.
      Vector normal_raw = face_normal(he_mesh, mesh, fi);
      Vector normal(0, 0, 0);
      if (dot(normal_raw, normal_raw) > math::EPS_NORMAL_SQ) {
        normal = normal_raw.normalized();
      } else if (dot(centroid, centroid) > math::EPS_LEN_SQ) {
        normal = centroid.normalized();
      }

      const bool do_twist = twist != 0.0f;
      Quaternion twist_q;
      if (do_twist)
        twist_q = make_rotation(normal, twist);

      // Emit the twisted primary face only when well-formed (>=3 sides); new
      // vertices and the he->vertex mapping are created unconditionally.
      const bool well_formed = count >= 3;
      if (well_formed)
        out_mesh.face_counts.push_back(narrow_face_count(count));
      he_idx = start;
      if (he_idx != HE_NONE) {
        int walked = 0; // anti-hang guard: face has `count` half-edges
        do {
          HS_CHECK(walked++ < count, "face re-walk overran side count");
          Vector v = mesh.vertices[he_mesh.half_edges[he_idx].vertex];
          Vector new_v = v + (centroid - v) * t;

          if (do_twist) {
            Vector local = new_v - centroid;
            new_v = centroid + rotate(local, twist_q);
          }

          out_mesh.vertices.push_back(new_v);
          uint16_t idx = narrow_index(out_mesh.vertices.size() - 1);
          he_to_vert_idx[he_idx] = idx;

          if (well_formed)
            out_mesh.faces.push_back(idx);

          he_idx = he_mesh.half_edges[he_idx].next;
        } while (he_idx != HE_NONE && he_idx != start);
      }
    }

    bool *visited_verts = static_cast<bool *>(
        temp.allocate(V * sizeof(bool), alignof(bool)));
    emit_vertex_orbit_faces<'N'>(
        he_mesh, out_mesh, visited_verts, orbit_buf, V, I, /*reverse=*/true,
        [&](uint16_t idx) {
          return he_to_vert_idx[he_mesh.half_edges[idx].prev];
        });

    bool *visited_edges = static_cast<bool *>(
        temp.allocate(I * sizeof(bool), alignof(bool)));

    for_each_edge(he_mesh, visited_edges, I,
                  [&](uint16_t he_idx, const HalfEdge &he) {
                    out_mesh.face_counts.push_back(3);
                    out_mesh.faces.push_back(he_to_vert_idx[he.prev]);
                    out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
                    out_mesh.faces.push_back(he_to_vert_idx[he_idx]);

                    out_mesh.face_counts.push_back(3);
                    out_mesh.faces.push_back(he_to_vert_idx[he.pair]);
                    out_mesh.faces.push_back(
                        he_to_vert_idx[he_mesh.half_edges[he.pair].prev]);
                    out_mesh.faces.push_back(he_to_vert_idx[he_idx]);
                  });
  }
  normalize(out_mesh);
  return out_mesh;
}

/**
 * @brief Gyro operator: dual of snub (pentagonal chiral subdivision).
 * @param mesh Source mesh.
 * @param target Arena used as the ping-pong source for the composition.
 * @param temp Arena used as the ping-pong destination for the composition.
 * @return Composed PolyMesh; the output lands in `temp`, not `target` (see
 *   COMPOSITION POLARITY at the top of the operator block).
 */
HS_COLD static PolyMesh gyro(const PolyMesh &mesh, Arena &target, Arena &temp) {
  return dual(snub(mesh, target, temp), temp, target);
}

// ---------------------------------------------------------------------------
// Compositional operators (Hart's notation)
//
// These are defined as compositions of primitive operators. Equivalences are
// from Hart's reference implementation. Memory-wise, each composition reuses
// the standard ping-pong of (target, temp) arenas with no extra allocation, and
// as a result returns its output in `temp`, not `target` (see COMPOSITION
// POLARITY at the top of the operator block).
//   meta   m = kj = kis of ambo (j = a)
//   needle n = kd = kis of dual
//   zip    z = dk = dual of kis (truncated dual)
//   bevel  b = ta = truncate of ambo (rectify-then-truncate)
// ---------------------------------------------------------------------------

/**
 * @brief Meta operator (Hart's `m`): kis of ambo (m = kj, j = a).
 * @param mesh Source mesh.
 * @param target Arena used as the ping-pong source for the composition.
 * @param temp Arena used as the ping-pong destination for the composition.
 * @return Composed PolyMesh; the output lands in `temp` (see COMPOSITION
 *   POLARITY at the top of the operator block).
 */
HS_COLD static PolyMesh meta(const PolyMesh &mesh, Arena &target, Arena &temp) {
  return kis(ambo(mesh, target, temp), temp, target);
}

/**
 * @brief Needle operator: kis of the dual (n = kd).
 * @param mesh Source mesh.
 * @param target Arena used as the ping-pong source for the composition.
 * @param temp Arena used as the ping-pong destination for the composition.
 * @return Composed PolyMesh; the output lands in `temp` (see COMPOSITION
 *   POLARITY at the top of the operator block).
 */
HS_COLD static PolyMesh needle(const PolyMesh &mesh, Arena &target, Arena &temp) {
  return kis(dual(mesh, target, temp), temp, target);
}

/**
 * @brief Zip operator: dual of kis, i.e. the truncated dual (z = dk).
 * @param mesh Source mesh.
 * @param target Arena used as the ping-pong source for the composition.
 * @param temp Arena used as the ping-pong destination for the composition.
 * @return Composed PolyMesh; the output lands in `temp` (see COMPOSITION
 *   POLARITY at the top of the operator block).
 */
HS_COLD static PolyMesh zip(const PolyMesh &mesh, Arena &target, Arena &temp) {
  return dual(kis(mesh, target, temp), temp, target);
}

/**
 * @brief Bevel operator: truncate of ambo (b = ta).
 * @param mesh Source mesh.
 * @param target Arena used as the ping-pong source for the composition.
 * @param temp Arena used as the ping-pong destination for the composition.
 * @param t Truncation depth forwarded to the truncate step, in [0..1].
 * @return Composed PolyMesh; the output lands in `temp` (see COMPOSITION
 *   POLARITY at the top of the operator block).
 */
HS_COLD static PolyMesh bevel(const PolyMesh &mesh, Arena &target, Arena &temp,
                               float t = 0.25f) {
  return truncate(ambo(mesh, target, temp), temp, target, t);
}

// TODO: Propeller (Hart's `p`) and whirl/loft are not implemented; their chiral
// blade winding on the sphere needs a dedicated test before adding.

} // namespace MeshOps
