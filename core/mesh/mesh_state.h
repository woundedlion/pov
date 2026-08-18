/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file mesh_state.h
 * @brief MeshState, the arena-backed flat mesh format the renderer and the mesh
 *        operators pass around, plus the shared arena copy helper.
 */

#include <cstdint>
#include <utility>

#include "math/3dmath.h"
#include "engine/memory.h"
#include "platform/platform.h"

namespace MeshOps {

/**
 * @brief Deep-copies n contiguous elements into a freshly-bound ArenaVector.
 * @tparam T Element type.
 * @param dst Destination vector, bound to exactly n elements from arena.
 * @param src Pointer to the first source element (ignored when n == 0).
 * @param n Number of elements to copy.
 * @param arena Arena supplying storage for dst.
 * @details Shared by MeshState::clone, MeshOps::clone and
 * CompiledHankin::clone.
 */
template <typename T>
inline void copy_vector(ArenaVector<T> &dst, const T *src, size_t n,
                        Arena &arena) {
  dst.bind(arena, n);
  dst.append_bulk(src, n);
}

} // namespace MeshOps

/**
 * @brief Represents the state of a mesh using arena storage to avoid heap
 * allocations.
 */
struct MeshState {
  ArenaVector<Vector> vertices; /**< Owned vertex positions, in world units. */
  ArenaVector<uint8_t> face_counts; /**< Owned per-face vertex counts. */
  ArenaVector<uint16_t> faces;      /**< Owned flattened face vertex indices. */
  ArenaVector<uint16_t>
      face_offsets; /**< Owned start offset of each face into faces. */
  ArenaVector<uint16_t> topology; /**< Optional owned per-face topology class
                   ids from classify_faces_by_topology; either empty or one
                   dense 16-bit-bounded id per face. */

  /**
   * @brief Constructs an empty, unbound mesh.
   */
  MeshState() = default;

  /**
   * @brief Move-constructs, transferring owned buffers and views.
   * @param other Source mesh; its borrowed views are cleared so it holds no
   * dangling borrows.
   */
  MeshState(MeshState &&other) noexcept
      : vertices(std::move(other.vertices)),
        face_counts(std::move(other.face_counts)),
        faces(std::move(other.faces)),
        face_offsets(std::move(other.face_offsets)),
        topology(std::move(other.topology)),
        face_counts_view(other.face_counts_view), faces_view(other.faces_view),
        face_offsets_view(other.face_offsets_view),
        topology_view(other.topology_view) {
    other.set_owned();
  }

  /**
   * @brief Move-assigns, transferring owned buffers and views.
   * @param other Source mesh; its borrowed views are cleared so the moved-from
   * mesh holds no dangling borrows.
   * @return Reference to this mesh.
   */
  MeshState &operator=(MeshState &&other) noexcept {
    if (this != &other) {
      vertices = std::move(other.vertices);
      face_counts = std::move(other.face_counts);
      faces = std::move(other.faces);
      face_offsets = std::move(other.face_offsets);
      topology = std::move(other.topology);
      face_counts_view = other.face_counts_view;
      faces_view = other.faces_view;
      face_offsets_view = other.face_offsets_view;
      topology_view = other.topology_view;
      other.set_owned();
    }
    return *this;
  }

  /**
   * @brief Resets to empty: clears owned buffers and drops the borrowed views.
   * @details Every arena binding survives, so this is the wrong reset for a
   * mesh whose arena has since been reset or rewound: bind()'s reuse path would
   * hand the mesh back blocks the arena has already reclaimed. Replace such a
   * mesh instead (`mesh = MeshState()`), as MeshOps::compile does.
   */
  void clear() {
    vertices.clear();
    face_counts.clear();
    faces.clear();
    face_offsets.clear();
    topology.clear();
    set_owned();
  }

  /**
   * @brief Checks whether the vertex buffer is bound (has been allocated).
   * @return True if the mesh owns vertex storage; says nothing about the
   *   topology arrays, which are unbound in borrowed mode.
   */
  bool is_bound() const { return vertices.is_bound(); }

  /**
   * @brief Returns the face-counts pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   * @details Discriminates on is_bound() (which mode this is), NOT on empty():
   * an owned-but-legitimately-empty mesh is bound with size 0, and gating on
   * empty() would wrongly fall through to a stale/unset borrowed view (same for
   * the sibling accessors below).
   */
  const uint8_t *get_face_counts_data() const {
    return face_counts.is_bound() ? face_counts.data()
                                  : face_counts_view.data();
  }
  /**
   * @brief Returns the number of face counts for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_face_counts_size() const {
    return face_counts.is_bound() ? face_counts.size()
                                  : face_counts_view.size();
  }

  /**
   * @brief Returns the faces pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   */
  const uint16_t *get_faces_data() const {
    return faces.is_bound() ? faces.data() : faces_view.data();
  }
  /**
   * @brief Returns the number of face indices for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_faces_size() const {
    return faces.is_bound() ? faces.size() : faces_view.size();
  }

  /**
   * @brief Returns the face-offsets pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   */
  const uint16_t *get_face_offsets_data() const {
    return face_offsets.is_bound() ? face_offsets.data()
                                   : face_offsets_view.data();
  }
  /**
   * @brief Returns the number of face offsets for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size.
   */
  size_t get_face_offsets_size() const {
    return face_offsets.is_bound() ? face_offsets.size()
                                   : face_offsets_view.size();
  }

  /**
   * @brief Returns the per-face topology-class pointer for the active mode.
   * @return Owned data in owned mode, otherwise the borrowed view pointer.
   *   Null when the mesh is unclassified.
   */
  const uint16_t *get_topology_data() const {
    return topology.is_bound() ? topology.data() : topology_view.data();
  }
  /**
   * @brief Returns the number of per-face topology classes for the active mode.
   * @return Owned size in owned mode, otherwise the borrowed view size; 0 when
   *   the mesh is unclassified.
   */
  size_t get_topology_size() const {
    return topology.is_bound() ? topology.size() : topology_view.size();
  }

  /**
   * @brief Returns the number of vertices in the mesh.
   * @return Vertex count.
   */
  size_t num_vertices() const { return vertices.size(); }
  /**
   * @brief Returns the number of faces in the mesh.
   * @return Face count (equal to the active face-counts size).
   */
  size_t num_faces() const { return get_face_counts_size(); }

  /**
   * @brief Deep-copies all owned data from src into dst using a target arena.
   * @param src Source mesh to copy from.
   * @param dst Destination mesh to populate.
   * @param arena Arena providing storage for the destination buffers.
   * @details Required by Cloneable. Traps if src aliases dst: the set_owned()
   * below would drop a borrowed src's views before they are read, yielding an
   * empty mesh.
   */
  HS_COLD_MEMBER static void clone(const MeshState &src, MeshState &dst,
                                   Arena &arena) {
    HS_CHECK(&src != &dst, "MeshState::clone src must not alias dst");
    // Reused dst may carry stale views; clone produces an owned-mode mesh.
    dst.set_owned();
    MeshOps::copy_vector(dst.vertices, src.vertices.data(), src.vertices.size(),
                         arena);
    MeshOps::copy_vector(dst.face_counts, src.get_face_counts_data(),
                         src.get_face_counts_size(), arena);
    MeshOps::copy_vector(dst.faces, src.get_faces_data(), src.get_faces_size(),
                         arena);
    MeshOps::copy_vector(dst.face_offsets, src.get_face_offsets_data(),
                         src.get_face_offsets_size(), arena);
    MeshOps::copy_vector(dst.topology, src.get_topology_data(),
                         src.get_topology_size(), arena);
  }

  /**
   * @brief Switches to owned mode by dropping the borrowed topology views so the
   *   is_bound() accessors read the owned buffers.
   * @details Call on a possibly-borrowed MeshState after (re)binding its owned
   *   topology, so a leftover view can't shadow the owned data.
   */
  void set_owned() {
    face_counts_view = {};
    faces_view = {};
    face_offsets_view = {};
    topology_view = {};
  }

  /**
   * @brief Reports whether every face offset is the running sum of the face
   *   counts before it.
   * @param face_counts_span Per-face vertex counts.
   * @param face_offsets_span Per-face start offsets into the flat faces list;
   *   any length other than one entry per face count reports false rather than
   *   indexing past the span.
   * @return True when offset[i] equals counts[0] + ... + counts[i-1] for every
   *   face.
   * @details Exact equality at each index, which also establishes that the
   *   offsets are non-decreasing and that no two face spans overlap. Endpoint
   *   agreement alone does not: interior offsets can be scrambled while the
   *   first and last still line up.
   */
  static bool offsets_are_prefix_sum(ArenaSpan<uint8_t> face_counts_span,
                                     ArenaSpan<uint16_t> face_offsets_span) {
    if (face_offsets_span.size() != face_counts_span.size())
      return false;
    size_t running = 0;
    for (size_t i = 0; i < face_counts_span.size(); ++i) {
      if (static_cast<size_t>(face_offsets_span[i]) != running)
        return false;
      running += face_counts_span[i];
    }
    return true;
  }

  /**
   * @brief Switches to borrowed mode: drops the owned face_counts/faces/
   *   face_offsets/topology arrays so they cannot shadow the views, then points
   *   the four views at the given spans.
   * @param face_counts_span Borrowed per-face vertex counts.
   * @param faces_span Borrowed flattened face vertex indices.
   * @param face_offsets_span Borrowed per-face start offsets into faces. Empty
   *   when the source mesh carries no offsets (only the solid scan path needs
   *   them); otherwise one entry per face.
   * @param topology_span Borrowed per-face topology class ids. Empty when the
   *   source mesh is unclassified; otherwise one entry per face.
   * @details Traps on inconsistent spans: a present offsets array must be one
   *   entry per face, and its last offset plus that face's count must cover the
   *   whole flat faces list. With no offsets the counts must sum to the flat
   *   faces length. A present topology array must be one entry per face. The
   *   interior offsets are audited against the prefix sum under HS_AUDIT_CHECK
   *   only: this runs per frame from MeshOps::transform, so the device pays
   *   nothing for the O(F) walk.
   */
  void set_view(ArenaSpan<uint8_t> face_counts_span,
                ArenaSpan<uint16_t> faces_span,
                ArenaSpan<uint16_t> face_offsets_span,
                ArenaSpan<uint16_t> topology_span = {}) {
    if (!face_offsets_span.is_empty()) {
      HS_CHECK(face_offsets_span.size() == face_counts_span.size(),
               "MeshState::set_view: one face offset per face count required");
      const size_t last = face_counts_span.size() - 1;
      HS_CHECK(static_cast<size_t>(face_offsets_span[last]) +
                       face_counts_span[last] ==
                   faces_span.size(),
               "MeshState::set_view: face offsets do not span faces");
      HS_AUDIT_CHECK(
          offsets_are_prefix_sum(face_counts_span, face_offsets_span),
          "MeshState::set_view: face offsets are not the prefix sum of the "
          "face counts");
    } else {
      size_t counted = 0;
      for (size_t i = 0; i < face_counts_span.size(); ++i)
        counted += face_counts_span[i];
      HS_CHECK(counted == faces_span.size(),
               "MeshState::set_view: face counts do not span faces");
    }
    HS_CHECK(topology_span.is_empty() ||
                 topology_span.size() == face_counts_span.size(),
             "MeshState::set_view: one topology class per face required");
    face_counts = {};
    faces = {};
    face_offsets = {};
    topology = {};
    face_counts_view = face_counts_span;
    faces_view = faces_span;
    face_offsets_view = face_offsets_span;
    topology_view = topology_span;
  }

private:
  // Private so set_view() stays the only way to enter borrowed mode: its
  // consistency traps are what downstream fi + fo[f] indexing rests on.
  ArenaSpan<uint8_t>
      face_counts_view; /**< Borrowed face-counts view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      faces_view; /**< Borrowed faces view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      face_offsets_view; /**< Borrowed face-offsets view, populated by MeshOps::transform. */
  ArenaSpan<uint16_t>
      topology_view; /**< Borrowed per-face topology view, populated by MeshOps::transform. */
};
