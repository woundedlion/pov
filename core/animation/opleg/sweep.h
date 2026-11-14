/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sweep.h
 * @brief Conway sweep preparation, sampling and endpoint provenance. */

#include "core/animation/opleg.h"

namespace Animation {

inline __attribute__((always_inline)) OpLeg::SweepPhase
OpLeg::sweep_phase(const Transients &tr, float progress) const {
  float k = progress;
  float settle_alpha = 0.0f;
  if (tr.settle_frames > 0) {
    float split = (tr.reverse ? tr.settle_frames : tr.sweep_frames) /
                  static_cast<float>(duration);
    if (!tr.reverse) {
      k = std::min(progress / split, 1.0f);
      settle_alpha = std::max(0.0f, (progress - split) / (1.0f - split));
    } else {
      settle_alpha = 1.0f - std::min(progress / split, 1.0f);
      k = std::max(0.0f, (progress - split) / (1.0f - split));
    }
  }
  return {k, settle_alpha};
}

inline __attribute__((always_inline)) OpLeg::SweepPhase
OpLeg::phase_at(const Transients &tr, float frame) const {
  // Clamp easing to the compiled topology-constant interval.
  const float progress =
      hs::clamp(easing_fn(frame / static_cast<float>(duration)), 0.0f, 1.0f);
  return sweep_phase(tr, progress);
}

inline HS_FLASH_MEMBER PolyMesh OpLeg::run_op(ConwayGraph::MorphOp op,
                                              const PolyMesh &seed,
                                              Arena &target, Arena &temp,
                                              float t, float twist) {
  switch (op) {
  case ConwayGraph::MorphOp::TRUNCATE:
    // Nudge t = 0.5 off the ambo short-circuit so a far-side leg keeps the
    // truncate topology.
    return MeshOps::truncate(seed, target, temp,
                             ConwayGraph::truncate_off_pinch(t));
  case ConwayGraph::MorphOp::EXPAND:
    return MeshOps::expand(seed, target, temp, t);
  case ConwayGraph::MorphOp::CHAMFER:
    return MeshOps::chamfer(seed, target, temp, t);
  default:
    return MeshOps::snub(seed, target, temp, t, twist);
  }
}

inline __attribute__((always_inline)) void
OpLeg::conway_at(const Transients &tr, PolyMesh &out, Arena &target,
                 Arena &temp, float t, float twist, float settle_alpha) {
  HS_CHECK(tr.seed_ref, "OpLeg: leg carries no seed mesh");
  out = run_op(tr.op, *tr.seed_ref, target, temp, t, twist);
  if (settle_alpha > 0.0f)
    settle_toward_relaxed(tr, out, settle_alpha);
}

inline HS_FLASH_MEMBER void OpLeg::swept_at(const Transients &tr, int frame,
                                            PolyMesh &out) const {
  const SweepPhase phase = phase_at(tr, static_cast<float>(frame));
  float tp = tr.t_start + (tr.t_end - tr.t_start) * phase.k;
  float tw = tr.twist_start + (tr.twist_end - tr.twist_start) * phase.k;
  if (tr.kind == LegKind::HANKIN_SWEEP)
    hankin_at(tr, out, scratch_arena_a, tp);
  else if (tr.kind == LegKind::RELAX_SLERP)
    relax_at(tr, out, scratch_arena_a, tp);
  else if (tr.kind == LegKind::MEDIAL_SLERP)
    medial_at(tr, out, scratch_arena_a, tp);
  else
    conway_at(tr, out, scratch_arena_a, scratch_arena_b, tp, tw,
              phase.settle_alpha);
}

inline HS_FLASH_MEMBER float OpLeg::sweep_parameter(float frame) const {
  const Transients &tr = *buf;
  return tr.t_start + (tr.t_end - tr.t_start) * phase_at(tr, frame).k;
}

inline HS_FLASH_MEMBER math::Vector OpLeg::swept_vertex_at(size_t i,
                                                           float k) const {
  const Transients &tr = *buf;
  if (tr.kind == LegKind::HANKIN_SWEEP) {
    const size_t statics = tr.hankin.static_vertices.size();
    if (i < statics)
      return tr.hankin.static_vertices[i];
    return k >= 1.0f ? hankin_vertex_at<true>(tr, i - statics, k)
                     : hankin_vertex_at<false>(tr, i - statics, k);
  }
  if (tr.kind == LegKind::RELAX_SLERP)
    return k >= 1.0f ? relax_vertex_at<true>(tr, i, k)
                     : relax_vertex_at<false>(tr, i, k);
  return k >= 1.0f ? medial_vertex_at<true>(tr, i, k)
                   : medial_vertex_at<false>(tr, i, k);
}

template <bool ARRIVAL>
inline __attribute__((always_inline)) math::Vector
OpLeg::hankin_vertex_at(const Transients &tr, size_t i, float k) {
  if constexpr (ARRIVAL)
    return tr.hk_final[i].decode().normalized();
  return math::slerp(
      tr.hankin.corner(tr.hankin.dynamic_instructions[i].v_corner).normalized(),
      tr.hk_final[i].decode(), k);
}

template <bool ARRIVAL>
inline __attribute__((always_inline)) math::Vector
OpLeg::relax_vertex_at(const Transients &tr, size_t i, float k) {
  if constexpr (ARRIVAL)
    return tr.relaxed[i];
  return math::slerp(tr.seed.vertices[i], tr.relaxed[i], k);
}

template <bool ARRIVAL>
inline __attribute__((always_inline)) math::Vector
OpLeg::medial_vertex_at(const Transients &tr, size_t i, float k) {
  if constexpr (ARRIVAL)
    return tr.medial_b[i].decode().normalized();
  return math::slerp(tr.medial_a[i].decode(), tr.medial_b[i].decode(), k);
}

inline HS_FLASH_MEMBER void
OpLeg::copy_topology(PolyMesh &out, Arena &arena,
                     const ArenaVector<uint8_t> &face_counts,
                     const ArenaVector<uint16_t> &faces) {
  out.face_counts.bind(arena, face_counts.size());
  out.face_counts.append_bulk(face_counts.data(), face_counts.size());
  out.faces.bind(arena, faces.size());
  out.faces.append_bulk(faces.data(), faces.size());
}

inline HS_FLASH_MEMBER void OpLeg::clone_geometry(const PolyMesh &src,
                                                  PolyMesh &dst, Arena &arena) {
  dst.vertices.bind(arena, src.vertices.size());
  dst.vertices.append_bulk(src.vertices.data(), src.vertices.size());
  copy_topology(dst, arena, src.face_counts, src.faces);
}

inline HS_FLASH_MEMBER void
OpLeg::hankin_at(const Transients &tr, PolyMesh &out, Arena &arena, float k) {
  const CompiledHankin &hk = tr.hankin;
  const size_t statics = hk.static_vertices.size();
  const size_t dyn = tr.hk_final.size();
  out.vertices.bind(arena, statics + dyn);
  out.vertices.append_bulk(hk.static_vertices.data(), statics);
  if (k >= 1.0f) {
    for (size_t i = 0; i < dyn; ++i)
      out.vertices.push_back(hankin_vertex_at<true>(tr, i, k));
  } else {
    for (size_t i = 0; i < dyn; ++i)
      out.vertices.push_back(hankin_vertex_at<false>(tr, i, k));
  }
  copy_topology(out, arena, hk.face_counts, hk.faces);
}

inline HS_FLASH_MEMBER void OpLeg::relax_at(const Transients &tr, PolyMesh &out,
                                            Arena &arena, float k) {
  const size_t n = tr.relaxed.size();
  out.vertices.bind(arena, n);
  if (k >= 1.0f)
    out.vertices.append_bulk(tr.relaxed.data(), n);
  else
    for (size_t i = 0; i < n; ++i)
      out.vertices.push_back(relax_vertex_at<false>(tr, i, k));
  copy_topology(out, arena, tr.seed.face_counts, tr.seed.faces);
}

inline HS_FLASH_MEMBER void
OpLeg::medial_at(const Transients &tr, PolyMesh &out, Arena &arena, float k) {
  const size_t n = tr.medial_a.size();
  out.vertices.bind(arena, n);
  if (k >= 1.0f) {
    for (size_t i = 0; i < n; ++i)
      out.vertices.push_back(medial_vertex_at<true>(tr, i, k));
  } else {
    for (size_t i = 0; i < n; ++i)
      out.vertices.push_back(medial_vertex_at<false>(tr, i, k));
  }
  copy_topology(out, arena, tr.seed.face_counts, tr.seed.faces);
}

inline __attribute__((always_inline)) void
OpLeg::settle_toward_relaxed(const Transients &tr, PolyMesh &swept,
                             float settle_alpha) {
  HS_CHECK(swept.vertices.size() == tr.relaxed.size(),
           "OpLeg relax: swept vertices %lu != relaxed vertices %lu",
           static_cast<unsigned long>(swept.vertices.size()),
           static_cast<unsigned long>(tr.relaxed.size()));
  // Alpha 1 copies the relaxed endpoint verbatim, so the settled bookend is
  // bitwise it.
  if (settle_alpha >= 1.0f)
    for (size_t i = 0; i < swept.vertices.size(); ++i)
      swept.vertices[i] = tr.relaxed[i];
  else
    for (size_t i = 0; i < swept.vertices.size(); ++i)
      swept.vertices[i] =
          math::slerp(swept.vertices[i], tr.relaxed[i], settle_alpha);
}

inline __attribute__((always_inline)) float
OpLeg::clamp_sweep_param(float t, bool truncate, bool far_side,
                         float trunc_floor) {
  t = std::max(t, truncate ? trunc_floor : ConwayGraph::T_EPS);
  if (truncate)
    t = std::min(t, far_side ? ConwayGraph::T_TRUNCATE_FAR_MAX
                             : ConwayGraph::T_TRUNCATE_PINCH_LOWER);
  return t;
}

inline HS_FLASH_MEMBER OpLeg::Transients &
OpLeg::init_transients(LegKind kind, int sweep_frames, Arena &arena,
                       const PaletteHandoff &handoff, BlendWeightFn blend_fn) {
  HS_CHECK(handoff.bank && handoff.prev_face_palette && handoff.prev_faces > 0,
           "OpLeg: incomplete palette handoff");
  Transients &tr = bind_transients(arena);
  tr.kind = kind;
  tr.sweep_frames = sweep_frames;
  tr.bank = handoff.bank;
  tr.blend_fn = blend_fn;
  return tr;
}

inline HS_FLASH_MEMBER void
OpLeg::bind_sweep_seed(Transients &tr, SweepSeed source, Arena &arena) {
  if (source.borrowed) {
    tr.seed_ref = &source.mesh;
  } else {
    clone_geometry(source.mesh, tr.seed, arena);
    tr.seed_ref = &tr.seed;
  }
  tr.seed_faces = source.mesh.face_counts.size();
}

template <typename Emit>
inline HS_FLASH_MEMBER void
OpLeg::first_seen_corners(const PolyMesh &dual, Arena &arena, Emit &&emit) {
  const size_t v = dual.vertices.size();
  bool *seen = arena.allocate_n<bool>(v);
  std::fill_n(seen, v, false);
  size_t out = 0, off = 0;
  for (size_t face = 0; face < dual.face_counts.size(); ++face) {
    for (int j = 0; j < dual.face_counts[face]; ++j) {
      const uint16_t vertex = dual.faces[off + j];
      HS_CHECK(vertex < v, "OpLeg: dual face index out of range");
      if (seen[vertex])
        continue;
      seen[vertex] = true;
      emit(vertex, static_cast<int>(face));
      ++out;
    }
    off += dual.face_counts[face];
  }
  HS_CHECK(out == v, "OpLeg: dual corners incomplete");
}

inline HS_FLASH_MEMBER const uint8_t *
OpLeg::dual_closing_palettes(const PolyMesh &dual,
                             const PaletteHandoff &handoff, Arena &scratch) {
  const size_t dual_faces = dual.face_counts.size();
  const size_t primal_faces = dual.vertices.size();
  HS_CHECK(!handoff.prev_face_centroid &&
               handoff.prev_faces == dual_faces + primal_faces,
           "OpLeg: invalid structural dual handoff");

  uint8_t *from = scratch.allocate_n<uint8_t>(handoff.prev_faces);
  for (size_t f = 0; f < dual_faces; ++f)
    from[f] = handoff.prev_face_palette[primal_faces + f];

  size_t out = dual_faces;
  first_seen_corners(dual, scratch, [&](uint16_t vertex, int) {
    from[out++] = handoff.prev_face_palette[vertex];
  });
  return from;
}

inline HS_FLASH_MEMBER void OpLeg::init_conway(const PaletteHandoff &handoff,
                                               const BookendClasses &bookend,
                                               Arena &arena,
                                               const ConwayLegFlags &flags) {
  const bool settle = flags.settle;
  const bool jitterbug = flags.jitterbug;
  const bool bridge_provenance = flags.bridge_provenance;
  Transients &tr = *buf;
  ScratchScope sa(scratch_arena_a);
  ScratchScope sb(scratch_arena_b);
  HS_CHECK(tr.seed_ref, "OpLeg: leg carries no seed mesh");
  const PolyMesh &seed = *tr.seed_ref;
  const uint8_t *forced_from = nullptr;
  if (handoff.correspondence == FaceCorrespondence::DUAL_CLOSING)
    forced_from = dual_closing_palettes(seed, handoff, scratch_arena_a);

  if (
#if HS_ENABLE_TEST_HOOKS
      !test_classify_partial_bookend &&
#endif
      forced_from && bookend.topology &&
      bookend.faces == seed.face_counts.size() && !settle) {
    tr.topo.bind(arena, seed.face_counts.size() + seed.vertices.size());
    for (size_t f = 0; f < seed.face_counts.size() + seed.vertices.size(); ++f)
      tr.topo.push_back(0);
    build_palette_mapping(tr, seed, handoff, bookend, arena, nullptr,
                          tr.seed_faces, forced_from);
    return;
  }

  // Closing-bridge face blocks transpose the handoff order; use the start
  // mesh's own centroids, built before arrival to avoid co-resident scratch.
  const math::Vector *start_centroid = nullptr;
  const bool closing_bridge =
      bridge_provenance && handoff.prev_face_centroid && tr.t_start > tr.t_end;
  if (closing_bridge)
    start_centroid = closing_bridge_centroids(tr, seed, handoff);

  PolyMesh arrival = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                            tr.t_end, tr.twist_end);
  PolyMesh *classified = &arrival;
  PolyMesh relaxed_mesh;
  if (settle)
    classified = hoist_settle_endpoint(tr, seed, arrival, arena, relaxed_mesh);
  tr.topo_is_bookend =
      hoist_arrival_topology(*classified, bookend, arena, tr.topo);

  PolyMesh start_mesh;
  if (handoff.prev_face_centroid) {
    if (bridge_provenance) {
      // Opening bridge leg: centroids are read only for its newborn faces,
      // where the arrival's own centroids stand in for the birth mesh.
      if (!start_centroid)
        start_centroid = face_centroids(*classified, scratch_arena_a);
    } else {
      start_centroid =
          opening_frame_centroids(tr, seed, relaxed_mesh, settle, start_mesh);
    }
  }

  build_palette_mapping(tr, *classified, handoff, bookend, arena,
                        start_centroid, survivor_faces(tr, seed, jitterbug),
                        forced_from);
}

inline __attribute__((always_inline)) const math::Vector *
OpLeg::closing_bridge_centroids(const Transients &tr, const PolyMesh &seed,
                                const PaletteHandoff &handoff) {
  math::Vector *cen =
      scratch_arena_a.allocate_n<math::Vector>(handoff.prev_faces);
  {
    ScratchScope ta(scratch_arena_a);
    ScratchScope tb(scratch_arena_b);
    PolyMesh start = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                            tr.t_start, tr.twist_start);
    HS_CHECK(start.face_counts.size() == handoff.prev_faces,
             "OpLeg: closing bridge start faces differ from the handoff");
    MeshOps::face_centroids_into(start, cen);
  }
  return cen;
}

inline __attribute__((always_inline)) PolyMesh *
OpLeg::hoist_settle_endpoint(Transients &tr, const PolyMesh &seed,
                             PolyMesh &arrival, Arena &arena,
                             PolyMesh &relaxed_mesh) {
  PolyMesh *classified = &arrival;
  if (!tr.reverse) {
    relaxed_mesh = MeshOps::relax(arrival, scratch_arena_b, scratch_arena_a,
                                  ConwayGraph::SETTLE_RELAX_ITERATIONS);
    classified = &relaxed_mesh;
  } else {
    PolyMesh start = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                            tr.t_start, tr.twist_start);
    relaxed_mesh = MeshOps::relax(start, scratch_arena_b, scratch_arena_a,
                                  ConwayGraph::SETTLE_RELAX_ITERATIONS);
  }
  tr.relaxed.bind(arena, relaxed_mesh.vertices.size());
  tr.relaxed.append_bulk(relaxed_mesh.vertices.data(),
                         relaxed_mesh.vertices.size());
  HS_CHECK(tr.relaxed.size() == arrival.vertices.size(),
           "OpLeg: relax changed the vertex count");
  return classified;
}

inline __attribute__((always_inline)) const math::Vector *
OpLeg::opening_frame_centroids(const Transients &tr, const PolyMesh &seed,
                               const PolyMesh &relaxed_mesh, bool settle,
                               PolyMesh &start_mesh) {
  const PolyMesh *start = &relaxed_mesh;
  if (!(settle && tr.reverse)) {
    start_mesh = run_op(tr.op, seed, scratch_arena_a, scratch_arena_b,
                        tr.t_start, tr.twist_start);
    start = &start_mesh;
  }
  HS_CHECK(start->face_counts.size() == tr.topo.size(),
           "OpLeg: start face count differs from arrival");
  return face_centroids(*start, scratch_arena_a);
}

} // namespace Animation
