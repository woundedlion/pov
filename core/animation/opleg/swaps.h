/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file swaps.h
 * @brief Gated morph swaps and topology classification. */

#include "core/animation/opleg.h"
#include "core/animation/opleg/sweep.h"

namespace Animation {

inline HS_FLASH_MEMBER void OpLeg::init_gated(const PaletteHandoff &handoff,
                                              const BookendClasses &bookend,
                                              Arena &arena) {
  Transients &tr = *buf;
  ScratchScope sa(scratch_arena_a);
  ScratchScope sb(scratch_arena_b);

  MeshOps::classify_faces_by_topology(tr.seed, scratch_arena_a, scratch_arena_b,
                                      arena);
  tr.seed_topo = std::move(tr.seed.topology);
  HS_CHECK(tr.seed_topo.size() == tr.seed.face_counts.size(),
           "OpLeg: gated seed classification is not one class per face");
  HS_CHECK(handoff.prev_faces == tr.seed.face_counts.size(),
           "OpLeg: gated swap departs from a mesh that is not its seed");

  PolyMesh arrival = swap_at(tr, scratch_arena_a, scratch_arena_b);
  MeshOps::classify_faces_by_topology(arrival, scratch_arena_a, scratch_arena_b,
                                      arena);
  tr.topo = std::move(arrival.topology);
  HS_CHECK(tr.topo.size() == arrival.face_counts.size(),
           "OpLeg: gated arrival classification is not one class per face");

  // Provenance across the swap: the children inherit their parent face's
  // (kis) or nearest orbit face's (dual) palette — colour locality.
  uint8_t *from = scratch_arena_a.allocate_n<uint8_t>(tr.topo.size());
  if (tr.swap_op == SwapOp::KIS)
    kis_provenance(tr, arrival, handoff, from);
  else
    dual_provenance(tr, arrival, handoff, from);

  build_palette_mapping(tr, arrival, handoff, bookend, arena,
                        /*start_centroid=*/nullptr, tr.seed.face_counts.size(),
                        from);

  // Seed-side table: the closing half draws at blend weight 0, so only each
  // pair's from palette is ever sampled and an identity pair suffices.
  const size_t seed_faces = tr.seed.face_counts.size();
  tr.seed_face_ramp.bind(arena, seed_faces);
  for (size_t f = 0; f < seed_faces; ++f) {
    const uint8_t pal = handoff.prev_face_palette[f];
    HS_CHECK(pal < PALETTES, "OpLeg palette: index out of range");
    tr.seed_face_ramp.push_back(intern_seed_ramp(tr, pal));
  }
}

inline __attribute__((always_inline)) uint8_t
OpLeg::intern_seed_ramp(Transients &tr, uint8_t pal) {
  int ramp = -1;
  for (int r = 0; r < tr.seed_num_ramps; ++r) {
    if (tr.seed_ramp_pal[r] == pal) {
      ramp = r;
      break;
    }
  }
  if (ramp < 0) {
    ramp = tr.seed_num_ramps++;
    tr.seed_ramp_pal[ramp] = pal;
  }
  return static_cast<uint8_t>(ramp);
}

inline HS_FLASH_MEMBER void OpLeg::step_gated(Canvas &canvas, int frame) {
  Transients &tr = *buf;
  const bool seed_side = frame <= tr.sweep_frames;

  ScratchScope sa(scratch_arena_a);
  ScratchScope sb(scratch_arena_b);

  PolyMesh swapped;
  {
    HS_PROFILE(hk_gate_swap);
    if (!seed_side)
      swapped = swap_at(tr, scratch_arena_a, scratch_arena_b);
  }
  const PolyMesh &mesh = seed_side ? tr.seed : swapped;

  // Colour holds at the departed palettes until the final frames.
  finish_frame(canvas, mesh,
               trailing_blend(frame, duration,
                              std::min(TRAILING_BLEND_FRAMES, tr.sweep_frames)),
               seed_side);
}

inline HS_FLASH_MEMBER PolyMesh OpLeg::swap_at(const Transients &tr,
                                               Arena &target, Arena &temp) {
  return tr.swap_op == SwapOp::KIS ? MeshOps::kis(tr.seed, target, temp)
                                   : MeshOps::dual(tr.seed, target, temp);
}

inline HS_FLASH_MEMBER void OpLeg::kis_provenance(const Transients &tr,
                                                  const PolyMesh &arrival,
                                                  const PaletteHandoff &handoff,
                                                  uint8_t *from) {
  size_t f = 0;
  for (size_t p = 0; p < tr.seed.face_counts.size(); ++p) {
    for (int k = 0; k < tr.seed.face_counts[p]; ++k) {
      HS_CHECK(f < arrival.face_counts.size(), "OpLeg: kis child overflow");
      from[f++] = handoff.prev_face_palette[p];
    }
  }
  HS_CHECK(f == arrival.face_counts.size(),
           "OpLeg: kis child count differs from the seed's side count");
}

inline HS_FLASH_MEMBER void
OpLeg::dual_provenance(const Transients &tr, const PolyMesh &arrival,
                       const PaletteHandoff &handoff, uint8_t *from) {
  size_t off = 0;
  for (size_t f = 0; f < arrival.face_counts.size(); ++f) {
    const int n = arrival.face_counts[f];
    const math::Vector c = math::normalized_or(
        MeshOps::face_vertex_sum(arrival.vertices.data(), arrival.faces.data(),
                                 off, n),
        arrival.vertices[arrival.faces[off]]);

    // The orbit's own source vertex: the dual face's centroid lies in that
    // vertex's cell, so the nearest seed vertex is it.
    const math::Vector v = tr.seed.vertices[nearest_point(
        c, tr.seed.vertices.data(), tr.seed.vertices.size())];
    const uint16_t best_face = nearest_orbit_vertex(arrival, off, n, v);
    HS_CHECK(best_face < handoff.prev_faces,
             "OpLeg: dual face index outside the departed face list");
    from[f] = handoff.prev_face_palette[best_face];
    off += n;
  }
}

inline __attribute__((always_inline)) uint16_t OpLeg::nearest_orbit_vertex(
    const PolyMesh &mesh, size_t off, int n, const math::Vector &v) {
  uint16_t best = mesh.faces[off];
  float best_sq = 1e9f;
  for (int k = 0; k < n; ++k) {
    const uint16_t j = mesh.faces[off + k];
    const math::Vector d = mesh.vertices[j] - v;
    const float dsq = math::dot(d, d);
    if (dsq < best_sq) {
      best_sq = dsq;
      best = j;
    }
  }
  return best;
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

inline HS_FLASH_MEMBER bool
OpLeg::hoist_arrival_topology(PolyMesh &arrival, const BookendClasses &bookend,
                              Arena &arena, ArenaVector<uint16_t> &topology) {
  const size_t faces = arrival.face_counts.size();
  bool reuse = bookend.topology != nullptr;
#if HS_ENABLE_TEST_HOOKS
  reuse = reuse && (!test_classify_partial_bookend || bookend.faces == faces);
#endif
  if (reuse) {
    topology.bind(arena, faces);
    if (bookend.faces == faces) {
      topology.append_bulk(bookend.topology, faces);
      return true;
    }
    // resolve_target_topology replaces every partial-bookend face slot.
    for (size_t f = 0; f < faces; ++f)
      topology.push_back(0);
    return false;
  }
#if HS_ENABLE_TEST_HOOKS
  test_classified_faces += faces;
  test_classified_corners += arrival.faces.size();
#endif
  MeshOps::classify_faces_by_topology(arrival, scratch_arena_a, scratch_arena_b,
                                      arena);
  topology = std::move(arrival.topology);
  HS_CHECK(topology.size() == faces,
           "OpLeg bookend: topology size %lu != face count %lu",
           static_cast<unsigned long>(topology.size()),
           static_cast<unsigned long>(faces));
  return false;
}

} // namespace Animation
