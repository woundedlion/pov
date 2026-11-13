/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file build.h
 * @brief Conway morph leg construction. */

#include "core/animation/opleg.h"

namespace Animation {

inline HS_FLASH_MEMBER OpLeg::BridgePreparation
OpLeg::prepare_bridge(const PolyMesh &seed, Arena &arena) {
  ScratchScope a(scratch_arena_a);
  ScratchScope b(scratch_arena_b);
  PolyMesh arrival = MeshOps::truncate(seed, scratch_arena_a, scratch_arena_b,
                                       ConwayGraph::T_TRUNCATE_PINCH_LOWER);
  MeshOps::classify_faces_by_topology(arrival, scratch_arena_a, scratch_arena_b,
                                      arena);
  BridgePreparation prepared;
  prepared.topology = std::move(arrival.topology);
  prepared.hosts.fill(UINT16_MAX);
  const auto *previous = face_centroids(seed, scratch_arena_a);
  const auto *centroids = face_centroids(arrival, scratch_arena_a);
  for (size_t f = seed.face_counts.size(); f < prepared.topology.size(); ++f) {
    const int slot = MeshPaletteBank::slot_of(prepared.topology[f]);
    if (prepared.hosts[slot] == UINT16_MAX)
      prepared.hosts[slot] = static_cast<uint16_t>(
          nearest_point(centroids[f], previous, seed.face_counts.size()));
  }
  return prepared;
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(const PolyMesh &seed, const EdgeSweepSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend, BlendWeightFn blend_fn,
             EasingFn easing_fn)
    : OpLeg(spec.sweep_frames + spec.settle_frames, easing_fn, draw) {
  HS_CHECK(spec.edge, "OpLeg: edge sweep carries no graph edge");
  const ConwayGraph::EdgeSpec &edge = *spec.edge;
  HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
  HS_CHECK(spec.settle_frames >= 0 && edge.settle == (spec.settle_frames > 0),
           "OpLeg: settle frames disagree with the edge");
  Transients &tr = init_transients(LegKind::CONWAY_SWEEP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  bind_sweep_seed(tr, seed, arena);
  tr.op = edge.op;
  tr.reverse = spec.reverse;
  tr.settle_frames = spec.settle_frames;

  tr.t_start = ConwayGraph::clamp_edge_endpoint(
      edge, spec.reverse ? edge.t_to : edge.t_from);
  tr.t_end = ConwayGraph::clamp_edge_endpoint(edge, spec.reverse ? edge.t_from
                                                                 : edge.t_to);
  tr.twist_start = spec.reverse ? edge.twist_to : edge.twist_from;
  tr.twist_end = spec.reverse ? edge.twist_from : edge.twist_to;

  init_conway(handoff, bookend, arena,
              {.settle = edge.settle,
               .jitterbug = ConwayGraph::is_jitterbug_edge(edge)});
  seal_transients();
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(SweepSeed source, const ParamSweepSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend, BlendWeightFn blend_fn,
             EasingFn easing_fn)
    : OpLeg(spec.sweep_frames, easing_fn, draw) {
  HS_CHECK(spec.sweep_frames >= 1,
           "OpLeg: parameter sweep needs a positive sweep length");
  Transients &tr = init_transients(LegKind::CONWAY_SWEEP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  bind_sweep_seed(tr, source, arena);
  tr.op = spec.op;

  // Truncate births below T_EPS when the larger endpoint is below
  // T_EPS / TRUNCATE_BIRTH_FRAC.
  const bool truncate = spec.op == ConwayGraph::MorphOp::TRUNCATE;
  // A far-side truncate leg (an endpoint > 0.5) sweeps through the ambo
  // pinch on the constant-topology truncate branch; a near-side leg is
  // capped below ambo.
  const bool far_side = truncate && std::max(spec.t_start, spec.t_end) > 0.5f;
  HS_CHECK(!truncate || std::max(spec.t_start, spec.t_end) > 0.0f,
           "OpLeg: truncate sweep needs a positive endpoint");
  const float trunc_floor =
      ConwayGraph::truncate_birth_floor(std::max(spec.t_start, spec.t_end));
  tr.t_start = clamp_sweep_param(spec.t_start, truncate, far_side, trunc_floor);
  tr.t_end = clamp_sweep_param(spec.t_end, truncate, far_side, trunc_floor);
  tr.twist_start = spec.twist_start;
  tr.twist_end = spec.twist_end;

  if (spec.prepared) {
    HS_CHECK(spec.bridge_provenance && truncate && spec.t_start == 0.0f &&
                 spec.t_end == 0.5f && !bookend.topology &&
                 handoff.prev_face_centroid != nullptr &&
                 handoff.correspondence == FaceCorrespondence::GEOMETRIC &&
                 handoff.prev_faces == source.mesh.face_counts.size(),
             "OpLeg: incompatible prepared bridge");
    tr.topo.bind(arena, spec.prepared->topology.size());
    tr.topo.append_bulk(spec.prepared->topology.data(),
                        spec.prepared->topology.size());
    ScratchScope a(scratch_arena_a);
    auto *from = scratch_arena_a.allocate_n<uint8_t>(tr.topo.size());
    for (size_t f = 0; f < tr.topo.size(); ++f) {
      const size_t host =
          f < handoff.prev_faces
              ? f
              : spec.prepared->hosts[MeshPaletteBank::slot_of(tr.topo[f])];
      HS_CHECK(host < handoff.prev_faces,
               "OpLeg: prepared bridge host out of range");
      from[f] = handoff.prev_face_palette[host];
    }
    build_palette_mapping(tr, source.mesh, handoff, bookend, arena, nullptr,
                          tr.seed_faces, from);
  } else {
    init_conway(handoff, bookend, arena,
                {.bridge_provenance = spec.bridge_provenance});
  }
  seal_transients();
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(const PolyMesh &seed, const HankinSweepSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend, BlendWeightFn blend_fn,
             EasingFn easing_fn)
    : OpLeg(spec.sweep_frames, easing_fn, draw) {
  HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
  Transients &tr = init_transients(LegKind::HANKIN_SWEEP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  tr.seed_faces = seed.face_counts.size();

  // Hankin legs sweep the slerp fraction from each star point's collapsed
  // corner, not the contact angle.
  HS_CHECK(spec.theta_start <= spec.theta_end,
           "OpLeg: hankin leg sweeps back to a smaller contact angle");
  const float theta_hi = std::max(spec.theta_end, THETA_EPS);
  tr.t_start = std::max(std::max(spec.theta_start, 0.0f) / theta_hi, K_EPS);
  tr.t_end = 1.0f;

  {
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    MeshOps::compile_hankin(seed, tr.hankin, arena, scratch_arena_a,
                            MeshOps::BaseVertices::BORROW);

    PolyMesh arrival;
    MeshOps::update_hankin(tr.hankin, arrival, scratch_arena_a, theta_hi);

    // Hoist classes before the snorm16 pack unless the exact bookend already
    // supplies them; quantization can split classifier angle buckets.
    tr.topo_is_bookend =
        hoist_arrival_topology(arrival, bookend, arena, tr.topo);

    // Only the star points are stored; the midpoint prefix is already in the
    // compiled topology.
    const size_t statics = tr.hankin.static_vertices.size();
    HS_CHECK(arrival.vertices.size() >= statics,
             "OpLeg hankin: arrival has fewer vertices than statics");
    const size_t dyn = arrival.vertices.size() - statics;
    tr.hk_final.bind(arena, dyn);
    for (size_t i = 0; i < dyn; ++i)
      tr.hk_final.push_back(
          math::Snorm3::encode(arrival.vertices[statics + i]));
    tr.landing.hankin = &tr.hankin;
    tr.landing.star_point = tr.hk_final.data();
    tr.landing.star_points = tr.hk_final.size();

    const math::Vector *start_centroid = nullptr;
    PolyMesh start_mesh;
    if (handoff.prev_face_centroid) {
      hankin_at(tr, start_mesh, scratch_arena_b, tr.t_start);
      HS_CHECK(start_mesh.face_counts.size() == tr.topo.size(),
               "OpLeg: start face count differs from arrival");
      start_centroid = face_centroids(start_mesh, scratch_arena_b);
    }

    // Star faces are emitted first, in base-face order, so the seed's face
    // count is the emission-order prefix corresponding 1:1 to it.
    const size_t survivors = tr.seed_faces;
    build_palette_mapping(tr, arrival, handoff, bookend, arena, start_centroid,
                          survivors);
  }
  seal_transients();
}

inline HS_FLASH_MEMBER OpLeg::OpLeg(const PolyMesh &seed, const RelaxSpec &spec,
                                    Arena &arena, MorphDrawFn draw,
                                    const PaletteHandoff &handoff,
                                    const BookendClasses &bookend,
                                    BlendWeightFn blend_fn, EasingFn easing_fn)
    : OpLeg(spec.sweep_frames, easing_fn, draw) {
  HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
  HS_CHECK(spec.bake || spec.iterations >= 1,
           "OpLeg: relax leg needs a positive iteration count");
  Transients &tr = init_transients(LegKind::RELAX_SLERP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  // relax_at slerps out of the seed vertices every frame, so the seed stays
  // in the leg arena; only its per-face class ids are dead here.
  clone_geometry(seed, tr.seed, arena);
  tr.seed_faces = seed.face_counts.size();
  tr.t_start = 0.0f;
  tr.t_end = 1.0f;

  {
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    PolyMesh arrival =
        spec.bake ? MeshOps::relax_baked(tr.seed, scratch_arena_a, *spec.bake)
                  : MeshOps::relax(tr.seed, scratch_arena_a, scratch_arena_b,
                                   spec.iterations);
    HS_CHECK(arrival.vertices.size() == tr.seed.vertices.size(),
             "OpLeg: relax changed the vertex count");
    tr.relaxed.bind(arena, arrival.vertices.size());
    tr.relaxed.append_bulk(arrival.vertices.data(), arrival.vertices.size());

    tr.topo_is_bookend =
        hoist_arrival_topology(arrival, bookend, arena, tr.topo);

    // The opening frame is the seed verbatim, so its face centroids are the
    // provenance source.
    const math::Vector *start_centroid = nullptr;
    if (handoff.prev_face_centroid)
      start_centroid = face_centroids(tr.seed, scratch_arena_a);

    build_palette_mapping(tr, arrival, handoff, bookend, arena, start_centroid,
                          tr.seed_faces);
  }
  seal_transients();
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(const PolyMesh &seed, const MedialSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend, BlendWeightFn blend_fn,
             EasingFn easing_fn)
    : OpLeg(spec.sweep_frames, easing_fn, draw) {
  HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
  Transients &tr = init_transients(LegKind::MEDIAL_SLERP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  tr.t_start = 0.0f;
  tr.t_end = 1.0f;

  {
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    // Classify full-precision scratch points before packing persistent endpoints;
    // snorm16 quantization can split symmetry orbits.
    PolyMesh med;
    ArenaVector<math::Vector> med_b;
    MeshOps::medial(seed, med, med_b, scratch_arena_a, scratch_arena_b);
    tr.seed_faces = med.face_counts.size();
    HS_CHECK(med_b.size() == med.vertices.size(),
             "OpLeg: medial vertex sets differ in size");

    const size_t medial_verts = med.vertices.size();
    copy_topology(tr.seed, arena, med.face_counts, med.faces);
    tr.medial_a.bind(arena, medial_verts);
    tr.medial_b.bind(arena, medial_verts);
    for (size_t i = 0; i < medial_verts; ++i) {
      tr.medial_a.push_back(math::Snorm3::encode(med.vertices[i]));
      tr.medial_b.push_back(math::Snorm3::encode(med_b[i]));
    }

    // The opening frame is ambo(P) verbatim (the packed a_e), so med's s=0
    // vertices are the geometric provenance source; snapshot their centroids
    // before the in-place overwrite below.
    const math::Vector *start_centroid = nullptr;
    if (handoff.prev_face_centroid)
      start_centroid = face_centroids(med, scratch_arena_a);

    // Prepare the s=1 arrival (ambo(dual(P))) in place. The exact bookend can
    // supply its classes; otherwise the full-precision mesh is classified.
    for (size_t i = 0; i < medial_verts; ++i)
      med.vertices[i] = med_b[i];

    tr.topo_is_bookend = hoist_arrival_topology(med, bookend, arena, tr.topo);
    tr.landing.arrival_topology = &tr.seed;
    tr.landing.arrival_point = tr.medial_b.data();
    tr.landing.arrival_points = tr.medial_b.size();

    // Every medial face survives the whole bridge.
    build_palette_mapping(tr, med, handoff, bookend, arena, start_centroid,
                          tr.seed_faces);
  }
  seal_transients();
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(const PolyMesh &from_mesh, const ReconcileSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend, BlendWeightFn blend_fn,
             EasingFn easing_fn)
    : OpLeg(spec.sweep_frames, easing_fn, draw) {
  HS_CHECK(spec.sweep_frames >= 1, "OpLeg needs a positive sweep length");
  HS_CHECK(spec.to_positions, "OpLeg: reconcile leg carries no endpoints");
  HS_CHECK(spec.to_count == from_mesh.vertices.size(),
           "OpLeg: reconcile endpoint count must match seed vertices");
  Transients &tr = init_transients(LegKind::MEDIAL_SLERP, spec.sweep_frames,
                                   arena, handoff, blend_fn);

  tr.t_start = 0.0f;
  tr.t_end = 1.0f;

  {
    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);

    const size_t n = from_mesh.vertices.size();
    tr.seed_faces = from_mesh.face_counts.size();
    copy_topology(tr.seed, arena, from_mesh.face_counts, from_mesh.faces);
    // Both endpoint sets are unit-sphere directions, snorm16-packed.
    tr.medial_a.bind(arena, n);
    tr.medial_b.bind(arena, n);
    for (size_t i = 0; i < n; ++i) {
      tr.medial_a.push_back(math::Snorm3::encode(from_mesh.vertices[i]));
      tr.medial_b.push_back(math::Snorm3::encode(spec.to_positions[i]));
    }

    // Arrival: the identity connectivity carrying the authored positions.
    PolyMesh arrival;
    arrival.vertices.bind(scratch_arena_a, n);
    arrival.vertices.append_bulk(spec.to_positions, n);
    copy_topology(arrival, scratch_arena_a, from_mesh.face_counts,
                  from_mesh.faces);

    tr.topo_is_bookend =
        hoist_arrival_topology(arrival, bookend, arena, tr.topo);

    // The opening frame is the identity mesh verbatim (the packed a_e), so its
    // face centroids are the geometric provenance source.
    const math::Vector *start_centroid = nullptr;
    if (handoff.prev_face_centroid)
      start_centroid = face_centroids(from_mesh, scratch_arena_a);

    build_palette_mapping(tr, arrival, handoff, bookend, arena, start_centroid,
                          tr.seed_faces);
  }
  seal_transients();
}

inline HS_FLASH_MEMBER
OpLeg::OpLeg(const PolyMesh &seed, const GatedSwapSpec &spec, Arena &arena,
             MorphDrawFn draw, const PaletteHandoff &handoff,
             const BookendClasses &bookend)
    : OpLeg(2 * spec.gate_frames + 1, math::ease_in_out_sin, draw) {
  HS_CHECK(spec.gate_frames >= 1, "OpLeg needs a positive gate length");
  Transients &tr = init_transients(LegKind::GATED_SWAP, spec.gate_frames, arena,
                                   handoff, classic_blend);

  MeshOps::clone(seed, tr.seed, arena);
  tr.seed_faces = seed.face_counts.size();
  tr.swap_op = spec.op;

  init_gated(handoff, bookend, arena);
  seal_transients();
}

} // namespace Animation
