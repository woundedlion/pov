/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file ramps.h
 * @brief Morph palette mapping, frame shading and arena lifetime. */

#include "core/animation/opleg.h"
#include "core/animation/opleg/swaps.h"

namespace Animation {

inline HS_FLASH_MEMBER void
OpLeg::finish_frame(Canvas &canvas, const PolyMesh &swept, float w,
                    bool seed_side, math::Vector *drawn, int previous_frame,
                    int current_frame) {
  Transients &tr = *buf;
  const ArenaVector<uint16_t> &topo = seed_side ? tr.seed_topo : tr.topo;
  const ArenaVector<uint8_t> &face_ramp =
      seed_side ? tr.seed_face_ramp : tr.face_ramp;

  MeshState compiled;
  {
    HS_PROFILE(hk_conway_compile);
    MeshOps::compile(swept, compiled, scratch_arena_a, scratch_arena_b);
  }
  HS_CHECK(compiled.face_counts.size() == topo.size(),
           "OpLeg: sweep changed the compiled face count");

  const float previous_w = tr.drawn ? tr.drawn_weight : w;
  tr.drawn = true;
  tr.drawn_weight = w;
  Shading sh{blend_ramps(tr, w, seed_side), face_ramp.data(), face_ramp.size()};
  if (!seed_side) {
    sh.bank = tr.bank;
    sh.ramp_from = tr.ramp_from;
    sh.ramp_to = tr.ramp_to;
  }
  sh.weight = w;
  sh.previous_weight = previous_w;
  if (compiled.vertices.size() == swept.vertices.size()) {
    sh.previous_vertices = drawn;
    if (previous_frame >= 0) {
      sh.leg = this;
      sh.drawn_frame = static_cast<float>(previous_frame);
      sh.frame = static_cast<float>(current_frame);
    }
  }
  draw_fn(canvas, compiled, sh);
}

inline __attribute__((always_inline)) BakedPalette *
OpLeg::blend_ramps(const Transients &tr, float w, bool seed_side) {
  const int num_ramps = seed_side ? tr.seed_num_ramps : tr.num_ramps;
  BakedPalette *ramps = scratch_arena_b.make_n<BakedPalette>(num_ramps);
  for (int r = 0; r < num_ramps; ++r) {
    if (seed_side) {
      ramps[r] = tr.bank->entries[tr.seed_ramp_pal[r]];
      continue;
    }
    const BakedPalette &from = tr.bank->entries[tr.ramp_from[r]];
    const BakedPalette &to = tr.bank->entries[tr.ramp_to[r]];
    if (tr.ramp_from[r] == tr.ramp_to[r]) {
      ramps[r] = to;
    } else {
      ramps[r] = bake_palette_blend(scratch_arena_b, from, to, w);
    }
  }
  return ramps;
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

inline HS_FLASH_MEMBER const math::Vector *
OpLeg::face_centroids(const PolyMesh &m, Arena &arena) {
  math::Vector *out = arena.allocate_n<math::Vector>(m.face_counts.size());
  MeshOps::face_centroids_into(m, out);
  return out;
}

inline HS_FLASH_MEMBER size_t
OpLeg::nearest_prev_face(const math::Vector &c, const PaletteHandoff &handoff) {
  return nearest_point(c, handoff.prev_face_centroid, handoff.prev_faces);
}

inline __attribute__((always_inline)) size_t OpLeg::nearest_point(
    const math::Vector &c, const math::Vector *points, size_t n) {
  size_t best = 0;
  float best_d = 1e9f;
  for (size_t j = 0; j < n; ++j) {
    const math::Vector d = c - points[j];
    const float dsq = math::dot(d, d);
    if (dsq < best_d) {
      best_d = dsq;
      best = j;
    }
  }
  return best;
}

inline float OpLeg::blend_weight(float p) {
  constexpr float IN = 0.2f, OUT = 0.8f;
  if (p <= IN)
    return 0.0f;
  if (p >= OUT)
    return 1.0f;
  return math::cubic_kernel((p - IN) / (OUT - IN));
}

inline __attribute__((always_inline)) const uint16_t *
OpLeg::resolve_target_topology(Transients &tr, const PolyMesh &arrival,
                               const PaletteHandoff &handoff,
                               const BookendClasses &bookend, Arena &arena,
                               size_t survivors) {
  const size_t total = tr.topo.size();
  HS_CHECK(!bookend.topology || bookend.faces == total ||
               bookend.faces == survivors,
           "OpLeg: bookend face count matches neither mapping");
  // A hoisted bookend classification already covers every face.
  if (!bookend.topology || tr.topo_is_bookend) {
    tr.target_topo.clear();
    return tr.topo.data();
  }

  const size_t landed = bookend.faces;
  const bool structural_closing =
      handoff.correspondence == FaceCorrespondence::DUAL_CLOSING;
  int *corner_host = nullptr;
  if (structural_closing && landed < total)
    corner_host = dual_corner_hosts(tr, landed, total);
  const math::Vector *arrival_centroid =
      landed < total && !structural_closing
          ? face_centroids(arrival, scratch_arena_a)
          : nullptr;
  tr.target_topo.bind(arena, total);
  for (size_t f = 0; f < total; ++f) {
    if (f < landed) {
      tr.target_topo.push_back(bookend.topology[f]);
      continue;
    }
    if (structural_closing) {
      tr.target_topo.push_back(bookend.topology[corner_host[f - landed]]);
      continue;
    }
    const size_t host =
        nearest_point(arrival_centroid[f], arrival_centroid, landed);
    tr.target_topo.push_back(bookend.topology[host]);
  }
  return tr.target_topo.data();
}

inline __attribute__((always_inline)) int *
OpLeg::dual_corner_hosts(const Transients &tr, size_t landed, size_t total) {
  HS_CHECK(tr.seed_ref, "OpLeg: leg carries no seed mesh");
  const PolyMesh &dual = *tr.seed_ref;
  const size_t corners = total - landed;
  HS_CHECK(landed == dual.face_counts.size() && corners == dual.vertices.size(),
           "OpLeg: closing-leg blocks differ from the dual seed");
  int *corner_host = scratch_arena_a.allocate_n<int>(corners);
  size_t out = 0;
  first_seen_corners(dual, scratch_arena_a,
                     [&](uint16_t, int host) { corner_host[out++] = host; });
  return corner_host;
}

inline __attribute__((always_inline)) uint8_t
OpLeg::DepartedPalettes::for_face(size_t f, uint8_t to) const {
  uint8_t from = to;
  if (forced_from) {
    from = forced_from[f];
  } else if (handoff.correspondence == FaceCorrespondence::IDENTITY) {
    HS_CHECK(handoff.prev_faces == total,
             "OpLeg: identity handoff face count differs");
    from = handoff.prev_face_palette[f];
  } else if (full_correspondence) {
    from = claim_nearest(f);
  } else if (start_centroid) { // prev_faces == survivors
    if (f < handoff.prev_faces)
      from = handoff.prev_face_palette[f];
    else
      from = newborn(f);
  } else if (f < handoff.prev_faces) {
    from = handoff.prev_face_palette[f];
  }
  return from;
}

inline __attribute__((always_inline)) uint8_t
OpLeg::DepartedPalettes::claim_nearest(size_t f) const {
  const size_t j = nearest_prev_face(start_centroid[f], handoff);
  const math::Vector d = start_centroid[f] - handoff.prev_face_centroid[j];
  HS_CHECK(!prev_used[j] && math::dot(d, d) < PROVENANCE_TOL_SQ,
           "OpLeg: start face has no unique departed counterpart");
  prev_used[j] = true;
  return handoff.prev_face_palette[j];
}

inline __attribute__((always_inline)) uint8_t
OpLeg::DepartedPalettes::newborn(size_t f) const {
  const int slot = MeshPaletteBank::slot_of(target_topo[f]);
  if (newborn_from[slot] < 0)
    newborn_from[slot] =
        handoff
            .prev_face_palette[nearest_prev_face(start_centroid[f], handoff)];
  return static_cast<uint8_t>(newborn_from[slot]);
}

inline __attribute__((always_inline)) uint8_t
OpLeg::intern_palette_ramp(Transients &tr, uint8_t from, uint8_t to) {
  for (int r = 0; r < tr.num_ramps; ++r)
    if (tr.ramp_from[r] == from && tr.ramp_to[r] == to)
      return static_cast<uint8_t>(r);

  HS_CHECK(tr.num_ramps < MAX_BLEND_PAIRS,
           "OpLeg: distinct palette pairs exceed the blend budget");
  const int ramp = tr.num_ramps++;
  tr.ramp_from[ramp] = from;
  tr.ramp_to[ramp] = to;
  return static_cast<uint8_t>(ramp);
}

inline HS_FLASH_MEMBER void
OpLeg::build_palette_mapping(Transients &tr, const PolyMesh &arrival,
                             const PaletteHandoff &handoff,
                             const BookendClasses &bookend, Arena &arena,
                             const math::Vector *start_centroid,
                             size_t survivors, const uint8_t *forced_from) {
  const size_t total = tr.topo.size();
  const size_t primary = tr.seed_faces;
  const bool preserve_palette = handoff.preserve_without_births &&
                                tr.kind != LegKind::GATED_SWAP &&
                                handoff.prev_faces == total;
  tr.landing.faces = total;
  tr.landing.primary_faces = primary;
  tr.landing.preserve_palette = preserve_palette;

  // Surviving faces keep bookend classes; collapsed slivers take their
  // arrival host's class.
  const uint16_t *target_topo =
      resolve_target_topology(tr, arrival, handoff, bookend, arena, survivors);
  tr.landing.topology = target_topo;

  for (int i = 0; i < PALETTES; ++i)
    tr.landing.to_palette[i] = static_cast<uint8_t>(i);
  if (!preserve_palette)
    hs::shuffle(tr.landing.to_palette.begin(), tr.landing.to_palette.end());

  // The whole face list corresponds, or an emission prefix (survivors or
  // seed primaries) does and the remaining faces are births.
  HS_CHECK(handoff.prev_faces == total || handoff.prev_faces == survivors ||
               handoff.prev_faces == primary,
           "OpLeg: handoff face count matches no mapping");

  // Full-correspondence geometric mapping: nearest departed centroid,
  // pinned as a bijection within tolerance.
  const bool full_correspondence =
      start_centroid && handoff.prev_faces == total;
  bool *prev_used = nullptr;
  if (full_correspondence) {
    prev_used = scratch_arena_a.allocate_n<bool>(handoff.prev_faces);
    std::fill_n(prev_used, handoff.prev_faces, false);
  }
  // Newborn palette slots share the from-palette of their first face.
  int newborn_from[PALETTES];
  for (int i = 0; i < PALETTES; ++i)
    newborn_from[i] = -1;
  const DepartedPalettes departed{handoff,     forced_from, start_centroid,
                                  target_topo, total,       full_correspondence,
                                  prev_used,   newborn_from};

  uint8_t *from_palette = arena.allocate_n<uint8_t>(total);
  tr.landing.from_palette = from_palette;

  tr.face_ramp.bind(arena, total);
  for (size_t f = 0; f < total; ++f) {
    uint8_t to =
        tr.landing.to_palette[MeshPaletteBank::slot_of(target_topo[f])];
    const uint8_t from = departed.for_face(f, to);
    if (preserve_palette)
      to = from;
    HS_CHECK(from < PALETTES && to < PALETTES,
             "OpLeg: crossfade ramp endpoint outside the palette bank");
    from_palette[f] = from;

    tr.face_ramp.push_back(intern_palette_ramp(tr, from, to));
  }
  tr.landing.blend_pairs = tr.num_ramps;

  HS_CHECK(ramps_fit_scratch(tr),
           "OpLeg: blended palette LUTs exceed scratch_arena_b");
}

inline __attribute__((always_inline)) bool
OpLeg::ramps_fit_scratch(const Transients &tr) {
  int blended = 0;
  for (int r = 0; r < tr.num_ramps; ++r)
    if (tr.ramp_from[r] != tr.ramp_to[r])
      ++blended;
  return static_cast<size_t>(tr.num_ramps) * sizeof(BakedPalette) +
             static_cast<size_t>(blended) *
                 BakedPalette::required_arena_bytes() <=
         scratch_arena_b.get_capacity();
}

inline OpLeg::Transients &OpLeg::bind_transients(Arena &arena) {
  buf = arena.make<Transients>();
  leg_arena = &arena;
  live_end = arena.get_offset();
#ifndef NDEBUG
  stamp.record(arena);
  live_bytes = sizeof(Transients);
#endif
  return *buf;
}

inline void OpLeg::seal_transients() {
#ifndef NDEBUG
  live_bytes += leg_arena->get_offset() - live_end;
#endif
  live_end = leg_arena->get_offset();
}

inline void OpLeg::check_alive() const {
  HS_CHECK(leg_arena->get_offset() >= live_end,
           "OpLeg: leg arena rewound under a live leg");
#ifndef NDEBUG
  HS_CHECK(stamp.block_alive(buf, live_bytes),
           "OpLeg: leg arena storage reclaimed under a live leg");
#endif
}

} // namespace Animation
