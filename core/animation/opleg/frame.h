/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file frame.h
 * @brief Conway morph frame drawing and sweep sampling. */

#include "core/animation/opleg.h"
#include "core/animation/opleg/build.h"
#include <cstring>

namespace Animation {

inline HS_FLASH_MEMBER void OpLeg::step(Canvas &canvas) {
  AnimationBase::step(canvas);
  draw_frame(canvas);
}

inline HS_FLASH_MEMBER void OpLeg::step_paused(Canvas &canvas) {
  draw_frame(canvas);
}

inline HS_FLASH_MEMBER void OpLeg::draw_frame(Canvas &canvas) {
  check_alive();
  Transients &tr = *buf;
  const int frame =
      static_cast<int>(std::min<uint32_t>(t, static_cast<uint32_t>(duration)));

  if (tr.kind == LegKind::GATED_SWAP) {
    step_gated(canvas, frame);
    return;
  }

  ScratchScope sa(scratch_arena_a);
  ScratchScope sb(scratch_arena_b);

  PolyMesh swept;
  const size_t a_before = scratch_arena_a.get_offset();
  const size_t b_before = scratch_arena_b.get_offset();
  {
    HS_PROFILE(hk_conway_sweep);
    swept_at(tr, frame, swept);
  }
  // The swept vertices of the frame drawn last, for the draw to move from;
  // a leg that sweeps vertex by vertex hands its path instead.
  const bool moved = tr.drawn && tr.drawn_frame != frame;
  const bool per_vertex = tr.kind != LegKind::CONWAY_SWEEP;
  math::Vector *drawn = nullptr;
  if (moved && !per_vertex && swept.vertices.size() == tr.drawn_vertices &&
      previous_fits(swept.vertices.size(),
                    scratch_arena_a.get_offset() - a_before,
                    scratch_arena_b.get_offset() - b_before)) {
    HS_PROFILE(hk_conway_previous);
    drawn = scratch_arena_a.allocate_n<math::Vector>(tr.drawn_vertices);
    ScratchScope prior_a(scratch_arena_a);
    ScratchScope prior_b(scratch_arena_b);
    PolyMesh prior;
    swept_at(tr, tr.drawn_frame, prior);
    if (prior.vertices.size() == tr.drawn_vertices)
      std::memcpy(drawn, prior.vertices.data(),
                  tr.drawn_vertices * sizeof(math::Vector));
    else
      drawn = nullptr;
  }
  const int previous_frame = tr.drawn_frame;
  tr.drawn_frame = frame;
  tr.drawn_vertices = swept.vertices.size();

  finish_frame(canvas, swept, tr.blend_fn(frame, duration), false, drawn,
               moved && per_vertex ? previous_frame : -1, frame);
}

inline HS_FLASH_MEMBER bool
OpLeg::previous_fits(size_t vertices, size_t sweep_a, size_t sweep_b) {
  const size_t need_a =
      vertices * sizeof(math::Vector) + 2 * sweep_a + PREVIOUS_SLACK;
  const size_t need_b = sweep_b + PREVIOUS_SLACK;
  return scratch_arena_a.get_capacity() - scratch_arena_a.get_offset() >=
             need_a &&
         scratch_arena_b.get_capacity() - scratch_arena_b.get_offset() >=
             need_b;
}

inline HS_FLASH_MEMBER float OpLeg::sweep_parameter(float frame) const {
  const Transients &tr = *buf;
  const float progress =
      hs::clamp(easing_fn(frame / static_cast<float>(duration)), 0.0f, 1.0f);
  return tr.t_start + (tr.t_end - tr.t_start) * sweep_phase(tr, progress).k;
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

inline HS_FLASH_MEMBER void OpLeg::arrival_mesh(const Landing &landing,
                                                PolyMesh &out, Arena &arena) {
  HS_CHECK(landing.hankin, "OpLeg: leg carries no baked arrival");
  out = PolyMesh();
  const CompiledHankin &hk = *landing.hankin;
  const size_t statics = hk.static_vertices.size();
  out.vertices.bind(arena, statics + landing.star_points);
  out.vertices.append_bulk(hk.static_vertices.data(), statics);
  for (size_t i = 0; i < landing.star_points; ++i)
    out.vertices.push_back(landing.star_point[i].decode().normalized());
  copy_topology(out, arena, hk.face_counts, hk.faces);
}

inline float OpLeg::trailing_blend(int frame, int duration, int window) {
  const int win = std::max(1, std::min(window, duration - 1));
  const int from_end = duration - frame;
  if (from_end >= win)
    return 0.0f;
  return math::cubic_kernel(static_cast<float>(win - from_end) /
                            static_cast<float>(win));
}

} // namespace Animation
