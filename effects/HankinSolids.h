/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/engine/engine.h"

// Unit-test accessor reaching the private graph-walk state (current node, held
// seed identity) so the §7.8 soak can pin node coverage and the per-state
// post-compaction arena footprint.
namespace hs_test {
namespace conway_soak_tests {
struct HankinWalkProbe;
} // namespace conway_soak_tests
} // namespace hs_test

/**
 * @brief Renders Hankin interlace patterns over Platonic/Archimedean solids.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Sweeps the interlace angle continuously, then transitions by
 * walking the Conway edge graph: each leg sweeps the destination solid's own
 * operator parameter, so faces visibly truncate, expand, and twist into the
 * next solid (docs/conway_morph_spec.md). Exactly one mesh is on screen at
 * all times; faces are colored by topology class via shuffled mesh palettes
 * with per-leg crossfades.
 */
template <int W, int H> class HankinSolids : public Effect {
public:
  /**
   * @brief Constructs the effect with a W x H canvas and empty filter pipeline.
   */
  HS_COLD_MEMBER HankinSolids()
      : Effect(W, H,
               {.strobe = true,
                .full_frame = decltype(filters)::any_crosses_segments}),
        filters() {}

  /**
   * @brief Sizes arenas, registers params, seeds the graph walk, and starts
   * the interlace sweep/morph cycle.
   */
  HS_COLD_MEMBER void init() override {
    // scratch_a (24 KB) covers its largest non-overlapping peak: the render
    // path (draw_mesh transforms the mesh into scratch_a, then Scan::Mesh::draw
    // stacks an SDF::FaceScratchBuffer on top), which for the heaviest hankin
    // mesh exceeds the generation/classify peak; morph frames stack the swept
    // mesh + compile output under the same render path. scratch_b (32 KB)
    // covers the largest of generation intermediates, the morph-cycle Persist
    // set, and a morph frame's blended palette LUTs. The per-solid H=144
    // scratch high-water is gated in CI by the arena-budget test.
    configure_arenas(GLOBAL_ARENA_SIZE - 24 * 1024 - 32 * 1024, 24 * 1024,
                     32 * 1024);
    register_param("Intensity", &params.intensity, 0.0f, 5.0f);
    register_animated_param("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);
    register_param("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    palette_bank_.bake_all(persistent_arena);

    // The walk starts at the tetrahedron with the registry seed; every family
    // seed is later derived from it (bridge ADOPTs), which is what lets a
    // reverse bridge regenerate the registry tetrahedron frame-exactly.
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed_base_ =
          Solids::finalize_solid(Solids::Platonic::tetrahedron(a, b), target);
    });
    node_ = ConwayGraph::TETRAHEDRON;
    seed_identity_ = ConwayGraph::TETRAHEDRON;
    ConwayGraph::record_visit(node_visits_, node_);

    MeshPaletteBank::shuffle_indices(palette_idx_);
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &) {
      compiled_hankin = CompiledHankin();
      MeshOps::compile_hankin(seed_base_, compiled_hankin, target, a);
      mesh_.clear();
      MeshOps::update_hankin(compiled_hankin, mesh_, target,
                             params.hankin_angle);
      record_node_faces(seed_base_);
    });
    classify_mesh_topology(mesh_);
    record_node_palettes();
    // First cycle has no predecessor to crossfade from: every slot opens on
    // its shuffled palette directly (strap_blend_mask_ stays 0).
    strap_from_ = palette_idx_;

    start_hankin_cycle();
  }

  /**
   * @brief Advances the timeline by one frame and renders into the canvas.
   */
  void draw_frame() override {
    // IIFE isolates the buffer_free() spin-wait in the Canvas ctor.
    Canvas canvas = [this]() -> Canvas {
      HS_PROFILE(hk_buffer_wait);
      return Canvas(*this);
    }();
    {
      HS_PROFILE(hk_timeline_step);
      timeline.step(canvas);
    }
  }

private:
  friend struct ::hs_test::conway_soak_tests::HankinWalkProbe;

  static constexpr int NUM_PALETTES = MeshPaletteBank::N;
  /** Largest node base-mesh face count (snubDodecahedron, F = 92). */
  static constexpr size_t MAX_NODE_FACES = 92;
  /** Largest hankin mesh face count (star faces plus the rosette/strap faces
   * born inside them); guarded by an HS_CHECK in resolve_host_faces. */
  static constexpr size_t MAX_HANKIN_FACES = 256;
  /** Frames over which the terminal rosette sliver fades out. Below about a
   * pixel wide the rosette is anti-aliasing only, so fading it uncovers
   * nothing — but it removes the one-frame snap as the last sliver and the
   * star-face boundary it holds off both disappear at once. */
  static constexpr int STRAP_TERMINAL_FRAMES = 3;
  /** Strap-crossfade window: frames of the 64-frame hankin sweep over which a
   * reborn strap slot glides from its previous color to the fresh target. The
   * sweep opens quadratically — straps reach visibility around frame 3-7 and
   * about half their peak area by frame 20 — so the turnover completes while
   * the straps are still small, and well before the closing bookend needs the
   * landed assignment displayed verbatim. */
  static constexpr int STRAP_BLEND_FRAMES = 20;

  /** Identity class-slot indexing for the per-slot resolved-LUT draw path. */
  static constexpr std::array<int, NUM_PALETTES> SLOT_IDENTITY = [] {
    std::array<int, NUM_PALETTES> a{};
    for (int i = 0; i < NUM_PALETTES; ++i)
      a[i] = i;
    return a;
  }();

  /** Indexable palette-bank view over per-slot resolved LUT pointers. */
  struct SlotLutView {
    const BakedPalette *const *slots; /**< One resolved LUT per class slot. */
    /**
     * @brief Returns the resolved LUT for a class slot.
     * @param i Class slot in [0, NUM_PALETTES).
     */
    const BakedPalette &operator[](int i) const { return *slots[i]; }
  };

  MeshPaletteBank palette_bank_;

  /**
   * @brief Classifies mesh faces into topology groups.
   * @param mesh Mesh whose faces are classified; the result is stored in
   * persistent_arena.
   * @details Saves/restores the scratch high-water marks rather than hard-
   * resetting to base, so a caller's prior allocations in the shared scratch
   * arenas survive.
   */
  HS_COLD_MEMBER void classify_mesh_topology(MeshState &mesh) {
    // ScratchScope frees only this call's own allocations, preserving prior
    // caller allocations in these shared arenas that a bare reset() would drop.
    ScratchScope a_guard(scratch_arena_a);
    ScratchScope b_guard(scratch_arena_b);
    MeshOps::classify_faces_by_topology(mesh, scratch_arena_a, scratch_arena_b,
                                        persistent_arena);
  }

  /**
   * @brief Records the node base mesh's face count, side counts, and unit
   * centroids for the next leg's palette handoff.
   * @param base The arrived node's base mesh, in emission order.
   */
  HS_COLD_MEMBER void record_node_faces(const PolyMesh &base) {
    node_faces_ = base.face_counts.size();
    HS_CHECK(node_faces_ <= MAX_NODE_FACES);
    size_t off = 0;
    for (size_t f = 0; f < node_faces_; ++f) {
      node_face_sides_[f] = base.face_counts[f];
      Vector c(0.0f, 0.0f, 0.0f);
      const int n = base.face_counts[f];
      for (int k = 0; k < n; ++k)
        c = c + base.vertices[base.faces[off + k]];
      node_face_centroid_[f] = c.normalized();
      off += n;
    }
  }

  /**
   * @brief Records the displayed palette of every node base face (via the
   * hankin star-face identity mapping) for the next leg's handoff.
   */
  HS_COLD_MEMBER void record_node_palettes() {
    HS_CHECK(node_faces_ <= mesh_.topology.size());
    for (size_t f = 0; f < node_faces_; ++f)
      node_face_palette_[f] = static_cast<uint8_t>(
          palette_idx_[wrap(mesh_.topology[f], NUM_PALETTES)]);
  }

  /**
   * @brief Strap-crossfade weight for a hankin-cycle frame: exactly 0 at the
   * opening bookend draw, smoothstep to exactly 1 by STRAP_BLEND_FRAMES.
   * @param cycle_frame Sprite draws since the cycle's opening bookend.
   */
  static float strap_blend_weight(int cycle_frame) {
    if (cycle_frame <= 0)
      return 0.0f;
    if (cycle_frame >= STRAP_BLEND_FRAMES)
      return 1.0f;
    const float u = static_cast<float>(cycle_frame) /
                    static_cast<float>(STRAP_BLEND_FRAMES);
    return u * u * (3.0f - 2.0f * u);
  }

  /**
   * @brief Resolves the per-slot palette LUTs one hankin-cycle frame shades
   * with, split by face role.
   * @param cycle_frame Sprite draws since the cycle's opening bookend.
   * @param blended Per-slot storage for this frame's blend results.
   * @param star_by_slot Receives the star-face LUT per class slot: always the
   * assignment's bank entry, bitwise (the star-face bookend exactness
   * contract).
   * @param strap_by_slot Receives the strap-face LUT per class slot: armed
   * slots run a (from, to) pre-blend over the opening window so reborn straps
   * open in their previous displayed color and glide to the fresh target
   * while still small; unarmed slots alias the star entry. The fragment path
   * stays a single LUT lookup either way.
   * @param scratch Arena receiving mid-blend LUT bakes (frame lifetime).
   */
  void resolve_hankin_slot_luts(
      int cycle_frame, BakedPalette (&blended)[NUM_PALETTES],
      const BakedPalette *(&star_by_slot)[NUM_PALETTES],
      const BakedPalette *(&strap_by_slot)[NUM_PALETTES], Arena &scratch) {
    const float w = strap_blend_weight(cycle_frame);
    for (int s = 0; s < NUM_PALETTES; ++s) {
      star_by_slot[s] = &palette_bank_.bank.entries[palette_idx_[s]];
      if (strap_blend_mask_ & (1u << s)) {
        bake_palette_blend(blended[s], scratch,
                           palette_bank_.bank.entries[strap_from_[s]],
                           palette_bank_.bank.entries[palette_idx_[s]], w);
        strap_by_slot[s] = &blended[s];
      } else {
        strap_by_slot[s] = star_by_slot[s];
      }
    }
  }

  /**
   * @brief Arms the opening-window strap crossfade for the new hankin cycle
   * (the spec-2.6 turnover applied to the rosette rebirths, which the
   * zero-area-birth assumption alone does not hide).
   * @param prev_idx Slot -> palette assignment the previous cycle displayed.
   * @param prev_used Slots the previous cycle put on screen at all.
   * @details A strap-bearing slot opens at its previous cycle's displayed
   * color — or, with no on-screen predecessor, at the palette of the nearest
   * star face by unit centroid (births inherit the underlying face's colors,
   * matching the leg-swap newborn rule) — and glides to its fresh target.
   * Slots a star class shares (the mod-NUM_PALETTES wrap can alias a rosette
   * class onto a star class's slot) arm like any other: the blend feeds only
   * the strap-face LUT, so the star faces stay on the exact bank entry.
   */
  HS_COLD_MEMBER void
  prepare_strap_crossfade(const std::array<int, NUM_PALETTES> &prev_idx,
                          const bool (&prev_used)[NUM_PALETTES]) {
    strap_blend_mask_ = 0;
    for (int s = 0; s < NUM_PALETTES; ++s)
      strap_from_[s] = palette_idx_[s];
    for (size_t f = node_faces_; f < mesh_.topology.size(); ++f) {
      const int slot = wrap(mesh_.topology[f], NUM_PALETTES);
      if (strap_blend_mask_ & (1u << slot))
        continue;
      int from;
      if (prev_used[slot]) {
        from = prev_idx[slot];
      } else {
        Vector c(0.0f, 0.0f, 0.0f);
        const size_t off = mesh_.face_offsets[f];
        const int n = mesh_.face_counts[f];
        for (int k = 0; k < n; ++k)
          c = c + mesh_.vertices[mesh_.faces[off + k]];
        c = c.normalized();
        size_t best = 0;
        float best_d = 1e9f;
        for (size_t j = 0; j < node_faces_; ++j) {
          const Vector d = c - node_face_centroid_[j];
          const float dsq = dot(d, d);
          if (dsq < best_d) {
            best_d = dsq;
            best = j;
          }
        }
        from = node_face_palette_[best];
      }
      if (from != palette_idx_[slot]) {
        strap_from_[slot] = from;
        strap_blend_mask_ |= static_cast<uint8_t>(1u << slot);
      }
    }
  }

  /**
   * @brief Builds the clean node mesh at one end of an edge from the held
   * seed — the registry chain, decomposed.
   * @param e Edge whose endpoint mesh is built.
   * @param to_end True for the t_to end, false for the t_from end.
   * @param a Output arena for even pipeline stages.
   * @param b Scratch arena for odd pipeline stages.
   * @return The endpoint mesh: t = 0 yields the leg seed itself, t = 0.5 the
   * clean ambo crossover form, a settled end the relax(50) canonical form,
   * the jitterbug bridge's octahedron end the clean ambo(seed) octahedron.
   */
  HS_COLD_MEMBER PolyMesh node_mesh_at(const ConwayGraph::EdgeSpec &e,
                                       bool to_end, Arena &a, Arena &b) {
    float t = to_end ? e.t_to : e.t_from;
    PolyMesh seed;
    MeshOps::clone(seed_base_, seed, a);
    Solids::SolidBuilder builder(std::move(seed), a, b);
    if (!ConwayGraph::is_platonic(e.seed_solid))
      builder.ambo();
    if (ConwayGraph::is_jitterbug_edge(e) && to_end) {
      // snub(tetra, 0.5, -pi/3) is the octahedron with 12 pairwise-coincident
      // vertices; ambo(tetra) is the same octahedron with clean topology.
      builder.ambo();
    } else if (t > 0.0f) {
      switch (e.op) {
      case ConwayGraph::MorphOp::TRUNCATE:
        builder.truncate(t);
        break;
      case ConwayGraph::MorphOp::EXPAND:
        builder.expand(t);
        break;
      case ConwayGraph::MorphOp::SNUB:
        builder.snub(t, to_end ? e.twist_to : e.twist_from);
        break;
      }
      if (e.settle && to_end)
        builder.relax(50);
    }
    return builder.build();
  }

  /**
   * @brief Camera-rotates and rasterizes one mesh, coloring each face by its
   * topology class and shading edges by distance.
   * @param canvas Target canvas to draw into.
   * @param mesh Source mesh in model space.
   * @param topology Per-face topology-class indices.
   * @param star_by_slot Resolved palette LUT per class slot for star faces
   * (emission index < node_faces_); see resolve_hankin_slot_luts.
   * @param strap_by_slot Resolved palette LUT per class slot for strap faces.
   * @param opacity Output alpha in [0, 1].
   * @param strap_open_fade Alpha multiplier for strap (emission index >=
   * node_faces_) faces over the opening window, in [0, 1]. Newborn straps split
   * former star interiors in one frame; fading their coverage in turns that
   * birth into a bounded reveal instead of a hard sliver pop.
   * @param strap_close_blend Blend of strap faces toward the rim color of the
   * base face they collapse onto, over the closing window, in [0, 1]. A strap
   * renders its full ramp even one pixel wide, so its interior differs from the
   * rim it closes onto and it winks out on vanishing; lerping toward that rim
   * first dissolves it into the face. resolve_host_faces resolves the host.
   * @param strap_terminal_fade Alpha multiplier over the last frames of the
   * close, in [0, 1]. The final sliver is anti-aliasing only, so fading it
   * uncovers nothing while removing the last-frame snap.
   * @details All three are 1.0 (no-op) through mid-cycle and the opening and
   * closing windows are disjoint. At either angle-0 bookend the straps are
   * zero-area, so a 0 changes no pixels and the star-face bookend stays
   * bitwise exact. When the LUT sets match and no shaping is active, the
   * shader skips the per-fragment role select.
   */
  void draw_mesh(Canvas &canvas, const MeshState &mesh,
                 const ArenaVector<int> &topology,
                 const BakedPalette *const (&star_by_slot)[NUM_PALETTES],
                 const BakedPalette *const (&strap_by_slot)[NUM_PALETTES],
                 float opacity, float strap_open_fade = 1.0f,
                 float strap_close_blend = 1.0f,
                 float strap_terminal_fade = 1.0f) {
    if (mesh.vertices.is_empty() || opacity < 0.01f)
      return;
    HS_PROFILE(hk_draw_mesh);

    ScratchScope scratch_a_guard(scratch_arena_a);
    MeshState rotated_mesh;
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(hk_mesh_transform);
      MeshOps::transform(mesh, rotated_mesh, scratch_arena_a, camera);
    }

    SlotLutView star_view{star_by_slot};
    SlotLutView strap_view{strap_by_slot};
    const bool fade_straps = strap_open_fade < 1.0f;
    const bool close_straps =
        strap_close_blend < 1.0f || strap_terminal_fade < 1.0f;
    bool split = fade_straps || close_straps;
    for (int s = 0; s < NUM_PALETTES; ++s)
      split |= star_by_slot[s] != strap_by_slot[s];
    const int star_faces = static_cast<int>(node_faces_);

    auto star_shader = [&](const Vector &, Fragment &f) {
      f.color = shade_mesh_topology(
          f, topology.data(), static_cast<int>(topology.size()), star_view,
          SLOT_IDENTITY, params.intensity, opacity);
    };
    auto split_shader = [&](const Vector &, Fragment &f) {
      const bool is_strap = static_cast<int>(f.v2) >= star_faces;
      const SlotLutView &view = is_strap ? strap_view : star_view;
      f.color = shade_mesh_topology(f, topology.data(),
                                    static_cast<int>(topology.size()), view,
                                    SLOT_IDENTITY, params.intensity, opacity);
      if (is_strap && fade_straps)
        f.color.alpha *= strap_open_fade;
      if (is_strap && close_straps) {
        const int fi = static_cast<int>(f.v2);
        if (fi >= 0 && fi < static_cast<int>(MAX_HANKIN_FACES)) {
          Color4 rim =
              palette_bank_.bank.entries[host_face_palette_[fi]].get(0.0f);
          rim.alpha = f.color.alpha;
          f.color = rim.lerp(f.color, strap_close_blend);
        }
        f.color.alpha *= strap_terminal_fade;
      }
    };

    {
      HS_PROFILE(hk_mesh_scan);
      if (split)
        Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, split_shader,
                               scratch_arena_a, params.debug_bb);
      else
        Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, star_shader,
                               scratch_arena_a, params.debug_bb);
    }
  }

  /**
   * @brief Camera-rotates and rasterizes a morph frame's swept mesh, shading
   * each face from its pre-blended palette ramp.
   * @param canvas Target canvas to draw into.
   * @param mesh Compiled swept mesh (scratch-backed, this frame only).
   * @param shading Per-face blended-ramp table from the ConwayMorph.
   */
  void draw_conway_mesh(Canvas &canvas, const MeshState &mesh,
                        const Animation::ConwayMorph::Shading &shading) {
    if (mesh.vertices.is_empty())
      return;
    HS_PROFILE(hk_draw_mesh);

    ScratchScope scratch_a_guard(scratch_arena_a);
    MeshState rotated_mesh;
    OrientTransformer camera(orientation);
    {
      HS_PROFILE(hk_mesh_transform);
      MeshOps::transform(mesh, rotated_mesh, scratch_arena_a, camera);
    }

    auto fragment_shader = [&](const Vector &, Fragment &f) {
      float t = hs::clamp(fragment_edge_dist(f) * params.intensity, 0.0f, 1.0f);
      int face = static_cast<int>(f.v2);
      int ramp = (face >= 0 && face < static_cast<int>(shading.faces))
                     ? shading.face_ramp[face]
                     : 0;
      Color4 c = shading.ramps[ramp].get(t);
      c.alpha = 1.0f;
      f.color = c;
    };

    {
      HS_PROFILE(hk_mesh_scan);
      Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, fragment_shader,
                             scratch_arena_a, params.debug_bb);
    }
  }

  /**
   * @brief Schedules one interlace-angle sweep plus the sprite that
   * re-evaluates and draws the mesh each frame.
   * @details The sweep starts one frame after the sprite and the sprite runs
   * one frame longer, so the first drawn frame renders the exact angle-0
   * bookend the leg completion pinned, and the sweep's completion pins the
   * closing bookend before the sprite's final draw. Both are gated on the
   * same pause flag so grabbing the slider holds the frame instead of
   * blanking it.
   */
  /**
   * @brief Resolves, per hankin-added face, the palette of the base face it
   * lives inside — the rim it opens from and collapses back onto.
   * @details A rosette is born on a base face's interior walls and opens inward
   * from them, so it belongs wholly to one base face; the host is found by
   * shared vertices, which is a topological fact and therefore angle
   * independent. A centroid test would misfile it at the closed end, where the
   * rosette has collapsed outward onto the rim and its centroid sits on the
   * boundary. Runs once per cycle; record_node_palettes must precede it.
   */
  HS_COLD_MEMBER void resolve_host_faces() {
    const size_t faces = compiled_hankin.face_counts.size();
    HS_CHECK(faces <= MAX_HANKIN_FACES,
             "HankinSolids: hankin face count exceeds the host-palette table");
    // face_offsets are not stored on CompiledHankin; walk the flat face list.
    size_t off = 0;
    size_t star_off[MAX_NODE_FACES];
    for (size_t f = 0; f < faces && f < node_faces_; ++f) {
      star_off[f] = off;
      off += compiled_hankin.face_counts[f];
    }
    for (size_t f = node_faces_; f < faces; ++f) {
      const int n = compiled_hankin.face_counts[f];
      size_t best = 0;
      int best_shared = -1;
      for (size_t j = 0; j < node_faces_; ++j) {
        const int nj = compiled_hankin.face_counts[j];
        int shared = 0;
        for (int a = 0; a < n; ++a)
          for (int b = 0; b < nj; ++b)
            if (compiled_hankin.faces[off + a] ==
                compiled_hankin.faces[star_off[j] + b])
              ++shared;
        if (shared > best_shared) {
          best_shared = shared;
          best = j;
        }
      }
      host_face_palette_[f] = node_face_palette_[best];
      off += n;
    }
  }

  HS_COLD_MEMBER void start_hankin_cycle() {
    constexpr int DURATION = 64;
    hankin_cycle_frame_ = 0;
    resolve_host_faces();
    timeline.add(2, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_linear, false,
                                        &anims_paused_)
                        .then([this]() {
                          // Bookend-in: the sweep's final sample lands ~0.002
                          // rad off the flat p_corner branch; force exact 0 so
                          // the sprite's last draw is the base solid.
                          params.hankin_angle = 0.0f;
                          this->start_morph_cycle();
                        }));

    // Snapshot the angle-independent counts for the per-frame HS_CHECK below.
    hankin_vertex_count_ = compiled_hankin.static_vertices.size() +
                           compiled_hankin.dynamic_vertices.size();
    hankin_face_count_ = compiled_hankin.face_counts.size();
    timeline.add(
        0, Animation::Sprite(
               [this](Canvas &c, float opacity) {
                 // update_hankin re-binds the mesh's vectors against
                 // persistent_arena every frame; the angle never changes the
                 // vertex/face counts, so bind reuses the blocks in place.
                 {
                   HS_PROFILE(hk_update_hankin);
                   MeshOps::update_hankin(compiled_hankin, mesh_,
                                          persistent_arena,
                                          params.hankin_angle);
                 }
                 // Always-on guard: grown counts would leak persistent_arena
                 // every frame on a permanent install.
                 HS_CHECK(mesh_.vertices.size() == hankin_vertex_count_ &&
                              mesh_.face_counts.size() == hankin_face_count_,
                          "HankinSolids: per-frame mesh counts changed; the "
                          "persistent re-bind would grow the arena");
                 const int cycle_frame = hankin_cycle_frame_;
                 if (!anims_paused_)
                   ++hankin_cycle_frame_;
                 // Blended strap LUTs live in scratch_b for this frame only;
                 // the draw path below uses scratch_a exclusively.
                 ScratchScope blend_guard(scratch_arena_b);
                 BakedPalette blended[NUM_PALETTES];
                 const BakedPalette *star_by_slot[NUM_PALETTES];
                 const BakedPalette *strap_by_slot[NUM_PALETTES];
                 resolve_hankin_slot_luts(cycle_frame, blended, star_by_slot,
                                          strap_by_slot, scratch_arena_b);
                 // Straps reveal over the opening window and dissolve into the
                 // host face's rim over the closing one, with the terminal
                 // sliver faded out. Every weight is exactly 0 at its bookend
                 // and 1 through mid-cycle, so mid-cycle draws are unchanged.
                 const int to_close = DURATION - cycle_frame;
                 draw_mesh(c, mesh_, mesh_.topology, star_by_slot,
                           strap_by_slot, opacity,
                           strap_blend_weight(cycle_frame),
                           strap_blend_weight(to_close),
                           hs::clamp(static_cast<float>(to_close) /
                                         STRAP_TERMINAL_FRAMES,
                                     0.0f, 1.0f));
               },
               DURATION + 1, 0, ease_linear, 0, ease_linear, &anims_paused_));
  }

  /**
   * @brief Picks the next graph edge, reconciles the held seed, and schedules
   * the ConwayMorph leg.
   */
  HS_COLD_MEMBER void start_morph_cycle() {
    using namespace ConwayGraph;
#ifdef HS_PROFILE_ORDERED_CYCLE
    cur_edge_ = pick_next_edge_ordered(node_, cur_edge_, leg_counter_++);
#else
    cur_edge_ = pick_next_edge(node_, cur_edge_, legs_in_family_, node_visits_,
                               static_cast<uint32_t>(hs::random()()));
#endif
    const EdgeSpec &e = EDGES[cur_edge_];
    HS_CHECK(edge_touches(cur_edge_, node_));
    reverse_ = (e.to_node == node_);

    SeedFix fix = seed_fix_at_start(cur_edge_, seed_identity_);
    HS_CHECK(fix != SeedFix::INVALID,
             "HankinSolids: no seed reconciliation for the picked edge");
    if (fix == SeedFix::DUAL_SWAP) {
      // Ambo crossover: ambo(dual(seed)) == ambo(seed), so the swap is
      // pixel-invisible at the displayed t = 0.5 form.
      HS_CHECK(node_ == CUBOCTAHEDRON || node_ == ICOSIDODECAHEDRON);
      generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        seed_base_ =
            Solids::finalize_solid(MeshOps::dual(seed_base_, a, b), target);
      });
      seed_identity_ = static_cast<uint8_t>(dual_platonic(seed_identity_));
    } else if (fix == SeedFix::REGEN_TETRA) {
      // Reverse family bridge: the held octa/icosa was derived from the
      // registry tetrahedron, so regenerating that tetrahedron is frame-exact.
      HS_CHECK(node_ == OCTAHEDRON || node_ == ICOSAHEDRON);
      generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
        seed_base_ =
            Solids::finalize_solid(Solids::Platonic::tetrahedron(a, b), target);
      });
      seed_identity_ = TETRAHEDRON;
    }
    HS_CHECK(fix == SeedFix::DERIVE_AMBO || seed_identity_ == e.seed_solid,
             "HankinSolids: leg seed identity mismatch");

    Animation::ConwayMorph::PaletteHandoff handoff{
        &palette_bank_.bank, node_face_palette_,        node_face_sides_,
        node_faces_,         fix == SeedFix::DUAL_SWAP, node_face_centroid_};

    // Bookend grouping of the arrival node: the closing bookend displays one
    // palette per hankin star-face class, so the leg's color targets key on
    // that classification — every face it merges converges to a single color
    // before the swap. Built from the same arrival mesh finish_morph_cycle
    // rebuilds, so the class ids match at completion.
    int arrival_topo[MAX_NODE_FACES];
    size_t arrival_faces = 0;
    {
      ScratchScope ba(scratch_arena_a);
      ScratchScope bb(scratch_arena_b);
      PolyMesh arrival =
          node_mesh_at(e, !reverse_, scratch_arena_b, scratch_arena_a);
      arrival_faces = arrival.face_counts.size();
      HS_CHECK(arrival_faces <= MAX_NODE_FACES);
      CompiledHankin ch;
      MeshOps::compile_hankin(arrival, ch, scratch_arena_b, scratch_arena_a);
      MeshState hk;
      MeshOps::update_hankin(ch, hk, scratch_arena_b, 0.0f);
      // One arena serves classify's scratch and output (LIFO-stacked scopes),
      // keeping scratch_b free for the compiled hankin + mesh peak.
      MeshOps::classify_faces_by_topology(hk, scratch_arena_a, scratch_arena_a,
                                          scratch_arena_a);
      for (size_t f = 0; f < arrival_faces; ++f)
        arrival_topo[f] = hk.topology[f];
    }
    Animation::ConwayMorph::BookendClasses bookend{arrival_topo, arrival_faces};

    ScratchScope sa(scratch_arena_a);
    ScratchScope sb(scratch_arena_b);
    PolyMesh derived;
    if (fix == SeedFix::DERIVE_AMBO)
      derived = MeshOps::ambo(seed_base_, scratch_arena_a, scratch_arena_b);
    const PolyMesh &leg_seed =
        fix == SeedFix::DERIVE_AMBO ? derived : seed_base_;

    hs::log("Conway leg: '%s' -> '%s'", Solids::simple_registry[node_].name,
            Solids::simple_registry[edge_other_end(cur_edge_, node_)].name);

    Animation::ConwayMorph anim(leg_seed, e, reverse_, persistent_arena,
                                draw_conway_fn_, handoff, SWEEP_FRAMES,
                                e.settle ? SETTLE_FRAMES : 0, bookend);
    pending_landing_ = &anim.landing();
    timeline.add(
        0, std::move(anim).then([this]() { this->finish_morph_cycle(); }));
  }

  /**
   * @brief Leg completion: clean-endpoint swap, reseed, hankin rebuild from
   * the arrived mesh, forward palette mapping, compaction, next cycle.
   */
  HS_COLD_MEMBER void finish_morph_cycle() {
    using namespace ConwayGraph;
    const EdgeSpec &e = EDGES[cur_edge_];
    const Animation::ConwayMorph::Landing &landing = *pending_landing_;
    const bool arrived_at_to = !reverse_;
    const uint8_t arrived = arrived_at_to ? e.to_node : e.from_node;

    if (family(arrived) != family(node_))
      legs_in_family_ = 0;
    else
      ++legs_in_family_;
    record_visit(node_visits_, arrived);
    node_ = arrived;
    hs::log("Loading shape: '%s'", Solids::simple_registry[node_].name);

    // Outgoing cycle's display state, consumed by the strap-crossfade prep
    // below: the per-slot colors, and which slots were on screen at all.
    const std::array<int, NUM_PALETTES> prev_idx = palette_idx_;
    bool prev_used[NUM_PALETTES] = {};
    for (size_t f = 0; f < mesh_.topology.size(); ++f)
      prev_used[wrap(mesh_.topology[f], NUM_PALETTES)] = true;

    // Build the arrived base mesh from the held seed and compile the hankin
    // pattern from that mesh — never a registry regenerate, so bridge
    // arrivals keep the orientation the walk produced.
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = node_mesh_at(e, arrived_at_to, a, b);
      if (e.reseed == Reseed::ADOPT && is_platonic(arrived)) {
        if (arrived_at_to) {
          // Family bridge: the arrived solid becomes the new family seed.
          seed_base_ = Solids::finalize_solid(base, target);
          seed_identity_ = node_;
        } else if (arrived == ICOSAHEDRON) {
          // Reverse jitterbug arrival: hold the icosahedron node's canonical
          // relax form (the mesh the tetra -> icosa bridge adopts), not the
          // unrelaxed jitterbug form the bookend displays.
          PolyMesh s;
          MeshOps::clone(seed_base_, s, a);
          seed_base_ =
              Solids::finalize_solid(Solids::SolidBuilder(std::move(s), a, b)
                                         .snub(0.5f, SNUB_BRIDGE_TWIST)
                                         .relax(50)
                                         .build(),
                                     target);
          seed_identity_ = node_;
        }
      }
      record_node_faces(base);

      compiled_hankin = CompiledHankin();
      MeshOps::compile_hankin(base, compiled_hankin, target, a);
      mesh_.clear();
      MeshOps::update_hankin(compiled_hankin, mesh_, target, 0.0f);
    });
    classify_mesh_topology(mesh_);

    // Forward palette mapping: base faces fill hankin star faces 1:1
    // (emission identity), so each populated slot carries the leg's landed
    // palette verbatim; only newborn (rosette-only) slots take fresh shuffles.
    HS_CHECK(node_faces_ <= landing.faces);
    MeshPaletteBank::shuffle_indices(palette_idx_);
    bool slot_mapped[NUM_PALETTES] = {};
    for (size_t f = 0; f < node_faces_; ++f) {
      // The leg's targets keyed on this same classification (computed at leg
      // start from the same arrival mesh); drift here would pop the bookend.
      HS_CHECK(landing.topology[f] == mesh_.topology[f],
               "HankinSolids: arrival classification drifted across the leg");
      int slot = wrap(mesh_.topology[f], NUM_PALETTES);
      if (!slot_mapped[slot]) {
        slot_mapped[slot] = true;
        palette_idx_[slot] =
            landing.to_palette[wrap(landing.topology[f], NUM_PALETTES)];
      }
    }
    record_node_palettes();
    prepare_strap_crossfade(prev_idx, prev_used);
    pending_landing_ = nullptr;

    // Bookend-out: the next cycle's first drawn sample is exactly angle 0.
    params.hankin_angle = 0.0f;

    {
      // Survivors split across the scratch pair: scratch_b holds the two
      // largest (compiled hankin + palette bank), scratch_a the rest.
      Persist<CompiledHankin> ph(compiled_hankin, scratch_arena_b,
                                 persistent_arena);
      Persist<MeshState> pf(mesh_, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank_, scratch_arena_b,
                                  persistent_arena);
      Persist<PolyMesh> ps(seed_base_, scratch_arena_a, persistent_arena);
      persistent_arena.reset();
    }

    MeshOps::update_hankin(compiled_hankin, mesh_, persistent_arena,
                           params.hankin_angle);
    start_hankin_cycle();
  }

  MeshState mesh_; /**< The single on-screen mesh (hankin form). */
  CompiledHankin compiled_hankin; /**< Active during the hankin cycle. */
  PolyMesh seed_base_;            /**< Held Platonic seed of the graph walk. */
  std::array<int, NUM_PALETTES> palette_idx_ =
      {}; /**< Class slot -> palette; value-init so a missed shuffle reads 0,
             not garbage. */
  std::array<int, NUM_PALETTES> strap_from_ =
      {}; /**< Per-slot crossfade origin palette for the current cycle's
             opening window (equals palette_idx_ on unarmed slots). */
  uint8_t strap_blend_mask_ = 0; /**< Bit s: slot s's strap faces crossfade
                                    this cycle. */
  int hankin_cycle_frame_ = 0;   /**< Sprite draws since the opening bookend
                                    of the active hankin cycle. */

  uint8_t node_face_palette_[MAX_NODE_FACES] =
      {}; /**< Displayed palette per node base face. */
  uint8_t host_face_palette_[MAX_HANKIN_FACES] =
      {}; /**< Per hankin-added face, the palette of the base face it lives
             inside — the rim color it collapses onto. See
             resolve_host_faces. */
  uint8_t node_face_sides_[MAX_NODE_FACES] =
      {}; /**< Clean side count per node base face. */
  Vector node_face_centroid_[MAX_NODE_FACES] =
      {};                 /**< Unit centroid per node base face (geometric
                             palette provenance). */
  size_t node_faces_ = 0; /**< Face count of the current node's base mesh. */

  uint8_t node_ = 0; /**< Current graph node (simple-registry index). */
  uint8_t seed_identity_ = 0; /**< Platonic solid seed_base_ represents. */
  int cur_edge_ = -1;         /**< Edge of the last (or in-flight) leg. */
  bool reverse_ = false;      /**< In-flight leg runs to_node -> from_node. */
  int legs_in_family_ =
      0; /**< Legs since the last family change (walk weighting). */
  uint8_t node_visits_[ConwayGraph::NUM_NODES] =
      {}; /**< Aged per-node visit counts (walk recency weighting). */
  uint32_t leg_counter_ = 0; /**< Monotonic leg count (ordered-cycle picks). */
  const Animation::ConwayMorph::Landing *pending_landing_ =
      nullptr; /**< In-flight leg's arrival data (leg-arena backed). */

  size_t hankin_vertex_count_ =
      0; /**< Expected per-frame vertex count for the active cycle. */
  size_t hankin_face_count_ =
      0; /**< Expected per-frame face count for the active cycle. */

  /**
   * @brief Draw callback for morph frames.
   * @details Held as a member for stable FunctionRef lifetime.
   */
  Fn<void(Canvas &, const MeshState &, const Animation::ConwayMorph::Shading &),
     8>
      draw_conway_fn_{[this](Canvas &c, const MeshState &m,
                             const Animation::ConwayMorph::Shading &sh) {
        draw_conway_mesh(c, m, sh);
      }};

  Orientation<> orientation; /**< Current camera orientation. */
  FastNoiseLite noise;       /**< Noise source driving the orientation walk. */
  Timeline timeline;         /**< Schedules sweeps, sprites, and morphs. */
  Pipeline<W, H> filters;    /**< Per-pixel filter pipeline applied on draw. */

  /**
   * @brief User-adjustable rendering parameters.
   */
  struct Params {
    float intensity = 1.2f;           /**< Edge-distance shading gain. */
    float hankin_angle = PI_F / 4.0f; /**< Interlace angle in radians. */
    bool debug_bb = false; /**< Draw face bounding boxes when true. */
  } params;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(HankinSolids)
