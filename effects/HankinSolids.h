/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include "core/effects_engine.h"

#include <algorithm>

#ifdef CORE_TEENSY
#define OS_CYCLES() ARM_DWT_CYCCNT
#else
#define OS_CYCLES() 0
#endif

template <int W, int H> class HankinSolids : public Effect {
public:
  FLASHMEM HankinSolids() : Effect(W, H), filters() {}

  void init() override {
    // scratch_b must hold CompiledHankin (~10KB) + BakedPaletteBank (~15KB)
    configure_arenas(GLOBAL_ARENA_SIZE - 16 * 1024 - 32 * 1024, 16 * 1024,
                     32 * 1024);
    registerParam("Intensity", &params.intensity, 0.0f, 5.0f);
    registerParam("Angle", &params.hankin_angle, 0.0f, PI_F / 2.0f);
    registerParam("Debug BB", &params.debug_bb);

    timeline.add(0, Animation::RandomWalk<W>(
                        orientation, Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));

    solid_idx = 0;

    for (int i = 0; i < NUM_PALETTES; ++i)
      baked_palette_bank_.entries[i].bake(persistent_arena,
                                          *source_palettes[i]);

    load_shape(carousel.current(), compiled_hankin,
               palettes_slots[carousel.front_index()], solid_idx,
               params.hankin_angle);

    start_hankin_cycle();
  }

  bool show_bg() const override { return false; }
  void draw_frame() override {
    uint32_t t_start = OS_CYCLES();

    Canvas canvas(*this);
    uint32_t t1 = OS_CYCLES();
    timeline.step(canvas);
    c_timeline += (OS_CYCLES() - t1);

    c_total_accum += (OS_CYCLES() - t_start);
    // Accumulate scan sub-metrics from this frame
    a_sdf_dist += hs::g_scan_metrics.sdf_dist;
    a_frag_shader += hs::g_scan_metrics.frag_shader;
    a_plot += hs::g_scan_metrics.plot;
    a_face_setup += hs::g_scan_metrics.face_setup;
    a_scan_loop += hs::g_scan_metrics.scan_loop;
    a_pixels_tested += hs::g_scan_metrics.pixels_tested;
    a_pixels_culled += hs::g_scan_metrics.pixels_culled;
    a_lut_hits += hs::g_scan_metrics.lut_hits;
    a_exact_hits += hs::g_scan_metrics.exact_hits;
    hs::g_scan_metrics.reset();
    current_cycle_frames++;
  }

private:
  static constexpr int NUM_PALETTES = 5;
  static constexpr std::array<const ProceduralPalette *, NUM_PALETTES>
      source_palettes = {&Palettes::embers, &Palettes::richSunset,
                         &Palettes::brightSunrise, &Palettes::bruisedMoss,
                         &Palettes::lavenderLake};
  BakedPaletteBank baked_palette_bank_;

  uint64_t c_raster = 0;
  uint64_t c_pipeline = 0;
  uint64_t c_timeline = 0;
  uint64_t c_total_accum = 0;
  uint64_t a_sdf_dist = 0;
  uint64_t a_frag_shader = 0;
  uint64_t a_plot = 0;
  uint64_t a_face_setup = 0;
  uint64_t a_scan_loop = 0;
  uint64_t a_pixels_tested = 0;
  uint64_t a_pixels_culled = 0;
  uint64_t a_lut_hits = 0;
  uint64_t a_exact_hits = 0;
  uint32_t current_cycle_frames = 0;

  void print_and_reset_metrics(const char* phase_name) {
    uint32_t f = current_cycle_frames ? current_cycle_frames : 1;
    printf("\n=== %s Complete (%lu frames) ===\n", phase_name, f);
    printf("Avg/frame: Total=%lu, Timeline=%lu\n", 
           (uint32_t)(c_total_accum / f), (uint32_t)(c_timeline / f));
    printf("  Pipeline(Logic)=%lu\n", (uint32_t)(c_pipeline / f));
    printf("  Raster Breakdown:\n");
    printf("    FaceSetup=%lu, ScanLoop=%lu\n",
           (uint32_t)(a_face_setup / f), (uint32_t)(a_scan_loop / f));
    printf("    SDF_Dist=%lu, FragShader=%lu, Plot=%lu\n",
           (uint32_t)(a_sdf_dist / f), (uint32_t)(a_frag_shader / f), (uint32_t)(a_plot / f));
    uint64_t accounted = a_sdf_dist + a_frag_shader + a_plot;
    uint64_t scan_other = a_scan_loop > accounted ? a_scan_loop - accounted : 0;
    printf("    ScanOther(bounds+iter)=%lu\n", (uint32_t)(scan_other / f));
    printf("  Pixels: tested=%lu, culled=%lu, interior=%lu\n",
           (uint32_t)(a_pixels_tested / f), (uint32_t)(a_pixels_culled / f),
           (uint32_t)((a_pixels_tested - a_pixels_culled) / f));
    printf("  LUT_hits=%lu, Exact_hits=%lu\n",
           (uint32_t)(a_lut_hits / f), (uint32_t)(a_exact_hits / f));
    
    c_raster = 0; c_pipeline = 0; c_timeline = 0;
    c_total_accum = 0;
    a_sdf_dist = 0; a_frag_shader = 0; a_plot = 0;
    a_face_setup = 0; a_scan_loop = 0;
    a_pixels_tested = 0; a_pixels_culled = 0;
    a_lut_hits = 0; a_exact_hits = 0;
    current_cycle_frames = 0;
    hs::g_scan_metrics.reset();
  }

  PolyMesh generate_base_solid(int idx, Arena &a, Arena &b) {
    auto solids = Solids::Collections::get_simple_solids();
    hs::log("Loading shape: '%s'", solids[idx].name);
    return Solids::finalize_solid(solids[idx].generate(a, b), a);
  }

  void classify_mesh_topology(MeshState &mesh) {
    scratch_arena_a.reset();
    scratch_arena_b.reset();
    MeshOps::classify_faces_by_topology(mesh, scratch_arena_a, scratch_arena_b,
                                        persistent_arena);
  }

  void shuffle_palette_indices(std::array<int, NUM_PALETTES> &out) {
    for (int i = 0; i < NUM_PALETTES; ++i)
      out[i] = i;
    std::shuffle(out.begin(), out.end(), hs::random());
  }

  /// Full pipeline: generate → compile → evaluate → classify → palette indices.
  void load_shape(MeshState &out_mesh, CompiledHankin &out_hankin,
                  std::array<int, NUM_PALETTES> &out_palette_idx, int idx,
                  float angle) {
    generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = generate_base_solid(idx, a, b);

      out_hankin = CompiledHankin();
      MeshOps::compile_hankin(base, out_hankin, target, a);

      out_mesh.clear();
      MeshOps::update_hankin(out_hankin, out_mesh, target, angle);
    });

    classify_mesh_topology(out_mesh);
    shuffle_palette_indices(out_palette_idx);
  }

  void draw_mesh(Canvas &canvas, const MeshState &mesh,
                 const ArenaVector<int> &topology,
                 const std::array<int, NUM_PALETTES> &palette_idx,
                 float opacity, const Quaternion &q) {
    if (mesh.vertices.empty() || opacity < 0.01f)
      return;

    ScratchScope _(scratch_arena_a);
    MeshState rotated_mesh;
    
    uint32_t t0 = OS_CYCLES();
    MeshOps::transform(mesh, rotated_mesh, scratch_arena_a);

    for (auto &v : rotated_mesh.vertices)
      v = rotate(v, q);
    
    c_pipeline += (OS_CYCLES() - t0);

    auto fragment_shader = [&](const Vector &p, Fragment &f) {
      int faceIdx = (int)std::round(f.v2);
      int topoIdx = (faceIdx >= 0 && faceIdx < (int)topology.size())
                        ? topology[faceIdx]
                        : 0;

      float distFromEdge = -f.v1;
      float normalizedDist =
          (f.size > 0.0001f) ? (distFromEdge / f.size) : 0.0f;
      float t = hs::clamp(normalizedDist * params.intensity, 0.0f, 1.0f);

      f.color =
          baked_palette_bank_.entries[palette_idx[topoIdx % NUM_PALETTES]].get(
              t);
      f.color.alpha *= opacity;
    };

    uint32_t t1 = OS_CYCLES();
    Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, fragment_shader,
                           scratch_arena_a, params.debug_bb);
    c_raster += (OS_CYCLES() - t1);
  }

  void promote_staged_hankin() {
    compiled_hankin = std::move(compiled_hankin_staging);
    compiled_hankin_staging = CompiledHankin();
  }

  void start_hankin_cycle() {
    constexpr int DURATION = 64;
    timeline.add(0, Animation::Mutation(params.hankin_angle,
                                        sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f),
                                        DURATION, ease_mid, false)
                        .then([this]() {
                          this->print_and_reset_metrics("Hankin Cycle");
                          this->start_morph_cycle(); 
                        }));

    int front = carousel.front_index();
    timeline.add(
        0, Animation::Sprite(
               [this, front](Canvas &c, float opacity) {
                 uint32_t t0 = OS_CYCLES();
                 MeshOps::update_hankin(compiled_hankin, carousel.slot(front),
                                        persistent_arena, params.hankin_angle);
                 c_pipeline += (OS_CYCLES() - t0);
                 draw_mesh(c, carousel.slot(front),
                           carousel.slot(front).topology, palettes_slots[front],
                           opacity, orientation.get());
               },
               DURATION));
  }

  FLASHMEM void start_morph_cycle() {
    constexpr int MORPH_FRAMES = 16;
    auto solids = Solids::Collections::get_simple_solids();
    int next_idx = (solid_idx + 1) % solids.size();

    int old_front = carousel.front_index();
    int new_slot = 1 - old_front;

    load_shape(carousel.slot(new_slot), compiled_hankin_staging,
               palettes_slots[new_slot], next_idx, params.hankin_angle);

    // Set slot indices before creating MeshMorph (draw callbacks reference
    // these)
    morph_old_slot_ = old_front;
    morph_new_slot_ = new_slot;

    timeline.add(
        0, Animation::MeshMorph(
               carousel.slot(old_front), carousel.slot(new_slot),
               persistent_arena, draw_morph_outgoing_fn_,
               draw_morph_incoming_fn_, MORPH_FRAMES, ease_in_out_sin)
               .then([this, next_idx, new_slot]() {
                 this->print_and_reset_metrics("Morph Cycle");
                 solid_idx = next_idx;
                 carousel.set_front(new_slot);
                 promote_staged_hankin();
                 // Manual compaction: preserve both carousel slots +
                 // compiled_hankin
                 {
                   Persist<CompiledHankin> ph(compiled_hankin, scratch_arena_b,
                                              persistent_arena);
                   Persist<MeshState> p0(carousel.slot(0), scratch_arena_a,
                                         persistent_arena);
                   Persist<MeshState> p1(carousel.slot(1), scratch_arena_a,
                                         persistent_arena);
                   Persist<BakedPaletteBank> pp(
                       baked_palette_bank_, scratch_arena_b, persistent_arena);
                   persistent_arena.reset();
                   hs::log("morph_cycle_then: finished arena compaction");
                 }

                 MeshOps::update_hankin(compiled_hankin, carousel.current(),
                                        persistent_arena, params.hankin_angle);
                 start_hankin_cycle();
               }));
  }

  MeshCarousel<W> carousel;
  CompiledHankin compiled_hankin;         // Active during hankin cycle
  CompiledHankin compiled_hankin_staging; // Built during morph cycle
  std::array<int, NUM_PALETTES> palettes_slots[2];

  // Slot indices for morph draw callbacks (set before each morph)
  int morph_old_slot_ = 0;
  int morph_new_slot_ = 1;

  /// Morph draw callbacks — members for stable FunctionRef lifetime
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_outgoing_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_old_slot_).topology,
                  palettes_slots[morph_old_slot_], o, orientation.get());
      }};
  Fn<void(Canvas &, const MeshState &, float), 8> draw_morph_incoming_fn_{
      [this](Canvas &c, const MeshState &m, float o) {
        draw_mesh(c, m, carousel.slot(morph_new_slot_).topology,
                  palettes_slots[morph_new_slot_], o, orientation.get());
      }};

  Orientation<W> orientation;
  FastNoiseLite noise;
  Timeline<W> timeline;
  Pipeline<W, H> filters;
  int solid_idx = 0;

  struct Params {
    float intensity = 1.2f;
    float hankin_angle = PI_F / 4.0f;
    bool debug_bb = false;
  } params;
};

#include "core/effect_registry.h"
REGISTER_EFFECT(HankinSolids)
