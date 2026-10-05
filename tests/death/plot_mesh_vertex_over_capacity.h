/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/**
 * @brief Death case: a face vertex index past the edge-dedup bitset must trap.
 * @details Plot surface — a vertex index beyond the TriangularBitset<128>
 *          capacity makes the face-walk draw() overload trap on the cold
 *          per-edge setup path instead of silently dropping the edge, which
 *          would leave a wireframe with missing lines and mask the sizing bug.
 */
inline void case_plot_mesh_vertex_over_capacity() {
  constexpr int W = 32, H = 16;
  OverCapacityMockMesh mesh;
  DeathEffect fx;
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::Mesh::draw<W, H>(pipe, c, mesh, [](const math::Vector &, Fragment &) {
  }); // index 130 -> trap
}

/** @brief An open face dual cannot supply a complete wireframe. */
inline void case_plot_four_regular_open_mesh() {
  configure_arenas_default();
  MeshState mesh;
  build_meshstate_solid<Solids::Octahedron>(mesh, persistent_arena,
                                            opaque(true));
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, mesh.faces.size());
  Plot::Mesh::extract_four_regular_edges(mesh, edges, scratch_arena_a);
}

/** @brief A tetrahedron's dual has odd cycles and cannot be two-colored. */
inline void case_plot_four_regular_non_bipartite() {
  configure_arenas_default();
  MeshState mesh;
  build_meshstate_solid<Solids::Tetrahedron>(mesh, persistent_arena);
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, mesh.faces.size());
  Plot::Mesh::extract_four_regular_edges(mesh, edges, scratch_arena_a);
}

/** @brief Face-edge lookup rejects an absent pair. */
inline void case_plot_find_missing_edge() {
  configure_arenas_default();
  ArenaVector<Plot::Mesh::Edge> edges(persistent_arena, 1);
  edges.push_back({opaque<uint16_t>(0), opaque<uint16_t>(1)});
  const auto index = Plot::Mesh::find_edge_index(edges, opaque<uint16_t>(1),
                                                 opaque<uint16_t>(2));
  if (index == opaque<uint16_t>(42))
    std::printf("x");
}

/**
 * @brief Death case: extract_edges with an over-capacity vertex index must trap.
 * @details Plot surface — the precomputed-edge path traps on the same cold setup
 *          path as the face-walk draw() overload, rather than silently filtering
 *          the edge out (which would produce an edge list with missing lines and
 *          mask the sizing bug).
 */
inline void case_plot_extract_edges_vertex_over_capacity() {
  OverCapacityMockMesh mesh;
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(scratch_arena_a, 8);
  Plot::Mesh::extract_edges(mesh, edges); // index 130 -> trap
}

/**
 * @brief Death case: a plot rejects a canvas that is not its <W, H>.
 * @details Plot surface -- the fragment walk projects onto the <W, H> grid and
 *          the sink strides the framebuffer by its own W, so a canvas of a
 *          different size writes past the row it means to.
 */
inline void case_plot_canvas_dim_mismatch() {
  constexpr int W = 32, H = 16;
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 2);
  Fragment f;
  f.pos = math::Vector(1, 0, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 1, 0);
  points.push_back(f);

  DeathEffect fx(W, opaque(H + 1));
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(pipe, c, points,
                        [](const math::Vector &, Fragment &) {});
}

/**
 * @brief Death case: a plot window over a multi-segment polyline must trap.
 * @details Plot surface -- the window narrows one segment's arc fraction, so a
 *          multi-segment polyline would apply the same [start, end] to every
 *          segment and silently drop whole edges rather than the intended
 *          out-of-band tail.
 */
inline void case_plot_window_multi_segment() {
  constexpr int W = 32, H = 16;
  ArenaVector<Fragment> points;
  points.bind(scratch_arena_a, 4);
  Fragment f;
  f.pos = math::Vector(1, 0, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 1, 0);
  points.push_back(f);
  f.pos = math::Vector(0, 0, 1);
  points.push_back(f); // 2 segments under a window -> HS_CHECK

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(
      pipe, c, points, [](const math::Vector &, Fragment &) {},
      {.plot_t_start = opaque(0.25f), .plot_t_end = opaque(0.75f)});
}

/**
 * @brief Death case: a feedback downsample that doesn't divide the resolution must trap.
 * @details Filter surface — Pixel::Feedback::flush traps rather than silently
 *          turning the whole feedback effect into a no-op; a cold
 *          authoring/config error the project routes to HS_CHECK (enabled
 *          remains the supported way to switch feedback off). The trap fires
 *          before any_pixel_lit / scratch allocation, so no buffers needed.
 */
inline void case_feedback_downsample_indivisible() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  Canvas c(fx);
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.downsample = opaque(5); // 32 % 5 != 0 -> HS_CHECK
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  pipe.begin_frame(c, 1.0f);
}

inline void case_feedback_uncached_scratch_budget() {
  DeathEffect fx;
  {
    Canvas first(fx);
    first(0, 0) = Pixel(65535, 65535, 65535);
  }
  fx.advance_display();
  Canvas canvas(fx);
  static uint8_t storage[16];
  scratch_arena_a.rebind(storage, sizeof(storage));
  ::Feedback::Style style{};
  Pipeline<32, 16, Filter::Pixel::Feedback<32, 16>> pipe{
      Filter::Pixel::Feedback<32, 16>(style)};
  pipe.begin_frame(canvas, 1.0f);
}

/**
 * @brief Death case: retuning a screen trail to a non-positive lifetime must trap.
 * @details Filter surface — Screen::Trails::set_lifetime carries the
 *          constructor's bound, so a slider that reaches zero traps here rather
 *          than dividing by it in flush()'s fade progress.
 */
inline void case_screen_trails_set_lifetime_nonpositive() {
  Filter::Screen::Trails<8> trails(4);
  trails.set_lifetime(opaque(0));
}

/**
 * @brief Death case: a world trail lifetime past the ttl byte must trap.
 * @details Filter surface — World::Trails packs the remaining lifetime into a
 *          uint8_t ttl, so a lifetime above 255 would wrap on seeding and give
 *          a near-dead trail instead of the long one asked for.
 */
inline void case_world_trails_lifetime_over_max() {
  Filter::World::Trails<8> trails(opaque(256)); // lifetime > 255 -> HS_CHECK
  (void)trails;
}

/**
 * @brief Death case: seeding a screen trail before init_storage() must trap.
 * @details Filter surface — plot() seeds the ring buffer, so without storage it
 *          would silently drop every trail point and an effect that never
 *          flushes would render trail-free instead of failing.
 */
inline void case_screen_trails_plot_without_storage() {
  Filter::Screen::Trails<8> trails(4); // no init_storage() -> HS_CHECK
  trails.plot(1.0f, 1.0f, Pixel(1, 1, 1), 0.0f, 1.0f,
              [](float, float, const Pixel &, float, float) {});
}

/**
 * @brief Death case: seeding a world trail before init_storage() must trap.
 * @details Filter surface — the 3D counterpart of the screen guard: plot()
 *          pushes into the ring buffer, which does not exist yet.
 */
inline void case_world_trails_plot_without_storage() {
  Filter::World::Trails<8> trails(4); // no init_storage() -> HS_CHECK
  trails.plot(math::Vector(0, 1, 0), Pixel(1, 1, 1), 0.0f, 1.0f,
              [](const math::Vector &, const Pixel &, float, float) {});
}

inline void case_plot_open_loop_seam() {
  constexpr int W = 32, H = 16;
  Fragment seam;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipeline;
  Plot::draw_fragments<W, H>(
      pipeline, canvas, nullptr, [](const math::Vector &, Fragment &) {},
      {.capacity = 2, .loop_seam = &seam}, [](Fragments &) {});
}

inline void case_raster_point_projection_pair_mismatch() {
  float rows[3] = {};
  float cols[3] = {};
  (void)Plot::PointProjections::paired({rows, opaque<size_t>(3)},
                                       {cols, opaque<size_t>(2)});
}

inline void case_raster_nonfinite_point() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope scope(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  Fragment point;
  point.pos = math::Vector(1, 0, 0);
  points.push_back(point);
  point.pos.x = opaque(std::numeric_limits<float>::quiet_NaN());
  points.push_back(point);
  DeathEffect effect(W, H);
  Canvas canvas(effect);
  Pipeline<W, H> pipeline;
  Plot::rasterize<W, H>(pipeline, canvas, points,
                        [](const math::Vector &, Fragment &) {});
}

inline void case_raster_edge_flags_short() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope scope(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 3);
  for (int i = 0; i < 3; ++i) {
    Fragment f;
    f.pos = math::Vector(1, 0, 0);
    points.push_back(f);
  }
  uint8_t flags[1] = {Plot::RasterOptions::EDGE_VISIBLE};
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  Pipeline<W, H> pipeline;
  Plot::rasterize<W, H>(
      pipeline, canvas, points, [](const math::Vector &, Fragment &) {},
      {.projection = Plot::RasterProjection::geodesic(flags)});
}

/**
 * @brief Death case: hoisted point projections shorter than the polyline must
 *        trap.
 * @details Plot surface — rasterize indexes point_rows/point_cols by point
 *          index, so an array sized to the EDGE count (as edge_flags is) reads
 *          one past the end on the last point.
 */
inline void case_raster_point_projections_short() {
  constexpr int W = 32, H = 16;
  configure_arenas_default();
  ScratchScope sc(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 3);
  for (int i = 0; i < 3; ++i) {
    Fragment f;
    f.pos = math::Vector(1, 0, 0);
    points.push_back(f);
  }
  float rows[3] = {0.0f, 0.0f, 0.0f};
  float cols[3] = {0.0f, 0.0f, 0.0f};
  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Plot::rasterize<W, H>(
      pipe, c, points, [](const math::Vector &, Fragment &) {},
      {.point_projections = Plot::PointProjections::paired(
           {rows, opaque<size_t>(2)}, {cols, opaque<size_t>(2)})});
}

/**
 * @brief Death case: a ring index past the last ring must trap.
 * @details SphericalFieldLayout surface — the chain walk saturates at row H-1
 *          while the offset keeps accumulating, so an out-of-range index would
 *          otherwise hand back a Ring pointing past the sample array.
 */
inline void case_spherical_field_ring_index_oob() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4);
  const auto ring = layout.ring(opaque(layout.ring_count()));
  if (ring.offset == 42)
    std::printf("x");
}

/**
 * @brief Death case: populating past the last ring must trap.
 * @details next_ring() saturates at the last ring, so an overrunning band would
 *          re-populate it and leave the caller believing it wrote fresh rings.
 */
inline void case_spherical_field_populate_ring_end_oob() {
  constexpr hs::SphericalFieldLayout<32, 16, 0> layout(4);
  static float values[layout.sample_count()];
  hs::SphericalField<float, 32, 16, 0> field(values, layout);
  field.populate(0, opaque(layout.ring_count()),
                 [](const math::Vector &v, const auto &) { return v.y; });
  if (values[0] == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: an order past the degree must trap.
 * @details reduced_legendre() has no term to recur on for |m| > l and returns
 *          0, so every sample of the mode comes back black.
 */
inline void case_spherical_harmonic_order_over_degree() {
  const float n = SHMath::normalization(opaque(2), 3);
  if (n == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a negative flat harmonic index must trap.
 * @details sqrtf of a negative argument is NaN, and the cast of a NaN to int
 *          is undefined, so the decoded level would be arbitrary.
 */
inline void case_spherical_harmonic_decode_negative_index() {
  auto [l, m] = SHMath::decode_lm(opaque(-1));
  if (l == 42 && m == 42)
    std::printf("x");
}

/**
 * @brief Death case: an infill band past the rendered domain must trap.
 * @details A south_infill wider than H puts every row at full longitude
 *          resolution, multiplying sample_count() by the spacing; the arena
 *          would then overflow at an unrelated call site.
 */
inline void case_spherical_field_infill_over_domain() {
  hs::SphericalFieldLayout<32, 16, 0> layout(4, 0, opaque(17));
  if (layout.sample_count() == 42)
    std::printf("x");
}

inline void case_latitude_geometry_degenerate_height() {
  math::LatitudeGeometry geometry(opaque(1), 0.1f, 3.0f);
  if (geometry.row_to_phi(0) == 42.0f)
    std::printf("x");
}

inline void case_latitude_geometry_reversed_span() {
  math::LatitudeGeometry geometry(16, opaque(2.0f), 1.0f);
  if (geometry.row_to_phi(0) == 42.0f)
    std::printf("x");
}

/**
 * @brief Death case: a negative feedback fade must trap in sync_hue.
 * @details Style is a public aggregate, so nothing but a slider bound keeps fade
 *          non-negative. logf of a negative yields NaN, the hue matrix carries it
 *          into every feedback pixel, and float_to_pixel16 clamps NaN to 65535 —
 *          a white buffer with no other symptom. The guard also catches a NaN
 *          fade, which compares false against zero.
 */
inline void case_feedback_negative_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(-0.25f); // logf(negative) -> HS_CHECK
  style.sync_hue();
}

/**
 * @brief Death case: an infinite feedback fade must trap in sync_hue.
 * @details A comparison against zero admits +INFINITY, whose -logf is -inf; the
 *          turn-to-cos/sin reduction of an infinite angle is NaN, which spreads
 *          through all nine hue matrix entries and whitens every feedback pixel.
 */
inline void case_feedback_infinite_fade() {
  ::Feedback::Style style = ::Feedback::Style::Smoke();
  style.fade = opaque(std::numeric_limits<float>::infinity());
  style.sync_hue();
}

/**
 * @brief Death case: Path::append_segment with zero samples must trap.
 * @details Animation surface — a zero sample count divides by zero in the
 *          t / samples term (easing(0/0) = NaN) and the loop would silently
 *          append a garbage point; the samples >= 1 guard traps the authoring
 *          error on the cold path-construction seam instead.
 */
inline void case_path_append_zero_samples() {
  Path<32> path;
  path.append_segment([](float s) { return math::Vector(s, 0.0f, 0.0f); }, 1.0f,
                      opaque(0),
                      [](float t) { return t; }); // samples < 1 -> HS_CHECK
}

/**
 * @brief Death case: a HueWobbleShade depth past the fast-trig argument range.
 */
inline void case_hue_wobble_depth_out_of_range() {
  static float phase = 0.0f;
  HueWobbleShade shade(&phase, 1.0f, opaque(1.0e6f));
  (void)shade;
}

/**
 * @brief Death case: a negative IridescentShade weight zeroes the overlay.
 */
inline void case_iridescent_weight_negative() {
  static float phase = 0.0f;
  IridescentShade shade(&phase, 3.0f, opaque(-0.5f));
  (void)shade;
}

/**
 * @brief Death case: AlphaFalloffShade requires a non-null callback.
 */
inline void case_alpha_falloff_null() {
  auto fn = opaque<AlphaFalloffShade::FalloffFunction>(nullptr);
  AlphaFalloffShade shade(fn);
  (void)shade;
}

inline void case_palette_cycler_mutated_policy() {
  static uint8_t storage[8192];
  Arena arena(storage, sizeof(storage));
  GenerativePalette first, second;
  const PaletteCycler::Entry entries[] = {first, second};
  PaletteCycler cycler;
  cycler.init(arena, entries, 2, 0, 4);
  PaletteRecipe changed;
  changed.chroma.headroom = 0.5f;
  second = GenerativePalette(changed);
  cycler.step();
  cycler.step();
}

/**
 * @brief Death case: GeneratedPaletteBank::palette rejects a mode outside the
 *        harmony enum.
 */
inline void case_generated_palette_bank_unknown_mode() {
  enum class Mode : uint8_t { TRIADIC, COMPLEMENTARY, ANALOGOUS };
  GeneratedPaletteBank bank;
  (void)bank.palette(opaque(static_cast<Mode>(3)));
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null palette source.
 */
inline void case_noise_hue_palette_direct_null_source() {
  static int8_t noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), noise_lut);
}

inline void case_noise_hue_palette_direct_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const int8_t *>(nullptr));
}

inline void case_noise_shimmer_palette_null_source() {
  static int8_t noise_lut[1];
  NoiseShimmerPalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), noise_lut);
}

inline void case_noise_shimmer_palette_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  NoiseShimmerPalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const int8_t *>(nullptr));
}

inline void case_noise_hue_palette_null_source() {
  static Pixel hue_rotation_lut[1];
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(opaque<const SolidColorPalette *>(nullptr), hue_rotation_lut,
               hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-rotation LUT.
 */
inline void case_noise_hue_palette_null_rotation_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static int8_t hue_noise_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, opaque<const Pixel *>(nullptr), hue_noise_lut);
}

/**
 * @brief Death case: NoiseHuePalette requires a non-null hue-noise LUT.
 */
inline void case_noise_hue_palette_null_noise_lut() {
  SolidColorPalette source(Color4(Pixel(255, 0, 0), 1.0f));
  static Pixel hue_rotation_lut[1];
  NoiseHuePalette<SolidColorPalette> palette;
  palette.bind(&source, hue_rotation_lut, opaque<const int8_t *>(nullptr));
}

/**
 * @brief Death case: cloning a BakedPalette from itself must trap.
 * @details Color surface — clone_from allocates fresh storage into this handle
 *          before reading @c src, so a self-clone memcpys uninitialized arena
 *          onto itself and leaves the LUT filled with garbage.
 */
inline void case_baked_palette_clone_from_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPaletteStorage lut;
  lut.bake(a, src);
  const BakedPalette &self = opaque(&lut)->view();
  lut.clone_from(self, a); // -> HS_CHECK
}

/**
 * @brief Death case: blending a BakedPalette with itself as an endpoint must
 *        trap.
 * @details Color surface — bake_blend reallocates this handle before walking
 *          the endpoints, so an endpoint that is the output reads the fresh
 *          uninitialized arena instead of the baked LUT.
 */
inline void case_baked_palette_bake_blend_self() {
  static uint8_t buf[4 * BakedPalette::required_arena_bytes()];
  Arena a(buf, sizeof(buf));
  SolidColorPalette src(Color4(Pixel(255, 0, 0), 1.0f));
  BakedPaletteStorage from, dst;
  from.bake(a, src);
  dst.bake(a, src);
  const BakedPalette &self = opaque(&dst)->view();
  dst.bake_blend(a, from, self, opaque(0.5f)); // -> HS_CHECK
}

/**
 * @brief Death case: a Gradient built from an empty stop list must trap.
 * @details The constructor requires at least one stop before reading the
 *          first stop's position and color.
 */
inline void case_gradient_no_stops() {
  Gradient grad({}); // empty stop list -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: a Gradient stop position outside [0,1] must trap.
 * @details Color surface — a stop position becomes a rounded LUT index via
 *          static_cast<int>(pos * 255 + 0.5f); sufficiently out-of-range
 *          positions can write beyond the table. The constructor traps the
 *          authoring error always-on at the cold literal-construction seam
 *          rather than corrupting memory.
 */
inline void case_gradient_stop_out_of_range() {
  Gradient grad{{0.0f, CPixel(0u, 0u, 0u)},
                {1.5f, CPixel(255u, 255u, 255u)}}; // pos > 1 -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: descending (unsorted) Gradient stops must trap.
 * @details Color surface — segments are only filled when end > start, so a
 *          transposed/unsorted pair would silently degenerate to wrong output.
 *          The constructor requires ascending positions and traps otherwise.
 */
inline void case_gradient_stops_unsorted() {
  Gradient grad{{0.6f, CPixel(0u, 0u, 0u)},
                {0.3f, CPixel(255u, 255u, 255u)}}; // descending -> HS_CHECK
  (void)grad;
}

/**
 * @brief Death case: a RandomTimer with min > max must trap.
 * @details Animation surface — reset() draws hs::rand_int(min, max + 1), a
 *          half-open range that is empty/inverted when min > max, giving an
 *          implementation-defined garbage delay. The constructor traps the
 *          inverted (or negative) range at the cold authoring seam.
 */
inline void case_random_timer_inverted_range() {
  Animation::RandomTimer timer({.min = opaque(5), .max = opaque(2)},
                               [](Canvas &) {}); // min > max -> HS_CHECK
  (void)timer;
}

/**
 * @brief Death case: calling an empty (default-constructed) Fn must trap.
 * @details Concepts surface — hs::inplace_function routes an empty-state call
 *          through ipf_empty_ops::invoke, which fail-fast traps via check_fail
 *          rather than dereferencing the empty buffer (std::function would throw
 *          bad_function_call; the engine builds without exceptions). The
 *          never-taken opaque(false) assignment keeps the optimizer from proving
 *          the function empty and folding the trap at compile time. The non-trap
 *          value semantics (copy/move/empty operator bool) are covered in-process
 *          by tests/test_concepts.h. Host/WASM only: the device Fn backend
 *          returns a zero-initialized R instead of trapping (row 9 of
 *          docs/ledgers/device_host_divergence_ledger.md).
 */
inline void case_empty_fn_call() {
  Fn<int(int), 16> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: invoking an empty FunctionRef must trap.
 * @details Concepts surface — the empty state's thunk diverges through
 *          function_ref_empty_call rather than calling through a null
 *          context. Unlike the Fn trap, this one ships to the device.
 */
inline void case_empty_function_ref_call() {
  FunctionRef<int(int)> f;
  if (opaque(false))
    f = [](int x) { return x; };
  int v = f(opaque(7)); // empty invoke -> check_fail -> trap
  if (v == 42)
    std::printf("x");
}

/**
 * @brief Death case: registering two effects under one name must trap.
 * @details Registry surface — the name keys the factory lookup and the
 *          lookup namespace, so duplicate names must be rejected.
 */
inline void case_effect_registry_duplicate_name() {
  EffectRegistration reg{};
  reg.name = "DeathDuplicate";
  reg.stable_id = "death-duplicate";
  validate_effect_registrations(std::array{reg, reg});
}

/** @brief Death case: two effects declaring the same stable ID must trap. */
inline void case_effect_registry_duplicate_stable_id() {
  EffectRegistration first{};
  first.name = "DeathStableA";
  first.stable_id = "death-stable";

  EffectRegistration second{};
  second.name = "DeathStableB";
  second.stable_id = "death-stable";
  validate_effect_registrations(std::array{first, second});
}

/** @brief Death case: a stable ID equal to another effect's name must trap. */
inline void case_effect_registry_stable_id_matches_name() {
  EffectRegistration first{};
  first.name = "DeathClassAlias";
  first.stable_id = "death-first";

  EffectRegistration second{};
  second.name = "DeathOther";
  second.stable_id = "DeathClassAlias";
  validate_effect_registrations(std::array{first, second});
}
