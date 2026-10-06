/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Gray-Scott reaction-diffusion: white-box dynamics coverage
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for GSReactionDiffusion's private fixed-point and
 *        physics internals.
 * @details The lattice is independent of <W,H>.
 */
struct GSWhiteBox {
  using GS = GSReactionDiffusion<SMALL_W, SMALL_H>;
  static constexpr int N = GS::RD_N;
  static constexpr float RENDER_MIN_WEIGHT = GS::KERNEL_MIN_TOTAL_WEIGHT;

  template <int Width, int Height> static float render_support_weight_floor() {
    using Grid = Scan::Shader::SsaaGrid<Width, Height>;
    math::PhiLUT<Height>::init();
    math::TrigLUT<Width, Height>::init();
    const float LIMIT = GS::template render_support_limit<Grid>();
    if (LIMIT <= 0.0f)
      return 0.0f;
    Grid grid;
    float minimum = 1.0f;
    for (int y = 0; y < Height; ++y) {
      grid.set_row(y);
      for (int x = 0; x < Width; x += Width / 8) {
        const math::Vector CENTER = math::pixel_to_vector<Width, Height>(x, y);
        for (int sample = 0; sample < Grid::SAMPLES; ++sample) {
          const float MAX_DISTANCE =
              LIMIT + (grid.at(x, sample) - CENTER).length();
          const float U =
              fmaxf(0.0f, 1.0f - MAX_DISTANCE * MAX_DISTANCE * GS::INV_R2);
          minimum = std::min(minimum, U * U);
        }
      }
    }
    return minimum;
  }

  static uint16_t to_q16(float v) { return GS::to_q16(v); }
  static float from_q16(uint16_t v) { return GS::from_q16(v); }
  static void fill_hot_flags(const uint16_t *b, uint8_t *hot1, uint8_t *hot2,
                             int count, uint16_t threshold) {
    GS::fill_hot_flags(b, hot1, hot2, count, threshold);
  }
  template <typename Compare>
  static void check_render_kernel(GS &gs, Compare &&compare) {
    ScratchScope guard(scratch_arena_a);
    auto lattice = gs.orient_lattice();
    using Grid = Scan::Shader::SsaaGrid<DEFAULT_W, DEFAULT_H>;
    if (!math::TrigLUT<DEFAULT_W, DEFAULT_H>::initialized)
      math::TrigLUT<DEFAULT_W, DEFAULT_H>::init();
    Grid grid;
    for (int y = 0; y < DEFAULT_H; ++y) {
      grid.set_row(y);
      for (int x = 0; x < DEFAULT_W; ++x) {
        Fragment frag;
        frag.pos = math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(x, y);
        gs.seed_face_lut(frag);
        const int center = gs.refine_render_center(frag.pos, lattice.get(),
                                                   static_cast<int>(frag.v0));
        float expected[7][4] = {}, actual[7][4] = {};
        gs.accumulate_stencil_ssaa2x2(
            grid, x, center, lattice.get(), [&](int slot, int) {
              return [&, slot](int sample, float weight) {
                expected[slot][sample] = weight;
              };
            });
        gs.accumulate_render_stencil(
            grid, x, center, lattice.get(), [&](int slot, int) {
              return [&, slot](int sample, float weight) {
                actual[slot][sample] = weight;
              };
            });
        for (int slot = 0; slot < 7; ++slot)
          for (int sample = 0; sample < 4; ++sample)
            compare(actual[slot][sample], expected[slot][sample]);
      }
    }
  }
  template <int Width, int Height, typename Compare>
  static void check_direct_draw(GSReactionDiffusion<Width, Height> &gs,
                                Compare &&compare) {
    using TestGS = GSReactionDiffusion<Width, Height>;
    ScratchScope guard(scratch_arena_a);
    auto lattice = gs.orient_lattice();
    std::vector<uint8_t> hot1(N), hot2(N);
    gs.fill_hot_flags(gs.state.B, hot1.data(), hot2.data(), N,
                      TestGS::to_q16(TestGS::B_CULL_THRESHOLD));
    std::vector<Pixel> expected(Width * Height);
    const int CLIPS[][4] = {{0, Height, 0, Width},
                            {Height / 4, Height / 2, Width / 4, Width / 2},
                            {0, Height, 0, Width / 8}};
    for (const auto &clip : CLIPS) {
      gs.set_clip(clip[0], clip[1], clip[2], clip[3]);
      gs.set_margin(3);
      {
        Canvas canvas(gs);
        auto clear = [&] {
          for (int y = 0; y < Height; ++y)
            for (int x = 0; x < Width; ++x)
              canvas(x, y) = Pixel(123, 456, 789);
        };
        clear();
        auto vertex = [&](Fragment &frag) {
          if (!hot2[static_cast<int>(frag.v0)])
            frag.v0 = -1.0f;
        };
        auto pixel = [&](Fragment &frag, const auto &grid, int x) {
          return gs.shade_pixel(static_cast<int>(frag.v0), frag.pos,
                                lattice.get(), grid, x, hot1.data());
        };
        gs.rasterize_lattice(canvas, vertex, pixel);
        for (int y = 0; y < Height; ++y)
          for (int x = 0; x < Width; ++x)
            expected[y * Width + x] = canvas(x, y);
        clear();
        gs.draw_lattice(canvas, lattice.get(), hot1.data(), hot2.data());
        for (int y = 0; y < Height; ++y)
          for (int x = 0; x < Width; ++x)
            compare(canvas(x, y), expected[y * Width + x]);
      }
      gs.advance_display();
    }
    gs.set_clip(0, Height, 0, Width);
  }
  template <typename Compare>
  static void check_staged_reseed(GS &gs, Compare &&compare) {
    const auto INITIAL_RNG = hs::random();
    gs.start_reaction();
    const std::vector<uint16_t> EXPECTED_A(gs.state.A, gs.state.A + N);
    const std::vector<uint16_t> EXPECTED_B(gs.state.B, gs.state.B + N);
    const std::vector<uint16_t> EXPECTED_PIGMENT(gs.state.pigment,
                                                 gs.state.pigment + N);
    const std::vector<Pixel> EXPECTED_PALETTES(
        gs.palettes, gs.palettes + GS::NUM_SEED_CLUSTERS * GS::PALETTE_SIZE);
    const uint32_t EXPECTED_NEXT_RANDOM = hs::random()();
    hs::random() = INITIAL_RNG;
    gs.transition.dissolve_frames = GS::DISSOLVE_FRAMES - 1;
    {
      Canvas canvas(gs);
      gs.render(canvas);
    }
    gs.advance_display();
    compare(gs.transition.next_seed, GS::NUM_SEED_CLUSTERS / 2,
            "partial reseed must place half the clusters");
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x)
        compare(gs.get_pixel(x, y), Pixel(0, 0, 0),
                "partial reseed must remain black");
    {
      Canvas canvas(gs);
      gs.render(canvas);
    }
    gs.advance_display();
    compare(gs.transition.next_seed, GS::NUM_SEED_CLUSTERS,
            "second reseed frame must complete all clusters");
    compare(gs.transition.grow_frames, 0,
            "reseed frames must not advance chemistry");
    compare(hs::random()(), EXPECTED_NEXT_RANDOM,
            "staged reseed must preserve random generator state");
    for (int i = 0; i < N; ++i) {
      compare(gs.state.A[i], EXPECTED_A[i],
              "staged reseed must preserve the A field");
      compare(gs.state.B[i], EXPECTED_B[i],
              "staged reseed must preserve the B field");
      compare(gs.state.pigment[i], EXPECTED_PIGMENT[i],
              "staged reseed must preserve pigment");
    }
    for (int i = 0; i < GS::NUM_SEED_CLUSTERS * GS::PALETTE_SIZE; ++i)
      compare(gs.palettes[i], EXPECTED_PALETTES[i],
              "staged reseed must preserve palette generation order");
    bool lit = false;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x)
        lit |= gs.get_pixel(x, y) != Pixel(0, 0, 0);
    compare(lit, true, "complete seed field was not displayed");
    gs.draw_frame();
    gs.advance_display();
    compare(gs.transition.grow_frames, 1,
            "chemistry must resume after the complete seed frame");
  }
  static auto source_palette() { return GS::make_palette(); }
  static Color4 palette_sample(const GS &gs, float t, int seed = 0) {
    return Color4(gs.palette_color(seed, t), 1.0f);
  }
  static Pixel modified_palette(const GS &gs, int seed, float t, float hue,
                                float shimmer) {
    return gs.modified_palette_color(seed, t, hue, shimmer);
  }
  static bool exact_color_path(GS &gs) {
    gs.refresh_color_palettes(true);
    return gs.color_noise_sample(0.0f).exact;
  }
  static void refresh_color_palettes(GS &gs, bool complete = false) {
    gs.refresh_color_palettes(complete);
  }
  static uint16_t color_palette_rows(const GS &gs) {
    return gs.color_palette_rows;
  }
  static bool exact_color_sample(const GS &gs, float noise) {
    return gs.color_noise_sample(noise).exact;
  }
  static Pixel staged_palette(const GS &gs, int seed, float t, float noise) {
    return gs.cached_palette_color(seed, t, noise);
  }
  static Pixel cached_palette(GS &gs, int seed, float t, float noise) {
    gs.refresh_color_palettes();
    return gs.cached_palette_color(seed, t, noise);
  }
  static Pixel wrapped_palette(const GS &gs, int seed, float t, float shift,
                               float lightness) {
    typename GS::SeedPalette source{gs.palettes + seed * GS::PALETTE_SIZE};
    NoiseHuePalette<typename GS::SeedPalette> hue(&source, gs.color_noise_lut);
    struct Shifted {
      const NoiseHuePalette<typename GS::SeedPalette> &palette;
      float shift;
      Color4 get(float value) const { return palette.get(value, shift); }
    } shifted{hue, shift};
    NoiseShimmerPalette<Shifted> shimmer(&shifted, gs.color_noise_lut);
    return shimmer.get(t, lightness).color;
  }
  static std::array<float, 4> color_params(const GS &gs) {
    return {gs.params.noise_speed, gs.params.noise_scale, gs.params.hue_shift,
            gs.params.shimmer};
  }
  static void set_color_params(GS &gs, float speed, float scale, float hue,
                               float shimmer) {
    gs.params.noise_speed = speed;
    gs.params.noise_scale = scale;
    gs.params.hue_shift = hue;
    gs.params.shimmer = shimmer;
  }
  static void advance_color_noise(GS &gs) { gs.advance_color_noise(); }
  static float projected_color_noise(const GS &gs,
                                     const math::Vector &direction) {
    return gs.sample_color_noise(ReactionGraph::CubemapLUT::project(direction));
  }
  static float noise_phase(const GS &gs) { return gs.color_noise_phase; }
  static void set_noise_phase(GS &gs, float phase) {
    gs.color_noise_phase = phase;
  }
  static float color_noise(const GS &gs, const math::Vector &direction) {
    return gs.sample_color_noise(direction);
  }
  static bool noise_bake_matches_reference(GS &gs, float scale, float phase) {
    std::array<int8_t, HueNoiseLutView::SIZE> expected;
    prepare_hue_noise_lut(expected, gs.color_noise, scale, phase);
    gs.params.noise_scale = scale;
    gs.color_noise_phase = phase;
    gs.refresh_color_noise();
    return std::equal(expected.begin(), expected.end(), gs.color_noise_lut);
  }
  static bool reaction_edited(GS &gs) { return gs.reaction_edited(); }
  static constexpr int SEEDS = GS::NUM_SEED_CLUSTERS;
  static void set_pigment(GS &gs, int node, int seed) {
    gs.state.pigment[node] = static_cast<uint16_t>(seed | (63u << 10));
  }
  static float pigment_weight(const GS &gs, int node, int seed) {
    float weights[SEEDS] = {};
    GS::add_pigment(weights, gs.state.pigment[node], 1.0f);
    return weights[seed];
  }
  template <int STEPS = 1>
  static void step_pigment(GS &gs, const float *a, const float *b) {
    std::vector<uint16_t> next(N);
    gs.template step_pigment<STEPS>(a, b, next.data());
  }
  static uint16_t *pigments(GS &gs) { return gs.state.pigment; }
  template <int STEPS = 1>
  static void pigment_reference(const GS &gs, const float *a, const float *b,
                                uint16_t *next) {
    const float DT = gs.params.dt * GS::STEP_DT_SCALE;
    const float DIFFUSION = gs.params.d_b * DT;
    for (int i = 0; i < N; ++i) {
      float weights[SEEDS] = {};
      float retained =
          fmaxf(0.0f, b[i] * (1.0f - GS::RD_K * DIFFUSION -
                              (gs.params.k + gs.params.feed) * DT) +
                          a[i] * b[i] * b[i] * DT);
      float diffusion = DIFFUSION;
      if constexpr (STEPS > 1) {
        float incoming = 0.0f;
        for (int nb : ReactionGraph::neighbors[i])
          incoming += b[nb] * DIFFUSION;
        float total = retained + incoming;
        float r = total > 0.0f ? retained / total : 0.0f;
        float power = 1.0f, scale = 1.0f;
        for (int step = 1; step < STEPS; ++step) {
          power *= r;
          scale += power;
        }
        retained *= power;
        diffusion *= scale;
      }
      GS::add_pigment(weights, gs.state.pigment[i], retained);
      for (int nb : ReactionGraph::neighbors[i])
        GS::add_pigment(weights, gs.state.pigment[nb], b[nb] * diffusion);
      next[i] = GS::pack_pigment(weights);
    }
  }
  static Pixel pigment_color(const GS &gs, int node, float t) {
    return gs.pigment_color(gs.state.pigment[node], t);
  }
  static void start_reaction(GS &gs) { gs.start_reaction(); }
  static void set_params(GS &gs, float feed, float k, float dA, float dB,
                         float dt) {
    gs.params.feed = feed;
    gs.params.k = k;
    gs.params.d_a = dA;
    gs.params.d_b = dB;
    gs.params.dt = dt;
  }
  static int dissolve_frame(const GS &gs) {
    return gs.transition.dissolve_frames;
  }
  static const uint16_t *b_field(const GS &gs) { return gs.state.B; }
  static const uint16_t *a_field(const GS &gs) { return gs.state.A; }
  static float dissolve_hash(int i, uint32_t seed) {
    return ::math::hash01(static_cast<uint32_t>(i), seed);
  }
  static void convert(GS &gs, float phase, uint32_t seed) {
    gs.transition.dissolve_seed = seed;
    gs.convert_below(phase);
  }
  static void set_node(GS &gs, int i, uint16_t a, uint16_t b) {
    gs.state.A[i] = a;
    gs.state.B[i] = b;
  }
  static constexpr float DISSOLVE_FADE_FRACTION = GS::DISSOLVE_FADE_FRACTION;
  static constexpr int DISSOLVE_FRAMES = GS::DISSOLVE_FRAMES;
  static constexpr int MIN_GROW_FRAMES = GS::MIN_GROW_FRAMES;
  static constexpr int STABLE_HOLD_FRAMES = GS::STABLE_HOLD_FRAMES;

  // refine_render_center's early-out certificates.
  static constexpr int POLE_BAND = GS::POLE_BAND;
  static constexpr float BULK_MIN_SPACING_FRAC = GS::BULK_MIN_SPACING_FRAC;
  static constexpr float BULK_CERTIFIED_D2 = GS::BULK_CERTIFIED_D2;
  static constexpr float POLE_CERTIFIED_D2 = GS::POLE_CERTIFIED_D2;

  static void step(GS &gs, const uint16_t *cA, const uint16_t *cB, uint16_t *nA,
                   uint16_t *nB) {
    // The production substep is float-resident; convert Q16 at the edges.
    std::vector<float> fA(N), fB(N), gA(N), gB(N);
    for (int i = 0; i < N; ++i) {
      fA[i] = GS::from_q16(cA[i]);
      fB[i] = GS::from_q16(cB[i]);
    }
    gs.step_physics(fA.data(), fB.data(), gA.data(), gB.data());
    for (int i = 0; i < N; ++i) {
      nA[i] = GS::to_q16(gA[i]);
      nB[i] = GS::to_q16(gB[i]);
    }
  }

  // Float-domain substep, without step()'s clamping Q16 edges.
  static void step_float(GS &gs, const float *cA, const float *cB, float *nA,
                         float *nB) {
    gs.step_physics(cA, cB, nA, nB);
  }

  static void step_float_inplace(GS &gs, float *a, float *b) {
    ScratchScope guard(scratch_arena_a);
    float *pending_a =
        scratch_arena_a.allocate_n<float>(GS::PHYSICS_HISTORY_SIZE);
    float *pending_b =
        scratch_arena_a.allocate_n<float>(GS::PHYSICS_HISTORY_SIZE);
    gs.step_physics_inplace(a, b, pending_a, pending_b);
  }

  static void validate_physics_neighbors() { GS::validate_physics_neighbors(); }
  static void validate_physics_neighbors(const ReactionGraph::NeighborRun *runs,
                                         unsigned count) {
    GS::validate_physics_neighbors(runs, count);
  }
  static constexpr int PHYSICS_NEIGHBOR_REACH = GS::PHYSICS_NEIGHBOR_REACH;
  static constexpr int STEPS_PER_FRAME = GS::STEPS_PER_FRAME;

  static void step_float_reference(const GS &gs, const float *cA,
                                   const float *cB, float *nA, float *nB) {
    const float dt = gs.params.dt * GS::STEP_DT_SCALE;
    for (int i = 0; i < N; ++i) {
      float a = cA[i];
      float b = cB[i];
      float sumA = 0.0f, sumB = 0.0f;
      for (int k = 0; k < ReactionGraph::RD_K; ++k) {
        int ni = ReactionGraph::neighbors[i][k];
        sumA += cA[ni];
        sumB += cB[ni];
      }
      float lA = sumA - ReactionGraph::RD_K * a;
      float lB = sumB - ReactionGraph::RD_K * b;
      float abb = a * b * b;
      nA[i] = hs::clamp(
          a + (gs.params.d_a * lA - abb + gs.params.feed * (1.0f - a)) * dt,
          0.0f, 1.0f);
      nB[i] = hs::clamp(
          b + (gs.params.d_b * lB + abb - (gs.params.k + gs.params.feed) * b) *
                  dt,
          0.0f, 1.0f);
    }
  }

  struct ShaderError {
    int different = 0;
    int above_rounding = 0;
    int lit = 0;
    int coverage = 0;
    int hard = 0;
    int stencil_samples = 0;
    int stencil_changes = 0;
    int center_mismatches = 0;
    int max_channel = 0;
    uint64_t total_channel = 0;
  };

  template <int W, int H>
  static ShaderError concentration_shader_error(GS &gs) {
    ScratchScope guard(scratch_arena_a);
    auto lattice = gs.orient_lattice();
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    Scan::Shader::SsaaGrid<W, H> grid;
    ShaderError error;
    for (int y = 0; y < H; ++y) {
      grid.set_row(y);
      for (int x = 0; x < W; ++x) {
        Fragment frag;
        frag.pos = math::pixel_to_vector<W, H>(x, y);
        gs.seed_face_lut(frag);
        const int SEED = static_cast<int>(frag.v0);
        const Pixel GOT =
            gs.shade_pixel(SEED, frag.pos, lattice.get(), grid, x);
        const Pixel EXPECTED =
            gs.shade_pixel_full(SEED, frag.pos, lattice.get(), grid, x);
        const int DR = std::abs(static_cast<int>(GOT.r) - EXPECTED.r);
        const int DG = std::abs(static_cast<int>(GOT.g) - EXPECTED.g);
        const int DB = std::abs(static_cast<int>(GOT.b) - EXPECTED.b);
        const bool GOT_LIT = GOT.r || GOT.g || GOT.b;
        const bool EXPECTED_LIT = EXPECTED.r || EXPECTED.g || EXPECTED.b;
        error.coverage += GOT_LIT != EXPECTED_LIT;
        error.lit += GOT_LIT || EXPECTED_LIT;
        error.max_channel =
            std::max(error.max_channel, std::max(DR, std::max(DG, DB)));
        error.total_channel += DR + DG + DB;
      }
    }
    return error;
  }

  template <int W, int H>
  static ShaderError
  shared_shader_error(GS &gs, bool fixed_stencil = false,
                      bool direct_color = false, bool aggregate_color = true,
                      bool center_pigment = false, bool full_kernel = false) {
    ScratchScope guard(scratch_arena_a);
    auto lattice = gs.orient_lattice();
    math::Vector *world_nodes = lattice.get();
    using Grid = Scan::Shader::SsaaGrid<W, H>;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    Grid grid;
    ShaderError error;
    std::vector<uint8_t> hot1(N), hot2(N);
    gs.fill_hot_flags(gs.state.B, hot1.data(), hot2.data(), N,
                      GS::to_q16(GS::B_CULL_THRESHOLD));
    for (int y = 0; y < H; ++y) {
      grid.set_row(y);
      for (int x = 0; x < W; ++x) {
        Fragment frag;
        frag.pos = math::pixel_to_vector<W, H>(x, y);
        gs.seed_face_lut(frag);
        int seed = static_cast<int>(frag.v0);
        int shared_center =
            gs.refine_render_center(frag.pos, world_nodes, seed);
        if (shared_center != gs.refine_center(frag.pos, world_nodes, seed))
          ++error.center_mismatches;
        Pixel got = full_kernel
                        ? gs.shade_pixel_full(seed, frag.pos, world_nodes, grid,
                                              x, hot1.data())
                        : gs.shade_pixel(seed, frag.pos, world_nodes, grid, x,
                                         hot1.data());
        float noise =
            gs.sample_color_noise(gs.inverse_orientation.apply(frag.pos));
        float shift = noise * gs.params.hue_shift;
        float lightness = std::max(0.0f, noise) * gs.params.shimmer;

        Pixel expected(0, 0, 0);
        if (!aggregate_color)
          for (int i = 0; i < 4; ++i) {
            math::Vector sample_v = grid.at(x, i);
            ++error.stencil_samples;
            if (gs.refine_center(sample_v, world_nodes, seed) != shared_center)
              ++error.stencil_changes;
            float b;
            if (fixed_stencil) {
              float tw = 0.0f, wb = 0.0f;
              gs.kernel_accumulate(sample_v, world_nodes, shared_center,
                                   [&](int ni, float w) {
                                     wb += gs.state.B[ni] * w;
                                     tw += w;
                                   });
              b = tw <= GS::KERNEL_MIN_TOTAL_WEIGHT ? 0.0f
                                                    : wb * (GS::Q16_INV / tw);
            } else {
              b = gs.interpolate_b(sample_v, seed, world_nodes);
            }
            if (b < GS::B_CULL_THRESHOLD)
              continue;
            float t = hs::clamp((b - GS::B_COLOR_FLOOR) * GS::B_COLOR_SCALE,
                                0.0f, 1.0f);
            float mass = 0.0f, rgb[3] = {};
            int center = fixed_stencil
                             ? shared_center
                             : gs.refine_center(sample_v, world_nodes, seed);
            gs.kernel_accumulate(
                sample_v, world_nodes, center, [&](int ni, float w) {
                  float weight = gs.state.B[ni] * w;
                  Pixel color =
                      GS::mix_pigment(gs.state.pigment[ni], [&](int id) {
                        return direct_color
                                   ? gs.modified_palette_color(id, t, shift,
                                                               lightness)
                                   : gs.cached_palette_color(id, t, noise);
                      });
                  mass += weight;
                  rgb[0] += color.r * weight;
                  rgb[1] += color.g * weight;
                  rgb[2] += color.b * weight;
                });
            if (mass > 0.0f)
              expected += Pixel(round_linear_channel(rgb[0] / (mass * 4)),
                                round_linear_channel(rgb[1] / (mass * 4)),
                                round_linear_channel(rgb[2] / (mass * 4)));
          }
        if (aggregate_color) {
          double palette_mass[GS::NUM_SEED_CLUSTERS] = {};
          double total_mass = 0.0, concentration = 0.0;
          int covered = 0;
          for (int i = 0; i < 4; ++i) {
            const math::Vector SAMPLE = grid.at(x, i);
            int center = fixed_stencil
                             ? shared_center
                             : gs.refine_center(SAMPLE, world_nodes, seed);
            ++error.stencil_samples;
            error.stencil_changes += center != shared_center;
            float tw = 0.0f, wb = 0.0f;
            gs.kernel_accumulate(
                SAMPLE, world_nodes, center, [&](int node, float w) {
                  const float MASS = gs.state.B[node] * w;
                  tw += w;
                  wb += MASS;
                  const uint16_t PIGMENT = gs.state.pigment[node];
                  const double FIRST_FRACTION =
                      ((PIGMENT >> 10) * 65535u + 31u) / 63u / 65535.0;
                  palette_mass[PIGMENT & 31] += MASS * FIRST_FRACTION;
                  palette_mass[(PIGMENT >> 5) & 31] +=
                      MASS * (1.0 - FIRST_FRACTION);
                  total_mass += MASS;
                });
            float b = tw <= GS::KERNEL_MIN_TOTAL_WEIGHT
                          ? 0.0f
                          : wb * (GS::Q16_INV / tw);
            if (b >= GS::B_CULL_THRESHOLD) {
              ++covered;
              concentration += hs::clamp(
                  (b - GS::B_COLOR_FLOOR) * GS::B_COLOR_SCALE, 0.0f, 1.0f);
            }
          }
          if (center_pigment) {
            std::fill(std::begin(palette_mass), std::end(palette_mass), 0.0);
            const uint16_t PIGMENT = gs.state.pigment[shared_center];
            const double FIRST_FRACTION =
                ((PIGMENT >> 10) * 65535u + 31u) / 63u / 65535.0;
            palette_mass[PIGMENT & 31] += FIRST_FRACTION;
            palette_mass[(PIGMENT >> 5) & 31] += 1.0 - FIRST_FRACTION;
            total_mass = 1.0;
          }
          if (covered && total_mass > 0.0) {
            float t = static_cast<float>(concentration / covered);
            double rgb[3] = {};
            for (int palette = 0; palette < GS::NUM_SEED_CLUSTERS; ++palette) {
              Pixel color =
                  direct_color
                      ? gs.modified_palette_color(palette, t, shift, lightness)
                      : gs.cached_palette_color(palette, t, noise);
              const double WEIGHT =
                  palette_mass[palette] * covered / (4 * total_mass);
              rgb[0] += color.r * WEIGHT;
              rgb[1] += color.g * WEIGHT;
              rgb[2] += color.b * WEIGHT;
            }
            expected = Pixel(round_linear_channel(rgb[0]),
                             round_linear_channel(rgb[1]),
                             round_linear_channel(rgb[2]));
          }
        }
        int dr = std::abs(static_cast<int>(got.r) - expected.r);
        int dg = std::abs(static_cast<int>(got.g) - expected.g);
        int db = std::abs(static_cast<int>(got.b) - expected.b);
        int pixel_max = std::max(dr, std::max(dg, db));
        bool got_lit = got.r || got.g || got.b;
        bool expected_lit = expected.r || expected.g || expected.b;
        if (got != expected)
          ++error.different;
        if (pixel_max > 2)
          ++error.above_rounding;
        if (got_lit || expected_lit)
          ++error.lit;
        if (got_lit != expected_lit)
          ++error.coverage;
        if (pixel_max > 4096)
          ++error.hard;
        error.max_channel = std::max(error.max_channel, pixel_max);
        error.total_channel += dr + dg + db;
      }
    }
    Pixel culled = gs.shade_pixel(-1, math::Vector(), world_nodes, grid, 0);
    if (culled != Pixel(0, 0, 0)) {
      ++error.different;
      ++error.above_rounding;
    }
    return error;
  }
};

/** @brief Verifies source palette recipes are opaque before the RGB bake. */
inline void test_gs_palette_is_opaque() {
  hs_test::reset_globals();
  for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed) {
    auto palette = GSWhiteBox::source_palette();
    for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
      HS_EXPECT_EQ(
          palette.get(static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1))
              .alpha,
          1.0f);
  }
}

/** @brief Verifies each reseeded reaction receives a freshly generated palette. */
inline void test_gs_reseed_generates_palette() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  std::array<Pixel, BakedPalette::LUT_SIZE> before;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i)
    before[i] = GSWhiteBox::palette_sample(gs, static_cast<float>(i) /
                                                   (BakedPalette::LUT_SIZE - 1))
                    .color;

  const size_t arena_offset = persistent_arena.get_offset();
  GSWhiteBox::start_reaction(gs);

  int changed = 0;
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const Color4 after = GSWhiteBox::palette_sample(
        gs, static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1));
    changed += after.color != before[i];
  }
  HS_EXPECT_GT(changed, 0);
  HS_EXPECT_EQ(persistent_arena.get_offset(), arena_offset);
}

/** @brief Shares cube projection without changing noise on faces, edges or corners. */
inline void test_gs_shared_cube_projection() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (float phase : {0.0f, 0.123f, 0.9999f}) {
    GSWhiteBox::set_noise_phase(gs, phase);
    GSWhiteBox::advance_color_noise(gs);
    for (int face = 0; face < HueNoiseLutView::FACE_COUNT; ++face)
      for (float u : {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f})
        for (float v : {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f}) {
          const auto DIRECTION = hue_noise_face_direction(face, u, v);
          const auto PROJECTION = ReactionGraph::CubemapLUT::project(DIRECTION);
          if (fabsf(u) < 1.0f && fabsf(v) < 1.0f) {
            HS_EXPECT_EQ(PROJECTION.face, face);
            HS_EXPECT_NEAR(PROJECTION.u, face < 2 ? -u : u, 2e-7f);
            HS_EXPECT_NEAR(PROJECTION.v, face >= 2 && face < 4 ? -v : v, 2e-7f);
          }
          HS_EXPECT_EQ(GSWhiteBox::projected_color_noise(gs, DIRECTION),
                       GSWhiteBox::color_noise(gs, DIRECTION));
        }
    for (const auto &direction : ReactionGraph::node_positions)
      HS_EXPECT_EQ(GSWhiteBox::projected_color_noise(gs, direction),
                   GSWhiteBox::color_noise(gs, direction));
  }
}

inline void test_gs_shared_noise_palette_modifiers() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  const size_t OFFSET = persistent_arena.get_offset();
  int rotated = 0;
  for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed) {
    for (float t : {0.0f, 0.37f, 1.0f})
      for (float shift : {-0.3f, 0.0f, 0.25f})
        for (float lift : {0.0f, 0.5f, 1.0f})
          HS_EXPECT_EQ(GSWhiteBox::modified_palette(gs, seed, t, shift, lift),
                       GSWhiteBox::wrapped_palette(gs, seed, t, shift, lift));
    Pixel base = GSWhiteBox::palette_sample(gs, 0.5f, seed).color;
    HS_EXPECT_EQ(GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.0f, 0.0f),
                 base);
    Pixel hue = GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.25f, 0.0f);
    Pixel expected = hue_rotate_lut_gamut(Color4(base, 1.0f), 0.25f).color;
    HS_EXPECT_EQ(hue, expected);
    rotated += hue != base;
    Pixel shimmer = GSWhiteBox::modified_palette(gs, seed, 0.5f, 0.0f, 0.5f);
    const LinRGB BEFORE = pixel_to_linrgb(base);
    const LinRGB AFTER = pixel_to_linrgb(shimmer);
    HS_EXPECT_GT(linear_rgb_to_oklab(AFTER.r, AFTER.g, AFTER.b).L,
                 linear_rgb_to_oklab(BEFORE.r, BEFORE.g, BEFORE.b).L);
  }
  HS_EXPECT_EQ(rotated, GSWhiteBox::SEEDS);

  GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  const float STILL = GSWhiteBox::color_noise(gs, math::X_AXIS);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_EQ(GSWhiteBox::noise_phase(gs), 0.0f);
  HS_EXPECT_EQ(GSWhiteBox::color_noise(gs, math::X_AXIS), STILL);
  GSWhiteBox::set_color_params(gs, 0.001f, 2.0f, 0.35f, 0.25f);
  for (int i = 0; i < 40; ++i)
    GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_GT(fabsf(GSWhiteBox::color_noise(gs, math::X_AXIS) - STILL), 1e-3f);
  const float MOVING = GSWhiteBox::color_noise(gs, math::X_AXIS);
  GSWhiteBox::set_color_params(gs, 0.0f, 4.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_GT(fabsf(GSWhiteBox::color_noise(gs, math::X_AXIS) - MOVING),
               1e-3f);
  GSWhiteBox::set_noise_phase(gs, 0.0005f);
  GSWhiteBox::set_color_params(gs, -0.001f, 4.0f, 0.35f, 0.25f);
  GSWhiteBox::advance_color_noise(gs);
  HS_EXPECT_NEAR(GSWhiteBox::noise_phase(gs), 0.9995f, 1e-6f);
  HS_EXPECT(!GSWhiteBox::reaction_edited(gs),
            "color controls must not restart the reaction");
  HS_EXPECT_EQ(persistent_arena.get_offset(), OFFSET);
}

inline void test_gs_noise_bake_matches_reference() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (float scale : {1.0f / 64.0f, 4.941984f, 8.0f}) {
    for (float phase : {0.0f, 1.0f / 65536.0f, 1.0f / 6.0f, 0.5f, 0.999f}) {
      HS_EXPECT_TRUE(
          GSWhiteBox::noise_bake_matches_reference(gs, scale, phase));
      HS_EXPECT_TRUE(
          GSWhiteBox::noise_bake_matches_reference(gs, scale, phase));
    }
  }
}

inline void test_gs_seed_palettes_and_pigment_transport() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int seed = 1; seed < GSWhiteBox::SEEDS; ++seed) {
    bool distinct = false;
    for (int j = 0; j < 9; ++j)
      distinct |= GSWhiteBox::palette_sample(gs, j / 8.0f, seed).color !=
                  GSWhiteBox::palette_sample(gs, j / 8.0f, seed - 1).color;
    HS_EXPECT(distinct, "each seed must receive an independent palette");
  }

  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f);
  constexpr int NODE = GSWhiteBox::N / 2;
  const int NEIGHBOR = ReactionGraph::neighbors[NODE][0];
  b[NODE] = 0.5f;
  GSWhiteBox::set_pigment(gs, NODE, 7);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NEIGHBOR, 7), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 7), 1.0f, 1e-6f);

  b[NEIGHBOR] = 0.5f;
  GSWhiteBox::set_pigment(gs, NEIGHBOR, 12);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  float mix = GSWhiteBox::pigment_weight(gs, NODE, 12);
  HS_EXPECT_GT(mix, 0.0f);
  HS_EXPECT_LT(mix, 1.0f);
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 7) + mix, 1.0f, 1e-6f);
  Pixel first = GSWhiteBox::palette_sample(gs, 0.5f, 7).color;
  Pixel second = GSWhiteBox::palette_sample(gs, 0.5f, 12).color;
  Pixel blended = GSWhiteBox::pigment_color(gs, NODE, 0.5f);
  HS_EXPECT_NEAR(blended.r, first.r * (1.0f - mix) + second.r * mix, 2.0f);
  HS_EXPECT_NEAR(blended.g, first.g * (1.0f - mix) + second.g * mix, 2.0f);
  HS_EXPECT_NEAR(blended.b, first.b * (1.0f - mix) + second.b * mix, 2.0f);

  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.0f, 2.5f);
  GSWhiteBox::step_pigment(gs, a.data(), b.data());
  HS_EXPECT_NEAR(GSWhiteBox::pigment_weight(gs, NODE, 12), mix, 1e-6f);

  const size_t OFFSET = persistent_arena.get_offset();
  GSWhiteBox::start_reaction(gs);
  HS_EXPECT_EQ(persistent_arena.get_offset(), OFFSET);
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    if (GSWhiteBox::b_field(gs)[i] != 0)
      continue;
    HS_EXPECT_EQ(GSWhiteBox::pigment_weight(gs, i, 0), 1.0f);
  }
}

inline void test_gs_dot_kernel_matches_squared_distance() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  const auto compare = [](float actual, float expected) {
    HS_EXPECT_LE(fabsf(actual - expected), 0.001f);
  };
  GSWhiteBox::check_render_kernel(gs, compare);
  for (int i = 0; i < 40; ++i) {
    gs.draw_frame();
    gs.advance_display();
  }
  GSWhiteBox::check_render_kernel(gs, compare);
}

inline void test_gs_direct_draw_matches_grid() {
  hs_test::reset_globals();
  GSReactionDiffusion<DEFAULT_W, DEFAULT_H> gs;
  gs.init();
  const auto compare = [](Pixel actual, Pixel expected) {
    HS_EXPECT(actual == expected, "GS direct traversal changes clipped pixels");
  };
  GSWhiteBox::check_direct_draw(gs, compare);
  for (int frame = 0; frame < 40; ++frame) {
    gs.draw_frame();
    gs.advance_display();
  }
  GSWhiteBox::check_direct_draw(gs, compare);
}

inline void test_gs_sparse_pigment_matches_dense() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
  std::vector<uint16_t> expected(GSWhiteBox::N);
  const float PARAMS[][5] = {{0.04f, 0.06f, 0.02f, 0.01f, 2.5f},
                             {0.1f, 0.1f, 0.05f, 0.05f, 3.0f},
                             {0.04f, 0.06f, 0.02f, 0.0f, 0.1f},
                             {0.04f, 0.06f, 0.02f, -0.0f, 0.1f}};
  for (const auto &params : PARAMS) {
    GSWhiteBox::set_params(gs, params[0], params[1], params[2], params[3],
                           params[4]);
    for (int mode = 0; mode < 5; ++mode) {
      for (int i = 0; i < GSWhiteBox::N; ++i) {
        uint32_t h = (static_cast<uint32_t>(i) + 173u) * 2654435761u;
        a[i] = static_cast<float>(h & 65535u) / 65535.0f;
        b[i] = mode == 0 ? 0.0f
                         : static_cast<float>((h >> 16) & 65535u) / 65535.0f;
        int first = mode == 1 ? 17 : static_cast<int>(h % GSWhiteBox::SEEDS);
        int second =
            mode == 1 ? 0 : static_cast<int>((h >> 8) % GSWhiteBox::SEEDS);
        int mix = mode <= 2 ? 63 : mode == 3 ? 32 : static_cast<int>(h & 63u);
        GSWhiteBox::pigments(gs)[i] =
            static_cast<uint16_t>(first | (second << 5) | (mix << 10));
      }
      for (int step = 0; step < GSWhiteBox::STEPS_PER_FRAME; ++step) {
        GSWhiteBox::pigment_reference(gs, a.data(), b.data(), expected.data());
        GSWhiteBox::step_pigment(gs, a.data(), b.data());
        for (int i = 0; i < GSWhiteBox::N; ++i)
          HS_EXPECT_EQ(GSWhiteBox::pigments(gs)[i], expected[i]);
        GSWhiteBox::pigment_reference<GSWhiteBox::STEPS_PER_FRAME>(
            gs, a.data(), b.data(), expected.data());
        GSWhiteBox::step_pigment<GSWhiteBox::STEPS_PER_FRAME>(gs, a.data(),
                                                              b.data());
        for (int i = 0; i < GSWhiteBox::N; ++i)
          HS_EXPECT_EQ(GSWhiteBox::pigments(gs)[i], expected[i]);
      }
    }
  }
}

/**
 * @brief Verifies the Q16 fixed-point round-trip and the +0.5 rounding/clamp
 *        boundaries.
 * @details to_q16(from_q16(v)) must be the identity over every representable
 *          value, and to_q16 must clamp out-of-range floats and round to nearest
 *          (so 1.0 tops out at 65535 with no overflow).
 */
inline void test_gs_q16_roundtrip() {
  HS_EXPECT_EQ(GSWhiteBox::to_q16(0.0f), (uint16_t)0);
  HS_EXPECT_EQ(GSWhiteBox::to_q16(1.0f), (uint16_t)65535);
  HS_EXPECT_EQ(GSWhiteBox::to_q16(2.0f), (uint16_t)65535); // clamp high
  HS_EXPECT_EQ(GSWhiteBox::to_q16(-0.5f), (uint16_t)0);    // clamp low
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(0), 0.0f, 1e-9f);
  HS_EXPECT_NEAR(GSWhiteBox::from_q16(65535), 1.0f, 1e-9f);
  int bad = 0;
  for (int v = 0; v <= 65535; ++v)
    if (GSWhiteBox::to_q16(GSWhiteBox::from_q16((uint16_t)v)) != (uint16_t)v)
      ++bad;
  HS_EXPECT_EQ(bad, 0);
}

/** @brief Compares both cull rings to the original directed adjacency table. */
inline void test_gs_hot_flags_match_directed_graph() {
  constexpr int N = GSWhiteBox::N;
  constexpr uint8_t GUARD = 0xa5;
  std::vector<uint16_t> b(N);
  std::vector<uint8_t> hot1(N + 2), hot2(N + 2), ref1(N), ref2(N);
  const uint16_t THRESHOLDS[] = {0, 1, GSWhiteBox::to_q16(0.1f), 65535};
  for (uint16_t threshold : THRESHOLDS) {
    const uint16_t BELOW = threshold == 0 ? 0 : threshold - 1;
    const uint16_t ABOVE = threshold == 65535 ? 65535 : threshold + 1;
    for (int pattern = 0; pattern < 9; ++pattern) {
      for (int i = 0; i < N; ++i) {
        const uint32_t HASH = static_cast<uint32_t>(i) * 2654435761u;
        switch (pattern) {
        case 0:
          b[i] = 0;
          break;
        case 1:
          b[i] = 65535;
          break;
        case 2:
          b[i] = BELOW;
          break;
        case 3:
          b[i] = threshold;
          break;
        case 4:
          b[i] = HASH % 257 == 0 ? ABOVE : BELOW;
          break;
        case 5:
          b[i] = HASH % 7 != 0 ? ABOVE : BELOW;
          break;
        case 6:
          b[i] = HASH % 3 == 0 ? BELOW : HASH % 3 == 1 ? threshold : ABOVE;
          break;
        case 7:
          b[i] = i == 0 ? ABOVE : BELOW;
          break;
        case 8:
          b[i] = i == N - 1 ? ABOVE : BELOW;
          break;
        }
      }
      for (int i = 0; i < N; ++i) {
        ref1[i] = b[i] >= threshold;
        for (int k = 0; k < ReactionGraph::RD_K; ++k)
          ref1[i] |= b[ReactionGraph::neighbors[i][k]] >= threshold;
      }
      for (int i = 0; i < N; ++i) {
        ref2[i] = ref1[i];
        for (int k = 0; k < ReactionGraph::RD_K; ++k)
          ref2[i] |= ref1[ReactionGraph::neighbors[i][k]];
      }
      std::fill(hot1.begin(), hot1.end(), GUARD);
      std::fill(hot2.begin(), hot2.end(), GUARD);
      GSWhiteBox::fill_hot_flags(b.data(), hot1.data() + 1, hot2.data() + 1, N,
                                 threshold);
      HS_EXPECT(std::equal(ref1.begin(), ref1.end(), hot1.begin() + 1),
                "GS one-ring flags match directed graph");
      HS_EXPECT(std::equal(ref2.begin(), ref2.end(), hot2.begin() + 1),
                "GS two-ring flags match directed graph");
      HS_EXPECT_EQ(hot1.front(), GUARD);
      HS_EXPECT_EQ(hot1.back(), GUARD);
      HS_EXPECT_EQ(hot2.front(), GUARD);
      HS_EXPECT_EQ(hot2.back(), GUARD);
    }
  }
}

/**
 * @brief Re-measures the lattice and requires refine_render_center's early-out
 *        certificates to bound the true minimum node spacing.
 * @details The early-out returns the seed unwalked whenever the query lies
 *          within sqrt(safe_d2) of it, which is the nearest node only while
 *          safe_d2 stays at or under a quarter of the seed's true
 *          nearest-neighbor distance squared. The certificates are hand-measured
 *          properties of the generated lattice; this recomputes them from
 *          node().
 *
 *          The per-node minimum is exact: node()'s y is strictly decreasing in
 *          the index and chord distance is at least |dy|, so the outward scan
 *          stops once |dy| reaches the running best.
 */
inline void test_gs_render_certificates_bound_lattice() {
  const int n = GSWhiteBox::N;
  std::vector<math::Vector> lattice(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i)
    lattice[static_cast<size_t>(i)] = ReactionGraph::node(i);

  const float bulk_floor =
      GSWhiteBox::BULK_MIN_SPACING_FRAC * ReactionGraph::D_AVG;
  float bulk_min_d2 = 4.0f; // antipodal chord squared: the loosest possible
  float pole_min_d2 = 4.0f;
  int tight_band = 0; // band width the sub-bulk-spacing nodes actually need
  for (int i = 0; i < n; ++i) {
    const math::Vector &p = lattice[static_cast<size_t>(i)];
    float best = 4.0f;
    for (int step = -1; step <= 1; step += 2) {
      for (int j = i + step; j >= 0 && j < n; j += step) {
        const math::Vector &q = lattice[static_cast<size_t>(j)];
        float dy = p.y - q.y;
        if (dy * dy >= best)
          break;
        float dx = p.x - q.x, dz = p.z - q.z;
        float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best)
          best = d2;
      }
    }
    if (i < GSWhiteBox::POLE_BAND || i >= n - GSWhiteBox::POLE_BAND) {
      pole_min_d2 = std::min(pole_min_d2, best);
    } else {
      bulk_min_d2 = std::min(bulk_min_d2, best);
    }
    if (best < bulk_floor * bulk_floor)
      tight_band = std::max(tight_band, std::min(i, n - 1 - i) + 1);
  }

  const float bulk_frac = std::sqrt(bulk_min_d2) / ReactionGraph::D_AVG;
  std::printf("  [info] GS lattice: bulk min spacing %.6f*D_AVG, bulk d2/4 "
              "%.9f, pole d2/4 %.9f, tight band %d\n",
              bulk_frac, 0.25f * bulk_min_d2, 0.25f * pole_min_d2, tight_band);

  HS_EXPECT_LE(GSWhiteBox::BULK_MIN_SPACING_FRAC, bulk_frac);
  HS_EXPECT_LE(GSWhiteBox::BULK_CERTIFIED_D2, 0.25f * bulk_min_d2);
  HS_EXPECT_LE(GSWhiteBox::POLE_CERTIFIED_D2, 0.25f * pole_min_d2);
  // Every node packed tighter than the bulk minimum must fall inside the band
  // that gets the smaller certificate.
  HS_EXPECT_LE(tight_band, GSWhiteBox::POLE_BAND);
}

/**
 * @brief Bounds the shared stencil against independent SSAA refinement.
 * @details Compares production-resolution pixels after 4, 16, and 40 rendered
 * frames. Near-tie Voronoi samples may change, but coverage and high-amplitude
 * errors remain confined to under 3% of lit pixels. The spatial-change count
 * excludes differences of at most two RGB16 quantization units.
 */
inline void test_gs_shared_stencil_error_is_bounded() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr int probe_frames[] = {4, 16, 40};
  int next_probe = 0;
  for (int frame = 1; frame <= probe_frames[2]; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != probe_frames[next_probe])
      continue;
    auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, false, false, true, true, true);
    std::printf("GS shared stencil frame=%d: different=%d above_rounding=%d "
                "lit=%d coverage=%d "
                "hard=%d center_changes=%d/%d center_mismatches=%d max=%d "
                "total=%llu\n",
                frame, error.different, error.above_rounding, error.lit,
                error.coverage, error.hard, error.stencil_changes,
                error.stencil_samples, error.center_mismatches,
                error.max_channel,
                static_cast<unsigned long long>(error.total_channel));
    HS_EXPECT_GT(error.lit, DEFAULT_W * DEFAULT_H / 100);
    HS_EXPECT_EQ(error.center_mismatches, 0);
    HS_EXPECT(error.above_rounding * 2 <= error.lit,
              "shared stencil changed over half of lit pixels beyond rounding");
    HS_EXPECT(error.coverage * 50 <= error.lit,
              "shared stencil moved over 2% of lit coverage");
    HS_EXPECT(error.hard * 100 <= error.lit * 3,
              "shared stencil put over 3% of lit pixels past 6.25% error");
    HS_EXPECT_LE(error.max_channel, 23000);
    HS_EXPECT(error.total_channel <=
                  static_cast<uint64_t>(error.lit) * 3u * 384u,
              "shared stencil mean channel error exceeds 3/512 full scale");
    if (next_probe < 2)
      ++next_probe;
  }
}

inline void test_gs_partial_color_palette_rows() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr uint16_t ALL_ROWS = (1u << 15) - 1;
  HS_EXPECT_EQ(GSWhiteBox::color_palette_rows(gs), ALL_ROWS);
  const size_t OFFSET = persistent_arena.get_offset();
  for (int generation = 0; generation < 3; ++generation) {
    if (generation == 1)
      GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, -0.2f, 0.3f);
    else
      GSWhiteBox::start_reaction(gs);
    GSWhiteBox::refresh_color_palettes(gs);
    HS_EXPECT_EQ(GSWhiteBox::color_palette_rows(gs),
                 static_cast<uint16_t>(31u << 5));
    HS_EXPECT_TRUE(GSWhiteBox::exact_color_sample(gs, -1.0f));
    HS_EXPECT_TRUE(GSWhiteBox::exact_color_sample(gs, 1.0f));
    float hue = generation == 0 ? 0.2f : -0.2f;
    float shimmer = generation == 0 ? 0.4f : 0.3f;
    for (float noise : {-1.0f + 12.0f / 14.0f, 0.0f}) {
      HS_EXPECT_FALSE(GSWhiteBox::exact_color_sample(gs, noise));
      for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed)
        for (float value : {0.0f, 1.0f})
          HS_EXPECT_EQ(
              GSWhiteBox::staged_palette(gs, seed, value, noise),
              GSWhiteBox::modified_palette(gs, seed, value, noise * hue,
                                           std::max(noise, 0.0f) * shimmer));
    }
    for (int frame = 0; frame < 7; ++frame) {
      for (float noise : {-1.0f, -0.8f, -0.5f, 0.0f, 0.5f, 0.8f, 1.0f}) {
        if (!GSWhiteBox::exact_color_sample(gs, noise))
          continue;
        for (int seed = 0; seed < GSWhiteBox::SEEDS; ++seed)
          for (float value : {0.0f, 0.37f, 1.0f})
            HS_EXPECT_EQ(
                GSWhiteBox::staged_palette(gs, seed, value, noise),
                GSWhiteBox::modified_palette(gs, seed, value, noise * hue,
                                             std::max(noise, 0.0f) * shimmer));
      }
      GSWhiteBox::refresh_color_palettes(gs);
    }
    HS_EXPECT_EQ(GSWhiteBox::color_palette_rows(gs), ALL_ROWS);
    HS_EXPECT_EQ(persistent_arena.get_offset(), OFFSET);
  }
}

/** @brief Checks certified support against the fixed-profile SSAA offsets. */
inline void test_gs_support_certificate_matches_display_geometry() {
  const auto check = [] {
    HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<96, 20>()),
                 GSWhiteBox::RENDER_MIN_WEIGHT);
    HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<288, 16>()),
                 GSWhiteBox::RENDER_MIN_WEIGHT);
    HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<288, 144>()),
                 GSWhiteBox::RENDER_MIN_WEIGHT);
    HS_EXPECT_GT((GSWhiteBox::render_support_weight_floor<4096, 16>()),
                 GSWhiteBox::RENDER_MIN_WEIGHT);
  };
  check();
}

/** @brief Preserves four-sample coverage while bounding reciprocal estimation error. */
inline void test_gs_signed_coverage_and_concentration() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int frame = 1; frame <= 650; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 40 && frame != 200 && frame != 400 &&
        frame != 450 && frame != 470 && frame != 600 && frame != 650)
      continue;
    const auto ERROR =
        GSWhiteBox::concentration_shader_error<DEFAULT_W, DEFAULT_H>(gs);
    HS_EXPECT_EQ(ERROR.coverage, 0);
    HS_EXPECT_LE(ERROR.max_channel, 6000);
    HS_EXPECT_LE(ERROR.total_channel,
                 static_cast<uint64_t>(ERROR.lit) * 3u * 125u);
  }
  for (int pattern = 0; pattern < 8; ++pattern) {
    for (int node = 0; node < GSWhiteBox::N; ++node) {
      const uint32_t HASH = static_cast<uint32_t>(node) * 2654435761u;
      const uint16_t B =
          pattern < 4
              ? static_cast<uint16_t>(6552 + pattern)
              : static_cast<uint16_t>(6552 + ((HASH >> (pattern - 4)) & 3));
      GSWhiteBox::set_node(gs, node, 65535, B);
    }
    const auto ERROR =
        GSWhiteBox::concentration_shader_error<DEFAULT_W, DEFAULT_H>(gs);
    HS_EXPECT_EQ(ERROR.coverage, 0);
    const auto FALLBACK = GSWhiteBox::concentration_shader_error<16, 8>(gs);
    HS_EXPECT_EQ(FALLBACK.max_channel, 0);
    HS_EXPECT_EQ(FALLBACK.total_channel, 0u);
  }
}

/** @brief Bounds nearest-pigment color fidelity against procedural and aggregate references. */
inline void test_gs_nearest_pigment_shader_fidelity() {
  // Float kernel rearrangement can move one pixel across the render floor.
  constexpr int MAX_COVERAGE_DIFFERENCES = 1;
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  for (int frame = 1; frame <= 600; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 40 && frame != 200 && frame != 400 &&
        frame != 500 && frame != 600)
      continue;
    const auto ERROR = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, true, true, false);
    HS_EXPECT_LE(ERROR.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(ERROR.max_channel, 22000);
    HS_EXPECT_LE(ERROR.total_channel,
                 static_cast<uint64_t>(ERROR.lit) * 3u * 850u);
    HS_EXPECT_LE(ERROR.hard * 10, ERROR.lit);
    const auto ROUNDING = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, true, false, true, true, true);
    HS_EXPECT_LE(ROUNDING.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(ROUNDING.max_channel, 2048);
    HS_EXPECT_LE(ROUNDING.total_channel,
                 static_cast<uint64_t>(ROUNDING.lit) * 6u);
    const auto AGGREGATE =
        GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(gs, true, false,
                                                              true, false);
    std::printf(
        "GS nearest pigment frame=%d procedural_mae=%.2f procedural_max=%d "
        "aggregate_mae=%.2f aggregate_max=%d rounding_max=%d\n",
        frame, static_cast<double>(ERROR.total_channel) / (ERROR.lit * 3),
        ERROR.max_channel,
        static_cast<double>(AGGREGATE.total_channel) / (AGGREGATE.lit * 3),
        AGGREGATE.max_channel, ROUNDING.max_channel);
    HS_EXPECT_LE(AGGREGATE.coverage, MAX_COVERAGE_DIFFERENCES);
    HS_EXPECT_LE(AGGREGATE.max_channel, 22000);
    HS_EXPECT_LE(AGGREGATE.total_channel,
                 static_cast<uint64_t>(AGGREGATE.lit) * 3u * 400u);
    if (frame != 200 && frame != 400 && frame != 600)
      continue;
    const auto CONTROLS = GSWhiteBox::color_params(gs);
    for (float hue : {-4.0f, -0.25f, 0.0f, 0.25f, 4.0f}) {
      for (float shimmer : {0.0f, 0.4f, 1.0f}) {
        GSWhiteBox::set_color_params(gs, CONTROLS[0], CONTROLS[1], hue,
                                     shimmer);
        GSWhiteBox::cached_palette(gs, 0, 0.0f, 0.0f);
        const auto EXTREME =
            GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(gs, true,
                                                                  true, false);
        if (EXTREME.max_channel > 30000)
          std::printf(
              "GS nearest pigment extreme frame=%d hue=%.2f shimmer=%.2f max=%d\n",
              frame, static_cast<double>(hue), static_cast<double>(shimmer),
              EXTREME.max_channel);
        HS_EXPECT_LE(EXTREME.coverage, MAX_COVERAGE_DIFFERENCES);
        HS_EXPECT_LE(EXTREME.max_channel, 35000);
        HS_EXPECT_LE(EXTREME.total_channel,
                     static_cast<uint64_t>(EXTREME.lit) * 3u * 1000u);
      }
    }
    GSWhiteBox::set_color_params(gs, CONTROLS[0], CONTROLS[1], CONTROLS[2],
                                 CONTROLS[3]);
  }
  for (float hue : {-4.0f, -0.3f, 0.0f, 0.3f, 4.0f}) {
    const float SHIMMER = hue == 0.0f ? 0.0f : 0.4f;
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, hue, SHIMMER);
    for (float value : {0.0f, 0.37f, 1.0f})
      for (float noise : {-1.0f, -0.25f, 0.0f, 0.8f, 1.0f})
        HS_EXPECT_EQ(
            GSWhiteBox::cached_palette(gs, 0, value, noise),
            GSWhiteBox::modified_palette(gs, 0, value, hue * noise,
                                         std::max(noise, 0.0f) * SHIMMER));
  }
  for (float shimmer : {0.0f, 0.4f, 0.4001f, 1.0f}) {
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, 0.2f, shimmer);
    HS_EXPECT_EQ(GSWhiteBox::exact_color_path(gs), shimmer > 0.4f);
    if (shimmer > 0.4f)
      for (float noise : {-1.0f, 0.0f, 0.31f, 1.0f})
        HS_EXPECT_EQ(
            GSWhiteBox::cached_palette(gs, 0, 0.37f, noise),
            GSWhiteBox::modified_palette(gs, 0, 0.37f, 0.2f * noise,
                                         std::max(noise, 0.0f) * shimmer));
  }
  for (float hue : {-0.2501f, -0.25f, 0.25f, 0.2501f}) {
    GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, hue, 0.4f);
    HS_EXPECT_EQ(GSWhiteBox::exact_color_path(gs), fabsf(hue) > 0.25f);
  }
}

/** @brief Checks nearest-pigment shading against independent scalar SSAA samples. */
inline void test_gs_nearest_pigment_shader_matches_scalar_reference() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::set_color_params(gs, 0.0f, 2.0f, 0.0f, 0.0f);
  for (int frame = 1; frame <= 40; ++frame) {
    gs.draw_frame();
    gs.advance_display();
    if (frame != 4 && frame != 16 && frame != 40)
      continue;
    auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
        gs, true, false, true, true, true);
    HS_EXPECT_EQ(error.coverage, 0);
    HS_EXPECT_LE(error.max_channel, 2048);
    HS_EXPECT_LE(error.total_channel, static_cast<uint64_t>(error.lit) * 3u);
  }
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    uint16_t b = i < GSWhiteBox::N / 2 ? 6554 : 32768;
    GSWhiteBox::set_node(gs, i, 65535, b);
  }
  auto error = GSWhiteBox::shared_shader_error<DEFAULT_W, DEFAULT_H>(
      gs, true, false, true, true, true);
  HS_EXPECT_EQ(error.coverage, 0);
  HS_EXPECT_LE(error.max_channel, 2048);
}

inline void test_gs_dissolve_frontier_fades_before_clear() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  constexpr float phase = 0.5f;
  constexpr uint32_t seed = 0x7b31d2a5u;
  int cleared = -1, fading = -1;
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    float h = GSWhiteBox::dissolve_hash(i, seed);
    if (cleared < 0 && h < phase)
      cleared = i;
    float fade = (h - phase) / GSWhiteBox::DISSOLVE_FADE_FRACTION;
    if (fading < 0 && fade > 0.25f && fade < 0.75f)
      fading = i;
  }
  HS_EXPECT(cleared >= 0 && fading >= 0,
            "dissolve probe did not find clear and fade nodes");
  if (cleared < 0 || fading < 0)
    return;
  GSWhiteBox::set_node(gs, cleared, 0, 65535);
  GSWhiteBox::set_node(gs, fading, 0, 65535);
  GSWhiteBox::convert(gs, phase, seed);
  HS_EXPECT_EQ(GSWhiteBox::a_field(gs)[cleared], (uint16_t)65535);
  HS_EXPECT_EQ(GSWhiteBox::b_field(gs)[cleared], (uint16_t)0);
  HS_EXPECT(GSWhiteBox::a_field(gs)[fading] > 0 &&
                GSWhiteBox::a_field(gs)[fading] < 65535,
            "dissolve fade did not ease A toward rest");
  HS_EXPECT(GSWhiteBox::b_field(gs)[fading] > 0 &&
                GSWhiteBox::b_field(gs)[fading] < 65535,
            "dissolve fade cleared B without an intermediate value");
}

/**
 * @brief The two-frame staged reseed matches start_reaction() fields, palettes
 *        and RNG position.
 */
inline void test_gs_staged_reseed_matches_synchronous() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  GSWhiteBox::check_staged_reseed(
      gs, [](const auto &actual, const auto &expected, const char *message) {
        HS_EXPECT(actual == expected, message);
      });
}

/**
 * @brief Verifies dissolve coverage falls below one eighth and reseeds.
 * @details On the last dissolve frame, coverage is below one eighth of the
 *          grown pattern. The closing frame restores the fresh seed pattern.
 */
inline void test_gs_dissolve_clears_and_reseeds() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();

  auto covered = [&]() {
    const uint16_t *b = GSWhiteBox::b_field(gs);
    int n = 0;
    for (int i = 0; i < GSWhiteBox::N; ++i)
      if (b[i] >= GSWhiteBox::to_q16(0.1f))
        ++n;
    return n;
  };

  // Run until the stalled field trips a dissolve; the default reaction settles
  // well inside this budget.
  const int budget =
      GSWhiteBox::MIN_GROW_FRAMES + GSWhiteBox::STABLE_HOLD_FRAMES + 400;
  int f = 0;
  for (; f < budget && GSWhiteBox::dissolve_frame(gs) < 0; ++f) {
    gs.draw_frame();
    gs.advance_display();
  }
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "stalled field did not start a dissolve");
  const int grown = covered();
  HS_EXPECT(grown > GSWhiteBox::N / 20,
            "default reaction grew no pattern to dissolve");

  // Sample the last frame before the dissolve window closes.
  for (int i = 0; i < GSWhiteBox::DISSOLVE_FRAMES - 1; ++i) {
    gs.draw_frame();
    gs.advance_display();
  }
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "dissolve ended before its window closed");
  HS_EXPECT(covered() < grown / 8,
            "dissolve left the sphere covered: cleared nodes healed");

  // The closing frame converts the remainder and reseeds.
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  HS_EXPECT(covered() > 0, "reseed left no nucleation sites");
}

/**
 * @brief Verifies editing the reaction constants starts a dissolve.
 * @details Feed/Kill/dA/dB define the reaction, so an edit clears and restarts
 *          it; Speed only sets the integration rate and must not.
 */
inline void test_gs_reaction_edit_starts_dissolve() {
  hs_test::reset_globals();
  GSWhiteBox::GS gs;
  gs.init();
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);

  for (float speed : {0.1f, 1.0f}) {
    GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, speed);
    gs.draw_frame();
    gs.advance_display();
    HS_EXPECT_EQ(GSWhiteBox::dissolve_frame(gs), -1);
  }

  GSWhiteBox::set_params(gs, 0.03f, 0.06f, 0.02f, 0.01f, 1.0f); // Feed moved
  gs.draw_frame();
  gs.advance_display();
  HS_EXPECT(GSWhiteBox::dissolve_frame(gs) >= 0,
            "a Feed edit did not start a dissolve");
}

/**
 * @brief Verifies the homogeneous rest state (A=1, B=0) is a Gray-Scott fixed
 *        point: one substep leaves it exactly unchanged.
 * @details With B=0 the reaction term A·B² and both Laplacians vanish and the
 *          feed term feed·(1-A) is zero at A=1, so the state must not move. A
 *          stray nonzero term in the update assembly perturbs it here.
 */
inline void test_gs_rest_state_is_fixed_point() {
  std::vector<uint16_t> cA(GSWhiteBox::N, 65535), cB(GSWhiteBox::N, 0),
      nA(GSWhiteBox::N), nB(GSWhiteBox::N);
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f); // defaults
  GSWhiteBox::step(gs, cA.data(), cB.data(), nA.data(), nB.data());
  int moved = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (nA[i] != 65535 || nB[i] != 0)
      ++moved;
  HS_EXPECT_EQ(moved, 0);
}

/**
 * @brief Pins the optimized float substep to the scalar equation bit-for-bit.
 */
inline void test_gs_substep_matches_scalar_reference() {
  std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
  std::vector<float> gotA(GSWhiteBox::N), gotB(GSWhiteBox::N);
  std::vector<float> refA(GSWhiteBox::N), refB(GSWhiteBox::N);
  for (int i = 0; i < GSWhiteBox::N; ++i) {
    a[i] = static_cast<float>((i * 73) & 65535) * (1.0f / 65535.0f);
    b[i] = static_cast<float>((i * 151 + 17) & 65535) * (1.0f / 65535.0f);
  }
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.037f, 0.061f, 0.019f, 0.013f, 2.3f);
  GSWhiteBox::step_float(gs, a.data(), b.data(), gotA.data(), gotB.data());
  GSWhiteBox::step_float_reference(gs, a.data(), b.data(), refA.data(),
                                   refB.data());
  HS_EXPECT(
      std::memcmp(gotA.data(), refA.data(), gotA.size() * sizeof(float)) == 0,
      "optimized A substep differs from scalar reference");
  HS_EXPECT(
      std::memcmp(gotB.data(), refB.data(), gotB.size() * sizeof(float)) == 0,
      "optimized B substep differs from scalar reference");
}

/** @brief Pins delayed writes to full float Jacobi across a complete frame. */
inline void test_gs_inplace_frame_matches_jacobi() {
  GSWhiteBox::validate_physics_neighbors();
  for (int i = 0; i < GSWhiteBox::N; ++i)
    for (int neighbor : ReactionGraph::neighbors[i])
      HS_EXPECT_LE(i - neighbor, GSWhiteBox::PHYSICS_NEIGHBOR_REACH);

  const float PARAMS[][5] = {{0.04f, 0.06f, 0.02f, 0.01f, 2.5f},
                             {0.1f, 0.1f, 0.05f, 0.05f, 3.0f},
                             {0.0f, 0.0f, 0.0f, 0.0f, 0.1f}};
  GSWhiteBox::GS gs;
  for (const auto &params : PARAMS) {
    GSWhiteBox::set_params(gs, params[0], params[1], params[2], params[3],
                           params[4]);
    std::vector<float> a(GSWhiteBox::N), b(GSWhiteBox::N);
    for (int i = 0; i < GSWhiteBox::N; ++i) {
      a[i] = static_cast<float>((i * 73) & 65535) * (1.0f / 65535.0f);
      b[i] = static_cast<float>((i * 151 + 17) & 65535) * (1.0f / 65535.0f);
    }
    std::vector<float> ref_a = a, ref_b = b;
    std::vector<float> next_a(GSWhiteBox::N), next_b(GSWhiteBox::N);
    for (int step = 0; step < GSWhiteBox::STEPS_PER_FRAME; ++step) {
      GSWhiteBox::step_float(gs, ref_a.data(), ref_b.data(), next_a.data(),
                             next_b.data());
      ref_a.swap(next_a);
      ref_b.swap(next_b);
      GSWhiteBox::step_float_inplace(gs, a.data(), b.data());
      HS_EXPECT(std::memcmp(a.data(), ref_a.data(), a.size() * sizeof(float)) ==
                    0,
                "delayed A writes differ from full Jacobi");
      HS_EXPECT(std::memcmp(b.data(), ref_b.data(), b.size() * sizeof(float)) ==
                    0,
                "delayed B writes differ from full Jacobi");
    }
  }
}

/**
 * @brief Verifies one substep has the right reaction/diffusion signs and that
 *        the Q16 clamp is actually applied.
 * @details Seed a single saturated-B nucleus on the otherwise-rest field. After
 *          one step at effective dt = 2.5 * (10 / 6) = 4.17: A at the seed is consumed (the
 *          1 - dt update underflows and must clamp to 0, not wrap), and B
 *          diffuses into at least one neighbor that started empty.
 */
inline void test_gs_substep_signs_and_clamp() {
  std::vector<uint16_t> cA(GSWhiteBox::N, 65535), cB(GSWhiteBox::N, 0),
      nA(GSWhiteBox::N), nB(GSWhiteBox::N);
  const int seed = 4000; // an interior lattice node with a full neighbor ring
  cB[seed] = 65535;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.02f, 0.01f, 2.5f);
  GSWhiteBox::step(gs, cA.data(), cB.data(), nA.data(), nB.data());

  // a + (dA·0 - 1 + feed·0)·dt = 1 - 4.17 < 0 → clamps to 0 (not an unclamped
  // negative-float-to-uint16 wrap).
  HS_EXPECT_EQ((int)nA[seed], 0);
  // B diffuses outward: at least one initially-empty neighbor is now lit.
  int spread = 0;
  for (int k = 0; k < ReactionGraph::RD_K; ++k) {
    int nb = ReactionGraph::neighbors[seed][k];
    if (nb >= 0 && nB[nb] > 0)
      ++spread;
  }
  HS_EXPECT_GT(spread, 0);
}

/**
 * @brief Verifies the explicit-Euler integrator does not diverge over many
 *        substeps at a high-diffusion stable setting.
 * @details The stability product dt·D·|λ|max = 5·0.03·12 = 1.8 ≤ 2, so the
 *          scheme must stay bounded: after 256 steps from seeded nuclei almost
 *          no node may sit at the upper rail. Whether B persists or decays is
 *          regime-dependent and not asserted.
 */
inline void test_gs_evolution_stays_bounded() {
  std::vector<uint16_t> a(GSWhiteBox::N, 65535), b(GSWhiteBox::N, 0),
      sa(GSWhiteBox::N), sb(GSWhiteBox::N);
  for (int s : {500, 2500, 4500, 6500})
    b[s] = 65535;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.04f, 0.06f, 0.03f, 0.03f, 3.0f);
  uint16_t *cA = a.data(), *cB = b.data(), *nA = sa.data(), *nB = sb.data();
  for (int s = 0; s < 256; ++s) {
    GSWhiteBox::step(gs, cA, cB, nA, nB);
    std::swap(cA, nA);
    std::swap(cB, nB);
  }
  // No blow-up: a stable run leaves at most a handful of nodes at the upper
  // rail; an unstable oscillation would clamp a large fraction there.
  int saturated = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (cB[i] == 65535)
      ++saturated;
  HS_EXPECT_LT(saturated, GSWhiteBox::N / 20);
}

/**
 * @brief Verifies the substep stays finite and inside [0, 1] at the joint
 *        feed/k corner of the slider box, past the Euler stability bound.
 * @details At top Speed/diffusion, 5 * 0.05 * 12 = 3 exceeds the Euler
 * bound of 2; the reaction term also exceeds it. Float fields remain finite
 * and in [0, 1] after 256 clamped substeps, before Q16 output conversion.
 */
inline void test_gs_reaction_corner_stays_bounded() {
  std::vector<float> a(GSWhiteBox::N, 1.0f), b(GSWhiteBox::N, 0.0f),
      sa(GSWhiteBox::N), sb(GSWhiteBox::N);
  for (int s : {500, 2500, 4500, 6500})
    b[s] = 1.0f;
  GSWhiteBox::GS gs;
  GSWhiteBox::set_params(gs, 0.1f, 0.1f, 0.05f, 0.05f,
                         3.0f); // joint slider-box corner
  float *cA = a.data(), *cB = b.data(), *nA = sa.data(), *nB = sb.data();
  for (int s = 0; s < 256; ++s) {
    GSWhiteBox::step_float(gs, cA, cB, nA, nB);
    std::swap(cA, nA);
    std::swap(cB, nB);
  }
  int escaped = 0;
  for (int i = 0; i < GSWhiteBox::N; ++i)
    if (!std::isfinite(cA[i]) || !std::isfinite(cB[i]) || cA[i] < 0.0f ||
        cA[i] > 1.0f || cB[i] < 0.0f || cB[i] > 1.0f)
      ++escaped;
  HS_EXPECT_EQ(escaped, 0);
}
