/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file
 * @brief Gray-Scott white-box accessors and reference shaders. */

/**
 * @brief White-box accessor for GSReactionDiffusion's private fixed-point and
 *        physics internals.
 * @details The lattice is independent of <W,H>.
 */
struct GSWhiteBox {
  using GS = GSReactionDiffusion<SMALL_W, SMALL_H>;
  static constexpr int N = GS::RD_N;
  static constexpr float RENDER_MIN_WEIGHT = GS::KERNEL_MIN_TOTAL_WEIGHT;

  /** @brief Packs byte-reference flags into a guarded node bit set. */
  static std::vector<uint32_t>
  pack_hot_flags(const std::vector<uint8_t> &flags) {
    std::vector<uint32_t> bits(GS::HOT_WORDS);
    for (size_t i = 0; i < flags.size(); ++i)
      if (flags[i])
        bits[GS::HOT_GUARD_WORDS + (i >> 5)] |= uint32_t{1} << (i & 31);
    return bits;
  }

  template <int Width, int Height> static float render_support_weight_floor() {
    using Grid = Scan::Shader::SsaaGrid<Width, Height>;
    math::PhiLUT<Height>::init();
    math::TrigLUT<Width, Height>::init();
    const float limit = GS::template render_support_limit<Grid>();
    if (limit <= 0.0f)
      return 0.0f;
    Grid grid;
    float minimum = 1.0f;
    for (int y = 0; y < Height; ++y) {
      grid.set_row(y);
      for (int x = 0; x < Width; x += Width / 8) {
        const math::Vector center = math::pixel_to_vector<Width, Height>(x, y);
        for (int sample = 0; sample < Grid::SAMPLES; ++sample) {
          const float max_distance =
              limit + (grid.at(x, sample) - center).length();
          const float u =
              fmaxf(0.0f, 1.0f - max_distance * max_distance * GS::INV_R2);
          minimum = std::min(minimum, u * u);
        }
      }
    }
    return minimum;
  }

  struct ColorNoiseSample {
    int row;
    uint16_t weight;
    float value;
    bool exact;
  };

  struct ColorValueSample {
    int column;
    uint16_t weight;
    float value;
  };

  static ColorNoiseSample color_noise_sample(const GS &gs, float noise) {
    float position =
        hs::clamp((noise + 1.0f) * (0.5f * (GS::COLOR_NOISE_STEPS - 1)), 0.0f,
                  static_cast<float>(GS::COLOR_NOISE_STEPS - 1));
    const auto selection = gs.shimmer_cache.selection(position);
    return {selection.row,
            static_cast<uint16_t>((position - selection.row) * 65535.0f), noise,
            selection.exact};
  }

  static ColorValueSample color_value_sample(float t) {
    float position = t * (GS::COLOR_VALUE_STEPS - 1);
    int column =
        std::min(static_cast<int>(position), GS::COLOR_VALUE_STEPS - 2);
    return {column, lut_index_weight(position, column), t};
  }

  /** @brief Bilinear read of the cached palette rows, or the exact color. */
  static Pixel cached_palette_color(const GS &gs, int seed, float t,
                                    float noise) {
    const ColorValueSample value = color_value_sample(t);
    const ColorNoiseSample noise_value = color_noise_sample(gs, noise);
    if (noise_value.exact)
      return gs.modified_palette_color(
          seed, value.value, noise_value.value * gs.params.hue_shift,
          fmaxf(noise_value.value, 0.0f) * gs.params.shimmer);
    const typename GS::FloatRgb *row =
        gs.shimmer_cache.row(seed, noise_value.row) + value.column;
    Pixel first = row[0].pixel().lerp16(row[1].pixel(), value.weight);
    Pixel second = row[GS::COLOR_VALUE_STEPS].pixel().lerp16(
        row[GS::COLOR_VALUE_STEPS + 1].pixel(), value.weight);
    return first.lerp16(second, noise_value.weight);
  }

  using Pigment = typename GS::Pigment;

  /** @brief A pigment mass as a 16-bit first-palette lerp weight. */
  static uint16_t first_weight16(int mass) {
    return static_cast<uint16_t>(
        (static_cast<unsigned>(mass) * 65535u + Pigment::FULL_MASS / 2) /
        Pigment::FULL_MASS);
  }

  /** @brief Blends a packed pigment's two palettes by its first-mass weight. */
  template <typename Sample>
  static Pixel mix_pigment(uint16_t pigment, Sample &&sample) {
    const auto pigment_value = Pigment::unpack(pigment);
    Pixel first = sample(pigment_value.first);
    if (pigment_value.mass == Pigment::FULL_MASS)
      return first;
    Pixel second = sample(pigment_value.second);
    return second.lerp16(first, first_weight16(pigment_value.mass));
  }

  /** @brief Jacobi reference substep for the in-place physics. */
  static void step_physics(GS &gs, const float *c_a, const float *c_b,
                           float *n_a, float *n_b) {
    gs.step_physics_nodes(c_a, c_b, [&](int i, float a, float b) {
      n_a[i] = a;
      n_b[i] = b;
    });
  }

  /**
   * @brief Kernel-weighted B at one point; 0 when no node is within the
   * support radius.
   */
  static float interpolate_b(const GS &gs, const math::Vector &p, int seed,
                             const math::Vector *nodes) {
    float tw = 0, wb = 0;
    GS::refine_and_accumulate(p, nodes, seed, [&](int i, float w) {
      wb += GS::from_q16(gs.state.B[i]) * w;
      tw += w;
    });
    if (tw <= GS::KERNEL_MIN_TOTAL_WEIGHT)
      return 0.0f;
    return wb / tw;
  }

  static uint16_t to_q16(float v) { return GS::to_q16(v); }
  static float from_q16(uint16_t v) { return GS::from_q16(v); }
  /** @brief shade_field() against the lattice rotated by the current
   *         orientation's matrix instead of its quaternion. */
  static void shade_matrix_oriented_field(GS &gs) {
    {
      ScratchScope guard(scratch_arena_a);
      Canvas canvas(gs);
      std::vector<uint8_t> hot1(N), hot2(N);
      gs.fill_hot_flags(gs.state.B, hot1.data(), hot2.data(), N,
                        GS::to_q16(GS::B_CULL_THRESHOLD));
      auto lattice = gs.orient_lattice(); // latches the seed lookup rotation
      const math::RotationMatrix m(gs.orientation.get());
      std::vector<math::Vector> world(N);
      for (int i = 0; i < N; ++i)
        world[i] = m.apply(gs.nodes[i]);
      const auto bits1 = pack_hot_flags(hot1), bits2 = pack_hot_flags(hot2);
      gs.draw_lattice(canvas, world.data(), gs.hot_bits(bits1.data()),
                      gs.hot_bits(bits2.data()));
    }
    gs.advance_display();
  }
  /** @brief Persistent per-node state and lifecycle of a GS instance. */
  struct Snapshot {
    std::vector<uint16_t> a, b, pigment;
    GS::Transition transition;

    bool operator==(const Snapshot &other) const {
      return a == other.a && b == other.b && pigment == other.pigment &&
             transition.next_seed == other.transition.next_seed &&
             transition.grow_frames == other.transition.grow_frames &&
             transition.stable_frames == other.transition.stable_frames &&
             transition.dissolve_frames == other.transition.dissolve_frames &&
             transition.dissolve_seed == other.transition.dissolve_seed;
    }
  };
  static Snapshot snapshot(const GS &gs) {
    return {std::vector<uint16_t>(gs.state.A, gs.state.A + N),
            std::vector<uint16_t>(gs.state.B, gs.state.B + N),
            std::vector<uint16_t>(gs.state.pigment, gs.state.pigment + N),
            gs.transition};
  }
  /** @brief Starts a dissolve that ends `frames` frames from now. */
  static void dissolve_within(GS &gs, int frames) {
    gs.transition.dissolve_frames = GS::DISSOLVE_FRAMES - frames;
    gs.transition.dissolve_seed = 0x5eed1234u;
  }
  /**
   * @brief Shades a Q16 B field over the whole frame through the lattice
   *        basis, converted by a dissolve at `phase` (when non-negative), with
   *        colour noise blended at `noise_weight` from the previous frame's
   *        table and the given pigment, then restores the state.
   */
  static std::vector<Pixel>
  shade_slab_field(GS &gs, const std::vector<uint16_t> &b, float phase,
                   float noise_weight, const std::vector<uint16_t> &pigment) {
    const Snapshot saved = snapshot(gs);
    std::copy(b.begin(), b.end(), gs.state.B);
    std::copy(pigment.begin(), pigment.end(), gs.state.pigment);
    if (phase >= 0.0f)
      gs.convert_below(phase);
    gs.slab_noise = {gs.noise_history ? gs.previous_noise_lut : nullptr,
                     noise_weight};
    shade_field(gs, true);
    gs.slab_noise = {};
    std::copy(saved.a.begin(), saved.a.end(), gs.state.A);
    std::copy(saved.b.begin(), saved.b.end(), gs.state.B);
    std::copy(saved.pigment.begin(), saved.pigment.end(), gs.state.pigment);
    std::vector<Pixel> frame;
    capture_frame<SMALL_W, SMALL_H>(gs, frame);
    return frame;
  }
  static bool noise_history(const GS &gs) { return gs.noise_history; }
  static float noise_phase(const GS &gs) { return gs.color_noise_phase; }
  static void set_noise_scale(GS &gs, float scale) {
    gs.params.noise_scale = scale;
  }
  static std::vector<int8_t> previous_noise(const GS &gs) {
    return {gs.previous_noise_lut,
            gs.previous_noise_lut + Pattern::SphereLutView::SIZE};
  }
  static std::vector<int8_t> current_noise(const GS &gs) {
    return {gs.color_noise_lut,
            gs.color_noise_lut + Pattern::SphereLutView::SIZE};
  }
  /** @brief The table the noise bake gives at `phase`, from the effect's
   *         noise source and scale. */
  static std::vector<int8_t> bake_noise(const GS &gs, float phase) {
    std::vector<int8_t> table(Pattern::SphereLutView::SIZE);
    Pattern::prepare_sphere_lut<Pattern::SphereNoise>(
        std::span<int8_t, Pattern::SphereLutView::SIZE>(
            table.data(), Pattern::SphereLutView::SIZE),
        Pattern::SphereNoise::prepare({.scale = gs.params.noise_scale}, {phase},
                                      gs.color_noise));
    return table;
  }
  /** @brief Replaces the view orientation's frame endpoints. */
  static void set_orientation(GS &gs, const math::Quaternion &previous,
                              const math::Quaternion &current) {
    gs.orientation.set(previous);
    gs.orientation.push(current);
  }
  /**
   * @brief Shades the current field over the whole frame, through the frame's
   *        lattice basis or against the lattice oriented at the current
   *        endpoint.
   */
  static void shade_field(GS &gs, bool columns) {
    {
      ScratchScope guard(scratch_arena_a);
      Canvas canvas(gs);
      if (columns) {
        std::vector<uint32_t> bits1(GS::HOT_WORDS), bits2(GS::HOT_WORDS);
        fill_hot_bits(gs.state.B, GS::to_q16(GS::B_CULL_THRESHOLD),
                      bits1.data(), bits2.data());
        gs.draw_lattice_columns(canvas, *gs.prepare_lattice_basis(),
                                GS::hot_bits(bits1.data()),
                                GS::hot_bits(bits2.data()));
      } else {
        std::vector<uint8_t> hot1(N), hot2(N);
        gs.fill_hot_flags(gs.state.B, hot1.data(), hot2.data(), N,
                          GS::to_q16(GS::B_CULL_THRESHOLD));
        auto lattice = gs.orient_lattice();
        const auto bits1 = pack_hot_flags(hot1), bits2 = pack_hot_flags(hot2);
        gs.draw_lattice(canvas, lattice.get(), gs.hot_bits(bits1.data()),
                        gs.hot_bits(bits2.data()));
      }
    }
    gs.advance_display();
  }
  static void fill_hot_flags(const uint16_t *b, uint8_t *hot1, uint8_t *hot2,
                             int count, uint16_t threshold) {
    GS::fill_hot_flags(b, hot1, hot2, count, threshold);
  }
  static constexpr int HOT_WORDS = GS::HOT_WORDS;
  static constexpr int HOT_GUARD_WORDS = GS::HOT_GUARD_WORDS;
  /**
   * @brief Builds the one- and two-ring bit sets the slab shader reads.
   * @param b Per-node B, Q16.
   * @param threshold Q16 render floor.
   * @param hot1 Receives the one-ring set, HOT_WORDS words.
   * @param hot2 Receives the two-ring set, HOT_WORDS words.
   */
  static void fill_hot_bits(const uint16_t *b, uint16_t threshold,
                            uint32_t *hot1, uint32_t *hot2) {
    std::vector<uint32_t> marked(HOT_WORDS);
    GS::mark_hot_bits(b, marked.data(), threshold);
    GS::dilate_hot_bits(marked.data(), hot1);
    GS::dilate_hot_bits(hot1, hot2);
  }
  static bool hot_bit(const uint32_t *bits, int node) {
    return GS::hot_bits(bits)[node];
  }
  /** @brief Returns a stencil callback storing slot and sample weights. */
  template <size_t Samples>
  static auto recorder(float (&weights)[GS::RD_K + 1][Samples]) {
    return [&](int slot, int) {
      return [&, slot](int sample, float weight) {
        HS_CONTEXT("slot, sample", slot, sample);
        HS_EXPECT_GE(slot, 0);
        HS_EXPECT_LT(slot, GS::RD_K + 1);
        HS_EXPECT_GE(sample, 0);
        HS_EXPECT_LT(sample, static_cast<int>(Samples));
        if (slot >= 0 && slot < GS::RD_K + 1 && sample >= 0 &&
            sample < static_cast<int>(Samples))
          weights[slot][sample] = weight;
      };
    };
  }

  /**
   * @brief Compares the render stencil against the SSAA reference stencil.
   * @return Count of non-zero reference weights compared.
   */
  static int check_render_kernel(GS &gs) {
    int weighted = 0;
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
        float expected[GS::RD_K + 1][Grid::SAMPLES] = {};
        float actual[GS::RD_K + 1][Grid::SAMPLES] = {};
        gs.accumulate_stencil_ssaa2x2(grid, x, center, lattice.get(),
                                      recorder(expected));
        gs.accumulate_render_stencil(grid, x, center, lattice.get(),
                                     recorder(actual));
        for (int slot = 0; slot < GS::RD_K + 1; ++slot)
          for (int sample = 0; sample < Grid::SAMPLES; ++sample) {
            weighted += expected[slot][sample] != 0.0f;
            HS_CONTEXT("px", x, y);
            HS_CONTEXT("slot, sample", slot, sample);
            HS_EXPECT_LE(fabsf(actual[slot][sample] - expected[slot][sample]),
                         0.001f);
          }
      }
    }
    return weighted;
  }
  /**
   * @brief Compares direct lattice drawing against the rasterizer per clip.
   * @return Fewest lit reference pixels drawn in any clip.
   */
  template <int Width, int Height>
  static size_t check_direct_draw(GSReactionDiffusion<Width, Height> &gs) {
    constexpr Pixel SENTINEL(123, 456, 789);
    size_t fewest_lit = SIZE_MAX;
    using TestGS = GSReactionDiffusion<Width, Height>;
    ScratchScope guard(scratch_arena_a);
    auto lattice = gs.orient_lattice();
    std::vector<uint8_t> hot1(N), hot2(N);
    gs.fill_hot_flags(gs.state.B, hot1.data(), hot2.data(), N,
                      TestGS::to_q16(TestGS::B_CULL_THRESHOLD));
    std::vector<Pixel> expected(Width * Height);
    const int clips[][4] = {{0, Height, 0, Width},
                            {Height / 4, Height / 2, Width / 4, Width / 2},
                            {0, Height, 0, Width / 8}};
    for (const auto &clip : clips) {
      gs.set_clip(clip[0], clip[1], clip[2], clip[3]);
      gs.set_margin(3);
      {
        Canvas canvas(gs);
        auto clear = [&] {
          for (int y = 0; y < Height; ++y)
            for (int x = 0; x < Width; ++x)
              canvas(x, y) = SENTINEL;
        };
        clear();
        auto pixel = [&](Fragment &frag, const auto &grid, int x) {
          const int seed = static_cast<int>(frag.v0);
          return hot2[seed] ? gs.shade_pixel(seed, frag.pos, lattice.get(),
                                             grid, x, hot1.data())
                            : Pixel(0, 0, 0);
        };
        gs.rasterize_lattice(canvas, pixel);
        size_t lit = 0;
        for (int y = 0; y < Height; ++y)
          for (int x = 0; x < Width; ++x) {
            expected[y * Width + x] = canvas(x, y);
            lit += canvas(x, y) != SENTINEL && !is_black(canvas(x, y));
          }
        fewest_lit = std::min(fewest_lit, lit);
        clear();
        const auto bits1 = pack_hot_flags(hot1), bits2 = pack_hot_flags(hot2);
        gs.draw_lattice(canvas, lattice.get(), gs.hot_bits(bits1.data()),
                        gs.hot_bits(bits2.data()));
        for (int y = 0; y < Height; ++y)
          for (int x = 0; x < Width; ++x) {
            HS_CONTEXT("px", x, y);
            HS_EXPECT_EQ(canvas(x, y), expected[y * Width + x]);
          }
      }
      gs.advance_display();
    }
    gs.set_clip(0, Height, 0, Width);
    return fewest_lit;
  }
  static void check_staged_reseed(GS &gs) {
    const auto initial_rng = hs::random();
    gs.start_reaction();
    const std::vector<uint16_t> expected_a(gs.state.A, gs.state.A + N);
    const std::vector<uint16_t> expected_b(gs.state.B, gs.state.B + N);
    const std::vector<uint16_t> expected_pigment(gs.state.pigment,
                                                 gs.state.pigment + N);
    const std::vector<Pixel> expected_palettes(
        gs.palettes, gs.palettes + GS::NUM_SEED_CLUSTERS * GS::PALETTE_SIZE);
    const uint32_t expected_next_random = hs::random()();
    hs::random() = initial_rng;
    gs.transition.dissolve_frames = GS::DISSOLVE_FRAMES - 1;
    {
      Canvas canvas(gs);
      gs.render(canvas);
    }
    gs.advance_display();
    HS_EXPECT_GT(gs.transition.next_seed, 0);
    HS_EXPECT_LT(gs.transition.next_seed, GS::NUM_SEED_CLUSTERS);
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x)
        HS_EXPECT_EQ(gs.get_pixel(x, y), Pixel(0, 0, 0));
    {
      Canvas canvas(gs);
      gs.render(canvas);
    }
    gs.advance_display();
    HS_EXPECT_EQ(gs.transition.next_seed, GS::NUM_SEED_CLUSTERS);
    HS_EXPECT_EQ(gs.transition.grow_frames, 0);
    HS_EXPECT_EQ(hs::random()(), expected_next_random);
    for (int i = 0; i < N; ++i) {
      HS_EXPECT_EQ(gs.state.A[i], expected_a[i]);
      HS_EXPECT_EQ(gs.state.B[i], expected_b[i]);
      HS_EXPECT_EQ(gs.state.pigment[i], expected_pigment[i]);
    }
    for (int i = 0; i < GS::NUM_SEED_CLUSTERS * GS::PALETTE_SIZE; ++i)
      HS_EXPECT_EQ(gs.palettes[i], expected_palettes[i]);
    const bool lit = count_lit_region<SMALL_W, SMALL_H>(gs) > 0;
    HS_EXPECT_EQ(lit, true);
    gs.draw_frame();
    gs.advance_display();
    HS_EXPECT_EQ(gs.transition.grow_frames, 1);
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
    return color_noise_sample(gs, 0.0f).exact;
  }
  static void refresh_color_palettes(GS &gs, bool complete = false) {
    gs.refresh_color_palettes(complete);
  }
  static uint16_t color_palette_rows(const GS &gs) {
    return gs.shimmer_cache.rows;
  }
  static bool exact_color_sample(const GS &gs, float noise) {
    return color_noise_sample(gs, noise).exact;
  }
  /** @brief `shade_pigment` of a solid `seed` pigment at unit scale. */
  static Pixel staged_palette(const GS &gs, int seed, float t, float noise) {
    return gs.shade_pigment(Pigment::pack(seed, seed, Pigment::FULL_MASS), t,
                            1.0f, noise);
  }
  static Pixel cached_palette(GS &gs, int seed, float t, float noise) {
    gs.refresh_color_palettes();
    return staged_palette(gs, seed, t, noise);
  }
  static Pixel wrapped_palette(const GS &gs, int seed, float t, float shift,
                               float lightness) {
    typename GS::SeedPalette source{gs.palettes + seed * GS::PALETTE_SIZE};
    const auto shifted =
        ColorModulation::HueRotation::apply(source.get(t), shift, {1.0f});
    return ColorModulation::PositiveLightnessLift::apply(shifted, lightness,
                                                         {1.0f})
        .color;
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
  static void set_noise_phase(GS &gs, float phase) {
    gs.color_noise_phase = phase;
  }
  static float color_noise(const GS &gs, const math::Vector &direction) {
    return gs.sample_color_noise(direction);
  }
  static bool noise_bake_matches_reference(GS &gs, float scale, float phase) {
    gs.params.noise_scale = scale;
    gs.color_noise_phase = phase;
    const auto expected = bake_noise(gs, phase);
    gs.refresh_color_noise();
    return std::equal(expected.begin(), expected.end(), gs.color_noise_lut);
  }
  static bool reaction_edited(GS &gs) { return gs.reaction_edited(); }
  static constexpr int SEEDS = GS::NUM_SEED_CLUSTERS;
  static constexpr int COLOR_NOISE_STEPS = GS::COLOR_NOISE_STEPS;
  static constexpr int COLOR_VALUE_STEPS = GS::COLOR_VALUE_STEPS;
  static constexpr float CACHED_HUE_LIMIT = GS::CACHED_HUE_LIMIT;
  static constexpr float CACHED_SHIMMER_LIMIT = GS::CACHED_SHIMMER_LIMIT;
  static void set_pigment(GS &gs, int node, int seed) {
    gs.state.pigment[node] = Pigment::pack(seed, 0, Pigment::FULL_MASS);
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
    std::copy(next.begin(), next.end(), gs.state.pigment);
  }
  static uint16_t *pigments(GS &gs) { return gs.state.pigment; }
  template <int STEPS = 1>
  static void pigment_reference(const GS &gs, const float *a, const float *b,
                                uint16_t *next) {
    const float dt = gs.params.dt * GS::STEP_DT_SCALE;
    const float diffusion_value = gs.params.d_b * dt;
    for (int i = 0; i < N; ++i) {
      float weights[SEEDS] = {};
      float retained =
          fmaxf(0.0f, b[i] * (1.0f - GS::RD_K * diffusion_value -
                              (gs.params.k + gs.params.feed) * dt) +
                          a[i] * b[i] * b[i] * dt);
      float diffusion = diffusion_value;
      if constexpr (STEPS > 1) {
        float incoming = 0.0f;
        for (int nb : ReactionGraph::neighbors[i])
          incoming += b[nb] * diffusion_value;
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
    return mix_pigment(gs.state.pigment[node],
                       [&](int seed) { return gs.palette_color(seed, t); });
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
  static void advance_transition(GS &gs, float mean_db) {
    gs.advance_transition(mean_db);
  }
  static bool reseed_pending(const GS &gs) {
    return gs.transition.reseed_pending();
  }
  static int grow_frames(const GS &gs) { return gs.transition.grow_frames; }
  static void finish_reseed(GS &gs) { gs.finish_reseed(); }
  /** @brief The stabilization floor on mean |dB| at Speed `dt`. */
  static float stable_floor(float dt) {
    return GS::MEAN_DB_STABLE * (dt * (1.0f / GS::DEFAULT_DT)) *
           (static_cast<float>(GS::EVOLUTION_STEPS_PER_FRAME) /
            GS::BASELINE_STEPS_PER_FRAME);
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
  static constexpr float B_CULL_THRESHOLD = GS::B_CULL_THRESHOLD;
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
        const int seed = static_cast<int>(frag.v0);
        const Pixel got =
            gs.shade_pixel(seed, frag.pos, lattice.get(), grid, x);
        const Pixel expected =
            gs.shade_pixel_full(seed, frag.pos, lattice.get(), grid, x);
        const int dr = channel_gap(got.r, expected.r);
        const int dg = channel_gap(got.g, expected.g);
        const int db = channel_gap(got.b, expected.b);
        const bool got_lit = !is_black(got);
        const bool expected_lit = !is_black(expected);
        error.coverage += got_lit != expected_lit;
        error.lit += got_lit || expected_lit;
        error.max_channel =
            std::max(error.max_channel, std::max(dr, std::max(dg, db)));
        error.total_channel += dr + dg + db;
      }
    }
    return error;
  }

  /** @brief Selects the reference shader `shared_shader_error` compares. */
  struct SharedShaderOptions {
    bool fixed_stencil = false;
    bool direct_color = false;
    bool aggregate_color = true;
    bool center_pigment = false;
    bool full_kernel = false;
  };

  template <typename Grid>
  static Pixel
  per_sample_reference(GS &gs, Grid &grid, int x, int seed, int shared_center,
                       math::Vector *world_nodes, float noise, float shift,
                       float lightness, ShaderError &error,
                       SharedShaderOptions options) {
    Pixel expected(0, 0, 0);
    for (int i = 0; i < Grid::SAMPLES; ++i) {
      math::Vector sample_v = grid.at(x, i);
      ++error.stencil_samples;
      if (gs.refine_center(sample_v, world_nodes, seed) != shared_center)
        ++error.stencil_changes;
      float b;
      if (options.fixed_stencil) {
        float tw = 0.0f, wb = 0.0f;
        gs.kernel_accumulate(sample_v, world_nodes, shared_center,
                             [&](int ni, float w) {
                               wb += gs.state.B[ni] * w;
                               tw += w;
                             });
        b = tw <= GS::KERNEL_MIN_TOTAL_WEIGHT ? 0.0f : wb * (GS::Q16_INV / tw);
      } else {
        b = interpolate_b(gs, sample_v, seed, world_nodes);
      }
      if (b < GS::B_CULL_THRESHOLD)
        continue;
      float t =
          hs::clamp((b - GS::B_COLOR_FLOOR) * GS::B_COLOR_SCALE, 0.0f, 1.0f);
      float mass = 0.0f, rgb[3] = {};
      int center = options.fixed_stencil
                       ? shared_center
                       : gs.refine_center(sample_v, world_nodes, seed);
      gs.kernel_accumulate(sample_v, world_nodes, center, [&](int ni, float w) {
        float weight = gs.state.B[ni] * w;
        Pixel color = mix_pigment(gs.state.pigment[ni], [&](int id) {
          return options.direct_color
                     ? gs.modified_palette_color(id, t, shift, lightness)
                     : cached_palette_color(gs, id, t, noise);
        });
        mass += weight;
        rgb[0] += color.r * weight;
        rgb[1] += color.g * weight;
        rgb[2] += color.b * weight;
      });
      if (mass > 0.0f)
        expected +=
            Pixel(round_linear_channel(rgb[0] / (mass * Grid::SAMPLES)),
                  round_linear_channel(rgb[1] / (mass * Grid::SAMPLES)),
                  round_linear_channel(rgb[2] / (mass * Grid::SAMPLES)));
    }
    return expected;
  }

  template <typename Grid>
  static Pixel aggregate_reference(GS &gs, Grid &grid, int x, int seed,
                                   int shared_center, math::Vector *world_nodes,
                                   float noise, float shift, float lightness,
                                   ShaderError &error,
                                   SharedShaderOptions options) {
    Pixel expected(0, 0, 0);
    double palette_mass[1 << Pigment::ID_BITS] = {};
    double total_mass = 0.0, concentration = 0.0;
    int covered = 0;
    for (int i = 0; i < Grid::SAMPLES; ++i) {
      const math::Vector sample = grid.at(x, i);
      int center = options.fixed_stencil
                       ? shared_center
                       : gs.refine_center(sample, world_nodes, seed);
      ++error.stencil_samples;
      error.stencil_changes += center != shared_center;
      float tw = 0.0f, wb = 0.0f;
      gs.kernel_accumulate(sample, world_nodes, center, [&](int node, float w) {
        const float mass = gs.state.B[node] * w;
        tw += w;
        wb += mass;
        const auto pigment = Pigment::unpack(gs.state.pigment[node]);
        const double first_fraction = first_weight16(pigment.mass) / 65535.0;
        palette_mass[pigment.first] += mass * first_fraction;
        palette_mass[pigment.second] += mass * (1.0 - first_fraction);
        total_mass += mass;
      });
      float b =
          tw <= GS::KERNEL_MIN_TOTAL_WEIGHT ? 0.0f : wb * (GS::Q16_INV / tw);
      if (b >= GS::B_CULL_THRESHOLD) {
        ++covered;
        concentration +=
            hs::clamp((b - GS::B_COLOR_FLOOR) * GS::B_COLOR_SCALE, 0.0f, 1.0f);
      }
    }
    if (options.center_pigment) {
      std::fill(std::begin(palette_mass), std::end(palette_mass), 0.0);
      const auto pigment = Pigment::unpack(gs.state.pigment[shared_center]);
      const double first_fraction = first_weight16(pigment.mass) / 65535.0;
      palette_mass[pigment.first] += first_fraction;
      palette_mass[pigment.second] += 1.0 - first_fraction;
      total_mass = 1.0;
    }
    if (covered && total_mass > 0.0) {
      float t = static_cast<float>(concentration / covered);
      double rgb[3] = {};
      for (int palette = 0; palette < GS::NUM_SEED_CLUSTERS; ++palette) {
        Pixel color =
            options.direct_color
                ? gs.modified_palette_color(palette, t, shift, lightness)
                : cached_palette_color(gs, palette, t, noise);
        const double weight =
            palette_mass[palette] * covered / (Grid::SAMPLES * total_mass);
        rgb[0] += color.r * weight;
        rgb[1] += color.g * weight;
        rgb[2] += color.b * weight;
      }
      expected =
          Pixel(round_linear_channel(rgb[0]), round_linear_channel(rgb[1]),
                round_linear_channel(rgb[2]));
    }
    return expected;
  }

  template <int W, int H>
  static ShaderError shared_shader_error(GS &gs, SharedShaderOptions options) {
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
        Pixel got = options.full_kernel
                        ? gs.shade_pixel_full(seed, frag.pos, world_nodes, grid,
                                              x, hot1.data())
                        : gs.shade_pixel(seed, frag.pos, world_nodes, grid, x,
                                         hot1.data());
        float noise =
            gs.sample_color_noise(gs.inverse_orientation.apply(frag.pos));
        float shift = noise * gs.params.hue_shift;
        float lightness = std::max(0.0f, noise) * gs.params.shimmer;

        const Pixel expected =
            options.aggregate_color
                ? aggregate_reference(gs, grid, x, seed, shared_center,
                                      world_nodes, noise, shift, lightness,
                                      error, options)
                : per_sample_reference(gs, grid, x, seed, shared_center,
                                       world_nodes, noise, shift, lightness,
                                       error, options);
        const int dr = channel_gap(got.r, expected.r);
        const int dg = channel_gap(got.g, expected.g);
        const int db = channel_gap(got.b, expected.b);
        const int pixel_max = std::max(dr, std::max(dg, db));
        const bool got_lit = !is_black(got);
        const bool expected_lit = !is_black(expected);
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
    return error;
  }
};
