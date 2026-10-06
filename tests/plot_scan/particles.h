/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ============================================================================
// Plot::ParticleSystem — trail rasterization
// ============================================================================

/** @brief Minimal particle for the ParticleSystem draw concept: trail + life. */
struct StubParticle {
  Animation::QuantizedVectorTrail<23>
      history;       /**< World-space trail positions. */
  uint16_t life = 0; /**< Remaining life (drives v3). */
};

/** @brief Minimal pool/active_count/max_life triple the draw concept reads. */
struct StubSystem {
  std::vector<StubParticle>
      pool;              /**< Pool; pool[i] read for i < active_count. */
  int active_count = 0;  /**< Live prefix length. */
  uint16_t max_life = 0; /**< Life normaliser for v3. */
  int active() const { return active_count; }
};

/** @brief Shader and output counts from drawing one particle. */
struct ParticleDrawCapture {
  int vertex_calls = 0;
  int deferred_calls = 0;
  int fragment_calls = 0;
  size_t plotted = 0;
};

/** @brief Draws one particle through every trail shader stage. */
inline ParticleDrawCapture capture_particle_draw(const StubParticle &particle) {
  constexpr int W = 96, H = 48;
  hs_test::StubEffect fx(W, H);
  StubSystem sys;
  sys.max_life = 100;
  sys.active_count = 1;
  sys.pool.push_back(particle);

  ParticleDrawCapture capture;
  CapturePipeline pipe;
  auto fragment_shader = [&](const math::Vector &, Fragment &) {
    ++capture.fragment_calls;
  };
  auto vertex_shader = [&](Fragment &) { ++capture.vertex_calls; };
  auto deferred_shader = [&](FragmentRegisters, const math::Vector &) {
    ++capture.deferred_calls;
  };
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(pipe, c, sys, fragment_shader,
                                     vertex_shader, deferred_shader);
  }
  capture.plotted = pipe.plotted.size();
  return capture;
}

/** @brief Captures the control-point fragments passed through the vertex stage.
 */
inline std::vector<Fragment>
capture_particle_vertices(const StubParticle &particle) {
  constexpr int W = 96, H = 48;
  hs_test::StubEffect fx(W, H);
  StubSystem sys;
  sys.max_life = 100;
  sys.active_count = 1;
  sys.pool.push_back(particle);

  std::vector<Fragment> vertices;
  CapturePipeline pipe;
  auto vertex_shader = [&](Fragment &f) { vertices.push_back(f); };
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(pipe, c, sys, noop_shader, vertex_shader);
  }
  return vertices;
}

/** @brief Captures control points emitted by an initialized particle system. */
template <typename SystemT>
inline std::vector<Fragment>
capture_particle_system_vertices(const SystemT &system) {
  constexpr int W = 96, H = 48;
  hs_test::StubEffect fx(W, H);
  std::vector<Fragment> vertices;
  CapturePipeline pipe;
  auto vertex_shader = [&](Fragment &f) { vertices.push_back(f); };
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(pipe, c, system, noop_shader,
                                     vertex_shader);
  }
  return vertices;
}

/** @brief Builds a quantized particle trail with deterministic spherical
 * points. */
inline StubParticle make_particle_trail(int samples) {
  StubParticle particle;
  particle.life = 60;
  for (int i = 0; i < samples; ++i) {
    float theta = 0.35f + 0.055f * i;
    float y = -0.3f + 0.02f * i;
    float radial = std::sqrt(1.0f - y * y);
    particle.history.record(
        math::Vector(radial * std::cos(theta), y, radial * std::sin(theta)));
  }
  return particle;
}

/** @brief Renders a particle through the current or callback materializer. */
inline std::vector<Pixel>
render_particle_materialization(const StubParticle &particle,
                                bool callback_reference) {
  constexpr int W = 96, H = 48;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> filters;
  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(50000, 30000, 10000),
                     hs::clamp(std::min(f.v0, f.v3), 0.0f, 1.0f));
  };
  {
    Canvas c(fx);
    if (callback_reference) {
      ScratchScope trail_guard(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a,
                 std::remove_cvref_t<decltype(particle.history)>::CAPACITY);
      const float inv_max_life = 1.0f / 100.0f;
      tween(particle.history, [&](const math::Vector &v, float t) {
        Fragment f;
        f.pos = v;
        f.v0 = t;
        f.v1 = 0.0f;
        f.v2 = 0.0f;
        f.v3 = static_cast<float>(particle.life) * inv_max_life;
        f.age = 0;
        f.color = Color4(0, 0, 0, 0);
        trail.push_back(f);
      });
      Plot::rasterize<W, H>(filters, c, trail, shade);
    } else {
      StubSystem sys;
      sys.max_life = 100;
      sys.active_count = 1;
      sys.pool.push_back(particle);
      Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade);
    }
  }
  fx.advance_display();

  std::vector<Pixel> pixels(static_cast<size_t>(W) * H);
  hs_test::capture_frame<W, H>(fx, pixels);
  return pixels;
}

/**
 * @brief Verifies ParticleSystem::draw rasterizes only the active prefix's trails
 *        and stamps the per-particle registers (v2 source index, v3 life ratio).
 * @details Pins three things the smoke loop never checks: (1) only particles in
 * [0, active_count) are drawn — an inactive particle parked on the ±Y poles never
 * contributes a point; (2) the drawn trail follows the particle's recorded
 * history (an equatorial +X→+Z arc); (3) every emitted fragment carries the
 * source index in v2 and life/max_life in v3, constant across the trail when
 * draw() uses its default register mappers.
 */
inline void test_particle_system_draws_active_trails_with_registers() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);

  StubSystem sys;
  sys.max_life = 100;
  sys.active_count = 1; // pool holds 2; only particle 0 is active

  StubParticle p0;
  p0.life = 60;
  const math::Vector t0[3] = {math::Vector(1, 0, 0),
                              math::Vector(0.7071f, 0.0f, 0.7071f),
                              math::Vector(0, 0, 1)}; // equatorial +X -> +Z arc
  for (const math::Vector &v : t0)
    p0.history.record(v);

  StubParticle p1; // inactive: parked on the poles, must NOT be drawn
  p1.life = 99;
  p1.history.record(math::Vector(0, 1, 0));
  p1.history.record(math::Vector(0, -1, 0));

  sys.pool.push_back(p0);
  sys.pool.push_back(p1);

  CapturePipeline pipe;
  float v2_lo = 1e9f, v2_hi = -1e9f, v3_lo = 1e9f, v3_hi = -1e9f;
  int nonfinite_registers = 0;
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(
        pipe, c, sys, [&](const math::Vector &, Fragment &f) {
          nonfinite_registers += !std::isfinite(f.v2) || !std::isfinite(f.v3);
          v2_lo = std::min(v2_lo, f.v2);
          v2_hi = std::max(v2_hi, f.v2);
          v3_lo = std::min(v3_lo, f.v3);
          v3_hi = std::max(v3_hi, f.v3);
        });
  }
  fx.advance_display();

  // (1)+(2) Active trail follows its recorded arc; the inactive ±Y particle is absent.
  HS_EXPECT_GT(pipe.plotted.size(), (size_t)2);
  for (const math::Vector &v : pipe.plotted) {
    HS_EXPECT_LE(arc_angular_distance(v, t0[0], t0[2]), 0.05f);
    HS_EXPECT_LT(std::fabs(v.y), 0.1f);
  }
  // (3) Registers: v2 == source index 0; v3 == life/max_life == 0.6, constant.
  HS_EXPECT_NEAR(v2_lo, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(v2_hi, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(v3_lo, 0.6f, 1e-3f);
  HS_EXPECT_NEAR(v3_hi, 0.6f, 1e-3f);
  HS_EXPECT_EQ(nonfinite_registers, 0);
}

/**
 * @brief Verifies an empty system needs no lifetime normalizer.
 */
inline void test_particle_system_empty_zero_lifetime_is_noop() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  StubSystem sys;
  CapturePipeline pipe;
  {
    Canvas c(fx);
    Plot::ParticleSystem::draw<W, H>(pipe, c, sys, noop_shader);
  }
  fx.advance_display();
  HS_EXPECT_TRUE(pipe.plotted.empty());
}

/**
 * @brief Dead and one-point particle trails bypass all render stages while a
 *        live two-point trail renders normally.
 */
inline void test_particle_system_skips_unrenderable_trails() {
  StubParticle dead;
  dead.history.record(math::Vector(1, 0, 0));
  dead.history.record(math::Vector(0, 0, 1));
  ParticleDrawCapture dead_capture = capture_particle_draw(dead);
  HS_EXPECT_EQ(dead_capture.vertex_calls, 0);
  HS_EXPECT_EQ(dead_capture.deferred_calls, 0);
  HS_EXPECT_EQ(dead_capture.fragment_calls, 0);
  HS_EXPECT_EQ(dead_capture.plotted, static_cast<size_t>(0));

  StubParticle one_point;
  one_point.life = 60;
  one_point.history.record(math::Vector(1, 0, 0));
  ParticleDrawCapture one_point_capture = capture_particle_draw(one_point);
  HS_EXPECT_EQ(one_point_capture.vertex_calls, 0);
  HS_EXPECT_EQ(one_point_capture.deferred_calls, 0);
  HS_EXPECT_EQ(one_point_capture.fragment_calls, 0);
  HS_EXPECT_EQ(one_point_capture.plotted, static_cast<size_t>(0));

  StubParticle trail;
  trail.life = 60;
  trail.history.record(math::Vector(1, 0, 0));
  trail.history.record(math::Vector(0, 0, 1));
  ParticleDrawCapture trail_capture = capture_particle_draw(trail);
  HS_EXPECT_EQ(trail_capture.vertex_calls, 2);
  HS_EXPECT_EQ(trail_capture.deferred_calls, 2);
  HS_EXPECT_GT(trail_capture.fragment_calls, 0);
  HS_EXPECT_GT(trail_capture.plotted, static_cast<size_t>(0));
}

/**
 * @brief Direct materialization preserves circular order, progress, and
 *        particle registers for partial, full, and wrapped histories.
 */
inline void test_particle_system_direct_trail_materialization_registers() {
  const int sample_counts[] = {7, 23, 30};
  for (int sample_count : sample_counts) {
    StubParticle particle = make_particle_trail(sample_count);
    std::vector<Fragment> vertices = capture_particle_vertices(particle);
    const size_t len = particle.history.length();
    HS_EXPECT_SIZE_OR_RETURN(vertices, len);
    for (size_t i = 0; i < len; ++i) {
      math::Vector expected = particle.history.get(i);
      HS_EXPECT_EQ(vertices[i].pos.x, expected.x);
      HS_EXPECT_EQ(vertices[i].pos.y, expected.y);
      HS_EXPECT_EQ(vertices[i].pos.z, expected.z);
      HS_EXPECT_NEAR(vertices[i].v0,
                     static_cast<float>(i) / static_cast<float>(len - 1),
                     std::numeric_limits<float>::epsilon());
      HS_EXPECT_EQ(vertices[i].v1, 0.0f);
      HS_EXPECT_EQ(vertices[i].v2, 0.0f);
      HS_EXPECT_EQ(vertices[i].v3, 60.0f * (1.0f / 100.0f));
    }
  }
}

/** @brief Sparse histories draw the exact live position between anchor frames. */
inline void test_particle_system_sparse_history_live_tip() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<96, 1, 4, 8, 8, false, 6> system;
  system.init(arena, /*friction=*/0.85f, /*gravity=*/0.0f,
              /*max_life=*/100.0f);
  system.spawn(math::Vector(0, 0, 1), math::Vector(0, 0, 0), 0);

  auto &particle = system.pool[0];
  particle.history.record(math::Vector(1, 0, 0));
  particle.history.record(math::Vector(0, 1, 0));
  particle.position = math::Vector(0, 0, 1);
  particle.life = 50;

  std::vector<Fragment> vertices = capture_particle_system_vertices(system);
  HS_EXPECT_SIZE_OR_RETURN(vertices, (size_t)3);
  HS_EXPECT_EQ(vertices[0].v0, 0.0f);
  HS_EXPECT_EQ(vertices[1].v0, 0.5f);
  HS_EXPECT_EQ(vertices[2].v0, 1.0f);
  HS_EXPECT_EQ(vertices[2].pos.x, particle.position.x);
  HS_EXPECT_EQ(vertices[2].pos.y, particle.position.y);
  HS_EXPECT_EQ(vertices[2].pos.z, particle.position.z);

  particle.life = 51;
  vertices = capture_particle_system_vertices(system);
  HS_EXPECT_SIZE_OR_RETURN(vertices, (size_t)2);
  HS_EXPECT_EQ(vertices.back().v0, 1.0f);
}

/**
 * @brief v0 ramps 0 at the OLDEST retained sample to 1 at the newest (head).
 * @details Pins the register's orientation, not just its spacing: consumers fade
 * on v0 (MindSplatter's min(v0, v3)), so a reversed ramp would fade the head.
 * Uses a wrapped history, where the oldest survivor is not the first record.
 */
inline void test_particle_system_v0_zero_at_oldest_sample() {
  StubParticle particle;
  particle.life = 60;
  constexpr size_t CAP =
      std::remove_cvref_t<decltype(particle.history)>::CAPACITY;
  const size_t recorded = CAP + 5;
  std::vector<math::Vector> order;
  for (size_t i = 0; i < recorded; ++i) {
    float theta = 0.2f + 0.04f * static_cast<float>(i);
    math::Vector v(std::cos(theta), 0.0f, std::sin(theta));
    particle.history.record(v);
    order.push_back(v);
  }

  std::vector<Fragment> vertices = capture_particle_vertices(particle);
  HS_EXPECT_SIZE_OR_RETURN(vertices, CAP);

  const math::Vector oldest = order[recorded - CAP];
  const math::Vector newest = order.back();
  HS_EXPECT_GT(math::angle_between(oldest, newest), 0.1f);
  HS_EXPECT_EQ(vertices.front().v0, 0.0f);
  HS_EXPECT_EQ(vertices.back().v0, 1.0f);
  HS_EXPECT_NEAR(math::angle_between(vertices.front().pos, oldest), 0.0f,
                 1e-3f);
  HS_EXPECT_NEAR(math::angle_between(vertices.back().pos, newest), 0.0f, 1e-3f);
}

/** @brief A custom v2 mapper runs once per materialized particle. */
inline void test_particle_system_custom_v2_mapper() {
  constexpr int W = 96, H = 48;
  hs_test::StubEffect fx(W, H);
  StubSystem sys;
  sys.max_life = 100;

  StubParticle first = make_particle_trail(3);
  StubParticle second = make_particle_trail(4);
  first.life = 25;
  second.life = 75;
  StubParticle dead = make_particle_trail(2);
  dead.life = 0;
  StubParticle one_point = make_particle_trail(1);
  sys.pool = {first, second, dead, one_point};
  sys.active_count = static_cast<int>(sys.pool.size());

  int mapper_calls = 0;
  int first_vertices = 0;
  int second_vertices = 0;
  int unclassified_vertices = 0;
  auto particle_v2 = [&](const StubParticle &p, int) {
    ++mapper_calls;
    return static_cast<float>(p.life) / 100.0f;
  };
  // Bucket by tolerance, not by ==: the shipping -ffast-math pair may turn the
  // mapper's division into a reciprocal multiply, and an exact compare would
  // then silently drop every vertex into neither bucket — failing the counts
  // below with a diagnostic that blames the vertex walk.
  auto vertex_shader = [&](Fragment &f) {
    if (hs_test::approx(f.v2, 0.25f, 1e-4f))
      ++first_vertices;
    else if (hs_test::approx(f.v2, 0.75f, 1e-4f))
      ++second_vertices;
    else
      ++unclassified_vertices;
  };
  CapturePipeline pipeline;
  {
    Canvas canvas(fx);
    Plot::ParticleSystem::draw<W, H>(pipeline, canvas, sys, noop_shader,
                                     vertex_shader, DeferredShaderRef{},
                                     particle_v2);
  }

  HS_EXPECT_EQ(mapper_calls, 2);
  HS_EXPECT_EQ(unclassified_vertices, 0);
  HS_EXPECT_EQ(first_vertices, 3);
  HS_EXPECT_EQ(second_vertices, 4);
}

/**
 * @brief Direct and callback materializers produce identical framebuffers for
 *        full linear and full wrapped quantized histories.
 */
inline void test_particle_system_direct_trail_materialization_output_parity() {
  const int sample_counts[] = {23, 30};
  for (int sample_count : sample_counts) {
    StubParticle particle = make_particle_trail(sample_count);
    std::vector<Pixel> direct =
        render_particle_materialization(particle, false);
    std::vector<Pixel> reference =
        render_particle_materialization(particle, true);
    int lit = 0;
    int diff = 0;
    for (size_t i = 0; i < direct.size(); ++i) {
      if (reference[i].r | reference[i].g | reference[i].b)
        ++lit;
      if (direct[i] != reference[i])
        ++diff;
    }
    HS_EXPECT_GT(lit, 0);
    HS_EXPECT_EQ(diff, 0);
  }
}

/** @brief Reference-lit and mismatching pixel counts over a canvas rectangle. */
struct BandDiff {
  int lit = 0;  /**< Reference pixels lit inside the rectangle. */
  int diff = 0; /**< Pixels differing from the reference. */
};

/**
 * @brief Compares a rendered frame against a reference over [y0,y1) x [x0,x1).
 * @tparam W Canvas width, the row stride of @p ref.
 * @param fx Effect holding the frame under test (already advance_display'd).
 * @param ref Full-canvas reference pixels.
 */
template <int W>
inline BandDiff band_diff(const hs_test::StubEffect &fx,
                          const std::vector<Pixel> &ref, int y0, int y1, int x0,
                          int x1) {
  BandDiff out;
  for (int y = y0; y < y1; ++y)
    for (int x = x0; x < x1; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      const Pixel &r = ref[static_cast<size_t>(y) * W + x];
      if (r.r | r.g | r.b)
        ++out.lit;
      if (p.r != r.r || p.g != r.g || p.b != r.b)
        ++out.diff;
    }
  return out;
}

/**
 * @brief Renders a system full-canvas through one combined vertex shader: the
 *        reference the band-clipped parity checks diff against.
 */
template <int W, int H, typename Shade, typename Combined>
inline std::vector<Pixel> particle_reference_frame(StubSystem &sys, Shade shade,
                                                   Combined combined) {
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> filters;
  {
    Canvas c(fx);
    initialize_parity_frame<W, H>(c);
    Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, combined);
  }
  fx.advance_display();
  std::vector<Pixel> ref(static_cast<size_t>(W) * H);
  hs_test::capture_frame<W, H>(fx, ref);
  return ref;
}

/** @brief Clip bands the 96x48 ParticleSystem parity sweeps run, as
 * {y0, y1, x0, x1}. */
inline constexpr int PARTICLE_PARITY_BANDS[4][4] = {
    {0, 24, 0, 48},   // quadrant; margin wraps rs past the seam
    {12, 36, 48, 96}, // opposite half
    {0, 48, 20, 70},  // interior wedge
    {30, 48, 0, 96},  // y-only clip (XClip inactive)
};

/**
 * @brief Deferred trail shader: bit-identical in-band pixels to an undeferred
 *        combined shader, skipped whole for trails whose every edge is culled,
 *        and handed the original pre-shader positions.
 * @details Two equatorial trails under an x-wedge clip: trail 0 crosses the
 *          band edge (mixed per-edge verdicts exercise the precomputed-bits
 *          path through rasterize), trail 1 lies wholly outside (its deferred
 *          pass must never run). The position pass negates the sphere, so the
 *          deferred shader can verify its `orig` argument is the pre-shader
 *          position, not the shaded one. Reference is a single combined vertex
 *          shader on a full canvas; the clipped deferred render must match it
 *          exactly inside the render band, including its margin ring.
 */
inline void test_particle_system_deferred_shader_parity_and_skip() {
  constexpr int W = 96, H = 48;
  const int band_x0 = 0, band_x1 = W / 2;

  // Shaded-space equatorial arcs (theta -> column x = theta * W / 2pi):
  // trail 0 spans columns ~40..58 (crosses band_x1 = 48), trail 1 ~64..76
  // (wholly outside [0,49) incl. the +-1 render margin, clear of the seam).
  auto shaded = [](float theta) {
    return math::Vector(cosf(theta), 0, sinf(theta));
  };
  const float T0[5] = {2.6f, 2.9f, 3.2f, 3.5f, 3.8f};
  const float T1[3] = {4.2f, 4.6f, 5.0f};

  StubSystem sys;
  sys.max_life = 100;
  sys.active_count = 2;
  StubParticle p0, p1;
  p0.life = 60;
  p1.life = 60;
  // History records ORIGINAL positions; the position pass negates them.
  for (float t : T0)
    p0.history.record(shaded(t) * -1.0f);
  for (float t : T1)
    p1.history.record(shaded(t) * -1.0f);
  sys.pool.push_back(p0);
  sys.pool.push_back(p1);

  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), hs::clamp(f.v3, 0.0f, 1.0f));
  };
  std::vector<math::Vector> shaded_positions;
  auto position_pass = [&](Fragment &f) {
    f.pos = f.pos * -1.0f;
    shaded_positions.push_back(f.pos);
  };
  int deferred_calls[2] = {0, 0};
  int orig_mismatches = 0;
  auto deferred_pass = [&](FragmentRegisters f, const math::Vector &orig) {
    // orig must be a pre-shader position: the negation of a shaded one.
    bool matched = false;
    for (const math::Vector &s : shaded_positions)
      if ((orig * -1.0f - s).length() <= 1e-4f) {
        matched = true;
        break;
      }
    if (!matched)
      orig_mismatches++;
    const float index = f.v2 + 0.5f;
    HS_EXPECT_TRUE(index >= 0.0f && index < 2.0f);
    if (!(index >= 0.0f && index < 2.0f))
      return;
    deferred_calls[static_cast<size_t>(index)]++;
    f.v3 *= 0.5f;
  };
  auto combined = [](Fragment &f) {
    f.pos = f.pos * -1.0f;
    f.v3 *= 0.5f;
  };

  const std::vector<Pixel> ref =
      particle_reference_frame<W, H>(sys, shade, combined);

  // Clipped, split shaders: in-band pixels identical; trail 1 never shaded.
  {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(0, H, band_x0, band_x1);
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, position_pass,
                                       deferred_pass);
    }
    fx.advance_display();
    const RenderBandDiff d = render_band_diff<W>(fx, ref);
    HS_EXPECT_GT(d.lit, 0);
    HS_EXPECT_GT(d.margin_lit, 0);
    expect_render_band_parity("particle deferred shader", d);
    HS_EXPECT_GT(deferred_calls[0], 0); // crossing trail: deferred pass ran
    HS_EXPECT_EQ(deferred_calls[1], 0); // fully-culled trail: skipped whole
    HS_EXPECT_EQ(orig_mismatches, 0);
  }

  // Full canvas, split shaders: identical everywhere, both trails shaded.
  {
    deferred_calls[0] = deferred_calls[1] = 0;
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, position_pass,
                                       deferred_pass);
    }
    fx.advance_display();
    HS_EXPECT_EQ((band_diff<W>(fx, ref, 0, H, 0, W).diff), 0);
    HS_EXPECT_GT(deferred_calls[0], 0);
    HS_EXPECT_GT(deferred_calls[1], 0);
    HS_EXPECT_EQ(orig_mismatches, 0);
  }
}

/**
 * @brief Randomized whole-trail gate conservativeness: a band-clipped
 *        ParticleSystem render is pixel-identical to the full render inside
 *        the render band, including its margin ring.
 * @details Random-walk trails, salted with pole-crossing and near-antipodal
 *          steps, drive the hoisted gate's coarse trail reject and per-edge
 *          bits through a hoistable pipeline; a cull false-negative drops
 *          in-band pixels. Bands cover seam-wrapping, interior, and y-only
 *          clips.
 */
inline void test_particle_system_gate_pixel_parity_random_trails() {
  constexpr int W = 96, H = 48;
  hs::random().seed(20260717);

  StubSystem sys;
  sys.max_life = 100;
  const int NT = 60;
  for (int t = 0; t < NT; ++t) {
    StubParticle p;
    p.life = static_cast<uint16_t>(40 + (t % 60));
    math::Vector v =
        (t % 5 == 0) ? math::Vector(0, t % 10 == 0 ? 1 : -1, 0) : rand_unit();
    for (int k = 0; k < 12; ++k) {
      p.history.record(v);
      const float step_x = hs::rand_f(-1, 1);
      const float step_y = hs::rand_f(-1, 1);
      const float step_z = hs::rand_f(-1, 1);
      math::Vector step(step_x, step_y, step_z);
      // Occasional huge step: a near-antipodal edge must trip the coarse
      // walk's half-sweep guard, not get mis-culled.
      float scale = (k == 7 && t % 9 == 0) ? 4.0f : 0.12f;
      v = (v + step * scale).normalized();
    }
    sys.pool.push_back(p);
  }
  sys.active_count = NT;

  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), hs::clamp(f.v3, 0.0f, 1.0f));
  };
  auto position_pass = [](Fragment &f) { f.pos = f.pos * -1.0f; };
  auto deferred_pass = [](FragmentRegisters f, const math::Vector &) {
    f.v3 *= 0.7f;
  };
  auto combined = [](Fragment &f) {
    f.pos = f.pos * -1.0f;
    f.v3 *= 0.7f;
  };

  const std::vector<Pixel> ref =
      particle_reference_frame<W, H>(sys, shade, combined);

  int margin_lit = 0;
  for (const auto &bd : PARTICLE_PARITY_BANDS) {
    // The non-deferred path must use the same precomputed gate bits.
    {
      hs_test::StubEffect fx(W, H);
      fx.set_clip(bd[0], bd[1], bd[2], bd[3]);
      Pipeline<W, H> filters;
      {
        Canvas c(fx);
        initialize_parity_frame<W, H>(c);
        Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, combined);
      }
      fx.advance_display();
      const RenderBandDiff diff = render_band_diff<W>(fx, ref);
      margin_lit += diff.margin_lit;
      expect_render_band_parity("particle combined shader", diff);
    }

    hs_test::StubEffect fx(W, H);
    fx.set_clip(bd[0], bd[1], bd[2], bd[3]);
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, position_pass,
                                       deferred_pass);
    }
    fx.advance_display();
    const RenderBandDiff d = render_band_diff<W>(fx, ref);
    margin_lit += d.margin_lit;
    HS_EXPECT_GT(d.lit, 20);
    expect_render_band_parity("particle split shader", d);
  }
  HS_EXPECT_GT(margin_lit, 20);
}

/**
 * @brief Verifies a clipped sub-pixel trail — the population the single-dot
 *        shortcut routes — is pixel-identical to the unclipped reference.
 * @details The gate-parity sweep above steps trails at ~2x base_step, so its
 *          edges mostly take the sampler. These trails step well under one
 *          screen step at several latitudes, so the shortcut (and its reuse of
 *          the gate's precomputed rows/columns) decides the pixels.
 */
inline void test_particle_system_subpixel_trail_dot_parity() {
  constexpr int W = 96, H = 48;
  constexpr float base_step = (2.0f * math::PI_F) / W;
  hs::random().seed(20260720);

  StubSystem sys;
  sys.max_life = 100;
  const int NT = 40;
  for (int t = 0; t < NT; ++t) {
    StubParticle p;
    p.life = static_cast<uint16_t>(50 + t % 40);
    float lat = -1.2f + 2.4f * (static_cast<float>(t) / (NT - 1));
    float az = hs::rand_f(0.0f, 2.0f * math::PI_F);
    math::Vector v(std::cos(lat) * std::cos(az), std::sin(lat),
                   std::cos(lat) * std::sin(az));
    for (int k = 0; k < 12; ++k) {
      p.history.record(v);
      const float step_x = hs::rand_f(-1, 1);
      const float step_y = hs::rand_f(-1, 1);
      const float step_z = hs::rand_f(-1, 1);
      math::Vector step(step_x, step_y, step_z);
      v = (v + step * (base_step * 0.2f)).normalized();
    }
    sys.pool.push_back(p);
  }
  sys.active_count = NT;

  auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 40000, 20000), hs::clamp(f.v3, 0.0f, 1.0f));
  };
  auto position_pass = [](Fragment &f) { f.pos = f.pos * -1.0f; };
  auto deferred_pass = [](FragmentRegisters f, const math::Vector &) {
    f.v3 *= 0.7f;
  };
  auto combined = [](Fragment &f) {
    f.pos = f.pos * -1.0f;
    f.v3 *= 0.7f;
  };

  const std::vector<Pixel> ref =
      particle_reference_frame<W, H>(sys, shade, combined);

  int margin_lit = 0;
  for (const auto &bd : PARTICLE_PARITY_BANDS) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(bd[0], bd[1], bd[2], bd[3]);
    Pipeline<W, H> filters;
    {
      Canvas c(fx);
      initialize_parity_frame<W, H>(c);
      Plot::ParticleSystem::draw<W, H>(filters, c, sys, shade, position_pass,
                                       deferred_pass);
    }
    fx.advance_display();
    const RenderBandDiff d = render_band_diff<W>(fx, ref);
    margin_lit += d.margin_lit;
    HS_EXPECT_GT(d.lit, 10);
    expect_render_band_parity("particle subpixel trail", d);
  }
  HS_EXPECT_GE(margin_lit, 10);
}
