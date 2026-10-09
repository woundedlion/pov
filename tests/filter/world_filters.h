/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// World filters — plot/cull taps, Pipeline cull routing and Mobius canvas parity
// ============================================================================

/**
 * @brief A small recorder for the (vector, color, age, alpha) tuples a World
 *        filter emits.
 */
struct Tap3D {
  math::Vector v;   /**< Emitted world-space position. */
  Pixel c;          /**< Emitted colour. */
  float age, alpha; /**< Emitted age and alpha. */
};

/**
 * @brief Verifies Hole masks a spherical cap: outside the radius the point
 *        passes through untouched; inside, alpha is scaled by
 *        quintic_kernel(d/r) and the very center emits nothing at all.
 */
inline void test_world_hole_masks_cap() {
  Filter::World::Hole hole(math::Vector(0, 1, 0),
                           0.5f); // cap at +Y, radius 0.5 rad

  // Far point (south pole) is well outside -> verbatim passthrough.
  int n = 0;
  Tap3D got{};
  hole.plot(math::Vector(0, -1, 0), Pixel(10000, 20000, 30000), 4.0f, 0.8f,
            [&](const math::Vector &v, const Pixel &c, float age, float a) {
              got = {v, c, age, a};
              ++n;
            });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_EQ((int)got.c.r, 10000);
  HS_EXPECT_EQ((int)got.c.g, 20000);
  HS_EXPECT_NEAR(got.age, 4.0f, 1e-6f);
  HS_EXPECT_NEAR(got.alpha, 0.8f, 1e-6f);
  HS_EXPECT_NEAR(got.v.y, -1.0f, 1e-6f);

  // Exact center: d=0 -> quintic_kernel(0)=0 -> the tap is dropped.
  n = 0;
  hole.plot(math::Vector(0, 1, 0), Pixel(10000, 20000, 30000), 0.0f, 1.0f,
            [&](const math::Vector &, const Pixel &, float, float) { ++n; });
  HS_EXPECT_EQ(n, 0);

  // Half-radius (d = 0.25 rad): alpha scaled by quintic_kernel(0.5) = 0.5,
  // colour untouched.
  math::Vector half(sinf(0.25f), cosf(0.25f), 0.0f); // 0.25 rad from +Y
  n = 0;
  hole.plot(half, Pixel(10000, 20000, 30000), 0.0f, 1.0f,
            [&](const math::Vector &, const Pixel &c, float, float a) {
              got.c = c;
              got.alpha = a;
              ++n;
            });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_EQ((int)got.c.r, 10000);
  HS_EXPECT_EQ((int)got.c.g, 20000);
  HS_EXPECT_NEAR_REL(got.alpha, 0.5f, 0.01f);
}

/**
 * @brief Verifies Hole's origin and radius can be retuned after construction.
 */
inline void test_world_hole_setters() {
  math::Vector center(0, 1, 0); // start at +Y
  Filter::World::Hole hole(center, 0.5f);

  // At the initial center the point is fully masked (quintic_kernel(0) = 0).
  Tap3D got{};
  int n = 0;
  auto capture = [&](const math::Vector &, const Pixel &c, float, float a) {
    got.c = c;
    got.alpha = a;
    ++n;
  };
  hole.plot(center, Pixel(10000, 20000, 30000), 0.0f, 1.0f, capture);
  HS_EXPECT_EQ(n, 0);

  // The old center now lies well outside -> verbatim passthrough.
  center = math::Vector(0, -1, 0);
  hole.set_origin(center);
  hole.plot(math::Vector(0, 1, 0), Pixel(10000, 20000, 30000), 0.0f, 1.0f,
            capture);
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_EQ((int)got.c.r, 10000);
  HS_EXPECT_EQ((int)got.c.g, 20000);
  HS_EXPECT_NEAR(got.alpha, 1.0f, 1e-6f);

  // The new center emits nothing, confirming the mask followed the origin.
  hole.plot(center, Pixel(10000, 20000, 30000), 0.0f, 1.0f, capture);
  HS_EXPECT_EQ(n, 1);

  math::Vector quarter(sinf(0.25f), -cosf(0.25f), 0.0f);
  hole.set_radius(0.1f);
  hole.plot(quarter, Pixel(10000, 20000, 30000), 0.0f, 1.0f, capture);
  HS_EXPECT_EQ(n, 2);
  HS_EXPECT_NEAR(got.alpha, 1.0f, 1e-6f);

  hole.set_radius(0.5f);
  hole.plot(quarter, Pixel(10000, 20000, 30000), 0.0f, 1.0f, capture);
  HS_EXPECT_EQ(n, 3);
  HS_EXPECT_EQ((int)got.c.r, 10000);
  HS_EXPECT_NEAR(got.alpha, 0.5f, 0.005f);
}

/**
 * @brief Verifies Orient rotates by the bound Orientation, and a single-frame
 *        (stationary) orientation emits one tap with age left untouched.
 * @details A lone snapshot is the newest sub-position (t = 1), so the (1 - t)
 *          age offset is zero.
 */
inline void test_world_orient_rotates_and_keeps_static_age() {
  math::Quaternion q =
      math::make_rotation(math::Y_AXIS, math::PI_F / 2); // 90 deg about +Y
  const math::Orientation<> ori(q);
  Filter::World::Orient orient(ori);

  int n = 0;
  Tap3D got{};
  orient.plot(math::X_AXIS, Pixel(1, 2, 3), 5.0f, 1.0f,
              [&](const math::Vector &v, const Pixel &c, float age, float a) {
                got = {v, c, age, a};
                ++n;
              });
  HS_EXPECT_EQ(n, 1);
  math::Vector expected = math::rotate(math::X_AXIS, q);
  HS_EXPECT_NEAR(got.v.x, expected.x, 1e-4f);
  HS_EXPECT_NEAR(got.v.y, expected.y, 1e-4f);
  HS_EXPECT_NEAR(got.v.z, expected.z, 1e-4f);
  HS_EXPECT_NEAR(got.age, 5.0f, 1e-4f); // age + (1 - t), t = 1 (age-neutral)
}

/**
 * @brief Verifies that with a 3-frame history the tween sweeps the trailing
 *        sub-positions (i = 1..n-1) and spreads age across one frame.
 * @details Age offsets are (1 - t) for t in {0.5, 1.0}.
 */
inline void test_world_orient_motion_blur_sweep_ages() {
  math::Orientation<> ori; // identity, 1 frame
  ori.push(math::make_rotation(math::Y_AXIS, math::PI_F / 4));
  ori.push(math::make_rotation(math::Y_AXIS, math::PI_F / 2)); // now 3 frames
  Filter::World::Orient orient(ori);

  int n = 0;
  float ages[4] = {0};
  orient.plot(math::X_AXIS, Pixel(1, 1, 1), 10.0f, 1.0f,
              [&](const math::Vector &, const Pixel &, float age, float) {
                if (n < 4)
                  ages[n] = age;
                ++n;
              });
  // tween emits len-1 = 2 taps (i=1,2), t in {0.5, 1.0} -> age + {0.5, 0.0}.
  HS_EXPECT_EQ(n, 2);
  HS_EXPECT_NEAR(ages[0], 10.5f, 1e-4f);
  HS_EXPECT_NEAR(ages[1], 10.0f, 1e-4f);
}

/**
 * @brief Verifies Orient's clip-cull re-emits the edge under the same tween
 *        rotations plot() applies, short-circuits on the first hit, rotates a
 *        planar basis alongside the endpoints, and keeps an edge the bound
 *        orientation carries into the band.
 */
inline void test_world_orient_cull_edge_mirrors_plot() {
  static_assert(Filter::has_cull_edge<Filter::World::Orient>);
  math::Orientation<> ori; // identity, 1 frame
  ori.push(math::make_rotation(math::Y_AXIS, math::PI_F / 4));
  ori.push(math::make_rotation(math::Y_AXIS, math::PI_F / 2)); // now 3 frames
  const Filter::World::Orient orient(ori);

  math::Vector seen[2];
  int n = 0;
  bool hit = orient.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &, const math::Basis *) {
        if (n < 2)
          seen[n] = a;
        ++n;
        return false;
      });
  HS_EXPECT_FALSE(hit);
  // tween skips index 0, so the sweep is frames 1 and 2 — the same copies plot()
  // draws.
  HS_EXPECT_EQ(n, 2);
  for (int i = 0; i < 2; ++i) {
    math::Vector expected = math::rotate(math::X_AXIS, ori.get(i + 1));
    HS_EXPECT_NEAR(seen[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(seen[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(seen[i].z, expected.z, 1e-4f);
  }

  n = 0;
  hit = orient.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &, const math::Vector &, const math::Basis *) {
        ++n;
        return true;
      });
  HS_EXPECT_TRUE(hit);
  HS_EXPECT_EQ(n, 1);

  // A planar edge rotates its basis alongside the endpoints.
  const math::Basis pb{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  math::Vector normals[2];
  n = 0;
  orient.cull_edge(
      math::X_AXIS, math::X_AXIS, &pb,
      [&](const math::Vector &, const math::Vector &, const math::Basis *bp) {
        if (bp && n < 2)
          normals[n] = bp->v;
        ++n;
        return false;
      });
  HS_EXPECT_EQ(n, 2);
  for (int i = 0; i < 2; ++i) {
    math::Vector expected = math::rotate(math::Y_AXIS, ori.get(i + 1));
    HS_EXPECT_NEAR(normals[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(normals[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(normals[i].z, expected.z, 1e-4f);
  }

  // An edge that sits outside a polar band until the orientation rotates it in:
  // the band predicate rejects the source endpoints, plot() lands the point
  // inside, and the cull must therefore keep the edge.
  math::Orientation<> tilt(math::make_rotation(math::Z_AXIS, math::PI_F / 2));
  Filter::World::Orient tilted(tilt);
  auto in_polar_band = [](const math::Vector &a, const math::Vector &b,
                          const math::Basis *) {
    return std::fabs(a.y) > 0.9f && std::fabs(b.y) > 0.9f;
  };
  HS_EXPECT_FALSE(in_polar_band(math::X_AXIS, math::X_AXIS, nullptr));
  math::Vector plotted{};
  tilted.plot(
      math::X_AXIS, Pixel(1, 1, 1), 0.0f, 1.0f,
      [&](const math::Vector &v, const Pixel &, float, float) { plotted = v; });
  HS_EXPECT_GT(std::fabs(plotted.y), 0.9f);
  HS_EXPECT_TRUE(
      tilted.cull_edge(math::X_AXIS, math::X_AXIS, nullptr, in_polar_band));
}

/**
 * @brief Verifies OrientSlice picks an orientation from a list by the point's
 *        projection onto an axis, then rotates by it.
 * @details With two distinct orientations a point near +axis selects the last,
 *          near -axis selects the first; disabled is a passthrough.
 */
inline void test_world_orient_slice_selects_by_projection() {
  math::Orientation<> oris[2];
  oris[0].set(math::make_rotation(math::X_AXIS, math::PI_F / 2)); // index 0
  oris[1].set(math::make_rotation(math::Z_AXIS, math::PI_F / 2)); // index 1
  std::span<const math::Orientation<>> span(oris, 2);
  Filter::World::OrientSlice slice(span, math::Y_AXIS);

  auto first_tap = [&](const math::Vector &probe) {
    math::Vector out{};
    slice.plot(
        probe, Pixel(1, 1, 1), 0.0f, 1.0f,
        [&](const math::Vector &v, const Pixel &, float, float) { out = v; });
    return out;
  };

  // Probe near +Y (projection ~ +1 -> t ~ 1 -> last index 1 = Z rotation).
  math::Vector near_pos = math::Vector(0.15f, 0.98f, 0.0f).normalized();
  math::Vector exp_pos = math::rotate(near_pos, oris[1].get());
  math::Vector got_pos = first_tap(near_pos);
  HS_EXPECT_NEAR(got_pos.x, exp_pos.x, 1e-3f);
  HS_EXPECT_NEAR(got_pos.y, exp_pos.y, 1e-3f);
  HS_EXPECT_NEAR(got_pos.z, exp_pos.z, 1e-3f);

  // Probe near -Y (projection ~ -1 -> t ~ 0 -> first index 0 = X rotation).
  math::Vector near_neg = math::Vector(0.15f, -0.98f, 0.0f).normalized();
  math::Vector exp_neg = math::rotate(near_neg, oris[0].get());
  math::Vector got_neg = first_tap(near_neg);
  HS_EXPECT_NEAR(got_neg.x, exp_neg.x, 1e-3f);
  HS_EXPECT_NEAR(got_neg.y, exp_neg.y, 1e-3f);

  // Disabled -> verbatim passthrough (no rotation).
  slice.set_enabled(false);
  math::Vector pass = first_tap(near_pos);
  HS_EXPECT_NEAR(pass.x, near_pos.x, 1e-6f);
  HS_EXPECT_NEAR(pass.y, near_pos.y, 1e-6f);
  HS_EXPECT_NEAR(pass.z, near_pos.z, 1e-6f);
}

/**
 * @brief Verifies OrientSlice's clip-cull bounds the edge over every candidate
 *        slice, short-circuits on the first hit, rotates a planar basis, passes
 *        through when disabled or empty, and keeps an edge the selected slice
 *        carries into the band.
 * @details The endpoints can fall in different slices, so the cull spans all
 *          candidates.
 */
inline void test_world_orient_slice_cull_edge_bounds_all_slices() {
  static_assert(Filter::has_cull_edge<Filter::World::OrientSlice>);
  math::Orientation<> oris[2];
  oris[0].set(math::make_rotation(math::X_AXIS,
                                  math::PI_F / 2)); // leaves +X where it is
  oris[1].set(math::make_rotation(math::Z_AXIS, math::PI_F / 2));
  std::span<const math::Orientation<>> span(oris, 2);
  Filter::World::OrientSlice slice(span, math::Y_AXIS);

  math::Vector seen[2];
  int n = 0;
  bool hit = slice.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &, const math::Basis *) {
        if (n < 2)
          seen[n] = a;
        ++n;
        return false;
      });
  HS_EXPECT_FALSE(hit);
  // One single-frame tween step per candidate slice.
  HS_EXPECT_EQ(n, 2);
  for (int i = 0; i < 2; ++i) {
    math::Vector expected = math::rotate(math::X_AXIS, oris[i].get());
    HS_EXPECT_NEAR(seen[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(seen[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(seen[i].z, expected.z, 1e-4f);
  }

  n = 0;
  hit = slice.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &, const math::Vector &, const math::Basis *) {
        ++n;
        return true;
      });
  HS_EXPECT_TRUE(hit);
  HS_EXPECT_EQ(n, 1);

  // A planar edge rotates its basis alongside the endpoints, per candidate.
  const math::Basis pb{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  math::Vector normals[2];
  n = 0;
  slice.cull_edge(
      math::X_AXIS, math::X_AXIS, &pb,
      [&](const math::Vector &, const math::Vector &, const math::Basis *bp) {
        if (bp && n < 2)
          normals[n] = bp->v;
        ++n;
        return false;
      });
  HS_EXPECT_EQ(n, 2);
  for (int i = 0; i < 2; ++i) {
    math::Vector expected = math::rotate(math::Y_AXIS, oris[i].get());
    HS_EXPECT_NEAR(normals[i].x, expected.x, 1e-4f);
    HS_EXPECT_NEAR(normals[i].y, expected.y, 1e-4f);
    HS_EXPECT_NEAR(normals[i].z, expected.z, 1e-4f);
  }

  // A probe near +Y selects the Z rotation, which swings it onto the |x| band:
  // the band predicate rejects the source endpoints, plot() lands the point
  // inside, and the cull must therefore keep the edge.
  const math::Vector probe = math::Vector(0.15f, 0.98f, 0.0f).normalized();
  auto in_x_band = [](const math::Vector &a, const math::Vector &b,
                      const math::Basis *) {
    return std::fabs(a.x) > 0.9f && std::fabs(b.x) > 0.9f;
  };
  HS_EXPECT_FALSE(in_x_band(probe, probe, nullptr));
  math::Vector plotted{};
  slice.plot(
      probe, Pixel(1, 1, 1), 0.0f, 1.0f,
      [&](const math::Vector &v, const Pixel &, float, float) { plotted = v; });
  HS_EXPECT_GT(std::fabs(plotted.x), 0.9f);
  HS_EXPECT_TRUE(slice.cull_edge(probe, probe, nullptr, in_x_band));

  // Disabled -> the edge reaches the tail once, unrotated.
  math::Vector passed{};
  auto record_once = [&](const math::Vector &a, const math::Vector &,
                         const math::Basis *) {
    passed = a;
    ++n;
    return false;
  };
  slice.set_enabled(false);
  n = 0;
  slice.cull_edge(math::X_AXIS, math::X_AXIS, nullptr, record_once);
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(passed.x, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(passed.y, 0.0f, 1e-6f);

  // An empty candidate list is the same passthrough.
  Filter::World::OrientSlice empty(std::span<const math::Orientation<>>{},
                                   math::Y_AXIS);
  n = 0;
  passed = math::Vector();
  empty.cull_edge(math::X_AXIS, math::X_AXIS, nullptr, record_once);
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(passed.x, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(passed.z, 0.0f, 1e-6f);
}

/**
 * @brief Verifies VertexReplicate fans a point onto N copies via rotations from
 *        vertices[0] to each vertex, every copy carrying the same source age.
 * @details Replication is spatial, not temporal. Probing with vertices[0] maps
 *          copy i back to vertices[i] exactly.
 */
inline void test_world_vertex_replicate_fanout_and_age() {
  constexpr int N = 3;
  std::array<math::Vector, N> verts = {math::X_AXIS, math::Y_AXIS,
                                       math::Z_AXIS};
  Filter::World::VertexReplicate<N> vr(verts);

  Tap3D taps[N]{};
  int n = 0;
  vr.plot(math::X_AXIS, Pixel(1, 1, 1), 7.0f, 1.0f,
          [&](const math::Vector &v, const Pixel &c, float age, float a) {
            if (n < N)
              taps[n] = {v, c, age, a};
            ++n;
          });
  HS_EXPECT_EQ(n, N);
  // rotate(vertices[0], make_rotation(v0, v_i)) == v_i.
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_NEAR(taps[i].v.x, verts[i].x, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.y, verts[i].y, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.z, verts[i].z, 1e-4f);
    // Age is unchanged for every copy (no per-copy +i offset).
    HS_EXPECT_NEAR(taps[i].age, 7.0f, 1e-6f);
  }

  std::array<math::Vector, N> updated = {math::X_AXIS, -math::Y_AXIS,
                                         -math::Z_AXIS};
  vr.set_vertices(updated);
  n = 0;
  vr.plot(math::X_AXIS, Pixel(1, 1, 1), 7.0f, 1.0f,
          [&](const math::Vector &v, const Pixel &c, float age, float a) {
            if (n < N)
              taps[n] = {v, c, age, a};
            ++n;
          });
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_NEAR(taps[i].v.x, updated[i].x, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.y, updated[i].y, 1e-4f);
    HS_EXPECT_NEAR(taps[i].v.z, updated[i].z, 1e-4f);
  }
}

/**
 * @brief Verifies VertexReplicate's clip-cull re-emits the edge under the same
 *        rotations plot() applies, and short-circuits on the first hit.
 */
inline void test_world_vertex_replicate_cull_edge_mirrors_plot() {
  constexpr int N = 3;
  static_assert(Filter::has_cull_edge<Filter::World::VertexReplicate<N>>);
  std::array<math::Vector, N> verts = {math::X_AXIS, math::Y_AXIS,
                                       math::Z_AXIS};
  const Filter::World::VertexReplicate<N> vr(verts);

  math::Vector seen[N];
  int n = 0;
  bool hit = vr.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &, const math::Basis *) {
        if (n < N)
          seen[n] = a;
        ++n;
        return false;
      });
  HS_EXPECT_FALSE(hit);
  HS_EXPECT_EQ(n, N);
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_NEAR(seen[i].x, verts[i].x, 1e-4f);
    HS_EXPECT_NEAR(seen[i].y, verts[i].y, 1e-4f);
    HS_EXPECT_NEAR(seen[i].z, verts[i].z, 1e-4f);
  }

  n = 0;
  hit = vr.cull_edge(
      math::X_AXIS, math::X_AXIS, nullptr,
      [&](const math::Vector &, const math::Vector &, const math::Basis *) {
        ++n;
        return true;
      });
  HS_EXPECT_TRUE(hit);
  HS_EXPECT_EQ(n, 1);

  // A planar edge rotates its basis alongside the endpoints.
  const math::Basis pb{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  math::Vector normals[N];
  n = 0;
  vr.cull_edge(
      math::X_AXIS, math::X_AXIS, &pb,
      [&](const math::Vector &, const math::Vector &, const math::Basis *bp) {
        if (bp && n < N)
          normals[n] = bp->u;
        ++n;
        return false;
      });
  HS_EXPECT_EQ(n, N);
  for (int i = 0; i < N; ++i) {
    HS_EXPECT_NEAR(normals[i].x, verts[i].x, 1e-4f);
    HS_EXPECT_NEAR(normals[i].y, verts[i].y, 1e-4f);
    HS_EXPECT_NEAR(normals[i].z, verts[i].z, 1e-4f);
  }
}

/**
 * @brief Verifies Pipeline::could_intersect_clip walks the whole stage chain:
 *        each world stage's cull_edge feeds the next, an identity stage without
 *        one forwards the edge unchanged, and the sink runs the predicate.
 * @details The tail stage sees the head's rotated copies.
 */
inline void test_pipeline_could_intersect_clip_forwards_through_stages() {
  constexpr int W = 32, H = 16;
  const math::Quaternion q = math::make_rotation(math::X_AXIS, math::PI_F / 2);
  math::Orientation<> ori(q); // one frame -> a single tween step
  Pipeline<W, H, Filter::World::Orient, Filter::World::Replicate<W>> pipe(ori,
                                                                          2);

  // Head image, then the tail's second Y-axis copy of that image.
  const math::Vector head_a = math::rotate(math::Y_AXIS, q);
  const math::Vector head_b = math::rotate(math::Z_AXIS, q);
  const math::Quaternion step = math::make_rotation(math::Y_AXIS, math::PI_F);
  const math::Vector tail_a = math::rotate(head_a, step).normalized();
  const math::Vector tail_b = math::rotate(head_b, step).normalized();

  math::Vector seen_a[2], seen_b[2];
  int n = 0;
  bool hit = pipe.could_intersect_clip(
      math::Y_AXIS, math::Z_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &b, const math::Basis *) {
        if (n < 2) {
          seen_a[n] = a;
          seen_b[n] = b;
        }
        ++n;
        return false;
      });
  HS_EXPECT_FALSE(hit);
  // Orient's one tween step x Replicate's two copies.
  HS_EXPECT_EQ(n, 2);
  HS_EXPECT_NEAR(seen_a[0].x, head_a.x, 1e-4f);
  HS_EXPECT_NEAR(seen_a[0].y, head_a.y, 1e-4f);
  HS_EXPECT_NEAR(seen_a[0].z, head_a.z, 1e-4f);
  HS_EXPECT_NEAR(seen_b[0].x, head_b.x, 1e-4f);
  HS_EXPECT_NEAR(seen_b[0].y, head_b.y, 1e-4f);
  HS_EXPECT_NEAR(seen_b[0].z, head_b.z, 1e-4f);
  HS_EXPECT_NEAR(seen_a[1].x, tail_a.x, 1e-4f);
  HS_EXPECT_NEAR(seen_a[1].y, tail_a.y, 1e-4f);
  HS_EXPECT_NEAR(seen_a[1].z, tail_a.z, 1e-4f);
  HS_EXPECT_NEAR(seen_b[1].x, tail_b.x, 1e-4f);
  HS_EXPECT_NEAR(seen_b[1].y, tail_b.y, 1e-4f);
  HS_EXPECT_NEAR(seen_b[1].z, tail_b.z, 1e-4f);

  // Only the composed transform reaches this target: neither the source
  // geometry nor the head-only image does.
  auto near_target = [&](const math::Vector &a, const math::Vector &,
                         const math::Basis *) {
    return math::distance_between(a, tail_a) < 1e-3f;
  };
  HS_EXPECT_FALSE(near_target(math::Y_AXIS, math::Z_AXIS, nullptr));
  HS_EXPECT_FALSE(near_target(head_a, head_b, nullptr));
  HS_EXPECT_TRUE(pipe.could_intersect_clip(math::Y_AXIS, math::Z_AXIS, nullptr,
                                           near_target));

  // First hit short-circuits the whole chain.
  n = 0;
  hit = pipe.could_intersect_clip(
      math::Y_AXIS, math::Z_AXIS, nullptr,
      [&](const math::Vector &, const math::Vector &, const math::Basis *) {
        ++n;
        return true;
      });
  HS_EXPECT_TRUE(hit);
  HS_EXPECT_EQ(n, 1);

  // An identity stage without cull_edge forwards the edge unchanged, so the
  // head's rotated copy still reaches the sink.
  Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>> plain(
      ori);
  math::Vector plain_a{};
  n = 0;
  plain.could_intersect_clip(
      math::Y_AXIS, math::Z_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &, const math::Basis *) {
        plain_a = a;
        ++n;
        return false;
      });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(plain_a.x, head_a.x, 1e-4f);
  HS_EXPECT_NEAR(plain_a.y, head_a.y, 1e-4f);
  HS_EXPECT_NEAR(plain_a.z, head_a.z, 1e-4f);

  // The filter-free sink is the terminal: it runs the predicate on the edge it
  // was handed and answers with it.
  Pipeline<W, H> sink;
  math::Vector sink_a{};
  n = 0;
  hit = sink.could_intersect_clip(
      math::Y_AXIS, math::Z_AXIS, nullptr,
      [&](const math::Vector &a, const math::Vector &, const math::Basis *) {
        sink_a = a;
        ++n;
        return true;
      });
  HS_EXPECT_TRUE(hit);
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(sink_a.y, 1.0f, 1e-6f);

  // The planar basis rides the same chain.
  const math::Basis pb{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  math::Vector tail_normal{};
  n = 0;
  pipe.could_intersect_clip(
      math::Y_AXIS, math::Z_AXIS, &pb,
      [&](const math::Vector &, const math::Vector &, const math::Basis *bp) {
        if (bp)
          tail_normal = bp->v;
        ++n;
        return false;
      });
  HS_EXPECT_EQ(n, 2);
  const math::Vector expected_normal =
      math::rotate(math::rotate(math::Y_AXIS, q), step);
  HS_EXPECT_NEAR(tail_normal.x, expected_normal.x, 1e-4f);
  HS_EXPECT_NEAR(tail_normal.y, expected_normal.y, 1e-4f);
  HS_EXPECT_NEAR(tail_normal.z, expected_normal.z, 1e-4f);
}

/**
 * @brief Verifies Mobius warps via stereographic -> Mobius -> inverse
 *        stereographic, with the default identity map and a non-identity map.
 * @details The default MobiusParams is the identity (a=1,b=0,c=0,d=1), so an
 *          interior point round-trips back to itself; a non-identity map
 *          actually moves it.
 */
inline void test_world_mobius_identity_and_transform() {
  const math::MobiusParams identity; // a=1,b=0,c=0,d=1
  Filter::World::Mobius mob(identity);

  const math::Vector v = math::Vector(0.4f, 0.3f, 0.86f).normalized();
  math::Vector out{};
  int n = 0;
  mob.plot(v, Pixel(1, 2, 3), 2.0f, 0.5f,
           [&](const math::Vector &o, const Pixel &, float, float) {
             out = o;
             ++n;
           });
  HS_EXPECT_EQ(n, 1);
  HS_EXPECT_NEAR(out.x, v.x, 1e-3f);
  HS_EXPECT_NEAR(out.y, v.y, 1e-3f);
  HS_EXPECT_NEAR(out.z, v.z, 1e-3f);

  // A translation map f(z) = z + 1 moves the point and keeps it on the sphere.
  math::MobiusParams shift(1, 0, 1, 0, 0, 0, 1, 0); // a=1, b=1, c=0, d=1
  Filter::World::Mobius mob2(shift);
  math::Vector out2{};
  mob2.plot(
      v, Pixel(1, 1, 1), 0.0f, 1.0f,
      [&](const math::Vector &o, const Pixel &, float, float) { out2 = o; });
  HS_EXPECT_NEAR(out2.length(), 1.0f, 1e-3f);
  HS_EXPECT_GT(math::distance_between(out2, v), 0.05f);
  const math::Vector want = math::mobius_transform(v, shift);
  HS_EXPECT_NEAR(out2.x, want.x, 1e-6f);
  HS_EXPECT_NEAR(out2.y, want.y, 1e-6f);
  HS_EXPECT_NEAR(out2.z, want.z, 1e-6f);

  // The filter reads its bound parameters live.
  shift.b.re = 2.0f;
  math::Vector out3{};
  mob2.plot(
      v, Pixel(1, 1, 1), 0.0f, 1.0f,
      [&](const math::Vector &o, const Pixel &, float, float) { out3 = o; });
  const math::Vector want3 = math::mobius_transform(v, shift);
  HS_EXPECT_GT(math::distance_between(want3, want), 0.05f);
  HS_EXPECT_NEAR(out3.x, want3.x, 1e-6f);
  HS_EXPECT_NEAR(out3.y, want3.y, 1e-6f);
  HS_EXPECT_NEAR(out3.z, want3.z, 1e-6f);
  shift.b.re = 1.0f;

  constexpr int W = 32, H = 16;
  std::array<Pixel, W * H> expected;
  {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipe;
    {
      Canvas canvas(fx);
      pipe.plot(canvas, out2, Pixel(40000, 20000, 10000), 0.0f, 0.7f);
    }
    fx.advance_display();
    HS_EXPECT_GT(count_lit_canvas(fx), size_t{0});
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        expected[y * W + x] = fx.get_pixel(x, y);
  }
  {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H, Filter::World::Mobius, Filter::Screen::AntiAlias<W, H>> pipe(
        mob2);
    {
      Canvas canvas(fx);
      pipe.plot(canvas, v, Pixel(40000, 20000, 10000), 0.0f, 0.7f);
    }
    fx.advance_display();
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const Pixel &pixel = expected[y * W + x];
        HS_EXPECT_PIXEL(fx.get_pixel(x, y), pixel.r, pixel.g, pixel.b);
      }
  }

  // Mobius offers no cull_edge bound, so the pipeline fold forces a
  // full-canvas render.
  static_assert(!Filter::has_cull_edge<Filter::World::Mobius>);
  HS_EXPECT_TRUE(Filter::World::Mobius::crosses_segments);
  HS_EXPECT_TRUE(
      (Pipeline<17, 9, Filter::World::Mobius>::any_crosses_segments));
}
