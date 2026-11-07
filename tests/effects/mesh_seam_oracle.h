/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Mesh seam geometry and colour oracles. */

/** @brief Adjacent-column colour steps at the arm boundaries and mid-arm. */
struct MeshSeamSteps {
  double mean_arm_ratios = 0.0; ///< Sum of mean arm-boundary ratios.
  double arm_ratios = 0.0; ///< Sum of the larger arm-boundary ratio per frame.
  double interior_ratios = 0.0; ///< Sum of the mid-arm ratios.
  int frames = 0;               ///< Consecutive frame pairs attempted.
  int arm_samples = 0; ///< Frame pairs with both arm boundaries measurable.
  int interior_samples = 0; ///< Frame pairs with measurable mid-arm neighbors.

  /**
   * @brief Accumulates one pair of consecutive frames.
   * @param previous Previous displayed frame.
   * @param current Current displayed frame.
   * @param width Frame width.
   * @param height Frame height.
   */
  void add(const std::vector<Pixel> &previous,
           const std::vector<Pixel> &current, int width, int height) {
    const auto column_difference = [&](const auto &a, int x, const auto &b,
                                       int z) {
      double sum = 0.0;
      for (int y = 0; y < height; ++y) {
        const Pixel &p = a[size_t(y) * width + (x + width) % width];
        const Pixel &q = b[size_t(y) * width + (z + width) % width];
        sum += std::abs(int(p.r) - q.r) + std::abs(int(p.g) - q.g) +
               std::abs(int(p.b) - q.b);
      }
      return sum;
    };
    double arm = 0.0, mean_arm = 0.0, interior = 0.0, interior_neighbors = 0.0;
    int arm_boundaries = 0;
    for (int boundary : {0, width / 2}) {
      double neighbors = 0.0;
      for (int j = 1; j <= 4; ++j) {
        neighbors += column_difference(previous, boundary - j - 1, previous,
                                       boundary - j);
        neighbors +=
            column_difference(current, boundary + j - 1, current, boundary + j);
      }
      const double seam =
          column_difference(previous, boundary - 1, current, boundary);
      if (neighbors > 0.0) {
        const double ratio = seam / (neighbors / 8.0);
        arm = std::max(arm, ratio);
        mean_arm += ratio;
        ++arm_boundaries;
      }
      const int middle = boundary + width / 4;
      interior += column_difference(current, middle - 1, current, middle);
      for (int j : {-3, -2, 1, 2})
        interior_neighbors +=
            column_difference(current, middle - 1 + j, current, middle + j);
    }
    if (arm_boundaries == 2) {
      arm_ratios += arm;
      mean_arm_ratios += mean_arm * 0.5;
      ++arm_samples;
    }
    if (interior_neighbors > 0.0) {
      interior_ratios += interior / (interior_neighbors / 4.0);
      ++interior_samples;
    }
    ++frames;
  }
};

/**
 * @brief Checks boundary steps against their neighboring steps.
 * @param steps Accumulated frame differences.
 * @param label Effect name for the result.
 * @param mean_bound Bound on the mean arm-boundary ratio.
 * @param arm_bound Bound on the mean of each frame's larger arm ratio.
 * @param interior_bound Bound on the mean mid-arm ratio.
 */
inline void expect_mesh_seam_continuous(const MeshSeamSteps &steps,
                                        const char *label,
                                        double mean_bound = 1.06,
                                        double arm_bound = 1.15,
                                        double interior_bound = 1.06) {
  HS_CONTEXT(label);
  HS_EXPECT_GT(steps.frames, 0);
  HS_EXPECT_GE(steps.arm_samples * 4, steps.frames * 3);
  HS_EXPECT_GE(steps.interior_samples * 4, steps.frames * 3);
  HS_EXPECT_GT(steps.arm_samples, 0);
  HS_EXPECT_GT(steps.interior_samples, 0);
  if (steps.arm_samples == 0 || steps.interior_samples == 0)
    return;
  std::printf(
      "  [%s] mean seam step %.3f x neighbors; mid-arm %.3f; samples %d/%d and %d/%d\n",
      label, steps.arm_ratios / steps.arm_samples,
      steps.interior_ratios / steps.interior_samples, steps.arm_samples,
      steps.frames, steps.interior_samples, steps.frames);
  std::printf("  mean across arms %.3f\n",
              steps.mean_arm_ratios / steps.arm_samples);
  HS_EXPECT_LE(steps.mean_arm_ratios / steps.arm_samples, mean_bound);
  HS_EXPECT_LE(steps.arm_ratios / steps.arm_samples, arm_bound);
  HS_EXPECT_LE(steps.interior_ratios / steps.interior_samples, interior_bound);
}

/** @brief One arm side's copy of a face in a mesh deformation scan. */
struct FaceCopy {
  std::vector<uint16_t> indices;     ///< The face's vertex indices.
  std::vector<math::Vector> polygon; ///< Its scanned vertices, model space.
  float size = 0.0f;                 ///< Inradius the scan shaded it with.
};

/** @brief A deformation scan's face copies: first arm, then second, by face. */
using FaceCopies = std::array<std::map<size_t, FaceCopy>, 2>;

/** @brief Measures successive frame seam steps for a pinned effect playback.
 * @tparam E Effect template.
 * @tparam W,H Display dimensions.
 * @param frames Number of frames to render.
 * @return Accumulated seam and interior steps. */
template <template <int, int> class E, int W, int H>
inline MeshSeamSteps mesh_seam_steps(int frames) {
  reset_globals();
  pin_frame_clock(0);
  E<W, H> effect;
  effect.init();
  std::vector<Pixel> previous, current;
  MeshSeamSteps steps;
  for (int f = 0; f < frames; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
    capture_frame<W, H>(effect, current);
    if (f > 0)
      steps.add(previous, current, W, H);
    previous.swap(current);
  }
  hs::clear_mock_time();
  return steps;
}

/// Largest extrapolation of a nonpolar face's contiguous capture branch.
inline constexpr float SEAM_EXTRAPOLATION = 1.0f;

/**
 * @brief Records the faces a deformation scan draws on one arm.
 * @tparam CAP Orientation history capacity.
 * @param canvas Canvas inside the scan.
 * @param mesh The scanned mesh, at the scan's captured vertices.
 * @param orientation Camera at the current endpoint.
 * @param copies Receives the arm's copies.
 * @return The arm, or -1 for a scan its render-x window does not hold to
 *         one arm.
 */
template <int CAP>
inline int record_face_copies(const Canvas &canvas, const MeshState &mesh,
                              const math::Orientation<CAP> &orientation,
                              FaceCopies &copies) {
  const Spacetime::TimeWindow &window = canvas.time_window();
  const Scan::ColumnBasis *time = canvas.column_basis();
  const int arm = window.arm(canvas.width());
  if (!window.active || arm < 0 || !time || !time->mesh.selected_faces)
    return -1;
  for (size_t f = 0; f < mesh.num_faces(); ++f) {
    if (time->mesh.selected_faces[f] != Scan::detail::FacePass::ARM_COPY)
      continue;
    FaceCopy &copy = copies[arm][f];
    copy = FaceCopy{};
    if (const Scan::FaceArmCopy *own = time->mesh.face_copy) {
      for (size_t k = 0; k < own->count; ++k) {
        copy.indices.push_back(own->sources[k]);
        copy.polygon.push_back(orientation.unorient(own->vertices[k]));
      }
      continue;
    }
    const uint16_t *face =
        mesh.get_faces_data() + mesh.get_face_offsets_data()[f];
    for (size_t k = 0; k < mesh.get_face_counts_data()[f]; ++k) {
      copy.indices.push_back(face[k]);
      copy.polygon.push_back(orientation.unorient(mesh.vertices[face[k]]));
    }
  }
  return arm;
}

/**
 * @brief Renders face copies as a still mesh through a camera.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam SetupFn Callable (size_t face, float size, FacePaletteShader &)
 *         shading a face as its effect does.
 * @param effect Effect whose display receives the render.
 * @param copies Face copies, by face.
 * @param faces The faces to draw.
 * @param camera Camera the copies are drawn through.
 * @param setup Shades each drawn face.
 * @param ids Whether to draw each face flat in its position in @p faces plus
 *        one, in the red channel, instead of shading it.
 * @return Displayed frame.
 */
template <int W, int H, typename EffectT, typename SetupFn>
inline std::vector<Pixel>
render_face_copies(EffectT &effect, const std::map<size_t, FaceCopy> &copies,
                   const std::vector<size_t> &faces,
                   const math::Quaternion &camera, SetupFn &&setup, bool ids) {
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(effect);
    ScratchScope guard(scratch_arena_a);
    size_t count = 0;
    for (size_t f : faces)
      count += copies.at(f).polygon.size();
    MeshState mesh;
    mesh.vertices.bind(scratch_arena_a, count);
    mesh.faces.bind(scratch_arena_a, count);
    mesh.face_counts.bind(scratch_arena_a, faces.size());
    mesh.face_offsets.bind(scratch_arena_a, faces.size());
    for (size_t f : faces) {
      const FaceCopy &copy = copies.at(f);
      mesh.face_offsets.push_back(static_cast<uint16_t>(mesh.faces.size()));
      mesh.face_counts.push_back(static_cast<uint8_t>(copy.polygon.size()));
      for (const math::Vector &v : copy.polygon) {
        mesh.faces.push_back(static_cast<uint16_t>(mesh.vertices.size()));
        mesh.vertices.push_back(math::rotate(v, camera));
      }
    }
    if (ids) {
      const auto id = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(static_cast<uint16_t>(f.v2 + 1.5f), 0, 0), 1);
      };
      Scan::Mesh::draw<W, H>(pipeline, canvas, mesh, id, scratch_arena_b);
    } else {
      FacePaletteShader shader;
      auto shade = [&](size_t j, float size) { setup(faces[j], size, shader); };
      Scan::Mesh::draw_specialized<W, H>(pipeline, canvas, mesh, shader,
                                         scratch_arena_b, nullptr, shade);
    }
  }
  effect.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(effect, out);
  return out;
}

/**
 * @brief Largest colour change a face's shading allows between two copies.
 * @tparam SetupFn As render_face_copies() takes it.
 * @param face Face index.
 * @param size The face's inradius in one copy.
 * @param shift Bound on how far the copies' inradius, and their edge depth
 *        at a pixel, differ.
 * @param setup Shades the face as its effect does.
 * @return 16-bit channel units: the shading at each edge depth against the
 *         shading within @p shift of that depth and of @p size.
 */
template <typename SetupFn>
inline int face_colour_bound(size_t face, float size, float shift,
                             SetupFn &&setup) {
  constexpr int DEPTHS = 64;
  constexpr float STEPS[] = {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f};
  FacePaletteShader one;
  setup(face, size, one);
  int bound = 0;
  for (float r : STEPS) {
    FacePaletteShader other;
    setup(face, fmaxf(size + r * shift, 0.0f), other);
    for (int k = 0; k <= DEPTHS; ++k) {
      Fragment a;
      a.v1 = -2.0f * size * k / DEPTHS;
      one(a.pos, a);
      for (float s : STEPS) {
        Fragment b;
        b.v1 = a.v1 + s * shift;
        other(b.pos, b);
        const Pixel &p = a.color.color, &q = b.color.color;
        bound = std::max<int>(bound, max_channel_gap(p, q));
      }
    }
  }
  return bound;
}

/// Columns either side of each arm boundary whose colours a seam check reads.
inline constexpr int SEAM_CHECK_COLUMNS = 3;

/**
 * @brief Checks that one frame's copies before an arm boundary and the next
 *        frame's copies after it are one polygon there that shades alike.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam ToleranceFn Callable (uint16_t vertex) returning the bound on that
 *         vertex's gap between the copies (radians of chord).
 * @tparam SetupFn As render_face_copies() takes it.
 * @param effect Effect whose display renders the copies.
 * @param label Draw named in failure context.
 * @param seam The boundary's first column after it: 0 or W / 2.
 * @param late Frame N's copies on the arm the boundary ends.
 * @param early Frame N+1's copies on the arm the boundary starts.
 * @param camera Camera the copies' colours are compared through.
 * @param tolerance Per-vertex gap bound.
 * @param setup Shades each face as its effect does.
 * @return The number of faces with interior colour samples compared.
 * @details Checks vertex gaps for all matching faces. A face's edge
 * depth and inradius move by at most twice its largest vertex gap, the
 * gnomonic stretch of a face within a quarter turn of its centre. Colours
 * are compared on the boundary columns' pixels a pixel clear of every edge.
 * Faces around a pole draw once and leave no copies.
 */
template <int W, int H, typename EffectT, typename ToleranceFn,
          typename SetupFn>
inline size_t
expect_seam_copies_match(EffectT &effect, const char *label, int seam,
                         const std::map<size_t, FaceCopy> &late,
                         const std::map<size_t, FaceCopy> &early,
                         const math::Quaternion &camera,
                         ToleranceFn &&tolerance, SetupFn &&setup) {
  HS_CONTEXT(label);
  HS_CONTEXT("seam column", seam);
  constexpr float PIXEL = math::TWO_PI_F / W;
  // Sign of x along the boundary's meridian.
  const float side = seam == 0 ? 1.0f : -1.0f;
  std::vector<size_t> faces;
  std::vector<int> colour_bounds;
  float worst_gap = 0.0f, worst_ratio = 0.0f, worst_size = 0.0f;
  for (const auto &[f, before] : late) {
    const auto after = early.find(f);
    if (after == early.end() || before.size <= 0.0f ||
        after->second.size <= 0.0f)
      continue;
    const auto seen = [&](const FaceCopy &copy) {
      std::vector<math::Vector> out;
      for (const auto &v : copy.polygon)
        out.push_back(math::rotate(v, camera));
      return out;
    };
    const std::vector<math::Vector> a = seen(before), b = seen(after->second);
    // A face meets the boundary when an edge crosses its meridian, or when
    // it lies within the checked columns on the boundary's side.
    const auto touches_seam = [&](const std::vector<math::Vector> &polygon) {
      float lo = 2.0f, hi = -2.0f, middle = 0.0f;
      for (size_t k = 0; k < polygon.size(); ++k) {
        const math::Vector &p = polygon[k];
        const math::Vector &q = polygon[(k + 1) % polygon.size()];
        lo = fminf(lo, p.z);
        hi = fmaxf(hi, p.z);
        middle += p.x;
        if (p.z * q.z <= 0.0f && p.z != q.z &&
            (p.x + (q.x - p.x) * p.z / (p.z - q.z)) * side > 0.0f)
          return true;
      }
      const float fringe = sinf(SEAM_CHECK_COLUMNS * PIXEL);
      return lo <= fringe && hi >= -fringe && middle * side > 0.0f;
    };
    if (!touches_seam(a) && !touches_seam(b))
      continue;
    // A face copy starts at its first vertex on its arm.
    const size_t n = before.indices.size();
    size_t offset = n;
    for (size_t o = 0; o < n && offset == n; ++o) {
      bool same = after->second.indices.size() == n;
      for (size_t k = 0; k < n && same; ++k)
        same = before.indices[k] == after->second.indices[(k + o) % n];
      if (same)
        offset = o;
    }
    HS_EXPECT_LT(offset, n);
    if (offset == n)
      continue;
    float face_tolerance = 0.0f;
    for (size_t k = 0; k < n; ++k) {
      const size_t j = (k + offset) % n;
      const float gap = sqrtf(
          math::distance_squared(before.polygon[k], after->second.polygon[j]));
      const float bound = tolerance(before.indices[k]);
      HS_EXPECT_LE(gap, bound);
      worst_gap = fmaxf(worst_gap, gap);
      worst_ratio = fmaxf(worst_ratio, gap / bound);
      face_tolerance = fmaxf(face_tolerance, bound);
    }
    const float size_gap = fabsf(before.size - after->second.size);
    HS_EXPECT_LE(size_gap, 2.0f * face_tolerance);
    worst_size = fmaxf(worst_size, size_gap / before.size);
    faces.push_back(f);
    colour_bounds.push_back(
        face_colour_bound(f, before.size, 2.0f * face_tolerance, setup));
  }
  const auto render = [&](const std::map<size_t, FaceCopy> &copies, bool ids) {
    return render_face_copies<W, H>(effect, copies, faces, camera, setup, ids);
  };
  const std::vector<Pixel> left = render(late, false);
  const std::vector<Pixel> right = render(early, false);
  const std::vector<Pixel> left_id = render(late, true);
  const std::vector<Pixel> right_id = render(early, true);
  const auto face_at = [&](int x, int y) {
    const int id = left_id[static_cast<size_t>(y) * W + x].r;
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx) {
        const int row = y + dy;
        if (row < 0 || row >= H)
          return 0;
        const size_t j = static_cast<size_t>(row) * W + (x + dx + W) % W;
        if (left_id[j].r != id || right_id[j].r != id)
          return 0;
      }
    return id;
  };
  int colour = 0, excess = 0, interior = 0;
  std::vector<bool> compared(faces.size() + 1, false);
  for (int y = 0; y < H; ++y)
    for (int c = -SEAM_CHECK_COLUMNS; c < SEAM_CHECK_COLUMNS; ++c) {
      const int x = (seam + c + W) % W;
      const int id = face_at(x, y);
      if (id == 0)
        continue;
      HS_EXPECT_LT(static_cast<size_t>(id), compared.size());
      if (static_cast<size_t>(id) >= compared.size())
        continue;
      compared[id] = true;
      ++interior;
      const size_t i = static_cast<size_t>(y) * W + x;
      const int difference = max_channel_gap(left[i], right[i]);
      colour = std::max(colour, difference);
      excess = std::max(excess, difference - colour_bounds[id - 1]);
    }
  const size_t compared_faces =
      static_cast<size_t>(std::count(compared.begin(), compared.end(), true));
  std::printf(
      "  [%s, seam %d] %zu colour faces, %zu geometry faces: vertex gap %.4f px "
      "(%.3f of bound), inradius %.5f, colour %.3f steps over %d "
      "interior pixels\n",
      label, seam, compared_faces, faces.size(), worst_gap / PIXEL, worst_ratio,
      worst_size, colour / 257.0, interior);
  HS_EXPECT_LE(excess, 0);
  return compared_faces;
}

/**
 * @brief expect_seam_copies_match() at both arm boundaries.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam EffectT Effect owning the display.
 * @tparam ToleranceFn As expect_seam_copies_match() takes it.
 * @tparam SetupFn As render_face_copies() takes it.
 * @param effect Effect whose display renders the copies.
 * @param label Draw named in failure context.
 * @param frames Frame N's and frame N+1's copies.
 * @param camera Camera at the frames' shared endpoint.
 * @param tolerance Per-vertex gap bound.
 * @param setup Shades each face as its effect does.
 * @return Faces with colour-compared interior pixels, over both boundaries.
 * @details Vertex-gap checks cover every recorded face.
 */
template <int W, int H, typename EffectT, typename ToleranceFn,
          typename SetupFn>
inline size_t expect_seams_match(EffectT &effect, const char *label,
                                 const std::array<FaceCopies, 2> &frames,
                                 const math::Quaternion &camera,
                                 ToleranceFn &&tolerance, SetupFn &&setup) {
  return expect_seam_copies_match<W, H>(effect, label, W / 2, frames[0][0],
                                        frames[1][1], camera, tolerance,
                                        setup) +
         expect_seam_copies_match<W, H>(effect, label, 0, frames[0][1],
                                        frames[1][0], camera, tolerance, setup);
}

/**
 * @brief Rotation angle between two unit quaternions.
 * @param a First rotation.
 * @param b Second rotation.
 * @return Radians in [0, PI].
 */
inline float quaternion_angle(const math::Quaternion &a,
                              const math::Quaternion &b) {
  const float cosine = fabsf(math::dot(a, b));
  return 2.0f * acosf(fminf(1.0f, cosine));
}

/**
 * @brief Bound on how far two copies' capture times of one vertex differ.
 * @param motion The vertex's motion over either frame (radians).
 * @param turn Largest camera rotation over either frame (radians).
 * @param horizontal The vertex's distance from the polar axis.
 * @return Frame fraction: math::fast_atan2's peak error and one correction
 *         step's residual in each copy, and the column basis the
 *         extrapolating copy holds at its arm's end column.
 */
inline float seam_capture_slack(float motion, float turn, float horizontal) {
  constexpr float FAST_ATAN2_ERROR =
      static_cast<float>(math::FAST_ATAN2_MAX_ERROR);
  const float ratio = (motion + turn) / (math::PI_F * horizontal);
  return 2.0f * FAST_ATAN2_ERROR / math::PI_F +
         2.0f * (1.0f + SEAM_EXTRAPOLATION) * ratio * ratio +
         SEAM_EXTRAPOLATION * turn / (math::PI_F * horizontal);
}

/**
 * @brief Bound on a vertex's gap between frame N's late copy and frame N+1's
 *        early copy.
 * @tparam W Canvas width.
 * @tparam LateFn Callable (float u) giving the vertex on frame N's path.
 * @tparam EarlyFn Callable (float u) giving the vertex on frame N+1's path.
 * @param late Frame N's path, model space.
 * @param early Frame N+1's path, model space.
 * @param camera Camera at the frames' shared endpoint.
 * @param turn Largest camera rotation over either frame (radians).
 * @return The paths' largest gap within SEAM_EXTRAPOLATION of the shared
 *         endpoint, which is second order in the motion, plus the vertex's
 *         motion over seam_capture_slack(), one column interval, and any
 *         endpoint-held motion.
 */
template <int W, typename LateFn, typename EarlyFn>
inline float seam_vertex_tolerance(LateFn &&late, EarlyFn &&early,
                                   const math::Quaternion &camera, float turn) {
  constexpr int SAMPLES = 8;
  constexpr float FLOAT_ROUNDING = 1e-6f;
  const auto gap = [](const math::Vector &a, const math::Vector &b) {
    return sqrtf(math::distance_squared(a, b));
  };
  float paths = 0.0f;
  for (int s = -SAMPLES; s <= SAMPLES; ++s) {
    const float d = SEAM_EXTRAPOLATION * s / SAMPLES;
    paths = fmaxf(paths, gap(late(1.0f + d), early(d)));
  }
  const float before = gap(late(0.0f), late(1.0f));
  const float after = gap(early(0.0f), early(1.0f));
  const float motion = fmaxf(before, after);
  const bool held = fminf(before, after) <= Scan::HELD_VERTEX_ARC<W>;
  const float held_slack = held ? motion : 0.0f;
  const math::Vector seen = math::rotate(late(1.0f), camera);
  const float horizontal = sqrtf(seen.x * seen.x + seen.z * seen.z);
  return held_slack + paths +
         (seam_capture_slack(motion, turn, horizontal) + 2.0f / W) * motion +
         FLOAT_ROUNDING;
}
