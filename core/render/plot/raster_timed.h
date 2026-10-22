/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file raster_timed.h
 * @brief Geometry callbacks and coherent column-time prediction.
 */

/** @brief Borrowed geometric control-point evaluator and birth registers. */
struct GeometryCallback {
  /** Evaluates one original control point at a frame fraction. */
  FunctionRef<math::Vector(size_t, float)> evaluate;
  /** Uses normalized linear interpolation between lattice positions. */
  bool linear_pair = false;
  /** Optional birth fraction per original control point. */
  std::span<const float> births{};
  bool azimuth_bound = false; /**< Callback bound measures azimuth directly. */
  Temporal::PathBound path{}; /**< Derivative bounds of every control point's
                                  path over the frame. */

  bool column_edges =
      false; /**< Rasterizes geodesics at destination column times. */
  bool source_positions =
      true; /**< Shaders read source-space positions under world stages;
                 false hands them world positions. */
  /** Bounds each vertex's path, in place of `path`, when set. */
  FunctionRef<Temporal::PathBound(size_t)> vertex_path{};
  bool chained = false; /**< Edges run along polylines, each vertex starting
                             at most one edge and ending at most one other:
                             the column walk plots each edge's start but not
                             its end, so a closed chain plots every vertex once
                             and an open one all but its last. */

  /** @brief Constructs an exact evaluator.
   * @param fn Borrowed control-point callback.
   * @return The exact geometry descriptor. */
  static GeometryCallback exact(FunctionRef<math::Vector(size_t, float)> fn) {
    return {fn, false, {}};
  }
  /** @brief Constructs a two-position interpolating evaluator.
   * @param fn Borrowed callback, read only at fractions zero and one.
   * @return The interpolating geometry descriptor. */
  static GeometryCallback
  interpolating(FunctionRef<math::Vector(size_t, float)> fn) {
    return {fn, true, {}};
  }
};

/** @brief A control point derived from one source edge. */
struct TimedControl {
  uint32_t edge;            /**< Original edge index. */
  float fraction;           /**< Arc fraction along the original edge. */
  math::PixelCoords screen; /**< Unwrapped endpoint projection. */
  float motion;             /**< Column motion bound. */
};

/** @brief Interpolates one geodesic control point.
 * @param a,b Endpoints.
 * @param s Arc fraction.
 * @return Unit point on the edge. */
HS_HOT_INLINE inline math::Vector
timed_edge_point(const math::Vector &a, const math::Vector &b, float s) {
  if (s <= 0.0f)
    return a;
  if (s >= 1.0f)
    return b;
  const auto span = make_geodesic_edge_span(a, b);
  if (!span.have_axis)
    return (a + (b - a) * s).normalized();
  return GeodesicEdgeSampler{a, math::cross(span.axis, a), span.total}
      .pos(s)
      .normalized();
}

/** @brief Lattice endpoint retained across adjacent source edges. */
struct TimedLatticeCarry {
  size_t vertex =
      std::numeric_limits<size_t>::max(); /**< Retained unwrapped vertex. */
  math::Vector position; /**< The preceding edge's final lattice point. */
  size_t origin =
      std::numeric_limits<size_t>::max(); /**< First vertex of a closed walk. */
  math::Vector
      origin_position; /**< Lattice value reused by the closing endpoint. */
};

/** @brief Scratch lattice for one source edge. */
class TimedGeometry {
public:
  /** @brief Borrows geometry and initializes the edge's lattice.
   * @param points Current original control points.
   * @param callback Geometry descriptor.
   * @param edge Unwrapped original edge index.
   * @param carry Lattice endpoint retained from the preceding edge. */
  TimedGeometry(const Fragments &points, const GeometryCallback &callback,
                size_t edge, TimedLatticeCarry &carry)
      : points(points), callback(callback), first(edge) {
    if (callback.evaluate && callback.linear_pair) {
      previous = scratch_arena_a.allocate_n<math::Vector>(2);
      for (size_t i = edge; i <= edge + 1; ++i) {
        if (i == carry.vertex) {
          previous[i - edge] = carry.position;
        } else if (carry.origin != std::numeric_limits<size_t>::max() &&
                   i % points.size() == carry.origin % points.size()) {
          previous[i - edge] = carry.origin_position;
        } else {
          Temporal::count_geometry<
              &Temporal::GeometryCounters::callback_evaluations>();
          previous[i - edge] = callback.evaluate(i % points.size(), 0.0f);
        }
        if (carry.origin == std::numeric_limits<size_t>::max()) {
          carry.origin = i;
          carry.origin_position = previous[i - edge];
        }
      }
      carry.vertex = edge + 1;
      carry.position = previous[1];
    }
  }

  /** @return Whether the callback bound measures azimuth directly. */
  bool azimuth_bound() const { return callback.azimuth_bound; }

  /** @brief Tests whether a cached linear edge is exactly stationary.
   * @param edge Original edge index.
   * @return True when its two lattice pairs are identical. */
  bool stationary(size_t edge) const {
    if (!callback.evaluate)
      return true;
    if (!callback.linear_pair)
      return false;
    for (size_t i = edge; i <= edge + 1; ++i) {
      const auto &a = previous[i - first];
      const auto &b = points[i % points.size()].pos;
      if (a.x != b.x || a.y != b.y || a.z != b.z)
        return false;
    }
    return true;
  }

  /** @brief Evaluates an original vertex.
   * @param i Unwrapped original vertex index.
   * @param t Frame fraction.
   * @return The vertex at that fraction. */
  HS_O3_FN HS_NOINLINE_NOCLONE math::Vector vertex(size_t i, float t) const {
    const auto &end = points[i % points.size()].pos;
    if (!callback.evaluate || t == 1.0f)
      return end;
    if (callback.linear_pair)
      return (previous[i - first] * (1.0f - t) + end * t).normalized();
    Temporal::count_geometry<
        &Temporal::GeometryCounters::callback_evaluations>();
    return callback.evaluate(i % points.size(), t);
  }

  /** @brief Evaluates an original or fractional control point.
   * @param p Control descriptor.
   * @param t Frame fraction.
   * @return Geometry at t. */
  HS_O3_FN HS_NOINLINE_NOCLONE math::Vector position(const TimedControl &p,
                                                     float t) const {
    if (p.fraction == 0.0f)
      return vertex(p.edge, t);
    if (p.fraction == 1.0f)
      return vertex(p.edge + 1, t);
    return timed_edge_point(vertex(p.edge, t), vertex(p.edge + 1, t),
                            p.fraction);
  }

  /** @brief Interpolates source shading registers.
   * @param p Control descriptor.
   * @return Source registers at the control point. */
  HS_HOT_INLINE Fragment registers(const TimedControl &p) const {
    return Fragment::lerp_registers(points[p.edge % points.size()],
                                    points[(p.edge + 1) % points.size()],
                                    p.fraction);
  }

  /** @brief Interpolates a control point's birth register.
   * @param p Control descriptor.
   * @return Birth fraction, or zero when absent. */
  HS_HOT_INLINE float birth(const TimedControl &p) const {
    if (callback.births.empty())
      return 0.0f;
    const float a = callback.births[p.edge % points.size()];
    return a + (callback.births[(p.edge + 1) % points.size()] - a) * p.fraction;
  }

private:
  const Fragments &points;          /**< Source fragments. */
  const GeometryCallback &callback; /**< Borrowed geometry descriptor. */
  size_t first;                     /**< First cached source vertex. */
  math::Vector *previous = nullptr; /**< Edge-local previous lattice. */
};

/** @brief Evaluated control point with time and birth registers. */
struct TimedVertex {
  Fragment fragment; /**< Geometry and shader registers. */
  float time;        /**< Capture fraction. */
  float birth;       /**< Appearance fraction. */
};

/** @brief Emits a still edge's ordinary samples at their column capture times.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Timed plotting pipeline.
 * @param canvas Destination canvas.
 * @param a,b Evaluated endpoints.
 * @param shader Erased time-aware shader.
 * @param terminal Include the final endpoint.
 * @param balanced Stretch the adaptive sample spacing.
 * @param screen_branch Optional branch selecting projected samples.
 */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void
timed_walk_edge(PipelineRef &pipeline, Canvas &canvas, const TimedVertex &a,
                const TimedVertex &b, TimedShaderFn shader, bool terminal,
                bool balanced, const int *screen_branch = nullptr) {
  Temporal::count_geometry<&Temporal::GeometryCounters::legacy_edge_walks>();
  const auto span = make_geodesic_edge_span(a.fragment.pos, b.fragment.pos);
  const GeodesicEdgeSampler sample{
      a.fragment.pos,
      span.have_axis ? math::cross(span.axis, a.fragment.pos) : math::Vector{},
      span.total};
  constexpr float BASE = 2.0f * math::PI_F / W;
  const auto axes = screen_step_axes(pipeline);
  const auto step_at = [&](const SamplePT &value) HS_HOT_INLINE {
    return pipeline_screen_step<W, H>(
        pipeline, value, pipeline.world_transform_is_identity, &axes);
  };
  const float branch_start = screen_branch ? *screen_branch * (W / 2.0f) : 0.0f;
  const Temporal::ColumnPrediction acceptance{
      branch_start + W / 4.0f, std::numeric_limits<float>::infinity(),
      branch_start - 1.0f, branch_start + W / 2.0f + 1.0f};
  const auto *gate = screen_branch ? &acceptance : nullptr;
  const auto emit = [&](float s, const math::Vector &position,
                        float alpha_scale = 1.0f) HS_HOT_INLINE {
    const auto projected = pipeline.probe(canvas, position, 1.0f);
    const float time =
        screen_branch
            ? hs::clamp(
                  0.5f + Temporal::centered_remainder(
                             projected.x - branch_start - W / 4.0f, float(W)) /
                             (W / 2.0f),
                  0.0f, 1.0f)
            : Temporal::column_fraction<W>(projected.x);
    const float birth = a.birth + (b.birth - a.birth) * s;
    Fragment registers;
    if constexpr (Cfg.interpolate_registers)
      registers = Fragment::lerp_registers(a.fragment, b.fragment, s);
    registers.pos = position;
    auto evaluate_boundary = [&](float at, math::PixelCoords &p, Pixel &color,
                                 float &age, float &alpha) HS_HOT_INLINE {
      if (birth > at)
        return false;
      Fragment f = registers;
      shader(f.pos, f, at);
      p = pipeline.probe(canvas, f.pos, at);
      color = f.color.color;
      age = f.age;
      alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
      return true;
    };
    Temporal::BoundaryEvaluator boundary(evaluate_boundary);
    auto prediction = gate ? *gate
                           : Temporal::ColumnPrediction{
                                 0.0f, std::numeric_limits<float>::infinity()};
    prediction.boundary = &boundary;
    prediction.record_sample = birth <= time;
    boundary.sample_exists = prediction.record_sample;
    Fragment f = registers;
    if (prediction.record_sample) {
      shader(f.pos, f, time);
      f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
    } else {
      f.color.alpha = 0.0f;
    }
    pipeline.plot_at(canvas, f.pos, f.color.color, f.age, f.color.alpha, time,
                     &prediction);
  };
  if (span.total < math::EPS_GEOMETRIC) {
    if (terminal)
      emit(0.0f, a.fragment.pos);
    return;
  }
  if constexpr (!Cfg.single_pass) {
    auto value = sample(0.0f);
    const float first_step = step_at(value);
    if (span.total <= first_step) {
      emit(0.0f, a.fragment.pos);
      if (terminal)
        emit(1.0f, b.fragment.pos);
      return;
    }
    ScratchScope guard(scratch_arena_a);
    auto *steps = scratch_arena_a.allocate_n<float>(rasterize_step_budget<W>());
    size_t size = 0;
    float simulated = 0.0f;
    while (simulated < span.total && size < rasterize_step_budget<W>()) {
      const float step = size == 0 ? first_step : step_at(value);
      steps[size++] = step;
      simulated += step;
      if (simulated < span.total)
        value = sample(simulated / span.total);
    }
    const float scale = span.total / simulated;
    emit(0.0f, newton_unit(sample.pos(0.0f)));
    float distance = 0.0f;
    const size_t limit = terminal ? size : size - 1;
    for (size_t i = 0; i < limit; ++i) {
      distance += steps[i] * scale;
      const bool end = terminal && i + 1 == limit;
      const float s = replay_t(end, distance, span.total);
      emit(s, end ? b.fragment.pos : newton_unit(sample.pos(s)));
    }
    return;
  }
  float distance = 0.0f;
  float gap = 0.0f;
  size_t steps = 0;
  do {
    const float s = distance / span.total;
    const auto value = sample(s);
    float step = step_at(value);
    const float default_step = step;
    if (balanced)
      step = balanced_step(step, BASE);
    const auto position = s == 0.0f ? a.fragment.pos : value.pos.normalized();
    emit(s, position, balanced ? step / default_step : 1.0f);
    advance_distance(distance, span.total, step, gap);
    if (distance >= span.total) {
      if (terminal)
        emit(1.0f, b.fragment.pos,
             balanced ? endpoint_alpha_scale(gap, default_step, step) : 1.0f);
      break;
    }
  } while (++steps < 2 * rasterize_step_budget<W>());
}

/** @brief Finds a cut edge's crossing of one column at one frame endpoint.
 * @tparam W Canvas width.
 * @param pipeline Projection pipeline.
 * @param canvas Frame context.
 * @param geometry Edge evaluator.
 * @param a,b Cut-edge controls.
 * @param t Frame endpoint, zero or one.
 * @param column Unwrapped boundary column including footprint extension.
 * @return Clamped arc fraction along the cut edge. */
template <int W>
HS_HOT_FLASH_MEMBER float
timed_crossing(PipelineRef &pipeline, Canvas &canvas,
               const TimedGeometry &geometry, const TimedControl &a,
               const TimedControl &b, float t, float column) {
  const auto va = geometry.position(a, t);
  const auto vb = geometry.position(b, t);
  const auto pa = pipeline.probe(canvas, va, t);
  const auto pb = pipeline.probe(canvas, vb, t);
  Temporal::count_geometry<&Temporal::GeometryCounters::solve_probes>(2);
  const float xa =
      a.screen.x + Temporal::centered_remainder(pa.x - a.screen.x, float(W));
  const float xb =
      b.screen.x + Temporal::centered_remainder(pb.x - b.screen.x, float(W));
  if (fabsf(xb - xa) < 1e-6f)
    return column <= xa ? 0.0f : 1.0f;
  float s = hs::clamp((column - xa) / (xb - xa), 0.0f, 1.0f);
  if (s > 0.0f && s < 1.0f) {
    const auto pc = pipeline.probe(canvas, timed_edge_point(va, vb, s), t);
    Temporal::count_geometry<&Temporal::GeometryCounters::solve_probes>();
    const float xc = xa + Temporal::centered_remainder(pc.x - pa.x, float(W));
    if ((xc < column) == (xa < column))
      s += (1.0f - s) * (column - xc) / (xb - xc);
    else
      s *= (column - xa) / (xc - xa);
  }
  return hs::clamp(s, 0.0f, 1.0f);
}

/** @brief Tests a still edge's footprint-padded span against arm boundaries.
 * @tparam W,H Canvas dimensions.
 * @param a,b Unwrapped endpoint projections.
 * @return Whether the edge requires a boundary cut. */
template <int W, int H>
HS_HOT_FLASH_MEMBER bool timed_boundary_edge(math::PixelCoords a,
                                             math::PixelCoords b) {
  b.x = a.x + Temporal::centered_remainder(b.x - a.x, float(W));
  const auto va = math::pixel_to_vector<W, H>(a.x, a.y);
  const auto vb = math::pixel_to_vector<W, H>(b.x, b.y);
  const auto span = make_geodesic_edge_span(va, vb);
  int col_start, col_length;
  if (!geodesic_col_span<W>(va, vb, span, col_start, col_length))
    return true;
  const float start =
      a.x + Temporal::centered_remainder(float(col_start) - a.x, float(W));
  const float lo = fminf(fminf(a.x, b.x), start) - 1.0f;
  const float hi = fmaxf(fmaxf(a.x, b.x), start + col_length) + 1.0f;
  return floorf(lo / (W / 2.0f)) != floorf(hi / (W / 2.0f));
}

/** @brief Actual column retained between successive samples of one part. */
struct TimedPredictionState {
  float column = 0.0f;  /**< Last projected column. */
  float pending = 0.0f; /**< Estimated advance to the next edge's start. */
  int branch = 0;       /**< Current half-turn branch. */
  bool valid = false;   /**< Whether a preceding sample is available. */
};

/** @brief Tests a world-space edge against the padded arm-meridian plane.
 * @param a,b World endpoints.
 * @param bound Angular displacement including the splat footprint.
 * @return Whether the edge requires separate arm parts. */
HS_HOT_FLASH_MEMBER inline bool timed_prediction_cut(const math::Vector &a,
                                                     const math::Vector &b,
                                                     float bound) {
  const float padding = sinf(fminf(bound, math::PI_F * 0.5f));
  return a.z * b.z <= 0.0f || fabsf(a.z) <= padding || fabsf(b.z) <= padding ||
         math::dot(a, b) < -0.999f;
}

/** @brief Converts an angular bound to columns at a world point.
 * @tparam W Canvas width.
 * @param world Unit world position.
 * @param radians Angular displacement bound.
 * @param sine_bound Sine of the angular displacement.
 * @param azimuth The bound measures azimuth directly.
 * @return Conservative longitudinal displacement in columns. */
template <int W>
HS_HOT_INLINE float timed_prediction_motion(const math::Vector &world,
                                            float radians, float sine_bound,
                                            bool azimuth = false) {
  if (azimuth)
    return radians * (W / (2.0f * math::PI_F));
  const float sine = sqrtf(fmaxf(0.0f, 1.0f - world.y * world.y));
  if (radians >= math::PI_F / 2.0f || sine <= sine_bound)
    return W / 4.0f;
  return asinf(sine_bound / sine) * (W / (2.0f * math::PI_F));
}

/** @brief Estimates the longitudinal advance along a world-space chord.
 * @tparam W Canvas width.
 * @param sampler World-space geodesic.
 * @param start,end Material fractions.
 * @param length Chord arc length.
 * @return Signed column advance. */
template <int W>
HS_HOT_INLINE float timed_prediction_advance(const GeodesicEdgeSampler &sampler,
                                             float start, float end,
                                             float length) {
  const auto middle = sampler((start + end) * 0.5f);
  const float sin2 = fmaxf(1e-7f, 1.0f - middle.pos.y * middle.pos.y);
  return (middle.pos.x * middle.tan.z - middle.pos.z * middle.tan.x) / sin2 *
         (W / (2.0f * math::PI_F)) * length * (end - start);
}

/** @brief Final transformed control retained across adjacent source edges. */
struct TimedRigidCarry {
  TimedVertex vertex{}; /**< Evaluated position, registers, time and birth. */
  math::Vector source;  /**< Evaluated point before world stages. */
  float column = 0.0f;  /**< Unwrapped actual control column. */
  uint32_t index = 0;   /**< Original source vertex index. */
  int branch = 0;       /**< Arm branch owning the cached evaluation. */
  bool valid = false;   /**< Whether the cached control is reusable. */
};

/** @brief Evaluates one rigid control and corrects its capture prediction once.
 * @tparam W,H Canvas dimensions.
 * @param pipeline World transformation chain.
 * @param canvas Frame context.
 * @param geometry Source evaluator.
 * @param point Source control.
 * @param expected Predicted unwrapped column.
 * @param branch Arm branch.
 * @param radians,sine_bound Frame motion bound and its sine.
 * @return Evaluated control and actual column. */
template <int W, int H>
HS_O3_FN HS_NOINLINE_NOCLONE TimedRigidCarry resolve_rigid_control(
    PipelineRef &pipeline, Canvas &canvas, const TimedGeometry &geometry,
    const TimedControl &point, float expected, int branch, float radians,
    float sine_bound) {
  constexpr float HALF = W / 2.0f;
  TimedRigidCarry result;
  result.index = point.edge + (point.fraction == 1.0f ? 1 : 0);
  result.branch = branch;
  result.valid = true;
  result.vertex.fragment = geometry.registers(point);
  result.vertex.birth = geometry.birth(point);
  result.vertex.time = hs::clamp(expected / HALF - branch, 0.0f, 1.0f);
  const auto evaluate = [&]() HS_HOT_INLINE {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::control_evaluations>();
    const auto local = geometry.position(point, result.vertex.time);
    result.source = local;
    result.vertex.fragment.pos =
        pipeline.world_at(canvas, local, result.vertex.time);
    const auto projected =
        math::vector_to_pixel<W, H>(result.vertex.fragment.pos);
    result.column = expected + Temporal::centered_remainder(
                                   projected.x - expected, float(W));
  };
  evaluate();
  const float motion =
      timed_prediction_motion<W>(result.vertex.fragment.pos, radians,
                                 sine_bound, geometry.azimuth_bound());
  const float corrected = hs::clamp(result.column / HALF - branch, 0.0f, 1.0f);
  if (fabsf(corrected - result.vertex.time) * HALF > 0.5f ||
      motion * fabsf(corrected - result.vertex.time) > 0.1f) {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::control_corrections>();
    result.vertex.time = corrected;
    evaluate();
  }
  return result;
}

/** @brief One adaptive rigid control and its cumulative world arc. */
struct TimedRigidNode {
  TimedRigidCarry control; /**< Evaluated source and world registers. */
  float fraction;          /**< Material fraction on the original edge. */
  float distance;          /**< Cumulative spacing arc. */
};

/** @brief Maximum adaptive rigid controls retained for one original edge. */
inline constexpr size_t MAX_RIGID_CONTROLS = 64;

/** @brief Ordinary geodesic samplers retained while traversing one control span. */
class TimedRigidSampler {
public:
  /** @brief Borrows the ordered adaptive controls.
   * @param nodes Control sequence. */
  explicit TimedRigidSampler(std::span<const TimedRigidNode> nodes)
      : nodes(nodes) {}

  /** @brief Selects the containing span without rebuilding it per sample.
   * @param distance Cumulative spacing arc. */
  HS_HOT_INLINE void seek(float distance) {
    size_t next = index == std::numeric_limits<size_t>::max() ? 0 : index;
    while (next > 0 && distance < nodes[next].distance)
      --next;
    while (next + 2 < nodes.size() && distance > nodes[next + 1].distance)
      ++next;
    if (next != index)
      prepare(next);
    const float length = nodes[index + 1].distance - nodes[index].distance;
    fraction =
        length > math::EPS_GEOMETRIC
            ? hs::clamp((distance - nodes[index].distance) / length, 0.0f, 1.0f)
            : 0.0f;
  }
  /** @brief Rewinds for the ordinary normalized replay. */
  void reset() { index = std::numeric_limits<size_t>::max(); }
  /** @brief Returns the current span's first control.
   * @return Beginning control. */
  HS_HOT_INLINE const TimedRigidNode &first() const { return nodes[index]; }
  /** @brief Returns the current span's final control.
   * @return Ending control. */
  HS_HOT_INLINE const TimedRigidNode &last() const { return nodes[index + 1]; }
  /** @brief Returns the local material interpolation fraction.
   * @return Fraction within the selected control span. */
  HS_HOT_INLINE float amount() const { return fraction; }
  /** @brief Evaluates the source point without calling moving stages.
   * @return Source-space position for the shader. */
  HS_HOT_INLINE math::Vector source_position() const {
    if (fraction == 0.0f)
      return first().control.source;
    if (fraction == 1.0f)
      return last().control.source;
    return newton_unit(source.pos(fraction));
  }
  /** @brief Evaluates the already-transformed point.
   * @return World-space position for projection. */
  HS_HOT_INLINE math::Vector world_position() const {
    if (fraction == 0.0f)
      return first().control.vertex.fragment.pos;
    if (fraction == 1.0f)
      return last().control.vertex.fragment.pos;
    return newton_unit(world.pos(fraction));
  }

private:
  std::span<const TimedRigidNode> nodes;             /**< Borrowed controls. */
  size_t index = std::numeric_limits<size_t>::max(); /**< Current span. */
  float fraction = 0.0f;        /**< Position within the span. */
  GeodesicEdgeSampler world{};  /**< World-space edge sampler. */
  GeodesicEdgeSampler source{}; /**< Source edge sampler. */
  /** @brief Rebuilds samplers only when crossing a geometric control.
   * @param next New control-span index. */
  HS_O3_FN HS_NOINLINE_NOCLONE void prepare(size_t next) {
    index = next;
    const auto &a = first().control;
    const auto &b = last().control;
    const auto ws =
        make_geodesic_edge_span(a.vertex.fragment.pos, b.vertex.fragment.pos);
    world = {a.vertex.fragment.pos,
             ws.have_axis ? math::cross(ws.axis, a.vertex.fragment.pos)
                          : math::Vector{},
             ws.total};
    const auto ss = make_geodesic_edge_span(a.source, b.source);
    source = {a.source,
              ss.have_axis ? math::cross(ss.axis, a.source) : math::Vector{},
              ss.total};
  }
};

/** @brief Frozen step schedules spanning a captured edge's time interval. */
struct TimedCadence {
  float *first = nullptr; /**< Normalized cumulative steps at the first time. */
  float *last = nullptr;  /**< Normalized cumulative steps at the last time. */
  size_t count = 0; /**< Maximum interval count; zero retains common spacing. */
  float begin = 0.0f; /**< First covered frame fraction. */
  float end = 1.0f;   /**< Last covered frame fraction. */

  /** @brief Material phase resolved at one frame fraction. */
  struct AtTime {
    const float *first; /**< First endpoint's cumulative steps. */
    const float *last;  /**< Last endpoint's cumulative steps. */
    float amount;       /**< Endpoint interpolation fraction. */
    float total;        /**< Cumulative extent at this time. */
    size_t intervals;   /**< Material interval count at this time. */

    /** @brief Resolves one ordinal without rebuilding the frame phase.
     * @param index Sample ordinal.
     * @param terminal Whether the trailing endpoint is included.
     * @param exists Receives whether the ordinal exists.
     * @return Material fraction on the edge. */
    HS_HOT_INLINE float fraction_at(size_t index, bool terminal,
                                    bool &exists) const {
      exists = index < intervals || (terminal && index == intervals);
      return index >= intervals
                 ? 1.0f
                 : (first[index] + (last[index] - first[index]) * amount) /
                       total;
    }
  };

  /** @brief Resolves the cadence scale once for a shared sample time.
   * @param time Capture fraction.
   * @return Reusable material phase. */
  HS_HOT_INLINE AtTime at_time(float time) const {
    const float amount = hs::clamp((time - begin) / (end - begin), 0.0f, 1.0f);
    const auto cumulative = [&](size_t i) HS_HOT_INLINE {
      return first[i] + (last[i] - first[i]) * amount;
    };
    size_t intervals = count;
    while (intervals > 1 && cumulative(intervals - 1) >= 1.0f)
      --intervals;
    return {first, last, amount, intervals ? cumulative(intervals) : 1.0f,
            intervals};
  }

  /** @brief Resolves an indexed sample at its capture time.
   * @param index Sample ordinal, including the leading endpoint.
   * @param time Capture fraction.
   * @param terminal Whether the trailing endpoint is included.
   * @param exists Receives whether the ordinal exists at this time.
   * @return Material fraction on the edge. */
  HS_HOT_INLINE float fraction_at(size_t index, float time, bool terminal,
                                  bool &exists) const {
    return at_time(time).fraction_at(index, terminal, exists);
  }
};

/** @brief Counts or caches the normalized cumulative steps of a frozen edge.
 * @tparam W,H Canvas dimensions.
 * @param sample Frozen world-space geodesic.
 * @param values Optional cumulative-step output, including the leading zero.
 * @param minimum Minimum count, extending the final step past the endpoint.
 * @return Number of intervals written or counted. */
template <int W, int H>
HS_HOT_FLASH_MEMBER size_t prepare_cadence_schedule(
    const GeodesicEdgeSampler &sample, float *values, size_t minimum) {
  constexpr float BASE = 2.0f * math::PI_F / W;
  size_t count = 0;
  float distance = 0.0f, step = BASE;
  if (values)
    values[0] = 0.0f;
  while ((distance < sample.total_dist || count < minimum) &&
         count < rasterize_step_budget<W>()) {
    if (distance < sample.total_dist) {
      const auto value = sample(distance / sample.total_dist);
      step = screen_step<W, H>(value.pos, value.tan, BASE);
    }
    distance += step;
    ++count;
    if (values)
      values[count] = distance / sample.total_dist;
  }
  return count;
}

/** @brief Prepares time-varying spacing when an edge gains or loses an interval.
 * @tparam W,H Canvas dimensions.
 * @param pipeline World and screen pipeline.
 * @param canvas Destination frame.
 * @param geometry Moving source geometry.
 * @param original Original edge registers.
 * @param nodes Captured controls delimiting the edge and its time interval.
 * @return Scratch-backed schedules, or an empty schedule for constant cadence. */
template <int W, int H>
HS_HOT_FLASH_MEMBER TimedCadence prepare_timed_cadence(
    PipelineRef &pipeline, Canvas &canvas, const TimedGeometry &geometry,
    const TimedControl &original, std::span<const TimedRigidNode> nodes) {
  constexpr float TAP_TIME = 2.0f / W;
  TimedCadence cadence;
  const float first_time = nodes.front().control.vertex.time;
  const float last_time = nodes.back().control.vertex.time;
  cadence.begin = fmaxf(0.0f, fminf(first_time, last_time) - TAP_TIME);
  cadence.end = fminf(1.0f, fmaxf(first_time, last_time) + TAP_TIME);
  GeodesicEdgeSampler samples[2];
  size_t counts[2];
  auto head = original, tail = original;
  head.fraction = nodes.front().fraction;
  tail.fraction = nodes.back().fraction;
  for (size_t i = 0; i < 2; ++i) {
    const float time = i == 0 ? cadence.begin : cadence.end;
    const auto a =
        pipeline.world_at(canvas, geometry.position(head, time), time);
    const auto b =
        pipeline.world_at(canvas, geometry.position(tail, time), time);
    Temporal::count_geometry<&Temporal::GeometryCounters::control_evaluations>(
        2);
    const auto span = make_geodesic_edge_span(a, b);
    if (span.total < math::EPS_GEOMETRIC)
      return {};
    samples[i] = {a,
                  span.have_axis ? math::cross(span.axis, a) : math::Vector{},
                  span.total};
    counts[i] = prepare_cadence_schedule<W, H>(samples[i], nullptr, 0);
  }
  if (counts[0] == counts[1])
    return {};
  cadence.count = std::max(counts[0], counts[1]);
  cadence.first = scratch_arena_a.allocate_n<float>(cadence.count + 1);
  cadence.last = scratch_arena_a.allocate_n<float>(cadence.count + 1);
  prepare_cadence_schedule<W, H>(samples[0], cadence.first, cadence.count);
  prepare_cadence_schedule<W, H>(samples[1], cadence.last, cadence.count);
  return cadence;
}

/** @brief Replays one original edge across its transformed control spans.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Screen emission pipeline.
 * @param canvas Frame context.
 * @param geometry Source evaluator for boundary taps.
 * @param original Original edge descriptor.
 * @param nodes Adaptive controls with cumulative material arcs.
 * @param spacing_sample Whole-edge sampler evaluated at a common frame fraction.
 * @param cadence Optional time-varying interval schedule.
 * @param shader Pure source-space shader.
 * @param state Actual last landing and current arm branch.
 * @param terminal Include the final endpoint.
 * @param balanced Apply balanced sampling.
 * @param bounded Restrict emission to the part footprint. */
template <int W, int H, RasterConfig Cfg>
HS_O3_FN HS_NOINLINE_NOCLONE void walk_rigid_samples(
    PipelineRef &pipeline, Canvas &canvas, const TimedGeometry &geometry,
    const TimedControl &original, std::span<const TimedRigidNode> nodes,
    const GeodesicEdgeSampler &spacing_sample, const TimedCadence &cadence,
    TimedShaderFn shader, TimedPredictionState &state, bool terminal,
    bool balanced, bool bounded) {
  constexpr float BASE = 2.0f * math::PI_F / W;
  constexpr float HALF = W / 2.0f;
  const float total = nodes.back().distance;
  TimedRigidSampler sampler(nodes);
  Temporal::ColumnPrediction prediction{(state.branch + 0.5f) * HALF,
                                        std::numeric_limits<float>::infinity()};
  if (bounded) {
    prediction.min_column = state.branch * HALF - 1.0f;
    prediction.max_column = (state.branch + 1) * HALF + 1.0f;
  }
  size_t current_index = 0;
  const auto emit = [&](float distance, float alpha_scale) HS_HOT_INLINE {
    sampler.seek(distance);
    bool cadence_exists = true;
    if (cadence.count) {
      for (int iteration = 0; iteration < 2; ++iteration) {
        const float guess = sampler.first().control.vertex.time +
                            (sampler.last().control.vertex.time -
                             sampler.first().control.vertex.time) *
                                sampler.amount();
        distance = cadence.fraction_at(current_index, guess, terminal,
                                       cadence_exists) *
                   total;
        sampler.seek(distance);
      }
    }
    const auto &a = sampler.first();
    const auto &b = sampler.last();
    const float amount = sampler.amount();
    const float time = a.control.vertex.time +
                       (b.control.vertex.time - a.control.vertex.time) * amount;
    const float birth =
        a.control.vertex.birth +
        (b.control.vertex.birth - a.control.vertex.birth) * amount;
    TimedControl material = original;
    material.fraction = a.fraction + (b.fraction - a.fraction) * amount;
    auto evaluate_boundary = [&](float at, math::PixelCoords &projected,
                                 Pixel &color, float &age,
                                 float &alpha) HS_HOT_INLINE {
      auto tap_material = material;
      if (cadence.count) {
        bool exists;
        const float fraction =
            cadence.fraction_at(current_index, at, terminal, exists);
        if (!exists)
          return false;
        tap_material.fraction =
            nodes.front().fraction +
            (nodes.back().fraction - nodes.front().fraction) * fraction;
      }
      if (geometry.birth(tap_material) > at)
        return false;
      Fragment f;
      if constexpr (Cfg.interpolate_registers)
        f = geometry.registers(tap_material);
      f.pos = geometry.position(tap_material, at);
      shader(f.pos, f, at);
      projected = pipeline.probe(canvas, f.pos, at);
      color = f.color.color;
      age = f.age;
      alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
      return true;
    };
    Temporal::BoundaryEvaluator boundary(evaluate_boundary);
    boundary.force_timing = original.motion > 0.0f;
    if (boundary.force_timing)
      Temporal::count_geometry<
          &Temporal::GeometryCounters::interior_timed_samples>();
    prediction.boundary = &boundary;
    prediction.record_sample = cadence_exists && birth <= time;
    boundary.sample_exists = prediction.record_sample;
    Fragment f;
    if constexpr (Cfg.interpolate_registers)
      f = geometry.registers(material);
    f.pos = sampler.source_position();
    if (prediction.record_sample) {
      shader(f.pos, f, time);
      f.color.alpha = balanced_sample_alpha(f.color.alpha, alpha_scale);
    } else {
      f.color.alpha = 0.0f;
    }
    Temporal::count_geometry<&Temporal::GeometryCounters::rigid_samples>();
    state.column = pipeline.plot_transformed_at(
        canvas, sampler.world_position(), f.color.color, f.age, f.color.alpha,
        time, &prediction);
    state.valid = true;
  };
  if (total < math::EPS_GEOMETRIC) {
    if (terminal)
      emit(0.0f, 1.0f);
    return;
  }
  if constexpr (!Cfg.single_pass) {
    ScratchScope guard(scratch_arena_a);
    auto *steps = scratch_arena_a.allocate_n<float>(rasterize_step_budget<W>());
    size_t count = 0;
    float distance = 0.0f;
    while (distance < total && count < rasterize_step_budget<W>()) {
      const auto value = spacing_sample(distance / total);
      const float step = screen_step<W, H>(value.pos, value.tan, BASE);
      steps[count++] = step;
      distance += step;
    }
    const float scale = total / distance;
    sampler.reset();
    emit(0.0f, 1.0f);
    distance = 0.0f;
    const size_t limit = cadence.count ? cadence.count
                         : terminal    ? count
                                       : count - 1;
    for (size_t i = 0; i < limit; ++i) {
      current_index = i + 1;
      distance += steps[std::min(i, count - 1)] * scale;
      emit(terminal && i + 1 == limit ? total : fminf(distance, total), 1.0f);
    }
  } else {
    float distance = 0.0f, gap = 0.0f;
    size_t count = 0;
    do {
      const auto value = spacing_sample(distance / total);
      const float normal = screen_step<W, H>(value.pos, value.tan, BASE);
      const float step = balanced ? balanced_step(normal, BASE) : normal;
      emit(distance, balanced ? step / normal : 1.0f);
      advance_distance(distance, total, step, gap);
      if (distance >= total) {
        if (terminal)
          emit(total,
               balanced ? endpoint_alpha_scale(gap, normal, step) : 1.0f);
        break;
      }
    } while (++count < 2 * rasterize_step_budget<W>());
  }
}

/** @brief Tests the capture and deformation bounds of a control interval.
 * @tparam W,H Canvas dimensions.
 * @param left,right Captured endpoints.
 * @param branch Owning arm pass.
 * @param radians,sine_bound Frame motion bound and its sine.
 * @param azimuth_only Whether motion preserves latitude.
 * @param column Predicted midpoint column output.
 * @return Whether the interval satisfies the adaptive tolerances. */
template <int W, int H>
HS_HOT_INLINE bool rigid_interval_is_resolved(
    const TimedRigidCarry &left, const TimedRigidCarry &right, int branch,
    float radians, float sine_bound, bool azimuth_only, float &column) {
  constexpr float HALF = W / 2.0f;
  const auto span = make_geodesic_edge_span(left.vertex.fragment.pos,
                                            right.vertex.fragment.pos);
  const auto middle = timed_edge_point(left.vertex.fragment.pos,
                                       right.vertex.fragment.pos, 0.5f);
  const auto projected = math::vector_to_pixel<W, H>(middle);
  const float mean_column = (left.column + right.column) * 0.5f;
  column = mean_column +
           Temporal::centered_remainder(projected.x - mean_column, float(W));
  const float time = (left.vertex.time + right.vertex.time) * 0.5f;
  const float captured = hs::clamp(column / HALF - branch, 0.0f, 1.0f);
  const float time_error = fabsf(captured - time);
  const float motion =
      timed_prediction_motion<W>(middle, radians, sine_bound, azimuth_only);
  const float sine = sqrtf(fmaxf(1e-8f, 1.0f - middle.y * middle.y));
  const float deformation = radians *
                            fabsf(right.vertex.time - left.vertex.time) *
                            span.total * (W / (2.0f * math::PI_F)) / sine;
  return time_error * HALF <= 0.5f && time_error * motion <= 0.1f &&
         deformation <= 0.1f;
}

/** @brief Captures adaptive controls and their common-time spacing.
 * @tparam W,H Canvas dimensions.
 * @param pipeline World and screen pipeline.
 * @param canvas Frame context.
 * @param geometry Source evaluator.
 * @param a,b Source controls.
 * @param reference_a,reference_b Endpoint world positions for the forecast.
 * @param state Previous sample landing.
 * @param carry Previous transformed control.
 * @param radians,sine_bound Frame motion bound and its sine.
 * @param initial_time Capture time at a new part.
 * @param nodes Storage for adaptive controls.
 * @return Control count, or zero when the scratch budget is exceeded. */
template <int W, int H>
HS_HOT_INLINE size_t prepare_rigid_controls(
    PipelineRef &pipeline, Canvas &canvas, const TimedGeometry &geometry,
    const TimedControl &a, const TimedControl &b,
    const math::Vector &reference_a, const math::Vector &reference_b,
    const TimedPredictionState &state, const TimedRigidCarry &carry,
    float radians, float sine_bound, float initial_time,
    TimedRigidNode *nodes) {
  constexpr float HALF = W / 2.0f;
  const float expected = state.valid ? state.column + state.pending
                                     : (state.branch + initial_time) * HALF;
  nodes[0].control =
      carry.valid && a.fraction == 0.0f && carry.index == a.edge &&
              carry.branch == state.branch
          ? carry
          : resolve_rigid_control<W, H>(pipeline, canvas, geometry, a, expected,
                                        state.branch, radians, sine_bound);
  const auto reference_span = make_geodesic_edge_span(reference_a, reference_b);
  const GeodesicEdgeSampler reference_sample{
      reference_a,
      reference_span.have_axis ? math::cross(reference_span.axis, reference_a)
                               : math::Vector{},
      reference_span.total};
  const float predicted_end =
      nodes[0].control.column +
      timed_prediction_advance<W>(reference_sample, 0.0f, 1.0f,
                                  reference_span.total);
  nodes[1].control =
      resolve_rigid_control<W, H>(pipeline, canvas, geometry, b, predicted_end,
                                  state.branch, radians, sine_bound);
  size_t count = 2;
  nodes[0].fraction = a.fraction;
  nodes[0].distance = 0.0f;
  nodes[1].fraction = b.edge == a.edge ? b.fraction : 1.0f;
  nodes[1].distance = 0.0f;
  for (size_t i = 0; i + 1 < count;) {
    float column;
    if (rigid_interval_is_resolved<W, H>(nodes[i].control, nodes[i + 1].control,
                                         state.branch, radians, sine_bound,
                                         geometry.azimuth_bound(), column)) {
      ++i;
      continue;
    }
    const float midpoint_fraction =
        (nodes[i].fraction + nodes[i + 1].fraction) * 0.5f;
    if (count == MAX_RIGID_CONTROLS || midpoint_fraction == nodes[i].fraction ||
        midpoint_fraction == nodes[i + 1].fraction) {
      Temporal::count_geometry<
          &Temporal::GeometryCounters::rigid_fallback_edges>();
      return 0;
    }
    TimedControl material = a;
    material.fraction = midpoint_fraction;
    for (size_t j = count; j > i + 1; --j)
      nodes[j] = nodes[j - 1];
    nodes[i + 1].control =
        resolve_rigid_control<W, H>(pipeline, canvas, geometry, material,
                                    column, state.branch, radians, sine_bound);
    nodes[i + 1].fraction = material.fraction;
    nodes[i + 1].distance = 0.0f;
    ++count;
    Temporal::count_geometry<&Temporal::GeometryCounters::refined_controls>();
  }
  return count;
}

/** @brief Assigns common-time arc distances to captured controls.
 * @tparam W,H Canvas dimensions.
 * @param pipeline World and screen pipeline.
 * @param canvas Frame context.
 * @param geometry Source evaluator.
 * @param a,b Source controls.
 * @param nodes Captured controls and output distances.
 * @param spacing_sample Common-time geodesic output. */
template <int W, int H>
HS_HOT_INLINE void
prepare_rigid_spacing(PipelineRef &pipeline, Canvas &canvas,
                      const TimedGeometry &geometry, const TimedControl &a,
                      const TimedControl &b, std::span<TimedRigidNode> nodes,
                      GeodesicEdgeSampler &spacing_sample) {
  const size_t count = nodes.size();
  const float spacing_time =
      (nodes[0].control.vertex.time + nodes[count - 1].control.vertex.time) *
      0.5f;
  const auto spacing_a = pipeline.world_at(
      canvas, geometry.position(a, spacing_time), spacing_time);
  const auto spacing_b = pipeline.world_at(
      canvas, geometry.position(b, spacing_time), spacing_time);
  Temporal::count_geometry<&Temporal::GeometryCounters::common_time_edges>();
  Temporal::count_geometry<&Temporal::GeometryCounters::control_evaluations>(2);
  const auto spacing_span = make_geodesic_edge_span(spacing_a, spacing_b);
  spacing_sample = GeodesicEdgeSampler{
      spacing_a,
      spacing_span.have_axis ? math::cross(spacing_span.axis, spacing_a)
                             : math::Vector{},
      spacing_span.total};
  const float fraction_span = nodes[count - 1].fraction - nodes[0].fraction;
  for (size_t i = 1; i < count; ++i)
    nodes[i].distance = fraction_span > 0.0f
                            ? spacing_span.total *
                                  (nodes[i].fraction - nodes[0].fraction) /
                                  fraction_span
                            : 0.0f;
}

/** @brief Walks a rigid edge between transformed, time-carrying controls.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline World and screen pipeline.
 * @param canvas Frame context.
 * @param geometry Source evaluator.
 * @param a,b Source controls.
 * @param reference_a,reference_b Endpoint world positions for the column forecast.
 * @param shader Pure time-aware shader.
 * @param state Previous sample landing.
 * @param carry Previous transformed control.
 * @param radians,sine_bound Frame motion bound and its sine.
 * @param terminal Include the final endpoint.
 * @param balanced Apply balanced sampling.
 * @param initial_time Capture time at a new part.
 * @param bounded Apply the arm part's footprint window.
 * @return False when adaptive controls exceed the scratch budget. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER bool
walk_rigid_edge(PipelineRef &pipeline, Canvas &canvas,
                const TimedGeometry &geometry, const TimedControl &a,
                const TimedControl &b, const math::Vector &reference_a,
                const math::Vector &reference_b, TimedShaderFn shader,
                TimedPredictionState &state, TimedRigidCarry &carry,
                float radians, float sine_bound, bool terminal, bool balanced,
                float initial_time, bool bounded) {
  ScratchScope control_scope(scratch_arena_a);
  auto *nodes = scratch_arena_a.allocate_n<TimedRigidNode>(MAX_RIGID_CONTROLS);
  GeodesicEdgeSampler spacing_sample;
  const size_t count = prepare_rigid_controls<W, H>(
      pipeline, canvas, geometry, a, b, reference_a, reference_b, state, carry,
      radians, sine_bound, initial_time, nodes);
  if (count == 0)
    return false;
  prepare_rigid_spacing<W, H>(pipeline, canvas, geometry, a, b,
                              std::span<TimedRigidNode>(nodes, count),
                              spacing_sample);
  TimedCadence cadence;
  if constexpr (!Cfg.single_pass) {
    if (a.motion > 0.0f)
      cadence = prepare_timed_cadence<W, H>(
          pipeline, canvas, geometry, a,
          std::span<const TimedRigidNode>(nodes, count));
  }
  walk_rigid_samples<W, H, Cfg>(pipeline, canvas, geometry, a,
                                std::span<const TimedRigidNode>(nodes, count),
                                spacing_sample, cadence, shader, state,
                                terminal, balanced, bounded);
  const auto &last = nodes[count - 1].control;
  state.pending = last.column - state.column;
  state.valid = true;
  carry = last;
  return true;
}

/** @brief Emits an edge with coherent capture prediction and one correction.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Timed pipeline.
 * @param canvas Target frame.
 * @param geometry Borrowed source geometry.
 * @param a,b Material endpoints.
 * @param shader Pure time-aware shader.
 * @param state Previous landed sample and branch.
 * @param radians Angular motion bound.
 * @param terminal Include the final endpoint.
 * @param balanced Apply balanced sampling.
 * @param initial_time Capture time when the part has no preceding sample.
 * @param bounded Restrict emitted taps to the branch's footprint interval. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void
timed_predict_edge(PipelineRef &pipeline, Canvas &canvas,
                   const TimedGeometry &geometry, const TimedControl &a,
                   const TimedControl &b, TimedShaderFn shader,
                   TimedPredictionState &state, float radians, bool terminal,
                   bool balanced, float initial_time, bool bounded) {
  constexpr float HALF = W / 2.0f;
  const float motion = radians > 0.0f ? W / 4.0f : 0.0f;
  constexpr float BASE = 2.0f * math::PI_F / W;
  const float start_column = state.valid ? state.column + state.pending
                                         : (state.branch + initial_time) * HALF;
  const float start_time =
      hs::clamp(start_column / HALF - state.branch, 0.0f, 1.0f);
  const auto world_a =
      pipeline.world_at(canvas, geometry.position(a, start_time), start_time);
  const auto same_time_b =
      pipeline.world_at(canvas, geometry.position(b, start_time), start_time);
  const auto initial_span = make_geodesic_edge_span(world_a, same_time_b);
  const GeodesicEdgeSampler initial_sample{
      world_a,
      initial_span.have_axis ? math::cross(initial_span.axis, world_a)
                             : math::Vector{},
      initial_span.total};
  const float end_time = hs::clamp(
      (start_column + timed_prediction_advance<W>(initial_sample, 0.0f, 1.0f,
                                                  initial_span.total)) /
              HALF -
          state.branch,
      0.0f, 1.0f);
  const auto world_b =
      pipeline.world_at(canvas, geometry.position(b, end_time), end_time);
  const auto span = make_geodesic_edge_span(world_a, world_b);
  const GeodesicEdgeSampler sample{
      world_a,
      span.have_axis ? math::cross(span.axis, world_a) : math::Vector{},
      span.total};
  float previous_fraction = 0.0f;
  bool emitted = false;
  const auto emit = [&](float fraction, float alpha_scale) HS_HOT_INLINE {
    TimedControl material = a;
    const float end_fraction = b.edge == a.edge ? b.fraction : 1.0f;
    material.fraction = a.fraction + (end_fraction - a.fraction) * fraction;
    float expected = emitted ? state.column + timed_prediction_advance<W>(
                                                  sample, previous_fraction,
                                                  fraction, span.total)
                             : start_column;
    if (bounded)
      expected = hs::clamp(expected, state.branch * HALF - 1.0f,
                           (state.branch + 1) * HALF + 1.0f);
    float time = hs::clamp(expected / HALF - state.branch, 0.0f, 1.0f);
    auto position = geometry.position(material, time);
    Temporal::ColumnPrediction prediction{
        expected, motion > 0.0f ? 0.1f * HALF / motion
                                : std::numeric_limits<float>::infinity()};
    if (bounded) {
      prediction.min_column = state.branch * HALF - 1.0f;
      prediction.max_column = (state.branch + 1) * HALF + 1.0f;
    }
    auto evaluate_boundary = [&](float at, math::PixelCoords &projected,
                                 Pixel &color, float &age,
                                 float &alpha) HS_HOT_INLINE {
      if (geometry.birth(material) > at)
        return false;
      Fragment fragment;
      if constexpr (Cfg.interpolate_registers)
        fragment = geometry.registers(material);
      fragment.pos = geometry.position(material, at);
      shader(fragment.pos, fragment, at);
      projected = pipeline.probe(canvas, fragment.pos, at);
      color = fragment.color.color;
      age = fragment.age;
      alpha = balanced_sample_alpha(fragment.color.alpha, alpha_scale);
      return true;
    };
    Temporal::BoundaryEvaluator boundary(evaluate_boundary);
    boundary.force_timing = radians > 0.0f;
    if (boundary.force_timing)
      Temporal::count_geometry<
          &Temporal::GeometryCounters::interior_timed_samples>();
    prediction.boundary = &boundary;
    const auto draw = [&](float at) HS_HOT_INLINE {
      boundary.sample_exists = geometry.birth(material) <= at;
      if (!boundary.sample_exists) {
        auto unborn = prediction;
        unborn.record_sample = false;
        return pipeline.plot_at(canvas, position, Pixel{}, 0.0f, 0.0f, at,
                                &unborn);
      }
      Fragment fragment;
      if constexpr (Cfg.interpolate_registers)
        fragment = geometry.registers(material);
      fragment.pos = position;
      shader(fragment.pos, fragment, at);
      const float alpha =
          balanced_sample_alpha(fragment.color.alpha, alpha_scale);
      return pipeline.plot_at(canvas, fragment.pos, fragment.color.color,
                              fragment.age, alpha, at, &prediction);
    };
    float landed = draw(time);
    if (prediction.needs_correction(landed)) {
      const float corrected_time =
          hs::clamp(landed / HALF - state.branch, 0.0f, 1.0f);
      if (corrected_time != time || (landed >= prediction.min_column &&
                                     landed <= prediction.max_column)) {
        Temporal::count_geometry<&Temporal::GeometryCounters::corrections>();
        time = corrected_time;
        prediction.tolerance = std::numeric_limits<float>::infinity();
        position = geometry.position(material, time);
        landed = draw(time);
      }
    }
    state.column = landed;
    state.valid = true;
    previous_fraction = fraction;
    emitted = true;
  };
  if (span.total < math::EPS_GEOMETRIC) {
    if (terminal)
      emit(0.0f, 1.0f);
    state.pending = 0.0f;
    return;
  }
  if constexpr (!Cfg.single_pass) {
    ScratchScope guard(scratch_arena_a);
    auto *steps = scratch_arena_a.allocate_n<float>(rasterize_step_budget<W>());
    size_t size = 0;
    float distance = 0.0f;
    while (distance < span.total && size < rasterize_step_budget<W>()) {
      const auto value = sample(distance / span.total);
      const float step = screen_step<W, H>(value.pos, value.tan, BASE);
      steps[size++] = step;
      distance += step;
    }
    const float scale = span.total / distance;
    emit(0.0f, 1.0f);
    distance = 0.0f;
    const size_t limit = terminal ? size : size - 1;
    for (size_t i = 0; i < limit; ++i) {
      distance += steps[i] * scale;
      emit(replay_t(terminal && i + 1 == limit, distance, span.total), 1.0f);
    }
  } else {
    float distance = 0.0f;
    float gap = 0.0f;
    size_t count = 0;
    do {
      const float fraction = distance / span.total;
      const auto value = sample(fraction);
      const float default_step = screen_step<W, H>(value.pos, value.tan, BASE);
      const float step =
          balanced ? balanced_step(default_step, BASE) : default_step;
      emit(fraction, balanced ? step / default_step : 1.0f);
      advance_distance(distance, span.total, step, gap);
      if (distance >= span.total) {
        if (terminal)
          emit(1.0f,
               balanced ? endpoint_alpha_scale(gap, default_step, step) : 1.0f);
        break;
      }
    } while (++count < 2 * rasterize_step_budget<W>());
  }
  state.pending = emitted ? timed_prediction_advance<W>(
                                sample, previous_fraction, 1.0f, span.total)
                          : 0.0f;
}

/** @brief Walks a stationary source edge with column-time shading.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline World and screen pipeline.
 * @param canvas Frame context.
 * @param points Endpoint source fragments.
 * @param geometry Source evaluator and birth fractions.
 * @param a,b Source controls.
 * @param reference_a,world_b Transformed endpoints.
 * @param shader Pure time-aware shader.
 * @param terminal Include the final endpoint.
 * @param balanced Apply balanced sampling. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_INLINE void walk_stationary_geometry_edge(
    PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
    const TimedGeometry &geometry, const TimedControl &a, const TimedControl &b,
    const math::Vector &reference_a, const math::Vector &world_b,
    TimedShaderFn shader, bool terminal, bool balanced) {
  constexpr float HALF = W / 2.0f;
  const auto pa = math::vector_to_pixel<W, H>(reference_a);
  auto pb = math::vector_to_pixel<W, H>(world_b);
  pb.x = pa.x + Temporal::centered_remainder(pb.x - pa.x, float(W));
  Temporal::count_geometry<&Temporal::GeometryCounters::end_probes>(2);
  const bool stationary_cut = timed_boundary_edge<W, H>(pa, pb);
  const int first_branch = static_cast<int>(
      floorf((stationary_cut ? fminf(pa.x, pb.x) - 1.0f : pa.x) / HALF));
  const int last_branch =
      stationary_cut
          ? static_cast<int>(floorf((fmaxf(pa.x, pb.x) + 1.0f) / HALF))
          : first_branch;
  for (int branch = first_branch; branch <= last_branch; ++branch) {
    TimedVertex va{points[a.edge % points.size()], pa.x / HALF - branch,
                   geometry.birth(a)};
    TimedVertex vb{points[(a.edge + 1) % points.size()], pb.x / HALF - branch,
                   geometry.birth(b)};
    Temporal::count_geometry<&Temporal::GeometryCounters::parts>();
    timed_walk_edge<W, H, Cfg>(pipeline, canvas, va, vb, shader, terminal,
                               balanced, stationary_cut ? &branch : nullptr);
  }
}

/** @brief Walks moving source edges using landed-column time prediction.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Timed pipeline.
 * @param canvas Target frame.
 * @param points Endpoint source geometry.
 * @param shader Pure time-aware shader.
 * @param opts Polyline options.
 * @param callback Source geometry descriptor.
 * @param radians Combined callback and world angular bound.
 * @param first_edge First original edge of a restricted run.
 * @param edge_count Restricted edge count, or zero for the whole polyline. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void rasterize_predicted_copy(
    PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
    TimedShaderFn shader, RasterOptions opts, const GeometryCallback &callback,
    float radians, size_t first_edge = 0, size_t edge_count = 0) {
  constexpr float HALF = W / 2.0f;
  const float sine_bound = sinf(radians);
  const size_t edges =
      edge_count ? edge_count
                 : (opts.loop.is_closed() ? points.size() : points.size() - 1);
  const float padding = radians + math::RADIANS_PER_COLUMN<W>;
  const CullRegion region(canvas.clip(), canvas.render_x_clip(), radians);
  const bool balanced =
      Cfg.sampling_policy == RasterSamplingPolicy::BALANCED ||
      (Cfg.sampling_policy == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling);
  size_t origin = first_edge;
  if (!edge_count && opts.loop.is_closed()) {
    auto previous = pipeline.world_at(canvas, points[0].pos, 1.0f);
    for (size_t i = 0; i < edges; ++i) {
      const auto next =
          pipeline.world_at(canvas, points[(i + 1) % points.size()].pos, 1.0f);
      if (timed_prediction_cut(previous, next, padding)) {
        origin = i;
        break;
      }
      previous = next;
    }
  }
  TimedLatticeCarry carry;
  TimedPredictionState state;
  TimedRigidCarry rigid_carry;
  auto world_a = pipeline.world_at(canvas, points[origin].pos, 1.0f);
  for (size_t offset = 0; offset < edges; ++offset) {
    const size_t edge = origin + offset;
    const auto &source_a = points[edge % points.size()].pos;
    const auto &source_b = points[(edge + 1) % points.size()].pos;
    const auto reference_a = world_a;
    const auto world_b = pipeline.world_at(canvas, source_b, 1.0f);
    const bool cut = timed_prediction_cut(world_a, world_b, padding);
    const float motion_a = timed_prediction_motion<W>(
        world_a, radians, sine_bound, callback.azimuth_bound);
    const float motion_b = timed_prediction_motion<W>(
        world_b, radians, sine_bound, callback.azimuth_bound);
    world_a = world_b;
    if (!edge_visible_in_region<W, H>(pipeline, region, source_a, source_b)) {
      state.valid = false;
      rigid_carry.valid = false;
      continue;
    }
    ScratchScope guard(scratch_arena_a);
    TimedGeometry geometry(points, callback, edge, carry);
    TimedControl a{static_cast<uint32_t>(edge), 0.0f, {}, motion_a};
    TimedControl b{static_cast<uint32_t>(edge), 1.0f, {}, motion_b};
    const bool terminal =
        offset + 1 == edges && !opts.loop.is_closed() && !opts.omit_end;
    if (pipeline.cull_bound() == 0.0f && geometry.stationary(edge)) {
      walk_stationary_geometry_edge<W, H, Cfg>(
          pipeline, canvas, points, geometry, a, b, reference_a, world_b,
          shader, terminal, balanced);
      state.valid = false;
      rigid_carry.valid = false;
      continue;
    }
    const auto endpoint_probe = [&](TimedControl &point) HS_HOT_INLINE {
      point.screen =
          pipeline.probe(canvas, geometry.position(point, 1.0f), 1.0f);
      point.motion =
          Temporal::MotionBound{radians, callback.azimuth_bound}.columns<W, H>(
              point.screen.y);
      Temporal::count_geometry<&Temporal::GeometryCounters::end_probes>();
    };
    const auto initial_capture = [&](const TimedControl &point,
                                     int branch) HS_HOT_INLINE {
      const auto probe = [&](float time) HS_HOT_INLINE {
        return pipeline.probe(canvas, geometry.position(point, time), time);
      };
      if (point.motion >= Temporal::MAX_CAPTURE_COLUMNS<W>)
        Temporal::count_geometry<&Temporal::GeometryCounters::held_anchors>();
      return Temporal::solve_capture_branch<W>(
          probe, point.screen, point.motion, callback.linear_pair, branch);
    };
    if (!cut) {
      float initial_time = 1.0f;
      if (!state.valid) {
        endpoint_probe(a);
        state.branch = static_cast<int>(floorf(a.screen.x / HALF));
        initial_time = initial_capture(a, state.branch);
        Temporal::count_geometry<&Temporal::GeometryCounters::parts>();
      }
      if (!(pipeline.rigid_world_stages() &&
            walk_rigid_edge<W, H, Cfg>(
                pipeline, canvas, geometry, a, b, reference_a, world_b, shader,
                state, rigid_carry, radians, sine_bound, terminal, balanced,
                initial_time, false)))
        timed_predict_edge<W, H, Cfg>(pipeline, canvas, geometry, a, b, shader,
                                      state, radians, terminal, balanced,
                                      initial_time, false);
      continue;
    }
    endpoint_probe(a);
    endpoint_probe(b);
    b.screen.x = a.screen.x + Temporal::centered_remainder(
                                  b.screen.x - a.screen.x, float(W));
    const float columns = fmaxf(a.motion, b.motion) + 1.0f;
    const int low_branch = static_cast<int>(
        floorf((fminf(a.screen.x, b.screen.x) - columns) / HALF));
    const int high_branch = static_cast<int>(
        floorf((fmaxf(a.screen.x, b.screen.x) + columns) / HALF));
    TimedPredictionState continuation;
    TimedRigidCarry rigid_continuation;
    for (int branch = low_branch; branch <= high_branch; ++branch) {
      const bool forward = b.screen.x >= a.screen.x;
      const float low = branch * HALF - 1.0f;
      const float high = (branch + 1) * HALF + 1.0f;
      const float start =
          timed_crossing<W>(pipeline, canvas, geometry, a, b,
                            forward ? 0.0f : 1.0f, forward ? low : high);
      const float end =
          timed_crossing<W>(pipeline, canvas, geometry, a, b,
                            forward ? 1.0f : 0.0f, forward ? high : low);
      if (end <= start)
        continue;
      TimedPredictionState part = state;
      part.valid = part.valid && part.branch == branch && start == 0.0f;
      part.branch = branch;
      float initial_time = 1.0f;
      if (!part.valid) {
        initial_time = initial_capture(a, branch);
        Temporal::count_geometry<&Temporal::GeometryCounters::parts>();
      }
      auto part_carry = rigid_carry;
      if (!(pipeline.rigid_world_stages() &&
            walk_rigid_edge<W, H, Cfg>(pipeline, canvas, geometry, a, b,
                                       reference_a, world_b, shader, part,
                                       part_carry, radians, sine_bound,
                                       terminal, balanced, initial_time, true)))
        timed_predict_edge<W, H, Cfg>(pipeline, canvas, geometry, a, b, shader,
                                      part, radians, terminal, balanced,
                                      initial_time, true);
      if (end == 1.0f) {
        continuation = part;
        rigid_continuation = part_carry;
      }
    }
    state = continuation;
    rigid_carry = rigid_continuation;
  }
}

#if HS_ENABLE_TEST_HOOKS
/** @brief Forces eligible dense geometry through the ordinary timed edge route. */
inline bool g_force_timed_edge_walk = false;
#endif

/** @brief Cached local and projected endpoints of one dense control. */
struct DenseControl {
  math::Vector previous;   /**< Source position at the previous endpoint. */
  math::PixelCoords start; /**< Previous projected endpoint. */
  math::PixelCoords end;   /**< Current projected endpoint. */
  bool safe = true;        /**< Endpoint stays within the affine motion band. */
  bool dense = false;      /**< Outgoing edge qualifies for dense sampling. */
  bool wide = false;       /**< Outgoing edge requires a curvature check. */
};

/** @brief One affine dense point at an arm pass. */
using DenseCapture = Temporal::LinearCapture;

/** @brief Emits one dense material point and its independent arm passes.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Prepared screen-entry pipeline.
 * @param canvas Destination frame.
 * @param a,b Original source shading controls.
 * @param da,db Cached endpoint trajectories.
 * @param fraction Material interpolation fraction.
 * @param birth_a,birth_b Endpoint appearance fractions.
 * @param shader Pure time-aware source shader.
 * @param interpolated Optional gap sample between captured controls. */
template <int W, int H, RasterConfig Cfg>
HS_O3_FN HS_NOINLINE_NOCLONE void
emit_dense_point(PipelineRef &pipeline, Canvas &canvas, const Fragment &a,
                 const Fragment &b, const DenseControl &da,
                 const DenseControl &db, float fraction, float birth_a,
                 float birth_b, TimedShaderFn shader,
                 const DenseCapture *interpolated = nullptr) {
  const auto mix = [fraction](float x, float y)
                       HS_HOT_INLINE { return x + (y - x) * fraction; };
  const math::PixelCoords start{mix(da.start.x, db.start.x),
                                mix(da.start.y, db.start.y)};
  const math::PixelCoords end{mix(da.end.x, db.end.x), mix(da.end.y, db.end.y)};
  const math::Vector previous =
      da.previous + (db.previous - da.previous) * fraction;
  const math::Vector current = a.pos + (b.pos - a.pos) * fraction;
  const float birth = mix(birth_a, birth_b);
  Fragment registers;
  if constexpr (Cfg.interpolate_registers)
    registers = Fragment::lerp_registers(a, b, fraction);
  const auto shade = [&](float t, Fragment &f) HS_HOT_INLINE {
    f = registers;
    f.pos = (previous + (current - previous) * t).normalized();
    shader(f.pos, f, t);
  };
  const auto emit_capture = [&](const DenseCapture &capture) HS_HOT_INLINE {
    auto boundary_at = [&](float t, math::PixelCoords &position, Pixel &color,
                           float &age, float &alpha) HS_HOT_INLINE {
      if (birth > t)
        return false;
      Fragment f;
      shade(t, f);
      position = {start.x + (end.x - start.x) * t,
                  start.y + (end.y - start.y) * t};
      if (position.x < -float(W) || position.x >= 2.0f * W)
        position.x = Temporal::centered_remainder(position.x, float(W));
      color = f.color.color;
      age = f.age;
      alpha = f.color.alpha;
      return true;
    };
    Temporal::BoundaryEvaluator boundary(boundary_at);
    const float dx = end.x - start.x, dy = end.y - start.y;
    const float residual =
        fabsf(capture.position.x - (start.x + dx * capture.time)) +
        fabsf(capture.position.y - (start.y + dy * capture.time));
    boundary.force_timing =
        residual + Temporal::tap_motion_error<W>(capture.position.x, dx, dy) >
        Temporal::MAX_SINGLE_TIME_TAP_ERROR;
    if (boundary.force_timing)
      Temporal::count_geometry<
          &Temporal::GeometryCounters::interior_timed_samples>();
    else
      Temporal::count_geometry<
          &Temporal::GeometryCounters::interior_single_time_samples>();
    Temporal::ColumnPrediction prediction{
        capture.position.x, std::numeric_limits<float>::infinity()};
    prediction.boundary = &boundary;
    prediction.record_sample = birth <= capture.time;
    boundary.sample_exists = prediction.record_sample;
    Fragment f = registers;
    if (prediction.record_sample)
      shade(capture.time, f);
    else
      f.color.alpha = 0.0f;
    auto position = capture.position;
    if (position.x < -float(W) || position.x >= 2.0f * W)
      position.x = Temporal::centered_remainder(position.x, float(W));
    pipeline.plot_projected_at(canvas, position, f.color.color, f.age,
                               f.color.alpha, capture.time, &prediction);
  };
  if (interpolated)
    emit_capture(*interpolated);
  else
    Temporal::linear_captures<W>(start, end, emit_capture);
}

/** @brief Emits one dense edge using cached endpoint trajectories.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Prepared screen-entry pipeline.
 * @param canvas Destination frame.
 * @param points Original source shading controls.
 * @param controls Cached local and projected endpoints.
 * @param i Original edge index.
 * @param birth_a,birth_b Endpoint appearance fractions.
 * @param shader Pure time-aware source shader.
 * @param emit_end Whether to include the final open endpoint. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_INLINE void
emit_dense_edge(PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
                const DenseControl *controls, size_t i, float birth_a,
                float birth_b, TimedShaderFn shader, bool emit_end) {
  Temporal::count_geometry<&Temporal::GeometryCounters::dense_controls>();
  Temporal::count_geometry<&Temporal::GeometryCounters::dense_run_edges>();
  emit_dense_point<W, H, Cfg>(pipeline, canvas, points[i], points[i],
                              controls[i], controls[i], 0.0f, birth_a, birth_a,
                              shader);
  const size_t next = (i + 1) % points.size();
  auto following = controls[next];
  const float shift = controls[i].end.x +
                      Temporal::centered_remainder(
                          following.end.x - controls[i].end.x, float(W)) -
                      following.end.x;
  following.end.x += shift;
  following.start.x += shift;
  constexpr float HALF = W / 2.0f;
  const float low_column = fminf(fminf(controls[i].start.x, controls[i].end.x),
                                 fminf(following.start.x, following.end.x));
  const float high_column = fmaxf(fmaxf(controls[i].start.x, controls[i].end.x),
                                  fmaxf(following.start.x, following.end.x));
  const int first_branch = static_cast<int>(floorf((low_column - 1.0f) / HALF));
  const int last_branch = static_cast<int>(floorf((high_column + 1.0f) / HALF));
  for (int branch = first_branch; branch <= last_branch; ++branch) {
    const auto capture_endpoint = [branch](
                                      const DenseControl &point) HS_HOT_INLINE {
      const float dx = point.end.x - point.start.x;
      const float t =
          hs::clamp((point.start.x - branch * HALF) / (HALF - dx), 0.0f, 1.0f);
      return DenseCapture{{point.start.x + dx * t,
                           point.start.y + (point.end.y - point.start.y) * t},
                          t};
    };
    const auto a = capture_endpoint(controls[i]);
    const auto b = capture_endpoint(following);
    const float dx = b.position.x - a.position.x;
    const float dy = b.position.y - a.position.y;
    const int steps =
        static_cast<int>(ceilf(sqrtf(dx * dx + dy * dy) / SCREEN_STEP_PX));
    for (int step = 1; step < steps; ++step) {
      const float amount = float(step) / steps;
      const DenseCapture capture{
          {a.position.x + dx * amount, a.position.y + dy * amount},
          a.time + (b.time - a.time) * amount};
      if (capture.position.x < branch * HALF - 1.0f ||
          capture.position.x > (branch + 1) * HALF + 1.0f)
        continue;
      Temporal::count_geometry<
          &Temporal::GeometryCounters::dense_gap_samples>();
      emit_dense_point<W, H, Cfg>(pipeline, canvas, points[i], points[next],
                                  controls[i], following, amount, birth_a,
                                  birth_b, shader, &capture);
    }
  }
  if (emit_end) {
    Temporal::count_geometry<&Temporal::GeometryCounters::dense_controls>();
    emit_dense_point<W, H, Cfg>(pipeline, canvas, points[next], points[next],
                                following, following, 0.0f, birth_b, birth_b,
                                shader);
  }
}

/** @brief Walks cached dense runs and dispatches sparse edge intervals.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Prepared timed pipeline.
 * @param canvas Destination frame.
 * @param points Original source shading controls.
 * @param controls Cached local and projected endpoints.
 * @param shader Pure time-aware source shader.
 * @param opts Polyline sampling options.
 * @param callback Original geometry callback and birth registers.
 * @param radians Combined angular motion bound.
 * @param visible Whether the motion-padded dense region intersects the clip. */
template <int W, int H, RasterConfig Cfg>
HS_O3_FN HS_NOINLINE_NOCLONE void
walk_dense_runs(PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
                const DenseControl *controls, TimedShaderFn shader,
                const RasterOptions &opts, const GeometryCallback &callback,
                float radians, bool visible) {
  const size_t edges =
      opts.loop.is_closed() ? points.size() : points.size() - 1;
  const auto birth = [&](size_t i) HS_HOT_INLINE {
    return callback.births.empty() ? 0.0f : callback.births[i];
  };
  for (size_t i = 0; i < edges;) {
    if (!controls[i].dense) {
      const size_t first = i;
      while (i < edges && !controls[i].dense)
        ++i;
      Temporal::count_geometry<&Temporal::GeometryCounters::dense_sparse_edges>(
          i - first);
      auto sparse_opts = opts;
      sparse_opts.loop = RasterLoop{};
      sparse_opts.omit_end =
          opts.loop.is_closed() || opts.omit_end || i < edges;
      rasterize_predicted_copy<W, H, Cfg>(pipeline, canvas, points, shader,
                                          sparse_opts, callback, radians, first,
                                          i - first);
      continue;
    }
    if (i == 0 || !controls[i - 1].dense)
      Temporal::count_geometry<&Temporal::GeometryCounters::dense_runs>();
    if (!visible) {
      ++i;
      continue;
    }
    emit_dense_edge<W, H, Cfg>(pipeline, canvas, points, controls, i, birth(i),
                               birth((i + 1) % points.size()), shader,
                               i + 1 == edges && !opts.loop.is_closed() &&
                                   !opts.omit_end);
    ++i;
  }
}

/** @brief Result of preparing cached dense controls. */
enum class DensePreparation {
  REJECTED,    /**< The ordinary timed route is required. */
  COMPLETE,    /**< The whole draw was culled. */
  VISIBLE,     /**< Dense and sparse runs may be visible. */
  SPARSE_ONLY, /**< Only sparse runs require visibility checks. */
};

/** @brief Prepares projected controls and tests dense-run eligibility.
 * @tparam W,H Canvas dimensions.
 * @param pipeline Prepared timed pipeline.
 * @param canvas Destination frame.
 * @param points Current source controls.
 * @param controls Storage for cached trajectories and run flags.
 * @param opts Polyline sampling options.
 * @param callback Original geometry callback.
 * @param radians Combined angular motion bound.
 * @return Visibility and dispatch state for the cached controls. */
template <int W, int H>
HS_O3_FN HS_NOINLINE_NOCLONE DensePreparation prepare_dense_controls(
    PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
    DenseControl *controls, const RasterOptions &opts,
    const GeometryCallback &callback, float radians) {
  const size_t edges =
      opts.loop.is_closed() ? points.size() : points.size() - 1;
  const float sine = sinf(radians), polar_limit = 2.0f * sine * sine;
  float min_x = std::numeric_limits<float>::infinity();
  float max_x = -min_x, min_y = min_x, max_y = max_x;
  const auto adjacent = [](const math::PixelCoords &a,
                           const math::PixelCoords &b) HS_HOT_INLINE {
    const float dx = Temporal::centered_remainder(b.x - a.x, float(W));
    const float dy = b.y - a.y;
    return dx * dx + dy * dy;
  };
  bool polar = false, sparse = false;
  for (size_t i = 0; i < points.size(); ++i) {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_candidate_controls>();
    const auto world = pipeline.world_at(canvas, points[i].pos, 1.0f);
    controls[i].safe = 1.0f - world.y * world.y > polar_limit;
    polar |= !controls[i].safe;
    controls[i].dense = false;
    auto projected = math::vector_to_pixel<W, H>(world);
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_endpoint_probes>();
    Temporal::count_geometry<&Temporal::GeometryCounters::end_probes>();
    if (i)
      projected.x = controls[i - 1].end.x +
                    Temporal::centered_remainder(
                        projected.x - controls[i - 1].end.x, float(W));
    controls[i].end = projected;
    min_x = fminf(min_x, projected.x);
    max_x = fmaxf(max_x, projected.x);
    min_y = fminf(min_y, projected.y);
    max_y = fmaxf(max_y, projected.y);
  }
  bool candidate = false;
  for (size_t i = 0; i < edges; ++i) {
    const size_t next = (i + 1) % points.size();
    const float separation = adjacent(controls[i].end, controls[next].end);
    const bool close = separation <= 6.25f;
    controls[i].wide = separation > 4.0f;
    sparse |= !close;
    if (!close)
      Temporal::count_geometry<
          &Temporal::GeometryCounters::dense_spacing_edges>();
    else if (!controls[i].safe || !controls[next].safe)
      Temporal::count_geometry<
          &Temporal::GeometryCounters::dense_polar_edges>();
    controls[i].dense = close && controls[i].safe && controls[next].safe;
    candidate |= controls[i].dense;
  }
  if (sparse)
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_reject_spacing>();
  if (polar)
    Temporal::count_geometry<&Temporal::GeometryCounters::dense_reject_polar>();
  if (!candidate)
    return DensePreparation::REJECTED;
  Temporal::count_geometry<&Temporal::GeometryCounters::dense_culls>();
  const CullRegion region(canvas.clip(), canvas.render_x_clip(), radians);
  const float pad = region.column_pad_rows<W, H>(min_y - 1.0f, max_y + 1.0f);
  const bool visible =
      region.could_intersect_rows<H>(min_y - 1.0f, max_y + 1.0f) &&
      region.arc_reaches_columns<W>(
          static_cast<int>(floorf(min_x - 1.0f)),
          static_cast<int>(ceilf(max_x - min_x + 3.0f)), pad);
  if (!visible && !sparse && !polar) {
    Temporal::count_geometry<&Temporal::GeometryCounters::dense_runs>();
    return DensePreparation::COMPLETE;
  }
  bool stationary = pipeline.cull_bound() == 0.0f;
  for (size_t i = 0; i < points.size(); ++i) {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::callback_evaluations>();
    controls[i].previous = callback.evaluate(i, 0.0f);
    const auto &a = controls[i].previous;
    const auto &b = points[i].pos;
    stationary = stationary && a.x == b.x && a.y == b.y && a.z == b.z;
    const auto world = pipeline.world_at(canvas, a, 0.0f);
    controls[i].safe &= 1.0f - world.y * world.y > polar_limit;
    auto projected = math::vector_to_pixel<W, H>(world);
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_endpoint_probes>();
    projected.x =
        controls[i].end.x +
        Temporal::centered_remainder(projected.x - controls[i].end.x, float(W));
    controls[i].safe &= fabsf(projected.x - controls[i].end.x) <
                        Temporal::MAX_CAPTURE_COLUMNS<W>;
    controls[i].start = projected;
  }
  if (stationary)
    return DensePreparation::REJECTED;
  candidate = false;
  for (size_t i = 0; i < edges; ++i) {
    const size_t next = (i + 1) % points.size();
    controls[i].dense &= controls[i].safe && controls[next].safe;
    if (controls[i].dense && controls[i].wide) {
      const auto sample_count = [](const math::PixelCoords &a,
                                   const math::PixelCoords &b) HS_HOT_INLINE {
        const float dx = Temporal::centered_remainder(b.x - a.x, float(W));
        const float dy = b.y - a.y;
        return ceilf(sqrtf(dx * dx + dy * dy) / SCREEN_STEP_PX);
      };
      if (sample_count(controls[i].start, controls[next].start) !=
          sample_count(controls[i].end, controls[next].end)) {
        controls[i].dense = false;
        Temporal::count_geometry<
            &Temporal::GeometryCounters::dense_cadence_edges>();
        continue;
      }
      const auto flat = [&](const math::Vector &a, const math::Vector &b,
                            const math::PixelCoords &pa,
                            const math::PixelCoords &pb,
                            float time) HS_HOT_INLINE {
        const auto middle = (a + b).normalized();
        const auto world = pipeline.world_at(canvas, middle, time);
        auto projected = math::vector_to_pixel<W, H>(world);
        Temporal::count_geometry<
            &Temporal::GeometryCounters::dense_midpoint_probes>();
        const float x =
            pa.x + 0.5f * Temporal::centered_remainder(pb.x - pa.x, float(W));
        const float dx =
            Temporal::centered_remainder(projected.x - x, float(W));
        const float dy = projected.y - 0.5f * (pa.y + pb.y);
        return dx * dx + dy * dy <= 0.000025f;
      };
      controls[i].dense = flat(points[i].pos, points[next].pos, controls[i].end,
                               controls[next].end, 1.0f) &&
                          flat(controls[i].previous, controls[next].previous,
                               controls[i].start, controls[next].start, 0.0f);
      if (controls[i].dense)
        Temporal::count_geometry<
            &Temporal::GeometryCounters::dense_wide_edges>();
      else
        Temporal::count_geometry<
            &Temporal::GeometryCounters::dense_curved_edges>();
    }
    candidate |= controls[i].dense;
  }
  if (!candidate)
    return DensePreparation::REJECTED;
  return visible ? DensePreparation::VISIBLE : DensePreparation::SPARSE_ONLY;
}

/** @brief Selects and emits a dense, interpolating rigid polyline.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Prepared timed pipeline.
 * @param canvas Destination frame.
 * @param points Current source controls.
 * @param shader Pure time-aware shader.
 * @param opts Polyline options.
 * @param callback Interpolating source trajectory.
 * @param radians Combined angular motion bound.
 * @return True when this draw was handled as a dense run. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER bool
rasterize_dense(PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
                TimedShaderFn shader, RasterOptions opts,
                const GeometryCallback &callback, float radians) {
#if HS_ENABLE_TEST_HOOKS
  if (g_force_timed_edge_walk)
    return false;
#endif
  if (Cfg.sampling_policy == RasterSamplingPolicy::BALANCED ||
      (Cfg.sampling_policy == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling)) {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_balanced_candidates>();
    return false;
  }
  if (!callback.evaluate || !callback.linear_pair ||
      !pipeline.rigid_world_stages() || !pipeline.dense_screen_safe() ||
      points.size() < 3)
    return false;
  Temporal::count_geometry<&Temporal::GeometryCounters::dense_candidates>();
  if (radians >= math::PI_F / 2.0f) {
    Temporal::count_geometry<&Temporal::GeometryCounters::dense_reject_polar>();
    return false;
  }
  const size_t available =
      scratch_arena_b.get_capacity() - scratch_arena_b.get_offset();
  if (available < alignof(DenseControl) ||
      points.size() >
          (available - alignof(DenseControl)) / sizeof(DenseControl)) {
    Temporal::count_geometry<
        &Temporal::GeometryCounters::dense_reject_scratch>();
    return false;
  }
  ScratchScope guard(scratch_arena_b);
  auto *controls = scratch_arena_b.allocate_n<DenseControl>(points.size());
  const auto preparation = prepare_dense_controls<W, H>(
      pipeline, canvas, points, controls, opts, callback, radians);
  if (preparation == DensePreparation::REJECTED)
    return false;
  if (preparation == DensePreparation::COMPLETE)
    return true;
  walk_dense_runs<W, H, Cfg>(pipeline, canvas, points, controls, shader, opts,
                             callback, radians,
                             preparation == DensePreparation::VISIBLE);
  return true;
}

/** @brief Samples zero-motion geometry with time-aware registers.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Prepared timed pipeline.
 * @param canvas Destination frame.
 * @param points Current source controls.
 * @param shader Pure time-aware shader.
 * @param opts Polyline sampling options.
 * @param callback Original birth registers. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_INLINE void
rasterize_stationary_timed(PipelineRef &pipeline, Canvas &canvas,
                           const Fragments &points, TimedShaderFn shader,
                           const RasterOptions &opts,
                           const GeometryCallback &callback) {
  ScratchScope guard(scratch_arena_a);
  const size_t edges =
      opts.loop.is_closed() ? points.size() : points.size() - 1;
  auto *projections =
      scratch_arena_a.allocate_n<math::PixelCoords>(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    projections[i] = pipeline.probe(canvas, points[i].pos, 1.0f);
    Temporal::count_geometry<&Temporal::GeometryCounters::end_probes>();
  }
  size_t origin = 0;
  if (opts.loop.is_closed()) {
    for (size_t i = 0; i < edges; ++i) {
      if (timed_boundary_edge<W, H>(projections[i],
                                    projections[(i + 1) % points.size()])) {
        origin = i;
        break;
      }
    }
  }
  constexpr float HALF = W / 2.0f;
  const CullRegion region(canvas.clip(), canvas.render_x_clip(), 0.0f);
  const bool balanced =
      Cfg.sampling_policy == RasterSamplingPolicy::BALANCED ||
      (Cfg.sampling_policy == RasterSamplingPolicy::SELECTABLE &&
       opts.balanced_sampling);
  for (size_t offset = 0; offset < edges; ++offset) {
    const size_t a = (origin + offset) % points.size();
    const size_t b = (a + 1) % points.size();
    if (!edge_visible_in_region<W, H>(pipeline, region, points[a].pos,
                                      points[b].pos))
      continue;
    const auto pa = projections[a];
    auto pb = projections[b];
    pb.x = pa.x + Temporal::centered_remainder(pb.x - pa.x, float(W));
    const bool cut = timed_boundary_edge<W, H>(pa, pb);
    const int first = static_cast<int>(
        floorf((cut ? fminf(pa.x, pb.x) - 1.0f : pa.x) / HALF));
    const int last =
        cut ? static_cast<int>(floorf((fmaxf(pa.x, pb.x) + 1.0f) / HALF))
            : first;
    const bool terminal =
        offset + 1 == edges && !opts.loop.is_closed() && !opts.omit_end;
    for (int branch = first; branch <= last; ++branch) {
      TimedVertex va{points[a], pa.x / HALF - branch,
                     callback.births.empty() ? 0.0f : callback.births[a]};
      TimedVertex vb{points[b], pb.x / HALF - branch,
                     callback.births.empty() ? 0.0f : callback.births[b]};
      if (offset == 0 || cut)
        Temporal::count_geometry<&Temporal::GeometryCounters::parts>();
      timed_walk_edge<W, H, Cfg>(pipeline, canvas, va, vb, shader, terminal,
                                 balanced, cut ? &branch : nullptr);
    }
  }
}

HS_O3_BEGIN

/** @brief Angular bound between a meridian crossing's arc fraction and the
 * column walk sampler's parameter, radians. */
inline constexpr float COLUMN_ARC_PARAMETER_PAD = 1e-5f;

/** @brief Largest projected motion between adjacent column times, pixels, at
 * which a column holding an edge end splats its taps at its own time. */
inline constexpr float COLUMN_SINGLE_TAP_ERROR = 0.02f;

/** @brief Intersects a geodesic basis with a forward meridian.
 * @param start_side,end_side Endpoint distances from the meridian plane.
 * @param tangent_side Start tangent distance from the meridian plane.
 * @param start_along,tangent_along Basis projections along the meridian.
 * @param inverse_length Reciprocal arc length, or zero for a collapsed edge.
 * @param fraction Receives the material arc fraction.
 * @return Whether the minor arc crosses the forward meridian. */
HS_HOT_INLINE bool column_basis_crossing(float start_side, float end_side,
                                         float tangent_side, float start_along,
                                         float tangent_along,
                                         float inverse_length,
                                         float &fraction) {
  if (start_side * end_side > 0.0f || fabsf(start_side - end_side) < 1e-9f)
    return false;
  float angle = math::precise_atan2(-start_side, tangent_side);
  const float forward = tangent_side * start_along - start_side * tangent_along;
  if (angle < 0.0f) {
    angle += math::PI_F;
    if (forward > 0.0f)
      return false;
  } else if (forward < 0.0f)
    return false;
  fraction = hs::clamp(angle * inverse_length, 0.0f, 1.0f);
  return true;
}

/** @brief Intersects a minor arc with a destination meridian.
 * @param a,b Arc endpoint directions; a is unit length.
 * @param direction Horizontal meridian direction.
 * @param length Minor arc length.
 * @param tangent Unit tangent at the first endpoint.
 * @param fraction Resulting arc fraction.
 * @return Whether the meridian crosses the minor arc. */
HS_HOT_INLINE bool
column_edge_crossing(const math::Vector &a, const math::Vector &b,
                     const math::Vector &direction, float length,
                     const math::Vector &tangent, float &fraction) {
  const math::Vector normal(-direction.z, 0.0f, direction.x);
  return column_basis_crossing(
      math::dot(a, normal), math::dot(b, normal), math::dot(tangent, normal),
      math::dot(a, direction), math::dot(tangent, direction),
      length > 1e-7f ? 1.0f / length : 0.0f, fraction);
}

/** @brief A mesh edge as two vertex indices. */
struct MeshEdge {
  uint16_t u; ///< First endpoint vertex index.
  uint16_t v; ///< Second endpoint vertex index.
};

/** @brief Mesh edges sharing one column walk. */
struct ColumnEdgeBatch {
  size_t count; /**< Edge count. */
  void (*prepare)(Fragments &, const math::Vector &, const math::Vector &,
                  size_t); /**< Writes an edge's endpoint registers into the
                                walk's fragments from its frame-end
                                vertices and index. */
  const math::Vector *vertices; /**< Frame-end vertex positions. */
  const Temporal::VertexMotion<>
      &motion;           /**< Source captures and interpolants of the
                              vertices this frame. */
  const MeshEdge *pairs; /**< Vertex indices of each edge. */
  bool chained = false;  /**< Edges leave their end vertex unplotted. */
};

/** @brief Rejects edges whose swept arcs miss a draw's rows and columns.
 * @tparam W Display width.
 * @details Rows use the column walk's own sample cut; columns use the render
 * band widened by two columns on each side. */
template <int W> struct ColumnEdgeCull {
  float top;           /**< Highest unit-sphere y the walk emits. */
  float bottom;        /**< Lowest unit-sphere y the walk emits. */
  math::Vector first;  /**< First longitude half-plane normal. */
  math::Vector second; /**< Second longitude half-plane normal. */
  int mode = 0; /**< 0: rows only; 1: hidden outside either half-plane; 2:
                   hidden inside both. */

  /** @brief Prepares the bounds for one draw.
   * @param columns Render-band column predicate.
   * @param top_y,bottom_y Unit-sphere y limits of emitted samples. */
  HS_HOT_FLASH_MEMBER ColumnEdgeCull(const ClipRegion::XClip &columns,
                                     float top_y, float bottom_y)
      : top(top_y), bottom(bottom_y) {
    if (!columns.active)
      return;
    constexpr float MARGIN = 2.0f;
    const float start =
        (float(columns.rs) - 0.5f - MARGIN) * math::RADIANS_PER_COLUMN<W>;
    const float width = (float(columns.length(W)) + 2.0f * MARGIN) *
                        math::RADIANS_PER_COLUMN<W>;
    if (width >= 2.0f * math::PI_F)
      return;
    const float end = start + width;
    if (width <= math::PI_F) {
      mode = 1;
      first = {-sinf(start), 0.0f, cosf(start)};
      second = {sinf(end), 0.0f, -cosf(end)};
    } else {
      mode = 2;
      first = {-sinf(end), 0.0f, cosf(end)};
      second = {sinf(start), 0.0f, -cosf(start)};
    }
  }

  /** @brief Highest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height of the row two above the clip. */
  template <int H> HS_HOT_INLINE static float top_of(const Canvas &canvas) {
    return cosf(math::DisplayGeometry<H>::row_to_phi(
        float(std::max(0, canvas.clip().y_start - 2))));
  }

  /** @brief Lowest unit-sphere y a canvas's column walk emits.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows bound the walk.
   * @return The height of the row two below the clip. */
  template <int H> HS_HOT_INLINE static float bottom_of(const Canvas &canvas) {
    return cosf(math::DisplayGeometry<H>::row_to_phi(
        float(std::min(H, canvas.clip().y_end + 2))));
  }

  /** @brief The bounds of a canvas's column walk.
   * @tparam H Display height.
   * @param canvas Frame whose clip rows and render band bound the walk.
   * @return The cull. */
  template <int H> static ColumnEdgeCull of(const Canvas &canvas) {
    return ColumnEdgeCull(canvas.render_x_clip(), top_of<H>(canvas),
                          bottom_of<H>(canvas));
  }

  /** @brief Tests whether every geodesic between points a frame sweeps
   * stays outside the bounds.
   * @param points Positions whose normalized convex hull, padded, holds every
   * drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_FLASH_MEMBER bool rejects(const math::Vector *points, int count,
                                   float pad) const {
    return rejects_hull(points, count, pad);
  }

  /** @brief Tests whether every arc an edge sweeps over the frame stays
   * outside the bounds.
   * @param points Both endpoints at the frame's start and end.
   * @param pad Bound on a drawn endpoint's distance from its chord.
   * @return True when no drawn sample of the edge can be emitted. */
  HS_HOT_FLASH_MEMBER bool rejects(const math::Vector (&points)[4],
                                   float pad) const {
    return rejects_hull(points, 4, pad);
  }

private:
  /** @brief rejects() for any count, inlined so a fixed count unrolls.
   * @param points Positions whose padded hull holds every drawn sample.
   * @param count Number of points.
   * @param pad Bound on a drawn sample's distance from that hull.
   * @return True when no drawn sample can be emitted. */
  HS_HOT_INLINE bool rejects_hull(const math::Vector *points, int count,
                                  float pad) const {
    float low, high;
    if (Temporal::swept_heights(points, count, pad, low, high) &&
        (high < bottom || low > top))
      return true;
    if (mode == 1)
      return Temporal::swept_beyond(points, count, pad, -first) ||
             Temporal::swept_beyond(points, count, pad, -second);
    return mode == 2 && Temporal::swept_beyond(points, count, pad, first) &&
           Temporal::swept_beyond(points, count, pad, second);
  }
};

/** @brief One arm pass's walk over the columns an edge can light. */
struct ColumnSegment {
  int start = 0;    /**< First column in walk order. */
  int count = 0;    /**< Columns to walk. */
  int deadline = 0; /**< Columns walked unlit before the walk gives up. */
};

/** @brief One moving edge's endpoint paths and frame-end world bases. */
struct ColumnEdgeFrame {
  math::Vector curves[2][3];  /**< Each endpoint's drawn path as
                                c0 + c1 t + c2 t^2: its source path, or its
                                world path under poses fixed over the frame. */
  math::Vector sources[2][3]; /**< Each endpoint's source path, when shaders
                                 read source positions. */
  const Temporal::PoseTable *poses =
      nullptr; /**< Poses by column time, or null when fixed. */
  GeodesicEdgeSampler bases[2]; /**< World bases at the frame's ends. */
  math::Vector ends[2];         /**< World end vertices at the frame's ends. */
  float pieces[2] = {}; /**< Edge length in screen steps at the frame's ends. */
  float deviation = 0.0f; /**< Bound on a drawn endpoint's distance from the
                             chord between its frame-end world positions. */
  Temporal::LongitudeTrack tracks[2]; /**< Longitude models of both ends. */
  bool known = false; /**< plan() bounded both tracks and the arc's side. */

  /** @brief Normalizes a near-unit vector.
   * @param value Nonzero vector.
   * @return Unit vector with a refined reciprocal-square-root estimate. */
  HS_HOT_INLINE static math::Vector normalize(const math::Vector &value) {
    const float norm = math::dot(value, value);
    const float estimate = screen_rsqrt(norm);
    return value * (estimate * (1.5f - 0.5f * norm * estimate * estimate));
  }

  /** @brief Builds the orthonormal basis of the arc between two directions.
   * @param a,b Nonzero endpoint directions.
   * @param total Arc length, radians.
   * @return Unit start, unit tangent toward @p b, and @p total; a collapsed
   * edge keeps a zero tangent. */
  HS_HOT_INLINE static GeodesicEdgeSampler
  basis(const math::Vector &a, const math::Vector &b, float total) {
    const math::Vector v1 = normalize(a);
    const math::Vector tangent = b - v1 * math::dot(v1, b);
    const float norm = math::dot(tangent, tangent);
    return {v1, norm > 1e-12f ? normalize(tangent) : math::Vector{}, total};
  }

  /** @brief Finds where an arc comes closest to a pole.
   * @param a,b Endpoint directions.
   * @param total Arc length, radians.
   * @param basis Receives the arc's orthonormal basis.
   * @param turn Receives the arc's point at the returned fraction.
   * @return The arc fraction of its height's extreme, or one when the extreme
   * lies outside the arc.
   * @details Longitude along a great circle is monotone and turns by half a
   * turn between its highest and lowest points, so the arc's two sides of
   * the extreme each span less than half the longitudes. The height has an
   * extreme inside the arc exactly when its slopes at the two ends,
   * b.y - a.y (a.b) and -(a.y - b.y (a.b)), have opposite signs. */
  HS_HOT_FLASH_MEMBER static float
  split_point(const math::Vector &a, const math::Vector &b, float total,
              GeodesicEdgeSampler &basis, math::Vector &turn) {
    basis = ColumnEdgeFrame::basis(a, b, total);
    const math::Vector end = normalize(b);
    const float cosine_ab = math::dot(basis.v1, end);
    float cut = 1.0f;
    if ((end.y - basis.v1.y * cosine_ab) * (basis.v1.y - end.y * cosine_ab) >
        0.0f) {
      float apex = atan2f(basis.v_perp.y, basis.v1.y);
      if (apex < 0.0f)
        apex += math::PI_F;
      if (apex > 0.0f && apex < total)
        cut = apex / total;
    }
    float sine, cosine;
    math::poly_sincosf_0_pi(cut * total, sine, cosine);
    turn = basis.v1 * cosine + basis.v_perp * sine;
    return cut;
  }

  /** @brief The geodesic basis of an arc between unit directions.
   * @param a,b Unit endpoints.
   * @return Start, unit tangent toward @p b and arc length; a collapsed arc
   * keeps a zero tangent. */
  HS_HOT_INLINE static GeodesicEdgeSampler arc(const math::Vector &a,
                                               const math::Vector &b) {
    const auto span = make_geodesic_edge_span(a, b);
    return {a, span.have_axis ? math::cross(span.axis, a) : math::Vector{},
            span.total};
  }

  /** @brief Lattice intervals the walk spaces along an edge.
   * @param pieces Edge length in screen steps.
   * @return The length rounded, at least one. */
  HS_HOT_INLINE static int lattice(float pieces) {
    return std::max(1, int(pieces + 0.5f));
  }

  /** @brief The lattice an edge keeps over the whole frame.
   * @tparam W,H Display dimensions.
   * @param starts Unit endpoints at the frame's start.
   * @param ends Unit endpoints at the frame's end.
   * @param most Largest count of interest.
   * @param arcs Receive the edge's arcs at the frame's end and, for a count
   * of at most @p most, its start.
   * @return The walk's lattice intervals when the edge's length rounds to
   * the same count, at most @p most, at both ends of the frame, and so at
   * every time between; zero otherwise. */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER static int
  steady_lattice(const math::Vector (&starts)[2], const math::Vector (&ends)[2],
                 int most, GeodesicEdgeSampler (&arcs)[2]) {
    arcs[1] = arc(ends[0], ends[1]);
    const int count = lattice(pieces_of<W, H>(arcs[1]));
    if (count > most)
      return 0;
    arcs[0] = arc(starts[0], starts[1]);
    return lattice(pieces_of<W, H>(arcs[0])) == count ? count : 0;
  }

  /** @brief A lattice point of an arc, as the walk places it.
   * @param basis The arc.
   * @param s Arc fraction in (0, 1).
   * @return The unit point at @p s. */
  HS_HOT_INLINE static math::Vector
  lattice_point(const GeodesicEdgeSampler &basis, float s) {
    float sine, cosine;
    math::poly_sincosf_0_pi(basis.total_dist * s, sine, cosine);
    return newton_unit(basis.v1 * cosine + basis.v_perp * sine);
  }

  /** @brief Measures an edge's length in screen steps.
   * @tparam W,H Display dimensions.
   * @param basis Orthonormal geodesic basis.
   * @return Arc length over the screen step at the arc's middle. */
  template <int W, int H>
  HS_HOT_INLINE static float pieces_of(const GeodesicEdgeSampler &basis) {
    const auto middle = basis.accurate(0.5f);
    return basis.total_dist / screen_step<W, H>(middle.pos, middle.tan,
                                                math::RADIANS_PER_COLUMN<W>);
  }

  /** @brief Evaluates a path held as polynomial coefficients.
   * @param curve Coefficients c0, c1, c2.
   * @param t Frame fraction.
   * @return c0 + c1 t + c2 t^2. */
  HS_HOT_INLINE static math::Vector at(const math::Vector (&curve)[3],
                                       float t) {
    return curve[0] + (curve[1] + curve[2] * t) * t;
  }

  /** @brief Prepares one edge's paths, frame-end bases and deviation bound.
   * @tparam W,H Display dimensions.
   * @param motion Vertex captures and interpolants.
   * @param from,to Start and end vertex.
   * @param table The world stages' poses by column time, or null without
   * world stages.
   * @param spin Rotation bound of the poses over the frame, radians.
   * @param source Whether to keep the source paths.
   * @details Each path is the parabola through its start, middle and end:
   * the middle capture for Temporal::Interpolant::QUADRATIC, the arc's middle
   * for Temporal::Interpolant::ARC, and the chord's for
   * Temporal::Interpolant::CHORD, which makes it the chord. It lies within
   * |middle - (start + end) / 2| of its chord. Under poses fixed over the
   * frame the paths are rotated once. Moving poses turning at most spin
   * radians per frame about a point moving along a path that departs from its
   * chord by d move a drawn endpoint off its world chord by d plus an eighth
   * of the second-derivative bound spin (spin + 2 travel), travel bounding
   * the path's speed by its chord plus 4 d; captures add their rounding. */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER void
  prepare(const Temporal::VertexMotion<> &motion, uint16_t from, uint16_t to,
          const Temporal::PoseTable *table, float spin, bool source) {
    poses = table && !table->stationary() ? table : nullptr;
    const bool fixed = table && !poses;
    const math::RotationMatrix pose =
        fixed ? table->at(table->times()) : math::RotationMatrix();
    float bend = 0.0f, travel = 0.0f;
    math::Vector world[2][2];
    for (int k = 0; k < 2; ++k) {
      const uint16_t vertex = k ? to : from;
      math::Vector start = motion.start(vertex), end = motion.end(vertex);
      const math::Vector chord_middle = (start + end) * 0.5f;
      const auto kind = motion.interpolant(vertex);
      math::Vector middle = kind == Temporal::Interpolant::CHORD ? chord_middle
                            : kind == Temporal::Interpolant::QUADRATIC
                                ? motion.middle(vertex)
                                : motion.path(vertex).middle;
      const math::Vector offset = middle - chord_middle;
      bend = fmaxf(bend, math::dot(offset, offset));
      travel = fmaxf(travel, math::distance_squared(start, end));
      if (source) {
        sources[k][0] = start;
        sources[k][1] = middle * 4.0f - start * 3.0f - end;
        sources[k][2] = (start + end - middle * 2.0f) * 2.0f;
      }
      if (fixed) {
        start = pose.apply(start);
        middle = pose.apply(middle);
        end = pose.apply(end);
      }
      curves[k][0] = start;
      curves[k][1] = middle * 4.0f - start * 3.0f - end;
      curves[k][2] = (start + end - middle * 2.0f) * 2.0f;
      if (poses) {
        start = poses->at(0).apply(start);
        end = poses->at(poses->times()).apply(end);
      }
      world[0][k] = normalize(start);
      world[1][k] = normalize(end);
    }
    for (int k = 0; k < 2; ++k) {
      const math::Vector &a = world[k][0], &b = world[k][1];
      const auto span = make_geodesic_edge_span(a, b);
      bases[k] = {a,
                  span.have_axis ? math::cross(span.axis, a) : math::Vector{},
                  span.total};
      ends[k] = b;
      pieces[k] = pieces_of<W, H>(bases[k]);
    }
    constexpr float CAPTURE_ROUNDING = 3e-5f;
    deviation = sqrtf(bend) + CAPTURE_ROUNDING;
    if (poses)
      deviation +=
          0.125f * spin * (spin + 2.0f * (sqrtf(travel) + 4.0f * deviation));
  }

  /** @brief Evaluates both endpoints at a column's time.
   * @param owner Column.
   * @param t The column's frame fraction.
   * @param a,b Receive the drawn world endpoints, unnormalized.
   * @param source Whether to evaluate the source endpoints too.
   * @param source_a,source_b Receive the source endpoints when asked. */
  HS_HOT_INLINE void at(int owner, float t, math::Vector &a, math::Vector &b,
                        bool source, math::Vector &source_a,
                        math::Vector &source_b) const {
    a = at(curves[0], t);
    b = at(curves[1], t);
    if (source) {
      source_a = at(sources[0], t);
      source_b = at(sources[1], t);
    }
    if (poses) {
      const math::RotationMatrix pose = poses->column(owner);
      a = pose.apply(a);
      b = pose.apply(b);
    }
  }

  /** @brief Projected endpoint motion over the frame.
   * @tparam W,H Display dimensions.
   * @param motion Receives start x, start y, end x and end y velocities in
   * pixels per frame. */
  template <int W, int H> HS_HOT_INLINE void velocities(float motion[4]) const {
    motion[0] = tracks[0].travel;
    motion[1] =
        y_to_screen_row<H>(bases[1].v1.y) - y_to_screen_row<H>(bases[0].v1.y);
    motion[2] = tracks[1].travel;
    motion[3] = y_to_screen_row<H>(ends[1].y) - y_to_screen_row<H>(ends[0].y);
  }

  /** @brief Bounds the columns each arm pass can light.
   * @tparam W Display width.
   * @param direction Walk direction.
   * @param segments Receives up to one walk per arm pass.
   * @return Number of walks.
   * @details Each column of an arm pass is drawn at its own time, so under the
   * linear track model of Temporal::model_track() the sweep crosses each
   * endpoint's longitude at most once per pass, at a column known within the
   * track's guard. A geodesic edge is at most a half circle and so covers the
   * shorter side between its endpoints' longitudes. A pass's lit columns
   * begin at its first column when that column lies on the arc's side and
   * otherwise at the earliest guarded crossing, and an unlit walk gives up
   * after the latest one; a pass with neither lights nothing. Tracks the
   * model cannot bound, and endpoints whose arc can change sides over the
   * frame, walk the whole pass. */
  template <int W>
  HS_HOT_FLASH_MEMBER int plan(int direction, ColumnSegment segments[2]) {
    constexpr int PASS = W / 2;
    constexpr float HALF = W / 2.0f;
    const bool start_known = Temporal::model_track<W>(
        bases[0].v1, bases[1].v1, deviation, 1.0f, tracks[0]);
    const bool end_known =
        Temporal::model_track<W>(ends[0], ends[1], deviation, 1.0f, tracks[1]);
    known = start_known && end_known;
    float span = 0.0f, widening = 0.0f;
    if (known) {
      span = Temporal::centered_remainder(tracks[1].origin - tracks[0].origin,
                                          float(W));
      widening = tracks[1].travel - tracks[0].travel;
      known = fmaxf(fabsf(span), fabsf(span + widening)) <
              HALF - tracks[0].guard - tracks[1].guard;
    }
    float entry[2] = {std::numeric_limits<float>::infinity(),
                      std::numeric_limits<float>::infinity()};
    float last[2] = {-entry[0], -entry[0]};
    if (known)
      for (const auto &track : tracks)
        Temporal::track_passes<W>(
            track, 0.0f, 1.0f, true, true,
            [&](const Temporal::ArmPass &pass) HS_HOT_INLINE {
              const int arm = pass.arm & 1;
              const float at = pass.time * HALF;
              const float walk = direction > 0 ? at : HALF - 1.0f - at;
              entry[arm] = fminf(entry[arm], walk - track.guard);
              last[arm] = fmaxf(last[arm], walk + track.guard);
            },
            track.guard);
    int count = 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int first = pass * PASS + (direction > 0 ? 0 : PASS - 1);
      if (!known) {
        segments[count++] = {first, PASS, PASS};
        continue;
      }
      int start = 0, deadline = PASS;
      if (entry[pass] > 0.0f) {
        const float time = direction > 0 ? 0.0f : (HALF - 1.0f) / HALF;
        const float a = tracks[0].origin + tracks[0].travel * time;
        const float width = span + widening * time;
        float into = width < 0.0f ? a - float(first) : float(first) - a;
        into -= float(W) * floorf(into * (1.0f / W));
        if (into <= fabsf(width)) {
          deadline = 1;
        } else {
          if (!(entry[pass] < HALF))
            continue;
          start = int(entry[pass]);
          deadline = int(ceilf(last[pass])) + 1 - start;
        }
      } else if (last[pass] < HALF) {
        deadline = int(ceilf(last[pass])) + 1;
      }
      segments[count++] = {math::fast_wrap(first + direction * start, W),
                           PASS - start, deadline};
    }
    return count;
  }
};

/** @brief Limits a walk to a column arc, keeping its hull when the arc's gap
 * lies inside the walk.
 * @tparam W Display width.
 * @param walk Walk to limit.
 * @param direction Walk direction.
 * @param begin,length Column arc, ascending from begin; length W or more
 * leaves the walk unchanged. */
template <int W>
HS_HOT_INLINE void clamp_walk(ColumnSegment &walk, int direction, int begin,
                              int length) {
  if (length >= W)
    return;
  const int head =
      direction > 0 ? begin : math::fast_wrap(begin + length - 1, W);
  const int into = math::fast_wrap(direction * (walk.start - head), W);
  if (into < length) {
    if (walk.count <= W - into)
      walk.count = std::min(walk.count, length - into);
    return;
  }
  const int gap = W - into;
  if (gap >= walk.count) {
    walk.count = 0;
    return;
  }
  walk.start = head;
  walk.count = std::min(walk.count - gap, length);
  walk.deadline -= gap;
}

/** @brief One edge's column walk: what its plain and split walks share.
 * @tparam W Display width. */
template <int W> struct ColumnWalk {
  PipelineRef &pipeline;         /**< Destination pipeline. */
  Canvas &canvas;                /**< Frame context. */
  const Fragments &points;       /**< The edge's endpoint registers. */
  TimedShaderFn shader;          /**< Time-aware source shader. */
  const ColumnEdgeFrame &frame;  /**< The edge's paths and bases. */
  const ColumnSegment *segments; /**< Walks per arm pass. */
  int segment_count;             /**< Number of walks. */
  const ClipRegion::XClip &band; /**< Render band. */
  int direction;                 /**< Walk direction. */
  float radians;                 /**< Maximum world displacement. */
  bool single_time;      /**< End columns splat at their own time too. */
  const float *velocity; /**< Endpoint motion, pixels per frame. */
  bool source_positions; /**< Shaders read source-space positions. */
  bool chained;          /**< The end vertex is left unplotted. */
  bool turning;    /**< The arc can turn its longitude direction this frame. */
  float top_y;     /**< Highest emitted unit-sphere y. */
  float bottom_y;  /**< Lowest emitted unit-sphere y. */
  float half_sine; /**< Sine of half a column. */
  float half_cosine; /**< Cosine of half a column. */
};

/** @brief Walks one edge's planned columns and plots their lattice samples.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @tparam SPLIT Walk the arc as its two sides of its closest approach to a
 * pole, each with its own column tests.
 * @param state The edge's walk.
 * @details A turning arc walks every planned column; otherwise each pass stops
 * where its lit run ends. */
template <int W, int H, RasterConfig Cfg, bool SPLIT>
HS_HOT_INLINE void walk_edge_columns(const ColumnWalk<W> &state) {
  auto &[pipeline, canvas, points, shader, frame, segments, segment_count, band,
         direction, radians, single_time, velocity, source_positions, chained,
         turning, top_y, bottom_y, half_sine, half_cosine] = state;
  const auto arm_change = [&](int from) HS_HOT_INLINE {
    return (direction > 0 ? math::fast_wrap(from + 1, W) : from) % (W / 2) == 0;
  };
  size_t sample_count = 0;
  int lattice_count = 0;
  float lattice_step = 1.0f;
  for (int segment = 0; segment < segment_count; ++segment) {
    const ColumnSegment &walk = segments[segment];
    int x = walk.start;
    bool seen = false;
    float carried = -1.0f;
    for (int column = 0; column < walk.count; ++column) {
      Temporal::count_geometry<&Temporal::GeometryCounters::columns_visited>();
      const int owner = x;
      if (radians != 0.0f && band.clipped(owner)) {
        const int skip = direction > 0
                             ? math::fast_wrap(band.rs - owner, W)
                             : math::fast_wrap(owner - band.re + 1, W);
        column += skip - 1;
        x = math::fast_wrap(owner + direction * skip, W);
        carried = -1.0f;
        continue;
      }
      x = math::fast_wrap(owner + direction, W);
      float t =
          radians == 0.0f ? 0.0f : Temporal::column_fraction<W>(float(owner));
      math::Vector a, b, source_a, source_b;
      GeodesicEdgeSampler sampler;
      float pieces, total;
      bool unit = false;
      if (radians != 0.0f) {
        frame.at(owner, t, a, b, source_positions, source_a, source_b);
        total = frame.bases[0].total_dist +
                (frame.bases[1].total_dist - frame.bases[0].total_dist) * t;
        pieces = frame.pieces[0] + (frame.pieces[1] - frame.pieces[0]) * t;
      } else {
        sampler = frame.bases[1];
        a = sampler.v1;
        b = frame.ends[1];
        total = sampler.total_dist;
        pieces = frame.pieces[1];
        unit = true;
        source_a = ColumnEdgeFrame::at(frame.sources[0], 1.0f);
        source_b = ColumnEdgeFrame::at(frame.sources[1], 1.0f);
      }
      Temporal::count_geometry<
          &Temporal::GeometryCounters::column_path_evaluations>();
      bool lit = radians == 0.0f, finished = false;
      size_t emitted = 0;
      float cut = 1.0f;
      math::Vector turn;
      if (SPLIT) {
        cut = ColumnEdgeFrame::split_point(a, b, total, sampler, turn);
        unit = true;
      }
      for (int part = 0; part < (SPLIT ? 2 : 1); ++part) {
        float s_in = 0.0f, s_out = 1.0f;
        bool start_column = true, end_column = true;
        math::Vector part_a = a, part_b = b;
        float from = 0.0f, to = 1.0f, length = total;
        GeodesicEdgeSampler piece;
        bool pieced = false;
        if (SPLIT) {
          from = part ? cut : 0.0f;
          to = part ? 1.0f : cut;
          if (!(to > from))
            break;
          part_a = part ? turn : a;
          part_b = part ? b : turn;
          length = (to - from) * total;
        }
        if (radians != 0.0f) {
          const float inverse_length = length > 1e-7f ? 1.0f / length : 0.0f;
          const float column_cosine = math::TrigLUT<W, H>::cos_theta(owner);
          const float column_sine = math::TrigLUT<W, H>::sin_theta[owner];
          const float sine = direction * half_sine;
          const math::Vector normal(
              -column_sine * half_cosine - column_cosine * sine, 0.0f,
              column_cosine * half_cosine - column_sine * sine);
          const math::Vector back(
              -column_sine * half_cosine + column_cosine * sine, 0.0f,
              column_cosine * half_cosine + column_sine * sine);
          const float end_side = math::dot(part_b, normal);
          const float back_side = math::dot(part_b, back);
          bool crossing = false;
          if (math::dot(part_a, normal) * end_side <= 0.0f) {
            if (!unit) {
              sampler = ColumnEdgeFrame::basis(a, b, total);
              unit = true;
            }
            if (SPLIT) {
              piece = ColumnEdgeFrame::basis(part_a, part_b, length);
              pieced = true;
            }
            const GeodesicEdgeSampler basis = SPLIT ? piece : sampler;
            const math::Vector along(
                column_cosine * half_cosine - column_sine * sine, 0.0f,
                column_sine * half_cosine + column_cosine * sine);
            crossing = column_basis_crossing(
                math::dot(basis.v1, normal), end_side,
                math::dot(basis.v_perp, normal), math::dot(basis.v1, along),
                math::dot(basis.v_perp, along), inverse_length, s_out);
          }
          end_column = !crossing && direction * end_side <= 0.0f &&
                       direction * back_side > 0.0f &&
                       part_b.x * column_cosine + part_b.z * column_sine > 0.0f;
          if (!crossing && !end_column)
            continue;
          if (!unit) {
            sampler = ColumnEdgeFrame::basis(a, b, total);
            unit = true;
          }
          if (SPLIT && !pieced)
            piece = ColumnEdgeFrame::basis(part_a, part_b, length);
          const GeodesicEdgeSampler basis = SPLIT ? piece : sampler;
          if (carried >= 0.0f) {
            start_column = false;
            s_in = fminf(carried, s_out);
          } else {
            const float start_back = math::dot(basis.v1, back);
            const math::Vector behind(
                column_cosine * half_cosine + column_sine * sine, 0.0f,
                column_sine * half_cosine - column_cosine * sine);
            start_column =
                start_back * back_side > 0.0f ||
                !column_basis_crossing(
                    start_back, back_side, math::dot(basis.v_perp, back),
                    math::dot(basis.v1, behind),
                    math::dot(basis.v_perp, behind), inverse_length, s_in);
            s_in = start_column ? 0.0f : fminf(s_in, s_out);
          }
          carried = crossing && !SPLIT && !arm_change(owner) ? s_out : -1.0f;
          lit = true;
          if (crossing)
            Temporal::count_geometry<
                &Temporal::GeometryCounters::column_crossings>();
        }
        const int count = ColumnEdgeFrame::lattice(pieces);
        if (count != lattice_count) {
          lattice_count = count;
          lattice_step = 1.0f / float(count);
        }
        bool edge_start = start_column, edge_end = end_column;
        if (SPLIT) {
          s_in = start_column ? from : from + (to - from) * s_in;
          s_out = end_column ? to : from + (to - from) * s_out;
          edge_start = start_column && from == 0.0f;
          edge_end = end_column && to == 1.0f;
        }
        const int first = edge_start ? 0 : int(ceilf(s_in * float(count)));
        const int last = edge_end ? count - int(chained)
                                  : int(ceilf(s_out * float(count))) - 1;
        const bool plain = single_time || (!edge_start && !edge_end);
        GeodesicEdgeSampler source;
        if (source_positions)
          source = ColumnEdgeFrame::basis(source_a, source_b, total);
        for (int k = first; k <= last; ++k) {
          Temporal::count_geometry<
              &Temporal::GeometryCounters::column_generated_samples>();
          if (++sample_count > 2 * rasterize_step_budget<W>())
            break;
          const float s = k >= count ? 1.0f : float(k) * lattice_step;
          float sine = 0.0f, cosine = 1.0f;
          if (k > 0 && k < count)
            math::poly_sincosf_0_pi(total * s, sine, cosine);
          const auto world =
              k == 0 ? sampler.v1
              : k >= count
                  ? ColumnEdgeFrame::normalize(b)
                  : newton_unit(sampler.v1 * cosine + sampler.v_perp * sine);
          if (world.y > top_y || world.y < bottom_y)
            continue;
          Fragment f;
          if constexpr (Cfg.interpolate_registers)
            f = Fragment::lerp_registers(points[0], points[1], s);
          Temporal::count_geometry<
              &Temporal::GeometryCounters::column_projections>();
          const auto p = math::vector_to_pixel<W, H>(world);
          const float arm_distance =
              fabsf(p.x - (W / 2.0f) * floorf(p.x * (2.0f / W) + 0.5f));
          if (radians == 0.0f)
            t = Temporal::column_fraction<W>(p.x);
          ++emitted;
          if (plain)
            Temporal::count_geometry<
                &Temporal::GeometryCounters::column_single_samples>();
          else
            Temporal::count_geometry<
                &Temporal::GeometryCounters::column_tap_samples>();
          if (!source_positions)
            f.pos = world;
          else
            f.pos =
                k == 0 ? source.v1
                : k >= count
                    ? ColumnEdgeFrame::normalize(source_b)
                    : newton_unit(source.v1 * cosine + source.v_perp * sine);
          shader(f.pos, f, t);
          Temporal::count_geometry<
              &Temporal::GeometryCounters::column_shader_calls>();
          Temporal::count_geometry<
              &Temporal::GeometryCounters::column_samples>();
          const bool arm_boundary = radians == 0.0f && arm_distance <= 1.0f;
          if (plain && !arm_boundary) {
            pipeline.plot_projected_at(canvas, p, f.color.color, f.age,
                                       f.color.alpha, t, nullptr);
            continue;
          }
          auto boundary_at = [&, s](float at, math::PixelCoords &projected,
                                    Pixel &color, float &age,
                                    float &alpha) HS_HOT_INLINE {
            Temporal::count_geometry<
                &Temporal::GeometryCounters::column_boundary_samples>();
            const float shift = at - t;
            const float x =
                p.x + (velocity[0] + (velocity[2] - velocity[0]) * s) * shift;
            projected = {x < 0.0f ? x + W
                         : x >= W ? x - W
                                  : x,
                         p.y + (velocity[1] + (velocity[3] - velocity[1]) * s) *
                                   shift};
            if (!arm_boundary) {
              color = f.color.color;
              age = f.age;
              alpha = f.color.alpha;
              return true;
            }
            Fragment fragment;
            if constexpr (Cfg.interpolate_registers)
              fragment = Fragment::lerp_registers(points[0], points[1], s);
            shader(f.pos, fragment, at);
            color = fragment.color.color;
            age = fragment.age;
            alpha = fragment.color.alpha;
            return true;
          };
          Temporal::BoundaryEvaluator boundary(boundary_at);
          boundary.force_timing = true;
          boundary.all_passes = arm_boundary;
          boundary.exact_at_time = true;
          Temporal::ColumnPrediction prediction{
              radians == 0.0f ? p.x : float(owner),
              std::numeric_limits<float>::infinity()};
          prediction.boundary = &boundary;
          pipeline.plot_projected_at(canvas, p, f.color.color, f.age,
                                     f.color.alpha, t, &prediction);
        }
        finished = finished || (end_column && to == 1.0f);
      }
      if (!lit) {
        carried = -1.0f;
        if (!turning && (seen || column + 1 >= walk.deadline))
          break;
        continue;
      }
      seen = true;
      if (emitted > 1)
        Temporal::count_geometry<
            &Temporal::GeometryCounters::column_span_fills>();
      if (sample_count > 2 * rasterize_step_budget<W>())
        break;
      if (finished && !turning)
        break;
    }
  }
}

/** @brief walk_edge_columns() for split arcs, out of line.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param state The edge's walk. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void walk_split_edge_columns(const ColumnWalk<W> &state) {
  walk_edge_columns<W, H, Cfg, true>(state);
}

/** @brief Rasterizes mesh edges one destination column at a time.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Register interpolation configuration.
 * @param pipeline Destination pipeline; its world stages are @p poses.
 * @param canvas Frame context.
 * @param points The walk's two endpoint fragments, written by
 * ColumnEdgeBatch::prepare.
 * @param shader Time-aware source shader.
 * @param batch Edges and their vertices' paths.
 * @param poses The world stages' poses by column time, or null without world
 * stages.
 * @param radians Maximum world displacement over the frame.
 * @param source_positions Whether shaders read source-space positions.
 * @details Each arm pass walks, in the edge's travel order, the render-band
 * columns its plan bounds, each with the edge's endpoints evaluated on their
 * paths at the column's time and rotated by its pose, and stops where the lit
 * run ends. Samples sit at arc fractions k/N, where N is the edge's length in
 * screen steps at that time, and a column plots the lattice points between
 * its two boundary meridians unweighted through the filter chain, so a moving
 * edge samples each column as the frozen edge at that column's time does. The
 * start vertex is lattice point zero; each edge also plots its end vertex
 * unless the batch is chained.
 * Columns strictly inside an edge splat at their own time; the end columns of
 * an edge moving faster than COLUMN_SINGLE_TAP_ERROR evaluate each neighbor
 * tap at that column's time. A stationary draw walks the lattice once. An
 * edge whose plan cannot bound it, such as one that can pass near a pole, is
 * walked as the two sides of its closest approach to the pole, out of line. */
template <int W, int H, RasterConfig Cfg>
FASTRUN HS_O3_FN HS_NOINLINE_NOCLONE void
rasterize_column_edges(PipelineRef &pipeline, Canvas &canvas, Fragments &points,
                       TimedShaderFn shader, const ColumnEdgeBatch &batch,
                       const Temporal::PoseTable *poses, float radians,
                       bool source_positions) {
  const float half_sine = sinf(0.5f * math::RADIANS_PER_COLUMN<W>);
  const float half_cosine = cosf(0.5f * math::RADIANS_PER_COLUMN<W>);

  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const float top_y = ColumnEdgeCull<W>::template top_of<H>(canvas);
  const float bottom_y = ColumnEdgeCull<W>::template bottom_of<H>(canvas);
  Temporal::count_geometry<&Temporal::GeometryCounters::column_draw_setups>();
  const ColumnEdgeCull<W> cull(canvas.render_x_clip(), top_y, bottom_y);
  const auto &band = canvas.render_x_clip();
  ScratchScope draw_scope(scratch_arena_a);
  const float spin = poses ? pipeline.cull_bound() : 0.0f;
  source_positions = source_positions && poses;
  auto &frame = *scratch_arena_a.make<ColumnEdgeFrame>();
  for (size_t edge = 0; edge < batch.count; ++edge) {
    const MeshEdge &pair = batch.pairs[edge];
    frame.prepare<W, H>(batch.motion, pair.u, pair.v, poses, spin,
                        source_positions);
    if (cull.rejects({frame.bases[0].v1, frame.bases[1].v1, frame.ends[0],
                      frame.ends[1]},
                     frame.deviation + 1e-4f)) {
      Temporal::count_geometry<
          &Temporal::GeometryCounters::column_edges_culled>();
      continue;
    }
    const int direction =
        frame.bases[1].v1.x * frame.bases[1].v_perp.z -
                    frame.bases[1].v1.z * frame.bases[1].v_perp.x >=
                0.0f
            ? 1
            : -1;
    ColumnSegment segments[2];
    int segment_count = 1;
    float velocity[4] = {};
    bool single_time = true;
    if (radians == 0.0f) {
      segments[0] = {0, 1, 1};
    } else {
      segment_count = frame.plan<W>(direction, segments);
      int kept = 0;
      for (int i = 0; i < segment_count; ++i) {
        ColumnSegment walk = segments[i];
        clamp_walk<W>(walk, direction, band.rs, band.length(W));
        if (walk.count > 0)
          segments[kept++] = walk;
      }
      segment_count = kept;
      if (segment_count == 0) {
        Temporal::count_geometry<
            &Temporal::GeometryCounters::column_edges_culled>();
        continue;
      }
      frame.velocities<W, H>(velocity);
      single_time = fmaxf(fmaxf(fabsf(velocity[0]), fabsf(velocity[1])),
                          fmaxf(fabsf(velocity[2]), fabsf(velocity[3]))) <=
                    COLUMN_SINGLE_TAP_ERROR * (W / 2);
    }
    Temporal::count_geometry<
        &Temporal::GeometryCounters::column_edges_prepared>();
    batch.prepare(points, batch.vertices[pair.u], batch.vertices[pair.v], edge);
    const auto sense = [&](const GeodesicEdgeSampler &basis) HS_HOT_INLINE {
      return basis.v1.x * basis.v_perp.z - basis.v1.z * basis.v_perp.x;
    };
    const bool split = radians != 0.0f && !frame.known;
    const ColumnWalk<W> walk{
        pipeline,
        canvas,
        points,
        shader,
        frame,
        segments,
        segment_count,
        band,
        direction,
        radians,
        single_time,
        velocity,
        source_positions,
        batch.chained,
        split && !(sense(frame.bases[0]) * sense(frame.bases[1]) >
                   4.0f * frame.deviation * frame.deviation),
        top_y,
        bottom_y,
        half_sine,
        half_cosine};
    if (split)
      walk_split_edge_columns<W, H, Cfg>(walk);
    else
      walk_edge_columns<W, H, Cfg, false>(walk);
  }
}

HS_O3_END

/** @brief Dispatches moving prediction or stationary time-aware sampling.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Timed pipeline.
 * @param canvas Destination frame.
 * @param points Original endpoint geometry and registers.
 * @param shader Erased time-aware shader.
 * @param opts Polyline options.
 * @param callback Borrowed geometry evaluator.
 * @param callback_bound Maximum angular callback displacement. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void
rasterize_timed_copy(PipelineRef &pipeline, Canvas &canvas,
                     const Fragments &points, TimedShaderFn shader,
                     RasterOptions opts, GeometryCallback callback,
                     float callback_bound) {
  if (points.size() < 2)
    return;
  if (callback_bound == 0.0f) {
    callback.evaluate = {};
    callback.linear_pair = false;
  }
  const float world_bound = pipeline.cull_bound();
  callback.azimuth_bound = callback.azimuth_bound && world_bound == 0.0f;
  const float bound = callback_bound + world_bound;
  if (bound == 0.0f) {
    rasterize_stationary_timed<W, H, Cfg>(pipeline, canvas, points, shader,
                                          opts, callback);
    return;
  }
  if (rasterize_dense<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                 callback, bound))
    return;
  rasterize_predicted_copy<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                      callback, bound);
}

/** @brief Walks a timed draw once per fixed pipeline copy.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Timed copy-selecting pipeline.
 * @param canvas Destination frame.
 * @param points Endpoint controls.
 * @param shader Time-aware shader.
 * @param opts Polyline options.
 * @param callback Geometry evaluator.
 * @param callback_bound Maximum callback displacement. */
template <int W, int H, RasterConfig Cfg>
HS_HOT_FLASH_MEMBER void
rasterize_timed(PipelineRef &pipeline, Canvas &canvas, const Fragments &points,
                TimedShaderFn shader, RasterOptions opts,
                const GeometryCallback &callback, float callback_bound) {
  const int copies = pipeline.copy_count();
  for (int copy = 0; copy < copies; ++copy) {
    pipeline.select_copy(copy);
    rasterize_timed_copy<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                    callback, callback_bound);
  }
  pipeline.select_copy(0);
}
