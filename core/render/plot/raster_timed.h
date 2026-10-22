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
  /** Uses two lattice positions and normalized linear interpolation. */
  bool linear_pair = false;
  /** Optional birth fraction per original control point. */
  std::span<const float> births{};

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
  TimedGeometry(const Fragments &points, GeometryCallback callback, size_t edge,
                TimedLatticeCarry &carry)
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
  GeometryCallback callback;        /**< Borrowed geometry descriptor. */
  size_t first;                     /**< First cached source vertex. */
  math::Vector *previous = nullptr; /**< Edge-local previous lattice. */
};

/** @brief Evaluated control point with time and birth registers. */
struct TimedVertex {
  Fragment fragment; /**< Geometry and shader registers. */
  float time;        /**< Capture fraction. */
  float birth;       /**< Appearance fraction. */
};

/** @brief Emits a still edge's ordinary samples at interpolated capture times.
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
    const float time = hs::clamp(a.time + (b.time - a.time) * s, 0.0f, 1.0f);
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
 * @return Conservative longitudinal displacement in columns. */
template <int W>
HS_HOT_INLINE float timed_prediction_motion(const math::Vector &world,
                                            float radians, float sine_bound) {
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
  const float motion = timed_prediction_motion<W>(result.vertex.fragment.pos,
                                                  radians, sine_bound);
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

/** @brief Replays one original edge across its transformed control spans.
 * @tparam W,H Canvas dimensions.
 * @tparam Cfg Raster configuration.
 * @param pipeline Screen emission pipeline.
 * @param canvas Frame context.
 * @param geometry Source evaluator for boundary taps.
 * @param original Original edge descriptor.
 * @param nodes Adaptive controls with cumulative material arcs.
 * @param spacing_sample Whole-edge sampler evaluated at a common frame fraction.
 * @param shader Pure source-space shader.
 * @param state Actual last landing and current arm branch.
 * @param terminal Include the final endpoint.
 * @param balanced Apply balanced sampling.
 * @param bounded Restrict emission to the part footprint. */
template <int W, int H, RasterConfig Cfg>
HS_O3_FN HS_NOINLINE_NOCLONE void walk_rigid_samples(
    PipelineRef &pipeline, Canvas &canvas, const TimedGeometry &geometry,
    const TimedControl &original, std::span<const TimedRigidNode> nodes,
    const GeodesicEdgeSampler &spacing_sample, TimedShaderFn shader,
    TimedPredictionState &state, bool terminal, bool balanced, bool bounded) {
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
  const auto emit = [&](float distance, float alpha_scale) HS_HOT_INLINE {
    sampler.seek(distance);
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
      if (geometry.birth(material) > at)
        return false;
      Fragment f;
      if constexpr (Cfg.interpolate_registers)
        f = geometry.registers(material);
      f.pos = geometry.position(material, at);
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
    prediction.record_sample = birth <= time;
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
    const size_t limit = terminal ? count : count - 1;
    for (size_t i = 0; i < limit; ++i) {
      distance += steps[i] * scale;
      emit(terminal && i + 1 == limit ? total : distance, 1.0f);
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
  constexpr float HALF = W / 2.0f;
  const float expected = state.valid ? state.column + state.pending
                                     : (state.branch + initial_time) * HALF;
  auto first =
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
      first.column + timed_prediction_advance<W>(reference_sample, 0.0f, 1.0f,
                                                 reference_span.total);
  auto last =
      resolve_rigid_control<W, H>(pipeline, canvas, geometry, b, predicted_end,
                                  state.branch, radians, sine_bound);
  ScratchScope control_scope(scratch_arena_a);
  auto *nodes = scratch_arena_a.allocate_n<TimedRigidNode>(MAX_RIGID_CONTROLS);
  size_t count = 2;
  nodes[0] = {first, a.fraction, 0.0f};
  nodes[1] = {last, b.edge == a.edge ? b.fraction : 1.0f, 0.0f};
  for (size_t i = 0; i + 1 < count;) {
    const auto &left = nodes[i].control;
    const auto &right = nodes[i + 1].control;
    const auto span = make_geodesic_edge_span(left.vertex.fragment.pos,
                                              right.vertex.fragment.pos);
    const auto middle = timed_edge_point(left.vertex.fragment.pos,
                                         right.vertex.fragment.pos, 0.5f);
    const auto projected = math::vector_to_pixel<W, H>(middle);
    const float mean_column = (left.column + right.column) * 0.5f;
    const float column = mean_column + Temporal::centered_remainder(
                                           projected.x - mean_column, float(W));
    const float time = (left.vertex.time + right.vertex.time) * 0.5f;
    const float captured = hs::clamp(column / HALF - state.branch, 0.0f, 1.0f);
    const float time_error = fabsf(captured - time);
    const float motion =
        timed_prediction_motion<W>(middle, radians, sine_bound);
    const float sine = sqrtf(fmaxf(1e-8f, 1.0f - middle.y * middle.y));
    const float deformation = radians *
                              fabsf(right.vertex.time - left.vertex.time) *
                              span.total * (W / (2.0f * math::PI_F)) / sine;
    if (time_error * HALF <= 0.5f && time_error * motion <= 0.1f &&
        deformation <= 0.1f) {
      ++i;
      continue;
    }
    const float midpoint_fraction =
        (nodes[i].fraction + nodes[i + 1].fraction) * 0.5f;
    if (count == MAX_RIGID_CONTROLS || midpoint_fraction == nodes[i].fraction ||
        midpoint_fraction == nodes[i + 1].fraction) {
      Temporal::count_geometry<
          &Temporal::GeometryCounters::rigid_fallback_edges>();
      return false;
    }
    TimedControl material = a;
    material.fraction = midpoint_fraction;
    auto control =
        resolve_rigid_control<W, H>(pipeline, canvas, geometry, material,
                                    column, state.branch, radians, sine_bound);
    for (size_t j = count; j > i + 1; --j)
      nodes[j] = nodes[j - 1];
    nodes[i + 1] = {control, material.fraction, 0.0f};
    ++count;
    Temporal::count_geometry<&Temporal::GeometryCounters::refined_controls>();
  }
  const float spacing_time = (first.vertex.time + last.vertex.time) * 0.5f;
  const auto spacing_a = pipeline.world_at(
      canvas, geometry.position(a, spacing_time), spacing_time);
  const auto spacing_b = pipeline.world_at(
      canvas, geometry.position(b, spacing_time), spacing_time);
  Temporal::count_geometry<&Temporal::GeometryCounters::common_time_edges>();
  Temporal::count_geometry<&Temporal::GeometryCounters::control_evaluations>(2);
  const auto spacing_span = make_geodesic_edge_span(spacing_a, spacing_b);
  const GeodesicEdgeSampler spacing_sample{
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
  walk_rigid_samples<W, H, Cfg>(pipeline, canvas, geometry, a,
                                std::span<const TimedRigidNode>(nodes, count),
                                spacing_sample, shader, state, terminal,
                                balanced, bounded);
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
    TimedShaderFn shader, RasterOptions opts, GeometryCallback callback,
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
    const float motion_a =
        timed_prediction_motion<W>(world_a, radians, sine_bound);
    const float motion_b =
        timed_prediction_motion<W>(world_b, radians, sine_bound);
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
        TimedVertex va{points[edge % points.size()], pa.x / HALF - branch,
                       geometry.birth(a)};
        TimedVertex vb{points[(edge + 1) % points.size()], pb.x / HALF - branch,
                       geometry.birth(b)};
        Temporal::count_geometry<&Temporal::GeometryCounters::parts>();
        timed_walk_edge<W, H, Cfg>(pipeline, canvas, va, vb, shader, terminal,
                                   balanced,
                                   stationary_cut ? &branch : nullptr);
      }
      state.valid = false;
      rigid_carry.valid = false;
      continue;
    }
    const auto endpoint_probe = [&](TimedControl &point) HS_HOT_INLINE {
      point.screen =
          pipeline.probe(canvas, geometry.position(point, 1.0f), 1.0f);
      point.motion =
          Temporal::MotionBound{radians}.columns<W, H>(point.screen.y);
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
                GeometryCallback callback, float radians) {
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
  const float bound = callback_bound + pipeline.cull_bound();
  if (bound > 0.0f) {
    if (rasterize_dense<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                   callback, bound))
      return;
    rasterize_predicted_copy<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                        callback, bound);
    return;
  }
  rasterize_stationary_timed<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                        callback);
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
                GeometryCallback callback, float callback_bound) {
  const int copies = pipeline.copy_count();
  for (int copy = 0; copy < copies; ++copy) {
    pipeline.select_copy(copy);
    rasterize_timed_copy<W, H, Cfg>(pipeline, canvas, points, shader, opts,
                                    callback, callback_bound);
  }
  pipeline.select_copy(0);
}
