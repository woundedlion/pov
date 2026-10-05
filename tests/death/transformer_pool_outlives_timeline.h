/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_death.h.

/**
 * @brief Death case: a pool outliving its Timeline must trap.
 * @details Transformer surface — the destructor reaches back into the timeline
 *          to drop the pool's clear hook, and the spawned completion callbacks
 *          reach back into the pool, so the two lifetimes are ordered. An owner
 *          that declares them the other way gets a dead reference here rather
 *          than at some later step().
 */
inline void case_transformer_pool_outlives_timeline() {
  configure_arenas_default();
  alignas(Timeline) static uint8_t tl_storage[sizeof(Timeline)];
  Timeline *tl = new (tl_storage) Timeline();
  RippleTransformer<2> rt(*tl);
  rt.init_storage(persistent_arena);
  tl->~Timeline();
  // ~RippleTransformer at scope exit -> HS_CHECK
}

/**
 * @brief Concrete Effect for the canvas death cases.
 * @details Defaults to 32x16 and exposes register_param via reg.
 */
struct DeathEffect : public Effect {
  using Effect::register_param;
  /**
   * @brief Constructs the effect at the requested resolution.
   * @param width Canvas width in pixels.
   * @param height Canvas height in pixels.
   */
  DeathEffect(int width = 32, int height = 16) : Effect(width, height) {}
  /**
   * @brief Draws one frame (no-op; the death cases never render).
   */
  void draw_frame() override {}
  /**
   * @brief Registers a parameter over the unit range, exposing register_param.
   * @param n Parameter name.
   * @param p Pointer to the backing float storage.
   */
  void reg(const char *n, float *p) { register_param(n, p, 0.0f, 1.0f); }
  void reg_float_options(float *p, const char *const *options, int count) {
    register_param("options", p, 0.0f, 1.0f, false, false, options, count);
  }
  /** @brief Registers a typed enum with the requested option count. */
  template <typename Enum> void reg_enum(Enum *p, int count) {
    static constexpr const char *OPTIONS[] = {"zero"};
    register_param("enum", p, OPTIONS, nullptr, count);
  }
  /**
   * @brief Registers an integer parameter, exposing register_int_param.
   * @param n Parameter name.
   * @param p Pointer to the backing integer storage.
   * @param min Minimum value, inclusive.
   * @param max Maximum value, inclusive.
   */
  template <typename Integer>
  void reg_int(const char *n, Integer *p, int min, int max) {
    register_int_param(n, p, min, max);
  }
};

/**
 * @brief Death case: a second simultaneously-live Effect must trap.
 * @details Canvas surface — the structural twin of case_timeline_double_construct.
 *          Every Effect aliases the same two static framebuffers and double-buffer
 *          indices, so a second live instance would scribble over the first's
 *          frames; the construction guard traps instead. The real app builds the
 *          next effect only after destroying the outgoing one.
 */
inline void case_effect_double_construct() {
  DeathEffect a;
  DeathEffect b; // second live Effect ctor -> HS_CHECK(!s_alive) -> trap
  if (opaque(a.strobe_columns()))
    std::printf("x");
}

/** @brief Death case: a zero Effect width must trap. */
inline void case_effect_width_zero() { DeathEffect fx(opaque(0), 16); }

/** @brief Death case: a zero Effect height must trap. */
inline void case_effect_height_zero() { DeathEffect fx(32, opaque(0)); }

/** @brief Death case: an Effect width above MAX_W must trap. */
inline void case_effect_width_over_max() {
  DeathEffect fx(opaque(MAX_W + 1), 16);
}

/** @brief Death case: an Effect height above MAX_H must trap. */
inline void case_effect_height_over_max() {
  DeathEffect fx(32, opaque(MAX_H + 1));
}

/**
 * @brief Initializes a particle system with the requested maximum lifetime.
 * @param max_life Maximum lifetime passed to ParticleSystem::init.
 */
inline void init_particle_system_with_lifetime(float max_life) {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, max_life);
}

/** @brief Death case: a zero particle lifetime must trap. */
inline void case_particle_lifetime_zero() {
  init_particle_system_with_lifetime(opaque(0.0f));
}

inline void case_particle_friction_nan() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, opaque(std::numeric_limits<float>::quiet_NaN()));
}

inline void case_particle_gravity_nan() {
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, opaque(std::numeric_limits<float>::quiet_NaN()));
}

/** @brief Death case: a NaN particle lifetime must trap. */
inline void case_particle_lifetime_nan() {
  init_particle_system_with_lifetime(
      opaque(std::numeric_limits<float>::quiet_NaN()));
}

inline void case_random_walk_nonfinite_options() {
  math::Orientation<> orientation;
  FastNoiseLite noise;
  Animation::RandomWalkOptions options;
  options.drift = opaque(std::numeric_limits<float>::quiet_NaN());
  Animation::RandomWalk<32> walk(orientation, math::Vector(0, 0, 1), noise,
                                 options);
}

/** @brief Death case: a particle lifetime above uint16_t must trap. */
inline void case_particle_lifetime_over_max() {
  init_particle_system_with_lifetime(opaque(65536.0f));
}

/** @brief Pipeline stub for the particle-render lifetime death case. */
struct DeathPlotPipeline {
  void plot(Canvas &, const math::Vector &, const Pixel &, float, float) {}
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/** @brief Death case: a nonempty render with zero max_life must trap. */
inline void case_particle_render_zero_lifetime() {
  configure_arenas_default();
  static uint8_t buf[4096];
  Arena arena(buf, sizeof(buf));
  Animation::ParticleSystem<32, 1> ps;
  ps.init(arena, 0.85f, 0.0f, 1.0f);
  ps.spawn(math::Vector(1, 0, 0), math::Vector(), 0);
  ps.max_life = 0;

  DeathEffect fx;
  Canvas canvas(fx);
  DeathPlotPipeline pipeline;
  Plot::ParticleSystem::draw<32, 16>(pipeline, canvas, ps,
                                     [](const math::Vector &, Fragment &) {});
}

/**
 * @brief Death case: a second simultaneously-live correction guard must trap.
 * @details LED surface — NoColorCorrection and NoTempCorrection share one
 *          liveness flag and set the global FastLED correction/temperature, so a second
 *          live guard of either type would leave the wrong baseline on the earlier
 *          guard's exit; the construction guard traps instead.
 */
inline void case_correction_guard_double_construct() {
  NoColorCorrection a;
  NoColorCorrection b; // second live guard -> liveness HS_CHECK -> trap
  if (correction_guard_live() == opaque(true))
    std::printf("x");
}

/**
 * @brief Death case: a live NoColorCorrection plus a NoTempCorrection must trap.
 * @details LED surface — the two guard types share the one liveness flag, so a
 *          second live guard of the OTHER type is as unsafe as a same-type
 *          double-construct; this is the case the shared "either type" contract
 *          exists to guarantee. The construction guard traps on either.
 */
inline void case_correction_guard_cross_type() {
  NoColorCorrection a;
  NoTempCorrection b; // second live guard of a different type -> trap
  if (correction_guard_live() == opaque(true))
    std::printf("x");
}

inline void case_float_options_missing_labels() {
  DeathEffect effect;
  float value = 0;
  effect.reg_float_options(&value, nullptr, 2);
}

inline void case_float_options_missing_count() {
  DeathEffect effect;
  float value = 0;
  const char *options[] = {"zero", "one"};
  effect.reg_float_options(&value, options, 0);
}

inline void case_float_options_wrong_range() {
  DeathEffect effect;
  float value = 0;
  const char *options[] = {"zero"};
  effect.reg_float_options(&value, options, 1);
}

/**
 * @brief Death case: overflowing the fixed ParamList must trap.
 * @details Canvas surface — register_param traps rather than silently dropping a
 *          registration, which would desync the GUI and, on WASM, break the
 *          no-realloc memory-view invariant.
 */
inline void case_register_param_overflow() {
  DeathEffect fx;
  static float slot = 0.0f;
  // Distinct names, so the capacity guard fires ahead of the duplicate guard.
  constexpr int CAPACITY = static_cast<int>(Effect::ParamList::FIXED_CAPACITY);
  static char names[CAPACITY + 1][8];
  for (int i = 0; i < opaque(CAPACITY + 1); ++i) {
    std::snprintf(names[i], sizeof(names[i]), "p%d", i);
    fx.reg(names[i], &slot);
  }
}

inline void case_register_param_duplicate() {
  DeathEffect fx;
  float value = 0.5f;
  fx.reg("duplicate", &value);
  fx.reg("duplicate", &value);
}

inline void case_register_param_default_outside_range() {
  DeathEffect fx;
  float value = 2.0f;
  fx.reg("outside", &value);
}

inline void case_restore_parameters_unknown_name() {
  DeathEffect effect;
  const std::array values{
      std::pair<std::string, float>{"unknown", opaque(0.5f)}};
  effect.replay_parameter_writes(values);
}

inline void case_restore_parameters_readonly_name() {
  DeathEffect effect;
  float value = 0.0f;
  effect.register_param("readonly", &value, ParamSpec<float>{.readonly = true});
  const std::array values{
      std::pair<std::string, float>{"readonly", opaque(0.5f)}};
  effect.replay_parameter_writes(values);
}

inline void case_restore_parameters_singular_mobius() {
  reset_globals();
  MobiusGrid<32, 16> effect;
  effect.init();
  const std::array values{
      std::pair<std::string, float>{"Mobius A Re", opaque(1.0f)},
      std::pair<std::string, float>{"Mobius A Im", opaque(0.0f)},
      std::pair<std::string, float>{"Mobius B Re", opaque(1.0f)},
      std::pair<std::string, float>{"Mobius B Im", opaque(0.0f)},
      std::pair<std::string, float>{"Mobius C Re", opaque(1.0f)},
      std::pair<std::string, float>{"Mobius C Im", opaque(0.0f)},
      std::pair<std::string, float>{"Mobius D Re", opaque(1.0f)},
      std::pair<std::string, float>{"Mobius D Im", opaque(0.0f)}};
  effect.replay_parameter_writes(values);
}

/**
 * @brief Death case: an integer param bound the target cannot store must trap.
 * @details Canvas surface — a value write narrows through
 *          static_cast<Integer>(float), which is undefined once the registered
 *          range leaves the storage type.
 */
inline void case_register_int_param_range() {
  DeathEffect fx;
  static uint8_t slot = 0;
  fx.reg_int("count", &slot, 0, opaque(256));
}

inline void case_register_enum_param_range() {
  enum class Mode : uint8_t { ZERO };
  DeathEffect fx;
  Mode slot = Mode::ZERO;
  fx.reg_enum(&slot, opaque(257));
}

inline void case_register_enum_param_bound_inexact() {
  enum class Mode : uint32_t { ZERO };
  DeathEffect fx;
  Mode slot = Mode::ZERO;
  fx.reg_enum(&slot, opaque(16777218));
}

inline void case_register_int_param_max_inexact() {
  DeathEffect fx;
  static int32_t slot = 0;
  fx.reg_int("count", &slot, 0, opaque(std::numeric_limits<int32_t>::max()));
}

inline void case_register_int_param_min_inexact() {
  DeathEffect fx;
  static int32_t slot = 0;
  fx.reg_int("count", &slot, opaque(-std::numeric_limits<int32_t>::max()), 0);
}

inline void case_param_spec_uint32_bound_outside_storage() {
  DeathEffect fx;
  uint32_t value = 0;
  fx.register_param("unsigned", &value,
                    ParamSpec<uint32_t>{.min = 0, .max = 4294967296LL});
}

inline void case_param_spec_uint32_bound_inexact() {
  DeathEffect fx;
  uint32_t value = 0;
  fx.register_param("unsigned", &value,
                    ParamSpec<uint32_t>{.min = 0, .max = 4294967295LL});
}

inline void case_param_spec_integer_preserve_policy() {
  DeathEffect fx;
  uint8_t value = 0;
  fx.register_param(
      "integer", &value,
      ParamSpec<uint8_t>{.initial_value =
                             ParamInitialValue::PRESERVE_REQUESTED_FLOAT});
}

inline void case_param_spec_float_bound_nonfinite() {
  DeathEffect fx;
  float value = 0.0f;
  fx.register_param(
      "float", &value,
      ParamSpec<float>{.max = std::numeric_limits<float>::infinity()});
}

/** @brief A null parameter name traps before lookup or diagnostic formatting. */
inline void case_param_spec_name_null() {
  DeathEffect fx;
  float value = 0.0f;
  fx.register_param(nullptr, &value, ParamSpec<float>{});
}

inline void case_param_spec_requested_nonfinite() {
  DeathEffect fx;
  float value = std::numeric_limits<float>::quiet_NaN();
  fx.register_param(
      "requested", &value,
      ParamSpec<float>{.initial_value =
                           ParamInitialValue::PRESERVE_REQUESTED_FLOAT});
}

inline void case_param_spec_option_label_null() {
  DeathEffect fx;
  uint8_t value = 0;
  const char *const options[] = {"Zero", nullptr};
  fx.register_param("labels", &value,
                    ParamSpec<uint8_t>::enumerated(options, 2));
}

inline void case_param_spec_export_label_null() {
  DeathEffect fx;
  enum class Mode : uint8_t { ZERO };
  Mode value = Mode::ZERO;
  const char *const options[] = {"Zero", "One"};
  const char *const exports[] = {"Mode::ZERO", nullptr};
  fx.register_param("exports", &value,
                    ParamSpec<Mode>::enumerated(options, 2, exports));
}

/**
 * @brief Death case: set_clip rejects x_end beyond the canvas width.
 */
inline void case_set_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  DeathEffect fx;
  fx.set_clip(0, H, 0, opaque(W + 1));
}

/** @brief Death case: clip state cannot change after a frame begins. */
inline void case_set_clip_mid_frame() {
  DeathEffect fx;
  Canvas canvas(fx);
  fx.set_clip(0, fx.height(), 0, fx.width());
}

/** @brief Death case: the publication envelope rejects values above one. */
inline void case_output_envelope_out_of_range() {
  DeathEffect fx;
  fx.set_output_envelope(opaque(1.01f));
}

/**
 * @brief Death case: set_clip_x rejects x_end beyond the canvas width.
 */
inline void case_set_clip_x_out_of_bounds() {
  constexpr int W = 32;
  DeathEffect fx;
  fx.set_clip_x(0, opaque(W + 1));
}

/**
 * @brief Death case: an arc start outside [0, w) must trap.
 * @details Clip surface — arcs_overlap wraps the seam-relative offset with one
 *          conditional add instead of a modulo, which only lands in range while
 *          both starts are already reduced; a start outside the cylinder would
 *          silently report the wrong overlap. Both lengths are positive and
 *          under w so the early-out branches do not preempt the guard.
 */
inline void case_arcs_overlap_start_out_of_range() {
  bool hit = ClipRegion::arcs_overlap(opaque(-1), opaque(2), opaque(0),
                                      opaque(2), opaque(8));
  if (hit)
    std::printf("x");
}

/** @brief Death case: the coherent shader rejects a zero block edge. */
inline void case_scan_block_coherent_zero_block() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas canvas(fx);
  const math::Vector position = math::UP;
  Scan::Shader::draw_block_coherent<W, H, 1>(
      canvas, opaque(0), &position, scratch_arena_a,
      [](const math::Vector &) { return Scan::Shader::BlockCell<1>{0}; },
      [](const math::Vector &, const Scan::Shader::BlockCandidates<1> &) {
        return Color4{};
      });
}

/**
 * @brief Death case: the shader rejects a clip beyond its LUT width.
 * @details Direct construction exercises the downstream guard independently
 *          of the Effect setters.
 */
inline void case_scan_clip_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = opaque(W + 1);
  cr.y_end = H;
  cr.margin = 0;
  cr.w = W;
  cr.h = H;
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: the shader rejects a render band past the canvas rows.
 * @details The rows the guard admits subscript the canvas, so a band reaching
 *          past H must be rejected even though the phi LUT holds those rows.
 */
inline void case_scan_clip_rows_out_of_bounds() {
  constexpr int W = 32, H = 16;
  ClipRegion cr;
  cr.x_end = W;
  cr.y_end = opaque(H + 1);
  cr.margin = 0;
  cr.w = W;
  cr.h = opaque(H + 1);
  Scan::Shader::check_lut_domain<W, H>(cr);
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later Face claimed.
 * @details The second build retargets the first face's spans, so the first no
 *          longer describes its own geometry.
 */
inline void case_face_scratch_retargeted() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[6];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 6; ++i) {
    float a = (2.0f * math::PI_F * i) / 6.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
  }
  for (int i = 0; i < 3; ++i) {
    idx_a[i] = static_cast<uint16_t>(i);
    idx_b[i] = static_cast<uint16_t>(i + 3);
  }
  static SDF::FaceScratchBuffer scratch;
  SDF::Face first(std::span<const math::Vector>(verts, 6),
                  std::span<const uint16_t>(idx_a, 3), scratch, HV, H);
  SDF::Face second(std::span<const math::Vector>(verts, 6),
                   std::span<const uint16_t>(idx_b, 3), scratch, HV, H);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/** @brief Death case: a virtual-height Face must use the active display grid. */
inline void case_face_virtual_height_mismatched_geometry() {
  constexpr int H = 16;
  const math::Vector verts[] = {math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const uint16_t indices[] = {0, 1, 2};
  static SDF::FaceScratchBuffer scratch;
  SDF::Face face(verts, indices, scratch, opaque(H + hs::H_OFFSET + 1), H);
  (void)face;
}

/**
 * @brief Death case: scanning a Face whose scratch buffer a later, culled Face
 *        claimed.
 * @details The second build fills the scratch buffer and only then culls on
 *          collapsed area, so the first face's spans are retargeted even though
 *          nothing will be drawn for the second.
 */
inline void case_face_scratch_retargeted_by_culled_face() {
  constexpr int H = 16, HV = H + hs::H_OFFSET;
  math::Basis basis =
      math::make_basis(math::Quaternion(), math::Vector(0, 1, 0));
  math::Vector verts[4];
  uint16_t idx_a[3], idx_b[3];
  for (int i = 0; i < 3; ++i) {
    float a = (2.0f * math::PI_F * i) / 3.0f;
    verts[i] = (basis.v * cosf(0.6f) +
                (basis.u * cosf(a) + basis.w * sinf(a)) * sinf(0.6f))
                   .normalized();
    idx_a[i] = static_cast<uint16_t>(i);
    // Coincident vertices enclose no area: culled after the scratch fill.
    idx_b[i] = 3;
  }
  verts[3] = basis.v;
  static SDF::FaceScratchBuffer scratch;
  SDF::Face first(std::span<const math::Vector>(verts, 4),
                  std::span<const uint16_t>(idx_a, 3), scratch, HV, H);
  SDF::Face second(std::span<const math::Vector>(verts, 4),
                   std::span<const uint16_t>(idx_b, 3), scratch, HV, H);
  (void)second;
  (void)first.get_vertical_bounds<H>();
}

/**
 * @brief Death case: a scan rejects a canvas that is not its <W, H>.
 * @details Direct construction exercises the guard independently of the draw
 *          primitives that call it.
 */
inline void case_scan_canvas_dim_mismatch() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, opaque(H + 1));
  Canvas c(fx);
  Scan::check_canvas_dims<W, H>(c);
}

/**
 * @brief Death case: a scan rejects a direct-raster sink prepared for no
 *        canvas.
 * @details The sink writes through its cached framebuffer base, which the
 *          double-buffered canvas leaves pointing at the buffer being scanned
 *          out until prepare() runs.
 */
inline void case_scan_pipeline_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  Scan::check_pipeline_prepared(sink, c);
}

/**
 * @brief Death case: erasing a direct-raster sink prepared for no canvas.
 * @details PipelineRef drops prepared_for(), so a draw taking the erased handle
 *          cannot check the cached base; the erasure checks it instead.
 */
inline void case_pipeline_ref_erase_not_prepared() {
  constexpr int W = 32, H = 16;
  DeathEffect fx(W, H);
  Canvas c(fx);
  static Filter::Screen::DirectAntiAliasSink<W, H> sink;
  PipelineRef erased(sink, c);
  (void)erased;
}

/**
 * @brief Builds an invalid scan mesh with either missing offsets or a bad index.
 * @param omit_offsets Selects the missing-offset guard before index validation.
 */
inline void scan_mesh_invalid_fixture(bool omit_offsets) {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[1024];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  mesh.vertices.bind(geom, 3);
  for (int i = 0; i < 3; ++i)
    mesh.vertices.push_back(math::Vector(0.0f, 1.0f, 0.0f));
  mesh.face_counts.bind(geom, 1);
  mesh.face_counts.push_back(opaque<uint8_t>(3));
  if (!omit_offsets) {
    mesh.face_offsets.bind(geom, 1);
    mesh.face_offsets.push_back(opaque<uint16_t>(0));
  }
  mesh.faces.bind(geom, 3);
  mesh.faces.push_back(opaque<uint16_t>(0));
  mesh.faces.push_back(opaque<uint16_t>(1));
  mesh.faces.push_back(opaque<uint16_t>(3));

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Scan::Mesh::draw<W, H>(
      pipe, c, mesh, [](const math::Vector &, Fragment &) {}, scratch);
}

/** @brief Death case: absent face offsets must trap. */
inline void case_scan_mesh_missing_offsets() {
  scan_mesh_invalid_fixture(true);
}

/** @brief Death case: a face index past the vertex pool must trap. */
inline void case_scan_mesh_face_index_out_of_range() {
  scan_mesh_invalid_fixture(false);
}

/**
 * @brief Death case: a class bake naming a class the class table does not hold.
 * @details ArenaVector::operator[] only asserts, so an out-of-range class id
 *          would read arbitrary memory as a CongruenceClass and hand its LUT to
 *          the per-pixel probe. Scan::Mesh bounds the id per face.
 */
inline void case_scan_mesh_class_id_out_of_range() {
  constexpr int W = 32, H = 16;
  static uint8_t geom_buf[4096];
  static uint8_t scratch_buf[64 * 1024];
  Arena geom(geom_buf, sizeof(geom_buf));
  Arena scratch(scratch_buf, sizeof(scratch_buf));

  MeshState mesh;
  mesh.vertices.bind(geom, 3);
  mesh.vertices.push_back(math::Vector(1.0f, 0.0f, 0.0f));
  mesh.vertices.push_back(math::Vector(0.0f, 1.0f, 0.0f));
  mesh.vertices.push_back(math::Vector(0.0f, 0.0f, 1.0f));
  mesh.face_counts.bind(geom, 1);
  mesh.face_counts.push_back(static_cast<uint8_t>(3));
  mesh.face_offsets.bind(geom, 1);
  mesh.face_offsets.push_back(static_cast<uint16_t>(0));
  mesh.faces.bind(geom, 3);
  for (uint16_t v = 0; v < 3; ++v)
    mesh.faces.push_back(v);

  MeshOps::MeshClassBake bake;
  bake.classes.bind(geom, 1);
  bake.face_recs.bind(geom, 1);
  bake.face_recs.push_back({opaque<uint8_t>(0), 0, 0}); // no class founded

  DeathEffect fx(W, H);
  Canvas c(fx);
  Pipeline<W, H> pipe;
  Scan::Mesh::draw<W, H>(
      pipe, c, mesh, [](const math::Vector &, Fragment &) {}, scratch, &bake);
}

/**
 * @brief Minimal duck-typed mesh: one 2-gon face whose second index (130)
 *        exceeds the TriangularBitset<128> capacity. Shared by both the
 *        face-walk draw() and the extract_edges over-capacity death cases so the
 *        mock interface is defined once, not kept in sync across two copies.
 * @details The trap fires before any vertex or pipeline access, so the vertex
 *          store only needs to satisfy the interface.
 */
struct OverCapacityMockMesh {
  /**
   * @brief Stand-in vertex store satisfying the mesh interface.
   */
  struct Verts {
    /**
     * @brief Returns a fixed vertex for any index.
     * @param Unused vertex index.
     * @return A constant Vector{0,1,0}.
     */
    math::Vector operator[](size_t) const {
      return math::Vector{0.0f, 1.0f, 0.0f};
    }
    /**
     * @brief Reports the vertex count.
     * @return Always 1.
     */
    size_t size() const { return 1; }
  } vertices;
  uint8_t fc[1];  /**< Face-counts data: a single 2-gon face. */
  uint16_t fi[2]; /**< Face-index data; second entry is over-capacity. */
  /**
   * @brief Builds the mock mesh with one over-capacity 2-gon face.
   * @details Stores the over-capacity index at runtime so the optimizer can't
   *          prove the trap at compile time and reshape the case (see opaque).
   */
  OverCapacityMockMesh() : fc{2}, fi{0, opaque<uint16_t>(130)} {}
  /**
   * @brief Returns the face-counts array.
   * @return Pointer to the face-counts data.
   */
  const uint8_t *get_face_counts_data() const { return fc; }
  /**
   * @brief Returns the number of faces.
   * @return Always 1.
   */
  size_t get_face_counts_size() const { return 1; }
  /**
   * @brief Returns the flat face-index array.
   * @return Pointer to the face-index data.
   */
  const uint16_t *get_faces_data() const { return fi; }
  /**
   * @brief Returns the flat face-index array length.
   * @return Always 2.
   */
  size_t get_faces_size() const { return 2; }
};
