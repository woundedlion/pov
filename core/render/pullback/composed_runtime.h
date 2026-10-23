/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/** @file composed_runtime.h
 * @brief `ComposedEffect`, the lifecycle base of composed pullback effects.
 */

/**
 * @brief A complete composed effect: the shared lifecycle plus a pipeline
 *        declared by a ranked stage Spec.
 * @details Parameter and runtime storage derive from the Spec's pipeline
 * providers. Required `Derived` members: `EFFECT_ID`, `PRESET_IDS` (preset
 * identities by preset number), `PARAMETER_SCHEMA_VERSION` (changes with the
 * Params layout) and `PRESET_DWELL_FRAMES`. Presets resolve through
 * the `PRESET_ROW` hook with `PRESET_IDS`, then `PRESETS`; only single-preset
 * effects may fall back to startup params. Optional hook bindings:
 * `INITIAL_PARAMS` and `AFTER_COMPOSED_INIT`. Optional members:
 * `ANIMATED_MOBIUS`, `CAMERA_SPIN_RATE`, and `WARP_NOISE_SEED` /
 * `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` shadows. A shade() shadow changes
 * the emitted body only by calling
 * RenderPipeline::evaluate(view, frame.ctx, frame.prepared).
 * `DESCRIPTOR_DIGEST` and `PRESET_BANK_DIGEST` pin the pattern document; the
 * runtime never reads them.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this base.
 * @tparam SpecT The effect's `Spec`.
 * @tparam Ownership Whether Derived declares hooks or inherits composed defaults.
 */
template <int W, int H, typename Derived, typename SpecT,
          Control::HookOwnership Ownership = Control::HookOwnership::AUTHORED>
class ComposedEffect
    : public ChoreographedEffect<
          Derived, ParamsFor<SpecT>,
          std::conditional_t<Ownership == Control::HookOwnership::INHERITED,
                             ComposedEffect<W, H, Derived, SpecT, Ownership>,
                             Derived>>,
      private ProjectionWalkState<SpecT::ANIMATED_PROJECTION> {
  using ParamsT = ParamsFor<SpecT>;
  using Choreography = ChoreographedEffect<
      Derived, ParamsT,
      std::conditional_t<Ownership == Control::HookOwnership::INHERITED,
                         ComposedEffect, Derived>>;
  friend Choreography;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  using Params = ParamsT; ///< Keyed parameter storage derived from `SpecT`.
  using Spec = SpecT;     ///< The effect's ranked stage Spec.
  using Metadata = MetadataFor<SpecT>; ///< Topology derived from typed stages.
  /// Whether the effect samples time across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = SpecT::TEMPORAL_SAMPLING;
  /// Whether hue tables retain exact float copies for interpolation.
  static constexpr bool FLOAT_HUE_TABLES =
      TEMPORAL_SAMPLING && Metadata::HUE == HueMode::NOISE;
  static_assert(!TEMPORAL_SAMPLING || (W >= 2 && W % 2 == 0),
                "Column-time sampling requires an even two-arm canvas width");
  /// Frame context, including retained color resources when sampling time.
  using FrameState =
      std::conditional_t<TEMPORAL_SAMPLING, FrameIntervalState<ParamsT>,
                         Pullback::FrameState<ParamsT>>;
  using Binding = Pullback::Binding<FrameState>; ///< Pipeline binding.
  /// Whether the projection orientation follows a random walk and spin.
  static constexpr bool ANIMATED_PROJECTION = SpecT::ANIMATED_PROJECTION;
  /**
   * @brief Whether any resource of a kind samples a noise field.
   * @tparam Kind Resource kind to test.
   * @return True when a `Kind` resource satisfies
   *         `ComposedDetail::RESOURCE_NOISE`.
   */
  template <ResourceKind Kind> static consteval bool has_noise() {
    bool result = false;
    Params{}.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Kind &&
                    ComposedDetail::RESOURCE_NOISE<Resource>)
        result = true;
    });
    return result;
  }
  /// Whether a surface resource samples noise.
  static constexpr bool HAS_SURFACE_NOISE = has_noise<ResourceKind::SURFACE>();
  /// Whether a warp resource samples noise.
  static constexpr bool HAS_WARP_NOISE = has_noise<ResourceKind::WARP>();
  /// Whether a source resource samples noise.
  static constexpr bool HAS_SOURCE_NOISE = has_noise<ResourceKind::SOURCE>();
  /// The Spec's pipeline bound to this effect's frame.
  using RenderPipeline = typename SpecT::template Pipeline<Binding>;
  static_assert(!SpecT::FACTORED_COLUMNS || RenderPipeline::LINEAR_COLUMNS,
                "factored columns require a rotation/projection prefix");
  static_assert(Metadata::COVERAGE != ProjectionCoverageMode::EDGE_FADE ||
                    Metadata::ProjectStage::EDGE_DISTANCE_AVAILABLE,
                "edge-fade coverage requires projection edge distance");
  static_assert(Metadata::PATH_TRACKED ==
                    (Metadata::HUE == HueMode::PATH_LENGTH),
                "path length hue must match pipeline tracking");
  /// Prepared per-frame state the pipeline shades from.
  using Frame = typename RenderPipeline::Frame;

  /** @brief Per-field noise seeds; an effect shadows one to decorrelate its
      field from the shared spatial phase. */
  static constexpr int32_t WARP_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SOURCE_NOISE_SEED = EFFECT_NOISE_SEED;  ///< Source.
  static constexpr int32_t SURFACE_NOISE_SEED = EFFECT_NOISE_SEED; ///< Surface.

  /** @brief Constructs the effect at W x H with its Spec column-strobe mode. */
  HS_COLD_MEMBER ComposedEffect()
      : Choreography(W, H, {.strobe = SpecT::STROBE_COLUMNS}) {}

  /**
   * @brief Claims the runtime's persistent storage and registers the effect.
   * @details Every allocation this makes comes from `persistent_arena`, so the
   * effect heap-allocates nothing after init. The `AFTER_COMPOSED_INIT` hook
   * runs last, once the parameters, palette cycler and camera walks are live.
   */
  HS_COLD_MEMBER void init() override {
    this->begin_choreography();
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena,
                          persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    if constexpr (Metadata::HUE == HueMode::NOISE)
      Pullback::init_effect_noise(state->color_noise, HUE_NOISE_SEED);
    params.visit([&]<typename Resource>(auto &) {
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>) {
        constexpr int32_t SEED = Resource::KIND == ResourceKind::WARP
                                     ? Derived::WARP_NOISE_SEED
                                 : Resource::KIND == ResourceKind::SOURCE
                                     ? Derived::SOURCE_NOISE_SEED
                                     : Derived::SURFACE_NOISE_SEED;
        constexpr int32_t INSTANCE_SEED = [] {
          if constexpr (ComposedDetail::qualified<
                            Resource, typename Params::ResourceTypes>()) {
            constexpr uint32_t HASH = fnv1a(Resource::KEY.view());
            return static_cast<int32_t>(static_cast<uint32_t>(SEED) ^ HASH);
          } else {
            return SEED;
          }
        }();
        Pullback::init_effect_noise(
            state->resources.template get<Resource::KEY>().noise,
            INSTANCE_SEED);
      }
    });
    palette_cycler.init_generated(persistent_arena, next_palette, this,
                                  GeneratedPaletteBank::DWELL_FRAMES,
                                  GeneratedPaletteBank::FADE_FRAMES,
                                  math::ease_in_out_sin);
    if constexpr (ANIMATED_PROJECTION)
      timeline.add(0, Animation::RandomWalk<W>(
                          this->projection_walk, math::UP,
                          state->projection_walk_noise,
                          typename Animation::RandomWalk<W>::Options{},
                          PROJECTION_WALK_SEED));
    timeline.add(0, Animation::RandomWalk<W>(
                        outer_walk, math::UP, state->outer_walk_noise,
                        typename Animation::RandomWalk<W>::Options{},
                        CAMERA_WALK_SEED));
    register_parameters();
    if constexpr (Choreography::template has_hook<
                      Control::ChoreographyHook::AFTER_COMPOSED_INIT>())
      this->template invoke_hook<
          Control::ChoreographyHook::AFTER_COMPOSED_INIT>();
  }

  /**
   * @brief Advances the frame clocks, resolves the frame state and shades.
   * @details Clocks, camera walks, palette and the preset lerp all step before
   * prepare_frame() snapshots them before shading or column-time interpolation.
   */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(fx_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(fx_advance);
      step_choreography();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const FrameState current = prepare_frame();
    if constexpr (TEMPORAL_SAMPLING) {
      auto &history = static_cast<Spacetime::FrameHistory<
          FrameState, Metadata::HUE, TEMPORAL_SAMPLING, FLOAT_HUE_TABLES> &>(
          *state);
      if (!history.initialized)
        history.capture(current, persistent_arena);
      if (!Spacetime::flat_surface()) {
        const FrameInterval sampler(history.previous, current);
        HS_PROFILE(fx_shader_draw);
        Scan::Shader::draw_column_times<W, H>(
            canvas,
            [&](int c) {
              HS_PROFILE(fx_prepare_column);
              return RenderPipeline::prepare(
                  sampler.at(Spacetime::SamplingSurface<W>::time(c)));
            },
            [](const Frame &frame, int x) {
              if constexpr (SpecT::FACTORED_COLUMNS &&
                            RenderPipeline::LINEAR_COLUMNS &&
                            &Derived::shade == &ComposedEffect::shade) {
                return [&frame, basis = RenderPipeline::column_basis(
                                    frame, math::TrigLUT<W, H>::cos_theta(x),
                                    math::TrigLUT<W, H>::sin_theta[x])](
                           float sine, float cosine) {
                  return RenderPipeline::shade_column(sine, cosine, frame,
                                                      basis);
                };
              } else {
                return [&frame](const math::Vector &view) {
                  return Derived::shade(view, frame);
                };
              }
            });
      } else {
        const Frame frame = RenderPipeline::prepare(current);
        HS_PROFILE(fx_shader_draw);
        draw_static_columns(canvas, frame);
      }
      history.capture(current, persistent_arena);
    } else {
      const Frame frame = RenderPipeline::prepare(current);
      HS_PROFILE(fx_shader_draw);
      draw_static_columns(canvas, frame);
    }
  }

  /**
   * @brief Shades one pixel; forwards to RenderPipeline::shade.
   * @param view Unit view direction for the pixel.
   * @param frame Per-frame transforms, params and LUTs from the runtime.
   * @return The pixel colour.
   */
  static HS_O3_FN Color4 shade(const math::Vector &view, const Frame &frame) {
    return RenderPipeline::shade(view, frame);
  }

  /** @brief Draws a fixed frame using a factored prefix when available.
   * @param canvas Destination canvas.
   * @param frame Prepared frame to shade. */
  static void draw_static_columns(Canvas &canvas, const Frame &frame) {
    if constexpr (SpecT::FACTORED_COLUMNS && RenderPipeline::LINEAR_COLUMNS &&
                  &Derived::shade == &ComposedEffect::shade) {
      Scan::Shader::draw_columns<W, H>(canvas, [&](int x) {
        return [basis = RenderPipeline::column_basis(
                    frame, math::TrigLUT<W, H>::cos_theta(x),
                    math::TrigLUT<W, H>::sin_theta[x]),
                &frame](float sine, float cosine) {
          return RenderPipeline::shade_column(sine, cosine, frame, basis);
        };
      });
    } else {
      Scan::Shader::draw_cached<W, H, 1>(canvas,
                                         [&frame](const math::Vector &view) {
                                           return Derived::shade(view, frame);
                                         });
    }
  }

  /** @brief Whether a parameter set is admissible, family by family.
   *  @param params Parameter set to check.
   *  @return True when every family is in range. */
  static bool valid_params(const Params &params) {
    return Pullback::valid(params);
  }
#if HS_ENABLE_PARAM_GUI_BRIDGE
  /**
   * @brief Warning text for a slider whose last write was refused.
   * @param name Parameter name.
   * @return Lens::MobiusLensParams::DEGENERATE_WARNING when `name` is the
   *         last refused parameter, else null.
   */
  const char *parameter_warning(const char *name) const override {
    return refused_name != nullptr && std::strcmp(name, refused_name) == 0
               ? Lens::MobiusLensParams::DEGENERATE_WARNING
               : nullptr;
  }
#endif
#if HS_ENABLE_TEST_HOOKS
  /** @brief Prepares and returns this draw's frame state.
      @return The prepare_frame() result. */
  FrameState frame_for_test() { return prepare_frame(); }
#endif

protected:
  /** Descriptors the arena-backed parameter array holds; every slider an effect
      registers has to fit. */
  static constexpr size_t PARAM_CAPACITY =
      ComposedDetail::parameter_count<SpecT>(typename Params::ResourceTypes{});

  using Choreography::anims_paused;
  using Choreography::params;

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /**
   * @brief Refuses a lens write that would make the lens invalid.
   * @details Applies the write to a copy of the targeted lens family and
   * checks it with Pullback::valid(); records the refused name for
   * parameter_warning().
   * @param parameter Descriptor being written.
   * @param value Proposed value.
   * @return False when the write is refused.
   */
  bool parameter_write_admitted(const ParamDef &parameter,
                                float value) override {
    bool admitted = true;
    params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::LENS) {
        const uintptr_t begin = reinterpret_cast<uintptr_t>(&family);
        const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
        if (target >= begin && target - begin < sizeof(family)) {
          auto candidate = family;
          ParamDef proposed = parameter;
          proposed.target =
              reinterpret_cast<unsigned char *>(&candidate) + (target - begin);
          this->write_parameter_unchecked(proposed, value);
          admitted = Pullback::valid(candidate);
        }
      }
    });
    if (!admitted) {
      refused_name = parameter.name;
      return false;
    }
    refused_name = nullptr;
    return true;
  }

  /// Name of the last refused parameter, or null after an admitted write.
  const char *refused_name = nullptr;
#endif

  using Choreography::step_choreography;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;
  using Choreography::use_parameter_storage;

  /**
   * @brief Slider name for a field of one resource instance.
   * @tparam Resource The instance the field belongs to.
   * @param name Unqualified field name.
   * @return `name`, or a copy prefixed with the instance key and a dot,
   *         allocated from `persistent_arena`, when the instance is
   *         qualified.
   */
  template <typename Resource>
  HS_COLD_MEMBER const char *resource_parameter_name(const char *name) {
    constexpr bool QUALIFY =
        ComposedDetail::qualified<Resource, typename Params::ResourceTypes>();
    if constexpr (!QUALIFY)
      return name;
    else {
      constexpr auto KEY = Resource::KEY.view();
      const size_t length = std::strlen(name);
      char *qualified =
          persistent_arena.allocate_n<char>(KEY.size() + length + 2);
      std::memcpy(qualified, KEY.data(), KEY.size());
      qualified[KEY.size()] = '.';
      std::memcpy(qualified + KEY.size() + 1, name, length + 1);
      return qualified;
    }
  }

  /** @brief Parameter folder positions in the pipeline. */
  enum class ParameterStage : uint8_t {
    CAMERA,           ///< Outer camera orientation.
    LENS,             ///< Sphere-space lens.
    SURFACE,          ///< Sphere displacement.
    PROJECTION_FRAME, ///< Projection orientation.
    PROJECTION,       ///< Sphere-to-plane map.
    WARP_1,           ///< Outer planar warp.
    WARP_2,           ///< Inner planar warp.
    FUNCTION,         ///< Scalar source.
    VALUE_TRANSFER,   ///< Scalar transfer.
    COVERAGE,         ///< Coverage function.
    COLORIZE          ///< Palette and hue controls.
  };
#if HS_ENABLE_PARAM_GUI_BRIDGE
  /// Group names, titles and pipeline positions for this spec.
  inline static constexpr auto PARAMETER_GROUPS = [] {
    std::array<ParameterGroup, 11> groups{
        { {"Camera", "Camera", 0},
          {"Lens", "Lens", 1},
          {"Surface Noise", "Surface Noise", 2},
          {"Projection Frame", "Projection Frame", 3},
          {"Projection", "Projection", 4},
          {"Planar Warp 1", "Planar Warp 1", 5},
          {"Planar Warp 2", "Planar Warp 2", 6},
          {"Function", "Function", 7},
          {"Value Transfer", "Value Transfer", 8},
          {"Coverage", "Coverage", 9},
          {"Colorize", "Colorize", 10} }};
    if constexpr (requires { SpecT::parameter_group_title(groups[0].name); })
      for (auto &group : groups)
        group.title = SpecT::parameter_group_title(group.name);
    return groups;
  }();
#endif

  /** @brief Assigns a pipeline group to the most recently registered field.
   *  @param stage Field's pipeline role. */
  void set_parameter_stage(ParameterStage stage) {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    this->set_last_parameter_group(
        &PARAMETER_GROUPS[static_cast<size_t>(stage)]);
#else
    (void)stage;
#endif
  }

  /** @brief Pipeline role of an independent parameter resource.
   *  @tparam Resource Resource supplying the field.
   *  @return Default stage for the resource's fields. */
  template <typename Resource>
  static constexpr ParameterStage resource_parameter_stage() {
    if constexpr (Resource::KIND == ResourceKind::SOURCE)
      return ParameterStage::FUNCTION;
    else if constexpr (Resource::KIND == ResourceKind::SURFACE)
      return ParameterStage::SURFACE;
    else if constexpr (Resource::KIND == ResourceKind::WARP)
      return Resource::KEY.view() == "outer_warp" ? ParameterStage::WARP_1
                                                  : ParameterStage::WARP_2;
    else if constexpr (Resource::KIND == ResourceKind::LENS)
      return ParameterStage::LENS;
    else if constexpr (Resource::KIND == ResourceKind::VALUE)
      return std::is_same_v<typename Resource::Family, IsoValueParams>
                 ? ParameterStage::VALUE_TRANSFER
                 : ParameterStage::COVERAGE;
    else
      return ParameterStage::PROJECTION;
  }

  /** @brief Registers the named descriptors of one independent instance.
   *  @tparam Resource The instance's parameter resource.
   *  @tparam T Parameter family type.
   *  @param family Live family the sliders write. */
  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_fields(T &family) {
    for (const auto &field : T::FIELDS)
      if (field.name != nullptr && field_gate_open(field.gate)) {
        this->register_param(resource_parameter_name<Resource>(field.name),
                             &(family.*field.member), field.description().spec);
        ParameterStage stage = resource_parameter_stage<Resource>();
        if constexpr (Resource::KIND == ResourceKind::PROJECTION) {
          if (field.member == &T::camera_wander)
            stage = ParameterStage::CAMERA;
          else if (field.member == &T::spin_rate || field.member == &T::wander)
            stage = ParameterStage::PROJECTION_FRAME;
        }
        set_parameter_stage(stage);
      }
  }

  /** @brief Adopts a snap target and re-derives the palette mapping weights.
   *  @param target Parameters to adopt. */
  HS_COLD_MEMBER void adopt_params(const Params &target) {
    params = target;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /** @brief Adopts an automatic target while retaining animated lens coefficients.
   *  @param target Parameters to adopt. */
  HS_COLD_MEMBER void finish_blend(const Params &target)
    requires(requires { Derived::ANIMATED_MOBIUS; } && Derived::ANIMATED_MOBIUS)
  {
    params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == ResourceKind::LENS) {
        const auto MOBIUS = family.mobius;
        family = target.template get<Resource::KEY>();
        family.mobius = MOBIUS;
      } else {
        family = target.template get<Resource::KEY>();
      }
    });
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /** @brief Re-derives the palette mapping weights after a slider write. */
  HS_COLD_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.template get<"color">().palette_mapping);
  }

  /** @brief Captures the palette-mapping endpoints of an arming crossfade.
   *  @param target Crossfade's destination parameters. */
  HS_COLD_MEMBER void transition_armed(const Params &target) {
    mapping_from = palette_mapping;
    mapping_to = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /**
   * @brief Writes the interpolated parameters of an in-flight transition.
   * @details Under `Derived::ANIMATED_MOBIUS` the live lens coefficients are
   * carried across the interpolation, so the transition does not overwrite the
   * timeline animation driving them.
   * @param progress Eased transition progress, 0 to 1.
   */
  HS_COLD_MEMBER void blend_params(float progress) {
    const auto before = params;
    params = Pullback::interpolate(transition.from, transition.to, progress);
    if constexpr (requires { Derived::ANIMATED_MOBIUS; })
      if constexpr (Derived::ANIMATED_MOBIUS)
        params.visit([&]<typename Resource>(auto &family) {
          if constexpr (Resource::KIND == ResourceKind::LENS)
            family.mobius = before.template get<Resource::KEY>().mobius;
        });
    palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
        mapping_from, mapping_to, progress);
  }

  /**
   * @brief Starts a repeating circular Mobius warp on the lens parameters.
   * @details Pausable, so the warp freezes with the rest of the parameter
   * animations. Compiles to nothing for an effect whose lens family carries no
   * Mobius coefficients.
   * @param scale Radius of the circular warp.
   * @param duration Frames per revolution.
   */
  template <ResourceKey Key = "lens">
  HS_COLD_MEMBER void start_mobius_animation(float scale, int duration) {
    if constexpr (Params::template HAS<Key>) {
      Animation::MobiusWarpCircular warp(params.template get<Key>().mobius,
                                         scale, duration);
      auto &phase = clocks.template get<Key>();
      warp.bind_phase(phase.angle, phase.magnitude);
      timeline.add_pausable(0, std::move(warp), &anims_paused);
    }
  }

public:
  /** @brief Default composed-effect choreography bindings.
   * @return Type-only hook list owned by this composed base. */
  static consteval auto choreography_hooks() {
    using Hooks = Control::HookSet<
        ComposedEffect,
        Control::HookBinding<Control::ChoreographyHook::VALID_PARAMS,
                             &ComposedEffect::valid_params>,
        Control::HookBinding<Control::ChoreographyHook::ADOPT_PARAMS,
                             &ComposedEffect::adopt_params>,
        Control::HookBinding<Control::ChoreographyHook::BLEND_PARAMS,
                             &ComposedEffect::blend_params>,
        Control::HookBinding<Control::ChoreographyHook::TRANSITION_ARMED,
                             &ComposedEffect::transition_armed>>;
    if constexpr (requires { Derived::ANIMATED_MOBIUS; }) {
      if constexpr (Derived::ANIMATED_MOBIUS)
        return Hooks::template extend<
            ComposedEffect,
            Control::HookBinding<Control::ChoreographyHook::FINISH_BLEND,
                                 &ComposedEffect::finish_blend>>();
      else
        return Hooks{};
    } else
      return Hooks{};
  }

private:
  struct State : ProjectionWalkNoise<ANIMATED_PROJECTION>,
                 OptionalHueRotationLut<Metadata::HUE != HueMode::NONE>,
                 OptionalLinearHueRotationLut<FLOAT_HUE_TABLES>,
                 OptionalHueNoiseLut<Metadata::HUE == HueMode::NOISE>,
                 Spacetime::FrameHistory<FrameState, Metadata::HUE,
                                         TEMPORAL_SAMPLING, FLOAT_HUE_TABLES> {
    ComposedDetail::ResourceStorage<ComposedDetail::ResourceNoise,
                                    typename Params::ResourceTypes>
        resources;
    FastNoiseLite outer_walk_noise;
  };

  // init() takes the palette cycler's generated arena, the external ParamDef
  // array, one State and qualified parameter names from the persistent arena.
  static constexpr size_t FOOTPRINT_BYTES =
      PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State) +
      (TEMPORAL_SAMPLING ? BakedPalette::required_arena_bytes() : 0) +
      alignof(ParamDef) +
      ComposedDetail::parameter_name_bytes<SpecT>(
          typename Params::ResourceTypes{});
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "Pullback::ComposedEffect persistent footprint exceeds the default "
      "partition");

  /** @brief Whether a gated field's slider exists for this effect. */
  static constexpr bool field_gate_open(Pullback::FieldGate gate) {
    return ComposedDetail::field_gate_open<SpecT>(gate);
  }

  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_warp_fields(T &warp) {
    static_assert(T::FIELDS[0].member == &T::speed &&
                      T::FIELDS[0].name == nullptr,
                  "warp speed must be the first unnamed descriptor");
    constexpr const char *SPEED_NAME =
        ComposedDetail::warp_speed_name(Resource::KEY.view()).data();
    this->register_param(resource_parameter_name<Resource>(SPEED_NAME),
                         &warp.speed, T::FIELDS[0].description().spec);
    set_parameter_stage(resource_parameter_stage<Resource>());
    register_fields<Resource>(warp);
  }

  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_lens_fields(T &lens) {
    constexpr math::Complex math::MobiusParams::*COEFFICIENTS[] = {
        &math::MobiusParams::a, &math::MobiusParams::b, &math::MobiusParams::c,
        &math::MobiusParams::d};
    constexpr float math::Complex::*CHANNELS[] = {&math::Complex::re,
                                                  &math::Complex::im};
    constexpr float LIMIT = MobiusLensParams::COEFFICIENT_LIMIT;
    for (size_t i = 0; i < std::size(COEFFICIENTS); ++i)
      for (size_t j = 0; j < std::size(CHANNELS); ++j) {
        register_animated_param(
            resource_parameter_name<Resource>(
                ComposedDetail::MOBIUS_PARAM_NAMES[i * 2 + j]),
            &((lens.mobius.*COEFFICIENTS[i]).*CHANNELS[j]), -LIMIT, LIMIT);
        set_parameter_stage(ParameterStage::LENS);
      }
  }

  template <float Color::ColorControls::*Member>
  static consteval const Field<ColorParams> &color_descriptor() {
    constexpr size_t INDEX = [] {
      for (size_t i = 0; i < ColorParams::FIELDS.size(); ++i)
        if (ColorParams::FIELDS[i].member == Member)
          return i;
      return ColorParams::FIELDS.size();
    }();
    static_assert(INDEX < ColorParams::FIELDS.size(),
                  "ColorParams member has no field descriptor");
    return ColorParams::FIELDS[INDEX];
  }

  /** @brief Registers a colour slider when its topology gate is live. */
  template <float Color::ColorControls::*Member>
  HS_COLD_MEMBER void register_color_field(const char *name) {
    constexpr const auto &FIELD = color_descriptor<Member>();
    if constexpr (ComposedDetail::color_gate_open<SpecT>(FIELD.topology_gate)) {
      register_animated_param(name, &(params.template get<"color">().*Member),
                              FIELD.min, FIELD.max);
      set_parameter_stage(ParameterStage::COLORIZE);
    }
  }

  template <ResourceKind Kind> HS_COLD_MEMBER void register_resource_kind() {
    params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == Kind) {
        if constexpr (Kind == ResourceKind::WARP)
          register_warp_fields<Resource>(family);
        else if constexpr (Kind == ResourceKind::LENS)
          register_lens_fields<Resource>(family);
        else
          register_fields<Resource>(family);
      }
    });
  }

  HS_COLD_MEMBER void register_parameters() {
    register_resource_kind<ResourceKind::SOURCE>();
    register_resource_kind<ResourceKind::PROJECTION>();
    register_resource_kind<ResourceKind::SURFACE>();
    register_resource_kind<ResourceKind::WARP>();
    register_resource_kind<ResourceKind::VALUE>();
    register_resource_kind<ResourceKind::LENS>();
    register_color_field<&ColorParams::palette_chroma>("Palette Chroma");
    register_animated_param(
        "Palette Mapping", &params.template get<"color">().palette_mapping,
        PALETTE_MAPPING_OPTIONS, PALETTE_MAPPING_EXPORT_OPTIONS,
        std::size(PALETTE_MAPPING_OPTIONS));
    set_parameter_stage(ParameterStage::COLORIZE);
    register_color_field<&ColorParams::mapping_frequency>("Mapping Frequency");
    register_color_field<&ColorParams::mapping_phase>("Mapping Phase");
    register_color_field<&ColorParams::phase_oscillation_depth>(
        "Phase Oscillation Depth");
    register_color_field<&ColorParams::phase_oscillation_speed>(
        "Phase Oscillation Speed");
    register_color_field<&ColorParams::brightness_bottom>("Brightness Bottom");
    register_color_field<&ColorParams::brightness_top>("Brightness Top");
    register_color_field<&ColorParams::opacity_low>("Opacity at Value 0");
    register_color_field<&ColorParams::opacity_high>("Opacity at Value 1");
    register_color_field<&ColorParams::hue_shift_amount>("Hue Shift Amount");
    register_color_field<&ColorParams::hue_noise_scale>("Hue Noise Scale");
    register_color_field<&ColorParams::hue_noise_speed>("Hue Noise Speed");
  }

  /**
   * @brief Steps every phase clock the effect's parameter families define.
   * @details Each clock is compiled in only when its field exists. Affine
   * rotations accumulate independently for each affine warp.
   */
  HS_COLD_MEMBER void advance_runtime() {
    params.visit([&]<typename Resource>(const auto &family) {
      auto &clock = clocks.template get<Resource::KEY>();
      if constexpr (Resource::KIND == ResourceKind::SOURCE) {
        Source::advance_clocks(family, clock.primary, clock.secondary,
                               clock.angle);
        if constexpr (requires { family.noise_time_rate; })
          clock.noise_time =
              math::wrap_t(clock.noise_time + family.noise_time_rate);
      } else if constexpr (Resource::KIND == ResourceKind::SURFACE) {
        if constexpr (std::is_same_v<typename Resource::Family,
                                     PeriodicRippleParams>)
          Surface::advance_ripple_phase(clock.phase, family);
        else
          clock.phase = math::wrap_t(clock.phase + family.speed);
      } else if constexpr (Resource::KIND == ResourceKind::WARP) {
        if constexpr (std::is_same_v<typename Resource::Family, AffineParams>)
          Warp::advance_affine_rotation(clock.rotation, family);
        clock.phase = math::wrap_t(clock.phase + family.speed);
      }
    });
    if constexpr (ANIMATED_PROJECTION)
      this->projection_spin = fmodf(
          this->projection_spin + params.template get<"projection">().spin_rate,
          math::TWO_PI_F);
    if constexpr (requires { Derived::CAMERA_SPIN_RATE; })
      camera_spin =
          fmodf(camera_spin + Derived::CAMERA_SPIN_RATE, math::TWO_PI_F);
    if constexpr (Metadata::HUE == HueMode::NOISE)
      hue_noise_phase = math::wrap_t(
          hue_noise_phase + params.template get<"color">().hue_noise_speed);
    palette_oscillation_phase =
        math::wrap_t(palette_oscillation_phase +
                     params.template get<"color">().phase_oscillation_speed);
  }

  HS_COLD_MEMBER void update_spatial_frames() {
    if constexpr (ANIMATED_PROJECTION) {
      const math::Quaternion projection = this->projection_walk.get();
      const math::Quaternion projection_delta =
          projection * this->projection_walk_previous.conjugate();
      this->projection_walk_previous = projection;
      this->projection_wander =
          (math::scaled_rotation_delta(
               projection_delta.normalized(),
               params.template get<"projection">().wander) *
           this->projection_wander)
              .normalized();
      this->projection_conjugate =
          (math::make_rotation(math::Y_AXIS, this->projection_spin) *
           this->base_orientation * this->projection_wander)
              .conjugate();
    }
    const math::Quaternion outer = outer_walk.get();
    const math::Quaternion outer_delta =
        outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander = (math::scaled_rotation_delta(
                        outer_delta.normalized(),
                        params.template get<"projection">().camera_wander) *
                    outer_wander)
                       .normalized();
    if constexpr (requires { Derived::CAMERA_SPIN_RATE; })
      outer_conjugate =
          (math::make_rotation(math::Y_AXIS, camera_spin) * outer_wander)
              .conjugate();
    else
      outer_conjugate = outer_wander.conjugate();
  }

  /** @brief The hue-rotation LUT base, or null when no hue mode reads one. */
  const Pixel *hue_rotation_lut_data() const {
    if constexpr (Metadata::HUE != HueMode::NONE)
      return state->hue_rotation_lut.data();
    else
      return nullptr;
  }

  /** @brief Optional preconverted hue table.
   * @return Float RGB entries, or null when disabled. */
  const math::Vector *hue_rotation_linear_data() const {
    if constexpr (FLOAT_HUE_TABLES)
      return static_cast<const OptionalLinearHueRotationLut<true> &>(*state)
          .hue_rotation_linear.data();
    else
      return nullptr;
  }

  /** @brief Optional preconverted noise table.
   * @return Quantized noise as floats, or null when disabled. */
  const float *hue_noise_linear_data() const {
    if constexpr (FLOAT_HUE_TABLES)
      return static_cast<const OptionalLinearHueRotationLut<true> &>(*state)
          .hue_noise_linear.data();
    else
      return nullptr;
  }

  /** @brief The hue-noise LUT base, or null outside HueMode::NOISE. */
  const int8_t *hue_noise_lut_data() const {
    if constexpr (Metadata::HUE == HueMode::NOISE)
      return state->hue_noise_lut.data();
    else
      return nullptr;
  }

  /**
   * @brief Bakes the frame's LUTs and snapshots everything the scan reads.
   * @details Column-time sampling retains valid hue tables even while hue rotation
   * is inactive. Other effects gate their builds on hue_rotation_active().
   * Hue noise is rebuilt when its scale or phase moves; hue rotation is rebuilt
   * when the palette cycler rebakes its display LUT.
   * @return The frame state for this draw, valid until the next draw_frame().
   */
  HS_COLD_MEMBER FrameState prepare_frame() {
    HS_PROFILE(fx_prepare_frame);
    if constexpr (Metadata::HUE == HueMode::NOISE) {
      if (TEMPORAL_SAMPLING ||
          hue_rotation_active<Metadata::HUE>(params.template get<"color">())) {
        [[maybe_unused]] const bool REBUILT = state->hue_noise_bake.refresh(
            state->hue_noise_lut, state->color_noise,
            params.template get<"color">().hue_noise_scale, hue_noise_phase);
        if constexpr (FLOAT_HUE_TABLES) {
          if (REBUILT) {
            auto &linear =
                static_cast<OptionalLinearHueRotationLut<true> &>(*state)
                    .hue_noise_linear;
            for (size_t index = 0; index < linear.size(); ++index)
              linear[index] = state->hue_noise_lut[index];
          }
        }
      }
    }
    if constexpr (Metadata::HUE != HueMode::NONE)
      if ((TEMPORAL_SAMPLING || hue_rotation_active<Metadata::HUE>(
                                    params.template get<"color">())) &&
          state->hue_rotation_lut_bake != palette_cycler.bake_generation()) {
        Pullback::Color::prepare_hue_rotation_lut(
            std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
                state->hue_rotation_lut),
            palette_cycler.palette());
        if constexpr (FLOAT_HUE_TABLES) {
          auto &linear =
              static_cast<OptionalLinearHueRotationLut<true> &>(*state)
                  .hue_rotation_linear;
          for (size_t index = 0; index < linear.size(); ++index) {
            const Pixel &RGB = state->hue_rotation_lut[index];
            linear[index] = math::Vector(RGB.r, RGB.g, RGB.b);
          }
        }
        state->hue_rotation_lut_bake = palette_cycler.bake_generation();
      }
    Pullback::FrameState<ParamsT> frame{
        .projection_conjugate = this->frame_conjugate(),
        .outer_conjugate = outer_conjugate,
        .palette = &palette_cycler.palette(),
        .hue_rotation_lut = hue_rotation_lut_data(),
        .hue_noise_lut = hue_noise_lut_data(),
        .params = params,
        .palette_mapping = palette_mapping,
        .resources = {},
        .palette_oscillation_phase = palette_oscillation_phase,
        .hue_rotation_linear = hue_rotation_linear_data(),
        .palette_alpha =
            palette_cycler.palette().constant_alpha().value_or(-1.0f),
        .hue_noise_linear = hue_noise_linear_data()};
    params.visit([&]<typename Resource>(const auto &) {
      auto &resource = frame.resources.template get<Resource::KEY>();
      static_cast<ComposedDetail::ResourceClock<Resource> &>(resource) =
          clocks.template get<Resource::KEY>();
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>)
        resource.noise = &state->resources.template get<Resource::KEY>().noise;
      if constexpr (Resource::KIND == ResourceKind::LENS)
        if (!Pullback::valid(frame.params.template get<Resource::KEY>()))
          frame.params.template get<Resource::KEY>().mobius = {};
      if constexpr (std::is_same_v<typename Resource::Family, MobiusLensParams>)
        resource.generated = lens_generated(
            frame.params.template get<Resource::KEY>(), resource);
    });
    if constexpr (TEMPORAL_SAMPLING)
      return FrameState{frame,
                        {frame.palette, frame.hue_rotation_lut,
                         frame.hue_noise_lut, 1.0f, frame.hue_rotation_linear,
                         frame.palette_alpha, frame.hue_noise_linear}};
    else
      return frame;
  }

  HS_COLD_MEMBER void update_palette_chroma() {
    if (palette_chroma == params.template get<"color">().palette_chroma)
      return;
    palette_chroma = params.template get<"color">().palette_chroma;
    palette_cycler.set_generated_chroma(palette_chroma);
  }

  /**
   * @brief Generator the palette cycler calls for each palette in the cycle.
   * @details Advances the hue on every palette after the first, so the cycle
   * opens on `palette_hue`, which starts at 0 and moves only here.
   * @param context The ComposedEffect instance, as registered with the cycler.
   * @param sequence Zero-based index of the palette being generated.
   * @param out Receives the recipe to bake.
   */
  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    ComposedEffect &effect = *static_cast<ComposedEffect *>(context);
    GeneratedPaletteBank::next_palette(
        effect.palette_hue, sequence, SpecT::HARMONY,
        effect.params.template get<"color">().palette_chroma, out);
  }

  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "Pullback::Color::PaletteMapping::CUP",
      "Pullback::Color::PaletteMapping::BELL",
      "Pullback::Color::PaletteMapping::LINEAR",
      "Pullback::Color::PaletteMapping::REVERSE"};
  static_assert(
      std::size(PALETTE_MAPPING_OPTIONS) ==
          static_cast<size_t>(Pullback::Color::PaletteMapping::REVERSE) + 1 &&
      std::size(PALETTE_MAPPING_EXPORT_OPTIONS) ==
          std::size(PALETTE_MAPPING_OPTIONS));

  State *state = nullptr;
  Pullback::Color::PaletteMappingWeights palette_mapping =
      Pullback::Color::PaletteMappingWeights::single(
          params.template get<"color">().palette_mapping);
  Pullback::Color::PaletteMappingWeights mapping_from;
  Pullback::Color::PaletteMappingWeights mapping_to;
  math::Orientation<> outer_walk;
  math::Quaternion outer_walk_previous;
  math::Quaternion outer_wander;
  math::Quaternion outer_conjugate;
  ComposedDetail::ResourceStorage<ComposedDetail::ResourceClock,
                                  typename Params::ResourceTypes>
      clocks;
  float camera_spin = 0.0f;
  float hue_noise_phase = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};
