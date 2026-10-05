/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/pullback/composed_effect.h.

/**
 * @brief A complete composed effect: the shared lifecycle plus a pipeline
 *        declared by a ranked stage Spec.
 * @details The effect states its Spec and identity constants; parameter and
 * runtime storage derive from its pipeline providers. The shared
 * lifecycle — parameter registration, preset choreography, palette cycling,
 * camera walks and noise clocks — are assembled here. Required `Derived`
 * members are EFFECT_ID, PRESET_IDS, PARAMETER_SCHEMA_VERSION and
 * PRESET_DWELL_FRAMES. Presets and their departures resolve
 * through `preset(index)`, then `PRESETS`; only single-preset effects may fall
 * back to startup params.
 * `initial_params` is optional. Other optional members are `ANIMATED_MOBIUS`,
 * `CAMERA_SPIN_RATE` and an `after_composed_init()` hook; `WARP_NOISE_SEED` /
 * `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` are inherited members an effect
 * shadows to decorrelate its warp, source, or surface noise fields. A shade() shadow that forwards to
 * RenderPipeline::shade changes only the entry trampoline's placement; the
 * pipeline body remains in hot flash. Different body emission requires calling
 * RenderPipeline::evaluate(view, frame.ctx, frame.prepared) from the shadow.
 * Surface-noise effects conventionally wrap the sphere run (displacement, lens
 * and projection) in Stage::Placed<CodeEmission::OUT_OF_LINE_FLASH, ...>.
 *
 * `EFFECT_ID` is the registry identity; `PRESET_IDS` lists immutable preset
 * identities indexed by preset number; `PARAMETER_SCHEMA_VERSION` changes
 * with the Params layout to reject stale snapshots; `PRESET_DWELL_FRAMES`
 * gives the frames held before the next transition.
 * `DESCRIPTOR_DIGEST` and `PRESET_BANK_DIGEST` pin the pattern document's
 * canonical descriptor (excluding parameter units) and preset bank. The
 * product-group generator parity tests check them; the runtime never reads
 * them. The browser computes its own matching digests from pattern documents.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this base.
 * @tparam SpecT The effect's `Spec`.
 */
template <int W, int H, typename Derived, typename SpecT>
class ComposedEffect : public ChoreographedEffect<Derived, ParamsFor<SpecT>>,
                       private ProjectionWalkState<SpecT::ANIMATED_PROJECTION> {
  using ParamsT = ParamsFor<SpecT>;
  using Choreography = ChoreographedEffect<Derived, ParamsT>;
  friend Choreography;
  static constexpr PaletteHarmony Harmony = SpecT::HARMONY;
  static constexpr HueMode HueV = SpecT::HUE;
  static constexpr Color::BrightnessEnvelope BrightnessV = SpecT::BRIGHTNESS;
  static constexpr bool AnimatedProjection = SpecT::ANIMATED_PROJECTION;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  using Params = ParamsT;
  using Spec = SpecT;
  static_assert(
      Spec::TRANSFER != TransferKind::ISO_CONTOUR || requires(Params p) {
        p.template get<"value">().iso_level;
        p.template get<"value">().iso_width;
      }, "iso-contour transfer requires iso-level and iso-width");
  static_assert(
      Spec::COVERAGE != ProjectionCoverageMode::EDGE_FADE ||
          requires(Params p) { p.template get<"value">().edge_width; },
      "edge-fade coverage requires edge-width");
  static_assert(
      Spec::FIELD_COVERAGE != FieldCoverageKind::VALUE_CUTOUT ||
          requires(Params p) {
            p.template get<"value">().cutout_threshold;
            p.template get<"value">().cutout_softness;
          },
      "value-cutout coverage requires threshold and softness");
  using FrameState = Pullback::FrameState<ParamsT>;
  using Binding = Pullback::Binding<FrameState>;
  static constexpr bool ANIMATED_PROJECTION = AnimatedProjection;
  template <ResourceKind Kind> static consteval bool has_noise() {
    bool result = false;
    Params{}.visit([&]<typename Resource>(const auto &) {
      if constexpr (Resource::KIND == Kind &&
                    ComposedDetail::RESOURCE_NOISE<Resource>)
        result = true;
    });
    return result;
  }
  static constexpr bool HAS_SURFACE_NOISE = has_noise<ResourceKind::SURFACE>();
  static constexpr bool HAS_WARP_NOISE = has_noise<ResourceKind::WARP>();
  static constexpr bool HAS_SOURCE_NOISE = has_noise<ResourceKind::SOURCE>();
  using RenderPipeline = typename SpecT::template Pipeline<Binding>;
  using Metadata =
      ComposedDetail::PipelineMetadata<Spec, Binding, RenderPipeline>;
  static_assert(Metadata::PATH_TRACKED == (Spec::HUE == HueMode::PATH_LENGTH),
                "path length hue metadata must match pipeline tracking");
  static_assert(Metadata::COVERAGE_MATCHES,
                "projection coverage metadata must match the pipeline");
  static_assert(Metadata::LENS_MATCHES,
                "lens policy metadata must match the pipeline");
  static_assert(Metadata::PROJECTION_MATCHES,
                "projection metadata must match the pipeline");
  static_assert(Metadata::SURFACE_PLACEMENT_MATCHES,
                "surface placement metadata must match the pipeline");
  static_assert(
      RenderPipeline::template any_stage<ComposedDetail::IsTransferStage> ==
          (Spec::TRANSFER != TransferKind::NONE),
      "transfer metadata must match the pipeline");
  static_assert(
      RenderPipeline::template any_stage<ComposedDetail::IsCoverageStage> ==
          (Spec::FIELD_COVERAGE != FieldCoverageKind::NONE),
      "field coverage metadata must match the pipeline");
  using Frame = typename RenderPipeline::Frame;

  /** @brief Per-field noise seeds; an effect shadows one to decorrelate its
      field from the shared spatial phase. */
  static constexpr int32_t WARP_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SOURCE_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SURFACE_NOISE_SEED = EFFECT_NOISE_SEED;

  /** @brief Constructs the effect at W x H with the POV column strobe on. */
  HS_COLD_MEMBER ComposedEffect() : Choreography(W, H, {.strobe = true}) {}

  /**
   * @brief Claims the runtime's persistent storage and registers the effect.
   * @details Every allocation this makes comes from `persistent_arena`, so the
   * effect heap-allocates nothing after init. `Derived::after_composed_init()`
   * runs last, once the parameters, palette cycler and camera walks are all
   * live, which is what lets it adjust the dwell or start a timeline
   * animation.
   */
  HS_COLD_MEMBER void init() override {
    this->begin_choreography();
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena,
                          persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    if constexpr (HueV == HueMode::NOISE)
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
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  math::ease_in_out_sin);
    if constexpr (AnimatedProjection)
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
    if constexpr (requires(Derived &effect) { effect.after_composed_init(); })
      static_cast<Derived &>(*this).after_composed_init();
  }

  /**
   * @brief Advances the frame clocks, resolves the frame state and shades.
   * @details The order is load-bearing: the runtime's clocks, camera walks and
   * palette all step before prepare_frame() snapshots them, so the whole scan
   * shades from one consistent frame. The preset-transition lerp steps with
   * the timeline, so its writes land before the snapshot too.
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
    const typename Derived::RenderPipeline::Frame frame =
        Derived::RenderPipeline::prepare(prepare_frame());
    {
      HS_PROFILE(fx_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas,
                                         [&frame](const math::Vector &view) {
                                           return Derived::shade(view, frame);
                                         });
    }
  }

  /**
   * @brief Shades one pixel; forwards to RenderPipeline::shade.
   * @param view Unit view direction for the pixel.
   * @param frame Per-frame transforms, params and LUTs from the runtime.
   */
  static HS_O3_FN Color4 shade(const math::Vector &view, const Frame &frame) {
    return RenderPipeline::shade(view, frame);
  }

  /** @brief Whether a parameter set is admissible, family by family. */
  static bool valid_params(const Params &params) {
    return Pullback::valid(params);
  }
#if HS_ENABLE_PARAM_GUI_BRIDGE
  const char *parameter_warning(const char *name) const override {
    return refused_name != nullptr && std::strcmp(name, refused_name) == 0
               ? Lens::MobiusLensParams::DEGENERATE_WARNING
               : nullptr;
  }
#endif
#if HS_ENABLE_TEST_HOOKS
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

  const char *refused_name = nullptr;
#endif

  using Choreography::step_choreography;
  using Choreography::register_animated_param;
  using Choreography::timeline;
  using Choreography::transition;
  using Choreography::use_parameter_storage;

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

  /** @brief Registers the named descriptors of one independent instance. */
  template <typename Resource, typename T>
  HS_COLD_MEMBER void register_fields(T &family) {
    for (const auto &field : T::FIELDS)
      if (field.name != nullptr && field_gate_open(field.gate))
        this->register_param(resource_parameter_name<Resource>(field.name),
                             &(family.*field.member), field.description().spec);
  }

  /** @brief Adopts a snap target and re-derives the palette mapping weights. */
  HS_COLD_MEMBER void adopt_params(const Params &target) {
    params = target;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        target.template get<"color">().palette_mapping);
  }

  /** @brief Adopts an automatic target while retaining animated lens coefficients. */
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

  HS_COLD_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.template get<"color">().palette_mapping);
  }

  /** @brief Captures the palette-mapping endpoints of an arming crossfade. */
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
    if constexpr (Params::template HAS<Key>)
      timeline.add_pausable(
          0,
          Animation::MobiusWarpCircular(params.template get<Key>().mobius,
                                        scale, duration, true),
          &anims_paused);
  }

private:
  struct State : ProjectionWalkNoise<AnimatedProjection>,
                 OptionalHueRotationLut<HueV != HueMode::NONE>,
                 OptionalHueNoiseLut<HueV == HueMode::NOISE> {
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
      for (size_t j = 0; j < std::size(CHANNELS); ++j)
        register_animated_param(
            resource_parameter_name<Resource>(
                ComposedDetail::MOBIUS_PARAM_NAMES[i * 2 + j]),
            &((lens.mobius.*COEFFICIENTS[i]).*CHANNELS[j]), -LIMIT, LIMIT);
  }

  template <float Color::ColorControls::*Member>
  static consteval const Field<ColorParams> &color_descriptor() {
    constexpr size_t index = [] {
      for (size_t i = 0; i < ColorParams::FIELDS.size(); ++i)
        if (ColorParams::FIELDS[i].member == Member)
          return i;
      return ColorParams::FIELDS.size();
    }();
    static_assert(index < ColorParams::FIELDS.size(),
                  "ColorParams member has no field descriptor");
    return ColorParams::FIELDS[index];
  }

  template <float Color::ColorControls::*Member>
  HS_COLD_MEMBER void register_color_field(const char *name) {
    constexpr const auto &field = color_descriptor<Member>();
    register_animated_param(name, &(params.template get<"color">().*Member),
                            field.min, field.max);
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
    register_color_field<&ColorParams::mapping_frequency>("Mapping Frequency");
    register_color_field<&ColorParams::mapping_phase>("Mapping Phase");
    register_color_field<&ColorParams::phase_oscillation_depth>(
        "Phase Oscillation Depth");
    register_color_field<&ColorParams::phase_oscillation_speed>(
        "Phase Oscillation Speed");
    if constexpr (BrightnessV != Pullback::Color::BrightnessEnvelope::NONE) {
      register_color_field<&ColorParams::brightness_bottom>(
          "Brightness Bottom");
      register_color_field<&ColorParams::brightness_top>("Brightness Top");
    }
    register_color_field<&ColorParams::opacity_low>("Opacity at Value 0");
    register_color_field<&ColorParams::opacity_high>("Opacity at Value 1");
    if constexpr (HueV != HueMode::NONE)
      register_color_field<&ColorParams::hue_shift_amount>("Hue Shift Amount");
    if constexpr (HueV == HueMode::NOISE) {
      register_color_field<&ColorParams::hue_noise_scale>("Hue Noise Scale");
      register_color_field<&ColorParams::hue_noise_speed>("Hue Noise Speed");
    }
  }

  /**
   * @brief Steps every phase clock the effect's parameter families define.
   * @details Each clock is compiled in only when its field exists, so an effect
   * pays for the clocks of its declared instances. Affine rotations are
   * accumulated independently for each affine warp.
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
    if constexpr (AnimatedProjection)
      this->projection_spin = fmodf(
          this->projection_spin + params.template get<"projection">().spin_rate,
          math::TWO_PI_F);
    if constexpr (requires { Derived::CAMERA_SPIN_RATE; })
      camera_spin =
          fmodf(camera_spin + Derived::CAMERA_SPIN_RATE, math::TWO_PI_F);
    if constexpr (HueV == HueMode::NOISE)
      hue_noise_phase = math::wrap_t(
          hue_noise_phase + params.template get<"color">().hue_noise_speed);
    palette_oscillation_phase =
        math::wrap_t(palette_oscillation_phase +
                     params.template get<"color">().phase_oscillation_speed);
  }

  // Rotation samples are eased within each walk step; chain walks apply the
  // recurrence directly and use different seeds. Nonzero wander differs.
  HS_COLD_MEMBER void update_spatial_frames() {
    // prepare_frame() reads projection_conjugate only for an animated
    // projection.
    if constexpr (AnimatedProjection) {
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
    if constexpr (HueV != HueMode::NONE)
      return state->hue_rotation_lut.data();
    else
      return nullptr;
  }

  /** @brief The hue-noise LUT base, or null outside HueMode::NOISE. */
  const int8_t *hue_noise_lut_data() const {
    if constexpr (HueV == HueMode::NOISE)
      return state->hue_noise_lut.data();
    else
      return nullptr;
  }

  /**
   * @brief Bakes the frame's LUTs and snapshots everything the scan reads.
   * @details Both LUT builds are gated on hue_rotation_active(), the flag the
   * returned frame hands the color stage; under it the hue-noise LUT is rebuilt
   * only when its scale or phase moved, and the hue-rotation LUT only when the
   * palette cycler rebaked its display LUT.
   * @return The frame state for this draw, valid until the next draw_frame().
   */
  HS_COLD_MEMBER FrameState prepare_frame() {
    HS_PROFILE(fx_prepare_frame);
    if constexpr (HueV == HueMode::NOISE) {
      if (hue_rotation_active<HueV>(params.template get<"color">())) {
        state->hue_noise_bake.refresh(
            state->hue_noise_lut, state->color_noise,
            params.template get<"color">().hue_noise_scale, hue_noise_phase);
      }
    }
    if constexpr (HueV != HueMode::NONE)
      if (hue_rotation_active<HueV>(params.template get<"color">()) &&
          state->hue_rotation_lut_bake != palette_cycler.bake_generation()) {
        Pullback::Color::prepare_hue_rotation_lut(
            std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
                state->hue_rotation_lut),
            palette_cycler.palette());
        state->hue_rotation_lut_bake = palette_cycler.bake_generation();
      }
    FrameState frame{.projection_conjugate = this->frame_conjugate(),
                     .outer_conjugate = outer_conjugate,
                     .palette = &palette_cycler.palette(),
                     .hue_rotation_lut = hue_rotation_lut_data(),
                     .hue_noise_lut = hue_noise_lut_data(),
                     .params = params,
                     .palette_mapping = palette_mapping,
                     .resources = {},
                     .palette_oscillation_phase = palette_oscillation_phase};
    params.visit([&]<typename Resource>(const auto &) {
      auto &resource = frame.resources.template get<Resource::KEY>();
      static_cast<ComposedDetail::ResourceClock<Resource> &>(resource) =
          clocks.template get<Resource::KEY>();
      if constexpr (ComposedDetail::RESOURCE_NOISE<Resource>)
        resource.noise = &state->resources.template get<Resource::KEY>().noise;
      if constexpr (Resource::KIND == ResourceKind::LENS)
        if (!Pullback::valid(frame.params.template get<Resource::KEY>()))
          frame.params.template get<Resource::KEY>().mobius = {};
    });
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
    if (sequence > 0)
      effect.palette_hue += 159;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, Harmony, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(effect.palette_hue),
        effect.params.template get<"color">().palette_chroma)};
  }

  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "Pullback::Color::PaletteMapping::CUP",
      "Pullback::Color::PaletteMapping::BELL",
      "Pullback::Color::PaletteMapping::LINEAR",
      "Pullback::Color::PaletteMapping::REVERSE"};

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
