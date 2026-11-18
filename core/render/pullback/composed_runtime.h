/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <span>
#include <type_traits>
#include "core/animation/orientation.h"
#include "core/color/authored_palette.h"
#include "core/color/palette_cycler.h"
#include "core/control/choreography.h"
#include "core/control/registry.h"
#include "core/math/mobius.h"
#include "core/math/noise_field.h"
#include "core/memory.h"
#include "core/render/pullback/composed_admission.h"
#include "core/render/pullback/composed_clocks.h"
#include "core/render/pullback/composed_descriptors.h"
#include "core/render/pullback/composed_interval.h"
#include "core/render/pullback/composed_palette.h"
#include "core/render/pullback/composed_registration.h"
#include "core/render/pullback/composed_traits.h"
#include "core/render/pullback/runtime_seeds.h"
#include "core/render/scan.h"
#include "core/render/spacetime.h"
#include "core/render/spacetime/field.h"

#if HS_ENABLE_TEST_HOOKS
namespace hs_test {
struct ComposedFrameWhiteBox;
}
#endif

/** @file composed_runtime.h
 * @brief `ComposedEffect`, the lifecycle base of composed pullback effects.
 */

namespace Pullback {

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
 * `CAMERA_SPIN_RATE`, and `WARP_NOISE_SEED` /
 * `SOURCE_NOISE_SEED` / `SURFACE_NOISE_SEED` shadows. A shade() shadow supplies
 * the per-pixel body. An optional `prepare_column_shader(frame, x)` factory
 * returns a column shader accepting latitude sine and cosine.
 * `DESCRIPTOR_DIGEST` and `PRESET_BANK_DIGEST` pin the pattern document; the
 * runtime never reads them.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Derived The effect class deriving from this base.
 * @tparam SpecT The effect's `Spec`.
 * @tparam Ownership Whether Derived declares hooks or inherits composed
 * defaults.
 */
template <int W, int H, typename Derived, typename SpecT,
          Control::HookOwnership Ownership = Control::HookOwnership::AUTHORED>
class ComposedEffect
    : public ChoreographedEffect<
          Derived, ParamsFor<SpecT>,
          std::conditional_t<Ownership == Control::HookOwnership::INHERITED,
                             ComposedEffect<W, H, Derived, SpecT, Ownership>,
                             Derived>>,
      public ComposedTraits<W, SpecT>,
      public ComposedAdmission<ComposedEffect<W, H, Derived, SpecT, Ownership>>,
      private ComposedRegistration<
          ComposedEffect<W, H, Derived, SpecT, Ownership>,
          ComposedTraits<W, SpecT>>,
      private ProjectionWalkState<
          ComposedTraits<W, SpecT>::ANIMATED_PROJECTION>,
      private OuterCameraState<MetadataFor<SpecT>::OUTER_CAMERA> {
  using Traits = ComposedTraits<W, SpecT>;
  using ParamsT = ParamsFor<SpecT>;
  using Choreography = ChoreographedEffect<
      Derived, ParamsT,
      std::conditional_t<Ownership == Control::HookOwnership::INHERITED,
                         ComposedEffect, Derived>>;
  using Registration = ComposedRegistration<ComposedEffect, Traits>;
  using IntervalCache = ComposedIntervalCache<Traits>;
  using Clocks = ComposedClocks<Traits>;
  friend Choreography;
  friend Registration;
  friend ComposedAdmission<ComposedEffect>;
#if HS_ENABLE_TEST_HOOKS
  friend struct hs_test::ComposedFrameWhiteBox;
#endif

public:
  using typename Traits::ColorFamily;
  using typename Traits::Frame;
  using typename Traits::FrameState;
  using typename Traits::Metadata;
  using typename Traits::Params;
  using typename Traits::RenderPipeline;
  /// Whether the projection orientation follows a random walk and spin.
  static constexpr bool ANIMATED_PROJECTION = Traits::ANIMATED_PROJECTION;
  /// Whether the color controls carry full authored sources.
  static constexpr bool AUTHORED_PALETTE = Traits::AUTHORED_PALETTE;
  /// Whether runtime identities follow authored nodes.
  static constexpr bool AUTHORED_TOPOLOGY = Traits::AUTHORED_TOPOLOGY;
  /// Whether the effect retains palette hue tables.
  static constexpr bool HUE_ROTATION_LUT = Traits::HUE_ROTATION_LUT;
  /// Whether the palette retains interval OKLab bases.
  static constexpr bool PALETTE_LAB_CACHE = Traits::PALETTE_LAB_CACHE;
  /// Whether the palette retains final hue colors over its interval bases.
  static constexpr bool PALETTE_RGB_CACHE = Traits::PALETTE_RGB_CACHE;
  /// Whether the effect samples time across each two-arm half-turn.
  static constexpr bool TEMPORAL_SAMPLING = Traits::TEMPORAL_SAMPLING;

  /** @brief Per-field noise seeds; an effect shadows one to decorrelate its
      field from the shared spatial phase. */
  static constexpr int32_t WARP_NOISE_SEED = EFFECT_NOISE_SEED;
  static constexpr int32_t SOURCE_NOISE_SEED = EFFECT_NOISE_SEED;  ///< Source.
  static constexpr int32_t SURFACE_NOISE_SEED = EFFECT_NOISE_SEED; ///< Surface.

  /** @brief Resolves a semantic resource key to its authored node label.
   * @param resource Parameter resource key.
   * @return The declared authored label, or the resource key. */
  static constexpr std::string_view
  authored_instance_label(std::string_view resource) {
    if constexpr (requires { Derived::authored_instance(resource); })
      return Derived::authored_instance(resource);
    else if constexpr (requires { SpecT::authored_instance(resource); })
      return SpecT::authored_instance(resource);
    else
      return resource;
  }

  /** @brief Constructs the effect at W x H with its lens column-strobe mode. */
  HS_SIZE_FLASH_MEMBER ComposedEffect()
      : Choreography(W, H, {.strobe = Metadata::STROBE_COLUMNS}) {}

  /**
   * @brief Claims the runtime's persistent storage and registers the effect.
   * @details Every allocation this makes comes from `persistent_arena`, so the
   * effect heap-allocates nothing after init. The `AFTER_COMPOSED_INIT` hook
   * runs last, once the parameters, palette cycler and camera walks are live.
   */
  HS_FLASH_MEMBER void init() override {
    static_assert(
        modulation_presets_are_supported(),
        "composed modulation presets must retain source and mechanic");
    this->begin_choreography();
    state = persistent_arena.make<State>();
    if constexpr (Traits::CACHE_DODECAHEDRAL_MATRICES)
      std::memcpy(state->dodecahedral_matrices.data(),
                  lenses::DodecahedralTable::GROUP.values.data(),
                  sizeof(state->dodecahedral_matrices));
    if constexpr (Traits::CACHE_DODECAHEDRAL_FACES)
      lenses::DodecahedralTable::unpack_faces(state->dodecahedral_faces);
    use_parameter_storage(persistent_arena,
                          persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    clocks.template seed<Derived>(state->resources);
    palette.init();
    if constexpr (ANIMATED_PROJECTION) {
      constexpr int32_t SEED = [] {
        if constexpr (AUTHORED_TOPOLOGY &&
                      Metadata::PROJECTION != ProjectionKind::NONE)
          return static_cast<int32_t>(
              instance_hash(authored_instance_label("projection"),
                            projection_operator_id(Metadata::PROJECTION)));
        return PROJECTION_WALK_SEED;
      }();
      auto walk = Animation::RandomWalk<W>(
          this->projection_walk, math::UP, state->projection_walk_noise,
          typename Animation::RandomWalk<W>::Options{}, SEED);
      if constexpr (AUTHORED_TOPOLOGY)
        state->projection_walk_noise.SetSeed(SEED);
      if constexpr (requires { params.template get<"projection">().frame; })
        timeline.add_pausable(0, std::move(walk), &projection_frame_paused);
      else
        timeline.add(0, std::move(walk));
    }
    if constexpr (Metadata::OUTER_CAMERA)
      timeline.add(0, Animation::RandomWalk<W>(
                          this->outer_walk, math::UP, state->outer_walk_noise,
                          typename Animation::RandomWalk<W>::Options{},
                          CAMERA_WALK_SEED));
    this->register_parameters();
#ifdef ARDUINO
    if constexpr (SpecT::CACHE_FULL_GAMUT_LUT)
      init_gamut_lut(persistent_arena, GAMUT_LUT_ANGLE_STEPS,
                     GAMUT_LUT_L_STEPS);
#endif
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
    if constexpr (ANIMATED_PROJECTION &&
                  requires { params.template get<"projection">().frame; })
      projection_frame_paused = !projection_frame_active();
    {
      HS_PROFILE(fx_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(fx_advance);
      step_choreography();
      advance_runtime();
      update_spatial_frames();
      palette.step();
    }
    prepare_frame(state->current);
    const FrameState &current = state->current;
    if constexpr (TEMPORAL_SAMPLING) {
      IntervalCache::begin(*state, current);
      if (!Spacetime::flat_surface()) {
        const auto &sampler = IntervalCache::bake(
            *state, state->current, [](const Params &value) {
              return palette_cache_controls_eligible(value);
            });
        HS_PROFILE(fx_shader_draw);
        Scan::Shader::draw_column_times<W, H>(
            canvas,
            [&](int c) -> decltype(auto) {
              HS_PROFILE(fx_prepare_column);
              sampler.at_into(state->column,
                              Spacetime::SamplingSurface<W>::time(c));
              if constexpr (SpecT::PERSISTENT_PREPARED_FRAME)
                return prepare_persistent_frame(state->column);
              else
                return RenderPipeline::prepare(state->column);
            },
            [](const Frame &frame, int x) {
              if constexpr (requires {
                              Derived::prepare_column_shader(frame, x);
                            }) {
                return Derived::prepare_column_shader(frame, x);
              } else if constexpr (factored_columns()) {
                return column_shader(frame, x);
              } else {
                return [&frame](const math::Vector &view) {
                  return Derived::shade(view, frame);
                };
              }
            });
      } else {
        draw_current_columns(canvas, current);
      }
      IntervalCache::capture(*state, current);
    } else {
      draw_current_columns(canvas, current);
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

  /** @brief Whether the default shader supports a factored column prefix.
   * @return True when the pipeline and effect permit factoring. */
  static constexpr bool factored_columns() {
    return SpecT::FACTORED_COLUMNS && RenderPipeline::LINEAR_COLUMNS &&
           &Derived::shade == &ComposedEffect::shade;
  }

  /** @brief Prepares the factored shader for one column.
   * @param frame Prepared frame to shade.
   * @param x Column index.
   * @return A latitude shader retaining the frame and column basis. */
  static HS_HOT_INLINE auto column_shader(const Frame &frame, int x) {
    return [&frame, basis = RenderPipeline::column_basis(
                        frame, math::TrigLUT<W, H>::cos_theta(x),
                        math::TrigLUT<W, H>::sin_theta[x])](float sine,
                                                            float cosine) {
      return RenderPipeline::shade_column(sine, cosine, frame, basis);
    };
  }

  /** @brief Draws a fixed frame using a factored prefix when available.
   * @param canvas Destination canvas.
   * @param frame Prepared frame to shade. */
  static void draw_static_columns(Canvas &canvas, const Frame &frame) {
    if constexpr (requires { Derived::prepare_column_shader(frame, 0); }) {
      Scan::Shader::draw_columns<W, H>(canvas, [&](int x) {
        return Derived::prepare_column_shader(frame, x);
      });
    } else if constexpr (factored_columns()) {
      Scan::Shader::draw_columns<W, H>(
          canvas, [&](int x) { return column_shader(frame, x); });
    } else {
      Scan::Shader::draw_cached<W, H, 1>(canvas,
                                         [&frame](const math::Vector &view) {
                                           return Derived::shade(view, frame);
                                         });
    }
  }

  /** @brief Replaces the persistent prepared frame from one context.
   * @param context Current column or fixed-frame context.
   * @return The current arena-owned prepared frame. */
  HS_FLASH_MEMBER const Frame &
  prepare_persistent_frame(const FrameState &context)
    requires(SpecT::PERSISTENT_PREPARED_FRAME)
  {
    return state->prepared.emplace(context);
  }

  /** @brief Prepares and shades one fixed context.
   * @param canvas Destination canvas.
   * @param current Frame context. */
  HS_FLASH_MEMBER void draw_current_columns(Canvas &canvas,
                                            const FrameState &current) {
    if constexpr (SpecT::PERSISTENT_PREPARED_FRAME) {
      const Frame &frame = prepare_persistent_frame(current);
      HS_PROFILE(fx_shader_draw);
      draw_static_columns(canvas, frame);
    } else {
      const Frame frame = RenderPipeline::prepare(current);
      HS_PROFILE(fx_shader_draw);
      draw_static_columns(canvas, frame);
    }
  }

  /** @brief Whether a parameter set is admissible, family by family.
   *  @param params Parameter set to check.
   *  @return True when every family is in range. */
  static bool valid_params(const Params &params) {
    bool admitted = Pullback::valid(params);
    if constexpr (!Metadata::EDGE_DISTANCE_AVAILABLE ||
                  Metadata::PROJECTION == ProjectionKind::BONNE ||
                  Metadata::PROJECTION == ProjectionKind::AIROCEAN)
      params.visit([&]<typename Resource>(const auto &family) {
        if constexpr (!Metadata::EDGE_DISTANCE_AVAILABLE &&
                      requires { family.coverage_mode; })
          admitted = admitted && family.coverage_mode <
                                     static_cast<uint8_t>(
                                         ProjectionCoverageMode::EDGE_FADE);
        if constexpr (Resource::KIND == ResourceKind::PROJECTION &&
                      std::is_base_of_v<AtlasProjectionParams,
                                        typename Resource::Family>)
          admitted = admitted &&
                     valid_projection_topology(Metadata::PROJECTION, family);
      });
    return admitted;
  }
#if HS_ENABLE_PARAM_GUI_BRIDGE
  /**
   * @brief Warning text for a slider whose last write was refused.
   * @param name Parameter name.
   * @return The family warning for the last refused parameter, else null.
   */
  const char *parameter_warning(const char *name) const override {
    return this->refused_parameter_warning(name);
  }
#endif
#if HS_ENABLE_TEST_HOOKS
  /** @brief Prepares and returns this draw's frame state.
      @return The prepare_frame() result. */
  FrameState frame_for_test() {
    FrameState frame{};
    prepare_frame(frame);
    return frame;
  }
#endif

protected:
  /// Descriptors the arena-backed parameter array holds.
  static constexpr size_t PARAM_CAPACITY = Traits::PARAM_CAPACITY;

  using Choreography::anims_paused;
  using Choreography::params;

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /**
   * @brief Refuses a write that would make its parameter family invalid.
   * @param parameter Descriptor being written.
   * @param value Proposed value.
   * @return False when the write is refused.
   */
  bool parameter_write_admitted(const ParamDef &parameter,
                                float value) override {
    return this->admit_parameter_write(parameter, value);
  }
#endif

  using Choreography::register_animated_param;
  using Choreography::step_choreography;
  using Choreography::timeline;
  using Choreography::transition;
  using Choreography::use_parameter_storage;

  /** @brief Adopts a snap target and re-derives the palette mapping weights.
   *  @param target Parameters to adopt. */
  HS_FLASH_MEMBER void adopt_params(const Params &target) {
    params = target;
    palette.adopt(target.template get<"color">());
    this->refresh_modulation_controls();
  }

  /** @brief Re-derives the palette mapping weights after a slider write. */
  HS_FLASH_MEMBER void parameter_written() override {
    Choreography::parameter_written();
    palette.adopt(this->params.template get<"color">());
    this->refresh_modulation_controls();
  }

  /** @brief Captures the palette-mapping endpoints of an arming crossfade.
   *  @param target Crossfade's destination parameters. */
  HS_FLASH_MEMBER void transition_armed(const Params &target) {
    palette.arm(target.template get<"color">());
  }

  /**
   * @brief Writes the interpolated parameters of an in-flight transition.
   * @param progress Eased transition progress, 0 to 1.
   */
  HS_FLASH_MEMBER void blend_params(float progress) {
    Pullback::interpolate_into(params, transition.from, transition.to,
                               progress);
    palette.blend(progress);
    this->refresh_modulation_controls();
  }

  /** @brief Starts the lens's repeating phase clock.
   * @param scale Horizontal generator multiplier.
   * @param duration Frames per revolution. */
  template <ResourceKey Key = "lens">
  HS_FLASH_MEMBER void start_lens_phase(float scale, int duration) {
    if constexpr (Params::template HAS<Key>) {
      auto &phase = clocks.template get<Key>();
      timeline.add_pausable(
          0,
          Animation::LensPhase(phase.angle, phase.magnitude, scale, duration),
          &anims_paused);
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
    return Hooks{};
  }

private:
  static consteval bool modulation_presets_are_supported() {
    bool supported = true;
    Choreography::preset_row(0).params.visit(
        [&]<typename Resource>(const auto &family) {
          if constexpr (Resource::KIND == ResourceKind::MODULATION)
            if constexpr (!Resource::Family::DYNAMIC)
              for (size_t index = 1;
                   index < Choreography::authored_preset_count(); ++index) {
                const auto config = Choreography::preset_row(index)
                                        .params.template get<Resource::KEY>()
                                        .modulation;
                supported =
                    supported && config.source == family.modulation.source &&
                    config.mechanic.kind == family.modulation.mechanic.kind;
              }
        });
    return supported;
  }

  struct State
      : ProjectionWalkNoise<ANIMATED_PROJECTION>,
        OuterCameraNoise<Metadata::OUTER_CAMERA>,
        OptionalHueRotationLut<HUE_ROTATION_LUT>,
        IntervalCache::LabLuts,
        OptionalDodecahedralMatrices<Traits::CACHE_DODECAHEDRAL_MATRICES>,
        OptionalDodecahedralFaces<Traits::CACHE_DODECAHEDRAL_FACES>,
        IntervalCache::RgbLuts,
        IntervalCache::History,
        Registration::Storage {
    typename Clocks::Storage resources;
    std::conditional_t<AUTHORED_PALETTE, AuthoredPaletteResources,
                       std::monostate>
        authored_palette;
    FrameState current{};
    std::conditional_t<SpecT::PERSISTENT_PREPARED_FRAME, std::optional<Frame>,
                       std::monostate>
        prepared;
    typename IntervalCache::Column column;
    typename IntervalCache::Interval interval;
  };

  /// Persistent-arena bytes one instance claims.
  static constexpr size_t FOOTPRINT_BYTES =
      Traits::template footprint_bytes<State>();
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "Pullback::ComposedEffect persistent footprint exceeds the default "
      "partition");

  /** @brief Registered flexible-modulation names.
   * @return Arena-owned names, or null before init(). */
  typename Registration::Storage *modulation_control_storage() { return state; }

  /** @brief Whether the current authored projection frame advances.
   * @return False for an identity projection frame. */
  bool projection_frame_active() const {
    if constexpr (requires { params.template get<"projection">().frame; })
      return params.template get<"projection">().frame != 0;
    return true;
  }

  /** @brief Whether the projection frame is the tilted spin-wander frame.
   * @return The authored frame value 2 when the frame is a parameter, else
   *         `TILTED_PROJECTION_FRAME`. */
  bool projection_frame_tilted() const {
    if constexpr (requires { params.template get<"projection">().frame; })
      return params.template get<"projection">().frame == 2;
    return Traits::TILTED_PROJECTION_FRAME;
  }

  /** @brief Steps the resource clocks, the projection and camera spins and
      the palette wobble. */
  HS_FLASH_MEMBER void advance_runtime() {
    clocks.advance();
    if constexpr (ANIMATED_PROJECTION)
      if (projection_frame_active())
        this->projection_spin = advance_frame_spin(
            this->projection_spin,
            this->params.template get<"projection">().spin_rate);
    if constexpr (Metadata::OUTER_CAMERA &&
                  requires { Derived::CAMERA_SPIN_RATE; })
      this->camera_spin =
          advance_frame_spin(this->camera_spin, Derived::CAMERA_SPIN_RATE);
    palette.advance();
  }

  HS_FLASH_MEMBER void update_spatial_frames() {
    if constexpr (ANIMATED_PROJECTION)
      if (projection_frame_active())
        this->update_projection_frame(
            this->params.template get<"projection">().wander,
            projection_frame_tilted());
    if constexpr (Metadata::OUTER_CAMERA) {
      constexpr bool SPIN = requires { Derived::CAMERA_SPIN_RATE; };
      this->template update_camera_frame<SPIN>(
          this->params.template get<"projection">().camera_wander);
    }
  }

  /** @brief Checks the authored cache qualification domain.
   * @param value Frame parameters.
   * @return Whether the color controls admit interval caching. */
  static bool palette_cache_controls_eligible(const Params &value) {
    if constexpr (SpecT::INITIAL_CACHE_CONTROLS_ONLY) {
      const auto initial = Derived::initial_params();
      const bool color_matches =
          value.template get<"color">() == initial.template get<"color">();
      if constexpr (Params::template HAS<"modulation">)
        return color_matches &&
               value.template get<"modulation">().modulation ==
                   initial.template get<"modulation">().modulation;
      else
        return color_matches;
    } else
      return true;
  }

  /**
   * @brief Bakes the frame's LUTs and snapshots everything the scan reads.
   * @details Temporal draws retain palette and source endpoints through the
   * scan.
   * @param frame Receives this draw's frame state.
   */
  HS_FLASH_MEMBER void prepare_frame(FrameState &frame) {
    HS_PROFILE(fx_prepare_frame);
    palette.prepare(frame, *state);
    frame.projection_conjugate = this->frame_conjugate();
    frame.outer_conjugate = this->camera_conjugate();
    frame.palette_lab_lut = nullptr;
    if constexpr (Traits::CACHE_DODECAHEDRAL_MATRICES)
      frame.dodecahedral_matrices = state->dodecahedral_matrices.data();
    else
      frame.dodecahedral_matrices = nullptr;
    if constexpr (Traits::CACHE_DODECAHEDRAL_FACES)
      frame.dodecahedral_faces = state->dodecahedral_faces.data();
    else
      frame.dodecahedral_faces = nullptr;
    frame.params = params;
    frame.resources = {};
    clocks.prepare_views(frame);
    if constexpr (TEMPORAL_SAMPLING) {
      frame.temporal = {frame.palette, frame.hue_rotation_lut, 1.0f};
      frame.temporal.hue_rotation_gamut = frame.hue_rotation_gamut;
    }
  }

  State *state = nullptr;
  Clocks clocks{this->params};
  bool projection_frame_paused = false;
  ComposedPalette<Traits> palette{this->params.template get<"color">()};
};

} // namespace Pullback
