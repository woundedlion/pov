/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <type_traits>
#include "core/control/params.h"
#include "core/memory.h"
#include "core/render/pullback/composed_descriptors.h"

/** @file composed_registration.h
 * @brief `ComposedRegistration`, the slider registration of a composed effect.
 */

namespace Pullback {

namespace ComposedDetail {

/** @brief Descriptor type retained by a parameter family's scalar field table.
 * @tparam Family Parameter family. */
template <typename Family>
using ScalarField =
    typename std::remove_cvref_t<decltype(Family::FIELDS)>::value_type;

/** @brief Visits named scalar registrations admitted by a family's gate mask.
 * @tparam Family Parameter family.
 * @param family Live parameter values.
 * @param gates Bit per admitted field gate.
 * @param context Registration callback context.
 * @param register_field Registers one admitted descriptor and its target. */
template <typename Family>
HS_SIZE_FLASH_MEMBER inline void register_scalar_family(
    Family &family, uint32_t gates, void *context,
    void (*register_field)(void *, const ScalarField<Family> &, float *,
                           const ParamSpec<float> &)) {
  static_assert(
      [] {
        for (const auto &field : Family::FIELDS)
          if (static_cast<unsigned>(field.gate) >=
              std::numeric_limits<uint32_t>::digits)
            return false;
        return true;
      }(),
      "family registration gates must fit the mask");
  for (const auto &field : Family::FIELDS)
    if (field.name != nullptr &&
        (gates & (uint32_t{1} << static_cast<unsigned>(field.gate)))) {
      const auto spec = field.description().spec;
      register_field(context, field, &(family.*field.member), spec);
    }
}

/// Whether the resource is a modulation recipe that permits topology changes.
template <typename Resource>
inline constexpr bool DYNAMIC_MODULATION = [] {
  if constexpr (Resource::KIND == ResourceKind::MODULATION)
    return Resource::Family::DYNAMIC;
  return false;
}();

/** @brief Registered names of one resource's flexible modulation controls.
 * @tparam Resource Keyed resource.
 * @tparam Enabled Whether the resource is a flexible modulation recipe. */
template <typename Resource, bool Enabled = DYNAMIC_MODULATION<Resource>>
struct ModulationControlNames {};
/** @brief Names used to refresh active-control flags.
 * @tparam Resource Flexible modulation resource. */
template <typename Resource> struct ModulationControlNames<Resource, true> {
  std::array<const char *, Modulation::Controls::EXTRA_SCHEMA.size()>
      modulation_control_names{}; ///< Persistent registered names.
};

} // namespace ComposedDetail

/**
 * @brief Registers a composed effect's sliders, folders and modulation
 *        controls.
 * @details CRTP base of `Effect`, which befriends it and supplies
 * `modulation_control_storage()`.
 * @tparam Effect The composed effect.
 * @tparam Traits The effect's `ComposedTraits`.
 */
template <typename Effect, typename Traits> class ComposedRegistration {
  using Params = typename Traits::Params;
  using SpecT = typename Traits::Spec;
  using Metadata = typename Traits::Metadata;
  using ColorFamily = typename Traits::ColorFamily;
  static constexpr bool AUTHORED_PALETTE = Traits::AUTHORED_PALETTE;

public:
  /// Arena-owned registered names of flexible modulation controls.
  using Storage =
      ComposedDetail::ResourceStorage<ComposedDetail::ModulationControlNames,
                                      typename Params::ResourceTypes>;

protected:
  /** @brief The effect this base registers for.
   * @return The derived effect. */
  Effect &effect() { return static_cast<Effect &>(*this); }

  /**
   * @brief Slider name for a field of one resource instance.
   * @tparam Resource The instance the field belongs to.
   * @param name Unqualified field name.
   * @return `name`, or a copy prefixed with the instance key and a dot,
   *         allocated from `persistent_arena`, when the instance is
   *         qualified.
   */
  template <typename Resource>
  HS_SIZE_FLASH_MEMBER const char *resource_parameter_name(const char *name) {
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
        {{"Camera", "Camera", 0},
         {"Lens", "Lens", 1},
         {"Surface Noise", "Surface Noise", 2},
         {"Projection Frame", "Projection Frame", 3},
         {"Projection", "Projection", 4},
         {"Planar Warp 1", "Planar Warp 1", 5},
         {"Planar Warp 2", "Planar Warp 2", 6},
         {"Function", "Function", 7},
         {"Value Transfer", "Value Transfer", 8},
         {"Coverage", "Coverage", 9},
         {"Colorize", "Colorize", 10}}};
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
    effect().set_last_parameter_group(
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
    if constexpr (Resource::KIND == ResourceKind::ROTATION)
      return ParameterStage::CAMERA;
    else if constexpr (Resource::KIND == ResourceKind::SOURCE)
      return ParameterStage::FUNCTION;
    else if constexpr (Resource::KIND == ResourceKind::SURFACE)
      return ParameterStage::SURFACE;
    else if constexpr (Resource::KIND == ResourceKind::WARP)
      return Resource::KEY.view() == "outer_warp" ? ParameterStage::WARP_1
                                                  : ParameterStage::WARP_2;
    else if constexpr (Resource::KIND == ResourceKind::LENS ||
                       Resource::KIND == ResourceKind::TWIST)
      return ParameterStage::LENS;
    else if constexpr (Resource::KIND == ResourceKind::VALUE)
      return (std::is_same_v<typename Resource::Family, IsoValueParams> ||
              std::is_same_v<typename Resource::Family, SmoothBandsParams>)
                 ? ParameterStage::VALUE_TRANSFER
                 : ParameterStage::COVERAGE;
    else
      return ParameterStage::PROJECTION;
  }

  /** @brief Registers one admitted scalar field and assigns its pipeline group.
   * @tparam Resource Parameter resource.
   * @tparam T Parameter family.
   * @param context Owning registration.
   * @param field Field descriptor.
   * @param target Live scalar.
   * @param spec Registration range and policies. */
  template <typename Resource, typename T>
  HS_SIZE_FLASH_MEMBER static void
  register_scalar_field(void *context,
                        const ComposedDetail::ScalarField<T> &field,
                        float *target, const ParamSpec<float> &spec) {
    auto &registration = *static_cast<ComposedRegistration *>(context);
    registration.effect().register_param(
        registration.template resource_parameter_name<Resource>(field.name),
        target, spec);
    ParameterStage stage = resource_parameter_stage<Resource>();
    if constexpr (Resource::KIND == ResourceKind::PROJECTION) {
      if (field.member == &T::camera_wander)
        stage = ParameterStage::CAMERA;
      else if (field.member == &T::spin_rate || field.member == &T::wander)
        stage = ParameterStage::PROJECTION_FRAME;
    }
    registration.set_parameter_stage(stage);
  }

  /** @brief Registers the named descriptors of one independent instance.
   *  @tparam Resource The instance's parameter resource.
   *  @tparam T Parameter family type.
   *  @param family Live family the sliders write. */
  template <typename Resource, typename T>
  HS_SIZE_FLASH_MEMBER void register_fields(T &family) {
    if constexpr (T::FIELDS.size() > 0) {
      constexpr uint32_t GATES = [] {
        uint32_t gates = 0;
        for (const auto &field : T::FIELDS)
          if (field_gate_open<Resource>(field.gate))
            gates |= uint32_t{1} << static_cast<unsigned>(field.gate);
        return gates;
      }();
      ComposedDetail::register_scalar_family(
          family, GATES, this, &register_scalar_field<Resource, T>);
    }
    if constexpr (requires { T::TOPOLOGY; })
      for (const auto &field : T::TOPOLOGY) {
        if (!ComposedDetail::resource_topology_open<SpecT, Resource>(field.id))
          continue;
        std::span<const char *const> options(field.value_ids,
                                             field.value_count);
        if constexpr (Resource::KIND == ResourceKind::PROJECTION)
          options = projection_topology_options(Metadata::PROJECTION, field.id,
                                                options);
        effect().register_animated_enum8_param(
            resource_parameter_name<Resource>(field.id),
            &(family.*field.member), options.data(), options.size());
        set_parameter_stage(resource_parameter_stage<Resource>());
      }
    if constexpr (requires {
                    T::EXTRA_SCHEMA;
                    family.extra_address(size_t{});
                  })
      for (size_t index = 0; index < T::EXTRA_SCHEMA.size(); ++index) {
        const auto &field = T::EXTRA_SCHEMA[index];
        const char *name = resource_parameter_name<Resource>(field.id);
        void *address = family.extra_address(index);
        if (field.enum_count)
          effect().register_animated_enum8_param(
              name, static_cast<uint8_t *>(address), field.enum_ids,
              field.enum_count);
        else
          effect().register_animated_param(name, static_cast<float *>(address),
                                           field.min, field.max);
        set_parameter_stage(resource_parameter_stage<Resource>());
      }
  }

  /** @brief Registers every slider of the effect's parameter families. */
  HS_SIZE_FLASH_MEMBER void register_parameters() {
    register_resource_kind<ResourceKind::ROTATION>();
    register_resource_kind<ResourceKind::SOURCE>();
    register_resource_kind<ResourceKind::PROJECTION>();
    register_resource_kind<ResourceKind::SURFACE>();
    register_resource_kind<ResourceKind::WARP>();
    register_resource_kind<ResourceKind::VALUE>();
    register_resource_kind<ResourceKind::LENS>();
    register_resource_kind<ResourceKind::TWIST>();
    register_color_field<&ColorParams::palette_chroma>("Palette Chroma");
    effect().register_animated_param(
        "Palette Mapping",
        &effect().params.template get<"color">().palette_mapping,
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
    if constexpr (AUTHORED_PALETTE) {
      auto &color = effect().params.template get<"color">();
      effect().register_animated_enum8_param(
          "Palette Mode", &color.palette_mode, Color::PALETTE_MODE_IDS,
          std::size(Color::PALETTE_MODE_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_enum8_param(
          "Palette Type", &color.palette_type, AUTHORED_PALETTE_TYPE_IDS,
          std::size(AUTHORED_PALETTE_TYPE_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_enum8_param(
          "Palette Harmony", &color.palette_harmony,
          AUTHORED_PALETTE_HARMONY_IDS,
          std::size(AUTHORED_PALETTE_HARMONY_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_enum8_param(
          "Palette Lightness Curve", &color.palette_lightness_curve,
          AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS,
          std::size(AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_enum8_param(
          "Palette Name", &color.palette_name, AUTHORED_PALETTE_NAME_IDS,
          std::size(AUTHORED_PALETTE_NAME_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_enum8_param(
          "Palette Name B", &color.palette_name_b, AUTHORED_PALETTE_NAME_IDS,
          std::size(AUTHORED_PALETTE_NAME_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      effect().register_animated_param(
          "Brightness Envelope", &color.brightness_envelope,
          Color::BRIGHTNESS_ENVELOPE_IDS, nullptr,
          std::size(Color::BRIGHTNESS_ENVELOPE_IDS));
      set_parameter_stage(ParameterStage::COLORIZE);
      for (size_t index = ColorParams::FIELDS.size();
           index < ColorFamily::FIELDS.size(); ++index) {
        const auto &field = ColorFamily::FIELDS[index];
        effect().register_param(field.name, &(color.*field.member),
                                field.description().spec);
        set_parameter_stage(ParameterStage::COLORIZE);
      }
    }
    effect().params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION)
        register_modulation<Resource>(family.modulation);
    });
  }

  /** @brief Refreshes the active source and variant controls after adoption. */
  HS_FLASH_MEMBER void refresh_modulation_controls() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const Storage *storage = effect().modulation_control_storage();
    if (!storage)
      return;
    effect().params.visit([&]<typename Resource>(const auto &family) {
      if constexpr (Resource::KIND == ResourceKind::MODULATION) {
        if constexpr (Resource::Family::DYNAMIC) {
          const auto &names =
              storage->template get<Resource::KEY>().modulation_control_names;
          if (!names[0])
            return;
          for (size_t index = 0; index < names.size(); ++index)
            effect().mark_readonly(names[index], !Modulation::Controls::active(
                                                     family.modulation, index));
        }
      }
    });
#endif
  }

private:
  /** @brief Whether a gated field's slider exists for one resource.
   *  @tparam Resource Parameter resource the field belongs to.
   *  @param gate Field gate.
   *  @return True when the field gets a slider. */
  template <typename Resource>
  static constexpr bool field_gate_open(Pullback::FieldGate gate) {
    return ComposedDetail::field_gate_open<SpecT, Resource>(gate);
  }

  template <typename Resource, typename T>
  HS_FLASH_MEMBER void register_warp_fields(T &warp) {
    static_assert(T::FIELDS[0].member == &T::speed &&
                      T::FIELDS[0].name == nullptr,
                  "warp speed must be the first unnamed descriptor");
    if constexpr (ComposedDetail::warp_speed_open<SpecT, Resource>()) {
      constexpr const char *SPEED_NAME =
          ComposedDetail::warp_speed_name(Resource::KEY.view()).data();
      effect().register_param(resource_parameter_name<Resource>(SPEED_NAME),
                              &warp.speed, T::FIELDS[0].description().spec);
      set_parameter_stage(resource_parameter_stage<Resource>());
    }
    register_fields<Resource>(warp);
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
  HS_SIZE_FLASH_MEMBER void register_color_field(const char *name) {
    constexpr const auto &FIELD = color_descriptor<Member>();
    if constexpr (AUTHORED_PALETTE ||
                  ComposedDetail::color_gate_open<SpecT>(FIELD.topology_gate)) {
      effect().register_animated_param(
          name, &(effect().params.template get<"color">().*Member), FIELD.min,
          FIELD.max);
      set_parameter_stage(ParameterStage::COLORIZE);
    }
  }

  template <ResourceKind Kind>
  HS_SIZE_FLASH_MEMBER void register_resource_kind() {
    effect().params.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == Kind) {
        if constexpr (Kind == ResourceKind::WARP)
          register_warp_fields<Resource>(family);
        else
          register_fields<Resource>(family);
      }
    });
  }

  template <typename Resource>
  HS_SIZE_FLASH_MEMBER void
  register_modulation(ColorModulation::Config &config) {
    if constexpr (Resource::Family::DYNAMIC) {
      auto &names = effect()
                        .modulation_control_storage()
                        ->template get<Resource::KEY>()
                        .modulation_control_names;
      for (size_t index = 0; index < Modulation::Controls::EXTRA_SCHEMA.size();
           ++index) {
        const auto &field = Modulation::Controls::EXTRA_SCHEMA[index];
        const char *name = resource_parameter_name<Resource>(field.id);
        names[index] = name;
        void *address = Modulation::Controls::address(config, index);
        if (field.enum_count)
          effect().register_animated_enum8_param(
              name, static_cast<uint8_t *>(address), field.enum_ids,
              field.enum_count);
        else
          effect().register_animated_param(name, static_cast<float *>(address),
                                           field.min, field.max);
        set_parameter_stage(ParameterStage::COLORIZE);
      }
      refresh_modulation_controls();
      return;
    }
    for (const auto &field : Modulation::control_fields(config)) {
      if (!field.value)
        continue;
      effect().register_animated_param(
          resource_parameter_name<Resource>(field.name), field.value, field.min,
          field.max);
      set_parameter_stage(ParameterStage::COLORIZE);
    }
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
};

} // namespace Pullback
