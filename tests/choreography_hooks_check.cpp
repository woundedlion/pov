/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/render/pullback/composed_effect.h"

struct HookParams {
  float value = 0.0f;
#ifdef HS_HOOK_DISCRETE_LERP
  uint8_t mode = 0;
  static constexpr auto FIELDS = std::tuple{
      Control::Field<HookParams, uint8_t>{.id = "mode",
                                          .member = &HookParams::mode,
                                          .spec = {.min = 0, .max = 1},
                                          .curve = Control::FieldCurve::LERP}};
#endif
};
struct HookFixture : ChoreographedEffect<HookFixture, HookParams> {
  using Base = ChoreographedEffect<HookFixture, HookParams>;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
#ifdef HS_HOOK_DEPARTURE_OUTLASTS
  static constexpr uint16_t PRESET_DWELL_FRAMES = 1;
#else
  static constexpr uint16_t PRESET_DWELL_FRAMES = 8;
#endif
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS{
      {{{0.0f}, Segue::Preset::Lerp{2, math::ease_linear}},
       {{1.0f}, Segue::Preset::Fade{2}}}};
#ifdef HS_HOOK_IDS_SIZE_MISMATCH
  static constexpr std::array<const char *, 1> PRESET_IDS{"a"};
#elif defined(HS_HOOK_PRESET_ROW) && !defined(HS_HOOK_PRESET_ROW_NO_IDS)
  static constexpr std::array<const char *, 2> PRESET_IDS{"a", "b"};
#endif
  static constexpr auto preset(size_t index) { return PRESETS[index]; }
  HookFixture() : Base(8, 8) {}
  static bool valid_params(const Params &) { return true; }
  static constexpr Params initial_params() { return {}; }
  Params preset_params(size_t index) const { return PRESETS[index].params; }
#ifndef HS_HOOK_INHERITED_DEFAULT
  void adopt_params(const Params &value) { params = value; }
#endif
  void blend_params(float) {}
  void transition_armed(const Params &) {}
#ifdef HS_HOOK_BAD_RETURN
  int finish_blend(const Params &) { return 0; }
#elif defined(HS_HOOK_BAD_ARGUMENT)
  void finish_blend(Params) {}
#else
  void finish_blend(const Params &) {}
#endif
  void set_preset_opacity(float) {}
  void after_profile_preset() {}
  void after_composed_init() {}
  void init() override { begin_choreography(); }
  void draw_frame() override { step_choreography(); }
  static consteval auto choreography_hooks() {
    using Hook = Control::ChoreographyHook;
    return Control::HookSet<
        HookFixture,
#ifdef HS_HOOK_PRESET_ROW
        Control::HookBinding<Hook::PRESET_ROW, &HookFixture::preset>,
#endif
#ifdef HS_HOOK_VALIDATION_DEFAULT
        Control::HookBinding<Hook::VALID_PARAMS, &Base::valid_params>,
#else
        Control::HookBinding<Hook::VALID_PARAMS, &HookFixture::valid_params>,
#endif
        Control::HookBinding<Hook::INITIAL_PARAMS,
                             &HookFixture::initial_params>,
        Control::HookBinding<Hook::PRESET_PARAMS, &HookFixture::preset_params>,
        Control::HookBinding<Hook::ADOPT_PARAMS, &HookFixture::adopt_params>,
#ifndef HS_HOOK_LERP_NO_BLEND
        Control::HookBinding<Hook::BLEND_PARAMS, &HookFixture::blend_params>,
#endif
        Control::HookBinding<Hook::TRANSITION_ARMED,
                             &HookFixture::transition_armed>,
#ifdef HS_HOOK_BOUND_TYPO
        Control::HookBinding<Hook::FINISH_BLEND, &HookFixture::blend_params>,
#else
        Control::HookBinding<Hook::FINISH_BLEND, &HookFixture::finish_blend>,
#endif
#ifdef HS_HOOK_DUPLICATE
        Control::HookBinding<Hook::FINISH_BLEND, &HookFixture::finish_blend>,
#endif
#ifndef HS_HOOK_FADE_NO_OPACITY
        Control::HookBinding<Hook::PRESET_OPACITY,
                             &HookFixture::set_preset_opacity>,
#endif
        Control::HookBinding<Hook::AFTER_PROFILE_PRESET,
                             &HookFixture::after_profile_preset>,
        Control::HookBinding<Hook::AFTER_COMPOSED_INIT,
                             &HookFixture::after_composed_init>>{};
  }
};

struct DefaultHookSpec : Pullback::Spec {
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Project<
          Pullback::Projection::Stereographic<Pullback::ProjectionProvider<B>>>,
      Pullback::Stage::Sample<Pullback::Source::Grid<
          Pullback::SourceProvider<B, Pullback::GridSourceParams>>>,
      Pullback::Stage::Colorize<
          Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
              B, Pullback::Color::BrightnessEnvelope::NONE>>>>;
};
struct AuthoredHookFixture
    : Pullback::ComposedEffect<8, 8, AuthoredHookFixture, DefaultHookSpec> {
  using Base =
      Pullback::ComposedEffect<8, 8, AuthoredHookFixture, DefaultHookSpec>;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 8;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
#ifdef HS_HOOK_FACTORY_TYPO
#define HS_HOOK_FACTORY_NAME choregraphy_hooks
#else
#define HS_HOOK_FACTORY_NAME choreography_hooks
#endif
  static consteval auto HS_HOOK_FACTORY_NAME() {
    return Base::choreography_hooks().template extend<AuthoredHookFixture>();
  }
#undef HS_HOOK_FACTORY_NAME
};
struct InheritedHookFixture
    : Pullback::ComposedEffect<8, 8, InheritedHookFixture, DefaultHookSpec,
                               Control::HookOwnership::INHERITED> {
  static constexpr uint16_t PRESET_DWELL_FRAMES = 8;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
};

static_assert(std::is_empty_v<decltype(HookFixture::choreography_hooks())>);

void construct_hook_companions() {
  HookFixture hooks;
  AuthoredHookFixture authored;
  InheritedHookFixture inherited;
}

#ifdef HS_HOOK_ROW_FIELDS
struct RowParams {
  float value = 0;
  static constexpr auto FIELDS = std::tuple{Control::Field<RowParams>{
      "value", &RowParams::value, "Value", {0.0f, 1.0f}}};
};
struct RowFixture : ChoreographedEffect<RowFixture, RowParams> {
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 8;
  static constexpr std::array<PresetEntry<Params>, 1> PRESETS{{
#ifdef HS_HOOK_ROW_OUT_OF_RANGE
      {{2.0f}},
#else
      {{0.5f}},
#endif
  }};
  RowFixture() : ChoreographedEffect(8, 8) {}
  void draw_frame() override {}
  static consteval auto choreography_hooks() {
    return Control::HookSet<RowFixture>{};
  }
};
[[maybe_unused]] RowFixture row_fixture;
#endif
