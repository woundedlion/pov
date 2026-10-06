/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// --- Individual death cases — each MUST trap (HS_CHECK / __builtin_trap) ------

// Params canvas death fixtures and guard cases.

/** @brief Rejects a margin below the filter pipeline requirement. */
inline void case_effect_margin_below_pipeline() {
  struct MarginEffect : Effect {
    MarginEffect() : Effect(32, 16, {.margin = 3, .required_margin = 3}) {}
    void draw_frame() override {}
  } effect;
  effect.set_margin(2);
}

/**
 * @brief Death case: reading a ParamDef with an unknown target type must trap.
 * @details An unknown tag has no supported value representation and traps
 *          before the descriptor reads the target.
 */
inline void case_param_def_unknown_get_target_type() {
  float storage = 0.5f;
  ParamDef def;
  def.target = &storage;
  def.target_type = static_cast<ParamDef::TargetType>(opaque<uint8_t>(9));
  if (def.get_from(&storage) == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: writing a ParamDef with an unknown target type must trap.
 * @details An unknown tag has no supported value representation and traps
 *          before the descriptor writes the target.
 */
inline void case_param_def_unknown_set_target_type() {
  float storage = 0.5f;
  ParamDef def;
  def.target = &storage;
  def.target_type = static_cast<ParamDef::TargetType>(opaque<uint8_t>(9));
  struct InternalWriter : ParamHost {
    using ParamHost::write_parameter_unchecked;
  };
  InternalWriter::write_parameter_unchecked(def, opaque(1.0f));
  if (storage == opaque(42.0f))
    std::printf("x");
}

/**
 * @brief Death case: a second simultaneously-live Effect must trap.
 * @details Every Effect aliases the same static framebuffers and double-buffer
 *          indices.
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
 * @brief Death case: a second simultaneously-live correction guard must trap.
 * @details NoColorCorrection and NoTempCorrection share one liveness flag and
 *          set the global FastLED correction/temperature.
 */
inline void case_correction_guard_double_construct() {
  NoColorCorrection a;
  NoColorCorrection b; // second live guard -> liveness HS_CHECK -> trap
  if (correction_guard_live() == opaque(true))
    std::printf("x");
}

/**
 * @brief Death case: a live NoColorCorrection plus a NoTempCorrection must trap.
 * @details The two guard types share one liveness flag.
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

/** @brief Death case: overflowing the fixed ParamList must trap. */
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

inline void case_param_spec_invalid_option_values() {
  DeathEffect fx;
  static constexpr const char *LABELS[] = {"Zero", "One", "Six"};
  static constexpr const int64_t IDS[] = {0, 1, 6};
  uint8_t valid = 6;
  uint8_t gap = opaque<uint8_t>(2);
  const ParamSpec<uint8_t> spec{.min = 0,
                                .max = 6,
                                .options = LABELS,
                                .option_count = 3,
                                .option_values = IDS};
  fx.register_param("valid", &valid, spec);
  fx.register_param("gap", &gap, spec);
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
 * @details arcs_overlap wraps the seam-relative offset with one conditional
 *          add, which assumes both starts are already reduced. Both lengths are
 *          positive and under w so the early-out branches do not preempt the
 *          guard.
 */
inline void case_arcs_overlap_start_out_of_range() {
  bool hit = ClipRegion::arcs_overlap(opaque(-1), opaque(2), opaque(0),
                                      opaque(2), opaque(8));
  if (hit)
    std::printf("x");
}

/** @brief Rejects a margin equal to the canvas width. */
inline void case_effect_margin_equal_width() {
  StubEffect effect(32, 16);
  effect.set_margin(opaque(32));
}
