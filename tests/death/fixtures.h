/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Death-test effect fixture.

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
   * @brief Constructs the effect with construction-time flags.
   * @param width Canvas width in pixels.
   * @param height Canvas height in pixels.
   * @param cfg Effect flags, as a pipeline's pipeline_config() folds them.
   */
  DeathEffect(int width, int height, EffectConfig cfg)
      : Effect(width, height, cfg) {}
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
    register_param("options", p,
                   ParamSpec<float>{.options = options, .option_count = count});
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
