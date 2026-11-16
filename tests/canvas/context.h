/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Unit tests for core/render/canvas.h and core/control/transition.h.
 *
 * A Canvas spins in its ctor while !buffer_free(), so a sequential drawing test
 * must advance_display() before the next Canvas.
 */
#pragma once

/** @file
 * @brief Shared canvas test fixtures. */

#include <cmath>
#include <tuple>
#include <type_traits>
#include <atomic>
#include <array>
#include <cstring>
#include <limits>
#include <string_view>
#include <thread>
#include <vector>

#include "core/control/transition.h"
#include "core/render/canvas.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::canvas_tests {

enum class TestMode { NONE, WARP, SPARKLE };
enum class TestModeI8 : int8_t { NONE, WARP, SPARKLE };
enum class TestModeU8 : uint8_t { NONE, WARP, SPARKLE };
enum class TestModeI16 : int16_t { NONE, WARP, SPARKLE };
enum class TestModeU16 : uint16_t { NONE, WARP, SPARKLE };
enum class TestModeU32 : uint32_t { NONE, WARP, SPARKLE };

/**
 * @brief Minimal concrete Effect for exercising the base-class machinery.
 * @details Exposes the protected hooks (register_param / persist_pixels) needed
 * by the tests.
 */
struct TestEffect : public Effect {
  using Effect::mark_global;
  using Effect::mark_readonly;
  using Effect::register_param;
  using Effect::register_int_param;

  float speed = 1.5f; /**< Sample float param backing store. */
  bool flag = false;  /**< Sample bool param backing store. */

  /**
   * @brief Constructs the test effect with the given canvas dimensions.
   * @param W Canvas width in pixels.
   * @param H Canvas height in pixels.
   * @param cfg Effect configuration forwarded to the base.
   */
  TestEffect(int W, int H, EffectConfig cfg = {}) : Effect(W, H, cfg) {}
  /**
   * @brief Per-frame draw hook; a no-op.
   */
  void draw_frame() override {}

  /**
   * @brief Toggles frame-to-frame pixel persistence.
   * @param b True to copy the previous frame into each new buffer.
   */
  void set_persist(bool b) { persist_pixels = b; }
  /**
   * @brief Registers a float parameter with the base Effect.
   * @param n Parameter name.
   * @param p Pointer to the float backing store.
   * @param mn Minimum allowed value.
   * @param mx Maximum allowed value.
   */
  void add_float(const char *n, float *p, float mn, float mx) {
    register_param(n, p, mn, mx);
  }
  /**
   * @brief Registers a bool parameter after setting its initial value.
   * @param n Parameter name.
   * @param p Pointer to the bool backing store.
   * @param d Initial value written into *p before registration.
   */
  void add_bool(const char *n, bool *p, bool d) {
    *p = d;
    register_param(n, p);
  }
  /**
   * @brief Registers an enumerated parameter with the base Effect.
   * @param n Parameter name.
   * @param p Pointer to the float backing store holding the option index.
   * @param options Option label array (static storage).
   * @param count Number of labels.
   */
  void add_enum(const char *n, float *p, const char *const *options,
                int count) {
    register_param(n, p, options, count);
  }
  template <typename Enum>
  void add_typed_enum(const char *n, Enum *p, const char *const *options,
                      const char *const *export_options, int count) {
    register_param(n, p, options, export_options, count);
  }
  template <typename State>
  void mirror_display(const State &requested, const State &displayed) {
    mirror_parameter_display_state(requested, displayed);
  }
  void clear_params() { reset_parameters(); }
};

template <typename Host>
concept named_registration_policies =
    requires(Host &host, float *f, bool *b, int *i, TestMode *mode) {
      host.register_param("float", f, ParamSpec<float>{.animated = true});
      host.register_param("bool", b, ParamSpec<bool>{.readonly = true});
      host.register_param("integer", i, ParamSpec<int>{.preset = false});
      host.register_param("enum", mode, ParamSpec<TestMode>{.animated = true});
    };
template <typename Host>
concept positional_float_policy = requires(Host &host, float *value) {
  host.register_param("float", value, 0.0f, 1.0f, true);
};
template <typename Host>
concept positional_bool_policy = requires(Host &host, bool *value) {
  host.register_param("bool", value, true);
};
template <typename Host>
concept positional_integer_policy = requires(Host &host, int *value) {
  host.register_int_param("integer", value, 0, 1, true);
};
template <typename Host>
concept positional_enum_policy =
    requires(Host &host, TestMode *value, const char *const *labels) {
      host.register_param("enum", value, labels, labels, 3, true);
    };

struct PresetTestEffect : public Effect {
  PresetTestEffect() : Effect(8, 8) { configure_presets(3); }

  void draw_frame() override {}

  bool advance() { return advance_preset(); }

  bool accept = true;
  size_t applied_from = 0;
  size_t applied_to = 0;
  bool applied_manually = false;
  bool applied_synchronized = false;
  size_t committed_index = 0;
  int commit_count = 0;

private:
  bool apply_preset(const PresetChange &change) override {
    applied_from = change.from;
    applied_to = change.to;
    applied_manually = change.origin == PresetChangeOrigin::MANUAL;
    applied_synchronized = change.origin == PresetChangeOrigin::SYNCHRONIZED;
    return accept;
  }

  void preset_changed(const PresetChange &) override {
    committed_index = getPresetIndex();
    ++commit_count;
  }
};

inline void prime_buffers(TestEffect &fx, Pixel a, Pixel b) {
  for (Pixel color : {a, b}) {
    {
      Canvas canvas(fx);
      for (int i = 0; i < fx.width() * fx.height(); ++i)
        canvas(i) = color;
    }
    fx.advance_display();
  }
}

struct StubPipe {
  int integer_hits = 0;
  int float_hits = 0;
  int integer_x = 0;
  int integer_y = 0;
  float float_x = 0.0f;
  float float_y = 0.0f;

  void plot(Canvas &, int x, int y, const Pixel &, float, float) {
    ++integer_hits;
    integer_x = x;
    integer_y = y;
  }
  void plot(Canvas &, float x, float y, const Pixel &, float, float) {
    ++float_hits;
    float_x = x;
    float_y = y;
  }
  void plot(Canvas &, const math::Vector &, const Pixel &, float, float) {}
};
} // namespace hs_test::canvas_tests
