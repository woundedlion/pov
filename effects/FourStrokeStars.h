/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/** @file FourStrokeStars.h
 * @brief Stationary Fibonacci stars sharing a four-stroke motion direction.
 */

#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "core/platform/effect_identity.h"
#include "core/render/scan/four_stroke.h"

/** @brief Scans five-point stars with a common wandering world-space light/dark edge cue.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels. */
template <int W, int H> class FourStrokeStars : public Effect {
public:
  /// Machine/display identity and compatibility keys.
  static constexpr hs::EffectIdentity EFFECT_IDENTITY =
      hs::effect_identities::FourStrokeStars::VALUE;
  /// Canonical machine ID.
  static constexpr std::string_view EFFECT_ID = EFFECT_IDENTITY.machine_id;

  /** @brief Constructs the W by H star field. */
  HS_FLASH_MEMBER FourStrokeStars()
      : Effect(W, H, {.replaces_frame = true, .can_skip_frames = true}) {}

  /** @brief Places the stationary stars and starts the direction walk. */
  HS_FLASH_MEMBER void init() override {
    register_described_params<&FourStrokeStars::parameter_fields>(params);
    for (int i = 0; i < STAR_COUNT; ++i)
      bases[i] = math::make_basis(math::Quaternion(1, 0, 0, 0),
                                  math::fib_spiral(STAR_COUNT, 0.5f, i));
    contours = persistent_arena.make<std::array<Contour, STAR_COUNT>>();
    prepare_contours();
    timeline.add(0, Animation::RandomWalk<W>(
                        direction, math::Y_AXIS, noise,
                        Animation::RandomWalk<W>::Options::Languid()));
  }

  /** @brief Shades the fixed footprints using the current global direction. */
  HS_HOT_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    const uint32_t due = canvas.scheduled_frame();
    for (uint32_t i = next_display_frame; i < due; ++i)
      timeline.step(canvas);
    next_display_frame = due + 1;
    timeline.step(canvas);
    if (north != math::DISPLAY_NORTH_PHI || south != math::DISPLAY_SOUTH_PHI)
      prepare_contours();
    static constexpr Pixel BASE_COLOR = FourStroke::linear_rgb(0xb5b5b5);
    const auto base = [](const math::Vector &, Fragment &fragment) {
      fragment.color = Color4(BASE_COLOR, 1);
    };
    const auto walk_arrow = direction.orient(math::Y_AXIS);
    const auto arrow = params.reverse ? -walk_arrow : walk_arrow;
    constexpr float LONGITUDE_STEP = math::TWO_PI_F / W;
    const float row_step = math::DisplayGeometry<H>::RADIANS_PER_ROW;
    const float pixel_angle = std::max(LONGITUDE_STEP, row_step);
    const Pixel background =
        params.black_background ? Pixel(0, 0, 0) : FourStroke::GRAY;
    const auto &clip = canvas.clip();
    for (int y = clip.y_start; y < clip.y_end; ++y)
      std::fill_n(canvas.data() + y * W + clip.x_start,
                  clip.x_end - clip.x_start, background);
    for (int i = 0; i < STAR_COUNT; ++i) {
      const SDF::Star shape(bases[i], RADIUS, 5, 0, false);
      constexpr float REFERENCE_TILE_PIXELS = 64;
      const float angular_width =
          params.reference ? 2 * shape.circumradius / REFERENCE_TILE_PIXELS
                           : 0.25f * shape.plane_d;
      const float rim_scale = std::min(1.f, angular_width / pixel_angle);
      const Scan::FourStrokeMaterial material{
          .controls = {.width = params.width * rim_scale,
                       .strength = params.strength,
                       .mode = params.reference
                                   ? FourStroke::Mode::REFERENCE_POLARITY
                                   : FourStroke::Mode::DIRECTIONAL_CONTRAST},
          .base = base,
          .arrow = arrow,
          .strip_direction =
              params.reference
                  ? std::optional{screen_direction(bases[i].v, arrow, row_step)}
                  : std::nullopt,
          .constant_base = BASE_COLOR};
      Scan::Star::draw<W, H>(pipeline, canvas, shape, (*contours)[i], material);
    }
  }

private:
  static constexpr int STAR_COUNT = 18; ///< Number of Fibonacci sites.
  static constexpr float RADIUS =
      0.20f; ///< Star radius in hemisphere fractions.
  using Contour =
      Scan::ProjectedStar<W, H>; ///< Static output-pixel star outline.
  std::array<Contour, STAR_COUNT> *contours =
      nullptr;                     ///< Prepared in arena storage.
  float north = 0, south = 0;      ///< Cached display latitude endpoints.
  uint32_t next_display_frame = 1; ///< Next one-based physical display window.

  static FourStroke::Direction
  screen_direction(math::Vector axis, math::Vector arrow, float row_step) {
    const float latitude_sine = std::hypot(axis.x, axis.z);
    if (latitude_sine == 0)
      return {};
    const math::Vector east{-axis.z / latitude_sine, 0, axis.x / latitude_sine};
    const float horizontal = math::dot(arrow, east);
    const float vertical = math::dot(arrow, math::cross(axis, east));
    const float x = horizontal / (latitude_sine * (math::TWO_PI_F / W));
    const float y = vertical / row_step;
    const float length = std::hypot(x, y);
    if (length == 0)
      return {};
    const float scale = std::hypot(horizontal, vertical) / length;
    return {x * scale, y * scale};
  }

  HS_FLASH_MEMBER void prepare_contours() {
    north = math::DISPLAY_NORTH_PHI;
    south = math::DISPLAY_SOUTH_PHI;
    for (int i = 0; i < STAR_COUNT; ++i)
      (*contours)[i].prepare(SDF::Star(bases[i], RADIUS, 5, 0, false));
  }

  struct Params {
    float strength = 1;            ///< Linear RGB mix.
    int width = 1;                 ///< Inner band width in output pixels.
    bool reference = true;         ///< Use the demo's opposing strips.
    bool reverse = false;          ///< Negate the common world-space direction.
    bool black_background = false; ///< Use a black backdrop.
  } params;                        ///< Live material controls.
  /** @brief Describes the exposed material controls.
   * @return Parameter field descriptors in display order. */
  static constexpr auto parameter_fields() {
    return std::tuple{
        Control::Field<Params, float>{.id = "strength",
                                      .member = &Params::strength,
                                      .name = "Strength",
                                      .spec = {.min = 0, .max = 1}},
        Control::Field<Params, int>{.id = "width",
                                    .member = &Params::width,
                                    .name = "Edge Width",
                                    .spec = {.min = 1, .max = 2}},
        Control::Field<Params, bool>{.id = "reference",
                                     .member = &Params::reference,
                                     .name = "Reference Polarity",
                                     .spec = {}},
        Control::Field<Params, bool>{.id = "reverse",
                                     .member = &Params::reverse,
                                     .name = "Reverse Direction",
                                     .spec = {}},
        Control::Field<Params, bool>{.id = "black_background",
                                     .member = &Params::black_background,
                                     .name = "Black Background",
                                     .spec = {}}};
  }
  std::array<math::Basis, STAR_COUNT>
      bases;                     ///< Fixed star frames in sphere coordinates.
  math::Orientation<> direction; ///< Pose rotating the shared cue direction.
  FastNoiseLite noise;           ///< Noise source for the direction walk.
  Timeline timeline;             ///< Advances the single direction walk.
  Pipeline<W, H> pipeline;       ///< Scan output pipeline.
};
