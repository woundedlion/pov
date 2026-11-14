/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file contour_layout.h
 * @brief Persistent radii, trig and waveform tables of a concentric contour
 *        stack, and the per-frame inputs its renderers share.
 */

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "core/animation/orientation.h"
#include "core/color/baked_palette.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/platform/platform.h"
#include "core/render/clip.h"
#include "core/render/plot/cull.h"
#include "core/render/plot/shapes.h"

namespace Plot {

/** @brief Plot primitives a contour stack strokes. */
enum class ContourShape : uint8_t {
  PLANAR_POLYGON,
  SPHERICAL_POLYGON,
  FLOWER,
  PLANAR_STAR,
  SPHERICAL_STAR
};

/** @brief Most sides, petals or star points a contour takes. */
inline constexpr int CONTOUR_MAX_SIDES = 16;
/** @brief Most vertices a contour takes. */
inline constexpr int CONTOUR_MAX_VERTICES = 2 * CONTOUR_MAX_SIDES;
/** @brief Contour count from which star edges switch to screen-step-balanced
 *  sampling: the dense planar-star path, and the balanced policy for
 *  spherical stars. */
inline constexpr float DENSE_CONTOUR_COUNT = 32.0f;

/** @brief Folds a pole-to-pole radius fraction onto its nearer pole.
 * @param radius_t Radius fraction in [0, 1].
 * @return Distance fraction from the nearer pole, in [0, 0.5]. */
constexpr float folded_radius_t(float radius_t) {
  return radius_t <= 0.5f ? radius_t : 1.0f - radius_t;
}

/** @brief Maps contour quantiles to the screen-space sampling envelope.
 * @param radius_t Uniform radius quantile in [0, 1].
 * @return Pole-to-pole radius fraction with screen-balanced density. */
inline HS_FLASH_MEMBER float screen_balanced_radius_t(float radius_t) {
  constexpr float DENSITY_FLOOR = 0.5f;
  constexpr float BREAK_ANGLE = math::PI_F / 6.0f;
  constexpr float DENSITY_INTEGRAL = 1.1278247916f;
  constexpr float BREAK_QUANTILE =
      DENSITY_FLOOR * BREAK_ANGLE / DENSITY_INTEGRAL;
  const bool far_side = radius_t > 0.5f;
  const float u = 2.0f * folded_radius_t(radius_t);
  const float theta = u < BREAK_QUANTILE ? u * DENSITY_INTEGRAL / DENSITY_FLOOR
                                         : acosf(DENSITY_INTEGRAL * (1.0f - u));
  const float folded = theta / math::PI_F;
  return far_side ? 1.0f - folded : folded;
}

/** @brief Tests a planar star contour's cap against the clip.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param clip Drawing clip.
 * @param near_cap,far_cap Cap centers at the two poles.
 * @param radius Contour radius.
 * @param pad Angular motion padding.
 * @return Whether the contour can touch the clip. */
template <int W, int H>
HS_HOT_INLINE bool
planar_star_cap_visible(const ClipRegion &clip, const CapCenter &near_cap,
                        const CapCenter &far_cap, float radius, float pad) {
  constexpr float AA_PAD = 2.0f * math::PI_F / W;
  const bool far_side = radius > 1.0f;
  const float cap_radius = far_side ? 2.0f - radius : radius;
  const float half_angle = cap_radius * (math::PI_F / 2.0f) + AA_PAD + pad;
  const float t2 = fminf(half_angle, math::PI_F);
  return cap_may_touch_clip<H>(clip, far_side ? far_cap : near_cap, t2,
                               sinf(t2));
}

/** @brief A contour's ordered slot and geometric radius. */
struct ContourSlot {
  int i;          ///< Contour index.
  float radius_t; ///< Palette coordinate.
  float radius;   ///< Geometric radius.
};

/** @brief A contour's phase and shaded color. */
struct ContourShading {
  float direction;     ///< Signed phase direction.
  float contour_phase; ///< Waveform displacement.
  Color4 shaded_color; ///< Stroke color and opacity.
};

/** @brief One frame's inputs to a contour stack renderer. */
struct ContourFrame {
  const BakedPalette &palette;            ///< Palette sampled by radius.
  const math::Orientation<> &orientation; ///< Stack orientation.
  ContourShape shape;                     ///< Primitive every contour takes.
  int count;                              ///< Contours drawn.
  int sides;        ///< Polygon side, flower petal or star point count.
  float amplitude;  ///< Peak contour rotation, radians.
  bool opposite;    ///< Toggles the rotation sense past the equator.
  bool sine;        ///< Whether the waveform is the sine.
  float phase;      ///< Waveform phase at the frame's end, turns.
  float twist_step; ///< Waveform phase advance over the frame, turns.
  float opacity;    ///< Frame opacity applied to every stroke.
};

/**
 * @brief Radii, waveform and draw order of up to @p MAX concentric contours.
 * @tparam MAX Contour capacity.
 * @details Contour i sits at radius fraction (i + 0.5) / count, remapped by
 * the spacing law. Tables are bound once from an arena and refilled by
 * prepare() when the count or spacing changes.
 */
template <int MAX> class ContourLayout {
public:
  /** @brief Planar-star radius trig cached per contour. */
  using RadiusTrig = Star<PlanarProjection>::RadiusTrig;
  /** @brief Planar-star step trig for the prepared side count. */
  using StepTrig = Star<PlanarProjection>::StepTrig;

  /** @brief Arena bytes init_storage() takes.
   *  @return Byte count. */
  static constexpr size_t storage_bytes() {
    return MAX * (4 * sizeof(float) + sizeof(uint16_t) + sizeof(RadiusTrig));
  }

  /** @brief Binds the tables.
   *  @param arena Arena the tables are allocated from. */
  HS_FLASH_MEMBER void init_storage(Arena &arena) {
    spaced = arena.allocate_n<float>(MAX);
    phase_sine = arena.allocate_n<float>(MAX);
    phase_cosine = arena.allocate_n<float>(MAX);
    wave = arena.allocate_n<float>(MAX);
    draw_order = arena.allocate_n<uint16_t>(MAX);
    trig = arena.allocate_n<RadiusTrig>(MAX);
  }

  /** @brief Draw-order position to contour index, alternating outward from
   *  the middle contour.
   *  @param count Contour count.
   *  @param ordinal Draw-order position.
   *  @return Contour index. */
  static constexpr int folded_draw_index(int count, int ordinal) {
    if ((count & 1) != 0) {
      if (ordinal == 0)
        return count / 2;
      const int offset = (ordinal + 1) / 2;
      return (ordinal & 1) != 0 ? count / 2 + offset : count / 2 - offset;
    }
    const int offset = ordinal / 2;
    return (ordinal & 1) == 0 ? count / 2 + offset : (count - 1) / 2 - offset;
  }

  /** @brief Phase sign of a polygon, flower or planar-star contour.
   *  @param opposite Toggles the rotation sense past the equator.
   *  @param radius Contour radius in [0, 2].
   *  @return +1 or -1. */
  static float phase_sign_contour(bool opposite, float radius) {
    return !opposite && radius > 1.0f ? -1.0f : 1.0f;
  }

  /** @brief Phase sign of a spherical-star contour.
   *  @param opposite Toggles the rotation sense past the equator.
   *  @param radius Contour radius in [0, 2].
   *  @return +1 or -1. */
  static float phase_sign_star(bool opposite, float radius) {
    return opposite && radius > 1.0f ? -1.0f : 1.0f;
  }

  /** @brief Fills the radius, trig and draw-order tables.
   *  @param count Contour count, in [1, MAX].
   *  @param screen_balanced Applies screen_balanced_radius_t(). */
  HS_FLASH_MEMBER void prepare(int count, bool screen_balanced) {
    HS_CHECK(count >= 1 && count <= MAX,
             "ContourLayout: contour count %d outside table capacity", count);
    for (int i = 0; i < count; ++i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      phase_sine[i] = sinf(2.0f * math::PI_F * radius_t);
      phase_cosine[i] = cosf(2.0f * math::PI_F * radius_t);
      spaced[i] =
          screen_balanced ? screen_balanced_radius_t(radius_t) : radius_t;
      trig[i] = Star<PlanarProjection>::radius_trig(2.0f * spaced[i]);
      draw_order[i] = static_cast<uint16_t>(folded_draw_index(count, i));
    }
  }

  /** @brief Caches the planar-star step trig for @p sides.
   *  @param sides Star point count. */
  void prepare_star_sides(int sides) {
    if (sides == prepared_star_sides)
      return;
    star_step = Star<PlanarProjection>::step_trig(sides);
    prepared_star_sides = sides;
  }

  /** @brief Samples the waveform at every contour's radius fraction.
   *  @tparam Wave Callable mapping a phase in turns to [-1, 1].
   *  @param phase Waveform phase, turns.
   *  @param count Contour count.
   *  @param sine Uses the cached sine table instead of @p evaluate.
   *  @param evaluate Waveform sampled at radius fraction plus @p phase. */
  template <typename Wave>
  HS_FLASH_MEMBER void prepare_waveform(float phase, int count, bool sine,
                                        const Wave &evaluate) {
    if (sine) {
      const float phase_angle = 2.0f * math::PI_F * phase;
      const float phase_sin = sinf(phase_angle);
      const float phase_cos = cosf(phase_angle);
      for (int i = 0; i < count; ++i)
        wave[i] = phase_sine[i] * phase_cos + phase_cosine[i] * phase_sin;
      return;
    }

    for (int i = 0; i < count; ++i) {
      const float radius_t =
          (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
      wave[i] = evaluate(radius_t + phase);
    }
  }

  /** @brief Selects a contour and refreshes its paired palette color.
   * @param ordinal Draw-order position.
   * @param count Number of contours.
   * @param continuous_star Whether contours have independent colors.
   * @param palette Palette sampled by radius.
   * @param pair_color Current paired palette color.
   * @return Contour index and radii. */
  HS_HOT_INLINE ContourSlot slot(int ordinal, int count, bool continuous_star,
                                 const BakedPalette &palette,
                                 Color4 &pair_color) const {
    const int i = continuous_star ? count - ordinal - 1 : draw_order[ordinal];
    const float radius_t =
        (static_cast<float>(i) + 0.5f) / static_cast<float>(count);
    const bool starts_pair = (count & 1) != 0
                                 ? ordinal == 0 || (ordinal & 1) != 0
                                 : (ordinal & 1) == 0;
    if (continuous_star || starts_pair)
      pair_color = palette.get(radius_t);
    return {i, radius_t, 2.0f * spaced[i]};
  }

  /** @brief Evaluates a contour's phase and opacity.
   * @param i Contour index.
   * @param radius Contour radius.
   * @param continuous_star Whether phase follows the spherical-star rule.
   * @param opposite Toggles the rotation sense past the equator.
   * @param amplitude Peak contour rotation, radians.
   * @param pair_color Selected palette color.
   * @param global_alpha Frame opacity.
   * @return Phase direction, displacement and shaded color. */
  HS_HOT_INLINE ContourShading shading(int i, float radius,
                                       bool continuous_star, bool opposite,
                                       float amplitude,
                                       const Color4 &pair_color,
                                       float global_alpha) const {
    const float direction = continuous_star
                                ? phase_sign_star(opposite, radius)
                                : phase_sign_contour(opposite, radius);
    const float contour_phase = direction * amplitude * wave[i];
    Color4 shaded_color = pair_color;
    shaded_color.alpha *= global_alpha;
    return {direction, contour_phase, shaded_color};
  }

  /** @brief Spaced radius fraction of contour @p i.
   *  @param i Contour index.
   *  @return Radius fraction in [0, 1]. */
  float spaced_radius_t(int i) const { return spaced[i]; }

  /** @brief Contour drawn at draw-order position @p ordinal.
   *  @param ordinal Draw-order position.
   *  @return Contour index. */
  int draw_index(int ordinal) const { return draw_order[ordinal]; }

  /** @brief Planar-star radius trig of contour @p i.
   *  @param i Contour index.
   *  @return Cached trig. */
  const RadiusTrig &radius_trig(int i) const { return trig[i]; }

  /** @brief Planar-star step trig of the prepared side count.
   *  @return Cached trig. */
  const StepTrig &star_step_trig() const { return star_step; }

private:
  float *spaced = nullptr;
  float *phase_sine = nullptr;
  float *phase_cosine = nullptr;
  float *wave = nullptr;
  uint16_t *draw_order = nullptr;
  RadiusTrig *trig = nullptr;
  StepTrig star_step{};
  int prepared_star_sides = 0;
};

} // namespace Plot
