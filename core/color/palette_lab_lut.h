/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/baked_palette.h"
#include "core/color/color_space.h"
#include "core/color/hue_rotation_lut.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>

/** @file palette_lab_lut.h
 * @brief Per-interval RGB-first palette blends cached in OKLab. */

/** @brief Borrowed neighboring time rows for opaque palettes. */
struct PaletteLabTimeView {
  const OKLab *low = nullptr;             ///< Lower time row.
  const OKLab *high = nullptr;            ///< Upper time row.
  float weight = 0.0f;                    ///< Upper time weight.
  int value_steps = 64;                   ///< Nodes in each borrowed row.
  const math::Vector *rgb_low = nullptr;  ///< Optional lower hue table.
  const math::Vector *rgb_high = nullptr; ///< Optional upper hue table.
};

/** @brief Borrowed palette-coordinate by frame-time base-color cache. */
struct PaletteLabLutBase {
  static constexpr int TIMES = 3; ///< Frame-time nodes.
  const OKLab *data = nullptr;    ///< Borrowed RGB-first blended colors.
  int value_steps = 64;           ///< Nodes in each time row.
  bool active = false;            ///< Whether both palettes are opaque.
  const math::Vector *rgb_nodes = nullptr; ///< Optional final hue table.
  bool rgb_midpoint = false; ///< Shares the interval midpoint for every column.
  uint8_t rgb_midpoint_row = 1; ///< Physical RGB row containing the midpoint.

  /** @brief Prepares the two time rows for one column.
   * @param fraction Frame progress in [0, 1].
   * @return Borrowed rows, or an inactive view. */
  HS_FLASH_INLINE PaletteLabTimeView view(float fraction) const {
    if (!active)
      return {};
    const float POSITION = hs::clamp(fraction, 0.0f, 1.0f) * (TIMES - 1);
    const int LOW_INDEX = std::min(int(POSITION), TIMES - 2);
    return {data + LOW_INDEX * value_steps,
            data + (LOW_INDEX + 1) * value_steps,
            POSITION - LOW_INDEX,
            value_steps,
            rgb_nodes
                ? rgb_nodes + (rgb_midpoint ? rgb_midpoint_row : LOW_INDEX) *
                                  HueRotationLutView::SIZE
                : nullptr,
            rgb_nodes ? rgb_nodes +
                            (rgb_midpoint ? rgb_midpoint_row : LOW_INDEX + 1) *
                                HueRotationLutView::SIZE
                      : nullptr};
  }

  /** @brief Samples neighboring palette nodes from both prepared time rows.
   * @param time Prepared column rows.
   * @param value Palette coordinate in [0, 1].
   * @return Bilinearly interpolated base OKLab. */
  __attribute__((always_inline)) static OKLab
  sample(const PaletteLabTimeView &time, float value) {
    assert(time.low && time.high && time.value_steps >= 2);
    const float POSITION =
        hs::clamp(value, 0.0f, 1.0f) * (time.value_steps - 1);
    const int LOW_INDEX = std::min(int(POSITION), time.value_steps - 2);
    const float WEIGHT = POSITION - LOW_INDEX;
    const auto lerp = [](const OKLab &a, const OKLab &b, float fraction) {
      return OKLab{hs::lerp(a.L, b.L, fraction), hs::lerp(a.a, b.a, fraction),
                   hs::lerp(a.b, b.b, fraction)};
    };
    return lerp(lerp(time.low[LOW_INDEX], time.low[LOW_INDEX + 1], WEIGHT),
                lerp(time.high[LOW_INDEX], time.high[LOW_INDEX + 1], WEIGHT),
                time.weight);
  }
};

/** @brief Owned interval palette bases.
 * @tparam ValueSteps Palette-coordinate nodes per time row. */
template <int ValueSteps> struct PaletteLabLutStorage : PaletteLabLutBase {
  static_assert(ValueSteps >= 2);
  static constexpr int VALUES = ValueSteps;  ///< Palette-coordinate nodes.
  std::array<OKLab, VALUES * TIMES> nodes{}; ///< RGB-first blended colors.

  /** @brief Binds the borrowed view to this cache's storage. */
  PaletteLabLutStorage() {
    data = nodes.data();
    value_steps = VALUES;
  }
  PaletteLabLutStorage(const PaletteLabLutStorage &other)
      : PaletteLabLutBase(other), nodes(other.nodes) {
    data = nodes.data();
  }
  PaletteLabLutStorage &operator=(const PaletteLabLutStorage &other) {
    if (this != &other) {
      static_cast<PaletteLabLutBase &>(*this) = other;
      nodes = other.nodes;
    }
    data = nodes.data();
    return *this;
  }

  /** @brief Bakes RGB16 endpoint blends before converting each node to OKLab.
   * @param from Previous immutable palette.
   * @param to Current immutable palette. */
  HS_COLD_MEMBER void bake(const BakedPalette &from, const BakedPalette &to) {
    rgb_nodes = nullptr;
    rgb_midpoint = false;
    rgb_midpoint_row = 1;
    const auto A = from.constant_alpha();
    const auto B = to.constant_alpha();
    active = A.has_value() && B.has_value() && *A == 1.0f && *B == 1.0f;
    if (!active)
      return;
    const auto convert = [](const Pixel &rgb) HS_HOT_INLINE {
      const LinRGB LINEAR = pixel_to_linrgb(rgb);
      return linear_rgb_to_oklab_fast(LINEAR.r, LINEAR.g, LINEAR.b);
    };
    for (int value = 0; value < VALUES; ++value) {
      const float COORDINATE = float(value) / (VALUES - 1);
      const Pixel FROM = from.get_color(COORDINATE);
      const Pixel TO = to.get_color(COORDINATE);
      const OKLab FROM_LAB = convert(FROM);
      const OKLab TO_LAB = FROM == TO ? FROM_LAB : convert(TO);
      const Pixel MIDPOINT = FROM.lerp16(TO, frac_to_q16(0.5f));
      nodes[value] = FROM_LAB;
      nodes[VALUES + value] = MIDPOINT == FROM ? FROM_LAB
                              : MIDPOINT == TO ? TO_LAB
                                               : convert(MIDPOINT);
      nodes[2 * VALUES + value] = TO_LAB;
    }
  }
};

/// Default interval palette cache resolution.
using PaletteLabLut = PaletteLabLutStorage<64>;

/** @brief Interval palette bases baked through the final hue and gamut mapping.
 * @tparam MidpointOnly Whether storage is limited to one admitted midpoint row. */
template <bool MidpointOnly> struct PaletteRgbHueLutImpl {
  static constexpr int RGB_TIMES =
      MidpointOnly ? 1 : PaletteLabLut::TIMES; ///< Stored time rows.
  static constexpr int MIDPOINT_ROW =
      MidpointOnly ? 0 : 1; ///< Stored midpoint row.
  static constexpr float MAX_INTERVAL_DELTA =
      .005f; ///< Qualified endpoint Lab step.
  static constexpr float MAX_REUSE_DELTA =
      .001f; ///< Qualified baked-base Lab step.
  std::array<OKLab, PaletteLabLut::VALUES>
      baked_midpoint{};         ///< Original baked midpoint key.
  bool midpoint_active = false; ///< Whether the middle RGB row has a baked key.
  static constexpr bool APPROXIMATE =
      true; ///< Samples a trilinear color table.
  std::array<math::Vector, RGB_TIMES * HueRotationLutView::SIZE>
      nodes{}; ///< Quantized RGB16 code units by time, value and hue.
  /** @brief Bakes final hue colors and reuses an exactly matching retained row.
   * @param base Interval base receiving this cache binding.
   * @param previous_base Prior interval's immutable base, when retained.
   * @param previous Prior interval's RGB owner, when retained. */
  HS_COLD_MEMBER void bake(PaletteLabLut &base,
                           const PaletteLabLut *previous_base = nullptr,
                           const PaletteRgbHueLutImpl *previous = nullptr) {
    base.rgb_nodes = nullptr;
    midpoint_active = false;
    if (!base.active)
      return;
    float difference = 0.0f;
    for (int value = 0; value < PaletteLabLut::VALUES; ++value) {
      const auto &a = base.nodes[value];
      const auto &b =
          base.nodes[(PaletteLabLut::TIMES - 1) * PaletteLabLut::VALUES +
                     value];
      difference = std::max(
          difference,
          std::max({fabsf(a.L - b.L), fabsf(a.a - b.a), fabsf(a.b - b.b)}));
    }
    base.rgb_midpoint = difference <= MAX_INTERVAL_DELTA;
    base.rgb_midpoint_row = MIDPOINT_ROW;
    if constexpr (MidpointOnly)
      if (!base.rgb_midpoint)
        return;
    if (base.rgb_midpoint && previous_base && previous && previous != this &&
        previous->midpoint_active &&
        previous_base->rgb_nodes == previous->nodes.data()) {
      float delta = 0.0f;
      for (int value = 0; value < PaletteLabLut::VALUES; ++value) {
        const auto &a = base.nodes[PaletteLabLut::VALUES + value];
        const auto &b = previous->baked_midpoint[value];
        delta = std::max(delta, std::max({fabsf(a.L - b.L), fabsf(a.a - b.a),
                                          fabsf(a.b - b.b)}));
      }
      if (delta <= MAX_REUSE_DELTA) {
        std::copy_n(previous->nodes.data() +
                        MIDPOINT_ROW * HueRotationLutView::SIZE,
                    HueRotationLutView::SIZE,
                    nodes.data() + MIDPOINT_ROW * HueRotationLutView::SIZE);
        baked_midpoint = previous->baked_midpoint;
        midpoint_active = true;
        base.rgb_nodes = nodes.data();
        return;
      }
    }
    std::array<float, HueRotationLutView::HUE_STEPS> cosines, sines;
    for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue)
      hue_sincos(float(hue) / HueRotationLutView::HUE_STEPS, cosines[hue],
                 sines[hue]);
    if constexpr (MidpointOnly) {
      for (int value = 0; value < PaletteLabLut::VALUES; ++value)
        for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue)
          nodes[value * HueRotationLutView::HUE_STEPS + hue] =
              rotate_node(base.nodes[PaletteLabLut::VALUES + value],
                          cosines[hue], sines[hue]);
    } else {
      const bool REUSE =
          !base.rgb_midpoint && previous_base && previous &&
          previous_base->active && !previous_base->rgb_midpoint &&
          previous_base->rgb_nodes == previous->nodes.data() &&
          std::memcmp(base.nodes.data(),
                      previous_base->nodes.data() +
                          (PaletteLabLut::TIMES - 1) * PaletteLabLut::VALUES,
                      PaletteLabLut::VALUES * sizeof(OKLab)) == 0;
      if (REUSE)
        std::copy_n(previous->nodes.data() +
                        (PaletteLabLut::TIMES - 1) * HueRotationLutView::SIZE,
                    HueRotationLutView::SIZE, nodes.data());
      for (int time = base.rgb_midpoint || REUSE ? 1 : 0;
           time < (base.rgb_midpoint ? 2 : PaletteLabLut::TIMES); ++time)
        for (int value = 0; value < PaletteLabLut::VALUES; ++value)
          for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue)
            nodes[time * HueRotationLutView::SIZE +
                  value * HueRotationLutView::HUE_STEPS + hue] =
                rotate_node(base.nodes[time * PaletteLabLut::VALUES + value],
                            cosines[hue], sines[hue]);
    }
    midpoint_active = base.rgb_midpoint;
    if (midpoint_active)
      std::copy_n(base.nodes.data() + PaletteLabLut::VALUES,
                  PaletteLabLut::VALUES, baked_midpoint.data());
    base.rgb_nodes = nodes.data();
  }
  /** @brief Final hue node with a prepared rotation.
   * @param lab Baked interval base.
   * @param cosine Hue cosine.
   * @param sine Hue sine.
   * @return Gamut-mapped RGB16 node. */
  HS_HOT_INLINE static math::Vector rotate_node(const OKLab &lab, float cosine,
                                                float sine) {
    const Pixel RGB_NODE = linrgb_to_pixel(oklab_to_linear_rgb_lut_gamut(
        {lab.L, lab.a * cosine - lab.b * sine, lab.a * sine + lab.b * cosine}));
    return math::Vector(float(RGB_NODE.r), float(RGB_NODE.g),
                        float(RGB_NODE.b));
  }
  /** @brief Trilinear palette, hue and time lookup.
   * @param time Prepared neighboring RGB rows.
   * @param value Palette coordinate in [0, 1].
   * @param hue Periodic hue turns.
   * @return Interpolated RGB16 color. */
  HS_HOT_INLINE static Pixel sample(const PaletteLabTimeView &time, float value,
                                    float hue) {
    const HueRotationValueSample coordinate(value);
    const HueRotationHueSample rotation(hue);
    const math::Vector LOW_RGB = sample_hue_rotation_linear(
        {nullptr, true, time.rgb_low}, coordinate, rotation);
    math::Vector result = LOW_RGB;
    if (time.rgb_low != time.rgb_high) {
      const math::Vector HIGH_RGB = sample_hue_rotation_linear(
          {nullptr, true, time.rgb_high}, coordinate, rotation);
      const float WEIGHT = frac_to_q16(time.weight) * Q16_TO_UNIT;
      result = LOW_RGB + WEIGHT * (HIGH_RGB - LOW_RGB);
    }
    return Pixel(round_linear_channel(result.x), round_linear_channel(result.y),
                 round_linear_channel(result.z));
  }
};

/** @brief Three-time-row interval hue cache. */
using PaletteRgbHueLut = PaletteRgbHueLutImpl<false>;
