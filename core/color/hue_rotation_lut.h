/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/color/color_space.h"
#include <array>
#include <span>

/** @file hue_rotation_lut.h
 * @brief Palette-coordinate by hue-rotation acceleration. */

/** @brief Largest float below 1: the top of the half-open palette coordinate
 *  the hue-rotation LUT is baked over and sampled at. */
inline constexpr float UNIT_OPEN_MAX = 0x1.fffffep-1f;

struct HueRotationGamutMask;

/** @brief View over palette-coordinate by hue-rotation colors. */
struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64; ///< Palette-coordinate samples.
  static constexpr int HUE_STEPS = 16;   ///< Hue-rotation samples per turn.
  static_assert(HUE_STEPS > 0 && (HUE_STEPS & (HUE_STEPS - 1)) == 0,
                "hue mask wrap requires a power of two");
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS; ///< Entries.

  /// Borrowed table, indexed `value_index * HUE_STEPS + hue_index`.
  const Pixel *data;
  bool active; ///< Whether `data` holds a baked table to sample.
  const Pixel *linear_data =
      nullptr; ///< Optional RGB16 nodes for linear interpolation.
  const HueRotationGamutMask *gamut_mask =
      nullptr; ///< Cells needing direct hue.
};

/** @brief Cells crossing gamut classes or the qualified RGB corner span. */
struct HueRotationGamutMask {
  static constexpr int MAX_CHANNEL_SPAN =
      14000; ///< Qualified RGB16 corner span.
  /** @brief Number of palette-coordinate by hue cells. */
  static constexpr int CELLS =
      (HueRotationLutView::VALUE_STEPS - 1) * HueRotationLutView::HUE_STEPS;
  std::array<uint8_t, (CELLS + 7) / 8> bits{}; ///< Value-major cell flags.

  /** @brief Clears all direct-sampling flags. */
  void clear() { bits.fill(0); }
  /** @brief Marks a cell for direct hue evaluation.
   * @param value Lower palette-coordinate node index.
   * @param hue Lower hue node index. */
  void mark_cell(int value, int hue) {
    const int CELL = value * HueRotationLutView::HUE_STEPS + hue;
    bits[CELL / 8] |= uint8_t(1u << (CELL % 8));
  }
  /** @brief Tests the cell containing a pair of lower node indices.
   * @param value Lower palette-coordinate node index.
   * @param hue Lower hue node index.
   * @return Whether direct hue evaluation is required. */
  bool requires_direct(int value, int hue) const {
    const int CELL = std::min(value, HueRotationLutView::VALUE_STEPS - 2) *
                         HueRotationLutView::HUE_STEPS +
                     hue;
    return (bits[CELL / 8] & uint8_t(1u << (CELL % 8))) != 0;
  }
};

/** @brief Accumulates differing corner classes in value-major node order. */
struct HueRotationGamutBuilder {
  HueRotationGamutMask *mask; ///< Destination cell flags.
  std::array<uint8_t, HueRotationLutView::HUE_STEPS>
      previous{}, ///< Previous row raw gamut classes.
      current{};  ///< Current row raw gamut classes.
  std::array<Pixel, HueRotationLutView::HUE_STEPS>
      previous_rgb{}, ///< Previous row mapped RGB16 nodes.
      current_rgb{};  ///< Current row mapped RGB16 nodes.

  /** @brief Binds the destination mask.
   * @param destination Cell flags populated by appended rows. */
  explicit HueRotationGamutBuilder(HueRotationGamutMask *destination)
      : mask(destination) {}

  /** @brief Appends a node in value-major order.
   * @param value Palette-coordinate node index.
   * @param hue Hue node index.
   * @param rgb Unmapped linear RGB.
   * @param node Gamut-mapped RGB16. */
  void append(int value, int hue, const LinRGB &rgb, const Pixel &node) {
    current_rgb[hue] = node;
    current[hue] = uint8_t((rgb.r < 0.0f ? 1 : 0) | (rgb.r > 1.0f ? 2 : 0) |
                           (rgb.g < 0.0f ? 4 : 0) | (rgb.g > 1.0f ? 8 : 0) |
                           (rgb.b < 0.0f ? 16 : 0) | (rgb.b > 1.0f ? 32 : 0));
    if (hue != HueRotationLutView::HUE_STEPS - 1)
      return;
    if (value != 0)
      for (int column = 0; column < HueRotationLutView::HUE_STEPS; ++column) {
        const int NEXT = (column + 1) & (HueRotationLutView::HUE_STEPS - 1);
        const uint8_t FIRST = previous[column];
        const auto span = [&](auto channel) {
          const int A = channel(previous_rgb[column]);
          const int B = channel(previous_rgb[NEXT]);
          const int C = channel(current_rgb[column]);
          const int D = channel(current_rgb[NEXT]);
          return std::max({A, B, C, D}) - std::min({A, B, C, D});
        };
        if (FIRST != previous[NEXT] || FIRST != current[column] ||
            FIRST != current[NEXT] ||
            (FIRST != 0 && (span([](Pixel p) { return p.r; }) >
                                HueRotationGamutMask::MAX_CHANNEL_SPAN ||
                            span([](Pixel p) { return p.g; }) >
                                HueRotationGamutMask::MAX_CHANNEL_SPAN ||
                            span([](Pixel p) { return p.b; }) >
                                HueRotationGamutMask::MAX_CHANNEL_SPAN)))
          mask->mark_cell(value - 1, column);
      }
    previous = current;
    previous_rgb = current_rgb;
  }
};

/** @brief Conditional retained hue-mask storage.
 * @tparam Enabled Whether retained hue sampling is enabled. */
template <bool Enabled> struct HueRotationGamutStorage {};
/** @brief Retained hue mask for enabled hue sampling. */
template <> struct HueRotationGamutStorage<true> {
  HueRotationGamutMask retained_hue_rotation_gamut; ///< Prior frame cell flags.
};

/** @brief Appends one raw hue node to a cell classifier.
 * @param builder Destination row classifier.
 * @param value Palette-coordinate node index.
 * @param hue Hue node index.
 * @param base Unrotated perceptual color.
 * @param cosine Prepared hue cosine.
 * @param sine Prepared hue sine.
 * @param node Gamut-mapped RGB16 node. */
inline void mark_hue_gamut_node(HueRotationGamutBuilder &builder, int value,
                                int hue, const OKLab &base, float cosine,
                                float sine, const Pixel &node) {
  const auto RAW = oklab_to_linear_rgb(rotate_chroma(base, cosine, sine));
  builder.append(value, hue, RAW, node);
}

/**
 * @brief Bakes all palette-coordinate and hue-rotation combinations.
 * @tparam Source Palette source exposing Color4 get(float) const.
 * @param output Destination LUT.
 * @param palette Palette source to hue-rotate.
 * @param gamut_mask Optional cell mask owned alongside the destination table.
 */
template <typename Source>
HS_FLASH_INLINE __attribute__((aligned(32))) inline void
prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE> output,
                         const Source &palette,
                         HueRotationGamutMask *gamut_mask = nullptr) {
  if (gamut_mask)
    gamut_mask->clear();
  HueRotationGamutBuilder builder(gamut_mask);
  for (int value_index = 0; value_index < HueRotationLutView::VALUE_STEPS;
       ++value_index) {
    const float value =
        UNIT_OPEN_MAX * value_index / (HueRotationLutView::VALUE_STEPS - 1);
    const HueRotateBase base = make_hue_rotate_base(palette.get(value));
    for (int hue_index = 0; hue_index < HueRotationLutView::HUE_STEPS;
         ++hue_index) {
      const float amount =
          static_cast<float>(hue_index) / HueRotationLutView::HUE_STEPS;
      const Pixel NODE = hue_rotate_lut_gamut(base, amount).color;
      if (gamut_mask) {
        float cosine, sine;
        hue_sincos(amount, cosine, sine);
        mark_hue_gamut_node(builder, value_index, hue_index, base.lab, cosine,
                            sine, NODE);
      }
      output[value_index * HueRotationLutView::HUE_STEPS + hue_index] = NODE;
    }
  }
}

/** @brief Palette-axis coordinates shared by hue-rotation LUT samples. */
struct HueRotationValueSample {
  int low;         ///< Lower palette row.
  int high;        ///< Upper palette row.
  uint16_t weight; ///< Upper row weight in Q16.

  /** @brief Prepares a palette-axis sample.
   * @param value Palette coordinate in [0, 1]. */
  explicit HueRotationValueSample(float value) {
    const float position =
        hs::clamp(value * (HueRotationLutView::VALUE_STEPS - 1), 0.0f,
                  static_cast<float>(HueRotationLutView::VALUE_STEPS - 1));
    low = lut_index_lo(position);
    high = std::min(low + 1, HueRotationLutView::VALUE_STEPS - 1);
    weight = lut_index_weight(position, low);
  }
};

/** @brief Hue-axis coordinates shared by palette endpoint samples. */
struct HueRotationHueSample {
  int low;             ///< Lower hue column, wrapped into the table.
  int high;            ///< Next hue column, wrapped into the table.
  uint16_t weight;     ///< Upper hue column weight in Q16.
  float linear_weight; ///< Q16 weight converted to a float fraction.

  /** @brief Prepares a periodic hue-axis sample.
   * @param amount Hue rotation in turns. */
  explicit HueRotationHueSample(float amount) {
    const float POSITION =
        hs::clamp((amount - floorf(amount)) * HueRotationLutView::HUE_STEPS,
                  0.0f, static_cast<float>(HueRotationLutView::HUE_STEPS));
    const int STEP = lut_index_lo(POSITION);
    low = STEP & (HueRotationLutView::HUE_STEPS - 1);
    high = (low + 1) & (HueRotationLutView::HUE_STEPS - 1);
    weight = lut_index_weight(POSITION, STEP);
    linear_weight = weight * Q16_TO_UNIT;
  }
};

/** @brief Samples hue rotation at prepared palette-axis coordinates.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param hue Prepared hue-axis coordinates.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view,
                        const HueRotationValueSample &value,
                        const HueRotationHueSample &hue) {
  const auto sample_row = [&](int row) {
    const int offset = row * HueRotationLutView::HUE_STEPS;
    return view.data[offset + hue.low].lerp16(view.data[offset + hue.high],
                                              hue.weight);
  };
  return sample_row(value.low).lerp16(sample_row(value.high), value.weight);
}

/** @brief Samples hue rotation at a prepared palette coordinate.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param amount Hue rotation in turns.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view,
                        const HueRotationValueSample &value, float amount) {
  return sample_hue_rotation_lut(view, value, HueRotationHueSample(amount));
}

/** @brief Samples a hue-rotation LUT with interpolation over both axes.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette coordinate in [0, 1].
 * @param amount Hue rotation in turns.
 * @return Hue-rotated linear RGB color. */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view, float value,
                        float amount) {
  return sample_hue_rotation_lut(view, HueRotationValueSample(value), amount);
}

/** @brief Samples hue rotation without rounding intermediate RGB blends.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param hue Prepared hue-axis coordinates.
 * @return Linear RGB channel values before quantization. */
__attribute__((always_inline)) inline math::Vector
sample_hue_rotation_linear(const HueRotationLutView &view,
                           const HueRotationValueSample &value,
                           const HueRotationHueSample &hue) {
  const auto row = [&](int index) {
    if (view.linear_data) {
      const Pixel &a =
          view.linear_data[index * HueRotationLutView::HUE_STEPS + hue.low];
      const Pixel &b =
          view.linear_data[index * HueRotationLutView::HUE_STEPS + hue.high];
      return math::Vector(hs::lerp(float(a.r), float(b.r), hue.linear_weight),
                          hs::lerp(float(a.g), float(b.g), hue.linear_weight),
                          hs::lerp(float(a.b), float(b.b), hue.linear_weight));
    }
    const Pixel &a = view.data[index * HueRotationLutView::HUE_STEPS + hue.low];
    const Pixel &b =
        view.data[index * HueRotationLutView::HUE_STEPS + hue.high];
    return math::Vector(hs::lerp(float(a.r), float(b.r), hue.linear_weight),
                        hs::lerp(float(a.g), float(b.g), hue.linear_weight),
                        hs::lerp(float(a.b), float(b.b), hue.linear_weight));
  };
  const math::Vector a = row(value.low);
  if (value.weight == 0)
    return a;
  const math::Vector b = row(value.high);
  return a + (value.weight * Q16_TO_UNIT) * (b - a);
}

/** @brief Samples hue rotation without intermediate RGB quantization.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette-axis coordinates.
 * @param amount Hue rotation in turns.
 * @return Linear RGB channel values before quantization. */
__attribute__((always_inline)) inline math::Vector
sample_hue_rotation_linear(const HueRotationLutView &view,
                           const HueRotationValueSample &value, float amount) {
  return sample_hue_rotation_linear(view, value, HueRotationHueSample(amount));
}
