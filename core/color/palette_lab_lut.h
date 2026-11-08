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
#include <cstring>

/** @file palette_lab_lut.h
 * @brief Per-interval RGB-first palette blends cached in OKLab. */

struct PaletteLabLutBase;

/** @brief Borrowed neighboring time rows for opaque palettes. */
struct PaletteLabTimeView {
  const OKLab *low = nullptr;      ///< Lower time row.
  const OKLab *high = nullptr;     ///< Upper time row.
  float weight = 0.0f;             ///< Upper time weight.
  int value_steps = 64;            ///< Nodes in each borrowed row.
  const Pixel *rgb_low = nullptr;  ///< Optional lower hue table.
  const Pixel *rgb_high = nullptr; ///< Optional upper hue table.
  const HueRotationGamutMask *rgb_mask_low = nullptr; ///< Lower RGB row flags.
  const HueRotationGamutMask *rgb_mask_high =
      nullptr; ///< Distinct upper row flags.
  const PaletteLabLutBase *canonical_base =
      nullptr;           ///< Borrowed endpoint owner.
  float fraction = 0.0f; ///< Frame progress within the interval.
};

/** @brief Borrowed palette-coordinate by frame-time base-color cache. */
struct PaletteLabLutBase {
  static constexpr int TIMES = 3;                    ///< Frame-time nodes.
  static constexpr float GAMUT_GUARD_WIDTH = 0.001f; ///< Linear-RGB gate width.
  const OKLab *data = nullptr;      ///< Borrowed RGB-first blended colors.
  int value_steps = 64;             ///< Nodes in each time row.
  bool active = false;              ///< Whether both palettes are opaque.
  const Pixel *rgb_nodes = nullptr; ///< Optional final hue table.
  const HueRotationGamutMask *rgb_masks = nullptr; ///< RGB row cell flags.
  bool rgb_midpoint = false; ///< Shares the interval midpoint for every column.
  uint8_t rgb_midpoint_row = 1; ///< Physical RGB row containing the midpoint.
  const BakedPalette *palette_from = nullptr; ///< Immutable interval start.
  const BakedPalette *palette_to = nullptr;   ///< Immutable interval end.

  /** @brief Whether raw RGB is near a discontinuous gamut gate.
   * @param raw Unmapped linear RGB.
   * @return Whether canonical-base restoration is required. */
  HS_HOT_INLINE static bool near_gamut_gate(const LinRGB &raw) {
    const float LEAST = __builtin_fminf(__builtin_fminf(raw.r, raw.g), raw.b);
    const float MOST = __builtin_fmaxf(__builtin_fmaxf(raw.r, raw.g), raw.b);
    return !(fabsf(LEAST - LINEAR_RGB_GAMUT_MIN) > GAMUT_GUARD_WIDTH &&
             fabsf(MOST - LINEAR_RGB_GAMUT_MAX) > GAMUT_GUARD_WIDTH);
  }

  /** @brief Finishes the immutable RGB-first base through a prepared hue.
   * @param time Borrowed endpoint palettes and interval progress.
   * @param value Palette coordinate.
   * @param cosine Prepared hue cosine.
   * @param sine Prepared hue sine.
   * @return Gamut-mapped linear RGB. */
  HS_HOT_FLASH_MEMBER __attribute__((flatten)) static LinRGB
  canonical_hue_rgb(const PaletteLabTimeView &time, float value, float cosine,
                    float sine) {
    const auto RGB = time.canonical_base->palette_from->get_color(value).lerp16(
        time.canonical_base->palette_to->get_color(value),
        frac_to_q16(time.fraction));
    const auto LINEAR_INPUT = pixel_to_linrgb(RGB);
    const auto LAB =
        rotate_chroma(linear_rgb_to_oklab_fast(LINEAR_INPUT.r, LINEAR_INPUT.g,
                                               LINEAR_INPUT.b),
                      cosine, sine);
    return oklab_to_linear_rgb_lut_gamut(LAB);
  }

  /** @brief Finishes the immutable RGB-first base through a prepared hue.
   * @param time Borrowed endpoint palettes and interval progress.
   * @param value Palette coordinate.
   * @param cosine Prepared hue cosine.
   * @param sine Prepared hue sine.
   * @return Gamut-mapped RGB16. */
  HS_HOT_INLINE static Pixel canonical_hue_color(const PaletteLabTimeView &time,
                                                 float value, float cosine,
                                                 float sine) {
    return linrgb_to_pixel(canonical_hue_rgb(time, value, cosine, sine));
  }

  /** @brief Prepares the two time rows for one column.
   * @param fraction Frame progress in [0, 1].
   * @return Borrowed rows, or an inactive view.
   * @pre Cache, rows, masks and endpoint palettes remain alive and unchanged
   * while the returned view is sampled. */
  HS_FLASH_INLINE PaletteLabTimeView view(float fraction) const {
    if (!active)
      return {};
    const float POSITION = hs::clamp(fraction, 0.0f, 1.0f) * (TIMES - 1);
    const int LOW_INDEX = std::min(int(POSITION), TIMES - 2);
    const auto *const MASK_LOW =
        rgb_masks ? rgb_masks + (rgb_midpoint ? rgb_midpoint_row : LOW_INDEX)
                  : nullptr;
    const auto *const MASK_HIGH =
        rgb_masks
            ? rgb_masks + (rgb_midpoint ? rgb_midpoint_row : LOW_INDEX + 1)
            : nullptr;
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
                      : nullptr,
            MASK_LOW,
            MASK_HIGH == MASK_LOW ? nullptr : MASK_HIGH,
            this,
            hs::clamp(fraction, 0.0f, 1.0f)};
  }

  /** @brief Samples neighboring palette nodes from both prepared time rows.
   * @param time Prepared column rows.
   * @param value Palette coordinate in [0, 1].
   * @return Bilinearly interpolated base OKLab. */
  __attribute__((always_inline)) static OKLab
  sample(const PaletteLabTimeView &time, float value) {
    HS_AUDIT_CHECK(time.low && time.high && time.value_steps >= 2);
    const float POSITION =
        hs::clamp(value, 0.0f, 1.0f) * (time.value_steps - 1);
    const int LOW_INDEX = std::min(int(POSITION), time.value_steps - 2);
    const float WEIGHT = POSITION - LOW_INDEX;
    const auto lerp = [](const OKLab &a, const OKLab &b,
                         float fraction) HS_HOT_INLINE {
      return OKLab{hs::lerp(a.L, b.L, fraction), hs::lerp(a.a, b.a, fraction),
                   hs::lerp(a.b, b.b, fraction)};
    };
    return lerp(lerp(time.low[LOW_INDEX], time.low[LOW_INDEX + 1], WEIGHT),
                lerp(time.high[LOW_INDEX], time.high[LOW_INDEX + 1], WEIGHT),
                time.weight);
  }

  /** @brief Finishes a marked hue cell from its perceptual base.
   * @param time Borrowed neighboring base rows.
   * @param value Palette coordinate.
   * @param hue Rotation in turns.
   * @return Gamut-mapped linear RGB. */
  HS_HOT_FLASH_MEMBER static LinRGB
  direct_hue_rgb(const PaletteLabTimeView &time, float value, float hue) {
    float cosine, sine;
    hue_sincos(hue, cosine, sine);
    OKLab lab =
        rotate_chroma(PaletteLabLutBase::sample(time, value), cosine, sine);
    LinRGB raw = oklab_to_linear_rgb(lab);
    if (PaletteLabLutBase::near_gamut_gate(raw) && time.canonical_base &&
        time.canonical_base->palette_from && time.canonical_base->palette_to) {
      return PaletteLabLutBase::canonical_hue_rgb(time, value, cosine, sine);
    }
    if (!linear_rgb_in_gamut(raw.r, raw.g, raw.b))
      raw = oklab_to_linear_rgb(gamut_scale_to_boundary_lut(lab));
    return raw;
  }

  /** @brief Finishes a marked hue cell from its perceptual base.
   * @param time Borrowed neighboring base rows.
   * @param value Palette coordinate.
   * @param hue Rotation in turns.
   * @return Gamut-mapped RGB16. */
  HS_HOT_INLINE static Pixel direct_hue_color(const PaletteLabTimeView &time,
                                              float value, float hue) {
    return linrgb_to_pixel(direct_hue_rgb(time, value, hue));
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
  /** @brief Copies the cache and binds the Lab view to its own nodes.
   * @param other Cache to copy, including borrowed RGB rows. */
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
  HS_FLASH_MEMBER void bake(const BakedPalette &from, const BakedPalette &to) {
    palette_from = &from;
    palette_to = &to;
    rgb_nodes = nullptr;
    rgb_masks = nullptr;
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

  void bake(const BakedPalette &&, const BakedPalette &) = delete;
  void bake(const BakedPalette &, const BakedPalette &&) = delete;
  void bake(const BakedPalette &&, const BakedPalette &&) = delete;
};

/// Default interval palette cache resolution.
using PaletteLabLut = PaletteLabLutStorage<64>;

/** @brief Interval palette bases baked through the final hue and gamut mapping.
 * @tparam MidpointOnly Whether storage is limited to one admitted midpoint row. */
template <bool MidpointOnly> struct PaletteRgbHueLutImpl {
  static_assert(PaletteLabLut::VALUES == HueRotationLutView::VALUE_STEPS,
                "palette value nodes must match hue rotation table rows");
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
  std::array<Pixel, RGB_TIMES * HueRotationLutView::SIZE>
      nodes{}; ///< Quantized RGB16 code units by time, value and hue.
  std::array<HueRotationGamutMask, RGB_TIMES>
      masks{}; ///< Direct cells by time row.
  /** @brief Bakes hue colors, reusing exact base rows or midpoint rows within
   * `MAX_REUSE_DELTA` of their original baked key.
   * @param base Interval base receiving this cache binding.
   * @param previous_base Prior interval's immutable base, when retained.
   * @param previous Prior interval's RGB owner, when retained. */
  HS_FLASH_MEMBER void bake(PaletteLabLut &base,
                            const PaletteLabLut *previous_base = nullptr,
                            const PaletteRgbHueLutImpl *previous = nullptr) {
    base.rgb_nodes = nullptr;
    base.rgb_masks = nullptr;
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
    std::array<float, HueRotationLutView::HUE_STEPS> cosines, sines;
    for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue)
      hue_sincos(float(hue) / HueRotationLutView::HUE_STEPS, cosines[hue],
                 sines[hue]);
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
        std::memcpy(nodes.data() + MIDPOINT_ROW * HueRotationLutView::SIZE,
                    previous->nodes.data() +
                        MIDPOINT_ROW * HueRotationLutView::SIZE,
                    HueRotationLutView::SIZE * sizeof(Pixel));
        baked_midpoint = previous->baked_midpoint;
        masks[MIDPOINT_ROW] = previous->masks[MIDPOINT_ROW];
        HueRotationGamutBuilder builder(&masks[MIDPOINT_ROW]);
        for (int value = 0; value < PaletteLabLut::VALUES; ++value)
          for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
            const auto &node =
                nodes[MIDPOINT_ROW * HueRotationLutView::SIZE +
                      value * HueRotationLutView::HUE_STEPS + hue];
            mark_hue_gamut_node(builder, value, hue,
                                base.nodes[PaletteLabLut::VALUES + value],
                                cosines[hue], sines[hue], node);
          }
        midpoint_active = true;
        base.rgb_nodes = nodes.data();
        base.rgb_masks = masks.data();
        return;
      }
    }
    const auto bake_row = [&](int row, int base_row) {
      masks[row].clear();
      HueRotationGamutBuilder builder(&masks[row]);
      for (int value = 0; value < PaletteLabLut::VALUES; ++value)
        for (int hue = 0; hue < HueRotationLutView::HUE_STEPS; ++hue) {
          const auto &lab =
              base.nodes[base_row * PaletteLabLut::VALUES + value];
          LinRGB raw;
          const auto NODE = rotate_node(lab, cosines[hue], sines[hue], raw);
          nodes[row * HueRotationLutView::SIZE +
                value * HueRotationLutView::HUE_STEPS + hue] = NODE;
          builder.append(value, hue, raw, NODE);
        }
    };
    if constexpr (MidpointOnly) {
      bake_row(0, 1);
    } else {
      const bool REUSE =
          !base.rgb_midpoint && previous_base && previous &&
          previous_base->active && !previous_base->rgb_midpoint &&
          previous_base->rgb_nodes == previous->nodes.data() &&
          std::memcmp(base.nodes.data(),
                      previous_base->nodes.data() +
                          (PaletteLabLut::TIMES - 1) * PaletteLabLut::VALUES,
                      PaletteLabLut::VALUES * sizeof(OKLab)) == 0;
      if (REUSE) {
        std::memcpy(nodes.data(),
                    previous->nodes.data() +
                        (PaletteLabLut::TIMES - 1) * HueRotationLutView::SIZE,
                    HueRotationLutView::SIZE * sizeof(Pixel));
        masks[0] = previous->masks[PaletteLabLut::TIMES - 1];
      }
      for (int time = base.rgb_midpoint || REUSE ? 1 : 0;
           time < (base.rgb_midpoint ? 2 : PaletteLabLut::TIMES); ++time)
        bake_row(time, time);
    }
    midpoint_active = base.rgb_midpoint;
    if (midpoint_active)
      std::memcpy(baked_midpoint.data(),
                  base.nodes.data() + PaletteLabLut::VALUES,
                  sizeof(baked_midpoint));
    base.rgb_nodes = nodes.data();
    base.rgb_masks = masks.data();
  }
  /** @brief Final hue node with a prepared rotation.
   * @param lab Baked interval base.
   * @param cosine Hue cosine.
   * @param sine Hue sine.
   * @param raw Receives rotated linear RGB before gamut mapping.
   * @return Gamut-mapped RGB16 node. */
  HS_HOT_INLINE static Pixel rotate_node(const OKLab &lab, float cosine,
                                         float sine, LinRGB &raw) {
    const auto ROTATED = rotate_chroma(lab, cosine, sine);
    raw = oklab_to_linear_rgb(ROTATED);
    LinRGB mapped = raw;
    if (!linear_rgb_in_gamut(raw.r, raw.g, raw.b)) {
      HS_PROFILE_DEEP(gamut_clip);
      mapped = oklab_to_linear_rgb(gamut_scale_to_boundary_lut(ROTATED));
    }
    return linrgb_to_pixel(mapped);
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
    if ((time.rgb_mask_low &&
         time.rgb_mask_low->requires_direct(coordinate.low, rotation.low)) ||
        (time.rgb_mask_high &&
         time.rgb_mask_high->requires_direct(coordinate.low, rotation.low))) {
      return PaletteLabLutBase::direct_hue_color(time, value, hue);
    }
    const Pixel LOW_RGB =
        sample_hue_rotation_lut({time.rgb_low, true}, coordinate, rotation);
    if (time.rgb_low == time.rgb_high)
      return LOW_RGB;
    const Pixel HIGH_RGB =
        sample_hue_rotation_lut({time.rgb_high, true}, coordinate, rotation);
    return LOW_RGB.lerp16(HIGH_RGB, frac_to_q16(time.weight));
  }
};

/** @brief Three-time-row interval hue cache. */
using PaletteRgbHueLut = PaletteRgbHueLutImpl<false>;
