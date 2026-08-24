/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <span>

#include "color/color.h"
#include "math/noise_field.h"

/**
 * @file noise_hue_palette.h
 * @brief Spatial noise-driven hue rotation for arbitrary palette sources.
 */

/** @brief View over palette-coordinate by hue-rotation colors. */
struct HueRotationLutView {
  static constexpr int VALUE_STEPS = 64;
  static constexpr int HUE_STEPS = 16;
  static constexpr size_t SIZE = VALUE_STEPS * HUE_STEPS;

  const Pixel *data;
  bool active;
};

/** @brief View over a cube-map noise field used to select hue rotation. */
struct HueNoiseLutView {
  static constexpr int FACE_COUNT = 6;
  static constexpr int FACE_STEPS = 24;
  static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS;
  static constexpr size_t SIZE = FACE_COUNT * FACE_SIZE;

  const int8_t *data;
  bool active;
};

/**
 * @brief Bakes all palette-coordinate and hue-rotation combinations.
 * @tparam Palette Palette source exposing Color4 get(float) const.
 * @param output Destination LUT.
 * @param palette Palette source to hue-rotate.
 */
template <typename Palette>
HS_FLASH_INLINE inline void
prepare_hue_rotation_lut(std::span<Pixel, HueRotationLutView::SIZE> output,
                         const Palette &palette) {
  constexpr float UNIT_OPEN_MAX = 0x1.fffffep-1f;
  for (int value_index = 0; value_index < HueRotationLutView::VALUE_STEPS;
       ++value_index) {
    const float value =
        UNIT_OPEN_MAX * value_index / (HueRotationLutView::VALUE_STEPS - 1);
    const HueRotateBase base = make_hue_rotate_base(palette.get(value));
    for (int hue_index = 0; hue_index < HueRotationLutView::HUE_STEPS;
         ++hue_index) {
      const float amount =
          static_cast<float>(hue_index) / HueRotationLutView::HUE_STEPS;
      output[value_index * HueRotationLutView::HUE_STEPS + hue_index] =
          hue_rotate_lut_gamut(base, amount).color;
    }
  }
}

__attribute__((always_inline)) inline Vector
hue_noise_face_direction(int face, float u, float v) {
  switch (face) {
  case 0:
    return Vector(1.0f, v, u).normalized();
  case 1:
    return Vector(-1.0f, v, -u).normalized();
  case 2:
    return Vector(u, 1.0f, v).normalized();
  case 3:
    return Vector(u, -1.0f, -v).normalized();
  case 4:
    return Vector(u, v, 1.0f).normalized();
  default:
    return Vector(-u, v, -1.0f).normalized();
  }
}

/**
 * @brief Bakes a seamless cube-map view of a sphere-domain noise field.
 * @param output Destination LUT.
 * @param noise Configured noise source.
 * @param scale Spatial frequency over the sphere.
 * @param phase Loop phase in turns.
 */
HS_FLASH_INLINE inline void
prepare_hue_noise_lut(std::span<int8_t, HueNoiseLutView::SIZE> output,
                      const FastNoiseLite &noise, float scale, float phase) {
  const float angle = TWO_PI_F * wrap_t(phase);
  const Vector loop_offset(NOISE_LOOP_RADIUS * cosf(angle),
                           NOISE_LOOP_RADIUS * sinf(angle), 0.0f);
  constexpr float STEP = 2.0f / (HueNoiseLutView::FACE_STEPS - 1);
  for (int face = 0; face < HueNoiseLutView::FACE_COUNT; ++face) {
    const int face_offset = face * HueNoiseLutView::FACE_SIZE;
    for (int y = 0; y < HueNoiseLutView::FACE_STEPS; ++y) {
      const float v = -1.0f + STEP * y;
      for (int x = 0; x < HueNoiseLutView::FACE_STEPS; ++x) {
        const float u = -1.0f + STEP * x;
        const Vector direction = hue_noise_face_direction(face, u, v);
        const Vector q = scale * direction + loop_offset;
        const float sample =
            hs::clamp(noise.GetNoiseSingle(q.x, q.y, q.z), -1.0f, 1.0f);
        const int quantized =
            static_cast<int>(sample * 127.0f + (sample < 0.0f ? -0.5f : 0.5f));
        output[face_offset + y * HueNoiseLutView::FACE_STEPS + x] =
            static_cast<int8_t>(quantized);
      }
    }
  }
}

/**
 * @brief Samples a hue-rotation LUT with interpolation over both axes.
 * @param view Prepared hue-rotation LUT.
 * @param value Palette coordinate in [0, 1].
 * @param amount Hue rotation in turns.
 * @return Hue-rotated linear RGB color.
 */
__attribute__((always_inline)) inline Pixel
sample_hue_rotation_lut(const HueRotationLutView &view, float value,
                        float amount) {
  const float value_position = value * (HueRotationLutView::VALUE_STEPS - 1);
  const int value_low = static_cast<int>(value_position);
  const int value_high =
      std::min(value_low + 1, HueRotationLutView::VALUE_STEPS - 1);
  const uint16_t value_weight =
      frac_to_q16(value_position - static_cast<float>(value_low));

  const float hue_position = amount * HueRotationLutView::HUE_STEPS;
  int hue_low = static_cast<int>(hue_position);
  if (hue_position < static_cast<float>(hue_low))
    --hue_low;
  const uint16_t hue_weight =
      frac_to_q16(hue_position - static_cast<float>(hue_low));
  const int hue_index_low = hue_low & (HueRotationLutView::HUE_STEPS - 1);
  const int hue_index_high =
      (hue_index_low + 1) & (HueRotationLutView::HUE_STEPS - 1);
  const auto sample_row = [&](int row) {
    const int offset = row * HueRotationLutView::HUE_STEPS;
    return view.data[offset + hue_index_low].lerp16(
        view.data[offset + hue_index_high], hue_weight);
  };
  return sample_row(value_low).lerp16(sample_row(value_high), value_weight);
}

/**
 * @brief Samples the cube-map noise field along a unit direction.
 * @param view Prepared hue-noise LUT.
 * @param direction Unit sphere direction.
 * @return Noise value in [-1, 1].
 */
HS_FLASH_INLINE inline float sample_hue_noise_lut(const HueNoiseLutView &view,
                                                  const Vector &direction) {
  const float ax = fabsf(direction.x);
  const float ay = fabsf(direction.y);
  const float az = fabsf(direction.z);
  int face;
  float u;
  float v;
  if (ax >= ay && ax >= az) {
    const float inverse = 1.0f / ax;
    face = direction.x >= 0.0f ? 0 : 1;
    u = (direction.x >= 0.0f ? direction.z : -direction.z) * inverse;
    v = direction.y * inverse;
  } else if (ay >= az) {
    const float inverse = 1.0f / ay;
    face = direction.y >= 0.0f ? 2 : 3;
    u = direction.x * inverse;
    v = (direction.y >= 0.0f ? direction.z : -direction.z) * inverse;
  } else {
    const float inverse = 1.0f / az;
    face = direction.z >= 0.0f ? 4 : 5;
    u = (direction.z >= 0.0f ? direction.x : -direction.x) * inverse;
    v = direction.y * inverse;
  }

  constexpr float SCALE =
      0.5f * static_cast<float>(HueNoiseLutView::FACE_STEPS - 1);
  const float x_position = (u + 1.0f) * SCALE;
  const float y_position = (v + 1.0f) * SCALE;
  const int x_low =
      std::min(static_cast<int>(x_position), HueNoiseLutView::FACE_STEPS - 2);
  const int y_low =
      std::min(static_cast<int>(y_position), HueNoiseLutView::FACE_STEPS - 2);
  const float x_fraction = x_position - x_low;
  const float y_fraction = y_position - y_low;
  const int offset = face * HueNoiseLutView::FACE_SIZE +
                     y_low * HueNoiseLutView::FACE_STEPS + x_low;
  const float row_low =
      hs::lerp(static_cast<float>(view.data[offset]),
               static_cast<float>(view.data[offset + 1]), x_fraction);
  const float row_high = hs::lerp(
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS]),
      static_cast<float>(view.data[offset + HueNoiseLutView::FACE_STEPS + 1]),
      x_fraction);
  return hs::lerp(row_low, row_high, y_fraction) * (1.0f / 127.0f);
}

/**
 * @brief Palette wrapper whose hue rotation varies over a sphere-domain noise
 *        field.
 * @tparam Source Palette source exposing Color4 get(float) const.
 */
template <typename Source> class NoiseHuePalette {
public:
  NoiseHuePalette() = default;

  /**
   * @brief Constructs a bound noise-hue palette.
   * @param source Base palette.
   * @param hue_rotation_lut Prepared palette-by-hue LUT.
   * @param hue_noise_lut Prepared cube-map noise LUT.
   */
  NoiseHuePalette(const Source *source, const Pixel *hue_rotation_lut,
                  const int8_t *hue_noise_lut) {
    bind(source, hue_rotation_lut, hue_noise_lut);
  }

  /**
   * @brief Binds the source and prepared LUT storage.
   * @param source Base palette.
   * @param hue_rotation_lut Prepared palette-by-hue LUT.
   * @param hue_noise_lut Prepared cube-map noise LUT.
   */
  void bind(const Source *source, const Pixel *hue_rotation_lut,
            const int8_t *hue_noise_lut) {
    HS_CHECK(source != nullptr, "NoiseHuePalette bound to null source");
    HS_CHECK(hue_rotation_lut != nullptr,
             "NoiseHuePalette bound to null hue-rotation LUT");
    HS_CHECK(hue_noise_lut != nullptr,
             "NoiseHuePalette bound to null hue-noise LUT");
    this->source = source;
    hue_rotation = {hue_rotation_lut, true};
    hue_noise = {hue_noise_lut, true};
  }

  /**
   * @brief Resolves the hue rotation at a sphere direction.
   * @param direction Unit sphere direction.
   * @param amount Hue rotation magnitude in turns.
   * @return Signed hue rotation in turns.
   */
  float hue_shift(const Vector &direction, float amount) const {
    assert(source != nullptr && "NoiseHuePalette used before bind()!");
    return amount * noise(direction);
  }

  /**
   * @brief Samples the prepared noise field at a sphere direction.
   * @param direction Non-zero direction.
   * @return Noise value in [-1, 1].
   */
  float noise(const Vector &direction) const {
    assert(source != nullptr && "NoiseHuePalette used before bind()!");
    return sample_hue_noise_lut(hue_noise, direction);
  }

  /**
   * @brief Samples a seamless two-axis UV noise field.
   * @param cos_u Cosine of the first periodic coordinate.
   * @param sin_u Sine of the first periodic coordinate.
   * @param cos_v Cosine of the second periodic coordinate.
   * @param sin_v Sine of the second periodic coordinate.
   * @return Noise value in [-1, 1].
   */
  float noise_uv(float cos_u, float sin_u, float cos_v, float sin_v) const {
    const float u_field = noise(Vector(cos_u, sin_u, cos_v));
    const float v_field = noise(Vector(cos_v, sin_v, sin_u));
    return 0.5f * (u_field + v_field);
  }

  /**
   * @brief Samples the palette with an already-resolved hue rotation.
   * @param value Palette coordinate in [0, 1].
   * @param hue_shift Hue rotation in turns.
   * @return Hue-rotated palette color with source alpha preserved.
   */
  Color4 get(float value, float hue_shift) const {
    assert(source != nullptr && "NoiseHuePalette used before bind()!");
    Color4 color = source->get(value);
    color.color = sample_hue_rotation_lut(hue_rotation, value, hue_shift);
    return color;
  }

  /**
   * @brief Samples the palette with a sphere-domain noise hue rotation.
   * @param value Palette coordinate in [0, 1].
   * @param direction Unit sphere direction.
   * @param amount Hue rotation magnitude in turns.
   * @return Noise-hue-rotated palette color.
   */
  Color4 get(float value, const Vector &direction, float amount) const {
    return get(value, hue_shift(direction, amount));
  }

private:
  const Source *source = nullptr;
  HueRotationLutView hue_rotation{nullptr, false};
  HueNoiseLutView hue_noise{nullptr, false};
};
