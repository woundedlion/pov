/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "effects/ShaderBall.h"
#include "tests/test_effects.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

constexpr std::array<const char *, 4> CASES{
    {"default", "endpoint_min", "endpoint_max", "interior"}};

void write_u16(FILE *output, uint16_t value) {
  const uint8_t bytes[] = {static_cast<uint8_t>(value),
                           static_cast<uint8_t>(value >> 8)};
  if (std::fwrite(bytes, sizeof(bytes), 1, output) != 1)
    std::abort();
}

void write_u32(FILE *output, uint32_t value) {
  write_u16(output, static_cast<uint16_t>(value));
  write_u16(output, static_cast<uint16_t>(value >> 16));
}

float case_value(const ParamDef &parameter, size_t parameter_index,
                 size_t case_index) {
  if (case_index == 1)
    return parameter.is_bool() ? 0.0f : parameter.min;
  if (case_index == 2)
    return parameter.is_bool() ? 1.0f : parameter.max;
  const float fractions[] = {0.25f, 0.5f, 0.75f};
  const float fraction = fractions[parameter_index % std::size(fractions)];
  if (parameter.is_bool())
    return fraction > 0.5f ? 1.0f : 0.0f;
  return parameter.min + fraction * (parameter.max - parameter.min);
}

template <int W, int H>
bool render_case(FILE *output, size_t preset, size_t case_index) {
  hs_test::effects_tests::reset_effect_globals();
  ShaderBall<W, H> effect;
  effect.init();
  if (!effect.selectPreset(preset))
    return false;
  effect.setAnimationsPaused(true);
  if (case_index != 0) {
    std::vector<std::string> parameters;
    for (const ParamDef &parameter : effect.getParameters())
      if (!parameter.readonly && !parameter.is_enum())
        parameters.emplace_back(parameter.name);
    for (size_t index = 0; index < parameters.size(); ++index) {
      const ParamDef *parameter =
          effect.getParameters().find(parameters[index].c_str());
      if (parameter == nullptr ||
          effect.updateParameter(parameters[index].c_str(),
                                 case_value(*parameter, index, case_index)) !=
              ParamSetResult::APPLIED)
        return false;
    }
  }
  effect.draw_frame();
  effect.advance_display();
  write_u16(output, static_cast<uint16_t>(preset));
  write_u16(output, static_cast<uint16_t>(case_index));
  const Pixel *pixels = effect.display_buffer();
  for (size_t index = 0; index < static_cast<size_t>(W) * H; ++index) {
    write_u16(output, pixels[index].r);
    write_u16(output, pixels[index].g);
    write_u16(output, pixels[index].b);
  }
  return true;
}

template <int W, int H> int capture(const char *path) {
  FILE *output = nullptr;
#ifdef _WIN32
  if (fopen_s(&output, path, "wb") != 0)
    output = nullptr;
#else
  output = std::fopen(path, "wb");
#endif
  if (output == nullptr)
    return 2;
  if (std::fwrite("HSPB", 4, 1, output) != 1)
    return 2;
  write_u16(output, 1);
  write_u16(output, W);
  write_u16(output, H);
  write_u16(output, 12);
  write_u16(output, static_cast<uint16_t>(CASES.size()));
  write_u32(output, 12U * static_cast<uint32_t>(CASES.size()));
  for (size_t preset = 0; preset < 12; ++preset)
    for (size_t case_index = 0; case_index < CASES.size(); ++case_index)
      if (!render_case<W, H>(output, preset, case_index)) {
        std::fclose(output);
        return 3;
      }
  return std::fclose(output) == 0 ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 4 || std::strcmp(argv[1], "--resolution") != 0) {
    std::fprintf(
        stderr, "usage: pullback_capture_native --resolution WxH output.bin\n");
    return 2;
  }
  if (std::strcmp(argv[2], "96x20") == 0)
    return capture<96, 20>(argv[3]);
  if (std::strcmp(argv[2], "288x144") == 0)
    return capture<288, 144>(argv[3]);
  std::fprintf(stderr, "unsupported pullback capture resolution: %s\n",
               argv[2]);
  return 2;
}
