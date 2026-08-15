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

namespace hs_test::shaderball_tests {

struct ShaderBallWhiteBox {
  template <int W, int H>
  static bool force_transition(ShaderBall<W, H> &effect, size_t source,
                               size_t destination, uint16_t elapsed,
                               uint16_t duration) {
    using SB = ShaderBall<W, H>;
    if (!effect.selectPreset(source))
      return false;
    effect.setAnimationsPaused(true);
    const auto &from = SB::PRESETS[source];
    const auto &to = SB::PRESETS[destination];
    if (!effect.prepare_resource_union(from, to))
      return false;
    effect.state->param_morph.active = false;
    effect.state->transition = {from,           to,      effect.runtime,
                                effect.runtime, elapsed, duration,
                                false,          true};
    effect.state->transition.from_pipeline = SB::resolve_pipeline_id(from);
    effect.state->transition.to_pipeline = SB::resolve_pipeline_id(to);
    return true;
  }
};

} // namespace hs_test::shaderball_tests

namespace {

enum class Operation : uint16_t {
  CASE_DEFAULT,
  CASE_ENDPOINT_MIN,
  CASE_ENDPOINT_MAX,
  CASE_INTERIOR,
  COLUMN_ZERO,
  WRAP_COLUMNS,
  NORTH_ROW,
  SOUTH_ROW,
  POLE_BANDS,
  HORIZON_COLUMNS,
  OCTANT_COLUMNS,
  MIRROR_GRID,
  CARDINAL_POINTS,
  FRAME_PERIMETER,
  EQUATOR_ROW,
  FRONT_AXIS_POINT,
  THROUGH_CLEAR_FROM,
  THROUGH_CLEAR_TO,
  COUNT,
};

struct Instruction {
  uint16_t preset;
  Operation operation;
  std::string name;
};

struct RecordMetadata {
  uint16_t source = UINT16_MAX;
  uint16_t destination = UINT16_MAX;
  uint16_t elapsed = 0;
  uint16_t duration = 0;
};

uint16_t read_u16(FILE *input) {
  uint8_t bytes[2];
  if (std::fread(bytes, sizeof(bytes), 1, input) != 1)
    std::abort();
  return static_cast<uint16_t>(bytes[0] | bytes[1] << 8);
}

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

FILE *open_file(const char *path, const char *mode) {
  FILE *file = nullptr;
#ifdef _WIN32
  if (fopen_s(&file, path, mode) != 0)
    return nullptr;
#else
  file = std::fopen(path, mode);
#endif
  return file;
}

std::vector<Instruction> read_instructions(const char *path) {
  FILE *input = open_file(path, "rb");
  if (input == nullptr)
    return {};
  char magic[4];
  if (std::fread(magic, sizeof(magic), 1, input) != 1 ||
      std::memcmp(magic, "HSPO", sizeof(magic)) != 0 || read_u16(input) != 1) {
    std::fclose(input);
    return {};
  }
  const uint16_t count = read_u16(input);
  std::vector<Instruction> instructions;
  instructions.reserve(count);
  for (uint16_t index = 0; index < count; ++index) {
    const uint16_t preset = read_u16(input);
    const auto operation = static_cast<Operation>(read_u16(input));
    const uint16_t length = read_u16(input);
    std::string name(length, '\0');
    if (length == 0 || std::fread(name.data(), length, 1, input) != 1 ||
        preset >= 12 || operation >= Operation::COUNT) {
      std::fclose(input);
      return {};
    }
    instructions.push_back({preset, operation, std::move(name)});
  }
  const bool complete = std::fgetc(input) == EOF;
  std::fclose(input);
  return complete ? instructions : std::vector<Instruction>{};
}

float case_value(const ParamDef &parameter, size_t parameter_index,
                 Operation operation) {
  if (operation == Operation::CASE_ENDPOINT_MIN)
    return parameter.is_bool() ? 0.0f : parameter.min;
  if (operation == Operation::CASE_ENDPOINT_MAX)
    return parameter.is_bool() ? 1.0f : parameter.max;
  constexpr std::array<float, 3> FRACTIONS{{0.25f, 0.5f, 0.75f}};
  const float fraction = FRACTIONS[parameter_index % FRACTIONS.size()];
  if (parameter.is_bool())
    return fraction > 0.5f ? 1.0f : 0.0f;
  return parameter.min + fraction * (parameter.max - parameter.min);
}

template <int W, int H> bool selected_pixel(Operation operation, int x, int y) {
  switch (operation) {
  case Operation::COLUMN_ZERO:
    return x == 0;
  case Operation::WRAP_COLUMNS:
    return x == 0 || x == W - 1;
  case Operation::NORTH_ROW:
    return y == 0;
  case Operation::SOUTH_ROW:
    return y == H - 1;
  case Operation::POLE_BANDS:
    return y < 2 || y >= H - 2;
  case Operation::HORIZON_COLUMNS:
    return x == W / 4 || x == 3 * W / 4;
  case Operation::OCTANT_COLUMNS:
    return x == W / 8 || x == 3 * W / 8 || x == 5 * W / 8 || x == 7 * W / 8;
  case Operation::MIRROR_GRID:
    return x % (W / 8) == 0 || y == H / 4 || y == H / 2 || y == 3 * H / 4;
  case Operation::CARDINAL_POINTS:
    return y == H / 2 && x % (W / 4) == 0;
  case Operation::FRAME_PERIMETER:
    return x == 0 || x == W - 1 || y == 0 || y == H - 1;
  case Operation::EQUATOR_ROW:
    return y == H / 2;
  case Operation::FRONT_AXIS_POINT:
    return x == 0 && y == H / 2;
  default:
    return true;
  }
}

template <int W, int H>
bool render_instruction(const Instruction &instruction,
                        std::vector<Pixel> &pixels, RecordMetadata &metadata) {
  hs_test::effects_tests::reset_effect_globals();
  ShaderBall<W, H> effect;
  effect.init();
  const bool transition_from =
      instruction.operation == Operation::THROUGH_CLEAR_FROM;
  const bool transition_to =
      instruction.operation == Operation::THROUGH_CLEAR_TO;
  if (transition_from || transition_to) {
    constexpr uint16_t DURATION = 60;
    const uint16_t source =
        transition_from ? instruction.preset : (instruction.preset + 11) % 12;
    const uint16_t destination =
        transition_from ? (instruction.preset + 1) % 12 : instruction.preset;
    const uint16_t elapsed = transition_from ? 0 : DURATION;
    if (!hs_test::shaderball_tests::ShaderBallWhiteBox::force_transition(
            effect, source, destination, elapsed, DURATION))
      return false;
    metadata = {source, destination, elapsed, DURATION};
  } else {
    if (!effect.selectPreset(instruction.preset))
      return false;
    effect.setAnimationsPaused(true);
    if (instruction.operation != Operation::CASE_DEFAULT &&
        instruction.operation <= Operation::CASE_INTERIOR) {
      std::vector<std::string> parameters;
      for (const ParamDef &parameter : effect.getParameters())
        if (!parameter.readonly && !parameter.is_enum())
          parameters.emplace_back(parameter.name);
      for (size_t index = 0; index < parameters.size(); ++index) {
        const ParamDef *parameter =
            effect.getParameters().find(parameters[index].c_str());
        if (parameter == nullptr ||
            effect.updateParameter(
                parameters[index].c_str(),
                case_value(*parameter, index, instruction.operation)) !=
                ParamSetResult::APPLIED)
          return false;
      }
    }
  }
  effect.draw_frame();
  effect.advance_display();
  const Pixel *display = effect.display_buffer();
  pixels.assign(display, display + static_cast<size_t>(W) * H);
  return true;
}

template <int W, int H>
bool write_record(FILE *output, const Instruction &instruction,
                  const std::vector<Pixel> &pixels,
                  const RecordMetadata &metadata) {
  write_u16(output, instruction.preset);
  write_u16(output, static_cast<uint16_t>(instruction.operation));
  write_u16(output, static_cast<uint16_t>(instruction.name.size()));
  if (std::fwrite(instruction.name.data(), instruction.name.size(), 1,
                  output) != 1)
    return false;
  write_u16(output, metadata.source);
  write_u16(output, metadata.destination);
  write_u16(output, metadata.elapsed);
  write_u16(output, metadata.duration);
  uint32_t selected = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      selected += selected_pixel<W, H>(instruction.operation, x, y);
  write_u32(output, selected);
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const bool selected = selected_pixel<W, H>(instruction.operation, x, y);
      const Pixel pixel =
          selected ? pixels[static_cast<size_t>(y) * W + x] : Pixel(0, 0, 0);
      write_u16(output, pixel.r);
      write_u16(output, pixel.g);
      write_u16(output, pixel.b);
    }
  }
  return true;
}

template <int W, int H>
int capture(const char *operations_path, const char *output_path) {
  const std::vector<Instruction> instructions =
      read_instructions(operations_path);
  if (instructions.empty())
    return 2;
  FILE *output = open_file(output_path, "wb");
  if (output == nullptr)
    return 2;
  if (std::fwrite("HSPB", 4, 1, output) != 1)
    return 2;
  write_u16(output, 2);
  write_u16(output, W);
  write_u16(output, H);
  write_u32(output, static_cast<uint32_t>(instructions.size()));
  for (const Instruction &instruction : instructions) {
    std::vector<Pixel> pixels;
    RecordMetadata metadata;
    if (!render_instruction<W, H>(instruction, pixels, metadata) ||
        !write_record<W, H>(output, instruction, pixels, metadata)) {
      std::fclose(output);
      return 3;
    }
  }
  return std::fclose(output) == 0 ? 0 : 2;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 6 || std::strcmp(argv[1], "--resolution") != 0 ||
      std::strcmp(argv[3], "--operations") != 0) {
    std::fprintf(stderr, "usage: pullback_capture_backend --resolution WxH "
                         "--operations operations.bin output.bin\n");
    return 2;
  }
  if (std::strcmp(argv[2], "96x20") == 0)
    return capture<96, 20>(argv[4], argv[5]);
  if (std::strcmp(argv[2], "288x144") == 0)
    return capture<288, 144>(argv[4], argv[5]);
  std::fprintf(stderr, "unsupported pullback capture resolution: %s\n",
               argv[2]);
  return 2;
}
