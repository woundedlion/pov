/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "tests/test_effects.h"
#include "workbench/shader/chain_host.h"
#include "chain_capture_fixtures.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

#define HS_PULLBACK_PRESET_COUNT(value) constexpr uint16_t PRESET_COUNT = value;
#define HS_PULLBACK_OPERATION(name)
#include "tools/pullback_operations.def"
#undef HS_PULLBACK_OPERATION
#undef HS_PULLBACK_PRESET_COUNT

enum class Operation : uint16_t {
#define HS_PULLBACK_PRESET_COUNT(value)
#define HS_PULLBACK_OPERATION(name) name,
#include "tools/pullback_operations.def"
#undef HS_PULLBACK_OPERATION
#undef HS_PULLBACK_PRESET_COUNT
  COUNT,
};

template <int W, int H> bool selected_pixel(Operation operation, int x, int y) {
  switch (operation) {
  case Operation::FULL_FRAME:
    return true;
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

} // namespace

namespace hs_test::shader_chain_tests {
struct ShaderChainWhiteBox {
  template <int W, int H> static auto context(ShaderChain<W, H> &effect) {
    const auto ctx = effect.make_frame_context(effect.colorize);
    effect.program.prepare(ctx);
    return ctx;
  }
  template <int W, int H>
  static const auto &program(const ShaderChain<W, H> &effect) {
    return effect.program;
  }
  template <int W, int H>
  static const auto &hue_noise(const ShaderChain<W, H> &effect) {
    return effect.resources->hue_noise;
  }
  template <int W, int H>
  static void oracle_phase(ShaderChain<W, H> &effect, float phase) {
    const_cast<Pullback::Interp::Op::ColorClockState *>(
        static_cast<const Pullback::Interp::Op::ColorClockState *>(
            effect.program.state_block(
                static_cast<size_t>(effect.colorize.index))))
        ->hue_noise_phase = phase;
  }
};
} // namespace hs_test::shader_chain_tests

namespace {

struct Instruction {
  uint16_t preset;
  Operation operation;
  std::string name;
};

struct OracleInstruction {
  std::string oracle;
  uint16_t preset;
  Operation operation;
  float hue_noise_phase;
};

struct OperationSet {
  std::vector<Instruction> frames;
  std::vector<OracleInstruction> oracles;
};

struct OracleMetric {
  std::string oracle;
  uint16_t maximum = 0;
  uint32_t samples = 0;
};

struct RecordMetadata {
  uint16_t source = UINT16_MAX;
  uint16_t destination = UINT16_MAX;
  uint16_t elapsed = 0;
  uint16_t duration = 0;
};

bool read_u16(FILE *input, uint16_t &value) {
  uint8_t bytes[2];
  if (std::fread(bytes, sizeof(bytes), 1, input) != 1)
    return false;
  value = static_cast<uint16_t>(bytes[0] | bytes[1] << 8);
  return true;
}

bool read_u32(FILE *input, uint32_t &value) {
  uint16_t low;
  uint16_t high;
  if (!read_u16(input, low) || !read_u16(input, high))
    return false;
  value = low | static_cast<uint32_t>(high) << 16;
  return true;
}

void write_u16(FILE *output, uint16_t value) {
  const uint8_t bytes[] = {static_cast<uint8_t>(value),
                           static_cast<uint8_t>(value >> 8)};
  if (std::fwrite(bytes, sizeof(bytes), 1, output) != 1) {
    std::fprintf(stderr, "pullback output: write failed\n");
    std::abort();
  }
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

OperationSet read_instructions(const char *path) {
  FILE *input = open_file(path, "rb");
  if (input == nullptr) {
    std::fprintf(stderr, "pullback read: cannot open %s\n", path);
    return {};
  }
  char magic[4];
  uint16_t version;
  if (std::fread(magic, sizeof(magic), 1, input) != 1 ||
      std::memcmp(magic, "HSPO", sizeof(magic)) != 0 ||
      !read_u16(input, version) || version != 2) {
    std::fprintf(stderr, "pullback read: invalid header in %s\n", path);
    std::fclose(input);
    return {};
  }
  uint16_t count;
  uint16_t oracle_count;
  if (!read_u16(input, count) || !read_u16(input, oracle_count)) {
    std::fprintf(stderr, "pullback read: missing counts in %s\n", path);
    std::fclose(input);
    return {};
  }
  OperationSet operations;
  operations.frames.reserve(count);
  operations.oracles.reserve(oracle_count);
  for (uint16_t index = 0; index < count; ++index) {
    uint16_t preset;
    uint16_t operation_value;
    uint16_t length;
    if (!read_u16(input, preset) || !read_u16(input, operation_value) ||
        !read_u16(input, length)) {
      std::fprintf(stderr, "pullback read: truncated frame %u in %s\n", index,
                   path);
      std::fclose(input);
      return {};
    }
    const auto operation = static_cast<Operation>(operation_value);
    std::string name(length, '\0');
    if (length == 0 || std::fread(name.data(), length, 1, input) != 1 ||
        preset >= PRESET_COUNT || operation >= Operation::COUNT) {
      std::fprintf(
          stderr,
          "pullback read: invalid frame %u preset=%u operation=%u name=%s in %s\n",
          index, preset, operation_value, name.c_str(), path);
      std::fclose(input);
      return {};
    }
    operations.frames.push_back({preset, operation, std::move(name)});
  }
  for (uint16_t index = 0; index < oracle_count; ++index) {
    uint16_t length;
    if (!read_u16(input, length)) {
      std::fprintf(stderr,
                   "pullback read: missing oracle %u name length in %s\n",
                   index, path);
      std::fclose(input);
      return {};
    }
    std::string oracle(length, '\0');
    if (length == 0 || std::fread(oracle.data(), length, 1, input) != 1) {
      std::fprintf(stderr, "pullback read: invalid oracle %u name in %s\n",
                   index, path);
      std::fclose(input);
      return {};
    }
    uint16_t preset;
    uint16_t operation_value;
    uint32_t phase_bits;
    if (!read_u16(input, preset) || !read_u16(input, operation_value) ||
        !read_u32(input, phase_bits)) {
      std::fprintf(stderr, "pullback read: truncated oracle %s in %s\n",
                   oracle.c_str(), path);
      std::fclose(input);
      return {};
    }
    const auto operation = static_cast<Operation>(operation_value);
    const float phase = std::bit_cast<float>(phase_bits);
    if (preset >= PRESET_COUNT || operation >= Operation::COUNT) {
      std::fprintf(
          stderr,
          "pullback read: invalid oracle %s preset=%u operation=%u in %s\n",
          oracle.c_str(), preset, operation_value, path);
      std::fclose(input);
      return {};
    }
    operations.oracles.push_back({std::move(oracle), preset, operation, phase});
  }
  const bool complete = std::fgetc(input) == EOF;
  std::fclose(input);
  if (!complete)
    std::fprintf(stderr, "pullback read: trailing data in %s\n", path);
  return complete ? operations : OperationSet{};
}

const ChainCaptureFixtures::Case *
find_fixture(const char *kind, const std::string &name, int width, int height,
             uint16_t preset, Operation operation, float phase = 0.0f) {
  for (const auto &fixture : ChainCaptureFixtures::CASES)
    if (std::string_view(fixture.kind) == kind && fixture.name == name &&
        fixture.width == width && fixture.height == height &&
        fixture.preset == preset &&
        fixture.operation == static_cast<uint16_t>(operation) &&
        std::bit_cast<uint32_t>(fixture.phase) ==
            std::bit_cast<uint32_t>(phase))
      return &fixture;
  return nullptr;
}

template <int W, int H>
bool render_instruction(const Instruction &instruction,
                        std::vector<Pixel> &pixels, RecordMetadata &metadata) {
  const auto *fixture = find_fixture("frame", instruction.name, W, H,
                                     instruction.preset, instruction.operation);
  if (!fixture)
    return false;
  hs_test::effects_tests::reset_effect_globals();
  ShaderChain<W, H> effect;
  effect.init();
  if (effect.restore_snapshot(ChainCaptureFixtures::snapshot(
          fixture->snapshot_index)) != ChainSnapshotRestoreResult::APPLIED)
    return false;
  metadata = {fixture->source, fixture->destination, fixture->elapsed,
              fixture->duration};
  effect.draw_frame();
  effect.advance_display();
  const Pixel *display = effect.display_buffer();
  pixels.assign(display, display + static_cast<size_t>(W) * H);
  return true;
}

Color4 exact_colorize(const Pullback::FieldSample &sample,
                      const Pullback::Color::GeneratedPaletteState &state,
                      const FastNoiseLite &noise, float noise_scale,
                      float noise_phase) {
  using namespace Pullback;
  const float value = Color::palette_mapping_coordinate(
      sample.value, state.mapping, state.mapping_frequency,
      state.mapping_offset);
  Color4 color = state.palette->get(value);
  if (state.hue_rotation.active && state.hue_mode == Color::HueMode::NOISE) {
    const auto q =
        math::noise_sphere_coordinate(sample.sphere, noise_scale, noise_phase);
    const float amount =
        state.hue_shift_amount * noise.GetNoiseSingle(q.x, q.y, q.z);
    color = hue_rotate_lut_gamut(make_hue_rotate_base(color), amount);
  } else if (state.hue_rotation.active) {
    const float amount =
        math::wrap_t(state.hue_shift_amount * sample.path_length);
    if (amount != 0.0f)
      color = hue_rotate_lut_gamut(make_hue_rotate_base(color), amount);
  }
  color.color =
      color.color * Color::brightness_envelope_gain(
                        sample.value, state.brightness_envelope,
                        state.brightness_bottom, state.brightness_top);
  color.alpha *= sample.coverage *
                 hs::lerp(state.opacity_low, state.opacity_high, sample.value);
  return color;
}

template <int W, int H>
Color4 exact_shade(const ShaderChain<W, H> &effect, const math::Vector &view,
                   const Pullback::Interp::FrameContext &ctx,
                   const std::string &oracle, size_t &substitutions) {
  using namespace Pullback;
  using namespace Pullback::Interp;
  using WB = hs_test::shader_chain_tests::ShaderChainWhiteBox;
  const auto &program = WB::program(effect);
  alignas(SLOT_ALIGN) uint8_t slot_a[SLOT_SIZE], slot_b[SLOT_SIZE];
  void *input = slot_a, *output = slot_b;
  ::new (input) SphereSample{view, 0.0f};
  const auto ops = program.ops();
  for (size_t index = 0; index < ops.size(); ++index) {
    const auto &op = *ops[index].op;
    const auto *params = program.param_block(index);
    const auto *prepared = program.prepared_block(index);
    if (oracle == "PEIRCE_FAST_SQUARE" &&
        std::string_view(op.operator_id) == Op::ProjectPeirceSquareFastV3::ID) {
      ++substitutions;
      const auto &p =
          *reinterpret_cast<const Op::ProjectPeirceSquareFastV3::Params *>(
              params);
      const auto &frame =
          *reinterpret_cast<const Op::ProjectOrientation *>(prepared);
      const auto &sphere = *static_cast<const SphereSample *>(input);
      const auto local = math::rotate(sphere.dir, frame.conjugate);
      ::new (output) PlaneSample{Kernel::project(
          sphere, local,
          Pullback::Projection::peirce(local, 0.0f, Op::PEIRCE_SQUARE_LAYOUT,
                                       0.0f, true, p.coordinate_scale,
                                       p.singularity_fade))};
    } else if (oracle == "HUE_ROTATION_AND_NOISE_LUTS" &&
               std::string_view(op.operator_id) ==
                   Op::ColorizeGeneratedPaletteV3::ID) {
      ++substitutions;
      const auto &p =
          *reinterpret_cast<const Op::GeneratedPaletteParams *>(params);
      const auto &clock =
          *static_cast<const Op::ColorClockState *>(program.state_block(index));
      ::new (output) Color4{exact_colorize(
          *static_cast<const FieldSample *>(input),
          *reinterpret_cast<const Color::GeneratedPaletteState *>(prepared),
          WB::hue_noise(effect), p.hue_noise_scale, clock.hue_noise_phase)};
    } else {
      op.runtime.run(input, output, ctx, params, prepared);
    }
    std::swap(input, output);
  }
  return *static_cast<const Color4 *>(input);
}

template <int W, int H>
bool measure_oracle(ShaderChain<W, H> &effect,
                    const OracleInstruction &instruction, uint16_t &maximum,
                    uint32_t &samples) {
  using WB = hs_test::shader_chain_tests::ShaderChainWhiteBox;
  const auto *fixture =
      find_fixture("oracle", instruction.oracle, W, H, instruction.preset,
                   instruction.operation, instruction.hue_noise_phase);
  if (!fixture ||
      effect.restore_snapshot(ChainCaptureFixtures::snapshot(
          fixture->snapshot_index)) != ChainSnapshotRestoreResult::APPLIED)
    return false;
  WB::oracle_phase(effect, instruction.hue_noise_phase);
  const auto ctx = WB::context(effect);
  const auto &program = WB::program(effect);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      if (!selected_pixel<W, H>(instruction.operation, x, y))
        continue;
      const auto view = math::pixel_to_vector<W, H>(x, y);
      const auto optimized = program.evaluate(view, ctx);
      size_t substitutions = 0;
      const auto exact =
          exact_shade(effect, view, ctx, instruction.oracle, substitutions);
      if (substitutions == 0)
        return false;
      const auto actual = optimized.color * optimized.alpha;
      const auto expected = exact.color * exact.alpha;
      const auto error = [](uint16_t a, uint16_t b) {
        return static_cast<uint16_t>(a > b ? a - b : b - a);
      };
      maximum = std::max(maximum, std::max({error(actual.r, expected.r),
                                            error(actual.g, expected.g),
                                            error(actual.b, expected.b)}));
      samples += 3;
    }
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
  uint32_t selected_count = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      selected_count += selected_pixel<W, H>(instruction.operation, x, y);
  write_u32(output, selected_count);
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
  const OperationSet operations = read_instructions(operations_path);
  if (operations.frames.empty() || operations.oracles.empty()) {
    std::fprintf(stderr,
                 "pullback capture: no frames or oracles loaded from %s\n",
                 operations_path);
    return 2;
  }
  std::vector<OracleMetric> metrics;
  for (const OracleInstruction &instruction : operations.oracles) {
    auto metric = std::find_if(metrics.begin(), metrics.end(),
                               [&](const OracleMetric &candidate) {
                                 return candidate.oracle == instruction.oracle;
                               });
    if (metric == metrics.end()) {
      metrics.push_back({instruction.oracle});
      metric = metrics.end() - 1;
    }
    hs_test::effects_tests::reset_effect_globals();
    ShaderChain<W, H> effect;
    effect.init();
    if (!measure_oracle<W, H>(effect, instruction, metric->maximum,
                              metric->samples)) {
      std::fprintf(
          stderr,
          "pullback oracle: measurement failed for %s preset=%u operation=%u\n",
          instruction.oracle.c_str(), instruction.preset,
          static_cast<unsigned>(instruction.operation));
      return 3;
    }
  }
  std::unique_ptr<FILE, decltype(&std::fclose)> output(
      open_file(output_path, "wb"), &std::fclose);
  if (output == nullptr) {
    std::fprintf(stderr, "pullback output: cannot open %s\n", output_path);
    return 2;
  }
  if (std::fwrite("HSPB", 4, 1, output.get()) != 1) {
    std::fprintf(stderr, "pullback output: cannot write header to %s\n",
                 output_path);
    return 2;
  }
  write_u16(output.get(), 3);
  write_u16(output.get(), W);
  write_u16(output.get(), H);
  write_u32(output.get(), static_cast<uint32_t>(operations.frames.size()));
  write_u16(output.get(), static_cast<uint16_t>(metrics.size()));
  for (const Instruction &instruction : operations.frames) {
    std::vector<Pixel> pixels;
    RecordMetadata metadata;
    if (!render_instruction<W, H>(instruction, pixels, metadata)) {
      std::fprintf(
          stderr,
          "pullback capture: no fixture or snapshot restore failed for "
          "%s preset=%u operation=%u\n",
          instruction.name.c_str(), instruction.preset,
          static_cast<unsigned>(instruction.operation));
      return 3;
    }
    if (!write_record<W, H>(output.get(), instruction, pixels, metadata)) {
      std::fprintf(
          stderr,
          "pullback output: cannot write %s preset=%u operation=%u to %s\n",
          instruction.name.c_str(), instruction.preset,
          static_cast<unsigned>(instruction.operation), output_path);
      return 3;
    }
  }
  for (const OracleMetric &metric : metrics) {
    write_u16(output.get(), static_cast<uint16_t>(metric.oracle.size()));
    if (std::fwrite(metric.oracle.data(), metric.oracle.size(), 1,
                    output.get()) != 1) {
      std::fprintf(stderr, "pullback output: cannot write oracle %s to %s\n",
                   metric.oracle.c_str(), output_path);
      return 2;
    }
    write_u16(output.get(), metric.maximum);
    write_u32(output.get(), metric.samples);
  }
  if (std::fclose(output.release()) != 0) {
    std::fprintf(stderr, "pullback output: cannot close %s\n", output_path);
    return 2;
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 6 || std::strcmp(argv[1], "--resolution") != 0 ||
      std::strcmp(argv[3], "--operations") != 0) {
    std::fprintf(stderr,
                 "usage: %s --resolution WxH "
                 "--operations operations.bin output.bin\n",
                 argv[0]);
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
