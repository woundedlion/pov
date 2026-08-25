/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

#include "tests/test_effects.h"
#include "workbench/ShaderWorkbench.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

static_assert(!HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND,
              "pullback capture must exercise compiled programs");

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

// The parameter cases are selected as an ordinal range below.
static_assert(static_cast<uint16_t>(Operation::CASE_DEFAULT) == 0 &&
                  static_cast<uint16_t>(Operation::CASE_INTERIOR) == 3,
              "CASE_DEFAULT..CASE_INTERIOR must lead pullback_operations.def");

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

namespace hs_test::shader_workbench_tests {

struct ShaderWorkbenchWhiteBox {
  template <typename Selected>
  static const auto &selected_config(const Selected &selected) {
    if constexpr (requires { selected.config; })
      return selected.config;
    else
      return selected;
  }

  template <typename SB, typename Selected>
  static auto selected_pipeline(const Selected &selected) {
    if constexpr (requires { selected.pipeline; })
      return selected.pipeline;
    else
      return SB::resolve_pipeline_id(selected);
  }

  template <typename SB, typename Selected>
  static const auto *selected_program(const Selected &selected) {
    if constexpr (requires { selected.pipeline; })
      return SB::get_inverse_program(selected.pipeline);
    else
      return SB::find_inverse_program(selected);
  }

  template <typename FieldSample>
  static float sample_path_length(const FieldSample &sample) {
    if constexpr (requires { sample.path_length; })
      return sample.path_length;
    else
      return sample.warp_displacement;
  }

  template <int W, int H>
  static bool force_transition(ShaderWorkbench<W, H> &effect, size_t source,
                               size_t destination, uint16_t elapsed,
                               uint16_t duration) {
    using SB = ShaderWorkbench<W, H>;
    if (!effect.selectPreset(source))
      return false;
    effect.setAnimationsPaused(true);
    const auto &from = SB::PRESETS[source];
    const auto &to = SB::PRESETS[destination];
    const auto &from_config = selected_config(from);
    const auto &to_config = selected_config(to);
    if (!effect.prepare_resource_union(from_config, to_config))
      return false;
    effect.state->param_morph.active = false;
    effect.state->transition = {from_config,
                                to_config,
                                effect.runtime,
                                effect.runtime,
                                elapsed,
                                duration,
                                false,
                                true,
                                selected_pipeline<SB>(from),
                                selected_pipeline<SB>(to)};
    return true;
  }

  template <int W, int H>
  static bool selected_pipeline_active(const ShaderWorkbench<W, H> &effect,
                                       size_t preset) {
    using SB = ShaderWorkbench<W, H>;
    return effect.active_pipeline == selected_pipeline<SB>(SB::PRESETS[preset]);
  }

  template <int W, int H>
  static bool measure_oracle(ShaderWorkbench<W, H> &effect,
                             const std::string &oracle, size_t preset,
                             float hue_noise_phase, Operation operation,
                             uint16_t &maximum, uint32_t &samples) {
    using SB = ShaderWorkbench<W, H>;
    if (preset >= SB::PRESETS.size())
      return false;
    effect.runtime.clocks.hue_noise_phase = hue_noise_phase;
    const auto &selected = SB::PRESETS[preset];
    const auto &config = selected_config(selected);
    if (!effect.prepare_resource_union(config, config))
      return false;
    const auto frame = effect.prepare_frame(config, effect.runtime);
    const auto *program = selected_program<SB>(selected);
    if (program == nullptr || !program->resources_ready(frame))
      return false;
    const auto optimized = program->shade;
    alignas(std::max_align_t)
        std::byte prepared_storage[SB::PREPARED_BLOB_BYTES];
    program->prepare(frame, prepared_storage);
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        if (!selected_pixel<W, H>(operation, x, y))
          continue;
        const Vector view = pixel_to_vector<W, H>(x, y);
        const Color4 optimized_color = optimized(view, frame, prepared_storage);
        Color4 exact_color;
        if (oracle == "PEIRCE_FAST_SQUARE")
          exact_color = exact_peirce_shade<SB>(view, frame);
        else if (oracle == "HUE_ROTATION_AND_NOISE_LUTS")
          exact_color = exact_hue_shade<SB>(view, frame);
        else
          return false;
        const Pixel actual = optimized_color.color * optimized_color.alpha;
        const Pixel expected = exact_color.color * exact_color.alpha;
        maximum =
            std::max(maximum, std::max({channel_error(actual.r, expected.r),
                                        channel_error(actual.g, expected.g),
                                        channel_error(actual.b, expected.b)}));
        samples += 3;
      }
    }
    return true;
  }

private:
  static uint16_t channel_error(uint16_t a, uint16_t b) {
    return a > b ? a - b : b - a;
  }

  template <typename SB>
  static Color4 exact_peirce_shade(const Vector &view,
                                   const typename SB::FrameState &frame) {
    const Vector outer_local = SB::outer_camera_lookup(view, frame);
    const Vector lensed = lenses::dodecahedral_kaleidoscope_lens(outer_local);
    const Vector local = rotate(lensed, frame.transforms.projection_conj);
    const Pullback::ProjectionResult result = Pullback::Projection::peirce(
        local, 0.0f, 1, 0.0f, true, frame.params.projection.coordinate_scale,
        frame.params.projection.singularity_fade);
    const typename SB::ProjectedLookup projected{
        result.coords, result.provenance, local, 0.0f};
    return SB::shade_projected(projected, frame);
  }

  template <typename SB>
  static Color4 exact_hue_shade(const Vector &view,
                                const typename SB::FrameState &frame) {
    const Vector outer_local = SB::outer_camera_lookup(view, frame);
    const auto projected = SB::surface_lens_project_lookup(outer_local, frame);
    const auto warped = SB::planar_warp_lookup(projected, frame);
    const Complex source_coords =
        SB::condition_source_coords(warped.coords, frame);
    const float field = SB::sample_source(source_coords, projected, frame);
    const auto material = SB::shape_material(field, projected, warped, frame);
    return exact_colorize<SB>(material, frame);
  }

  template <typename SB>
  static Color4 exact_colorize(const typename SB::FieldSample &sample,
                               const typename SB::FrameState &frame) {
    const float oscillation =
        frame.params.color.phase_oscillation_depth *
        fast_sinf(TWO_PI_F * frame.clocks.palette_oscillation_phase);
    const float palette_value = SB::palette_mapping_coordinate(
        sample.value, frame.slots.palette_mapping,
        frame.params.color.mapping_frequency,
        frame.params.color.mapping_phase + oscillation);
    Color4 color = frame.resources.generated_palette->get(palette_value);
    if (frame.prepared_hue_rotation.active &&
        frame.slots.hue_shift == SB::HueShiftMode::NOISE) {
      const Vector q = noise_sphere_coordinate(
          sample.sphere, frame.params.color.hue_noise_scale,
          frame.clocks.hue_noise_phase);
      const float noise =
          frame.resources.color_noise->GetNoiseSingle(q.x, q.y, q.z);
      const float amount = frame.params.color.hue_shift_amount * noise;
      const HueRotateBase base = make_hue_rotate_base(color);
      color = hue_rotate_lut_gamut(base, amount);
    } else if (frame.prepared_hue_rotation.active) {
      const float amount = wrap_t(frame.params.color.hue_shift_amount *
                                  sample_path_length(sample));
      if (amount != 0.0f)
        color = hue_rotate_lut_gamut(make_hue_rotate_base(color), amount);
    }
    color.color =
        color.color * SB::brightness_envelope_gain(
                          sample.value, frame.slots.brightness_envelope,
                          frame.params.color.brightness_depth);
    color.alpha *= sample.coverage * hs::lerp(frame.params.color.opacity_low,
                                              frame.params.color.opacity_high,
                                              sample.value);
    return color;
  }
};

} // namespace hs_test::shader_workbench_tests

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

OperationSet read_instructions(const char *path) {
  FILE *input = open_file(path, "rb");
  if (input == nullptr)
    return {};
  char magic[4];
  uint16_t version;
  if (std::fread(magic, sizeof(magic), 1, input) != 1 ||
      std::memcmp(magic, "HSPO", sizeof(magic)) != 0 ||
      !read_u16(input, version) || version != 2) {
    std::fclose(input);
    return {};
  }
  uint16_t count;
  uint16_t oracle_count;
  if (!read_u16(input, count) || !read_u16(input, oracle_count)) {
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
      std::fclose(input);
      return {};
    }
    const auto operation = static_cast<Operation>(operation_value);
    std::string name(length, '\0');
    if (length == 0 || std::fread(name.data(), length, 1, input) != 1 ||
        preset >= PRESET_COUNT || operation >= Operation::COUNT) {
      std::fclose(input);
      return {};
    }
    operations.frames.push_back({preset, operation, std::move(name)});
  }
  for (uint16_t index = 0; index < oracle_count; ++index) {
    uint16_t length;
    if (!read_u16(input, length)) {
      std::fclose(input);
      return {};
    }
    std::string oracle(length, '\0');
    if (length == 0 || std::fread(oracle.data(), length, 1, input) != 1) {
      std::fclose(input);
      return {};
    }
    uint16_t preset;
    uint16_t operation_value;
    uint32_t phase_bits;
    if (!read_u16(input, preset) || !read_u16(input, operation_value) ||
        !read_u32(input, phase_bits)) {
      std::fclose(input);
      return {};
    }
    const auto operation = static_cast<Operation>(operation_value);
    const float phase = std::bit_cast<float>(phase_bits);
    if (preset >= PRESET_COUNT || operation >= Operation::COUNT) {
      std::fclose(input);
      return {};
    }
    operations.oracles.push_back({std::move(oracle), preset, operation, phase});
  }
  const bool complete = std::fgetc(input) == EOF;
  std::fclose(input);
  return complete ? operations : OperationSet{};
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

template <int W, int H>
bool render_instruction(const Instruction &instruction,
                        std::vector<Pixel> &pixels, RecordMetadata &metadata) {
  hs_test::effects_tests::reset_effect_globals();
  ShaderWorkbench<W, H> effect;
  effect.init();
  const bool transition_from =
      instruction.operation == Operation::THROUGH_CLEAR_FROM;
  const bool transition_to =
      instruction.operation == Operation::THROUGH_CLEAR_TO;
  if (transition_from || transition_to) {
    constexpr uint16_t DURATION = 60;
    const uint16_t source =
        transition_from
            ? instruction.preset
            : (instruction.preset + PRESET_COUNT - 1) % PRESET_COUNT;
    const uint16_t destination = transition_from
                                     ? (instruction.preset + 1) % PRESET_COUNT
                                     : instruction.preset;
    const uint16_t elapsed = transition_from ? 0 : DURATION;
    if (!hs_test::shader_workbench_tests::ShaderWorkbenchWhiteBox::
            force_transition(effect, source, destination, elapsed, DURATION))
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
  if (!transition_from && !transition_to &&
      !hs_test::shader_workbench_tests::ShaderWorkbenchWhiteBox::
          selected_pipeline_active(effect, instruction.preset))
    return false;
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
  const OperationSet operations = read_instructions(operations_path);
  if (operations.frames.empty() || operations.oracles.empty())
    return 2;
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
    ShaderWorkbench<W, H> effect;
    effect.init();
    if (!hs_test::shader_workbench_tests::ShaderWorkbenchWhiteBox::
            measure_oracle(effect, instruction.oracle, instruction.preset,
                           instruction.hue_noise_phase, instruction.operation,
                           metric->maximum, metric->samples))
      return 3;
  }
  std::unique_ptr<FILE, decltype(&std::fclose)> output(
      open_file(output_path, "wb"), &std::fclose);
  if (output == nullptr)
    return 2;
  if (std::fwrite("HSPB", 4, 1, output.get()) != 1)
    return 2;
  write_u16(output.get(), 3);
  write_u16(output.get(), W);
  write_u16(output.get(), H);
  write_u32(output.get(), static_cast<uint32_t>(operations.frames.size()));
  write_u16(output.get(), static_cast<uint16_t>(metrics.size()));
  for (const Instruction &instruction : operations.frames) {
    std::vector<Pixel> pixels;
    RecordMetadata metadata;
    if (!render_instruction<W, H>(instruction, pixels, metadata) ||
        !write_record<W, H>(output.get(), instruction, pixels, metadata)) {
      return 3;
    }
  }
  for (const OracleMetric &metric : metrics) {
    write_u16(output.get(), static_cast<uint16_t>(metric.oracle.size()));
    if (std::fwrite(metric.oracle.data(), metric.oracle.size(), 1,
                    output.get()) != 1)
      return 2;
    write_u16(output.get(), metric.maximum);
    write_u32(output.get(), metric.samples);
  }
  return std::fclose(output.release()) == 0 ? 0 : 2;
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
