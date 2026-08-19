/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Chain-interpreter core: operator-table integrity, catalog golden pin,
 * static-versus-erased parity over every slice operator and topology variant,
 * transactional refusals, and instance-state identity/migration.
 *
 * Catalog regen: the golden at tests/data/shader_chain_catalog.json is
 * rewritten (instead of compared) by running the module with
 * HS_SHADER_CHAIN_CATALOG_REGEN=1, e.g.
 *   HS_SHADER_CHAIN_CATALOG_REGEN=1 ./run_tests shader_chain
 * then committing the file.
 */
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "core/render/pullback/catalog_export.h"
#include "core/render/pullback/interpreter.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shader_chain_tests {

namespace PB = Pullback;
namespace In = Pullback::Interp;

// --- fixtures -------------------------------------------------------------

/** Two program arenas plus the program bound over them. */
struct ProgramFixture {
  alignas(std::max_align_t) uint8_t block_a[In::CHAIN_ARENA_BYTES];
  alignas(std::max_align_t) uint8_t block_b[In::CHAIN_ARENA_BYTES];
  In::ChainProgram program;

  explicit ProgramFixture(
      size_t capacity = In::CHAIN_ARENA_BYTES,
      std::span<const In::OperatorDescriptor> table =
          std::span<const In::OperatorDescriptor>(In::OPERATOR_TABLE)) {
    program.bind_storage(block_a, block_b, capacity, table);
  }
};

struct GradientSource {
  uint16_t base;
  Color4 get(float t) const {
    return Color4(Pixel(static_cast<uint16_t>(base + t * 40000.0f),
                        static_cast<uint16_t>(t * 30000.0f),
                        static_cast<uint16_t>(50000.0f - t * 40000.0f)),
                  1.0f);
  }
};

/** Engine-owned shared color resources the FrameContext borrows. */
struct ColorResources {
  alignas(std::max_align_t) uint8_t
      palette_bytes[3 * BakedPalette::required_arena_bytes() + 64];
  BakedPalette palettes[3];
  std::array<Pixel, PB::Color::HueRotationLutView::SIZE> hue_rotation{};
  std::array<int8_t, PB::Color::HueNoiseLutView::SIZE> hue_noise{};

  ColorResources() {
    Arena arena(palette_bytes, sizeof(palette_bytes));
    palettes[0].bake(arena, GradientSource{1000});
    palettes[1].bake(arena, GradientSource{9000});
    palettes[2].bake(arena, GradientSource{21000});
    PB::Color::prepare_hue_rotation_lut(
        std::span<Pixel, PB::Color::HueRotationLutView::SIZE>(hue_rotation),
        palettes[0]);
    FastNoiseLite noise;
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(4211);
    noise.SetFrequency(1.0f);
    PB::Color::prepare_hue_noise_lut(
        std::span<int8_t, PB::Color::HueNoiseLutView::SIZE>(hue_noise), noise,
        1.5f, 0.25f);
  }

  In::FrameContext context() const {
    In::FrameContext ctx;
    ctx.frame = 7;
    ctx.time = 7.0f / 30.0f;
    ctx.projection_base = make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
    ctx.palettes = {&palettes[0], &palettes[1], &palettes[2]};
    ctx.hue_rotation_lut = hue_rotation.data();
    ctx.hue_noise_lut = hue_noise.data();
    return ctx;
  }
};

inline const ColorResources &shared_resources() {
  static const ColorResources RESOURCES;
  return RESOURCES;
}

inline constexpr In::ChainEntryRequest DEFAULT_CHAIN[] = {
    {"camera", "sphere.rotate.v2"},
    {"project", "project.stereographic.v2"},
    {"sample", "sample.grid.v2"},
    {"colorize", "colorize.generated-palette.v2"},
};

inline std::array<Vector, 14> sweep_views() {
  std::array<Vector, 14> views = {
      Vector(0, 1, 0),         Vector(0, -1, 0),
      Vector(1, 0, 0),         Vector(-1, 0, 0),
      Vector(0, 0, 1),         Vector(0, 0, -1),
      Vector(1, 1, 1),         Vector(-1, 1, -1),
      Vector(1, -2, 0.5f),     Vector(-0.3f, 0.9f, 0.6f),
      Vector(0.05f, 0.99f, 0), Vector(0.05f, -0.99f, 0),
      Vector(2, 0.1f, -1),     Vector(-1, -1, 2),
  };
  for (Vector &view : views)
    view = view.normalized();
  return views;
}

template <typename T> T &param_as(In::ChainProgram &program, size_t index) {
  return *reinterpret_cast<T *>(program.param_block(index));
}

template <typename T>
const T &state_as(const In::ChainProgram &program, size_t index) {
  return *static_cast<const T *>(program.state_block(index));
}

inline bool color4_identical(const Color4 &a, const Color4 &b) {
  return a.color.r == b.color.r && a.color.g == b.color.g &&
         a.color.b == b.color.b &&
         std::memcmp(&a.alpha, &b.alpha, sizeof(float)) == 0;
}

inline bool float_identical(float a, float b) {
  return std::memcmp(&a, &b, sizeof(float)) == 0;
}

// Member-wise: struct padding bytes are not compared.
inline bool sphere_identical(const PB::SphereSample &a,
                             const PB::SphereSample &b) {
  return float_identical(a.dir.x, b.dir.x) &&
         float_identical(a.dir.y, b.dir.y) &&
         float_identical(a.dir.z, b.dir.z) &&
         float_identical(a.path_length, b.path_length);
}

inline bool plane_identical(const PB::PlaneSample &a,
                            const PB::PlaneSample &b) {
  return float_identical(a.coords.re, b.coords.re) &&
         float_identical(a.coords.im, b.coords.im) &&
         a.provenance.region_id == b.provenance.region_id &&
         a.provenance.component_id == b.provenance.component_id &&
         a.provenance.boundary_flags == b.provenance.boundary_flags &&
         float_identical(a.provenance.fade_edge_distance,
                         b.provenance.fade_edge_distance) &&
         float_identical(a.provenance.value_weight,
                         b.provenance.value_weight) &&
         a.provenance.flags == b.provenance.flags &&
         a.provenance.traits == b.provenance.traits &&
         a.provenance.edge_class == b.provenance.edge_class &&
         float_identical(a.provenance.domain_coverage,
                         b.provenance.domain_coverage) &&
         float_identical(a.sphere.x, b.sphere.x) &&
         float_identical(a.sphere.y, b.sphere.y) &&
         float_identical(a.sphere.z, b.sphere.z) &&
         float_identical(a.path_length, b.path_length);
}

/** Writes every tabled float of a family to its low or high endpoint. */
enum class ValueSet { DEFAULTS, MINIMUMS, MAXIMUMS };

template <typename T> void apply_value_set(T &params, ValueSet set) {
  if (set == ValueSet::DEFAULTS)
    return;
  for (const auto &field : T::FIELDS)
    params.*(field.member) = set == ValueSet::MINIMUMS ? field.min : field.max;
}

// --- template mirror of the slice chain -----------------------------------

struct MirrorFrame {
  Quaternion camera_conjugate;
  Quaternion projection_conjugate;
  In::Op::ProjectChainParams projection;
  In::Op::GridSampleParams sample;
  In::Op::GeneratedPaletteParams color;
  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float oscillation_phase = 0.0f;
  const BakedPalette *palette = nullptr;
  const Pixel *hue_rotation_lut = nullptr;
  const int8_t *hue_noise_lut = nullptr;
};

struct MirrorBinding {
  using FrameState = MirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct MirrorCamera {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const Quaternion &conjugate(const MirrorFrame &frame) {
    return frame.camera_conjugate;
  }
};

struct MirrorProjection {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const Quaternion &conjugate(const MirrorFrame &frame) {
    return frame.projection_conjugate;
  }
  static float pole_fade(const MirrorFrame &frame) {
    return frame.projection.pole_fade;
  }
};

struct MirrorSource {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static const PB::Source::GridSourceParams &params(const MirrorFrame &frame) {
    return frame.sample;
  }
  static PB::Source::PreparedSource prepare(const MirrorFrame &frame) {
    return PB::Source::prepare(frame.source_primary, frame.source_secondary,
                               frame.source_angle);
  }
};

struct MirrorValue {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static float edge_width(const MirrorFrame &frame) {
    return frame.sample.edge_width;
  }
};

template <PB::Color::HueMode HueV, PB::Color::BrightnessEnvelope BrightV>
struct MirrorColor {
  using Binding = MirrorBinding;
  using FrameState = MirrorFrame;
  static PB::Color::PaletteMappingWeights
  mapping_weights(const MirrorFrame &frame) {
    return PB::Color::PaletteMappingWeights::single(
        static_cast<PB::Color::PaletteMapping>(frame.color.mapping_mode));
  }
  static float mapping_frequency(const MirrorFrame &frame) {
    return frame.color.mapping_frequency;
  }
  static float mapping_phase(const MirrorFrame &frame) {
    return frame.color.mapping_phase;
  }
  static float oscillation_depth(const MirrorFrame &frame) {
    return frame.color.phase_oscillation_depth;
  }
  static float oscillation_phase(const MirrorFrame &frame) {
    return frame.oscillation_phase;
  }
  static const BakedPalette &palette(const MirrorFrame &frame) {
    return *frame.palette;
  }
  static PB::Color::HueMode hue_mode(const MirrorFrame &) { return HueV; }
  static float hue_shift_amount(const MirrorFrame &frame) {
    return frame.color.hue_shift_amount;
  }
  static PB::Color::HueRotationLutView hue_rotation(const MirrorFrame &frame) {
    return {frame.hue_rotation_lut, frame.color.hue_shift_amount != 0.0f};
  }
  static PB::Color::HueNoiseLutView hue_noise(const MirrorFrame &frame) {
    return {frame.hue_noise_lut, HueV == PB::Color::HueMode::NOISE &&
                                     frame.color.hue_shift_amount != 0.0f};
  }
  static PB::Color::BrightnessEnvelope
  brightness_envelope(const MirrorFrame &) {
    return BrightV;
  }
  static float brightness_depth(const MirrorFrame &frame) {
    return frame.color.brightness_depth;
  }
  static float opacity_low(const MirrorFrame &frame) {
    return frame.color.opacity_low;
  }
  static float opacity_high(const MirrorFrame &frame) {
    return frame.color.opacity_high;
  }
};

template <typename WeightP, typename CoverageP, PB::Color::HueMode HueV,
          PB::Color::BrightnessEnvelope BrightV>
using MirrorPipeline = PB::Pipeline<
    MirrorBinding, PB::Stage::Rotate<MirrorCamera>,
    PB::Stage::Project<PB::Projection::Stereographic<MirrorProjection>>,
    PB::Stage::Sample<PB::Source::Grid<MirrorSource>, WeightP, CoverageP>,
    PB::Stage::Colorize<
        PB::Color::GeneratedPalette<MirrorColor<HueV, BrightV>>>>;

/** Mirror frame populated from the erased program's blocks: providers read
    the same raw state the erased prepare reads. */
inline MirrorFrame mirror_from(In::ChainProgram &program,
                               const In::FrameContext &ctx) {
  MirrorFrame frame;
  const auto &camera = state_as<In::Op::SpatialWalkState>(program, 0);
  const auto &projection = state_as<In::Op::SpatialWalkState>(program, 1);
  const auto &source = state_as<In::Op::SourceClockState>(program, 2);
  const auto &color = state_as<In::Op::ColorClockState>(program, 3);
  frame.camera_conjugate =
      (make_rotation(Y_AXIS, camera.spin_phase) * camera.wander).conjugate();
  frame.projection_conjugate = (make_rotation(Y_AXIS, projection.spin_phase) *
                                ctx.projection_base * projection.wander)
                                   .conjugate();
  frame.projection = param_as<In::Op::ProjectChainParams>(program, 1);
  frame.sample = param_as<In::Op::GridSampleParams>(program, 2);
  frame.color = param_as<In::Op::GeneratedPaletteParams>(program, 3);
  frame.source_primary = source.primary;
  frame.source_secondary = source.secondary;
  frame.source_angle = source.angle;
  frame.oscillation_phase = color.oscillation_phase;
  frame.palette = ctx.palettes[frame.color.palette_mode];
  frame.hue_rotation_lut = ctx.hue_rotation_lut;
  frame.hue_noise_lut = ctx.hue_noise_lut;
  return frame;
}

template <typename WeightP, typename CoverageP, PB::Color::HueMode HueV,
          PB::Color::BrightnessEnvelope BrightV>
void expect_parity(In::ChainProgram &program, const In::FrameContext &ctx) {
  using Pipe = MirrorPipeline<WeightP, CoverageP, HueV, BrightV>;
  program.prepare(ctx);
  const MirrorFrame mirror = mirror_from(program, ctx);
  const typename Pipe::Frame reference_frame = Pipe::prepare(mirror);
  for (const Vector &view : sweep_views()) {
    const Color4 erased = program.evaluate(view, ctx);
    const Color4 reference = Pipe::shade(view, reference_frame);
    HS_EXPECT_TRUE(color4_identical(erased, reference));
  }
}

/** Compiles the default chain and steps it @p frames times. */
inline void arm_default_chain(In::ChainProgram &program, int frames,
                              ValueSet set) {
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

// --- lifecycle test models ------------------------------------------------

struct CountLifecycle {
  static inline int inits = 0;
  static inline int migrates = 0;
  static inline int destroys = 0;
  static inline bool fail_migrate = false;
  static void reset() {
    inits = 0;
    migrates = 0;
    destroys = 0;
    fail_migrate = false;
  }
};

struct CountParams {
  float gain = 0.5f;
  static constexpr auto FIELDS = std::array{
      PB::Field<CountParams>{"gain", &CountParams::gain, "Gain", 0.0f, 1.0f,
                             PB::FieldCurve::LERP},
  };
};

struct CountingState {
  float accumulator = 0.0f;
  bool live = false;
  CountingState() = default;
  CountingState(const CountingState &) = delete;
  CountingState &operator=(const CountingState &) = delete;
  ~CountingState() {
    if (live)
      ++CountLifecycle::destroys;
  }
};

template <int Tag> struct CountModel {
  static constexpr const char *ID =
      Tag == 0 ? "test.count-a.v2" : "test.count-b.v2";
  static constexpr const char *NAME = "Count";
  using Input = PB::SphereSample;
  using Output = PB::SphereSample;
  using Params = CountParams;
  using State = CountingState;
  struct Prepared {};

  static void init(State &state, In::InstanceId) {
    state.live = true;
    ++CountLifecycle::inits;
  }
  static In::Status migrate(State &dst, const State &src, In::InstanceId) {
    ++CountLifecycle::migrates;
    if (CountLifecycle::fail_migrate)
      return In::Status::FAILED;
    dst.accumulator = src.accumulator;
    dst.live = true;
    return In::Status::OK;
  }
  static void advance(State &state, const Params &) {
    state.accumulator += 1.0f;
  }
  static Prepared prepare(const In::FrameContext &, const Params &,
                          const State &) {
    return {};
  }
  static PB::SphereSample run(const PB::SphereSample &input,
                              const In::FrameContext &, const Params &,
                              const Prepared &) {
    return input;
  }
};

inline constexpr auto FAT_SCHEMA = [] {
  std::array<In::ParamFieldInfo, In::MAX_CHAIN_PARAMS + 1> out{};
  for (auto &field : out)
    field = In::ParamFieldInfo{
        "fat", nullptr, 0.0f, 1.0f,    0.0f,   PB::FieldCurve::LERP,
        false, 0,       0,    nullptr, nullptr};
  return out;
}();

/** Never run; exists to trip the schema-field budget alone. */
inline constexpr In::OperatorDescriptor make_fat_descriptor() {
  In::OperatorDescriptor descriptor =
      In::make_operator_descriptor<CountModel<0>>();
  descriptor.operator_id = "test.fat.v2";
  descriptor.schema = FAT_SCHEMA.data();
  descriptor.schema_count = static_cast<uint16_t>(FAT_SCHEMA.size());
  return descriptor;
}

inline std::span<const In::OperatorDescriptor> extended_table() {
  static constexpr std::array<In::OperatorDescriptor, 7> TABLE{
      In::OPERATOR_TABLE[0],
      In::OPERATOR_TABLE[1],
      In::OPERATOR_TABLE[2],
      In::OPERATOR_TABLE[3],
      In::make_operator_descriptor<CountModel<0>>(),
      In::make_operator_descriptor<CountModel<1>>(),
      make_fat_descriptor(),
  };
  return TABLE;
}

// --- cases ----------------------------------------------------------------

inline void test_shader_chain_table_integrity() {
  static_assert(In::operator_ids_unique());
  static_assert(In::operator_table_monotone());
  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 4u);
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_EXPECT_TRUE(op.operator_id != nullptr && op.display_name != nullptr);
    // Monotonicity pin: proven at table construction, never re-walked.
    HS_EXPECT_LE(static_cast<int>(op.input), static_cast<int>(op.output));
    HS_EXPECT_GT(op.schema_count, 0);
    HS_EXPECT_TRUE(op.runtime.construct_params != nullptr);
    HS_EXPECT_TRUE(op.runtime.init != nullptr);
    HS_EXPECT_TRUE(op.runtime.migrate != nullptr);
    HS_EXPECT_TRUE(op.runtime.destroy != nullptr);
    HS_EXPECT_TRUE(op.runtime.advance != nullptr);
    HS_EXPECT_TRUE(op.runtime.prepare != nullptr);
    HS_EXPECT_TRUE(op.runtime.run != nullptr);
    HS_EXPECT_TRUE(op.runtime.param_address != nullptr);
    HS_EXPECT_GT(op.runtime.param.size, 0u);
    HS_EXPECT_GT(op.runtime.state.size, 0u);
    for (uint16_t index = 0; index < op.schema_count; ++index) {
      const In::ParamFieldInfo &field = op.schema[index];
      HS_EXPECT_TRUE(field.id != nullptr);
      if (field.topology)
        HS_EXPECT_GE(field.enum_count, 2);
      else
        HS_EXPECT_LE(field.min, field.max);
    }
    const bool wants_oracle = op.approximate;
    HS_EXPECT_EQ(op.oracle != PB::ApproximationOracleId::NONE, wants_oracle);
    HS_EXPECT_EQ(op.metric_count > 0, wants_oracle);
  }
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v2") == &In::OPERATOR_TABLE[2]);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v1") == nullptr);
  // The colorize entry owns the chain's only approximation.
  HS_EXPECT_FALSE(In::OPERATOR_TABLE[0].approximate);
  HS_EXPECT_TRUE(In::OPERATOR_TABLE[3].approximate);
  HS_EXPECT_EQ(
      static_cast<int>(In::OPERATOR_TABLE[3].oracle),
      static_cast<int>(PB::ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS));
}

inline void test_shader_chain_schema_and_field_ids() {
  static_assert(PB::field_ids_unique<In::Op::RotateChainParams>());
  static_assert(PB::field_ids_unique<In::Op::ProjectChainParams>());
  static_assert(PB::field_ids_unique<In::Op::GridSampleParams>());
  static_assert(PB::field_ids_unique<In::Op::GeneratedPaletteParams>());
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::Rotate>));
  static_assert(
      In::schema_ids_unique(In::SCHEMA<In::Op::ProjectStereographic>));
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::SampleGrid>));
  static_assert(
      In::schema_ids_unique(In::SCHEMA<In::Op::ColorizeGeneratedPalette>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::SampleGrid>));
  static_assert(
      In::topology_wellformed(In::SCHEMA<In::Op::ColorizeGeneratedPalette>));

  // Schema order is the family table then the topology enum8s; defaults come
  // from the default-constructed family.
  const In::OperatorDescriptor &sample = In::OPERATOR_TABLE[2];
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  HS_EXPECT_EQ(sample.schema_count, GRID_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(sample.schema[0].id) == "pattern-freq");
  HS_EXPECT_TRUE(std::string_view(sample.schema[GRID_FIELDS - 1].id) ==
                 "edge-width");
  const In::ParamFieldInfo &weight = sample.schema[GRID_FIELDS];
  HS_EXPECT_TRUE(std::string_view(weight.id) == "weight-mode");
  HS_EXPECT_TRUE(weight.topology);
  HS_EXPECT_EQ(weight.enum_count, 2);
  HS_EXPECT_EQ(weight.enum_def, 1);
  HS_EXPECT_TRUE(std::string_view(weight.enum_ids[1]) == "projection");
  const In::ParamFieldInfo &coverage = sample.schema[GRID_FIELDS + 1];
  HS_EXPECT_EQ(coverage.enum_count, 4);
  HS_EXPECT_TRUE(std::string_view(coverage.enum_ids[3]) == "edge-fade");

  const In::OperatorDescriptor &colorize = In::OPERATOR_TABLE[3];
  constexpr size_t COLOR_FIELDS = PB::Color::ColorParams::FIELDS.size();
  HS_EXPECT_EQ(colorize.schema_count, COLOR_FIELDS + 4);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS].id) ==
                 "palette-mode");
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 1].id) ==
                 "palette-mapping");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_def, 2);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 3].id) ==
                 "brightness-envelope");

  const In::OperatorDescriptor &rotate = In::OPERATOR_TABLE[0];
  HS_EXPECT_EQ(rotate.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(rotate.schema[0].id) == "wander");
  HS_EXPECT_TRUE(std::string_view(rotate.schema[1].id) == "spin-speed");
  HS_EXPECT_EQ(rotate.schema[1].max, 0.05f);
}

inline void test_shader_chain_instance_id_wellformed() {
  HS_EXPECT_TRUE(In::instance_id_wellformed("camera"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("warp1"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("a"));
  HS_EXPECT_TRUE(In::instance_id_wellformed("wave-shear-2"));
  HS_EXPECT_FALSE(In::instance_id_wellformed(""));
  HS_EXPECT_FALSE(In::instance_id_wellformed("Camera"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("1cam"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam.era"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam--x"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam-"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("-cam"));
  HS_EXPECT_FALSE(In::instance_id_wellformed("cam era"));
  const std::string at_cap(In::MAX_INSTANCE_ID, 'a');
  HS_EXPECT_TRUE(In::instance_id_wellformed(at_cap));
  HS_EXPECT_FALSE(In::instance_id_wellformed(at_cap + "a"));
}

inline void test_shader_chain_slot_and_hash_contract() {
  static_assert(In::SLOT_SIZE == sizeof(PB::PlaneSample));
  static_assert(In::SLOT_ALIGN == alignof(PB::PlaneSample));
  static_assert(In::SLOT_SIZE >= sizeof(Color4));
  HS_EXPECT_EQ(In::SLOT_SIZE, 44u);
  HS_EXPECT_EQ(In::SLOT_ALIGN, 4u);
  HS_EXPECT_EQ(static_cast<int>(In::carrier_id_of<PB::SphereSample>()),
               static_cast<int>(In::CarrierId::SPHERE));
  HS_EXPECT_EQ(static_cast<int>(In::carrier_id_of<Color4>()),
               static_cast<int>(In::CarrierId::COLOR));
  constexpr uint32_t HASH = In::instance_hash("camera", "sphere.rotate.v2");
  static_assert(HASH == In::instance_hash("camera", "sphere.rotate.v2"));
  HS_EXPECT_NE(HASH, In::instance_hash("camera2", "sphere.rotate.v2"));
  HS_EXPECT_NE(HASH, In::instance_hash("camera", "project.stereographic.v2"));
  // The separator keeps ("ab", "c") and ("a", "bc") distinct.
  HS_EXPECT_NE(In::instance_hash("ab", "c"), In::instance_hash("a", "bc"));
}

inline std::string read_file(const char *path) {
  std::string content;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *file = std::fopen(path, "rb");
#pragma clang diagnostic pop
  if (file == nullptr)
    return content;
  char buffer[4096];
  size_t bytes;
  while ((bytes = std::fread(buffer, 1, sizeof(buffer), file)) > 0)
    content.append(buffer, bytes);
  std::fclose(file);
  return content;
}

inline void test_shader_chain_catalog_golden() {
  std::string catalog;
  In::append_catalog_json(catalog);
  catalog += '\n';
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  const bool regen = std::getenv("HS_SHADER_CHAIN_CATALOG_REGEN") != nullptr;
#pragma clang diagnostic pop
  if (regen) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    std::FILE *file = std::fopen(HS_SHADER_CHAIN_CATALOG_PATH, "wb");
#pragma clang diagnostic pop
    HS_EXPECT_TRUE(file != nullptr);
    if (file != nullptr) {
      std::fwrite(catalog.data(), 1, catalog.size(), file);
      std::fclose(file);
    }
  }
  const std::string golden = read_file(HS_SHADER_CHAIN_CATALOG_PATH);
  HS_EXPECT_FALSE(golden.empty());
  HS_EXPECT_TRUE(catalog == golden);
  if (catalog != golden)
    std::printf("shader_chain catalog drift: %zu generated vs %zu golden "
                "bytes; regen with HS_SHADER_CHAIN_CATALOG_REGEN=1\n",
                catalog.size(), golden.size());
}

inline void test_shader_chain_catalog_shape() {
  std::string catalog;
  In::append_catalog_json(catalog);
  const auto contains = [&catalog](const char *needle) {
    return catalog.find(needle) != std::string::npos;
  };
  HS_EXPECT_TRUE(contains("\"catalog_version\":2"));
  HS_EXPECT_TRUE(contains("\"budgets\":{\"max_chain_ops\":32,"
                          "\"arena_bytes\":98304,\"max_params\":224,"
                          "\"max_instance_id_length\":48,"
                          "\"per_op_overhead_bytes\":49}"));
  HS_EXPECT_TRUE(
      contains("\"carriers\":[\"sphere\",\"plane\",\"field\",\"color\"]"));
  HS_EXPECT_TRUE(contains("\"id\":\"sphere.rotate.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"project.stereographic.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.grid.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"colorize.generated-palette.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"weight-mode\",\"topology\":true,"
                          "\"values\":[\"none\",\"projection\"],"
                          "\"default\":\"projection\""));
  HS_EXPECT_TRUE(contains("\"values\":[\"none\",\"weight\",\"weight-squared\","
                          "\"edge-fade\"]"));
  HS_EXPECT_TRUE(contains("\"approximate\":true"));
  // The cataloged block layouts are the real ABI figures.
  std::string sample_blocks = "\"blocks\":{\"param\":{\"size\":";
  sample_blocks += std::to_string(sizeof(In::Op::GridSampleParams));
  HS_EXPECT_TRUE(contains(sample_blocks.c_str()));
}

inline void test_shader_chain_default_chain_renders() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  HS_EXPECT_FALSE(program.compiled());
  arm_default_chain(program, 3, ValueSet::DEFAULTS);
  HS_EXPECT_TRUE(program.compiled());
  HS_EXPECT_EQ(program.ops().size(), 4u);
  HS_EXPECT_TRUE(std::string_view(program.ops()[0].instance) == "camera");
  HS_EXPECT_TRUE(std::string_view(program.ops()[3].instance) == "colorize");
  HS_EXPECT_EQ(program.ops()[0].stable_hash,
               In::instance_hash("camera", "sphere.rotate.v2"));
  HS_EXPECT_GT(program.used_bytes(), 0u);
  HS_EXPECT_LE(program.used_bytes(), In::CHAIN_ARENA_BYTES);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  bool any_opaque = false;
  for (const Vector &view : sweep_views()) {
    const Color4 color = program.evaluate(view, ctx);
    HS_EXPECT_TRUE(std::isfinite(color.alpha));
    HS_EXPECT_GE(color.alpha, 0.0f);
    if (color.alpha > 0.0f &&
        (color.color.r | color.color.g | color.color.b) != 0)
      any_opaque = true;
  }
  HS_EXPECT_TRUE(any_opaque);
  program.clear();
  HS_EXPECT_FALSE(program.compiled());
}

inline void test_shader_chain_param_address_channel() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 1, ValueSet::DEFAULTS);
  const In::OperatorDescriptor &sample_op = *program.ops()[2].op;
  uint8_t *block = program.param_block(2);
  auto &typed = *reinterpret_cast<In::Op::GridSampleParams *>(block);
  // Every schema index resolves inside the block, floats first.
  for (uint16_t index = 0; index < sample_op.schema_count; ++index) {
    const uint8_t *address = static_cast<const uint8_t *>(
        sample_op.runtime.param_address(block, index));
    HS_EXPECT_TRUE(address >= block &&
                   address < block + sample_op.runtime.param.size);
  }
  auto *freq = static_cast<float *>(sample_op.runtime.param_address(block, 0));
  *freq = 5.5f;
  HS_EXPECT_EQ(typed.pattern_freq, 5.5f);
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  auto *coverage_mode = static_cast<uint8_t *>(
      sample_op.runtime.param_address(block, GRID_FIELDS + 1));
  *coverage_mode = static_cast<uint8_t>(In::Op::CoverageMode::EDGE_FADE);
  HS_EXPECT_EQ(typed.coverage_mode,
               static_cast<uint8_t>(In::Op::CoverageMode::EDGE_FADE));
  // The write is visible to the render: identical view, different coverage.
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  const Vector view = Vector(1, 1, 1).normalized();
  const Color4 faded = program.evaluate(view, ctx);
  *coverage_mode = static_cast<uint8_t>(In::Op::CoverageMode::NONE);
  program.prepare(ctx);
  const Color4 full = program.evaluate(view, ctx);
  HS_EXPECT_GE(full.alpha, faded.alpha);
  program.clear();
}

inline void test_shader_chain_parity_rotate_project() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 5, set);
    const In::FrameContext ctx = shared_resources().context();
    program.prepare(ctx);
    // The erased prepared conjugates equal the same composition recomputed
    // from the state blocks.
    const MirrorFrame mirror = mirror_from(program, ctx);
    const auto &camera_prepared =
        *reinterpret_cast<const In::Op::Rotate::Prepared *>(
            program.prepared_block(0));
    const auto &project_prepared =
        *reinterpret_cast<const In::Op::ProjectStereographic::Prepared *>(
            program.prepared_block(1));
    HS_EXPECT_EQ(std::memcmp(&camera_prepared.conjugate,
                             &mirror.camera_conjugate, sizeof(Quaternion)),
                 0);
    HS_EXPECT_EQ(std::memcmp(&project_prepared.conjugate,
                             &mirror.projection_conjugate, sizeof(Quaternion)),
                 0);
    // Op-level kernel parity for the two sphere-family ops.
    using BoundRotate = PB::Stage::Rotate<MirrorCamera>::Bind<MirrorBinding>;
    using BoundProject = PB::Stage::Project<
        PB::Projection::Stereographic<MirrorProjection>>::Bind<MirrorBinding>;
    const In::OperatorDescriptor &rotate_op = *program.ops()[0].op;
    const In::OperatorDescriptor &project_op = *program.ops()[1].op;
    for (const Vector &view : sweep_views()) {
      const PB::SphereSample seed{view, 0.0f};
      alignas(In::SLOT_ALIGN) uint8_t erased_out[In::SLOT_SIZE];
      rotate_op.runtime.run(&seed, erased_out, ctx, program.param_block(0),
                            program.prepared_block(0));
      const auto &erased_rotated =
          *std::launder(reinterpret_cast<PB::SphereSample *>(erased_out));
      const PB::SphereSample reference_rotated =
          BoundRotate::run(seed, mirror, {});
      HS_EXPECT_TRUE(sphere_identical(erased_rotated, reference_rotated));
      alignas(In::SLOT_ALIGN) uint8_t projected_out[In::SLOT_SIZE];
      project_op.runtime.run(&erased_rotated, projected_out, ctx,
                             program.param_block(1), program.prepared_block(1));
      const auto &erased_projected =
          *std::launder(reinterpret_cast<PB::PlaneSample *>(projected_out));
      const PB::PlaneSample reference_projected =
          BoundProject::run(reference_rotated, mirror, {});
      HS_EXPECT_TRUE(plane_identical(erased_projected, reference_projected));
    }
    program.clear();
  }
}

inline void run_sample_variant(In::ChainProgram &program,
                               const In::FrameContext &ctx,
                               In::Op::WeightMode weight,
                               In::Op::CoverageMode coverage) {
  auto &params = param_as<In::Op::GridSampleParams>(program, 2);
  params.weight_mode = static_cast<uint8_t>(weight);
  params.coverage_mode = static_cast<uint8_t>(coverage);
  constexpr auto HUE = PB::Color::HueMode::NOISE;
  constexpr auto ENV = PB::Color::BrightnessEnvelope::NONE;
  if (weight == In::Op::WeightMode::NONE) {
    switch (coverage) {
    case In::Op::CoverageMode::NONE:
      expect_parity<PB::Weight::None, PB::ProjectedCoverage::None, HUE, ENV>(
          program, ctx);
      break;
    case In::Op::CoverageMode::WEIGHT:
      expect_parity<PB::Weight::None, PB::ProjectedCoverage::Weight, HUE, ENV>(
          program, ctx);
      break;
    case In::Op::CoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::None, PB::ProjectedCoverage::WeightSquared, HUE,
                    ENV>(program, ctx);
      break;
    case In::Op::CoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::None,
                    PB::ProjectedCoverage::EdgeFade<MirrorValue>, HUE, ENV>(
          program, ctx);
      break;
    }
  } else {
    switch (coverage) {
    case In::Op::CoverageMode::NONE:
      expect_parity<PB::Weight::Projection, PB::ProjectedCoverage::None, HUE,
                    ENV>(program, ctx);
      break;
    case In::Op::CoverageMode::WEIGHT:
      expect_parity<PB::Weight::Projection, PB::ProjectedCoverage::Weight, HUE,
                    ENV>(program, ctx);
      break;
    case In::Op::CoverageMode::WEIGHT_SQUARED:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectedCoverage::WeightSquared, HUE, ENV>(program,
                                                                    ctx);
      break;
    case In::Op::CoverageMode::EDGE_FADE:
      expect_parity<PB::Weight::Projection,
                    PB::ProjectedCoverage::EdgeFade<MirrorValue>, HUE, ENV>(
          program, ctx);
      break;
    }
  }
}

inline void test_shader_chain_parity_sample_variants() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    for (const In::Op::WeightMode weight :
         {In::Op::WeightMode::NONE, In::Op::WeightMode::PROJECTION})
      for (const In::Op::CoverageMode coverage :
           {In::Op::CoverageMode::NONE, In::Op::CoverageMode::WEIGHT,
            In::Op::CoverageMode::WEIGHT_SQUARED,
            In::Op::CoverageMode::EDGE_FADE})
        run_sample_variant(program, ctx, weight, coverage);
    program.clear();
  }
}

inline void run_colorize_variant(In::ChainProgram &program,
                                 const In::FrameContext &ctx,
                                 In::Op::HueShiftMode hue,
                                 In::Op::EnvelopeMode envelope) {
  auto &params = param_as<In::Op::GeneratedPaletteParams>(program, 3);
  params.hue_mode = static_cast<uint8_t>(hue);
  params.envelope_mode = static_cast<uint8_t>(envelope);
  for (uint8_t palette_mode = 0; palette_mode < 3; ++palette_mode) {
    params.palette_mode = palette_mode;
    for (uint8_t mapping = 0; mapping < 4; ++mapping) {
      params.mapping_mode = mapping;
      using WeightP = PB::Weight::Projection;
      using CoverageP = PB::ProjectedCoverage::Weight;
      if (hue == In::Op::HueShiftMode::NOISE) {
        if (envelope == In::Op::EnvelopeMode::NONE)
          expect_parity<WeightP, CoverageP, PB::Color::HueMode::NOISE,
                        PB::Color::BrightnessEnvelope::NONE>(program, ctx);
        else
          expect_parity<WeightP, CoverageP, PB::Color::HueMode::NOISE,
                        PB::Color::BrightnessEnvelope::CUP>(program, ctx);
      } else {
        if (envelope == In::Op::EnvelopeMode::NONE)
          expect_parity<WeightP, CoverageP, PB::Color::HueMode::PATH_LENGTH,
                        PB::Color::BrightnessEnvelope::NONE>(program, ctx);
        else
          expect_parity<WeightP, CoverageP, PB::Color::HueMode::PATH_LENGTH,
                        PB::Color::BrightnessEnvelope::CUP>(program, ctx);
      }
    }
  }
}

inline void test_shader_chain_parity_colorize_variants() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_default_chain(program, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    for (const In::Op::HueShiftMode hue :
         {In::Op::HueShiftMode::NOISE, In::Op::HueShiftMode::PATH_LENGTH})
      for (const In::Op::EnvelopeMode envelope :
           {In::Op::EnvelopeMode::NONE, In::Op::EnvelopeMode::CUP})
        run_colorize_variant(program, ctx, hue, envelope);
    program.clear();
  }
}

/** Renders the sweep into @p out for a byte-identity comparison. */
inline void snapshot_render(In::ChainProgram &program,
                            const In::FrameContext &ctx,
                            std::array<Color4, 14> &out) {
  const auto views = sweep_views();
  for (size_t index = 0; index < views.size(); ++index)
    out[index] = program.evaluate(views[index], ctx);
}

inline void expect_refusal(In::ChainProgram &program,
                           std::span<const In::ChainEntryRequest> request,
                           In::ChainStatus expected, int16_t expected_index,
                           const In::FrameContext &ctx,
                           const std::array<Color4, 14> &baseline) {
  const In::ChainRefusal refusal = program.compile(request);
  HS_EXPECT_EQ(static_cast<int>(refusal.code), static_cast<int>(expected));
  HS_EXPECT_EQ(refusal.entry_index, expected_index);
  // Transactional: the live program, its params, prepared state and render
  // are all byte-identical after the refusal.
  HS_EXPECT_EQ(program.ops().size(), 4u);
  std::array<Color4, 14> after;
  snapshot_render(program, ctx, after);
  for (size_t index = 0; index < after.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after[index], baseline[index]));
}

inline void test_shader_chain_refusal_shape() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 2, ValueSet::MAXIMUMS);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  std::array<Color4, 14> baseline;
  snapshot_render(program, ctx, baseline);

  expect_refusal(program, {}, In::ChainStatus::EMPTY, -1, ctx, baseline);

  std::array<In::ChainEntryRequest, In::MAX_CHAIN_OPS + 1> too_long;
  std::array<std::string, In::MAX_CHAIN_OPS + 1> labels;
  for (size_t index = 0; index < too_long.size(); ++index) {
    labels[index] = "cam" + std::to_string(index);
    too_long[index] = {labels[index], "sphere.rotate.v2"};
  }
  expect_refusal(program, too_long, In::ChainStatus::TOO_LONG, -1, ctx,
                 baseline);

  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.mirror-tile.v2"},
  };
  expect_refusal(program, unknown, In::ChainStatus::UNKNOWN_OPERATOR, 1, ctx,
                 baseline);

  const In::ChainEntryRequest duplicate[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera", "project.stereographic.v2"},
  };
  expect_refusal(program, duplicate, In::ChainStatus::DUPLICATE_INSTANCE, 1,
                 ctx, baseline);

  const In::ChainEntryRequest malformed[] = {
      {"camera", "sphere.rotate.v2"},
      {"Bad.Label", "project.stereographic.v2"},
  };
  expect_refusal(program, malformed, In::ChainStatus::MALFORMED_INSTANCE, 1,
                 ctx, baseline);

  const In::ChainEntryRequest bad_entry[] = {
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  expect_refusal(program, bad_entry, In::ChainStatus::ENTRY_FAMILY, 0, ctx,
                 baseline);

  const In::ChainEntryRequest bad_exit[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
  };
  expect_refusal(program, bad_exit, In::ChainStatus::EXIT_FAMILY, 2, ctx,
                 baseline);

  const In::ChainEntryRequest mismatch[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  expect_refusal(program, mismatch, In::ChainStatus::CARRIER_MISMATCH, 2, ctx,
                 baseline);

  // State continuity: the shape refusals above left every clock untouched.
  const auto &source = state_as<In::Op::SourceClockState>(program, 2);
  const float speed = param_as<In::Op::GridSampleParams>(program, 2).speed;
  HS_EXPECT_GT(speed, 0.0f);
  HS_EXPECT_EQ(source.primary, fmodf(fmodf(speed, TWO_PI_F) + speed, TWO_PI_F));
  program.clear();
}

inline void test_shader_chain_refusal_budget_overflows() {
  // Exact-fit boundary: capacity == used commits, capacity - 1 refuses.
  auto measured = std::make_unique<ProgramFixture>();
  arm_default_chain(measured->program, 0, ValueSet::DEFAULTS);
  const size_t needed = measured->program.used_bytes();
  measured->program.clear();

  auto exact = std::make_unique<ProgramFixture>(needed);
  const In::ChainRefusal fits = exact->program.compile(
      std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(fits.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(exact->program.used_bytes(), needed);

  auto small = std::make_unique<ProgramFixture>(needed - 1);
  const In::ChainRefusal overflow = small->program.compile(
      std::span<const In::ChainEntryRequest>(DEFAULT_CHAIN));
  HS_EXPECT_EQ(static_cast<int>(overflow.code),
               static_cast<int>(In::ChainStatus::ARENA_OVERFLOW));
  HS_EXPECT_EQ(overflow.entry_index, -1);
  HS_EXPECT_FALSE(small->program.compiled());

  // A committed program survives a later over-budget recompile.
  const In::ChainEntryRequest longer[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal rejected = exact->program.compile(longer);
  HS_EXPECT_EQ(static_cast<int>(rejected.code),
               static_cast<int>(In::ChainStatus::ARENA_OVERFLOW));
  HS_EXPECT_EQ(exact->program.ops().size(), 4u);
  const In::FrameContext ctx = shared_resources().context();
  exact->program.prepare(ctx);
  const Color4 color =
      exact->program.evaluate(Vector(1, 1, 1).normalized(), ctx);
  HS_EXPECT_TRUE(std::isfinite(color.alpha));
  exact->program.clear();

  // Schema-field budget: one fat operator overflows MAX_CHAIN_PARAMS.
  CountLifecycle::reset();
  auto fat =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  const In::ChainEntryRequest fat_chain[] = {
      {"fat", "test.fat.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal params_overflow = fat->program.compile(fat_chain);
  HS_EXPECT_EQ(static_cast<int>(params_overflow.code),
               static_cast<int>(In::ChainStatus::PARAM_OVERFLOW));
  HS_EXPECT_EQ(params_overflow.entry_index, -1);
  HS_EXPECT_FALSE(fat->program.compiled());
  // Refused before layout: no lifecycle callback ran.
  HS_EXPECT_EQ(CountLifecycle::inits, 0);
}

inline void test_shader_chain_refusal_migrate_failed() {
  CountLifecycle::reset();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"counter", "test.count-a.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal first = program.compile(chain);
  HS_EXPECT_EQ(static_cast<int>(first.code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::inits, 1);
  program.advance();
  program.advance();
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  const In::FrameContext ctx = shared_resources().context();
  program.prepare(ctx);
  std::array<Color4, 14> baseline;
  snapshot_render(program, ctx, baseline);

  CountLifecycle::fail_migrate = true;
  const int destroys_before = CountLifecycle::destroys;
  const In::ChainRefusal failed = program.compile(chain);
  HS_EXPECT_EQ(static_cast<int>(failed.code),
               static_cast<int>(In::ChainStatus::MIGRATE_FAILED));
  HS_EXPECT_EQ(failed.entry_index, 0);
  // The candidate tore down as a unit and the live program is untouched: the
  // failed dst was rolled back unconstructed, no live state was destroyed,
  // and the accumulated phase survives.
  HS_EXPECT_EQ(CountLifecycle::destroys, destroys_before);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  std::array<Color4, 14> after;
  snapshot_render(program, ctx, after);
  for (size_t index = 0; index < after.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after[index], baseline[index]));

  // A failing migrate deeper in the chain destroys the candidate states
  // constructed before it — and only those. "fresh" is a new pair (init at
  // entry 0), "counter" a surviving pair whose migrate fails at entry 1.
  const In::ChainEntryRequest deep_chain[] = {
      {"fresh", "test.count-b.v2"},
      {"counter", "test.count-a.v2"},
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const int inits_before = CountLifecycle::inits;
  const int migrates_before = CountLifecycle::migrates;
  const int deep_destroys_before = CountLifecycle::destroys;
  const In::ChainRefusal deep = program.compile(deep_chain);
  HS_EXPECT_EQ(static_cast<int>(deep.code),
               static_cast<int>(In::ChainStatus::MIGRATE_FAILED));
  HS_EXPECT_EQ(deep.entry_index, 1);
  // The fresh candidate was constructed, then torn down with the candidate
  // arena; the live program keeps its state.
  HS_EXPECT_EQ(CountLifecycle::inits, inits_before + 1);
  HS_EXPECT_EQ(CountLifecycle::migrates, migrates_before + 1);
  HS_EXPECT_EQ(CountLifecycle::destroys, deep_destroys_before + 1);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 2.0f);
  std::array<Color4, 14> after_deep;
  snapshot_render(program, ctx, after_deep);
  for (size_t index = 0; index < after_deep.size(); ++index)
    HS_EXPECT_TRUE(color4_identical(after_deep[index], baseline[index]));
  CountLifecycle::fail_migrate = false;
  program.clear();
}

inline void test_shader_chain_state_identity_migration() {
  CountLifecycle::reset();
  auto fixture =
      std::make_unique<ProgramFixture>(In::CHAIN_ARENA_BYTES, extended_table());
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest first[] = {
      {"a", "test.count-a.v2"},
      {"b", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  HS_EXPECT_EQ(static_cast<int>(program.compile(first).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::inits, 2);
  program.advance();
  program.advance();
  program.advance();

  // Surviving pair migrates and keeps its phase; the removed instance is
  // destroyed; the new label gets a fresh init.
  const In::ChainEntryRequest second[] = {
      {"a", "test.count-a.v2"},
      {"c", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const int destroys_before = CountLifecycle::destroys;
  HS_EXPECT_EQ(static_cast<int>(program.compile(second).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(CountLifecycle::migrates, 1);
  HS_EXPECT_EQ(CountLifecycle::inits, 3);
  // Commit destroyed the loser arena: old "a" and old "b".
  HS_EXPECT_EQ(CountLifecycle::destroys, destroys_before + 2);
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 3.0f);
  HS_EXPECT_EQ(state_as<CountingState>(program, 1).accumulator, 0.0f);

  // Same label, different operator: a fresh init, never a migration.
  const In::ChainEntryRequest third[] = {
      {"a", "test.count-b.v2"},
      {"c", "test.count-a.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const int migrates_before = CountLifecycle::migrates;
  HS_EXPECT_EQ(static_cast<int>(program.compile(third).code),
               static_cast<int>(In::ChainStatus::OK));
  HS_EXPECT_EQ(state_as<CountingState>(program, 0).accumulator, 0.0f);
  // "a" re-inits under the new operator; "c" migrates.
  HS_EXPECT_EQ(CountLifecycle::inits, 4);
  HS_EXPECT_EQ(CountLifecycle::migrates, migrates_before + 1);
  program.clear();
  // Every live construction (init or successful migrate clone) was destroyed.
  HS_EXPECT_EQ(CountLifecycle::destroys,
               CountLifecycle::inits + CountLifecycle::migrates);
}

inline void test_shader_chain_state_continuity_slice() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 6, ValueSet::MAXIMUMS);
  const auto camera_before = state_as<In::Op::SpatialWalkState>(program, 0);
  const auto source_before = state_as<In::Op::SourceClockState>(program, 2);
  HS_EXPECT_GT(source_before.primary, 0.0f);

  const In::ChainEntryRequest edited[] = {
      {"camera", "sphere.rotate.v2"},
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  HS_EXPECT_EQ(static_cast<int>(program.compile(edited).code),
               static_cast<int>(In::ChainStatus::OK));
  // An unchanged pair keeps its accumulated walk and clocks across an edit
  // elsewhere; the new instance starts fresh.
  const auto &camera_after = state_as<In::Op::SpatialWalkState>(program, 0);
  HS_EXPECT_EQ(std::memcmp(&camera_after.wander, &camera_before.wander,
                           sizeof(Quaternion)),
               0);
  HS_EXPECT_EQ(camera_after.spin_phase, camera_before.spin_phase);
  HS_EXPECT_EQ(camera_after.walk_time, camera_before.walk_time);
  const auto &fresh = state_as<In::Op::SpatialWalkState>(program, 1);
  HS_EXPECT_EQ(fresh.spin_phase, 0.0f);
  HS_EXPECT_EQ(fresh.walk_time, 0.0f);
  const auto &source_after = state_as<In::Op::SourceClockState>(program, 3);
  HS_EXPECT_EQ(source_after.primary, source_before.primary);
  HS_EXPECT_EQ(source_after.secondary, source_before.secondary);
  HS_EXPECT_EQ(source_after.angle, source_before.angle);
  program.clear();
}

inline void test_shader_chain_determinism() {
  auto first = std::make_unique<ProgramFixture>();
  auto second = std::make_unique<ProgramFixture>();
  arm_default_chain(first->program, 5, ValueSet::MAXIMUMS);
  arm_default_chain(second->program, 5, ValueSet::MAXIMUMS);
  const In::FrameContext ctx = shared_resources().context();
  first->program.prepare(ctx);
  second->program.prepare(ctx);
  for (const Vector &view : sweep_views()) {
    const Color4 a = first->program.evaluate(view, ctx);
    const Color4 b = second->program.evaluate(view, ctx);
    HS_EXPECT_TRUE(color4_identical(a, b));
  }
  // The walk seed is the instance hash: a different label walks differently.
  auto relabeled = std::make_unique<ProgramFixture>();
  const In::ChainEntryRequest renamed[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  HS_EXPECT_EQ(static_cast<int>(relabeled->program.compile(renamed).code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(relabeled->program, 0),
                  ValueSet::MAXIMUMS);
  for (int frame = 0; frame < 5; ++frame)
    relabeled->program.advance();
  const auto &walk_a = state_as<In::Op::SpatialWalkState>(first->program, 0);
  const auto &walk_b =
      state_as<In::Op::SpatialWalkState>(relabeled->program, 0);
  HS_EXPECT_NE(std::memcmp(&walk_a.wander, &walk_b.wander, sizeof(Quaternion)),
               0);
  first->program.clear();
  second->program.clear();
  relabeled->program.clear();
}

inline int run_shader_chain_tests() {
  ModuleFixture fixture("shader_chain");
  test_shader_chain_table_integrity();
  test_shader_chain_schema_and_field_ids();
  test_shader_chain_instance_id_wellformed();
  test_shader_chain_slot_and_hash_contract();
  test_shader_chain_catalog_golden();
  test_shader_chain_catalog_shape();
  test_shader_chain_default_chain_renders();
  test_shader_chain_param_address_channel();
  test_shader_chain_parity_rotate_project();
  test_shader_chain_parity_sample_variants();
  test_shader_chain_parity_colorize_variants();
  test_shader_chain_refusal_shape();
  test_shader_chain_refusal_budget_overflows();
  test_shader_chain_refusal_migrate_failed();
  test_shader_chain_state_identity_migration();
  test_shader_chain_state_continuity_slice();
  test_shader_chain_determinism();
  return fixture.result();
}

} // namespace shader_chain_tests
} // namespace hs_test
