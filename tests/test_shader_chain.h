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
 *
 * The golden's block sizes and alignments are the native ABI, emitted by the
 * host build this suite runs in. scripts/engine_catalog.json is a second
 * catalog stating the wasm32 ABI the browser workbench budgets against. The
 * two disagree in every prepared block by construction: they are not to be
 * reconciled, and copying either over the other retargets a consumer's budget
 * math.
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
#include "workbench/ShaderChain.h"

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

inline bool field_identical(const PB::FieldSample &a,
                            const PB::FieldSample &b) {
  return float_identical(a.value, b.value) &&
         float_identical(a.coverage, b.coverage) &&
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
  static float singularity_fade(const MirrorFrame &frame) {
    return frame.projection.singularity_fade;
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

consteval In::OperatorDescriptor table_entry(std::string_view operator_id) {
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE)
    if (operator_id == op.operator_id)
      return op;
  throw "unknown operator id";
}

inline std::span<const In::OperatorDescriptor> extended_table() {
  static constexpr std::array<In::OperatorDescriptor, 7> TABLE{
      table_entry("sphere.rotate.v2"),
      table_entry("project.stereographic.v2"),
      table_entry("sample.grid.v2"),
      table_entry("colorize.generated-palette.v2"),
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
  HS_EXPECT_EQ(In::OPERATOR_TABLE.size(), 37u);
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE) {
    HS_EXPECT_TRUE(op.operator_id != nullptr && op.display_name != nullptr);
    // Monotonicity pin: proven at table construction, never re-walked.
    HS_EXPECT_LE(static_cast<int>(op.input), static_cast<int>(op.output));
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
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sample.grid.v1") == nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.curl.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.displace.ripple.v2") != nullptr);
  HS_EXPECT_TRUE(In::find_operator("sphere.lens.kaleidoscope.v2") != nullptr);
  // Exactly two approximations ship: the colorizer's LUT path and the fast
  // square Peirce projection, each with its own oracle.
  HS_EXPECT_FALSE(In::find_operator("sphere.rotate.v2")->approximate);
  HS_EXPECT_FALSE(In::find_operator("project.peirce.v2")->approximate);
  const In::OperatorDescriptor &colorize =
      *In::find_operator("colorize.generated-palette.v2");
  HS_EXPECT_TRUE(colorize.approximate);
  HS_EXPECT_EQ(
      static_cast<int>(colorize.oracle),
      static_cast<int>(PB::ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS));
  const In::OperatorDescriptor &peirce_fast =
      *In::find_operator("project.peirce-square-fast.v2");
  HS_EXPECT_TRUE(peirce_fast.approximate);
  HS_EXPECT_EQ(static_cast<int>(peirce_fast.oracle),
               static_cast<int>(PB::ApproximationOracleId::PEIRCE_FAST_SQUARE));
  HS_EXPECT_EQ(peirce_fast.metric_count, 3);
  size_t approximate_count = 0;
  for (const In::OperatorDescriptor &op : In::OPERATOR_TABLE)
    approximate_count += op.approximate ? 1 : 0;
  HS_EXPECT_EQ(approximate_count, 2u);
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

  static_assert(PB::field_ids_unique<In::Op::CurlDisplaceParams>());
  static_assert(PB::field_ids_unique<In::Op::DirectDisplaceParams>());
  static_assert(PB::field_ids_unique<PB::Surface::PeriodicRippleParams>());
  static_assert(PB::field_ids_unique<In::Op::MobiusChainParams>());
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::schema_ids_unique(In::SCHEMA<In::Op::DisplaceRipple>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::DisplaceCurl>));
  static_assert(In::topology_wellformed(In::SCHEMA<In::Op::LensKaleidoscope>));

  // Schema order is the family table then the topology enum8s; defaults come
  // from the default-constructed family.
  const In::OperatorDescriptor &sample = *In::find_operator("sample.grid.v2");
  constexpr size_t GRID_FIELDS = In::Op::GridSampleParams::FIELDS.size();
  HS_EXPECT_EQ(sample.schema_count, GRID_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(sample.schema[0].id) == "pattern-freq");
  HS_EXPECT_EQ(sample.schema[0].max, 64.0f);
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

  const In::OperatorDescriptor &colorize =
      *In::find_operator("colorize.generated-palette.v2");
  constexpr size_t COLOR_FIELDS = PB::Color::ColorParams::FIELDS.size();
  HS_EXPECT_EQ(colorize.schema_count, COLOR_FIELDS + 4);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS].id) ==
                 "palette-mode");
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 1].id) ==
                 "palette-mapping");
  HS_EXPECT_EQ(colorize.schema[COLOR_FIELDS + 1].enum_def, 2);
  HS_EXPECT_TRUE(std::string_view(colorize.schema[COLOR_FIELDS + 3].id) ==
                 "brightness-envelope");

  const In::OperatorDescriptor &rotate = *In::find_operator("sphere.rotate.v2");
  HS_EXPECT_EQ(rotate.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(rotate.schema[0].id) == "wander");
  HS_EXPECT_TRUE(std::string_view(rotate.schema[1].id) == "spin-speed");
  HS_EXPECT_EQ(rotate.schema[1].max, 0.05f);

  // Sphere batch: the parameterless lenses register empty schemas; the
  // kaleidoscope is topology-only; the mobius family carries the v1 ids.
  HS_EXPECT_EQ(In::find_operator("sphere.lens.glitch.v2")->schema_count, 0);
  HS_EXPECT_EQ(In::find_operator("sphere.lens.twist.v2")->schema_count, 0);
  const In::OperatorDescriptor &kaleidoscope =
      *In::find_operator("sphere.lens.kaleidoscope.v2");
  HS_EXPECT_EQ(kaleidoscope.schema_count, 1);
  HS_EXPECT_TRUE(kaleidoscope.schema[0].topology);
  HS_EXPECT_EQ(kaleidoscope.schema[0].enum_count, 9);
  HS_EXPECT_TRUE(std::string_view(kaleidoscope.schema[0].enum_ids[0]) ==
                 "azimuthal");
  HS_EXPECT_TRUE(std::string_view(kaleidoscope.schema[0].enum_ids[8]) ==
                 "octagonal-prism");
  const In::OperatorDescriptor &mobius =
      *In::find_operator("sphere.lens.mobius.v2");
  HS_EXPECT_EQ(mobius.schema_count, 8);
  HS_EXPECT_TRUE(std::string_view(mobius.schema[0].id) == "mobius-a-re");
  HS_EXPECT_TRUE(std::string_view(mobius.schema[7].id) == "mobius-d-im");
  HS_EXPECT_EQ(mobius.schema[0].min, -4.0f);
  HS_EXPECT_EQ(mobius.schema[0].max, 4.0f);
  HS_EXPECT_EQ(mobius.schema[0].def, 0.7071067811865475f);
  const In::OperatorDescriptor &curl =
      *In::find_operator("sphere.displace.curl.v2");
  constexpr size_t CURL_FIELDS = In::Op::CurlDisplaceParams::FIELDS.size();
  HS_EXPECT_EQ(curl.schema_count, CURL_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS].id) == "basis");
  HS_EXPECT_EQ(curl.schema[CURL_FIELDS].enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS + 1].id) ==
                 "integrator");
  HS_EXPECT_TRUE(std::string_view(curl.schema[CURL_FIELDS + 1].enum_ids[2]) ==
                 "midpoint-2x");
  const In::OperatorDescriptor &direct =
      *In::find_operator("sphere.displace.direct.v2");
  constexpr size_t DIRECT_FIELDS = In::Op::DirectDisplaceParams::FIELDS.size();
  HS_EXPECT_EQ(direct.schema_count, DIRECT_FIELDS + 1);
  HS_EXPECT_TRUE(std::string_view(direct.schema[DIRECT_FIELDS].id) == "basis");
  const In::OperatorDescriptor &ripple =
      *In::find_operator("sphere.displace.ripple.v2");
  HS_EXPECT_EQ(ripple.schema_count,
               PB::Surface::PeriodicRippleParams::FIELDS.size());
  HS_EXPECT_TRUE(std::string_view(ripple.schema[0].id) == "period");
  HS_EXPECT_EQ(ripple.schema[0].min, 30.0f);
  HS_EXPECT_EQ(ripple.schema[0].max, 143.0f);
  HS_EXPECT_EQ(ripple.schema[0].def, 80.0f);
  HS_EXPECT_EQ(ripple.schema[1].max, 0.15f);
  HS_EXPECT_EQ(ripple.schema[1].def, 0.15f);
  HS_EXPECT_EQ(ripple.schema[2].max, 5.0f);
  HS_EXPECT_EQ(ripple.schema[2].def, 0.1f);
  HS_EXPECT_EQ(ripple.schema[3].min, 0.7f);
  HS_EXPECT_EQ(ripple.schema[3].def, 0.7f);
  HS_EXPECT_TRUE(std::string_view(ripple.schema[5].id) == "center-polar");
  In::Op::RipplePhaseState ripple_state;
  PB::Surface::PeriodicRippleParams ripple_params;
  for (int frame = 0; frame < 40; ++frame)
    In::Op::DisplaceRipple::advance(ripple_state, ripple_params);
  HS_EXPECT_EQ(ripple_state.phase, 40.0f);
  for (int frame = 0; frame < 40; ++frame)
    In::Op::DisplaceRipple::advance(ripple_state, ripple_params);
  HS_EXPECT_EQ(ripple_state.phase, 0.0f);

  // Warp batch: every op is PLANE->PLANE with "speed" first; the polar chart
  // carries the full sixteen-harmonic list; curl-flow is basis-only.
  for (const char *id :
       {"warp.affine.v2", "warp.wave-shear.v2", "warp.vector-noise.v2",
        "warp.mirror-tile.v2", "warp.polar-chart.v2", "warp.curl-flow.v2"}) {
    const In::OperatorDescriptor &warp = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(warp.input),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_EQ(static_cast<int>(warp.output),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_TRUE(std::string_view(warp.schema[0].id) == "speed");
  }
  const In::OperatorDescriptor &polar =
      *In::find_operator("warp.polar-chart.v2");
  constexpr size_t POLAR_FIELDS = In::Op::PolarChartParams::FIELDS.size();
  HS_EXPECT_EQ(polar.schema_count, POLAR_FIELDS + 2);
  HS_EXPECT_TRUE(std::string_view(polar.schema[POLAR_FIELDS].id) == "mode");
  const In::ParamFieldInfo &harmonic = polar.schema[POLAR_FIELDS + 1];
  HS_EXPECT_TRUE(std::string_view(harmonic.id) == "harmonic");
  HS_EXPECT_EQ(harmonic.enum_count, PB::Warp::MAX_POLAR_HARMONIC);
  HS_EXPECT_TRUE(std::string_view(harmonic.enum_ids[0]) == "h1");
  HS_EXPECT_TRUE(std::string_view(harmonic.enum_ids[15]) == "h16");
  const In::OperatorDescriptor &shear =
      *In::find_operator("warp.wave-shear.v2");
  const In::ParamFieldInfo &shear_envelope =
      shear.schema[shear.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(shear_envelope.id) == "envelope");
  HS_EXPECT_EQ(shear_envelope.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(shear_envelope.enum_ids[1]) ==
                 "projection-weight");
  const In::OperatorDescriptor &curl_flow =
      *In::find_operator("warp.curl-flow.v2");
  constexpr size_t CURL_FLOW_FIELDS = In::Op::CurlFlowParams::FIELDS.size();
  HS_EXPECT_EQ(curl_flow.schema_count, CURL_FLOW_FIELDS + 1);
  HS_EXPECT_TRUE(std::string_view(curl_flow.schema[CURL_FLOW_FIELDS].id) ==
                 "basis");

  // Projected sample batch: every crossing carries the union edge-width field
  // and the weight/coverage topology pair; only projected-noise adds a basis.
  for (const char *id :
       {"sample.grid.v2", "sample.twin-wave.v2", "sample.rings.v2",
        "sample.spiral.v2", "sample.lattice.v2", "sample.fractal.v2",
        "sample.tessellation.v2", "sample.projected-noise.v2"}) {
    const In::OperatorDescriptor &crossing = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(crossing.input),
                 static_cast<int>(In::CarrierId::PLANE));
    HS_EXPECT_EQ(static_cast<int>(crossing.output),
                 static_cast<int>(In::CarrierId::FIELD));
    bool has_edge_width = false;
    bool has_weight = false;
    bool has_coverage = false;
    bool has_basis = false;
    for (uint16_t field = 0; field < crossing.schema_count; ++field) {
      const std::string_view field_id{crossing.schema[field].id};
      has_edge_width |= field_id == "edge-width";
      has_weight |= field_id == "weight-mode";
      has_coverage |= field_id == "coverage-mode";
      has_basis |= field_id == "basis";
    }
    HS_EXPECT_TRUE(has_edge_width && has_weight && has_coverage);
    HS_EXPECT_EQ(has_basis,
                 std::string_view(id) == "sample.projected-noise.v2");
  }

  for (const char *id :
       {"sample.spherical-rings.v3", "sample.spherical-noise.v3"}) {
    const In::OperatorDescriptor &crossing = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(crossing.input),
                 static_cast<int>(In::CarrierId::SPHERE));
    HS_EXPECT_EQ(static_cast<int>(crossing.output),
                 static_cast<int>(In::CarrierId::FIELD));
    for (uint16_t field = 0; field < crossing.schema_count; ++field) {
      const std::string_view field_id{crossing.schema[field].id};
      HS_EXPECT_FALSE(field_id == "edge-width" || field_id == "weight-mode" ||
                      field_id == "coverage-mode");
    }
  }

  // Projection batch: the meridian-consuming projections extend the shared
  // family with central-meridian; gnomonic and the fast square Peirce do not.
  for (const char *id :
       {"project.stereographic.v2", "project.folded-sinusoidal.v2",
        "project.equirectangular.v2", "project.gnomonic.v2",
        "project.peirce.v2", "project.peirce-square-fast.v2",
        "project.bonne.v2", "project.airocean.v2"}) {
    const In::OperatorDescriptor &projection = *In::find_operator(id);
    HS_EXPECT_EQ(static_cast<int>(projection.input),
                 static_cast<int>(In::CarrierId::SPHERE));
    HS_EXPECT_EQ(static_cast<int>(projection.output),
                 static_cast<int>(In::CarrierId::PLANE));
    bool has_meridian = false;
    for (uint16_t field = 0; field < projection.schema_count; ++field)
      has_meridian |=
          std::string_view(projection.schema[field].id) == "central-meridian";
    const std::string_view id_view{id};
    HS_EXPECT_EQ(has_meridian, id_view != "project.stereographic.v2" &&
                                   id_view != "project.gnomonic.v2" &&
                                   id_view != "project.peirce-square-fast.v2");
  }
  const In::OperatorDescriptor &gnomonic =
      *In::find_operator("project.gnomonic.v2");
  const In::ParamFieldInfo &gnomonic_hemisphere =
      gnomonic.schema[gnomonic.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(gnomonic_hemisphere.id) == "hemisphere");
  HS_EXPECT_EQ(gnomonic_hemisphere.enum_count, 3);
  HS_EXPECT_TRUE(std::string_view(gnomonic_hemisphere.enum_ids[0]) == "folded");
  const In::OperatorDescriptor &bonne = *In::find_operator("project.bonne.v2");
  const In::ParamFieldInfo &bonne_hemisphere =
      bonne.schema[bonne.schema_count - 1];
  HS_EXPECT_TRUE(std::string_view(bonne_hemisphere.id) == "hemisphere");
  HS_EXPECT_EQ(bonne_hemisphere.enum_count, 2);
  HS_EXPECT_TRUE(std::string_view(bonne_hemisphere.enum_ids[0]) == "north");

  // Field batch: ridge is schema-free; iso-contour, smooth-bands and the value
  // cutout carry their own families.
  HS_EXPECT_EQ(In::find_operator("field.transfer.ridge.v2")->schema_count, 0);
  const In::OperatorDescriptor &iso =
      *In::find_operator("field.transfer.iso-contour.v2");
  HS_EXPECT_EQ(iso.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(iso.schema[0].id) == "iso-level");
  HS_EXPECT_TRUE(std::string_view(iso.schema[1].id) == "iso-width");
  const In::OperatorDescriptor &bands =
      *In::find_operator("field.transfer.smooth-bands.v2");
  HS_EXPECT_EQ(bands.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(bands.schema[0].id) == "band-count");
  HS_EXPECT_EQ(bands.schema[0].def, 4.0f);
  HS_EXPECT_TRUE(std::string_view(bands.schema[1].id) == "band-phase");
  const In::OperatorDescriptor &cutout =
      *In::find_operator("field.coverage.value-cutout.v2");
  HS_EXPECT_EQ(cutout.schema_count, 2);
  HS_EXPECT_TRUE(std::string_view(cutout.schema[0].id) == "cutout-threshold");
  HS_EXPECT_TRUE(std::string_view(cutout.schema[1].id) == "cutout-softness");
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
                          "\"per_op_overhead_bytes\":49,"
                          "\"per_param_name_bytes\":81}"));
  HS_EXPECT_TRUE(
      contains("\"carriers\":[\"sphere\",\"plane\",\"field\",\"color\"]"));
  HS_EXPECT_TRUE(contains("\"id\":\"sphere.rotate.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"project.stereographic.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.grid.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.spherical-rings.v3\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.fractal.v2\""));
  HS_EXPECT_TRUE(contains("\"id\":\"sample.tessellation.v2\""));
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

// --- sphere-op mirrors ----------------------------------------------------

/** Mirror frame for the sphere endomorphism batch; populated from the erased
    program's own blocks so both sides read identical raw state. */
struct SphereOpMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::CurlDisplaceParams curl;
  In::Op::DirectDisplaceParams direct;
  PB::Surface::PeriodicRippleParams ripple;
  MobiusParams mobius;
  float phase = 0.0f;
};

struct SphereOpMirrorBinding {
  using FrameState = SphereOpMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct CurlMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static float scale(const FrameState &frame) { return frame.curl.scale; }
  static float strength(const FrameState &frame) { return frame.curl.strength; }
  static PB::Surface::PreparedLoop prepare(const FrameState &frame) {
    return PB::Surface::prepare(frame.phase);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct DirectMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static float scale(const FrameState &frame) { return frame.direct.scale; }
  static float strength(const FrameState &frame) {
    return frame.direct.strength;
  }
  static PB::Surface::PreparedDirect prepare(const FrameState &frame) {
    return PB::Surface::prepare_direct(frame.phase, frame.direct.direction);
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct RippleMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const PB::Surface::PeriodicRippleParams &
  params(const FrameState &frame) {
    return frame.ripple;
  }
  static float phase(const FrameState &frame) {
    return frame.phase / frame.ripple.period;
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct MobiusMirrorProvider {
  using Binding = SphereOpMirrorBinding;
  using FrameState = SphereOpMirrorFrame;
  static const MobiusParams &params(const FrameState &frame) {
    return frame.mobius;
  }
};

/** Compiles a chain with one sphere endomorphism at entry 1 and applies the
    value set to every block. */
template <typename OpParams>
inline void arm_sphere_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"op", op_id},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 3), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 4), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** Erased-vs-bound parity of the sphere endomorphism at entry 1: identical
    seeds, bit-identical outputs across the sweep. */
template <typename BoundStage>
inline void expect_sphere_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx,
                                    const SphereOpMirrorFrame &mirror) {
  program.prepare(ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  for (const Vector &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.25f};
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&seed, out, ctx, program.param_block(1),
                   program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::SphereSample *>(out));
    const PB::SphereSample reference = BoundStage::run(seed, mirror, prepared);
    HS_EXPECT_TRUE(sphere_identical(erased, reference));
  }
}

inline SphereOpMirrorFrame sphere_op_mirror(In::ChainProgram &program) {
  SphereOpMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::DisplaceCurl::ID || id == In::Op::DisplaceDirect::ID) {
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 1);
    mirror.noise = &state.noise;
    mirror.phase = state.phase;
  }
  if (id == In::Op::DisplaceRipple::ID) {
    const auto &state = state_as<In::Op::RipplePhaseState>(program, 1);
    mirror.ripple = param_as<PB::Surface::PeriodicRippleParams>(program, 1);
    mirror.phase = state.phase;
  }
  if (id == In::Op::DisplaceCurl::ID)
    mirror.curl = param_as<In::Op::CurlDisplaceParams>(program, 1);
  if (id == In::Op::DisplaceDirect::ID)
    mirror.direct = param_as<In::Op::DirectDisplaceParams>(program, 1);
  if (id == In::Op::LensMobius::ID) {
    const auto &params = param_as<In::Op::MobiusChainParams>(program, 1);
    mirror.mobius =
        MobiusParams{params.a_re, params.a_im, params.b_re, params.b_im,
                     params.c_re, params.c_im, params.d_re, params.d_im};
  }
  return mirror;
}

template <::NoiseBasis Basis, typename Integrator>
inline void run_curl_displace_variant(In::ChainProgram &program,
                                      const In::FrameContext &ctx) {
  using Bound = typename PB::Stage::Displace<
      PB::Surface::CurlNoise<CurlMirrorProvider, Basis,
                             Integrator>>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

template <::NoiseBasis Basis>
inline void run_curl_displace_basis(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::CurlDisplaceParams>(program, 1);
  params.basis = static_cast<uint8_t>(Basis);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::EULER);
  run_curl_displace_variant<Basis, PB::Surface::Euler>(program, ctx);
  params.integrator = static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint>(program, ctx);
  params.integrator =
      static_cast<uint8_t>(PB::Surface::Integrator::MIDPOINT_2X);
  run_curl_displace_variant<Basis, PB::Surface::Midpoint2x>(program, ctx);
}

inline void test_shader_chain_parity_displace_curl() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<In::Op::CurlDisplaceParams>(
        program, In::Op::DisplaceCurl::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_curl_displace_basis<::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_displace_basis<::NoiseBasis::FBM3>(program, ctx);
    run_curl_displace_basis<::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

template <::NoiseBasis Basis>
inline void run_direct_displace_variant(In::ChainProgram &program,
                                        const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::DirectDisplaceParams>(program, 1);
  params.basis = static_cast<uint8_t>(Basis);
  using Bound = typename PB::Stage::Displace<PB::Surface::DirectNoise<
      DirectMirrorProvider, Basis>>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

inline void test_shader_chain_parity_displace_direct() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<In::Op::DirectDisplaceParams>(
        program, In::Op::DisplaceDirect::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_direct_displace_variant<::NoiseBasis::SIMPLEX>(program, ctx);
    run_direct_displace_variant<::NoiseBasis::FBM3>(program, ctx);
    run_direct_displace_variant<::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_displace_ripple() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sphere_op_chain<PB::Surface::PeriodicRippleParams>(
        program, In::Op::DisplaceRipple::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    using Bound = typename PB::Stage::Displace<PB::Surface::PeriodicRipple<
        RippleMirrorProvider>>::template Bind<SphereOpMirrorBinding>;
    expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
    program.clear();
  }
}

template <typename LensPolicy, typename OpParams>
inline void run_lens_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_sphere_op_chain<OpParams>(program, op_id, 3, set);
  const In::FrameContext ctx = shared_resources().context();
  using Bound = typename PB::Stage::Lens<LensPolicy>::template Bind<
      SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
  program.clear();
}

template <typename Policy>
inline void run_kaleidoscope_variant(In::ChainProgram &program,
                                     const In::FrameContext &ctx,
                                     In::Op::KaleidoscopeSymmetry symmetry) {
  param_as<In::Op::KaleidoscopeChainParams>(program, 1).symmetry =
      static_cast<uint8_t>(symmetry);
  using Bound =
      typename PB::Stage::Lens<Policy>::template Bind<SphereOpMirrorBinding>;
  expect_sphere_op_parity<Bound>(program, ctx, sphere_op_mirror(program));
}

inline void test_shader_chain_parity_lens_ops() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    run_lens_parity<PB::Lens::Glitch, PB::Lens::NoLensParams>(
        In::Op::LensGlitch::ID, set);
    run_lens_parity<PB::Lens::Twist, PB::Lens::NoLensParams>(
        In::Op::LensTwist::ID, set);
    run_lens_parity<PB::Lens::Mobius<MobiusMirrorProvider>,
                    In::Op::MobiusChainParams>(In::Op::LensMobius::ID, set);
  }

  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_sphere_op_chain<In::Op::KaleidoscopeChainParams>(
      program, In::Op::LensKaleidoscope::ID, 3, ValueSet::DEFAULTS);
  const In::FrameContext ctx = shared_resources().context();
  using Symmetry = In::Op::KaleidoscopeSymmetry;
  run_kaleidoscope_variant<PB::Lens::Kaleidoscope>(program, ctx,
                                                   Symmetry::AZIMUTHAL);
  run_kaleidoscope_variant<PB::Lens::TetrahedralKaleidoscope>(
      program, ctx, Symmetry::TETRAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::OctahedralKaleidoscope>(
      program, ctx, Symmetry::OCTAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::DodecahedralKaleidoscope>(
      program, ctx, Symmetry::DODECAHEDRAL);
  run_kaleidoscope_variant<PB::Lens::TriangularPrismKaleidoscope>(
      program, ctx, Symmetry::TRIANGULAR_PRISM);
  run_kaleidoscope_variant<PB::Lens::SquarePrismKaleidoscope>(
      program, ctx, Symmetry::SQUARE_PRISM);
  run_kaleidoscope_variant<PB::Lens::PentagonalPrismKaleidoscope>(
      program, ctx, Symmetry::PENTAGONAL_PRISM);
  run_kaleidoscope_variant<PB::Lens::HexagonalPrismKaleidoscope>(
      program, ctx, Symmetry::HEXAGONAL_PRISM);
  run_kaleidoscope_variant<PB::Lens::OctagonalPrismKaleidoscope>(
      program, ctx, Symmetry::OCTAGONAL_PRISM);
  program.clear();
}

// --- warp-op mirrors ------------------------------------------------------

/** Mirror frame for the planar warp batch. */
struct WarpMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::AffineWarpParams affine;
  In::Op::WaveShearWarpParams wave_shear;
  In::Op::VectorNoiseWarpParams vector_noise;
  In::Op::MirrorWarpParams mirror;
  In::Op::PolarChartParams polar;
  In::Op::CurlFlowParams curl;
  float phase = 0.0f;
  float rotation = 0.0f;
};

struct WarpMirrorBinding {
  using FrameState = WarpMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct AffineMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::AffineWarpParams &params(const FrameState &frame) {
    return frame.affine;
  }
  static PB::Warp::PreparedAffineSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.affine, frame.phase, frame.rotation, 1.0f);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct WaveShearMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::WaveShearWarpParams &params(const FrameState &frame) {
    return frame.wave_shear;
  }
  static PB::Warp::PreparedRotation prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.wave_shear, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct VectorNoiseMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::VectorNoiseWarpParams &params(const FrameState &frame) {
    return frame.vector_noise;
  }
  static PB::Warp::PreparedVectorNoiseSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.vector_noise, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static bool path_length_required(const FrameState &) { return true; }
};

struct MirrorTileMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::MirrorWarpParams &params(const FrameState &frame) {
    return frame.mirror;
  }
  static PB::Warp::PreparedMirrorSlot prepare(const FrameState &frame) {
    return PB::Warp::prepare(frame.mirror, frame.phase);
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct PolarChartMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::PolarChartParams &params(const FrameState &frame) {
    return frame.polar;
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static bool path_length_required(const FrameState &) { return true; }
};

struct CurlFlowMirrorProvider {
  using Binding = WarpMirrorBinding;
  using FrameState = WarpMirrorFrame;
  static const In::Op::CurlFlowParams &params(const FrameState &frame) {
    return frame.curl;
  }
  static float phase(const FrameState &frame) { return frame.phase; }
  static const FastNoiseLite &noise(const FrameState &frame) {
    return *frame.noise;
  }
  static bool path_length_required(const FrameState &) { return true; }
};

/** Compiles a chain with one planar warp at entry 2 and applies the value set
    to every block. */
template <typename OpParams>
inline void arm_warp_op_chain(In::ChainProgram &program, const char *op_id,
                              int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"warp", op_id},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<OpParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 3), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 4), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** The warp's input carrier for @p view: the erased camera and projection
    entries run over the seed, so both parity sides read one real sample. */
inline PB::PlaneSample warp_input(In::ChainProgram &program,
                                  const In::FrameContext &ctx,
                                  const Vector &view) {
  const PB::SphereSample seed{view, 0.0f};
  alignas(In::SLOT_ALIGN) uint8_t rotated[In::SLOT_SIZE];
  program.ops()[0].op->runtime.run(&seed, rotated, ctx, program.param_block(0),
                                   program.prepared_block(0));
  alignas(In::SLOT_ALIGN) uint8_t projected[In::SLOT_SIZE];
  program.ops()[1].op->runtime.run(rotated, projected, ctx,
                                   program.param_block(1),
                                   program.prepared_block(1));
  return *std::launder(reinterpret_cast<PB::PlaneSample *>(projected));
}

/** Erased-vs-bound parity of the warp at entry 2. */
template <typename BoundStage>
inline void expect_warp_op_parity(In::ChainProgram &program,
                                  const In::FrameContext &ctx,
                                  const WarpMirrorFrame &mirror) {
  program.prepare(ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  for (const Vector &view : sweep_views()) {
    const PB::PlaneSample input = warp_input(program, ctx, view);
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&input, out, ctx, program.param_block(2),
                   program.prepared_block(2));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::PlaneSample *>(out));
    const PB::PlaneSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(plane_identical(erased, reference));
  }
}

inline WarpMirrorFrame warp_mirror(In::ChainProgram &program) {
  WarpMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::WarpAffine::ID) {
    const auto &state = state_as<In::Op::AffineClockState>(program, 2);
    mirror.phase = state.phase;
    mirror.rotation = state.rotation;
    mirror.affine = param_as<In::Op::AffineWarpParams>(program, 2);
  } else if (id == In::Op::WarpVectorNoise::ID ||
             id == In::Op::WarpCurlFlow::ID) {
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 2);
    mirror.noise = &state.noise;
    mirror.phase = state.phase;
    if (id == In::Op::WarpVectorNoise::ID)
      mirror.vector_noise = param_as<In::Op::VectorNoiseWarpParams>(program, 2);
    else
      mirror.curl = param_as<In::Op::CurlFlowParams>(program, 2);
  } else {
    mirror.phase = state_as<In::Op::WarpPhaseState>(program, 2).phase;
    if (id == In::Op::WarpWaveShear::ID)
      mirror.wave_shear = param_as<In::Op::WaveShearWarpParams>(program, 2);
    else if (id == In::Op::WarpMirrorTile::ID)
      mirror.mirror = param_as<In::Op::MirrorWarpParams>(program, 2);
    else if (id == In::Op::WarpPolarChart::ID)
      mirror.polar = param_as<In::Op::PolarChartParams>(program, 2);
  }
  return mirror;
}

inline void test_shader_chain_parity_warp_affine_mirror() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    const In::FrameContext ctx = shared_resources().context();
    arm_warp_op_chain<In::Op::AffineWarpParams>(program, In::Op::WarpAffine::ID,
                                                4, set);
    using BoundAffine =
        typename PB::Stage::Warp<PB::Warp::AffineFrame<AffineMirrorProvider>>::
            template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundAffine>(program, ctx, warp_mirror(program));
    program.clear();

    arm_warp_op_chain<In::Op::MirrorWarpParams>(
        program, In::Op::WarpMirrorTile::ID, 4, set);
    using BoundMirror = typename PB::Stage::Warp<PB::Warp::MirrorTile<
        MirrorTileMirrorProvider>>::template Bind<WarpMirrorBinding>;
    expect_warp_op_parity<BoundMirror>(program, ctx, warp_mirror(program));
    program.clear();
  }
}

template <typename Envelope>
inline void run_wave_shear_variant(In::ChainProgram &program,
                                   const In::FrameContext &ctx,
                                   In::Op::WarpEnvelope envelope) {
  param_as<In::Op::WaveShearWarpParams>(program, 2).envelope =
      static_cast<uint8_t>(envelope);
  using Bound = typename PB::Stage::Warp<PB::Warp::WaveShear<
      WaveShearMirrorProvider, Envelope>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

inline void test_shader_chain_parity_warp_wave_shear() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::WaveShearWarpParams>(
        program, In::Op::WarpWaveShear::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_wave_shear_variant<PB::Warp::FlatEnvelope>(program, ctx,
                                                   In::Op::WarpEnvelope::FLAT);
    run_wave_shear_variant<PB::Warp::ProjectionWeightEnvelope>(
        program, ctx, In::Op::WarpEnvelope::PROJECTION_WEIGHT);
    run_wave_shear_variant<PB::Warp::EdgeFadeEnvelope>(
        program, ctx, In::Op::WarpEnvelope::EDGE_FADE);
    program.clear();
  }
}

template <::NoiseBasis Basis, typename Envelope>
inline void run_vector_noise_variant(In::ChainProgram &program,
                                     const In::FrameContext &ctx,
                                     In::Op::WarpEnvelope envelope) {
  auto &params = param_as<In::Op::VectorNoiseWarpParams>(program, 2);
  params.basis = static_cast<uint8_t>(Basis);
  params.envelope = static_cast<uint8_t>(envelope);
  using Bound = typename PB::Stage::Warp<
      PB::Warp::VectorNoise<VectorNoiseMirrorProvider, Basis,
                            Envelope>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

template <::NoiseBasis Basis>
inline void run_vector_noise_basis(In::ChainProgram &program,
                                   const In::FrameContext &ctx) {
  run_vector_noise_variant<Basis, PB::Warp::FlatEnvelope>(
      program, ctx, In::Op::WarpEnvelope::FLAT);
  run_vector_noise_variant<Basis, PB::Warp::ProjectionWeightEnvelope>(
      program, ctx, In::Op::WarpEnvelope::PROJECTION_WEIGHT);
  run_vector_noise_variant<Basis, PB::Warp::EdgeFadeEnvelope>(
      program, ctx, In::Op::WarpEnvelope::EDGE_FADE);
}

inline void test_shader_chain_parity_warp_vector_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::VectorNoiseWarpParams>(
        program, In::Op::WarpVectorNoise::ID, 4, set);
    FastNoiseLite authored_noise;
    authored_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    authored_noise.SetSeed(static_cast<int32_t>(
        In::instance_hash("warp", In::Op::WarpVectorNoise::ID)));
    authored_noise.SetFrequency(1.0f);
    const auto &state = state_as<In::Op::NoisePhaseState>(program, 2);
    HS_EXPECT_EQ(state.noise.GetNoise(0.25f, -0.5f, 0.75f),
                 authored_noise.GetNoise(0.25f, -0.5f, 0.75f));
    const In::FrameContext ctx = shared_resources().context();
    run_vector_noise_basis<::NoiseBasis::SIMPLEX>(program, ctx);
    run_vector_noise_basis<::NoiseBasis::FBM3>(program, ctx);
    run_vector_noise_basis<::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

/** Two instances of one noise operator own decorrelated fields: each is seeded
    from its own instance identity, not from a shared constant. */
inline void test_shader_chain_noise_instances_decorrelate() {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"warp-a", In::Op::WarpVectorNoise::ID},
      {"warp-b", In::Op::WarpVectorNoise::ID},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  const auto &first = state_as<In::Op::NoisePhaseState>(program, 2);
  const auto &second = state_as<In::Op::NoisePhaseState>(program, 3);
  bool differs = false;
  for (const Vector &view : sweep_views())
    if (first.noise.GetNoise(view.x, view.y, view.z) !=
        second.noise.GetNoise(view.x, view.y, view.z))
      differs = true;
  HS_EXPECT_TRUE(differs);
  // Deterministic: the seed is the instance's stable hash, so a recompile of
  // the same shape reconstructs the same field.
  FastNoiseLite authored;
  authored.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  authored.SetSeed(static_cast<int32_t>(
      In::instance_hash("warp-b", In::Op::WarpVectorNoise::ID)));
  authored.SetFrequency(1.0f);
  HS_EXPECT_EQ(second.noise.GetNoise(0.25f, -0.5f, 0.75f),
               authored.GetNoise(0.25f, -0.5f, 0.75f));
  program.clear();
}

template <typename Mode, uint8_t Harmonic>
inline void run_polar_variant(In::ChainProgram &program,
                              const In::FrameContext &ctx) {
  auto &params = param_as<In::Op::PolarChartParams>(program, 2);
  params.mode =
      static_cast<uint8_t>(std::is_same_v<Mode, PB::Warp::LogarithmicPolar>
                               ? In::Op::PolarMode::LOGARITHMIC
                               : In::Op::PolarMode::LINEAR);
  params.harmonic = Harmonic - 1;
  using Bound = typename PB::Stage::Warp<
      PB::Warp::PolarChart<PolarChartMirrorProvider, Mode,
                           Harmonic>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

template <typename Mode>
inline void run_polar_harmonics(In::ChainProgram &program,
                                const In::FrameContext &ctx) {
  [&]<size_t... I>(std::index_sequence<I...>) {
    (run_polar_variant<Mode, static_cast<uint8_t>(I + 1)>(program, ctx), ...);
  }(std::make_index_sequence<PB::Warp::MAX_POLAR_HARMONIC>{});
}

inline void test_shader_chain_parity_warp_polar_chart() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::PolarChartParams>(
        program, In::Op::WarpPolarChart::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_polar_harmonics<PB::Warp::LinearPolar>(program, ctx);
    run_polar_harmonics<PB::Warp::LogarithmicPolar>(program, ctx);
    program.clear();
  }
}

template <::NoiseBasis Basis>
inline void run_curl_flow_variant(In::ChainProgram &program,
                                  const In::FrameContext &ctx) {
  param_as<In::Op::CurlFlowParams>(program, 2).basis =
      static_cast<uint8_t>(Basis);
  using Bound = typename PB::Stage::Warp<
      PB::Warp::CurlFlow<CurlFlowMirrorProvider, Basis,
                         PB::Warp::Euler1>>::template Bind<WarpMirrorBinding>;
  expect_warp_op_parity<Bound>(program, ctx, warp_mirror(program));
}

inline void test_shader_chain_parity_warp_curl_flow() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_warp_op_chain<In::Op::CurlFlowParams>(program, In::Op::WarpCurlFlow::ID,
                                              4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_curl_flow_variant<::NoiseBasis::SIMPLEX>(program, ctx);
    run_curl_flow_variant<::NoiseBasis::FBM3>(program, ctx);
    run_curl_flow_variant<::NoiseBasis::RIDGED3>(program, ctx);
    program.clear();
  }
}

// --- projection-op mirrors ------------------------------------------------

/** Mirror frame for the projection batch. */
struct ProjMirrorFrame {
  Quaternion conjugate;
  In::Op::MeridianProjectChainParams meridian;
  In::Op::GnomonicChainParams gnomonic;
  In::Op::BonneChainParams bonne;
};

struct ProjMirrorBinding {
  using FrameState = ProjMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

template <typename Derived> struct ProjMirrorBase {
  using Binding = ProjMirrorBinding;
  using FrameState = ProjMirrorFrame;
  static const Quaternion &conjugate(const ProjMirrorFrame &frame) {
    return frame.conjugate;
  }
};

struct MeridianProjMirror : ProjMirrorBase<MeridianProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.meridian.central_meridian;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.meridian.singularity_fade;
  }
};

struct GnomonicProjMirror : ProjMirrorBase<GnomonicProjMirror> {
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.gnomonic.singularity_fade;
  }
};

struct PeirceProjMirror : ProjMirrorBase<PeirceProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.meridian.central_meridian;
  }
  static float layout_scroll(const ProjMirrorFrame &) { return 0.0f; }
  static float coordinate_scale(const ProjMirrorFrame &) {
    return In::Op::PROJECT_COORDINATE_SCALE;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.meridian.singularity_fade;
  }
};

struct PeirceFastProjMirror : ProjMirrorBase<PeirceFastProjMirror> {
  static constexpr bool ZERO_CENTRAL_MERIDIAN = true;
  static float coordinate_scale(const ProjMirrorFrame &) {
    return In::Op::PROJECT_COORDINATE_SCALE;
  }
  static float singularity_fade(const ProjMirrorFrame &frame) {
    return frame.meridian.singularity_fade;
  }
};

struct BonneProjMirror : ProjMirrorBase<BonneProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.bonne.central_meridian;
  }
  static float standard_parallel(const ProjMirrorFrame &) {
    return In::Op::BONNE_STANDARD_PARALLEL;
  }
  static float coordinate_scale(const ProjMirrorFrame &) {
    return In::Op::PROJECT_COORDINATE_SCALE;
  }
};

struct AiroceanProjMirror : ProjMirrorBase<AiroceanProjMirror> {
  static float central_meridian(const ProjMirrorFrame &frame) {
    return frame.meridian.central_meridian;
  }
  static float coordinate_scale(const ProjMirrorFrame &) {
    return In::Op::PROJECT_COORDINATE_SCALE;
  }
};

/** Compiles a chain with one projection at entry 1 and applies the value set
    to every block. */
template <typename OpParams>
inline void arm_project_op_chain(In::ChainProgram &program, const char *op_id,
                                 int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", op_id},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

inline ProjMirrorFrame project_mirror(In::ChainProgram &program,
                                      const In::FrameContext &ctx) {
  ProjMirrorFrame mirror;
  const auto &state = state_as<In::Op::SpatialWalkState>(program, 1);
  mirror.conjugate = (make_rotation(Y_AXIS, state.spin_phase) *
                      ctx.projection_base * state.wander)
                         .conjugate();
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::ProjectGnomonic::ID) {
    mirror.gnomonic = param_as<In::Op::GnomonicChainParams>(program, 1);
  } else if (id == In::Op::ProjectBonne::ID) {
    mirror.bonne = param_as<In::Op::BonneChainParams>(program, 1);
  } else if (id == In::Op::ProjectPeirceSquareFast::ID) {
    static_cast<In::Op::ProjectChainParams &>(mirror.meridian) =
        param_as<In::Op::ProjectChainParams>(program, 1);
  } else {
    mirror.meridian = param_as<In::Op::MeridianProjectChainParams>(program, 1);
  }
  return mirror;
}

/** Erased-vs-bound parity of the projection at entry 1. */
template <typename BoundStage>
inline void expect_project_op_parity(In::ChainProgram &program,
                                     const In::FrameContext &ctx) {
  program.prepare(ctx);
  const ProjMirrorFrame mirror = project_mirror(program, ctx);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[1].op;
  for (const Vector &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.25f};
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&seed, out, ctx, program.param_block(1),
                   program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::PlaneSample *>(out));
    const PB::PlaneSample reference = BoundStage::run(seed, mirror, prepared);
    HS_EXPECT_TRUE(plane_identical(erased, reference));
  }
}

template <typename Policy, typename OpParams>
inline void run_project_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_project_op_chain<OpParams>(program, op_id, 4, set);
  const In::FrameContext ctx = shared_resources().context();
  using Bound =
      typename PB::Stage::Project<Policy>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
  program.clear();
}

inline void test_shader_chain_parity_project_ops() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    run_project_parity<PB::Projection::FoldedSinusoidal<MeridianProjMirror>,
                       In::Op::MeridianProjectChainParams>(
        In::Op::ProjectFoldedSinusoidal::ID, set);
    run_project_parity<PB::Projection::Equirectangular<MeridianProjMirror>,
                       In::Op::MeridianProjectChainParams>(
        In::Op::ProjectEquirectangular::ID, set);
    run_project_parity<
        PB::Projection::Peirce<PeirceProjMirror, In::Op::PEIRCE_SQUARE_LAYOUT,
                               true>,
        In::Op::MeridianProjectChainParams>(In::Op::ProjectPeirce::ID, set);
    run_project_parity<PB::Projection::PeirceFastSquare<PeirceFastProjMirror>,
                       In::Op::ProjectChainParams>(
        In::Op::ProjectPeirceSquareFast::ID, set);
    run_project_parity<
        PB::Projection::Airocean<AiroceanProjMirror, false, true>,
        In::Op::MeridianProjectChainParams>(In::Op::ProjectAirocean::ID, set);
  }
}

template <PB::Projection::GnomonicHemisphere Hemisphere>
inline void run_gnomonic_variant(In::ChainProgram &program,
                                 const In::FrameContext &ctx) {
  param_as<In::Op::GnomonicChainParams>(program, 1).hemisphere =
      static_cast<uint8_t>(Hemisphere);
  using Bound = typename PB::Stage::Project<PB::Projection::Gnomonic<
      GnomonicProjMirror, Hemisphere>>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
}

template <bool North>
inline void run_bonne_variant(In::ChainProgram &program,
                              const In::FrameContext &ctx) {
  param_as<In::Op::BonneChainParams>(program, 1).hemisphere = North ? 0 : 1;
  using Bound = typename PB::Stage::Project<PB::Projection::Bonne<
      BonneProjMirror, North>>::template Bind<ProjMirrorBinding>;
  expect_project_op_parity<Bound>(program, ctx);
}

inline void test_shader_chain_parity_project_hemispheres() {
  using Hemisphere = PB::Projection::GnomonicHemisphere;
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_project_op_chain<In::Op::GnomonicChainParams>(
        program, In::Op::ProjectGnomonic::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_gnomonic_variant<Hemisphere::FOLDED>(program, ctx);
    run_gnomonic_variant<Hemisphere::FRONT>(program, ctx);
    run_gnomonic_variant<Hemisphere::BACK>(program, ctx);
    program.clear();

    arm_project_op_chain<In::Op::BonneChainParams>(
        program, In::Op::ProjectBonne::ID, 4, set);
    run_bonne_variant<true>(program, ctx);
    run_bonne_variant<false>(program, ctx);
    program.clear();
  }
}

// --- field-op mirrors -----------------------------------------------------

/** Mirror frame for the FIELD endomorphism batch. */
struct FieldMirrorFrame {
  In::Op::IsoContourChainParams iso;
  In::Op::SmoothBandsChainParams bands;
  In::Op::ValueCutoutChainParams cutout;
};

struct FieldMirrorBinding {
  using FrameState = FieldMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

struct IsoContourFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float iso_level(const FieldMirrorFrame &frame) {
    return frame.iso.iso_level;
  }
  static float iso_width(const FieldMirrorFrame &frame) {
    return frame.iso.iso_width;
  }
};

struct SmoothBandsFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float band_count(const FieldMirrorFrame &frame) {
    return frame.bands.band_count;
  }
  static float band_phase(const FieldMirrorFrame &frame) {
    return frame.bands.band_phase;
  }
};

struct ValueCutoutFieldMirror {
  using Binding = FieldMirrorBinding;
  using FrameState = FieldMirrorFrame;
  static float cutout_threshold(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_threshold;
  }
  static float cutout_softness(const FieldMirrorFrame &frame) {
    return frame.cutout.cutout_softness;
  }
};

/** Compiles a chain with one FIELD endomorphism at entry 3 and applies the
    value set to every block. */
template <typename OpParams>
inline void arm_field_op_chain(In::ChainProgram &program, const char *op_id,
                               int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"op", op_id},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GridSampleParams>(program, 2), set);
  apply_value_set(param_as<OpParams>(program, 3), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 4), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

/** The FIELD op's input carrier for @p view: the erased camera, projection
    and sample entries run over the seed. */
inline PB::FieldSample field_input(In::ChainProgram &program,
                                   const In::FrameContext &ctx,
                                   const Vector &view) {
  const PB::PlaneSample plane = warp_input(program, ctx, view);
  alignas(In::SLOT_ALIGN) uint8_t sampled[In::SLOT_SIZE];
  program.ops()[2].op->runtime.run(&plane, sampled, ctx, program.param_block(2),
                                   program.prepared_block(2));
  return *std::launder(reinterpret_cast<PB::FieldSample *>(sampled));
}

inline FieldMirrorFrame field_mirror(In::ChainProgram &program) {
  FieldMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[3].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::TransferIsoContour::ID)
    mirror.iso = param_as<In::Op::IsoContourChainParams>(program, 3);
  else if (id == In::Op::TransferSmoothBands::ID)
    mirror.bands = param_as<In::Op::SmoothBandsChainParams>(program, 3);
  else if (id == In::Op::CoverageValueCutout::ID)
    mirror.cutout = param_as<In::Op::ValueCutoutChainParams>(program, 3);
  return mirror;
}

/** Erased-vs-bound parity of the FIELD endomorphism at entry 3. */
template <typename BoundStage>
inline void expect_field_op_parity(In::ChainProgram &program,
                                   const In::FrameContext &ctx) {
  program.prepare(ctx);
  const FieldMirrorFrame mirror = field_mirror(program);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[3].op;
  for (const Vector &view : sweep_views()) {
    const PB::FieldSample input = field_input(program, ctx, view);
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&input, out, ctx, program.param_block(3),
                   program.prepared_block(3));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::FieldSample *>(out));
    const PB::FieldSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(field_identical(erased, reference));
  }
}

template <typename BoundStage, typename OpParams>
inline void run_field_parity(const char *op_id, ValueSet set) {
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_field_op_chain<OpParams>(program, op_id, 4, set);
  const In::FrameContext ctx = shared_resources().context();
  expect_field_op_parity<BoundStage>(program, ctx);
  program.clear();
}

inline void test_shader_chain_parity_field_ops() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    run_field_parity<typename PB::Stage::Transfer<PB::Transfer::Ridge>::
                         template Bind<FieldMirrorBinding>,
                     PB::Transfer::NoValueParams>(In::Op::TransferRidge::ID,
                                                  set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::IsoContour<
            IsoContourFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::IsoContourChainParams>(In::Op::TransferIsoContour::ID, set);
    run_field_parity<
        typename PB::Stage::Transfer<PB::Transfer::SmoothBands<
            SmoothBandsFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::SmoothBandsChainParams>(In::Op::TransferSmoothBands::ID, set);
    run_field_parity<
        typename PB::Stage::Coverage<PB::Coverage::ValueCutout<
            ValueCutoutFieldMirror>>::template Bind<FieldMirrorBinding>,
        In::Op::ValueCutoutChainParams>(In::Op::CoverageValueCutout::ID, set);
  }
}

// --- sample-op mirrors ----------------------------------------------------

/** Mirror frame for the sample source batch. */
struct SampleMirrorFrame {
  const FastNoiseLite *noise = nullptr;
  In::Op::TwinWaveSampleParams twin_wave;
  In::Op::RingsSampleParams rings;
  In::Op::SphericalRingsSampleParams spherical_rings;
  In::Op::SpiralSampleParams spiral;
  In::Op::LatticeSampleParams lattice;
  In::Op::FractalSampleParams fractal;
  In::Op::TessellationSampleParams tessellation;
  In::Op::ProjectedNoiseSampleParams projected;
  In::Op::SphericalNoiseSampleParams spherical;
  Vector ring_axis = Y_AXIS;
  float ring_phase = 0.0f;
  float primary = 0.0f;
  float secondary = 0.0f;
  float angle = 0.0f;
  float noise_time = 0.0f;
  float edge_width = 0.1f;
};

struct SampleMirrorBinding {
  using FrameState = SampleMirrorFrame;
  using Instrumentation = PB::NoInstrumentation;
};

template <typename Derived> struct ClockedSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static PB::Source::PreparedSource prepare(const FrameState &frame) {
    return PB::Source::prepare(frame.primary, frame.secondary, frame.angle);
  }
};

struct TwinWaveSampleMirror : ClockedSampleMirror<TwinWaveSampleMirror> {
  static const In::Op::TwinWaveSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.twin_wave;
  }
};

struct RingsSampleMirror : ClockedSampleMirror<RingsSampleMirror> {
  static const In::Op::RingsSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.rings;
  }
};

struct SphericalRingsSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::SphericalRingsSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.spherical_rings;
  }
  static PB::Source::PreparedSphericalRings
  prepare(const SampleMirrorFrame &frame) {
    return {frame.ring_axis, frame.ring_phase};
  }
};

struct SpiralSampleMirror : ClockedSampleMirror<SpiralSampleMirror> {
  static const In::Op::SpiralSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.spiral;
  }
};

struct LatticeSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const In::Op::LatticeSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.lattice;
  }
};

struct FractalSampleMirror : ClockedSampleMirror<FractalSampleMirror> {
  static const In::Op::FractalSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.fractal;
  }
};

struct TessellationSampleMirror
    : ClockedSampleMirror<TessellationSampleMirror> {
  static const In::Op::TessellationSampleParams &
  params(const SampleMirrorFrame &frame) {
    return frame.tessellation;
  }
};

struct ProjectedNoiseSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const FastNoiseLite &noise(const SampleMirrorFrame &frame) {
    return *frame.noise;
  }
  static float noise_scale(const SampleMirrorFrame &frame) {
    return frame.projected.noise_scale;
  }
  static float noise_time(const SampleMirrorFrame &frame) {
    return frame.noise_time;
  }
  static float noise_contrast(const SampleMirrorFrame &frame) {
    return frame.projected.noise_contrast;
  }
};

struct SphericalNoiseSampleMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static const FastNoiseLite &noise(const SampleMirrorFrame &frame) {
    return *frame.noise;
  }
  static float noise_scale(const SampleMirrorFrame &frame) {
    return frame.spherical.noise_scale;
  }
  static float noise_time(const SampleMirrorFrame &frame) {
    return frame.noise_time;
  }
  static float noise_contrast(const SampleMirrorFrame &frame) {
    return frame.spherical.noise_contrast;
  }
};

/** Edge-width provider for the mirror's edge-fade coverage policy. */
struct SampleEdgeMirror {
  using Binding = SampleMirrorBinding;
  using FrameState = SampleMirrorFrame;
  static float edge_width(const SampleMirrorFrame &frame) {
    return frame.edge_width;
  }
};

/** Compiles a chain with one sample crossing at entry 2 and applies the value
    set to every block. */
template <typename OpParams>
inline void arm_sample_op_chain(In::ChainProgram &program, const char *op_id,
                                int frames, ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<In::Op::ProjectChainParams>(program, 1), set);
  apply_value_set(param_as<OpParams>(program, 2), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 3), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

template <typename OpParams>
inline void arm_spherical_sample_op_chain(In::ChainProgram &program,
                                          const char *op_id, int frames,
                                          ValueSet set) {
  const In::ChainEntryRequest chain[] = {
      {"camera", "sphere.rotate.v2"},
      {"source", op_id},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal refusal =
      program.compile(std::span<const In::ChainEntryRequest>(chain));
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::OK));
  apply_value_set(param_as<In::Op::RotateChainParams>(program, 0), set);
  apply_value_set(param_as<OpParams>(program, 1), set);
  apply_value_set(param_as<In::Op::GeneratedPaletteParams>(program, 2), set);
  for (int frame = 0; frame < frames; ++frame)
    program.advance();
}

inline SampleMirrorFrame sample_mirror(In::ChainProgram &program,
                                       size_t source_index = 2) {
  SampleMirrorFrame mirror;
  const In::OperatorDescriptor &op = *program.ops()[source_index].op;
  const std::string_view id{op.operator_id};
  if (id == In::Op::SampleProjectedNoise::ID ||
      id == In::Op::SampleSphericalNoise::ID) {
    const auto &state =
        state_as<In::Op::NoisePhaseState>(program, source_index);
    mirror.noise = &state.noise;
    mirror.noise_time = state.phase;
    if (id == In::Op::SampleProjectedNoise::ID) {
      mirror.projected =
          param_as<In::Op::ProjectedNoiseSampleParams>(program, source_index);
      mirror.edge_width = mirror.projected.edge_width;
    } else {
      mirror.spherical =
          param_as<In::Op::SphericalNoiseSampleParams>(program, source_index);
    }
  } else if (id == In::Op::SampleLattice::ID) {
    mirror.lattice =
        param_as<In::Op::LatticeSampleParams>(program, source_index);
    mirror.edge_width = mirror.lattice.edge_width;
  } else if (id == In::Op::SampleSphericalRings::ID) {
    const auto &state =
        state_as<In::Op::SphericalRingsState>(program, source_index);
    mirror.spherical_rings =
        param_as<In::Op::SphericalRingsSampleParams>(program, source_index);
    const Quaternion orientation =
        make_rotation(X_AXIS, state.walk.spin_phase) * state.walk.wander;
    mirror.ring_axis = rotate(Y_AXIS, orientation);
    mirror.ring_phase = state.phase;
  } else {
    const auto &state =
        state_as<In::Op::SourceClockState>(program, source_index);
    mirror.primary = state.primary;
    mirror.secondary = state.secondary;
    mirror.angle = state.angle;
    if (id == In::Op::SampleTwinWave::ID) {
      mirror.twin_wave =
          param_as<In::Op::TwinWaveSampleParams>(program, source_index);
      mirror.edge_width = mirror.twin_wave.edge_width;
    } else if (id == In::Op::SampleRings::ID) {
      mirror.rings = param_as<In::Op::RingsSampleParams>(program, source_index);
      mirror.edge_width = mirror.rings.edge_width;
    } else if (id == In::Op::SampleSpiral::ID) {
      mirror.spiral =
          param_as<In::Op::SpiralSampleParams>(program, source_index);
      mirror.edge_width = mirror.spiral.edge_width;
    } else if (id == In::Op::SampleFractal::ID) {
      mirror.fractal =
          param_as<In::Op::FractalSampleParams>(program, source_index);
      mirror.edge_width = mirror.fractal.edge_width;
    } else if (id == In::Op::SampleTessellation::ID) {
      mirror.tessellation =
          param_as<In::Op::TessellationSampleParams>(program, source_index);
      mirror.edge_width = mirror.tessellation.edge_width;
    }
  }
  return mirror;
}

template <typename SourceP>
inline void expect_spherical_sample_op_parity(In::ChainProgram &program,
                                              const In::FrameContext &ctx) {
  using BoundStage = typename PB::Stage::SampleSphere<SourceP>::template Bind<
      SampleMirrorBinding>;
  program.prepare(ctx);
  const SampleMirrorFrame mirror = sample_mirror(program, 1);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &rotate_op = *program.ops()[0].op;
  const In::OperatorDescriptor &sample_op = *program.ops()[1].op;
  for (const Vector &view : sweep_views()) {
    const PB::SphereSample seed{view, 0.0f};
    alignas(In::SLOT_ALIGN) uint8_t rotated_out[In::SLOT_SIZE];
    rotate_op.runtime.run(&seed, rotated_out, ctx, program.param_block(0),
                          program.prepared_block(0));
    const auto &input =
        *std::launder(reinterpret_cast<PB::SphereSample *>(rotated_out));
    alignas(In::SLOT_ALIGN) uint8_t sampled_out[In::SLOT_SIZE];
    sample_op.runtime.run(&input, sampled_out, ctx, program.param_block(1),
                          program.prepared_block(1));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::FieldSample *>(sampled_out));
    const PB::FieldSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(field_identical(erased, reference));
  }
}

/** Erased-vs-bound parity of the sample crossing at entry 2. */
template <typename BoundStage>
inline void expect_sample_op_parity(In::ChainProgram &program,
                                    const In::FrameContext &ctx) {
  program.prepare(ctx);
  const SampleMirrorFrame mirror = sample_mirror(program);
  const typename BoundStage::Prepared prepared = BoundStage::prepare(mirror);
  const In::OperatorDescriptor &op = *program.ops()[2].op;
  for (const Vector &view : sweep_views()) {
    const PB::PlaneSample input = warp_input(program, ctx, view);
    alignas(In::SLOT_ALIGN) uint8_t out[In::SLOT_SIZE];
    op.runtime.run(&input, out, ctx, program.param_block(2),
                   program.prepared_block(2));
    const auto &erased =
        *std::launder(reinterpret_cast<PB::FieldSample *>(out));
    const PB::FieldSample reference = BoundStage::run(input, mirror, prepared);
    HS_EXPECT_TRUE(field_identical(erased, reference));
  }
}

/** Address of a topology enum8 inside entry @p index's param block, resolved
    through the erased param-address channel. */
inline uint8_t *topology_byte(In::ChainProgram &program, size_t index,
                              const char *field_id) {
  const In::OperatorDescriptor &op = *program.ops()[index].op;
  for (uint16_t field = 0; field < op.schema_count; ++field)
    if (std::string_view(op.schema[field].id) == field_id)
      return static_cast<uint8_t *>(
          op.runtime.param_address(program.param_block(index), field));
  return nullptr;
}

template <typename SourceP, typename WeightP>
inline void run_sample_coverage_cases(In::ChainProgram &program,
                                      const In::FrameContext &ctx,
                                      uint8_t coverage) {
  using EdgeP = PB::ProjectedCoverage::EdgeFade<SampleEdgeMirror>;
  switch (static_cast<In::Op::CoverageMode>(coverage)) {
  case In::Op::CoverageMode::NONE:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP,
        PB::ProjectedCoverage::None>::template Bind<SampleMirrorBinding>>(
        program, ctx);
    break;
  case In::Op::CoverageMode::WEIGHT:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP,
        PB::ProjectedCoverage::Weight>::template Bind<SampleMirrorBinding>>(
        program, ctx);
    break;
  case In::Op::CoverageMode::WEIGHT_SQUARED:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP, PB::ProjectedCoverage::WeightSquared>::
                                template Bind<SampleMirrorBinding>>(program,
                                                                    ctx);
    break;
  case In::Op::CoverageMode::EDGE_FADE:
    expect_sample_op_parity<typename PB::Stage::Sample<
        SourceP, WeightP, EdgeP>::template Bind<SampleMirrorBinding>>(program,
                                                                      ctx);
    break;
  }
}

/** Runs the full weight-by-coverage topology matrix of the sample crossing at
    entry 2 against @p SourceP. */
template <typename SourceP>
inline void run_sample_op_matrix(In::ChainProgram &program,
                                 const In::FrameContext &ctx) {
  uint8_t *weight = topology_byte(program, 2, "weight-mode");
  uint8_t *coverage = topology_byte(program, 2, "coverage-mode");
  HS_EXPECT_TRUE(weight != nullptr && coverage != nullptr);
  for (uint8_t w = 0; w < 2; ++w)
    for (uint8_t c = 0; c < 4; ++c) {
      *weight = w;
      *coverage = c;
      if (static_cast<In::Op::WeightMode>(w) == In::Op::WeightMode::NONE)
        run_sample_coverage_cases<SourceP, PB::Weight::None>(program, ctx, c);
      else
        run_sample_coverage_cases<SourceP, PB::Weight::Projection>(program, ctx,
                                                                   c);
    }
}

inline void test_shader_chain_parity_sample_twin_wave() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::TwinWaveSampleParams>(
        program, In::Op::SampleTwinWave::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::TwinWave<TwinWaveSampleMirror>>(program,
                                                                     ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_rings() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::RingsSampleParams>(
        program, In::Op::SampleRings::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Rings<RingsSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spherical_rings() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_spherical_sample_op_chain<In::Op::SphericalRingsSampleParams>(
        program, In::Op::SampleSphericalRings::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    expect_spherical_sample_op_parity<
        PB::Source::SphericalRings<SphericalRingsSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spiral() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::SpiralSampleParams>(
        program, In::Op::SampleSpiral::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::Spiral<SpiralSampleMirror>>(program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_lattice() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::LatticeSampleParams>(
        program, In::Op::SampleLattice::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::PrimitiveLattice<LatticeSampleMirror>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_fractal() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::FractalSampleParams>(
        program, In::Op::SampleFractal::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    run_sample_op_matrix<PB::Source::EscapeFractal<FractalSampleMirror>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_tessellation() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::TessellationSampleParams>(
        program, In::Op::SampleTessellation::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    auto &params = param_as<In::Op::TessellationSampleParams>(program, 2);
    params.kind =
        static_cast<uint8_t>(PB::Source::TessellationKind::TRIANGULAR);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::TRIANGULAR>>(
        program, ctx);
    params.kind = static_cast<uint8_t>(PB::Source::TessellationKind::SQUARE);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::SQUARE>>(
        program, ctx);
    params.kind = static_cast<uint8_t>(PB::Source::TessellationKind::HEXAGONAL);
    run_sample_op_matrix<PB::Source::Tessellation<
        TessellationSampleMirror, PB::Source::TessellationKind::HEXAGONAL>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_projected_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_sample_op_chain<In::Op::ProjectedNoiseSampleParams>(
        program, In::Op::SampleProjectedNoise::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    uint8_t *basis = topology_byte(program, 2, "basis");
    HS_EXPECT_TRUE(basis != nullptr);
    *basis = static_cast<uint8_t>(::NoiseBasis::SIMPLEX);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    ::NoiseBasis::SIMPLEX>>(
        program, ctx);
    *basis = static_cast<uint8_t>(::NoiseBasis::FBM3);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    ::NoiseBasis::FBM3>>(
        program, ctx);
    *basis = static_cast<uint8_t>(::NoiseBasis::RIDGED3);
    run_sample_op_matrix<PB::Source::ProjectedNoise<ProjectedNoiseSampleMirror,
                                                    ::NoiseBasis::RIDGED3>>(
        program, ctx);
    program.clear();
  }
}

inline void test_shader_chain_parity_sample_spherical_noise() {
  for (const ValueSet set :
       {ValueSet::DEFAULTS, ValueSet::MINIMUMS, ValueSet::MAXIMUMS}) {
    auto fixture = std::make_unique<ProgramFixture>();
    In::ChainProgram &program = fixture->program;
    arm_spherical_sample_op_chain<In::Op::SphericalNoiseSampleParams>(
        program, In::Op::SampleSphericalNoise::ID, 4, set);
    const In::FrameContext ctx = shared_resources().context();
    expect_spherical_sample_op_parity<PB::Source::SphericalNoise<
        SphericalNoiseSampleMirror, ::NoiseBasis::SIMPLEX>>(program, ctx);
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

  // Warp::Vortex is deferred (no shipped FIELDS family), so the id stays
  // unknown as the operator set builds out.
  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.vortex.v2"},
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
  FastNoiseLite authored_walk_noise;
  authored_walk_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  authored_walk_noise.SetSeed(PB::CAMERA_WALK_SEED);
  authored_walk_noise.SetFrequency(In::Op::WALK_NOISE_SCALE);
  const auto &seeded_walk =
      state_as<In::Op::SpatialWalkState>(first->program, 0);
  HS_EXPECT_EQ(seeded_walk.walk_noise.GetNoise(0.25f, -0.5f, 0.75f),
               authored_walk_noise.GetNoise(0.25f, -0.5f, 0.75f));
  const In::FrameContext ctx = shared_resources().context();
  first->program.prepare(ctx);
  second->program.prepare(ctx);
  for (const Vector &view : sweep_views()) {
    const Color4 a = first->program.evaluate(view, ctx);
    const Color4 b = second->program.evaluate(view, ctx);
    HS_EXPECT_TRUE(color4_identical(a, b));
  }
  // Renaming does not change effect-authored walk resources.
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
  HS_EXPECT_EQ(std::memcmp(&walk_a.wander, &walk_b.wander, sizeof(Quaternion)),
               0);
  first->program.clear();
  second->program.clear();
  relabeled->program.clear();
}

inline void test_shader_chain_param_names_and_budget() {
  static_assert(In::PER_PARAM_NAME_BYTES ==
                In::MAX_INSTANCE_ID + 1 + In::MAX_FIELD_ID + 1);
  static_assert(In::operator_schema_ids_fit_names());
  auto fixture = std::make_unique<ProgramFixture>();
  In::ChainProgram &program = fixture->program;
  arm_default_chain(program, 0, ValueSet::DEFAULTS);

  // The names are exactly "{instance}.{field-id}", per schema entry.
  const auto ops = program.ops();
  for (size_t index = 0; index < ops.size(); ++index) {
    const In::OperatorDescriptor &op = *ops[index].op;
    for (uint16_t field = 0; field < op.schema_count; ++field) {
      std::string expected{ops[index].instance};
      expected += '.';
      expected += op.schema[field].id;
      HS_EXPECT_TRUE(expected == program.param_name(index, field));
    }
  }
  HS_EXPECT_TRUE(std::string_view(program.param_name(0, 0)) == "camera.wander");
  HS_EXPECT_TRUE(std::string_view(program.param_name(3, 0)) ==
                 "colorize.hue-shift-amount");

  // The committed footprint equals the accounted layout: cataloged blocks,
  // the per-op overhead, and the fixed per-field name reservation.
  size_t expected_bytes = 0;
  const auto align_to = [](size_t offset, size_t alignment) {
    return (offset + alignment - 1) & ~(alignment - 1);
  };
  for (const In::ChainProgram::ChainOp &op : ops) {
    const In::OperatorRuntime &runtime = op.op->runtime;
    expected_bytes = align_to(expected_bytes, runtime.param.align);
    expected_bytes += runtime.param.size;
    expected_bytes = align_to(expected_bytes, runtime.prepared.align);
    expected_bytes += runtime.prepared.size;
    expected_bytes = align_to(expected_bytes, runtime.state.align);
    expected_bytes += runtime.state.size;
    expected_bytes += In::PER_OP_OVERHEAD_BYTES;
    expected_bytes +=
        static_cast<size_t>(op.op->schema_count) * In::PER_PARAM_NAME_BYTES;
  }
  HS_EXPECT_EQ(program.used_bytes(), expected_bytes);
  program.clear();
}

/** Reaches the effect's committed program and palette state. */
struct ShaderChainWhiteBox {
  using FX = ShaderChain<96, 20>;

  static const In::ChainProgram &program(const FX &effect) {
    return effect.program;
  }
  static Pixel palette_color(const FX &effect, float value) {
    return effect.generated_palettes.palette(In::Op::PaletteMode::TRIADIC)
        .get(value)
        .color;
  }
};

inline void test_shader_chain_effect_registers_params() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  size_t expected = 0;
  for (const In::ChainEntryRequest &entry : DEFAULT_CHAIN)
    expected += In::find_operator(entry.operator_id)->schema_count;
  const ParamList &params = effect.getParameters();
  HS_EXPECT_EQ(params.size(), expected);
  HS_EXPECT_TRUE(params.find("camera.wander") != nullptr);
  HS_EXPECT_TRUE(params.find("project.singularity-fade") != nullptr);
  HS_EXPECT_TRUE(params.find("sample.pattern-freq") != nullptr);
  HS_EXPECT_TRUE(params.find("colorize.palette-chroma") != nullptr);
  const ParamDef *coverage = params.find("sample.coverage-mode");
  HS_EXPECT_TRUE(coverage != nullptr);
  HS_EXPECT_TRUE(coverage->is_enum());
  HS_EXPECT_EQ(coverage->option_count, 4);
  HS_EXPECT_TRUE(std::string_view(coverage->options[3]) == "edge-fade");
  // The value channel reaches the committed param block through the
  // registered target.
  HS_EXPECT_EQ(
      static_cast<int>(effect.updateParameter("sample.coverage-mode", 3.0f)),
      static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(params.find("sample.coverage-mode")->get_requested(), 3.0f);

  uint64_t lit = 0;
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}

inline void test_shader_chain_effect_rebind_generation() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t initial = effect.getParameterSchemaGeneration();
  const In::ChainEntryRequest edited[] = {
      {"camera2", "sphere.rotate.v2"},
      {"project", "project.stereographic.v2"},
      {"sample", "sample.grid.v2"},
      {"colorize", "colorize.generated-palette.v2"},
  };
  const In::ChainRefusal committed = effect.set_chain(edited);
  HS_EXPECT_EQ(static_cast<int>(committed.code),
               static_cast<int>(In::ChainStatus::OK));
  // set_chain rebinds before returning: fresh generation, fresh names.
  HS_EXPECT_NE(effect.getParameterSchemaGeneration(), initial);
  HS_EXPECT_TRUE(effect.getParameters().find("camera2.wander") != nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") == nullptr);
  effect.draw_frame();
  effect.advance_display();
}

inline void test_shader_chain_effect_refusal_keeps_schema() {
  reset_globals();
  ShaderChain<96, 20> effect;
  effect.init();
  const uint32_t committed = effect.getParameterSchemaGeneration();
  const size_t param_count = effect.getParameters().size();
  const In::ChainEntryRequest unknown[] = {
      {"camera", "sphere.rotate.v2"},
      {"warp", "warp.vortex.v2"},
  };
  const In::ChainRefusal refusal = effect.set_chain(unknown);
  HS_EXPECT_EQ(static_cast<int>(refusal.code),
               static_cast<int>(In::ChainStatus::UNKNOWN_OPERATOR));
  HS_EXPECT_EQ(refusal.entry_index, 1);
  // Transactional at the effect layer too: definitions, generation, and the
  // committed program all survive, and the effect still renders.
  HS_EXPECT_EQ(effect.getParameterSchemaGeneration(), committed);
  HS_EXPECT_EQ(effect.getParameters().size(), param_count);
  HS_EXPECT_TRUE(effect.getParameters().find("camera.wander") != nullptr);
  uint64_t lit = 0;
  for (int frame = 0; frame < 2; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
}

/** @brief Pause gates preset selection only: chain clocks and the visible
    palette keep advancing. */
inline void test_shader_chain_pause_semantics() {
  using WB = ShaderChainWhiteBox;
  reset_globals();
  WB::FX effect;
  effect.init();
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("sample.speed", 0.05f)),
               static_cast<int>(ParamSetResult::APPLIED));
  HS_EXPECT_EQ(static_cast<int>(effect.updateParameter("camera.wander", 1.0f)),
               static_cast<int>(ParamSetResult::APPLIED));
  // An animated-param write engages the pause by itself.
  effect.setAnimationsPaused(true);
  const In::Op::SourceClockState source_before =
      state_as<In::Op::SourceClockState>(WB::program(effect), 2);
  const Quaternion wander_before =
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander;
  const Pixel color_before = WB::palette_color(effect, 0.25f);

  for (int frame = 0; frame < 60; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }

  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_NE(
      state_as<In::Op::SourceClockState>(WB::program(effect), 2).primary,
      source_before.primary);
  HS_EXPECT_TRUE(
      state_as<In::Op::SpatialWalkState>(WB::program(effect), 0).wander !=
      wander_before);
  const Pixel color_after = WB::palette_color(effect, 0.25f);
  HS_EXPECT_TRUE(color_after.r != color_before.r ||
                 color_after.g != color_before.g ||
                 color_after.b != color_before.b);
  uint64_t lit = 0;
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      lit += static_cast<uint64_t>(pixel.r) + pixel.g + pixel.b;
    }
  HS_EXPECT_GT(lit, 0u);
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
  test_shader_chain_parity_displace_curl();
  test_shader_chain_parity_displace_direct();
  test_shader_chain_parity_displace_ripple();
  test_shader_chain_parity_lens_ops();
  test_shader_chain_parity_project_ops();
  test_shader_chain_parity_project_hemispheres();
  test_shader_chain_parity_field_ops();
  test_shader_chain_parity_warp_affine_mirror();
  test_shader_chain_parity_warp_wave_shear();
  test_shader_chain_parity_warp_vector_noise();
  test_shader_chain_noise_instances_decorrelate();
  test_shader_chain_parity_warp_polar_chart();
  test_shader_chain_parity_warp_curl_flow();
  test_shader_chain_parity_sample_variants();
  test_shader_chain_parity_sample_twin_wave();
  test_shader_chain_parity_sample_rings();
  test_shader_chain_parity_sample_spherical_rings();
  test_shader_chain_parity_sample_spiral();
  test_shader_chain_parity_sample_lattice();
  test_shader_chain_parity_sample_fractal();
  test_shader_chain_parity_sample_tessellation();
  test_shader_chain_parity_sample_projected_noise();
  test_shader_chain_parity_sample_spherical_noise();
  test_shader_chain_parity_colorize_variants();
  test_shader_chain_refusal_shape();
  test_shader_chain_refusal_budget_overflows();
  test_shader_chain_refusal_migrate_failed();
  test_shader_chain_state_identity_migration();
  test_shader_chain_state_continuity_slice();
  test_shader_chain_determinism();
  test_shader_chain_param_names_and_budget();
  test_shader_chain_effect_registers_params();
  test_shader_chain_effect_rebind_generation();
  test_shader_chain_effect_refusal_keeps_schema();
  test_shader_chain_pause_semantics();
  return fixture.result();
}

} // namespace shader_chain_tests
} // namespace hs_test
