/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_shader_chain.h.

// --- fixtures -------------------------------------------------------------

inline constexpr const char *MISMATCHED_TOPOLOGY_IDS[] = {"first", "second"};

struct MismatchedTopologyDefaults {
  uint8_t mode = 1;
  static constexpr auto TOPOLOGY = std::array{
      In::TopologyField<MismatchedTopologyDefaults>{
          "mode", &MismatchedTopologyDefaults::mode, MISMATCHED_TOPOLOGY_IDS,
          0},
  };
};

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
  BakedPaletteStorage palettes[3];
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
    ctx.projection_base = Pullback::projection_base_orientation();
    ctx.palettes = {&palettes[0].view(), &palettes[1].view(),
                    &palettes[2].view()};
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
    {"sample", "sample.grid.v3"},
    {"colorize", "colorize.generated-palette.v3"},
};

inline std::array<math::Vector, 14> sweep_views() {
  std::array<math::Vector, 14> views = {
      math::Vector(0, 1, 0),         math::Vector(0, -1, 0),
      math::Vector(1, 0, 0),         math::Vector(-1, 0, 0),
      math::Vector(0, 0, 1),         math::Vector(0, 0, -1),
      math::Vector(1, 1, 1),         math::Vector(-1, 1, -1),
      math::Vector(1, -2, 0.5f),     math::Vector(-0.3f, 0.9f, 0.6f),
      math::Vector(0.05f, 0.99f, 0), math::Vector(0.05f, -0.99f, 0),
      math::Vector(2, 0.1f, -1),     math::Vector(-1, -1, 2),
  };
  for (math::Vector &view : views)
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

/** @brief Selects defaults or a tabled parameter endpoint. */
enum class ValueSet { DEFAULTS, MINIMUMS, MAXIMUMS };

/** Failure-context label for a value set. */
inline const char *value_set_name(ValueSet set) {
  switch (set) {
  case ValueSet::MINIMUMS:
    return "minimums";
  case ValueSet::MAXIMUMS:
    return "maximums";
  default:
    return "defaults";
  }
}

/** @brief Applies tabled endpoints; DEFAULTS leaves the family unchanged. */
template <typename T> void apply_value_set(T &params, ValueSet set) {
  if (set == ValueSet::DEFAULTS)
    return;
  for (const auto &field : T::FIELDS)
    params.*(field.member) = set == ValueSet::MINIMUMS ? field.min : field.max;
}
