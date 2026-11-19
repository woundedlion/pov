/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Parameter schema test cases and fixtures. */

#include "tests/composed_effect/context.h"

namespace hs_test::composed_effect_tests {

// ============================================================================
// Parameter layout hash and schema pins
// ============================================================================

/** @brief
 * FNV-1a accumulator over parameter family and member layout metadata.
 */
struct ParameterLayoutHash {
  /** @brief Current FNV-1a accumulator. */
  uint64_t value = 14695981039346656037ULL;

  /** @brief Hashes an integer in least-significant-byte order.
   * @param number Unsigned integer added to the layout hash.
   */
  void number(uint64_t number) {
    for (unsigned i = 0; i < 8; ++i) {
      value = fnv1a64_byte(value, static_cast<uint8_t>(number));
      number >>= 8;
    }
  }

  /** @brief Hashes the byte length and content of a string.
   * @param text Text bytes included in the layout hash.
   */
  void text(std::string_view text) {
    number(text.size());
    for (const unsigned char ch : text)
      value = fnv1a64_byte(value, ch);
  }

  /** @brief Hashes a named member, its type and its offset from the root.
   * @tparam Member Member type.
   * @tparam Root Containing parameter type.
   * @param root Containing parameter object.
   * @param member Member whose type and offset are hashed.
   * @param name Stable name of the parameter member.
   */
  template <typename Root, typename Member>
  void member(const Root &root, const Member &member, std::string_view name) {
    text(name);
    number(reinterpret_cast<uintptr_t>(&member) -
           reinterpret_cast<uintptr_t>(&root));
    number(sizeof(Member));
    number(alignof(Member));
    if constexpr (std::is_same_v<Member, float>) {
      text("float32");
    } else if constexpr (std::is_same_v<Member, bool>) {
      text("bool");
    } else if constexpr (std::is_same_v<Member, uint8_t>) {
      text("uint8");
    } else if constexpr (std::is_same_v<Member, uint16_t>) {
      text("uint16");
    } else if constexpr (std::is_same_v<Member, int32_t>) {
      text("int32");
    } else if constexpr (std::is_same_v<Member, Pattern::SourceKind>) {
      text("pattern-source-enum8");
    } else if constexpr (std::is_same_v<Member, Pullback::Source::PlanarKind>) {
      text("planar-source-enum8");
    } else if constexpr (std::is_same_v<Member, Pattern::Evaluation>) {
      text("pattern-evaluation-enum8");
    } else if constexpr (std::is_same_v<Member, Pattern::NoiseBasis>) {
      text("noise-basis-enum8");
    } else if constexpr (std::is_same_v<Member, Pattern::NoiseVariant>) {
      text("noise-variant-enum8");
    } else if constexpr (std::is_same_v<Member, ColorModulation::Mechanic>) {
      text("modulation-mechanic-enum8");
    } else if constexpr (std::is_same_v<Member,
                                        Pullback::Color::PaletteMapping>) {
      text("palette-mapping-enum8");
    } else if constexpr (std::is_same_v<Member,
                                        Pullback::Color::BrightnessEnvelope>) {
      text("brightness-envelope-enum8");
    } else if constexpr (requires {
                           typename Member::value_type;
                           member.size();
                           member.data();
                         }) {
      text("array");
      number(member.size());
      for (size_t index = 0; index < member.size(); ++index)
        this->member(root, member[index], std::to_string(index));
    } else if constexpr (requires {
                           member.type;
                           member.kind;
                           member.l;
                           member.m;
                           member.amplitude;
                           member.phase;
                         }) {
      text("potential-term");
      this->member(root, member.type, "type");
      this->member(root, member.kind, "kind");
      this->member(root, member.l, "l");
      this->member(root, member.m, "m");
      this->member(root, member.amplitude, "amplitude");
      this->member(root, member.phase, "phase");
    } else if constexpr (Pullback::HasFields<Member>) {
      text("fields");
      number(Member::FIELDS.size());
      for (const auto &field : Member::FIELDS)
        this->member(root, member.*field.member, field.id);
      visit_untabled_members(
          [&](std::string_view name, const auto &value) {
            this->member(root, value, name);
          },
          member);
    } else if constexpr (std::is_same_v<Member, ColorModulation::Config>) {
      text("modulation-config");
      this->member(root, member.enabled, "enabled");
      this->member(root, member.source, "source");
      this->member(root, member.version, "version");
      this->member(root, member.evaluation, "evaluation");
      this->member(root, member.noise, "noise");
      this->member(root, member.sine, "sine");
      this->member(root, member.interference, "interference");
      this->member(root, member.grid, "grid");
      this->member(root, member.seed, "seed");
      this->member(root, member.mechanic, "mechanic");
    } else if constexpr (std::is_same_v<Member, Pattern::NoiseParams>) {
      text("noise-params");
      this->member(root, member.scale, "scale");
      this->member(root, member.speed, "speed");
      this->member(root, member.basis, "basis");
      this->member(root, member.contrast, "contrast");
      this->member(root, member.variant, "variant");
    } else if constexpr (std::is_same_v<Member, Pattern::SineParams>) {
      text("sine-params");
      this->member(root, member.axis, "axis");
      this->member(root, member.frequency, "frequency");
      this->member(root, member.offset, "offset");
      this->member(root, member.speed, "speed");
    } else if constexpr (std::is_same_v<Member, Pattern::InterferenceParams> ||
                         std::is_same_v<Member, Pattern::GridParams>) {
      text(std::is_same_v<Member, Pattern::GridParams> ? "grid-params"
                                                       : "interference-params");
      this->member(root, static_cast<const Pattern::SineParams &>(member),
                   "primary");
      this->member(root, member.axis_b, "axis-b");
      this->member(root, member.frequency_b, "frequency-b");
      this->member(root, member.offset_b, "offset-b");
    } else if constexpr (std::is_same_v<Member, ColorModulation::Params>) {
      text("modulation-mechanic-params");
      this->member(root, member.kind, "kind");
      this->member(root, member.amount, "amount");
    } else if constexpr (std::is_same_v<Member, math::Vector>) {
      text("vector3");
      this->member(root, member.x, "x");
      this->member(root, member.y, "y");
      this->member(root, member.z, "z");

    } else {
      static_assert(std::is_same_v<Member, void>, "unhashed parameter type");
    }
  }
};

/** @brief Numeric layout fixture storing a before b. */
struct SchemaFieldsAB {
  /** @brief Numeric fields whose member offsets enter the layout hash. */
  float a, /** @brief Numeric field b in the layout fixture. */ b;
  /** @brief Named fields included in the layout hash. */
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SchemaFieldsAB>{"a", &SchemaFieldsAB::a, nullptr, 0, 1},
      Pullback::Field<SchemaFieldsAB>{"b", &SchemaFieldsAB::b, nullptr, 0, 1}};
};
/** @brief Numeric layout fixture storing b before a. */
struct SchemaFieldsBA {
  /** @brief Numeric fields whose member offsets enter the layout hash. */
  float b, /** @brief Numeric field a in the reordered fixture. */ a;
  /** @brief Named fields included in the layout hash. */
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SchemaFieldsBA>{"a", &SchemaFieldsBA::a, nullptr, 0, 1},
      Pullback::Field<SchemaFieldsBA>{"b", &SchemaFieldsBA::b, nullptr, 0, 1}};
};

/** @brief Modulation header with the enable flag before the selector. */
struct SchemaModulationHeaderAB {
  /** @brief Whether the modulation fixture is enabled. */
  bool enabled;
  /** @brief Modulation source selector. */
  Pattern::SourceKind source;
};
/** @brief Modulation header with the selector before the enable flag. */
struct SchemaModulationHeaderBA {
  /** @brief Modulation source selector. */
  Pattern::SourceKind source;
  /** @brief Whether the modulation fixture is enabled. */
  bool enabled;
};

/** @brief
 * Field reorders and type changes affect the hash even at equal layout size.
 */
inline void test_parameter_layout_reorder() {
  const SchemaFieldsAB ORIGINAL{};
  const SchemaFieldsBA REORDERED{};
  static_assert(sizeof(ORIGINAL) == sizeof(REORDERED));
  ParameterLayoutHash before, after;
  before.member(ORIGINAL, ORIGINAL, "family");
  after.member(REORDERED, REORDERED, "family");
  HS_EXPECT_NE(before.value, after.value);
  const SchemaModulationHeaderAB ORIGINAL_HEADER{};
  const SchemaModulationHeaderBA REORDERED_HEADER{};
  static_assert(sizeof(ORIGINAL_HEADER) == sizeof(REORDERED_HEADER));
  ParameterLayoutHash header_before, header_after;
  header_before.member(ORIGINAL_HEADER, ORIGINAL_HEADER.enabled, "enabled");
  header_before.member(ORIGINAL_HEADER, ORIGINAL_HEADER.source, "source");
  header_after.member(REORDERED_HEADER, REORDERED_HEADER.enabled, "enabled");
  header_after.member(REORDERED_HEADER, REORDERED_HEADER.source, "source");
  HS_EXPECT_NE(header_before.value, header_after.value);
  ParameterLayoutHash boolean_type, enum_type;
  boolean_type.member(ORIGINAL_HEADER, ORIGINAL_HEADER.enabled, "tag");
  enum_type.member(REORDERED_HEADER, REORDERED_HEADER.source, "tag");
  HS_EXPECT_NE(boolean_type.value, enum_type.value);
}

/** @brief Hashes named parameter fields and their offsets in the snapshot.
 * @tparam Params Composed parameter block type.
 * @return Hash of the parameter block size and named member layout.
 */
template <typename Params> uint64_t parameter_layout_hash() {
  const Params PARAMS{};
  ParameterLayoutHash hash;
  hash.number(sizeof(Params));
  hash.number(alignof(Params));
  PARAMS.visit([&]<typename Resource>(const auto &family) {
    hash.member(PARAMS, family, Resource::KEY.view());
  });
  return hash.value;
}

/** @brief Schema version, byte size and field layout of a persisted snapshot. */
struct ParameterSchemaPin {
  /** @brief Effect identity of the schema pin. */
  const char *effect;
  /** @brief Version declared for the serialized parameter block. */
  uint32_t schema_version;
  /** @brief Byte size of the serialized parameter block. */
  size_t params_bytes;
  /** @brief Expected hash of named member types and offsets. */
  uint64_t layout_hash;
};

/** @brief Persisted schema versions and native parameter-layout hashes. */
constexpr ParameterSchemaPin PARAMETER_SCHEMA_PINS[] = {
    {"Involution", 5, 520, 14928048846491712631ULL},
    {"AlienBrain", 5, 284, 7611044476778441529ULL},
    {"AlienOcean", 5, 292, 15246186803803213796ULL},
    {"AlienCore", 5, 292, 15246186803803213796ULL},
    {"KaleidoscopeMandala", 5, 308, 14118180186280879034ULL},
    {"Millefiori", 2, 436, 4151912069677071222ULL},
    {"GridSpace", 6, 272, 9211873771733847186ULL},
    {"GridWorms", 3, 252, 11193120348943479710ULL},
    {"GridWarp", 3, 480, 6265964079350224296ULL},
    {"LiquidWindows", 2, 296, 8838112055952055764ULL},
    {"Iridescence", 2, 516, 4473797754982781858ULL},
    {"RiverStones", 2, 668, 7241475773858249476ULL},
    {"KaleidoscopeMelt", 4, 240, 12469853473749147524ULL},
    {"ChromaticLichen", 5, 240, 2714246693381544909ULL},
    {"MermaidSkin", 5, 240, 2714246693381544909ULL},
    {"JewelMelt", 4, 240, 12469853473749147524ULL},
    {"AshCloud", 6, 432, 12146385198927097535ULL},
    {"KaleidoscopeHexOil", 6, 284, 1921882905749031131ULL},
    {"KaleidoscopeStainedGlass", 5, 308, 15089019455005349960ULL},
    {"KaleidoscopeSmooth", 7, 288, 607714572854722266ULL},
    {"KaleidoscopeHexBright", 4, 252, 12248859291337775303ULL},
    {"KaleidoscopeFlowers", 5, 288, 607714572854722266ULL},
    {"CosmicEyeball", 5, 292, 13478471776509544105ULL},
    {"MobiusGrid", 5, 276, 15629938187438781434ULL},

};

/** @brief Pins one specialization's schema version to its field layout.
 * @tparam E Effect template instantiated for the schema check.
 * @param name Effect identity used to locate the schema pin.
 */
template <template <int, int> class E>
inline void check_parameter_schema_pin(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  HS_CONTEXT(name);
  const ParameterSchemaPin *pin = nullptr;
  for (const ParameterSchemaPin &row : PARAMETER_SCHEMA_PINS)
    if (std::string_view(row.effect) == name)
      pin = &row;
  HS_EXPECT_TRUE(pin != nullptr);
  if (pin == nullptr)
    return;
  HS_EXPECT_EQ(FX::PARAMETER_SCHEMA_VERSION, pin->schema_version);
  HS_EXPECT_EQ(sizeof(typename FX::Params), pin->params_bytes);
  const uint64_t layout = parameter_layout_hash<typename FX::Params>();
  HS_EXPECT_EQ(layout, pin->layout_hash);
}

/**
 * @brief Sweeps the schema-version pin over every specialization.
 * @details One row per roster entry, so an effect the roster gains reds here
 * until it is pinned.
 */
inline void test_composed_parameter_schema_pins() {
#define HS_COMPOSED_SCHEMA_PIN(name) check_parameter_schema_pin<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SCHEMA_PIN)
#undef HS_COMPOSED_SCHEMA_PIN
#define HS_COMPOSED_ROSTER_ONE(name) +1
  constexpr size_t ROSTER = 0 HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_ROSTER_ONE);
#undef HS_COMPOSED_ROSTER_ONE
  HS_EXPECT_EQ(std::size(PARAMETER_SCHEMA_PINS), ROSTER);
}

} // namespace hs_test::composed_effect_tests
