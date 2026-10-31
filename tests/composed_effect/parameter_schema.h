/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Parameter layout hash and schema pins
// ============================================================================

struct ParameterLayoutHash {
  uint64_t value = 14695981039346656037ULL;

  void number(uint64_t number) {
    for (unsigned i = 0; i < 8; ++i) {
      value = fnv1a64_byte(value, static_cast<uint8_t>(number));
      number >>= 8;
    }
  }

  void text(std::string_view text) {
    number(text.size());
    for (const unsigned char ch : text)
      value = fnv1a64_byte(value, ch);
  }

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
    } else if constexpr (Pullback::HasFields<Member>) {
      text("fields");
      number(Member::FIELDS.size());
      for (const auto &field : Member::FIELDS)
        this->member(root, member.*field.member, field.id);
      if constexpr (std::is_same_v<Member, Pullback::ColorParams>)
        this->member(root, member.palette_mapping, "palette-mapping");
      if constexpr (requires { member.modulation; })
        this->member(root, member.modulation, "modulation");
      if constexpr (requires { member.burst_size; }) {
        this->member(root, member.burst_size, "burst-size");
        this->member(root, member.random_center, "random-center");
      }
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
    } else if constexpr (std::is_same_v<Member, Pullback::MobiusLensParams>) {
      this->member(root, member.mobius, "mobius");
    } else if constexpr (std::is_same_v<Member, math::MobiusParams>) {
      this->member(root, member.a, "a");
      this->member(root, member.b, "b");
      this->member(root, member.c, "c");
      this->member(root, member.d, "d");
    } else if constexpr (std::is_same_v<Member, math::Complex>) {
      this->member(root, member.re, "re");
      this->member(root, member.im, "im");
    } else {
      static_assert(std::is_same_v<Member, void>, "unhashed parameter type");
    }
  }
};

struct SchemaFieldsAB {
  float a, b;
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SchemaFieldsAB>{"a", &SchemaFieldsAB::a, nullptr, 0, 1},
      Pullback::Field<SchemaFieldsAB>{"b", &SchemaFieldsAB::b, nullptr, 0, 1}};
};
struct SchemaFieldsBA {
  float b, a;
  static constexpr auto FIELDS = std::array{
      Pullback::Field<SchemaFieldsBA>{"a", &SchemaFieldsBA::a, nullptr, 0, 1},
      Pullback::Field<SchemaFieldsBA>{"b", &SchemaFieldsBA::b, nullptr, 0, 1}};
};

struct SchemaModulationHeaderAB {
  bool enabled;
  Pattern::SourceKind source;
};
struct SchemaModulationHeaderBA {
  Pattern::SourceKind source;
  bool enabled;
};

inline void test_parameter_layout_reorder() {
  const SchemaFieldsAB original{};
  const SchemaFieldsBA reordered{};
  static_assert(sizeof(original) == sizeof(reordered));
  ParameterLayoutHash before, after;
  before.member(original, original, "family");
  after.member(reordered, reordered, "family");
  HS_EXPECT_NE(before.value, after.value);
  const SchemaModulationHeaderAB original_header{};
  const SchemaModulationHeaderBA reordered_header{};
  static_assert(sizeof(original_header) == sizeof(reordered_header));
  ParameterLayoutHash header_before, header_after;
  header_before.member(original_header, original_header.enabled, "enabled");
  header_before.member(original_header, original_header.source, "source");
  header_after.member(reordered_header, reordered_header.enabled, "enabled");
  header_after.member(reordered_header, reordered_header.source, "source");
  HS_EXPECT_NE(header_before.value, header_after.value);
  ParameterLayoutHash boolean_type, enum_type;
  boolean_type.member(original_header, original_header.enabled, "tag");
  enum_type.member(reordered_header, reordered_header.source, "tag");
  HS_EXPECT_NE(boolean_type.value, enum_type.value);
}

/** @brief Hashes named parameter fields and their offsets in the snapshot. */
template <typename Params> uint64_t parameter_layout_hash() {
  const Params params{};
  ParameterLayoutHash hash;
  hash.number(sizeof(Params));
  hash.number(alignof(Params));
  params.visit([&]<typename Resource>(const auto &family) {
    hash.member(params, family, Resource::KEY.view());
  });
  return hash.value;
}

/** @brief Schema version, byte size and field layout of a persisted snapshot. */
struct ParameterSchemaPin {
  const char *effect;
  uint32_t schema_version;
  size_t params_bytes;
  uint64_t layout_hash;
};

constexpr ParameterSchemaPin PARAMETER_SCHEMA_PINS[] = {
    {"AlienBrain", 4, 256, 1520615164818923855ULL},
    {"KaleidoscopeHexSoft", 4, 252, 12248859291337775303ULL},
    {"AlienOcean", 4, 264, 10147967365089435384ULL},
    {"AlienCore", 4, 264, 10147967365089435384ULL},
    {"KaleidoscopeMandala", 4, 280, 13820729389425608030ULL},
    {"Millefiori", 1, 436, 1470632974023726218ULL},
    {"GridSpace", 4, 264, 11830402334610884471ULL},
    {"KaleidoscopeMelt", 3, 240, 738385364115786086ULL},
    {"ChromaticLichen", 4, 248, 3008862667822062567ULL},
    {"MermaidSkin", 4, 248, 3008862667822062567ULL},
    {"JewelMelt", 3, 240, 738385364115786086ULL},
    {"AshCloud", 5, 432, 556905902123095251ULL},
    {"KaleidoscopeHexOil", 4, 240, 8963401164026260972ULL},
    {"KaleidoscopeStainedGlass", 4, 280, 5435562790506759308ULL},
    {"KaleidoscopeSmooth", 6, 260, 16824253992083180273ULL},
    {"KaleidoscopeHexBright", 4, 252, 12248859291337775303ULL},
    {"KaleidoscopeFlowers", 4, 260, 16824253992083180273ULL},
    {"CosmicEyeball", 4, 264, 15015288303892513197ULL},
    {"MobiusGrid", 4, 284, 16410136890467484544ULL},

};

/** @brief Pins one specialization's schema version to its field layout. */
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
#define HS_COMPOSED_SCHEMA_PIN(name, seconds)                                  \
  check_parameter_schema_pin<name>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_SCHEMA_PIN)
#undef HS_COMPOSED_SCHEMA_PIN
#define HS_COMPOSED_ROSTER_ONE(name, seconds) +1
  constexpr size_t roster = 0 HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_ROSTER_ONE);
#undef HS_COMPOSED_ROSTER_ONE
  HS_EXPECT_EQ(std::size(PARAMETER_SCHEMA_PINS), roster);
}
