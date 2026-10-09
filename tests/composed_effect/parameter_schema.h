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
      value = (value ^ (number & 0xff)) * 1099511628211ULL;
      number >>= 8;
    }
  }

  void text(std::string_view text) {
    number(text.size());
    for (const unsigned char ch : text)
      value = (value ^ ch) * 1099511628211ULL;
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

inline void test_parameter_layout_reorder() {
  const SchemaFieldsAB original{};
  const SchemaFieldsBA reordered{};
  static_assert(sizeof(original) == sizeof(reordered));
  ParameterLayoutHash before, after;
  before.member(original, original, "family");
  after.member(reordered, reordered, "family");
  HS_EXPECT_NE(before.value, after.value);
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
    {"AlienBrain", 2, 116, 2450429636986970733ULL},
    {"KaleidoscopeHexSoft", 2, 112, 8936720859510048520ULL},
    {"AlienOcean", 2, 124, 5732924260995353001ULL},
    {"AlienCore", 2, 124, 5732924260995353001ULL},
    {"KaleidoscopeMandala", 2, 140, 13066144578015210057ULL},
    {"GridSpace", 2, 124, 1602663564915579084ULL},
    {"LatticeMelt", 6, 100, 6646121068342335885ULL},
    {"ChromaticLichen", 2, 108, 149313388846503452ULL},
    {"MermaidSkin", 2, 108, 149313388846503452ULL},
    {"JewelMelt", 1, 100, 6646121068342335885ULL},
    {"AshCloud", 2, 108, 6209082325371804633ULL},
    {"KaleidoscopePentBright", 2, 124, 3683580502139334648ULL},
    {"KaleidoscopeHexOil", 2, 100, 9129498782342370022ULL},
    {"KaleidoscopeStainedGlass", 2, 140, 1167435404544472545ULL},
    {"KaleidoscopeSmooth", 4, 120, 8155564830961471120ULL},
    {"KaleidoscopeHexBright", 2, 112, 8936720859510048520ULL},
    {"KaleidoscopeFlowers", 2, 120, 8155564830961471120ULL},
    {"CosmicEyeball", 2, 124, 5732924260995353001ULL},
    {"MobiusGrid", 2, 144, 15460173752614586553ULL},
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
