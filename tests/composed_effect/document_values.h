/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Promoted shader documents: JSON reader and document-value parity
// ============================================================================

/** @brief Minimal JSON value for reading the promoted shader documents. */
struct JsonValue {
  enum class Kind : uint8_t { NUL, BOOLEAN, NUMBER, STRING, ARRAY, OBJECT };
  Kind kind = Kind::NUL;
  bool boolean = false;
  double number = 0.0;
  std::string text;
  std::vector<JsonValue> items;
  std::vector<std::string> member_keys;
  std::vector<JsonValue> member_values;

  const JsonValue *find(std::string_view key) const {
    if (kind != Kind::OBJECT)
      return nullptr;
    for (size_t index = 0; index < member_keys.size(); ++index)
      if (member_keys[index] == key)
        return &member_values[index];
    return nullptr;
  }
};

/** @brief Fail-stop recursive-descent JSON parser over a document string. */
struct JsonParser {
  std::string_view source;
  size_t position = 0;
  bool failed = false;

  void skip_space() {
    while (position < source.size() &&
           (source[position] == ' ' || source[position] == '\t' ||
            source[position] == '\n' || source[position] == '\r'))
      ++position;
  }

  bool consume(char expected) {
    skip_space();
    if (position < source.size() && source[position] == expected) {
      ++position;
      return true;
    }
    failed = true;
    return false;
  }

  std::string parse_string() {
    std::string result;
    if (!consume('"'))
      return result;
    while (position < source.size() && source[position] != '"') {
      char c = source[position++];
      if (c != '\\') {
        result.push_back(c);
        continue;
      }
      if (position >= source.size())
        break;
      const char escape = source[position++];
      switch (escape) {
      case 'b':
        result.push_back('\b');
        break;
      case 'f':
        result.push_back('\f');
        break;
      case 'n':
        result.push_back('\n');
        break;
      case 'r':
        result.push_back('\r');
        break;
      case 't':
        result.push_back('\t');
        break;
      case 'u': {
        if (position + 4 > source.size()) {
          failed = true;
          return result;
        }
        unsigned code = 0;
        for (int digit = 0; digit < 4; ++digit) {
          const char h = source[position++];
          code <<= 4;
          if (h >= '0' && h <= '9')
            code |= static_cast<unsigned>(h - '0');
          else if (h >= 'a' && h <= 'f')
            code |= static_cast<unsigned>(h - 'a' + 10);
          else if (h >= 'A' && h <= 'F')
            code |= static_cast<unsigned>(h - 'A' + 10);
          else
            failed = true;
        }
        if (code < 0x80) {
          result.push_back(static_cast<char>(code));
        } else if (code < 0x800) {
          result.push_back(static_cast<char>(0xC0 | (code >> 6)));
          result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        } else {
          result.push_back(static_cast<char>(0xE0 | (code >> 12)));
          result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
          result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
        }
        break;
      }
      default:
        result.push_back(escape);
        break;
      }
    }
    if (!consume('"'))
      failed = true;
    return result;
  }

  JsonValue parse_value() {
    JsonValue value;
    skip_space();
    if (failed || position >= source.size()) {
      failed = true;
      return value;
    }
    const char c = source[position];
    if (c == '{') {
      ++position;
      value.kind = JsonValue::Kind::OBJECT;
      skip_space();
      if (position < source.size() && source[position] == '}') {
        ++position;
        return value;
      }
      do {
        skip_space();
        std::string key = parse_string();
        if (!consume(':'))
          return value;
        JsonValue member = parse_value();
        value.member_keys.push_back(std::move(key));
        value.member_values.push_back(std::move(member));
        skip_space();
      } while (!failed && position < source.size() && source[position] == ',' &&
               ++position);
      consume('}');
      return value;
    }
    if (c == '[') {
      ++position;
      value.kind = JsonValue::Kind::ARRAY;
      skip_space();
      if (position < source.size() && source[position] == ']') {
        ++position;
        return value;
      }
      do {
        value.items.push_back(parse_value());
        skip_space();
      } while (!failed && position < source.size() && source[position] == ',' &&
               ++position);
      consume(']');
      return value;
    }
    if (c == '"') {
      value.kind = JsonValue::Kind::STRING;
      value.text = parse_string();
      return value;
    }
    if (source.substr(position, 4) == "true") {
      position += 4;
      value.kind = JsonValue::Kind::BOOLEAN;
      value.boolean = true;
      return value;
    }
    if (source.substr(position, 5) == "false") {
      position += 5;
      value.kind = JsonValue::Kind::BOOLEAN;
      return value;
    }
    if (source.substr(position, 4) == "null") {
      position += 4;
      return value;
    }
    value.kind = JsonValue::Kind::NUMBER;
    char *end = nullptr;
    const std::string number(source.substr(position, 64));
    value.number = std::strtod(number.c_str(), &end);
    if (end == number.c_str())
      failed = true;
    position += static_cast<size_t>(end - number.c_str());
    return value;
  }
};

template <typename FX>
inline std::optional<JsonValue> load_promoted_document() {
  std::string name(FX::EFFECT_ID);
  std::replace(name.begin(), name.end(), '-', '_');
  const std::string text = read_text_file(
      std::string(HS_PROMOTED_PATTERNS_DIR "/") + name + ".shader.json");
  HS_EXPECT(!text.empty(), "promoted shader document is readable");
  if (text.empty())
    return std::nullopt;
  JsonParser parser{text};
  JsonValue document = parser.parse_value();
  parser.skip_space();
  HS_EXPECT_EQ(parser.position, text.size());
  HS_EXPECT_FALSE(parser.failed);
  if (parser.failed || parser.position != text.size())
    return std::nullopt;
  return document;
}

/** @brief Returns a promoted document's `descriptor.chain`, or null. */
inline const JsonValue *document_chain(const JsonValue &document) {
  const JsonValue *descriptor = document.find("descriptor");
  return descriptor != nullptr ? descriptor->find("chain") : nullptr;
}

/** @brief Returns a promoted document's `preset_bank.presets`, or null. */
inline const JsonValue *document_presets(const JsonValue &document) {
  const JsonValue *bank = document.find("preset_bank");
  return bank != nullptr ? bank->find("presets") : nullptr;
}

using DocumentReachValue = std::array<std::string, 3>;

/** @brief Records shipped topology values, including omitted schema defaults. */
template <template <int, int> class E>
inline void record_document_reach(std::vector<DocumentReachValue> &observed) {
  using FX = E<SMALL_W, SMALL_H>;
  const auto loaded = load_promoted_document<FX>();
  if (!loaded)
    return;
  const auto *chain = document_chain(*loaded);
  const auto *presets = document_presets(*loaded);
  HS_EXPECT_TRUE(chain != nullptr && presets != nullptr);
  if (chain == nullptr || presets == nullptr)
    return;
  for (const auto &entry : chain->items) {
    const auto *label = entry.find("label");
    const auto *id = entry.find("operator");
    HS_EXPECT_TRUE(label != nullptr && id != nullptr);
    if (label == nullptr || id == nullptr)
      continue;
    const auto *op = In::find_operator(id->text.c_str());
    HS_EXPECT_TRUE(op != nullptr);
    if (op == nullptr)
      continue;
    for (const auto &preset : presets->items) {
      const auto *values = preset.find("values");
      HS_EXPECT_TRUE(values != nullptr);
      if (values == nullptr)
        continue;
      for (const auto &field : op->schema_span()) {
        if (!field.topology)
          continue;
        const auto *value = values->find(label->text + "." + field.id);
        const char *selected = field.enum_ids[field.enum_def];
        observed.push_back(
            {id->text, field.id, value != nullptr ? value->text : selected});
      }
    }
  }
}

/** @brief Every reachable table value is emitted by a shipped document. */
inline void test_composed_derivation_values_are_emitted() {
  std::vector<DocumentReachValue> observed;
#define HS_COMPOSED_REACH(name) record_document_reach<name>(observed);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_REACH)
#undef HS_COMPOSED_REACH
  for (const auto &row : DERIVATION_REACH) {
    if (row.topology_id == nullptr)
      continue;
    HS_CONTEXT(row.operator_id);
    HS_CONTEXT(row.topology_id);
    for (const char *value : row.reachable) {
      if (value == nullptr)
        break;
      HS_CONTEXT(value);
      HS_EXPECT_TRUE(std::find(observed.begin(), observed.end(),
                               DocumentReachValue{row.operator_id,
                                                  row.topology_id, value}) !=
                     observed.end());
    }
  }
}

/**
 * @brief Finds the values object of the last preset entry named @p preset_id.
 * @param presets The `preset_bank.presets` array.
 * @param preset_id Preset to look up.
 * @return The entry's `values` member, or null when absent.
 */
inline const JsonValue *preset_values(const JsonValue &presets,
                                      std::string_view preset_id) {
  const JsonValue *preset = nullptr;
  for (const JsonValue &candidate : presets.items) {
    const JsonValue *id = candidate.find("preset_id");
    if (id != nullptr && id->text == preset_id)
      preset = &candidate;
  }
  return preset != nullptr ? preset->find("values") : nullptr;
}

/** @brief Appended modulation schema defaults match canonical controls. */
inline void test_modulation_schema_defaults_match_controls() {
  using Controls = Pullback::Modulation::Controls;
  Controls controls;
  for (size_t index = 0; index < Controls::EXTRA_SCHEMA.size(); ++index) {
    const auto &field = Controls::EXTRA_SCHEMA[index];
    HS_CONTEXT(field.id);
    const void *address = Controls::address(controls.modulation, index);
    HS_EXPECT_TRUE(address != nullptr);
    if (address == nullptr)
      continue;
    if (field.topology)
      HS_EXPECT_EQ(*static_cast<const uint8_t *>(address), field.enum_def);
    else
      HS_EXPECT_NEAR(*static_cast<const float *>(address), field.def, 0.0f);
  }
}

inline void test_catalog_semantic_export() {
  namespace Catalog = Pullback::Interp;
  std::string text;
  Catalog::append_catalog_json(text);
  JsonParser parser{text};
  const JsonValue root = parser.parse_value();
  parser.skip_space();
  HS_EXPECT_FALSE(parser.failed);
  HS_EXPECT_EQ(parser.position, text.size());
  const auto member = [](const JsonValue &object,
                         const char *name) -> const JsonValue & {
    const JsonValue *value = object.find(name);
    HS_EXPECT(value != nullptr, "catalog member is present");
    static const JsonValue EMPTY;
    return value != nullptr ? *value : EMPTY;
  };
  const auto &operators = member(root, "operators").items;
  HS_EXPECT_EQ(operators.size(), Catalog::OPERATOR_TABLE.size());
  if (operators.size() != Catalog::OPERATOR_TABLE.size())
    return;
  for (size_t index = 0; index < operators.size(); ++index) {
    const auto &expected = Catalog::OPERATOR_TABLE[index];
    const auto &actual = operators[index];
    HS_CONTEXT(expected.operator_id);
    HS_EXPECT_TRUE(member(actual, "id").text == expected.operator_id);
    HS_EXPECT_TRUE(member(actual, "input").text ==
                   Catalog::CARRIER_NAMES[static_cast<size_t>(expected.input)]);
    HS_EXPECT_TRUE(
        member(actual, "output").text ==
        Catalog::CARRIER_NAMES[static_cast<size_t>(expected.output)]);
    const auto &params = member(actual, "params").items;
    HS_EXPECT_EQ(params.size(), expected.schema_count);
    if (params.size() != expected.schema_count)
      continue;
    for (size_t field = 0; field < params.size(); ++field) {
      const auto &schema = expected.schema[field];
      const auto &param = params[field];
      HS_CONTEXT(schema.id);
      HS_EXPECT_TRUE(member(param, "id").text == schema.id);
      if (schema.topology) {
        HS_EXPECT_TRUE(member(param, "topology").boolean);
        const auto &values = member(param, "values").items;
        HS_EXPECT_EQ(values.size(), schema.enum_count);
        for (size_t value = 0;
             value < values.size() && value < schema.enum_count; ++value)
          HS_EXPECT_TRUE(values[value].text == schema.enum_ids[value]);
        HS_EXPECT_TRUE(member(param, "default").text ==
                       schema.enum_ids[schema.enum_def]);
      } else {
        HS_EXPECT_EQ(static_cast<float>(member(param, "min").number),
                     schema.min);
        HS_EXPECT_EQ(static_cast<float>(member(param, "max").number),
                     schema.max);
        HS_EXPECT_EQ(static_cast<float>(member(param, "default").number),
                     schema.def);
      }
    }
  }
}

/** @brief The parameter family a document chain instance addresses. */
enum class SlotRole : uint8_t {
  CAMERA,
  LENS,
  SURFACE,
  PROJECT,
  WARP,
  SAMPLE,
  FIELD,
  COLORIZE,
  MODULATION,
  UNKNOWN
};

/** @brief One document chain entry, classified by its operator id. */
struct DocumentSlot {
  std::string label;
  std::string resource_key;
  std::string operator_id;
  SlotRole role = SlotRole::UNKNOWN;
  int warp_side = 0; /**< 0 = outer warp family, 1 = inner warp family,
                          2 = affine source-motion family. */
};

/** @brief Native resource addressed by a document chain entry. */
inline std::string_view document_resource_key(const DocumentSlot &slot) {
  if (!slot.resource_key.empty())
    return slot.resource_key;
  switch (slot.role) {
  case SlotRole::CAMERA:
  case SlotRole::PROJECT:
    return "projection";
  case SlotRole::WARP:
    return slot.warp_side == 0   ? "outer_warp"
           : slot.warp_side == 1 ? "inner_warp"
                                 : "affine";
  case SlotRole::SAMPLE:
    return "source";
  case SlotRole::FIELD:
    return "value";
  case SlotRole::COLORIZE:
    return "color";
  default:
    return slot.label;
  }
}

inline SlotRole classify_operator(std::string_view operator_id) {
  if (operator_id.starts_with("sphere.rotate."))
    return SlotRole::CAMERA;
  if (operator_id.starts_with("sphere.lens."))
    return SlotRole::LENS;
  if (operator_id.starts_with("sphere.displace."))
    return SlotRole::SURFACE;
  if (operator_id.starts_with("project."))
    return SlotRole::PROJECT;
  if (operator_id.starts_with("warp."))
    return SlotRole::WARP;
  if (operator_id.starts_with("sample."))
    return SlotRole::SAMPLE;
  if (operator_id.starts_with("field."))
    return SlotRole::FIELD;
  if (operator_id == "color.modulate.v1")
    return SlotRole::MODULATION;
  if (operator_id.starts_with("colorize."))
    return SlotRole::COLORIZE;
  return SlotRole::UNKNOWN;
}

/** @brief Writes @p value into the family field tabled under @p id, if any. */
template <typename Family>
inline bool assign_field(Family &family, std::string_view id, float value) {
  if constexpr (Pullback::HasFields<Family>) {
    for (const auto &field : Family::FIELDS)
      if (std::string_view(field.id) == id) {
        family.*(field.member) = value;
        return true;
      }
  }
  if constexpr (requires { family.extra_address(size_t{}); })
    for (size_t index = 0; index < Family::EXTRA_SCHEMA.size(); ++index)
      if (id == Family::EXTRA_SCHEMA[index].id &&
          !Family::EXTRA_SCHEMA[index].topology) {
        *static_cast<float *>(family.extra_address(index)) = value;
        return true;
      }
  if constexpr (requires { family.modulation; }) {
    Pullback::Interp::Op::ModulationParams bridge;
    bridge.modulation = family.modulation;
    for (size_t index = 0; index < bridge.EXTRA_SCHEMA.size(); ++index)
      if (id == bridge.EXTRA_SCHEMA[index].id &&
          !bridge.EXTRA_SCHEMA[index].topology) {
        *static_cast<float *>(bridge.extra_address(index)) = value;
        family.modulation = bridge.modulation;
        return true;
      }
  }
  return false;
}

/** @brief Assigns a generator field by chain ID. */
template <typename Params>
inline bool assign_mobius(Params &params, std::string_view id, float value) {
  if constexpr (std::is_same_v<typename Params::template Family<"lens">,
                               Pullback::MobiusGeneratorParams>)
    return assign_field(params.template get<"lens">(), id, value);
  return false;
}

namespace ChainOp = Pullback::Interp::Op;

/** @brief Affine parts one preset's values move away from identity.
 * @param values The preset's `values` object.
 * @param label Label of the affine chain instance.
 * @return ROTATES for a nonzero rotation rate, TRANSLATES for a nonzero
 * rounded winding, SCALES_SHEARS for a non-unit scale or nonzero shear. */
inline Pullback::AffineForm document_affine_form(const JsonValue &values,
                                                 const std::string &label) {
  const auto number = [&](const char *id, float identity) {
    const JsonValue *value = values.find(label + "." + id);
    if (value == nullptr)
      return identity;
    HS_EXPECT(value->kind == JsonValue::Kind::NUMBER,
              "affine value is a number");
    return static_cast<float>(value->number);
  };
  using Form = Pullback::AffineForm;
  Form form = Form::IDENTITY;
  if (number("rotation-rate", 0.0f) != 0.0f)
    form = form | Form::ROTATES;
  if (std::round(number("translation-x", 0.0f)) != 0.0f ||
      std::round(number("translation-y", 0.0f)) != 0.0f)
    form = form | Form::TRANSLATES;
  if (number("scale-x", 1.0f) != 1.0f || number("scale-y", 1.0f) != 1.0f ||
      number("shear", 0.0f) != 0.0f)
    form = form | Form::SCALES_SHEARS;
  return form;
}

/** @brief PALETTE_MODE_IDS entry for @p harmony, or null. */
constexpr const char *palette_mode_id(PaletteHarmony harmony) {
  switch (harmony) {
  case PaletteHarmony::TRIADIC:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::TRIADIC)];
  case PaletteHarmony::COMPLEMENTARY:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::COMPLEMENTARY)];
  case PaletteHarmony::ANALOGOUS:
    return ChainOp::PALETTE_MODE_IDS[static_cast<uint8_t>(
        ChainOp::PaletteMode::ANALOGOUS)];
  default:
    return nullptr;
  }
}

template <typename Lens> constexpr std::string_view lens_symmetry_id() {
  if constexpr (std::is_same_v<Lens, Pullback::Lens::Kaleidoscope>)
    return "azimuthal";
  if constexpr (std::is_base_of_v<Pullback::Lens::TetrahedralKaleidoscope,
                                  Lens>)
    return "tetrahedral";
  if constexpr (std::is_same_v<Lens, Pullback::Lens::OctahedralKaleidoscope>)
    return "octahedral";
  if constexpr (std::is_base_of_v<Pullback::Lens::DodecahedralKaleidoscope,
                                  Lens>)
    return "dodecahedral";
  if constexpr (std::is_same_v<Lens,
                               Pullback::Lens::TriangularPrismKaleidoscope>)
    return "triangular-prism";
  if constexpr (std::is_same_v<Lens, Pullback::Lens::SquarePrismKaleidoscope>)
    return "square-prism";
  if constexpr (std::is_same_v<Lens,
                               Pullback::Lens::PentagonalPrismKaleidoscope>)
    return "pentagonal-prism";
  if constexpr (std::is_same_v<Lens,
                               Pullback::Lens::HexagonalPrismKaleidoscope>)
    return "hexagonal-prism";
  if constexpr (std::is_same_v<Lens,
                               Pullback::Lens::OctagonalPrismKaleidoscope>)
    return "octagonal-prism";
  return {};
}

/**
 * @brief Applies one document preset entry onto @p built, or verifies it
 *        against the effect's compile-time constants.
 * @return False when the key addresses nothing this effect owns.
 * @details Numeric entries write through the owning family's field table, so a
 * document key the engine does not table fails loudly. String entries
 * are chain topology: the ones with a composed-effect equivalent are checked
 * against the derived topology and DERIVATION_REACH.
 */
template <typename FX>
inline bool
apply_document_value(typename FX::Params &built, const DocumentSlot &slot,
                     std::string_view field_id, const JsonValue &value) {
  using Traits = TraitsOf<FX>;
  if (value.kind == JsonValue::Kind::NUMBER) {
    const float number = static_cast<float>(value.number);
    bool assigned = false;
    built.visit([&]<typename Resource>(auto &family) {
      if (Resource::KEY.view() == document_resource_key(slot))
        assigned =
            assign_field(family,
                         slot.role == SlotRole::CAMERA &&
                                 Resource::KIND ==
                                     Pullback::ResourceKind::PROJECTION &&
                                 field_id == "wander"
                             ? "camera-wander"
                             : field_id,
                         number) ||
            assigned;
    });
    if (assigned)
      return true;
    switch (slot.role) {
    case SlotRole::CAMERA:
      if (field_id == "wander")
        return assign_field(built.template get<"projection">(), "camera-wander",
                            number);
      if (field_id == "spin-speed") {
        if constexpr (requires { FX::CAMERA_SPIN_RATE; }) {
          HS_EXPECT_EQ(bits(number), bits(FX::CAMERA_SPIN_RATE));
          return true;
        }
        return false;
      }
      return false;
    case SlotRole::PROJECT:
      return assign_field(built.template get<"projection">(), field_id, number);
    case SlotRole::WARP:
      if (slot.warp_side == 0) {
        if constexpr (FX::Params::template HAS<"outer_warp">)
          return assign_field(built.template get<"outer_warp">(), field_id,
                              number);
      } else if (slot.warp_side == 1) {
        if constexpr (FX::Params::template HAS<"inner_warp">)
          return assign_field(built.template get<"inner_warp">(), field_id,
                              number);
      } else {
        if constexpr (FX::Params::template HAS<"affine">)
          return assign_field(built.template get<"affine">(), field_id, number);
      }
      return false;
    case SlotRole::SURFACE: {
      bool assigned = false;
      built.visit([&]<typename Resource>(auto &family) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE)
          if (Resource::KEY.view() == slot.label)
            assigned = assign_field(family, field_id, number);
      });
      return assigned;
    }
    case SlotRole::SAMPLE:
      if (assign_field(built.template get<"source">(), field_id, number))
        return true;
      if constexpr (FX::Params::template HAS<"value">)
        if (assign_field(built.template get<"value">(), field_id, number))
          return true;
      // The sample operators always carry an edge-width; without edge-fade
      // coverage it is inert in the chain and has no composed-effect field.
      return field_id == "edge-width" &&
             FX::Metadata::COVERAGE !=
                 Pullback::ProjectionCoverageMode::EDGE_FADE;
    case SlotRole::FIELD:
      if constexpr (FX::Params::template HAS<"value">)
        return assign_field(built.template get<"value">(), field_id, number);
      return false;
    case SlotRole::LENS: {
      bool assigned = assign_mobius(built, field_id, number);
      built.visit([&]<typename Resource>(auto &family) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::LENS ||
                      Resource::KIND == Pullback::ResourceKind::TWIST)
          if (Resource::KEY.view() == slot.label)
            assigned = assign_field(family, field_id, number) || assigned;
      });
      return assigned;
    }
    case SlotRole::COLORIZE:
      return assign_field(built.template get<"color">(), field_id, number);
    case SlotRole::MODULATION: {
      bool assigned = false;
      built.visit([&]<typename Resource>(auto &family) {
        if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION)
          if (Resource::KEY.view() == slot.label)
            assigned = assign_field(family, field_id, number);
      });
      return assigned;
    }
    default:
      return false;
    }
  }
  if (value.kind != JsonValue::Kind::STRING)
    return false;
  const std::string_view text = value.text;
  if (slot.role != SlotRole::COLORIZE) {
    bool assigned = false;
    built.visit([&]<typename Resource>(auto &family) {
      using Family = typename Resource::Family;
      if constexpr (requires { Family::TOPOLOGY; }) {
        const bool matches =
            Resource::KEY.view() == document_resource_key(slot);
        if (!matches)
          return;
        for (const auto &field : Family::TOPOLOGY)
          if (field_id == field.id)
            for (const auto &schema :
                 In::find_operator(slot.operator_id.c_str())->schema_span())
              if (field_id == schema.id && schema.topology)
                for (uint8_t index = 0; index < schema.enum_count; ++index)
                  if (text == schema.enum_ids[index]) {
                    family.*field.member = index;
                    assigned = true;
                  }
      }
      if constexpr (requires { family.extra_address(size_t{}); }) {
        if (Resource::KEY.view() != document_resource_key(slot))
          return;
        for (size_t index = 0; index < Family::EXTRA_SCHEMA.size(); ++index) {
          const auto &field = Family::EXTRA_SCHEMA[index];
          if (field_id != field.id || !field.topology)
            continue;
          for (uint8_t option = 0; option < field.enum_count; ++option)
            if (text == field.enum_ids[option]) {
              *static_cast<uint8_t *>(family.extra_address(index)) = option;
              assigned = true;
            }
        }
      }
    });
    if (assigned)
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
    if (assigned)
      return true;
  }
  switch (slot.role) {
  case SlotRole::MODULATION: {
    bool assigned = false;
    built.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
        if (Resource::KEY.view() != slot.label)
          return;
        ChainOp::ModulationParams bridge;
        bridge.modulation = family.modulation;
        for (size_t index = 0; index < bridge.EXTRA_SCHEMA.size(); ++index) {
          const auto &field = bridge.EXTRA_SCHEMA[index];
          if (field_id != field.id || !field.topology)
            continue;
          for (uint8_t option = 0; option < field.enum_count; ++option)
            if (text == field.enum_ids[option]) {
              *static_cast<uint8_t *>(bridge.extra_address(index)) = option;
              family.modulation = bridge.modulation;
              assigned = true;
            }
        }
      }
    });
    if (assigned)
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
    return assigned;
  }
  case SlotRole::COLORIZE:
    if constexpr (FX::AUTHORED_PALETTE) {
      auto &color = built.template get<"color">();
      const auto assign = [&](uint8_t &member, const auto &ids) {
        for (uint8_t index = 0; index < std::size(ids); ++index)
          if (text == ids[index]) {
            HS_EXPECT_TRUE(
                derivation_value_reachable(slot.operator_id, field_id, text));
            member = index;
            return true;
          }
        return false;
      };
      if (field_id == "palette-mode")
        return assign(color.palette_mode, ChainOp::PALETTE_MODE_IDS);
      if (field_id == "palette-type")
        return assign(color.palette_type, AUTHORED_PALETTE_TYPE_IDS);
      if (field_id == "palette-harmony")
        return assign(color.palette_harmony, AUTHORED_PALETTE_HARMONY_IDS);
      if (field_id == "palette-lightness-curve")
        return assign(color.palette_lightness_curve,
                      AUTHORED_PALETTE_LIGHTNESS_CURVE_IDS);
      if (field_id == "palette-name")
        return assign(color.palette_name, AUTHORED_PALETTE_NAME_IDS);
      if (field_id == "palette-name-b")
        return assign(color.palette_name_b, AUTHORED_PALETTE_NAME_IDS);
      if (field_id == "brightness-envelope") {
        uint8_t index = 0;
        if (!assign(index, ChainOp::BRIGHTNESS_ENVELOPE_IDS))
          return false;
        color.brightness_envelope =
            static_cast<Pullback::Color::BrightnessEnvelope>(index);
        return true;
      }
    }
    if (field_id == "palette-type") {
      HS_EXPECT_TRUE(text == "cycling-generated");
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    if (field_id == "palette-mapping") {
      for (uint8_t index = 0; index < std::size(ChainOp::PALETTE_MAPPING_IDS);
           ++index)
        if (text == ChainOp::PALETTE_MAPPING_IDS[index]) {
          built.template get<"color">().palette_mapping =
              static_cast<Pullback::Color::PaletteMapping>(index);
          HS_EXPECT_TRUE(
              derivation_value_reachable(slot.operator_id, field_id, text));
          return true;
        }
      return false;
    }
    if (field_id == "palette-mode") {
      const char *expected = palette_mode_id(Traits::HARMONY);
      HS_EXPECT_TRUE(expected != nullptr && text == expected);
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    if (field_id == "brightness-envelope") {
      HS_EXPECT_TRUE(text ==
                     ChainOp::BRIGHTNESS_ENVELOPE_IDS[static_cast<uint8_t>(
                         Traits::BRIGHTNESS)]);
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  case SlotRole::SAMPLE:
    if (field_id == "coverage-mode") {
      if constexpr (static_cast<size_t>(FX::Metadata::COVERAGE) <
                    std::size(ChainOp::COVERAGE_MODE_IDS))
        HS_EXPECT_TRUE(text == ChainOp::COVERAGE_MODE_IDS[static_cast<uint8_t>(
                                   FX::Metadata::COVERAGE)]);
      else
        return false;
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    if (field_id == "weight-mode") {
      HS_EXPECT_TRUE(text == "projection");
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    if (field_id == "kind" || field_id == "fill") {
      using Source = typename FX::Params::template Family<"source">;
      if constexpr (requires {
                      Source::KIND;
                      Source::FILL;
                    }) {
        HS_EXPECT_TRUE(
            text == (field_id == "kind"
                         ? ChainOp::TESSELLATION_KIND_IDS[static_cast<uint8_t>(
                               Source::KIND)]
                         : ChainOp::TESSELLATION_FILL_IDS[static_cast<uint8_t>(
                               Source::FILL)]));
        HS_EXPECT_TRUE(
            derivation_value_reachable(slot.operator_id, field_id, text));
        return true;
      }
      return false;
    }
    return false;
  case SlotRole::PROJECT:
    if (field_id == "frame") {
      HS_EXPECT_TRUE(text ==
                     (FX::ANIMATED_PROJECTION ? "spin-wander" : "identity"));
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    if (field_id == "hemisphere") {
      HS_EXPECT_TRUE(FX::Metadata::PROJECTION ==
                         Pullback::ProjectionKind::GNOMONIC_FOLDED &&
                     text == "folded");
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  case SlotRole::LENS:
    if (field_id == "symmetry") {
      using Leaves = Pullback::ComposedDetail::PipelineLeaves<
          typename FX::Spec::template Pipeline<
              Pullback::ComposedDetail::DiscoveryBinding>>;
      const bool found = [&]<size_t... I>(std::index_sequence<I...>) {
        return (
            ([]<typename Stage>(std::string_view symmetry) {
              if constexpr (Pullback::ComposedDetail::IsLensStage<Stage>::value)
                return symmetry ==
                       lens_symmetry_id<typename Stage::LensPolicy>();
              return false;
            }.template
             operator()<typename Leaves::template StageAt<I>>(text)) ||
            ...);
      }(std::make_index_sequence<Leaves::STAGE_COUNT>{});
      HS_EXPECT_TRUE(found);
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  case SlotRole::SURFACE:
    if (field_id == "burst-size" || field_id == "random-center") {
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return text == (field_id == "burst-size" ? "four" : "on");
    }
    [[fallthrough]];
  case SlotRole::WARP:
    if (field_id == "path-length") {
      using Leaves = Pullback::ComposedDetail::PipelineLeaves<
          typename FX::Spec::template Pipeline<
              Pullback::ComposedDetail::DiscoveryBinding>>;
      const std::string_view key = document_resource_key(slot);
      const int tracked = [&]<size_t... I>(std::index_sequence<I...>) {
        int found = -1;
        (
            [&]<typename Stage>() {
              using Tracking = AffineStageTracking<Stage>;
              if constexpr (Tracking::AFFINE)
                if (Tracking::KEY == key)
                  found = Tracking::TRACK ? 1 : 0;
            }.template operator()<typename Leaves::template StageAt<I>>(),
            ...);
        return found;
      }(std::make_index_sequence<Leaves::STAGE_COUNT>{});
      HS_EXPECT(tracked >= 0, "path-length names a compiled affine stage");
      HS_EXPECT(text == (tracked == 1 ? "on" : "off"),
                "path-length matches the compiled affine TrackPath");
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return tracked >= 0;
    }
    if (field_id == "basis" || field_id == "integrator" ||
        field_id == "limit" || field_id == "mode" || field_id == "envelope" ||
        field_id == "harmonic") {
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  default:
    return false;
  }
}

/**
 * @brief Pins one effect's authored parameter values to its shader document.
 * @tparam E Composed effect class template.
 * @param name Effect name, for the failure context.
 * @details Every preset in the effect's shader document is rebuilt into a
 * Params through the engine's field tables and compared bit-exactly both
 * ways.
 */
template <template <int, int> class E>
inline void check_document_values(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  constexpr bool AUTHORED = [] {
    if constexpr (requires { FX::Spec::AUTHORED_TOPOLOGY; })
      return FX::Spec::AUTHORED_TOPOLOGY;
    return false;
  }();
  HS_CONTEXT(name);

  const auto loaded = load_promoted_document<FX>();
  if (!loaded)
    return;
  const JsonValue &document = *loaded;

  const JsonValue *effect_id = document.find("effect_id");
  HS_EXPECT_TRUE(effect_id != nullptr && effect_id->text == FX::EFFECT_ID);

  const JsonValue *chain = document_chain(document);
  HS_EXPECT_TRUE(chain != nullptr);
  if (chain == nullptr)
    return;
  std::vector<DocumentSlot> slots;
  int warps_seen = 0;
  int surface_at = -1, lens_at = -1;
  for (const JsonValue &entry : chain->items) {
    const JsonValue *label = entry.find("label");
    const JsonValue *operator_id = entry.find("operator");
    HS_EXPECT_TRUE(label != nullptr && operator_id != nullptr);
    if (label == nullptr || operator_id == nullptr)
      return;
    DocumentSlot slot;
    slot.label = label->text;
    if (const auto *resource = entry.find("resource_key"))
      slot.resource_key = resource->text;
    slot.operator_id = operator_id->text;
    slot.role = classify_operator(operator_id->text);
    if constexpr (AUTHORED)
      if (slot.role == SlotRole::CAMERA && slot.resource_key.empty())
        slot.resource_key = slot.label;
    HS_EXPECT(slot.role != SlotRole::UNKNOWN, "chain operator classified");
    if (slot.role == SlotRole::PROJECT) {
      const char *expected = nullptr;
      switch (FX::Metadata::PROJECTION) {
      case Pullback::ProjectionKind::GNOMONIC:
        expected = "project.gnomonic.v2";
        break;
      case Pullback::ProjectionKind::PEIRCE:
        expected = "project.peirce.v3";
        break;
      case Pullback::ProjectionKind::PEIRCE_SQUARE_FAST:
        expected = "project.peirce-square-fast.v3";
        break;
      case Pullback::ProjectionKind::BONNE:
        expected = "project.bonne.v3";
        break;
      case Pullback::ProjectionKind::AIROCEAN:
        expected = "project.airocean.v3";
        break;
      case Pullback::ProjectionKind::NONE:
        break;
      case Pullback::ProjectionKind::TETRAHEDRAL:
        expected = "project.tetrahedral.v1";
        break;
      case Pullback::ProjectionKind::CUBIC:
        expected = "project.cubic.v1";
        break;
      case Pullback::ProjectionKind::OCTAHEDRAL:
        expected = "project.octahedral.v1";
        break;
      case Pullback::ProjectionKind::DODECAHEDRAL:
        expected = "project.dodecahedral.v1";
        break;
      case Pullback::ProjectionKind::ICOSAHEDRAL:
        expected = "project.icosahedral.v1";
        break;
      case Pullback::ProjectionKind::STEREOGRAPHIC:
        expected = "project.stereographic.v2";
        break;
      case Pullback::ProjectionKind::GNOMONIC_FOLDED:
        expected = "project.gnomonic.v2";
        break;
      case Pullback::ProjectionKind::EQUIRECTANGULAR:
        expected = "project.equirectangular.v2";
        break;
      case Pullback::ProjectionKind::FOLDED_SINUSOIDAL:
        expected = "project.folded-sinusoidal.v2";
        break;
      }
      HS_EXPECT_TRUE(expected != nullptr && operator_id->text == expected);
    }
    if (slot.role == SlotRole::WARP) {
      // Labels warp1 and warp2 retain their slot when the other slot is absent.
      if (slot.label == "warp1")
        slot.warp_side = 0;
      else if (slot.label == "warp2")
        slot.warp_side = 1;
      else if (slot.label == "affine")
        slot.warp_side = 2;
      else
        slot.warp_side = warps_seen;
      ++warps_seen;
    }
    if (slot.role == SlotRole::SURFACE)
      surface_at = static_cast<int>(slots.size());
    if (slot.role == SlotRole::LENS && lens_at < 0)
      lens_at = static_cast<int>(slots.size());
    slots.push_back(std::move(slot));
  }

  if (surface_at >= 0 && lens_at >= 0) {
    HS_EXPECT_EQ(surface_at > lens_at,
                 TraitsOf<FX>::SURFACE_PLACEMENT ==
                     Pullback::SurfacePlacement::AFTER_LENS);
  }

  const JsonValue *presets = document_presets(document);
  HS_EXPECT_TRUE(presets != nullptr);
  if (presets == nullptr)
    return;
  const JsonValue *bank = document.find("preset_bank");
  HS_EXPECT_EQ(presets->items.size(), FX::PRESET_IDS.size());

  const JsonValue *choreography = bank->find("choreography");
  const JsonValue *order =
      choreography != nullptr ? choreography->find("generated_order") : nullptr;
  const JsonValue *dwell =
      choreography != nullptr ? choreography->find("dwell") : nullptr;
  HS_EXPECT_TRUE(order != nullptr && dwell != nullptr);
  if (order != nullptr && order->items.size() == FX::PRESET_IDS.size())
    for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index)
      HS_EXPECT_TRUE(order->items[index].text == FX::PRESET_IDS[index]);
  else
    HS_EXPECT_TRUE(order != nullptr &&
                   order->items.size() == FX::PRESET_IDS.size());
  if (dwell != nullptr) {
    HS_EXPECT_EQ(dwell->member_keys.size(), FX::PRESET_IDS.size());
    for (size_t member = 0; member < dwell->member_keys.size(); ++member) {
      HS_CONTEXT(dwell->member_keys[member].c_str());
      size_t preset = 0;
      while (preset < FX::PRESET_IDS.size() &&
             FX::PRESET_IDS[preset] != dwell->member_keys[member])
        ++preset;
      HS_EXPECT_LT(preset, FX::PRESET_IDS.size());
      HS_EXPECT_EQ(dwell->member_values[member].number,
                   double{FX::PRESET_DWELL_FRAMES});
    }
  }
  const JsonValue *edges = bank->find("edges");
  if constexpr (FX::PRESET_IDS.size() > 1)
    HS_EXPECT_TRUE(edges != nullptr &&
                   edges->items.size() == FX::PRESET_IDS.size());
  std::array<bool, FX::PRESET_IDS.size()> departures{};
  if (edges != nullptr)
    for (const JsonValue &edge : edges->items) {
      const JsonValue *from = edge.find("from");
      const JsonValue *to = edge.find("to");
      const JsonValue *duration = edge.find("duration");
      HS_EXPECT_TRUE(from != nullptr && duration != nullptr);
      if (from == nullptr || duration == nullptr)
        continue;
      size_t departing = 0;
      while (departing < FX::PRESET_IDS.size() &&
             FX::PRESET_IDS[departing] != from->text)
        ++departing;
      HS_EXPECT_TRUE(departing < FX::PRESET_IDS.size());
      if (departing < FX::PRESET_IDS.size()) {
        HS_EXPECT_FALSE(departures[departing]);
        departures[departing] = true;
        HS_EXPECT_TRUE(
            to != nullptr &&
            to->text ==
                FX::PRESET_IDS[(departing + 1) % FX::PRESET_IDS.size()]);
        HS_EXPECT_EQ(duration->number,
                     static_cast<double>(Segue::Preset::frames(
                         FX::preset_departure(departing))));
      }
    }

  if constexpr (FX::PRESET_IDS.size() > 1)
    for (bool seen : departures)
      HS_EXPECT_TRUE(seen);

  std::vector<std::pair<std::string, Pullback::AffineForm>> affine_forms;
  for (const DocumentSlot &slot : slots)
    if (slot.operator_id.starts_with("warp.affine."))
      affine_forms.emplace_back(std::string(document_resource_key(slot)),
                                Pullback::AffineForm::IDENTITY);

  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    const JsonValue *values = preset_values(*presets, FX::PRESET_IDS[index]);
    HS_EXPECT_TRUE(values != nullptr);
    if (values == nullptr)
      continue;
    for (auto &[key, form] : affine_forms)
      for (const DocumentSlot &slot : slots)
        if (slot.operator_id.starts_with("warp.affine.") &&
            document_resource_key(slot) == key)
          form = form | document_affine_form(*values, slot.label);

    if constexpr (!FX::ANIMATED_PROJECTION)
      for (const DocumentSlot &slot : slots)
        if (slot.role == SlotRole::PROJECT) {
          const JsonValue *frame = values->find(slot.label + ".frame");
          HS_EXPECT_TRUE(frame != nullptr && frame->text == "identity");
        }

    Params built{};
    const float MISSING = std::bit_cast<float>(uint32_t{0x7fc00001});
    const auto poison = [MISSING, values, &slots]<typename Resource>(
                            typename Resource::Family &family) {
      (void)values;
      (void)slots;
      using Family = typename Resource::Family;
      if constexpr (Pullback::HasFields<Family> && !std::is_void_v<Family>)
        for (const auto &field : Family::FIELDS) {
          if (!gate_open<FX>(field.gate))
            continue;
          if (AUTHORED || field.topology_gate.field != nullptr) {
            bool addressed = false;
            for (const auto &slot : slots)
              if (document_resource_key(slot) == Resource::KEY.view()) {
                const auto id =
                    slot.role == SlotRole::CAMERA &&
                            std::string_view(field.id) == "camera-wander"
                        ? "wander"
                        : field.id;
                addressed =
                    values->find(slot.label + "." + id) != nullptr || addressed;
              }
            if (!addressed)
              continue;
          }
          if (std::string_view(field.id) == "edge-width" &&
              !requires { family.coverage_mode; } &&
              FX::Metadata::COVERAGE !=
                  Pullback::ProjectionCoverageMode::EDGE_FADE)
            continue;
          if constexpr (std::is_same_v<Family, Pullback::ColorParams>) {
            if ((field.member == &Family::brightness_bottom ||
                 field.member == &Family::brightness_top) &&
                TraitsOf<FX>::BRIGHTNESS ==
                    Pullback::Color::BrightnessEnvelope::NONE)
              continue;
          }
          family.*(field.member) = MISSING;
        }
    };
    built.visit([&]<typename Resource>(auto &family) {
      if constexpr (Pullback::HasFields<typename Resource::Family>)
        poison.template operator()<Resource>(family);
    });
    built.visit([&]<typename Resource>(auto &family) {
      if constexpr (Resource::KIND == Pullback::ResourceKind::MODULATION) {
        family.modulation.mechanic.amount = MISSING;
        if (family.modulation.source == Pattern::SourceKind::NOISE) {
          family.modulation.noise.scale = MISSING;
          family.modulation.noise.speed = MISSING;
        }
      }
    });
    built.template get<"color">().palette_mapping =
        static_cast<Pullback::Color::PaletteMapping>(255);
    for (size_t member = 0; member < values->member_keys.size(); ++member) {
      const std::string &key = values->member_keys[member];
      const JsonValue &value = values->member_values[member];
      HS_CONTEXT(key.c_str());
      const size_t dot = key.find('.');
      HS_EXPECT_TRUE(dot != std::string::npos);
      if (dot == std::string::npos)
        continue;
      const std::string_view label = std::string_view(key).substr(0, dot);
      const std::string_view field_id = std::string_view(key).substr(dot + 1);
      const DocumentSlot *slot = nullptr;
      for (const DocumentSlot &candidate : slots)
        if (candidate.label == label)
          slot = &candidate;
      HS_EXPECT(slot != nullptr, "value key names a chain instance");
      if (slot == nullptr)
        continue;
      HS_EXPECT(apply_document_value<FX>(built, *slot, field_id, value),
                "document value maps onto the effect");
    }
    if constexpr (requires { FX::CAMERA_SPIN_RATE; }) {
      bool camera_spin_present = false;
      for (const DocumentSlot &slot : slots)
        if (slot.role == SlotRole::CAMERA)
          camera_spin_present =
              values->find(slot.label + ".spin-speed") != nullptr;
      HS_EXPECT_TRUE(camera_spin_present);
    }
    const Params expected = preset_params<FX>(index);
    verify_params_equal(built, expected);
  }

  using Leaves = Pullback::ComposedDetail::PipelineLeaves<
      typename FX::Spec::template Pipeline<
          Pullback::ComposedDetail::DiscoveryBinding>>;
  size_t compiled_affine = 0;
  [&]<size_t... I>(std::index_sequence<I...>) {
    (
        [&]<typename Stage>() {
          using Tracking = AffineStageTracking<Stage>;
          if constexpr (Tracking::AFFINE) {
            ++compiled_affine;
            HS_CONTEXT(Tracking::KEY.data());
            const auto document = std::find_if(
                affine_forms.begin(), affine_forms.end(),
                [](const auto &entry) { return entry.first == Tracking::KEY; });
            HS_EXPECT(document != affine_forms.end(),
                      "compiled affine stage has a document instance");
            if (document != affine_forms.end())
              HS_EXPECT(document->second == Tracking::FORM,
                        "compiled affine form is the document's minimal form");
          }
        }.template operator()<typename Leaves::template StageAt<I>>(),
        ...);
  }(std::make_index_sequence<Leaves::STAGE_COUNT>{});
  HS_EXPECT_EQ(compiled_affine, affine_forms.size());
}
