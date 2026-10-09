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

/** @brief Reads available bytes, or returns empty if the file cannot be opened. */
inline std::string read_document(const std::string &path) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  std::FILE *file = std::fopen(path.c_str(), "rb");
#pragma clang diagnostic pop
  if (file == nullptr)
    return {};
  std::string text;
  char buffer[4096];
  for (size_t n; (n = std::fread(buffer, 1, sizeof buffer, file)) > 0;)
    text.append(buffer, n);
  std::fclose(file);
  return text;
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
  UNKNOWN
};

/** @brief One document chain entry, classified by its operator id. */
struct DocumentSlot {
  std::string label;
  std::string operator_id;
  SlotRole role = SlotRole::UNKNOWN;
  int warp_side = 0; /**< 0 = outer warp family, 1 = inner warp family. */
};

inline bool derivation_value_reachable(std::string_view operator_id,
                                       std::string_view field_id,
                                       std::string_view value);

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
  return false;
}

/** @brief Writes one Mobius coefficient addressed by its chain field id. */
template <typename Params>
inline bool assign_mobius(Params &params, std::string_view id, float value) {
  if constexpr (Params::template HAS<"lens">) {
    auto &mobius = params.template get<"lens">().mobius;
    const struct {
      std::string_view id;
      float *slot;
    } coefficients[] = {
        {"mobius-a-re", &mobius.a.re}, {"mobius-a-im", &mobius.a.im},
        {"mobius-b-re", &mobius.b.re}, {"mobius-b-im", &mobius.b.im},
        {"mobius-c-re", &mobius.c.re}, {"mobius-c-im", &mobius.c.im},
        {"mobius-d-re", &mobius.d.re}, {"mobius-d-im", &mobius.d.im},
    };
    for (const auto &coefficient : coefficients)
      if (coefficient.id == id) {
        *coefficient.slot = value;
        return true;
      }
  }
  return false;
}

namespace ChainOp = Pullback::Interp::Op;

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
  if constexpr (std::is_same_v<Lens, Pullback::Lens::TetrahedralKaleidoscope>)
    return "tetrahedral";
  if constexpr (std::is_same_v<Lens, Pullback::Lens::OctahedralKaleidoscope>)
    return "octahedral";
  if constexpr (std::is_same_v<Lens, Pullback::Lens::DodecahedralKaleidoscope>)
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
 * against the effect's Spec, base-template arguments and DERIVATION_REACH.
 */
template <typename FX>
inline bool
apply_document_value(typename FX::Params &built, const DocumentSlot &slot,
                     std::string_view field_id, const JsonValue &value) {
  using Traits = TraitsOf<FX>;
  using Spec = typename Traits::Spec;
  if (value.kind == JsonValue::Kind::NUMBER) {
    const float number = static_cast<float>(value.number);
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
      } else {
        if constexpr (FX::Params::template HAS<"inner_warp">)
          return assign_field(built.template get<"inner_warp">(), field_id,
                              number);
      }
      return false;
    case SlotRole::SURFACE:
      if constexpr (FX::Params::template HAS<"surface">)
        return assign_field(built.template get<"surface">(), field_id, number);
      return false;
    case SlotRole::SAMPLE:
      if (assign_field(built.template get<"source">(), field_id, number))
        return true;
      if constexpr (FX::Params::template HAS<"value">)
        if (assign_field(built.template get<"value">(), field_id, number))
          return true;
      // The sample operators always carry an edge-width; without edge-fade
      // coverage it is inert in the chain and has no composed-effect field.
      return field_id == "edge-width" &&
             Spec::COVERAGE != Pullback::ProjectionCoverageMode::EDGE_FADE;
    case SlotRole::FIELD:
      if constexpr (FX::Params::template HAS<"value">)
        return assign_field(built.template get<"value">(), field_id, number);
      return false;
    case SlotRole::LENS:
      return assign_mobius(built, field_id, number);
    case SlotRole::COLORIZE:
      return assign_field(built.template get<"color">(), field_id, number);
    default:
      return false;
    }
  }
  if (value.kind != JsonValue::Kind::STRING)
    return false;
  const std::string_view text = value.text;
  switch (slot.role) {
  case SlotRole::COLORIZE:
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
    if (field_id == "hue-shift-mode") {
      HS_EXPECT_TRUE(
          text ==
          ChainOp::HUE_SHIFT_MODE_IDS[static_cast<uint8_t>(Traits::HUE)]);
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
      HS_EXPECT_TRUE(
          text ==
          ChainOp::COVERAGE_MODE_IDS[static_cast<uint8_t>(Spec::COVERAGE)]);
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
      HS_EXPECT_TRUE(Spec::PROJECTION ==
                         Pullback::ProjectionKind::GNOMONIC_FOLDED &&
                     text == "folded");
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  case SlotRole::LENS:
    if (field_id == "symmetry") {
      HS_EXPECT_EQ(text, lens_symmetry_id<typename Spec::LensPolicy>());
      HS_EXPECT_TRUE(
          derivation_value_reachable(slot.operator_id, field_id, text));
      return true;
    }
    return false;
  case SlotRole::SURFACE:
  case SlotRole::WARP:
    if (field_id == "basis" || field_id == "integrator" || field_id == "mode" ||
        field_id == "envelope" || field_id == "harmonic") {
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
 * Params through the engine's field tables and compared both ways. Reciprocal
 * lattice cell scale allows rounding error; other values compare bit-exactly.
 */
template <template <int, int> class E>
inline void check_document_values(const char *name) {
  using FX = E<SMALL_W, SMALL_H>;
  using Params = typename FX::Params;
  HS_CONTEXT(name);

  std::string file_name(FX::EFFECT_ID);
  for (char &c : file_name)
    if (c == '-')
      c = '_';
  const std::string path =
      std::string(HS_PROMOTED_PATTERNS_DIR "/") + file_name + ".shader.json";
  const std::string text = read_document(path);
  HS_EXPECT(!text.empty(), "promoted shader document is readable");
  if (text.empty())
    return;

  JsonParser parser{text};
  const JsonValue document = parser.parse_value();
  HS_EXPECT(!parser.failed, "promoted shader document parses");
  if (parser.failed)
    return;

  const JsonValue *effect_id = document.find("effect_id");
  HS_EXPECT_TRUE(effect_id != nullptr && effect_id->text == FX::EFFECT_ID);

  const JsonValue *descriptor = document.find("descriptor");
  const JsonValue *chain =
      descriptor != nullptr ? descriptor->find("chain") : nullptr;
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
    slot.operator_id = operator_id->text;
    slot.role = classify_operator(operator_id->text);
    HS_EXPECT(slot.role != SlotRole::UNKNOWN, "chain operator classified");
    if (slot.role == SlotRole::PROJECT) {
      using Spec = typename TraitsOf<FX>::Spec;
      const char *expected = nullptr;
      switch (Spec::PROJECTION) {
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
      // The v1 expansion keeps the warp's slot position in its label even
      // when the other slot's identity op is omitted from the chain.
      if (slot.label == "warp1")
        slot.warp_side = 0;
      else if (slot.label == "warp2")
        slot.warp_side = 1;
      else
        slot.warp_side = warps_seen;
      ++warps_seen;
    }
    if (slot.role == SlotRole::SURFACE)
      surface_at = static_cast<int>(slots.size());
    if (slot.role == SlotRole::LENS)
      lens_at = static_cast<int>(slots.size());
    slots.push_back(std::move(slot));
  }

  if (surface_at >= 0 && lens_at >= 0) {
    HS_EXPECT_EQ(surface_at > lens_at,
                 TraitsOf<FX>::SURFACE_PLACEMENT ==
                     Pullback::SurfacePlacement::AFTER_LENS);
  }

  const JsonValue *bank = document.find("preset_bank");
  const JsonValue *presets = bank != nullptr ? bank->find("presets") : nullptr;
  HS_EXPECT_TRUE(presets != nullptr);
  if (presets == nullptr)
    return;
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

  for (size_t index = 0; index < FX::PRESET_IDS.size(); ++index) {
    HS_CONTEXT("preset", static_cast<int>(index));
    const JsonValue *preset = nullptr;
    for (const JsonValue &candidate : presets->items) {
      const JsonValue *preset_id = candidate.find("preset_id");
      if (preset_id != nullptr && preset_id->text == FX::PRESET_IDS[index])
        preset = &candidate;
    }
    HS_EXPECT_TRUE(preset != nullptr);
    if (preset == nullptr)
      continue;
    const JsonValue *values = preset->find("values");
    HS_EXPECT_TRUE(values != nullptr);
    if (values == nullptr)
      continue;

    if constexpr (!FX::ANIMATED_PROJECTION)
      for (const DocumentSlot &slot : slots)
        if (slot.role == SlotRole::PROJECT) {
          const JsonValue *frame = values->find(slot.label + ".frame");
          HS_EXPECT_TRUE(frame != nullptr && frame->text == "identity");
        }

    Params built{};
    const float MISSING = std::bit_cast<float>(uint32_t{0x7fc00001});
    const auto poison = [MISSING]<typename Family>(Family &family) {
      if constexpr (Pullback::HasFields<Family> && !std::is_void_v<Family>)
        for (const auto &field : Family::FIELDS) {
          if (!gate_open<FX>(field.gate))
            continue;
          if (std::string_view(field.id) == "edge-width" &&
              FX::Spec::COVERAGE != Pullback::ProjectionCoverageMode::EDGE_FADE)
            continue;
          if constexpr (std::is_same_v<Family, Pullback::ColorParams>) {
            if ((field.member == &Family::brightness_bottom ||
                 field.member == &Family::brightness_top) &&
                TraitsOf<FX>::BRIGHTNESS ==
                    Pullback::Color::BrightnessEnvelope::NONE)
              continue;
            if (field.member == &Family::hue_shift_amount &&
                TraitsOf<FX>::HUE == Pullback::HueMode::NONE)
              continue;
            if ((field.member == &Family::hue_noise_scale ||
                 field.member == &Family::hue_noise_speed) &&
                TraitsOf<FX>::HUE != Pullback::HueMode::NOISE)
              continue;
          }
          family.*(field.member) = MISSING;
        }
    };
    built.visit([&]<typename Resource>(auto &family) {
      if constexpr (Pullback::HasFields<typename Resource::Family>)
        poison(family);
    });
    built.template get<"color">().palette_mapping =
        static_cast<Pullback::Color::PaletteMapping>(255);
    if constexpr (FX::Params::template HAS<"lens">)
      built.template get<"lens">().mobius = {MISSING, MISSING, MISSING,
                                             MISSING, MISSING, MISSING,
                                             MISSING, MISSING};
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
      if constexpr (requires {
                      built.template get<"source">().lattice_cell_scale;
                    }) {
        if (slot->role == SlotRole::WARP && field_id == "lattice-period") {
          HS_EXPECT_NEAR(static_cast<float>(value.number),
                         1.0f / preset_params_or_initial<FX>(index)
                                    .template get<"source">()
                                    .lattice_cell_scale,
                         1e-6f);
          built.template get<"source">().lattice_cell_scale =
              1.0f / static_cast<float>(value.number);
          continue;
        }
      }
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
    const Params expected = preset_params_or_initial<FX>(index);
    if constexpr (requires {
                    built.template get<"source">().lattice_cell_scale;
                  }) {
      HS_EXPECT_NEAR(built.template get<"source">().lattice_cell_scale,
                     expected.template get<"source">().lattice_cell_scale,
                     1e-6f);
      built.template get<"source">().lattice_cell_scale =
          expected.template get<"source">().lattice_cell_scale;
    }
    verify_params_equal(built, expected);
  }
}
