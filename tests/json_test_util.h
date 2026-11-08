/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Minimal JSON reader for test documents and catalogs.
 */
#pragma once

#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

namespace hs_test {

/** @brief Minimal parsed JSON value. */
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

} // namespace hs_test
