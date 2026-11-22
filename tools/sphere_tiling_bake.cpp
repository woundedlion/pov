/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <vector>

#include "core/render/pullback/sphere_tiling.h"

namespace {

constexpr struct BakeSpec {
  const char *name;
  const char *symbol;
  int steps;
  bool narrow = false;
} BAKES[] = {
    {"truncatedIcosahedron", "TRUNCATED_ICOSAHEDRON", 16},
    {"truncatedIcosahedron_truncate50d_ambo_dual",
     "TRUNCATED_ICOSAHEDRON_TRUNCATE50D_AMBO_DUAL", 20, true},
    {"dodecahedron_hk54_ambo_hk72", "DODECAHEDRON_HK54_AMBO_HK72", 24, true},
};

size_t solid_index(std::string_view name) {
  size_t index = 0;
  for (auto registry : Solids::all_registries())
    for (const auto &entry : registry) {
      if (entry.name == name)
        return index;
      ++index;
    }
  return index;
}

void vector_literal(FILE *output, const math::Vector &v) {
  std::fprintf(output, "{%af, %af, %af}", v.x, v.y, v.z);
}

void table_open(FILE *output, const char *type, const char *name,
                const char *suffix) {
  std::fprintf(output, "inline const %s %s_%s[] HS_PROGMEM_UNIQUE(%s_%s) = {\n",
               type, name, suffix, name, suffix);
}

template <class T>
void integer_table(FILE *output, std::span<const T> values, const char *type,
                   const char *name, const char *suffix) {
  table_open(output, type, name, suffix);
  for (size_t i = 0; i < values.size(); ++i) {
    if (i % 16 == 0)
      std::fputs("    ", output);
    std::fprintf(output, "%u,", static_cast<unsigned>(values[i]));
    std::fputc(i % 16 == 15 || i + 1 == values.size() ? '\n' : ' ', output);
  }
  std::fputs("};\n\n", output);
}

void emit_narrow(FILE *output, const char *symbol,
                 const Pullback::Source::SphereTilingEdgeView &view) {
  using namespace Pullback::Source;
  using namespace tiling_edges_detail;
  constexpr int STEPS = SphereTilingNarrowIndex::STEPS;
  constexpr float STEP = 2.0f / STEPS;
  constexpr float PAD = 1.0e-4f;
  if (view.edges.size() > SphereTilingNarrowIndex::INTERIOR) {
    std::fputs("narrow candidate IDs overflow\n", stderr);
    std::exit(1);
  }
  const float limit_sine = sinf(SphereTilingNarrowIndex::OUTER + PAD);
  std::vector<uint16_t> cells;
  cells.reserve(6 * STEPS * STEPS);
  std::vector<uint16_t> candidates;
  for (int face = 0; face < 6; ++face) {
    for (int y = 0; y < STEPS; ++y) {
      for (int x = 0; x < STEPS; ++x) {
        const float u = -1.0f + x * STEP;
        const float v = -1.0f + y * STEP;
        const math::Vector center = Pattern::sphere_face_direction(
            face, u + 0.5f * STEP, v + 0.5f * STEP);
        const math::Vector corners[] = {
            Pattern::sphere_face_direction(face, u, v),
            Pattern::sphere_face_direction(face, u + STEP, v),
            Pattern::sphere_face_direction(face, u, v + STEP),
            Pattern::sphere_face_direction(face, u + STEP, v + STEP)};
        float radius_chord = 0.0f;
        for (const auto &corner : corners)
          radius_chord =
              std::max(radius_chord, math::distance_squared(center, corner));
        const float limit_chord = chord_squared(
            chord_angle(radius_chord) + SphereTilingNarrowIndex::OUTER + PAD);
        const size_t begin = candidates.size();
        for (size_t id = 0; id < view.edges.size(); ++id) {
          const auto &edge = view.edges[id];
          if (edge_chord_squared(edge, center, view.vertices) > limit_chord)
            continue;
          bool front = true, back = true, interior = edge.inv_length > 0.0f;
          for (const auto &corner : corners) {
            const float sine = math::dot(edge.normal, corner);
            front &= sine > limit_sine;
            back &= sine < -limit_sine;
            interior &= math::dot(edge.from_a, corner) > PAD &&
                        math::dot(edge.from_b, corner) > PAD;
          }
          if (front || back)
            continue;
          candidates.push_back(static_cast<uint16_t>(
              id | (interior ? SphereTilingNarrowIndex::INTERIOR : 0)));
        }
        const size_t count = candidates.size() - begin;
        if (count == 0)
          cells.push_back(0);
        else if (count == 1) {
          cells.push_back(candidates.back() | SphereTilingNarrowIndex::LAST);
          candidates.pop_back();
        } else {
          if (begin + 1 >= SphereTilingNarrowIndex::LAST) {
            std::fputs("narrow candidate offsets overflow\n", stderr);
            std::exit(1);
          }
          cells.push_back(static_cast<uint16_t>(begin + 1));
          candidates.back() |= SphereTilingNarrowIndex::LAST;
        }
      }
    }
  }
  std::fprintf(stderr, "  narrow: %zu cells, %zu multi-edge candidates\n",
               cells.size(), candidates.size());
  integer_table(output, std::span<const uint16_t>(cells), "uint16_t", symbol,
                "NARROW_CELLS");
  integer_table(output, std::span<const uint16_t>(candidates), "uint16_t",
                symbol, "NARROW_CANDIDATES");
}

void emit(FILE *output, const BakeSpec &spec, size_t solid) {
  const Pullback::Source::SphereTilingGeometry geometry{solid};
  const auto index = geometry.build_edge_index(spec.steps);
  const auto view = geometry.edge_view(index);
  std::fprintf(stderr,
               "%s: %zu vertices, %zu edges, %zu cells, %zu candidates\n",
               spec.name, view.vertices.size(), view.edges.size(),
               view.offsets.size() - 1, view.candidate_ids.size());
  std::fprintf(output, "inline constexpr unsigned %s_SOLID = %zu;\n",
               spec.symbol, solid);
  std::fprintf(output, "inline constexpr int %s_STEPS = %d;\n\n", spec.symbol,
               spec.steps);

  table_open(output, "math::Vector", spec.symbol, "VERTICES");
  for (const auto &vertex : view.vertices) {
    std::fputs("    ", output);
    vector_literal(output, vertex);
    std::fputs(",\n", output);
  }
  std::fputs("};\n\n", output);

  table_open(output, "SphereTilingEdge", spec.symbol, "EDGES");
  for (const auto &edge : view.edges) {
    std::fprintf(output, "    {%u, %u, ", static_cast<unsigned>(edge.a),
                 static_cast<unsigned>(edge.b));
    vector_literal(output, edge.normal);
    std::fputs(", ", output);
    vector_literal(output, edge.from_a);
    std::fputs(", ", output);
    vector_literal(output, edge.from_b);
    std::fprintf(output, ", %af},\n", edge.inv_length);
  }
  std::fputs("};\n\n", output);
  integer_table(output, view.offsets, "uint32_t", spec.symbol, "OFFSETS");
  integer_table(output, view.candidate_ids, "uint16_t", spec.symbol,
                "CANDIDATES");
  if (spec.narrow)
    emit_narrow(output, spec.symbol, view);
}

void write_header(FILE *output) {
  std::fputs(R"(/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file sphere_tiling_edges_baked.h @brief Flash-resident sphere tiling arcs. */

#include "core/platform/attributes.h"
#include "core/render/pullback/sphere_tiling_edges.h"

namespace Pullback::Source {

/// @cond
namespace baked_tiling_detail {

// Generated by the sphere_tiling_bake target.
// clang-format off

)",
             output);
  for (const auto &spec : BAKES)
    emit(output, spec, solid_index(spec.name));
  std::fputs("// clang-format on\n\n} // namespace baked_tiling_detail\n"
             "/// @endcond\n\n",
             output);

  std::fputs(R"(/** @brief Checks whether a registry solid has baked edges.
 * @param solid Solid registry ordinal.
 * @return Whether an edge index is available. */
constexpr bool supports_baked_sphere_tiling_edges(unsigned solid) {
  using namespace baked_tiling_detail;
  switch (solid) {
)",
             output);
  for (const auto &spec : BAKES)
    std::fprintf(output, "  case %s_SOLID:\n", spec.symbol);
  std::fputs(R"(    return true;
  default:
    return false;
  }
}

/** @brief Selects the baked edge geometry for a registry solid.
 * @param solid Solid registry ordinal.
 * @return View of constant edge data, or an empty view. */
inline SphereTilingEdgeView baked_sphere_tiling_edges(unsigned solid) {
  using namespace baked_tiling_detail;
  switch (solid) {
)",
             output);
  for (const auto &spec : BAKES) {
    const auto *name = spec.symbol;
    std::fprintf(output, "  case %s_SOLID:\n", name);
    if (spec.narrow) {
      std::fprintf(output, "    return {%s_VERTICES,\n", name);
      std::fprintf(output, "            %s_EDGES,\n", name);
      std::fprintf(output, "            %s_OFFSETS,\n", name);
      std::fprintf(output, "            %s_CANDIDATES,\n", name);
      std::fprintf(output, "            %s_STEPS,\n", name);
      std::fprintf(output, "            {%s_NARROW_CELLS,\n", name);
      std::fprintf(output, "             %s_NARROW_CANDIDATES}};\n", name);
    } else {
      std::fprintf(output, "    return {%s_VERTICES, %s_EDGES,\n", name, name);
      std::fprintf(output, "            %s_OFFSETS, %s_CANDIDATES,\n", name,
                   name);
      std::fprintf(output, "            %s_STEPS};\n", name);
    }
  }
  std::fputs(R"(  default:
    return {};
  }
}

} // namespace Pullback::Source
)",
             output);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fputs("usage: sphere_tiling_bake <output-header>\n", stderr);
    return 1;
  }
  for (const auto &spec : BAKES) {
    if (solid_index(spec.name) == Solids::NUM_ENTRIES) {
      std::fprintf(stderr, "unknown solid: %s\n", spec.name);
      return 1;
    }
  }
#ifdef _MSC_VER
  FILE *output = nullptr;
  fopen_s(&output, argv[1], "wb");
#else
  FILE *output = std::fopen(argv[1], "wb");
#endif
  if (!output) {
    std::perror(argv[1]);
    return 1;
  }
  write_header(output);
  return std::fclose(output) == 0 ? 0 : 1;
}
