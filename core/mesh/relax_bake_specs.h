/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "mesh/relax_bake.h"

/** @brief Authored relax-bake names and iteration budgets. */
namespace Solids::RelaxBakeSpecs {
constexpr MeshOps::RelaxBake make_spec(const char *name, uint16_t iterations) {
  return {name, nullptr, 0, 0, 0, iterations, 0, 0, 0};
}
inline constexpr MeshOps::RelaxBake truncated_cuboctahedron_converged =
    make_spec("truncated_cuboctahedron_converged", 4096);
inline constexpr MeshOps::RelaxBake snub_cube_converged =
    make_spec("snub_cube_converged", 4096);
inline constexpr MeshOps::RelaxBake rhombicosidodecahedron_converged =
    make_spec("rhombicosidodecahedron_converged", 4096);
inline constexpr MeshOps::RelaxBake truncated_icosidodecahedron_converged =
    make_spec("truncated_icosidodecahedron_converged", 4096);
inline constexpr MeshOps::RelaxBake snub_dodecahedron_converged =
    make_spec("snub_dodecahedron_converged", 4096);
inline constexpr MeshOps::RelaxBake dodecahedron_ambo_bevel33_converged =
    make_spec("dodecahedron_ambo_bevel33_converged", 4096);
inline constexpr MeshOps::RelaxBake truncated_icosahedron_ambo_converged =
    make_spec("truncated_icosahedron_ambo_converged", 4096);
inline constexpr MeshOps::RelaxBake dodecahedron_bevel20_converged =
    make_spec("dodecahedron_bevel20_converged", 4096);
inline constexpr MeshOps::RelaxBake
    truncated_icosidodecahedron_bevel50_relax100 =
        make_spec("truncated_icosidodecahedron_bevel50_relax100", 100);
inline constexpr MeshOps::RelaxBake
    dodecahedron_hankin_ambo_hankin_ambo_converged =
        make_spec("dodecahedron_hankin_ambo_hankin_ambo_converged", 4096);
inline constexpr MeshOps::RelaxBake icosahedron_snub_converged =
    make_spec("icosahedron_snub_converged", 4096);
} // namespace Solids::RelaxBakeSpecs
