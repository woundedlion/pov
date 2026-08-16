/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <string_view>

#include "effects/ShaderBall.h"

#define HS_SHADER_DERIVED_EFFECT_LIST(X)                                       \
  X(Shader)                                                                    \
  X(SignalWeave)                                                               \
  X(KaleidoWave)                                                               \
  X(AlienOcean)                                                                \
  X(GlitchGrid)                                                                \
  X(FacetWave)                                                                 \
  X(ContourLattice)                                                            \
  X(PrismLattice)                                                              \
  X(VectorFacets)                                                              \
  X(HexWave)                                                                   \
  X(EquatorGrid)                                                               \
  X(CosmicEyeball)                                                             \
  X(MobiusGrid)

template <int W, int H> class Shader : public ShaderBall<W, H> {
public:
  static constexpr std::string_view EFFECT_ID = "shader";
};

#define HS_SINGLE_PROMOTED_SHADER_LOOK(ClassName, EffectId, PresetId, Index,   \
                                       Program)                                \
  template <int W, int H> class ClassName : public ShaderBall<W, H> {          \
  public:                                                                      \
    static constexpr std::string_view EFFECT_ID = EffectId;                    \
    static constexpr std::array<std::string_view, 1> PRESET_IDS{PresetId};     \
    static constexpr std::array<uint8_t, 1> SOURCE_PRESET_INDICES{Index};      \
    HS_COLD_MEMBER void init() override {                                      \
      this->set_fixed_preset_view(SOURCE_PRESET_INDICES);                      \
      ShaderBall<W, H>::init();                                                \
    }                                                                          \
    HS_FLASH_MEMBER void draw_frame() override {                               \
      if constexpr (Program >= 0)                                              \
        this->template draw_fixed_program_frame<Program>();                    \
      else                                                                     \
        ShaderBall<W, H>::draw_frame();                                        \
    }                                                                          \
  }

HS_SINGLE_PROMOTED_SHADER_LOOK(SignalWeave, "signal-weave", "signal-weave", 0,
                               0);
HS_SINGLE_PROMOTED_SHADER_LOOK(KaleidoWave, "kaleido-wave", "twin-wave", 1, 1);
HS_SINGLE_PROMOTED_SHADER_LOOK(AlienOcean, "alien-ocean", "folded-grid", 2, 2);
HS_SINGLE_PROMOTED_SHADER_LOOK(GlitchGrid, "glitch-grid", "folded-glitch", 3,
                               3);
HS_SINGLE_PROMOTED_SHADER_LOOK(FacetWave, "facet-wave", "wave-mirror", 5, 5);
HS_SINGLE_PROMOTED_SHADER_LOOK(ContourLattice, "contour-lattice",
                               "affine-contour", 6, 6);
HS_SINGLE_PROMOTED_SHADER_LOOK(PrismLattice, "prism-lattice", "polar-wave", 9,
                               8);
HS_SINGLE_PROMOTED_SHADER_LOOK(VectorFacets, "vector-facets", "vector-mirror",
                               10, 9);
HS_SINGLE_PROMOTED_SHADER_LOOK(HexWave, "hex-wave", "hex-twin-wave", 12, 11);
HS_SINGLE_PROMOTED_SHADER_LOOK(CosmicEyeball, "cosmic-eyeball", "mirrored-grid",
                               18, 13);

#undef HS_SINGLE_PROMOTED_SHADER_LOOK

template <int W, int H> class MobiusGrid : public ShaderBall<W, H> {
public:
  static constexpr std::string_view EFFECT_ID = "mobius-grid";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"mobius-grid"};
  static constexpr std::array<uint8_t, 1> SOURCE_PRESET_INDICES{19};

  HS_COLD_MEMBER void init() override {
    this->set_fixed_preset_view(SOURCE_PRESET_INDICES);
    ShaderBall<W, H>::init();
    this->start_circular_mobius_animation(1.0f, 160);
  }

  HS_FLASH_MEMBER void draw_frame() override {
    this->template draw_fixed_program_frame<14>();
  }
};

template <int W, int H> class EquatorGrid : public ShaderBall<W, H> {
public:
  static constexpr std::string_view EFFECT_ID = "equator-grid";
  static constexpr std::array<std::string_view, 3> PRESET_IDS{
      "double-map", "open-grid", "fine-grid"};
  static constexpr std::array<uint8_t, 3> SOURCE_PRESET_INDICES{15, 16, 17};

  HS_COLD_MEMBER void init() override {
    this->set_fixed_preset_view(SOURCE_PRESET_INDICES);
    ShaderBall<W, H>::init();
  }

  HS_FLASH_MEMBER void draw_frame() override {
    this->template draw_fixed_program_frame<12>();
  }
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(Shader)
REGISTER_EFFECT(SignalWeave)
REGISTER_EFFECT(KaleidoWave)
REGISTER_EFFECT(AlienOcean)
REGISTER_EFFECT(GlitchGrid)
REGISTER_EFFECT(FacetWave)
REGISTER_EFFECT(ContourLattice)
REGISTER_EFFECT(PrismLattice)
REGISTER_EFFECT(VectorFacets)
REGISTER_EFFECT(HexWave)
REGISTER_EFFECT(EquatorGrid)
REGISTER_EFFECT(CosmicEyeball)
REGISTER_EFFECT(MobiusGrid)
