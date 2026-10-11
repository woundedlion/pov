/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_identity.h
 * @brief Canonical effect identities, compatibility aliases and RNG keys.
 */

#include <array>
#include <string_view>

namespace hs {
/** @brief Resolution-independent machine, display and compatibility identity. */
struct EffectIdentity {
  std::string_view machine_id;   ///< Canonical kebab-case machine ID.
  std::string_view display_name; ///< Class/browser name retained for display.
  std::string_view
      legacy_id; ///< Compatibility ID; equals machine_id when there is no distinct alias.
  std::string_view rng_seed_key; ///< Immutable key for the effect's RNG stream.

  /** @return Canonical, display and legacy lookup keys. */
  constexpr std::array<std::string_view, 3> accepted_ids() const {
    return {machine_id, display_name, legacy_id};
  }
  /**
   * @brief Matches an exact nonempty machine ID or compatibility alias.
   * @param candidate Lookup key.
   * @return True when the key names this effect.
   */
  constexpr bool matches(std::string_view candidate) const {
    return !candidate.empty() &&
           (candidate == machine_id || candidate == display_name ||
            candidate == legacy_id);
  }
};

/**
 * @brief Checks a lowercase, hyphen-separated machine ID.
 * @param id Machine ID.
 * @return True for nonempty alphanumeric words starting with a lowercase letter.
 */
constexpr bool canonical_effect_id(std::string_view id) {
  if (id.empty() || id.front() < 'a' || id.front() > 'z' || id.back() == '-')
    return false;
  char previous = 0;
  for (char value : id) {
    if (!(value >= 'a' && value <= 'z') && !(value >= '0' && value <= '9') &&
        value != '-')
      return false;
    if (value == '-' && previous == '-')
      return false;
    previous = value;
  }
  return true;
}

/**
 * @brief Returns an effect's canonical machine ID, or an unregistered class name.
 * @tparam EffectType Concrete effect type.
 * @param class_name Fallback class name.
 * @return Declared EFFECT_ID when present, else class_name.
 */
template <typename EffectType>
constexpr std::string_view stable_effect_id(std::string_view class_name) {
  if constexpr (requires { EffectType::EFFECT_ID; })
    return EffectType::EFFECT_ID;
  return class_name;
}

/**
 * @brief Returns the immutable RNG key for an effect type.
 * @tparam EffectType Concrete effect type.
 * @param class_name Fallback for an unregistered type.
 * @return The declared identity's RNG key, or the type's original ID fallback.
 */
template <typename EffectType>
constexpr std::string_view effect_seed_key(std::string_view class_name) {
  if constexpr (requires { EffectType::EFFECT_IDENTITY; })
    return EffectType::EFFECT_IDENTITY.rng_seed_key;
  return stable_effect_id<EffectType>(class_name);
}

namespace effect_identities {
/**
 * @brief Applies X to each implementation class and its identity fields.
 * @param X Macro invoked as X(ClassName, canonical_id, display_name, legacy_id, rng_key).
 */
#define HS_EFFECT_IDENTITY_LIST(X)                                             \
  X(AlienBrain, "alien-brain", "AlienBrain", "alien-brain", "alien-brain")     \
  X(AlienCore, "alien-core", "AlienCore", "alien-core", "alien-core")          \
  X(AlienOcean, "alien-ocean", "AlienOcean", "alien-ocean", "alien-ocean")     \
  X(AshCloud, "ash-cloud", "AshCloud", "ash-cloud", "ash-cloud")               \
  X(BZReactionDiffusion, "bz-reaction-diffusion", "BZReactionDiffusion",       \
    "BZReactionDiffusion", "BZReactionDiffusion")                              \
  X(ChromaticLichen, "chromatic-lichen", "ChromaticLichen",                    \
    "chromatic-lichen", "chromatic-lichen")                                    \
  X(Comets, "comets", "Comets", "Comets", "Comets")                            \
  X(CosmicEyeball, "cosmic-eyeball", "CosmicEyeball", "cosmic-eyeball",        \
    "cosmic-eyeball")                                                          \
  X(DisplacementField, "displacement-field", "DisplacementField",              \
    "DisplacementField", "DisplacementField")                                  \
  X(DreamBalls, "dream-balls", "DreamBalls", "DreamBalls", "DreamBalls")       \
  X(Dynamo, "dynamo", "Dynamo", "Dynamo", "Dynamo")                            \
  X(Fishbowl, "fishbowl", "Fishbowl", "Fishbowl", "Fishbowl")                  \
  X(Galaxies, "galaxies", "Galaxies", "Galaxies", "Galaxies")                  \
  X(GnomonicStars, "gnomonic-stars", "GnomonicStars", "GnomonicStars",         \
    "GnomonicStars")                                                           \
  X(GridSpace, "grid-space", "GridSpace", "grid-space", "grid-space")          \
  X(GSReactionDiffusion, "gs-reaction-diffusion", "GSReactionDiffusion",       \
    "GSReactionDiffusion", "GSReactionDiffusion")                              \
  X(HankinSolids, "hankin-solids", "HankinSolids", "HankinSolids",             \
    "HankinSolids")                                                            \
  X(HopfFibration, "hopf-fibration", "HopfFibration", "HopfFibration",         \
    "HopfFibration")                                                           \
  X(HyperLattice, "hyper-lattice", "HyperLattice", "HyperLattice",             \
    "HyperLattice")                                                            \
  X(IslamicStars, "islamic-stars", "IslamicStars", "IslamicStars",             \
    "IslamicStars")                                                            \
  X(JewelMelt, "jewel-melt", "JewelMelt", "jewel-melt", "jewel-melt")          \
  X(KaleidoscopeFlowers, "kaleidoscope-flowers", "KaleidoscopeFlowers",        \
    "kaleidoscope-flowers", "kaleidoscope-flowers")                            \
  X(KaleidoscopeHexBright, "kaleidoscope-hex-bright", "KaleidoscopeHexBright", \
    "kaleidoscope-hex-bright", "kaleidoscope-hex-bright")                      \
  X(KaleidoscopeHexOil, "kaleidoscope-hex-oil", "KaleidoscopeHexOil",          \
    "kaleidoscope-hex-oil", "kaleidoscope-hex-oil")                            \
  X(KaleidoscopeHexSoft, "kaleidoscope-hex-soft", "KaleidoscopeHexSoft",       \
    "kaleidoscope-hex-soft", "kaleidoscope-hex-soft")                          \
  X(KaleidoscopeMandala, "kaleidoscope-mandala", "KaleidoscopeMandala",        \
    "kaleidoscope-mandala", "kaleidoscope-mandala")                            \
  X(KaleidoscopePentBright, "kaleidoscope-pent-bright",                        \
    "KaleidoscopePentBright", "kaleidoscope-pent-bright",                      \
    "kaleidoscope-pent-bright")                                                \
  X(KaleidoscopeSmooth, "kaleidoscope-smooth", "KaleidoscopeSmooth",           \
    "kaleidoscope-smooth", "kaleidoscope-smooth")                              \
  X(KaleidoscopeStainedGlass, "kaleidoscope-stained-glass",                    \
    "KaleidoscopeStainedGlass", "kaleidoscope-stained-glass",                  \
    "kaleidoscope-stained-glass")                                              \
  X(LatticeMelt, "lattice-melt", "LatticeMelt", "lattice-melt",                \
    "lattice-melt")                                                            \
  X(MermaidSkin, "mermaid-skin", "MermaidSkin", "mermaid-skin",                \
    "mermaid-skin")                                                            \
  X(MeshFeedback, "mesh-feedback", "MeshFeedback", "MeshFeedback",             \
    "MeshFeedback")                                                            \
  X(MindSplatter, "mind-splatter", "MindSplatter", "MindSplatter",             \
    "MindSplatter")                                                            \
  X(MobiusGrid, "mobius-grid", "MobiusGrid", "mobius-grid", "mobius-grid")     \
  X(MobiusRings, "mobius-rings", "MobiusRings", "MobiusRings", "MobiusRings")  \
  X(PetalFlow, "petal-flow", "PetalFlow", "PetalFlow", "PetalFlow")            \
  X(Raymarch, "raymarch", "Raymarch", "Raymarch", "Raymarch")                  \
  X(RingShower, "ring-shower", "RingShower", "RingShower", "RingShower")       \
  X(RingSpin, "ring-spin", "RingSpin", "RingSpin", "RingSpin")                 \
  X(ShapeShifter, "shape-shifter", "ShapeShifter", "ShapeShifter",             \
    "ShapeShifter")                                                            \
  X(SphericalHarmonics, "spherical-harmonics", "SphericalHarmonics",           \
    "SphericalHarmonics", "SphericalHarmonics")                                \
  X(Thrusters, "thrusters", "Thrusters", "Thrusters", "Thrusters")             \
  X(Voronoi, "voronoi", "Voronoi", "Voronoi", "Voronoi")                       \
  X(ShaderChain, "shader-chain", "ShaderChain", "shader-chain", "shader-chain")

/**
 * @brief Declares one immutable identity and checks its machine ID.
 * @param Name Class/browser name.
 * @param Id Canonical machine ID.
 * @param Display Retained class/browser lookup name.
 * @param Legacy Previously persisted ID.
 * @param Seed Immutable RNG compatibility key.
 */
#define HS_DECLARE_EFFECT_IDENTITY(Name, Id, Display, Legacy, Seed)            \
  struct Name {                                                                \
    static constexpr EffectIdentity VALUE{Id, Display, Legacy, Seed};          \
  };                                                                           \
  static_assert(canonical_effect_id(Id),                                       \
                "effect machine IDs must be kebab-case");
HS_EFFECT_IDENTITY_LIST(HS_DECLARE_EFFECT_IDENTITY)
#undef HS_DECLARE_EFFECT_IDENTITY
} // namespace effect_identities
} // namespace hs
