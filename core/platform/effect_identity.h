/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file effect_identity.h
 * @brief Canonical effect identities, compatibility aliases and RNG keys.
 */

#include "core/platform/rng.h"
#include <array>
#include <string_view>

namespace hs {
/** @brief Resolution-independent machine, display and compatibility identity. */
struct EffectIdentity {
  std::string_view machine_id;   ///< Canonical kebab-case machine ID.
  std::string_view display_name; ///< Class/browser name, also a lookup alias.
  std::string_view rng_seed_key; ///< Immutable key for the effect's RNG stream.

  /** @return Canonical and display lookup keys. */
  constexpr std::array<std::string_view, 2> accepted_ids() const {
    return {machine_id, display_name};
  }
  /**
   * @brief Matches an exact nonempty machine ID or display name.
   * @param candidate Lookup key.
   * @return True when the key names this effect.
   */
  constexpr bool matches(std::string_view candidate) const {
    return !candidate.empty() &&
           (candidate == machine_id || candidate == display_name);
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
 * @param X Macro invoked as X(ClassName, canonical_id, rng_key).
 */
#define HS_EFFECT_IDENTITY_LIST(X)                                             \
  X(AlienBrain, "alien-brain", "alien-brain")                                  \
  X(AlienCore, "alien-core", "alien-core")                                     \
  X(AlienOcean, "alien-ocean", "alien-ocean")                                  \
  X(AshCloud, "ash-cloud", "ash-cloud")                                        \
  X(BZReactionDiffusion, "bz-reaction-diffusion", "BZReactionDiffusion")       \
  X(ChromaticLichen, "chromatic-lichen", "chromatic-lichen")                   \
  X(Comets, "comets", "Comets")                                                \
  X(CosmicEyeball, "cosmic-eyeball", "cosmic-eyeball")                         \
  X(DisplacementField, "displacement-field", "DisplacementField")              \
  X(DreamBalls, "dream-balls", "DreamBalls")                                   \
  X(Dynamo, "dynamo", "Dynamo")                                                \
  X(Fishbowl, "fishbowl", "Fishbowl")                                          \
  X(Galaxies, "galaxies", "Galaxies")                                          \
  X(GnomonicStars, "gnomonic-stars", "GnomonicStars")                          \
  X(GridSpace, "grid-space", "grid-space")                                     \
  X(GridWarp, "grid-warp", "grid-warp")                                        \
  X(GridWorms, "grid-worms", "grid-worms")                                     \
  X(GSReactionDiffusion, "gs-reaction-diffusion", "GSReactionDiffusion")       \
  X(HankinSolids, "hankin-solids", "HankinSolids")                             \
  X(HopfFibration, "hopf-fibration", "HopfFibration")                          \
  X(HyperLattice, "hyper-lattice", "HyperLattice")                             \
  X(Involution, "involution", "involution")                                    \
  X(Iridescence, "iridescence", "iridescence")                                 \
  X(IslamicStars, "islamic-stars", "IslamicStars")                             \
  X(JewelMelt, "jewel-melt", "jewel-melt")                                     \
  X(KaleidoscopeFlowers, "kaleidoscope-flowers", "kaleidoscope-flowers")       \
  X(KaleidoscopeHexBright, "kaleidoscope-hex-bright",                          \
    "kaleidoscope-hex-bright")                                                 \
  X(KaleidoscopeHexOil, "kaleidoscope-hex-oil", "kaleidoscope-hex-oil")        \
  X(KaleidoscopeMandala, "kaleidoscope-mandala", "kaleidoscope-mandala")       \
  X(KaleidoscopeSmooth, "kaleidoscope-smooth", "kaleidoscope-smooth")          \
  X(KaleidoscopeStainedGlass, "kaleidoscope-stained-glass",                    \
    "kaleidoscope-stained-glass")                                              \
  X(KaleidoscopeMelt, "kaleidoscope-melt", "kaleidoscope-melt")                \
  X(LiquidWindows, "liquid-windows", "liquid-windows")                         \
  X(MermaidSkin, "mermaid-skin", "mermaid-skin")                               \
  X(MeshFeedback, "mesh-feedback", "MeshFeedback")                             \
  X(Millefiori, "millefiori", "millefiori")                                    \
  X(MindSplatter, "mind-splatter", "MindSplatter")                             \
  X(MobiusGrid, "mobius-grid", "mobius-grid")                                  \
  X(MobiusRings, "mobius-rings", "MobiusRings")                                \
  X(PetalFlow, "petal-flow", "PetalFlow")                                      \
  X(Raymarch, "raymarch", "Raymarch")                                          \
  X(RingShower, "ring-shower", "RingShower")                                   \
  X(RingSpin, "ring-spin", "RingSpin")                                         \
  X(ShapeShifter, "shape-shifter", "ShapeShifter")                             \
  X(SphericalHarmonics, "spherical-harmonics", "SphericalHarmonics")           \
  X(Thrusters, "thrusters", "Thrusters")                                       \
  X(ShaderChain, "shader-chain", "shader-chain")

/**
 * @brief Declares one immutable identity and checks its machine ID.
 * @param Name Class/browser name.
 * @param Id Canonical machine ID.
 * @param Seed Immutable RNG compatibility key.
 */
#define HS_DECLARE_EFFECT_IDENTITY(Name, Id, Seed)                             \
  struct Name {                                                                \
    static constexpr EffectIdentity VALUE{Id, #Name, Seed};                    \
  };                                                                           \
  static_assert(canonical_effect_id(Id),                                       \
                "effect machine IDs must be kebab-case");
HS_EFFECT_IDENTITY_LIST(HS_DECLARE_EFFECT_IDENTITY)
#undef HS_DECLARE_EFFECT_IDENTITY
} // namespace effect_identities
/**
 * @brief Resets an effect's RNG without changing its historical seed key.
 * @param identity Canonical effect and compatibility identity.
 */
inline void seed_effect_stream(const EffectIdentity &identity) {
  seed_effect_stream(identity.rng_seed_key);
}

/**
 * @brief Resets the RNG stream for a registered effect type.
 * @tparam EffectType Concrete effect with a declared effect identity.
 */
template <typename EffectType> inline void seed_effect_stream() {
  seed_effect_stream(effect_seed_key<EffectType>(EffectType::EFFECT_ID));
}

} // namespace hs
