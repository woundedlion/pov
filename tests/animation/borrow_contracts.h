/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Borrow-contract guards (compile-time)
// ----------------------------------------------------------------------------
// Animations that borrow effect-owned state accept an lvalue and reject a
// temporary, which would dangle.
// ============================================================================
namespace borrow_guard {
/** @brief Orientation alias used by the borrow-contract static_asserts. */
using Ori = math::Orientation<16>;

static_assert(std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                      ProceduralPath &, int>,
              "Motion must accept an lvalue (effect-owned) path");
static_assert(!std::is_constructible_v<Animation::Motion<288, 16>, Ori &,
                                       ProceduralPath &&, int>,
              "Motion must REJECT a temporary path (would dangle)");

static_assert(
    std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &,
                            const Lerpable &, int, EasingFn>,
    "Lerp must accept lvalue (effect-owned) start/target");
static_assert(
    !std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &&,
                             const Lerpable &, int, EasingFn>,
    "Lerp must REJECT a temporary start (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::Lerp, Lerpable &, const Lerpable &,
                             const Lerpable &&, int, EasingFn>,
    "Lerp must REJECT a temporary target (would dangle)");

static_assert(
    std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                            const GenerativePalette::Snapshot &,
                            const GenerativePalette::Snapshot &, int, EasingFn>,
    "ColorWipe must accept effect-owned snapshots");
static_assert(
    !std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                             GenerativePalette::Snapshot &&,
                             const GenerativePalette::Snapshot &, int,
                             EasingFn>,
    "ColorWipe must REJECT a temporary start snapshot (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::ColorWipe, GenerativePalette &,
                             const GenerativePalette::Snapshot &,
                             GenerativePalette::Snapshot &&, int, EasingFn>,
    "ColorWipe must REJECT a temporary target snapshot (would dangle)");

static_assert(
    std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                            const float &, const float &, int>,
    "MobiusFlow must accept lvalue (effect-owned) scalars");
static_assert(
    !std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                             const float &&, const float &, int>,
    "MobiusFlow must REJECT a temporary num_rings (would dangle)");
static_assert(
    !std::is_constructible_v<Animation::MobiusFlow, math::MobiusParams &,
                             const float &, const float &&, int>,
    "MobiusFlow must REJECT a temporary num_lines (would dangle)");
} // namespace borrow_guard
