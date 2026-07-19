/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>

#include "engine/platform.h"
#include "engine/memory.h"
#include "engine/util.h"
#include "math/3dmath.h"
#include "color/color.h"

/**
 * @brief Represents a "Fragment" or a potential pixel/vertex with associated
 * data registers. Mirrors the JS Fragment structure for shader compatibility.
 */
struct Fragment {
  Vector pos;        /**< Position (typically a unit vector on the sphere). */
  float v0 = 0.0f;   /**< Register 0 (usually normalized progress t) */
  float v1 = 0.0f;   /**< Register 1 (usually arc length/distance) */
  float v2 = 0.0f;   /**< Register 2 (stroke coverage or index/id) */
  float v3 = 0.0f;   /**< Register 3 (auxiliary) */
  float size = 1.0f; /**< Size metric (e.g. radius/apothem) for normalization */
  float age = 0.0f;  /**< Age of the operation/trail */
  Color4 color = Color4(0, 0, 0, 0); /**< Output Color (RGBA) */

  /**
   * @brief Linear interpolation between two fragments.
   * @param a Start fragment.
   * @param b End fragment.
   * @param t Interpolation factor (0.0 to 1.0).
   * @return The interpolated fragment; pos, v0-v3, age, size, and color are all
   * interpolated, so no register (notably size, the fragment_edge_dist
   * denominator) resets to its struct default.
   */
  static Fragment lerp(const Fragment &a, const Fragment &b, float t) {
    Fragment f;
    f.pos =
        a.pos + (b.pos - a.pos) * t;
    f.v0 = a.v0 + (b.v0 - a.v0) * t;
    f.v1 = a.v1 + (b.v1 - a.v1) * t;
    f.v2 = a.v2 + (b.v2 - a.v2) * t;
    f.v3 = a.v3 + (b.v3 - a.v3) * t;
    f.age = a.age + (b.age - a.age) * t;
    f.size = a.size + (b.size - a.size) * t;
    f.color = a.color.lerp(b.color, t);
    return f;
  }
};

/**
 * @brief Normalized inward depth from the nearest face edge for a rasterized
 * fragment, in face-relative units.
 * @param f Rasterized fragment; v1 holds the signed edge distance (negative
 * inside the face) and size holds the face's reference size.
 * @return `-v1 / size` (inward depth in face-relative units), or 0 for
 * degenerate (near-zero-size) faces.
 * @details Shared by the topology shaders (HankinSolids/IslamicStars) which
 * both gradient-map this depth.
 */
inline float fragment_edge_dist(const Fragment &f) {
  return (f.size > math::TOLERANCE) ? (-f.v1 / f.size) : 0.0f;
}

/**
 * @brief Resolves a fragment's palette slot from its face's topology class.
 * @tparam NumPalettes Palette count; the class index wraps modulo this.
 * @param f Rasterized fragment; v2 carries the integer face index.
 * @param topology Per-face topology-class indices.
 * @param num_faces Length of `topology`; an out-of-range face index falls back
 * to class 0 rather than reading out of bounds.
 * @return The palette slot in [0, NumPalettes).
 */
template <size_t NumPalettes>
inline int mesh_topology_slot(const Fragment &f, const int *topology,
                              int num_faces) {
  int faceIdx = static_cast<int>(f.v2);
  int topoIdx = (faceIdx >= 0 && faceIdx < num_faces) ? topology[faceIdx] : 0;
  return wrap(topoIdx, static_cast<int>(NumPalettes));
}

/**
 * @brief Shared face-topology fragment shading for the mesh effects.
 * @tparam PaletteBank Indexable bank of palettes exposing `bank[i].get(t)`.
 * @tparam NumPalettes Palette count (deduced from `palette_idx`).
 * @param f Rasterized fragment; v2 carries the integer face index.
 * @param topology Per-face topology-class indices.
 * @param num_faces Length of `topology`; an out-of-range face index falls back
 * to class 0 rather than reading out of bounds.
 * @param palette_bank Bank of per-class palettes.
 * @param palette_idx Maps a topology class to a palette slot in the bank.
 * @param gain Multiplier on the edge-distance gradient before clamping to [0,1].
 * @param opacity Output alpha.
 * @return The face's palette color shaded by edge distance, at `opacity`.
 * @details Single home for the face-color/edge-shade policy shared by
 * IslamicStars (gain 1.0) and HankinSolids (gain = intensity), so the two
 * cannot drift. The class index wraps modulo NumPalettes.
 */
template <typename PaletteBank, size_t NumPalettes>
inline Color4 shade_mesh_topology(const Fragment &f, const int *topology,
                                  int num_faces, PaletteBank &palette_bank,
                                  const std::array<int, NumPalettes> &palette_idx,
                                  float gain, float opacity) {
  float t = hs::clamp(fragment_edge_dist(f) * gain, 0.0f, 1.0f);
  int slot = mesh_topology_slot<NumPalettes>(f, topology, num_faces);
  Color4 c = palette_bank[palette_idx[slot]].get(t);
  c.alpha = opacity;
  return c;
}

/**
 * @brief Segue-aware variant of shade_mesh_topology: routes the edge distance
 * and resulting color through a segue policy's shading hooks.
 * @tparam SegueT Segue policy type (see namespace Segue in
 * animation/mesh.h); duck-typed, so this header needs no dependency on it.
 * @param f Rasterized fragment; v2 carries the integer face index.
 * @param topology Per-face topology-class indices.
 * @param num_faces Length of `topology`.
 * @param palette_bank Bank of per-class palettes.
 * @param palette_idx Maps a topology class to a palette slot in the bank.
 * @param gain Multiplier on the edge-distance gradient before clamping to [0,1].
 * @param segue Policy whose fill/grade/opacity hooks shape the fragment.
 * @param phase Segue phase for this fragment (face-local for per-face segues).
 * @return The shaded color; fully transparent when the segue's fill culls the
 * fragment.
 */
template <typename PaletteBank, size_t NumPalettes, typename SegueT>
inline Color4 shade_mesh_topology(const Fragment &f, const int *topology,
                                  int num_faces, PaletteBank &palette_bank,
                                  const std::array<int, NumPalettes> &palette_idx,
                                  float gain, const SegueT &segue, float phase) {
  float t = hs::clamp(fragment_edge_dist(f) * gain, 0.0f, 1.0f);
  float cover = segue.fill(t, phase);
  if (cover <= 0.0f)
    return Color4();
  int slot = mesh_topology_slot<NumPalettes>(f, topology, num_faces);
  Color4 c = segue.grade(palette_bank[palette_idx[slot]].get(t), phase);
  c.alpha = cover * segue.opacity(phase);
  return c;
}

/**
 * @brief Face-hoisted segue shade: the caller resolves the fragment's palette
 * once per face and passes it directly, so the per-fragment path skips the
 * topology-slot lookup and palette-bank indirection.
 * @tparam Palette Baked palette type exposing `Color4 get(float) const`.
 * @tparam SegueT Segue policy (see shade_mesh_topology's segue overload).
 * @param f Rasterized fragment.
 * @param palette The face's already-resolved palette.
 * @param gain Multiplier on the edge-distance gradient before clamping to [0,1].
 * @param segue Policy whose fill/grade/opacity hooks shape the fragment.
 * @param phase Segue phase for this fragment (face-local for per-face segues).
 * @return The shaded color; fully transparent when the segue's fill culls it.
 */
template <typename Palette, typename SegueT>
inline Color4 shade_mesh_topology(const Fragment &f, const Palette &palette,
                                  float gain, const SegueT &segue, float phase) {
  float t = hs::clamp(fragment_edge_dist(f) * gain, 0.0f, 1.0f);
  float cover = segue.fill(t, phase);
  if (cover <= 0.0f)
    return Color4();
  Color4 c = segue.grade(palette.get(t), phase);
  c.alpha = cover * segue.opacity(phase);
  return c;
}

/**
 * @brief Unit half-vector for the metallic Blinn-Phong specular lobe, with the
 *        light tilted off-axis along the surface tangent.
 * @param light_dir Direction toward the light in world space (unit length).
 * @param view_dir Direction toward the viewer in world space (unit length).
 * @param tangent Surface tangent used to tilt the specular highlight off-axis.
 * @return The unit half-vector, or the un-normalized sum when either
 *         intermediate degenerates to near-zero length.
 * @details Depends on no per-pixel input, so a shader over a fixed light/view/
 *          tangent frame computes it once and passes it to every
 *          shade_blinn_phong call.
 */
HS_O3_BEGIN
inline Vector blinn_phong_half(const Vector &light_dir, const Vector &view_dir,
                               const Vector &tangent) {
  Vector light = light_dir + tangent * 0.3f;
  float ll = light.length();
  if (ll > math::TOLERANCE)
    light /= ll;
  Vector half = light + view_dir;
  float hl = half.length();
  if (hl > math::TOLERANCE)
    half /= hl;
  return half;
}

/**
 * @brief Metallic headlight shade factor: half-Lambert diffuse, tight
 *        specular, Fresnel rim.
 * @param normal_w Surface normal in world space (unit length).
 * @param view_dir Direction toward the viewer in world space (unit length);
 *        the light is coincident with the viewer, so this is the light
 *        direction too and one dot product feeds both diffuse and Fresnel.
 * @param half_w Unit half-vector from blinn_phong_half().
 * @param diffuse_w Weight of the half-Lambert diffuse term.
 * @param specular_w Weight of the specular term.
 * @param fresnel_w Weight of the Fresnel rim term.
 * @return Scalar shade factor (ambient base plus weighted diffuse/specular/
 *         Fresnel terms); the caller multiplies it by the surface color for
 *         the metallic look.
 */
inline float shade_blinn_phong(const Vector &normal_w, const Vector &view_dir,
                               const Vector &half_w, float diffuse_w,
                               float specular_w, float fresnel_w) {
  float ndotv = dot(normal_w, view_dir);
  float half_lam = ndotv * 0.5f + 0.5f;
  float diffuse = half_lam * half_lam;

  float ndoth = std::max(0.0f, dot(normal_w, half_w));
  // ndoth^32 via repeated squaring
  float spec = ndoth * ndoth; // ^2
  spec *= spec;               // ^4
  spec *= spec;               // ^8
  spec *= spec;               // ^16
  spec *= spec;               // ^32

  float fresnel = 1.0f - hs::clamp(ndotv, 0.0f, 1.0f);
  fresnel = fresnel * fresnel * fresnel;

  return 0.05f + diffuse * diffuse_w + spec * specular_w + fresnel * fresnel_w;
}
HS_O3_END

/**
 * @brief A list of fragments, equivalent to 'Points' in the JS context but with
 * full register support.
 */
using Fragments = ArenaVector<Fragment>;
