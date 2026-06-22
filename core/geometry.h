/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <algorithm>
#include <array>
#include "platform.h"
#include "constants.h"
#include "3dmath.h"
#include "color.h"
#include "memory.h" // for ArenaVector (Fragments); geometry must not depend on mesh.h
#include "static_circular_buffer.h" // for StaticCircularBuffer (Dots/Points)
#include "util.h"                    // for wrap()

/**
 * @brief Represents a "Fragment" or a potential pixel/vertex with associated
 * data registers. Mirrors the JS Fragment structure for shader compatibility.
 */
struct Fragment {
  Vector pos;        /**< Position (typically a unit vector on the sphere). */
  float v0 = 0.0f;   /**< Register 0 (usually normalized progress t) */
  float v1 = 0.0f;   /**< Register 1 (usually arc length/distance) */
  float v2 = 0.0f;   /**< Register 2 (usually index/id) */
  float v3 = 0.0f;   /**< Register 3 (auxiliary) */
  float size = 1.0f; /**< Size metric (e.g. radius/apothem) for normalization */
  float age = 0.0f;  /**< Age of the operation/trail */
  Color4 color = Color4(0, 0, 0, 0); /**< Output Color (RGBA) */

  /**
   * @brief Linear interpolation between two fragments.
   * @param a Start fragment.
   * @param b End fragment.
   * @param t Interpolation factor (0.0 to 1.0).
   * @return The interpolated fragment. Interpolates pos, v0-v3, age, size, and
   * color, so a sample between two control points carries every register rather
   * than resetting size (the fragment_edge_dist denominator) and color to their
   * struct defaults.
   */
  static Fragment lerp(const Fragment &a, const Fragment &b, float t) {
    Fragment f;
    f.pos =
        a.pos + (b.pos - a.pos) * t; // Linear, not slerp: any slerp of pos is
                                     // applied by the caller.
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
  return (f.size > 0.0001f) ? (-f.v1 / f.size) : 0.0f;
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
  int faceIdx = static_cast<int>(f.v2);
  int topoIdx = (faceIdx >= 0 && faceIdx < num_faces) ? topology[faceIdx] : 0;
  float t = hs::clamp(fragment_edge_dist(f) * gain, 0.0f, 1.0f);
  // topology[] can hold a raw (possibly negative) class id, and C++ % keeps the
  // sign of the dividend — a bare `topoIdx % NumPalettes` could index palette_idx
  // out of bounds. Floor it into [0, NumPalettes); branchless-cheap on this
  // per-pixel path, so no HS_CHECK trap here.
  int slot = topoIdx % static_cast<int>(NumPalettes);
  if (slot < 0)
    slot += static_cast<int>(NumPalettes);
  Color4 c = palette_bank[palette_idx[slot]].get(t);
  c.alpha = opacity;
  return c;
}

/**
 * @brief A list of fragments, equivalent to 'Points' in the JS context but with
 * full register support.
 */
using Fragments = ArenaVector<Fragment>;

/**
 * @brief No-op vertex shader; leaves every fragment unchanged.
 */
struct NullVertexShader {
  /**
   * @brief Does nothing to the fragment (ignored).
   */
  void operator()(Fragment &) const {}
};

/**
 * @brief No-op fragment shader; emits a fully transparent color.
 * @details Used where a shader slot is required but no shading is wanted.
 */
struct NullFragmentShader {
  /**
   * @brief Returns a transparent color regardless of input.
   * @details The sample position and source fragment are ignored.
   * @return Transparent black Color4(0, 0, 0, 0).
   */
  Color4 operator()(const Vector &, const Fragment &) const {
    return Color4(0, 0, 0, 0);
  }
};

/**
 * @brief Unit vector along the Cartesian X-axis.
 */
static constexpr Vector X_AXIS(1, 0, 0);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
static constexpr Vector Y_AXIS(0, 1, 0);
/**
 * @brief Unit vector along the Cartesian Z-axis.
 */
static constexpr Vector Z_AXIS(0, 0, 1);
/**
 * @brief Unit vector along the Cartesian Y-axis.
 */
static constexpr Vector UP = Y_AXIS;

/**
 * @brief Structure representing 2D floating-point pixel coordinates.
 */
struct PixelCoords {
  float x; /**< Horizontal coordinate. */
  float y; /**< Vertical coordinate. */
};

/**
 * @brief Structure representing a single rendered point (dot) in the scene.
 * @details Stores the 3D position and the color to be plotted.
 */
struct Dot {
  /**
   * @brief Constructs a Dot.
   * @param v The 3D position vector.
   * @param color The pixel color with alpha.
   */
  Dot(const Vector &v, const Color4 &color) : position(v), color(color) {}

  /**
   * @brief Copy constructor — defaulted to keep Dot trivially copyable so the
   * 1024-deep StaticCircularBuffer can memcpy/vectorize copies. A hand-written
   * member-wise body is identical in effect but defeats that (see Vector's
   * defaulted-copy rationale in 3dmath.h).
   */
  Dot(const Dot &d) = default;

  Vector position; /**< The 3D position (unit vector). */
  Color4 color;    /**< The color of the dot. */
};

/**
 * @brief Type alias for a circular buffer used to store active dots/fragments.
 * @details Capacity is set to 1024.
 */
using Dots = StaticCircularBuffer<Dot, 1024>;

/**
 * @brief Type alias for a circular buffer used to store geometry points
 * (Vectors).
 * @details Capacity is set to 1024.
 */
using Points = StaticCircularBuffer<Vector, 1024>;

/**
 * @brief Struct to hold Log-Polar coordinates.
 */
struct LogPolar {
  float rho;   /**< Log-radius (natural log of the complex-plane radius). */
  float theta; /**< Angle in radians. */
};

/**
 * @brief Converts a pixel y-coordinate to a spherical phi angle.
 * @param y The pixel y-coordinate [0, h_virt - 1].
 * @param h_virt The virtual height.
 * @return The spherical phi angle in radians.
 */
inline float y_to_phi(float y, int h_virt) {
  // Mirrors the templated y_to_phi<H>'s static_assert(H_VIRT > 1): the mapping
  // degenerates (divide-by-zero at h_virt == 1, garbage for h_virt <= 0). Cold
  // setup/conversion helper, not the per-pixel loop; every real caller passes a
  // compile-time H_VIRT, so the check constant-folds away under optimization.
  HS_CHECK(h_virt > 1, "y_to_phi: h_virt must be > 1");
  return (y * PI_F) / (h_virt - 1);
}

/**
 * @brief Converts a spherical phi angle to a pixel y-coordinate.
 * @param phi The spherical phi angle in radians.
 * @param h_virt The virtual height.
 * @return The pixel y-coordinate [0, h_virt - 1].
 */
inline float phi_to_y(float phi, int h_virt) {
  // See y_to_phi(float,int): same h_virt > 1 precondition as the templated
  // overload's static_assert. Constant-folds away for compile-time H_VIRT.
  HS_CHECK(h_virt > 1, "phi_to_y: h_virt must be > 1");
  return (phi * (h_virt - 1)) / PI_F;
}

/**
 * @brief phi -> pixel-y for a compile-time logical height H.
 * @tparam H Logical (not virtual) height; H_VIRT is derived as H + hs::H_OFFSET.
 * @param phi The spherical phi angle in radians.
 * @return The pixel y-coordinate.
 * @details Derives H_VIRT from H plus hs::H_OFFSET so callers pass the logical
 * height, not the virtual one.
 */
template <int H> inline float phi_to_y(float phi) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  static_assert(H_VIRT > 1, "phi<->y mapping degenerates when H_VIRT <= 1");
  return (phi * (H_VIRT - 1)) / PI_F;
}

/**
 * @brief Precomputed lookup table for scanline phi angles.
 * @tparam H Logical height; the table has H_VIRT = H + hs::H_OFFSET entries.
 */
template <int H> struct PhiLUT {
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static std::array<float, H_VIRT> data; /**< phi per virtual row, radians. */
  // Lazy-init guard. THREAD-SAFETY CONTRACT: the `if (!initialized) init()`
  // pattern at the call sites is a non-atomic check-then-set and is safe ONLY
  // because rendering is single-threaded — every consumer (the scanline
  // rasterizers, y_to_phi) runs on the render thread, and engine setup calls
  // init_geometry_luts() eagerly before the first frame, so on hardware the
  // column-sweep ISR never observes a half-filled table. It is NOT a
  // concurrency safeguard; a second writer would race. See init_geometry_luts.
  static bool initialized; /**< Lazy-init guard; true once data is filled. */
  /**
   * @brief Fills the phi table for every virtual row and marks it initialized.
   */
  static void init() {
    for (int y = 0; y < H_VIRT; y++) {
      data[y] = y_to_phi(static_cast<float>(y), H_VIRT);
    }
    initialized = true;
  }
};

template <int H> std::array<float, PhiLUT<H>::H_VIRT> PhiLUT<H>::data;
template <int H> bool PhiLUT<H>::initialized = false;

/**
 * @brief LUT-backed pixel-y -> phi for integer rows at compile-time height H.
 * @tparam H Logical height selecting the PhiLUT<H> table.
 * @param y Integer pixel row in [0, H_VIRT).
 * @return The spherical phi angle in radians for that row.
 * @details Lazily fills PhiLUT on first touch; traps an out-of-range row via
 * HS_CHECK.
 */
template <int H> inline float y_to_phi(int y) {
  if (!PhiLUT<H>::initialized) {
    PhiLUT<H>::init();
  }
  HS_CHECK(y >= 0 && y < PhiLUT<H>::H_VIRT);
  return PhiLUT<H>::data[y];
}

/**
 * @brief Pixel-y -> phi for fractional rows at compile-time height H.
 * @tparam H Logical height; H_VIRT is H + hs::H_OFFSET.
 * @param y Fractional pixel row.
 * @return The spherical phi angle in radians.
 * @details Snaps to the LUT for near-integer y; otherwise computes the angle
 * analytically.
 */
template <int H> inline float y_to_phi(float y) {
  if (std::abs(y - std::floor(y)) < TOLERANCE) {
    int iy = static_cast<int>(y);
    if (iy >= 0 && iy < PhiLUT<H>::H_VIRT) {
      return y_to_phi<H>(iy);
    }
  }
  constexpr int H_VIRT = H + hs::H_OFFSET;
  static_assert(H_VIRT > 1, "phi<->y mapping degenerates when H_VIRT <= 1");
  return (y * PI_F) / (H_VIRT - 1);
}

/**
 * @brief Split trig lookup tables for efficient vector reconstruction.
 * @tparam W Width (column count).
 * @tparam H Logical height; phi tables have H_VIRT = H + hs::H_OFFSET entries.
 * @details Caches sin/cos for theta (per column) and phi (per row) separately.
 * Reconstructs vectors with 3 multiplies instead of storing full Vectors.
 * Memory: ~(4*W + 4*H_VIRT) floats vs W*H_VIRT Vectors — a ~145x reduction.
 */
template <int W, int H> struct TrigLUT {
  static constexpr int H_VIRT = H + hs::H_OFFSET;
  static std::array<float, W> sin_theta;     /**< sin(theta) per column. */
  static std::array<float, W> cos_theta;     /**< cos(theta) per column. */
  static std::array<float, H_VIRT> sin_phi;  /**< sin(phi) per virtual row. */
  static std::array<float, H_VIRT> cos_phi;  /**< cos(phi) per virtual row. */
  // Lazy-init guard. Same non-atomic check-then-set thread-safety contract as
  // PhiLUT::initialized above: single-render-thread only, with eager
  // init_geometry_luts() at engine setup as the production first-touch.
  static bool initialized; /**< Lazy-init guard; true once tables are filled. */
  /**
   * @brief Fills the theta and phi sin/cos tables and marks them initialized.
   * @details Ensures PhiLUT<H> is populated first to source the phi angles.
   */
  static void init() {
    if (!PhiLUT<H>::initialized) {
      PhiLUT<H>::init();
    }
    for (int x = 0; x < W; x++) {
      float theta = (x * 2 * PI_F) / W;
      sin_theta[x] = sinf(theta);
      cos_theta[x] = cosf(theta);
    }
    for (int y = 0; y < H_VIRT; y++) {
      float phi = PhiLUT<H>::data[y];
      sin_phi[y] = sinf(phi);
      cos_phi[y] = cosf(phi);
    }
    initialized = true;
  }
};

template <int W, int H> std::array<float, W> TrigLUT<W, H>::sin_theta;
template <int W, int H> std::array<float, W> TrigLUT<W, H>::cos_theta;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::sin_phi;
template <int W, int H>
std::array<float, TrigLUT<W, H>::H_VIRT> TrigLUT<W, H>::cos_phi;
template <int W, int H> bool TrigLUT<W, H>::initialized = false;

/**
 * @brief Eagerly fill the scanline LUTs for resolution <W, H>.
 * @tparam W Width (column count).
 * @tparam H Logical height.
 * @details Engine setup calls this once, before the first frame, so the tables
 * are fully populated before any rendering — and, on hardware, before the
 * column-sweep ISR could ever observe a partially-filled table. With eager
 * init in place the per-call `if (!initialized) init()` guards scattered
 * through the scanline rasterizers (Scan/Plot/SDF/Filter, `pixel_to_vector`)
 * are never the *first* touch in production; they remain only as a
 * lazy fallback for unit tests and offline tools that render at other
 * resolutions without going through engine setup. Those guards are a
 * non-atomic check-then-set, NOT a concurrency safeguard: their correctness
 * rests on this eager call and on the single-render-thread assumption (see the
 * THREAD-SAFETY CONTRACT on PhiLUT/TrigLUT::initialized). Idempotent.
 */
template <int W, int H> inline void init_geometry_luts() {
  PhiLUT<H>::init();
  TrigLUT<W, H>::init();
}

/**
 * @brief Recovers an effect's compile-time <W, H> from its type so a driver's
 * `show<E>()` can eager-init the LUTs without the caller restating the
 * resolution.
 * @tparam E The effect type, of the form `Eff<W, H>`.
 * @details Every effect is `template <int W, int H> class E`, so the partial
 * specialization matches them all.
 */
template <typename E> struct GeometryResolution;
/**
 * @brief Partial specialization that destructures an effect's <W, H>.
 * @tparam Eff The effect class template.
 * @tparam W Width recovered from the effect type.
 * @tparam H Height recovered from the effect type.
 */
template <template <int, int> class Eff, int W, int H>
struct GeometryResolution<Eff<W, H>> {
  /**
   * @brief Eager-inits the geometry LUTs for the recovered <W, H>.
   */
  static void init() { init_geometry_luts<W, H>(); }
};

/**
 * @brief Reconstruct a vector from pixel coordinates using split trig LUTs.
 * @tparam W Width.
 * @tparam H Height.
 * @param x X coordinate (column).
 * @param y Y coordinate (row).
 * @return Unit vector on the sphere.
 */
template <int W, int H> Vector pixel_to_vector(int x, int y) {
  if (!TrigLUT<W, H>::initialized) {
    TrigLUT<W, H>::init();
  }
  // Per-pixel hot path: use a stripped assert (compiled out under NDEBUG on
  // device), not an always-on HS_CHECK. The LUT-domain invariant is trapped
  // once per draw on the cold path by Scan::check_lut_domain; this assert only
  // backs the native/WASM-debug builds. Local avoids the comma in TrigLUT<W, H>.
  constexpr int kHVirt = TrigLUT<W, H>::H_VIRT;
  assert(x >= 0 && x < W && y >= 0 && y < kHVirt);
  float sp = TrigLUT<W, H>::sin_phi[y];
  return Vector(sp * TrigLUT<W, H>::cos_theta[x], TrigLUT<W, H>::cos_phi[y],
                sp * TrigLUT<W, H>::sin_theta[x]);
}

/**
 * @brief Reconstruct a unit vector from fractional pixel coordinates.
 * @tparam W Width.
 * @tparam H Height.
 * @param x Fractional X coordinate (column).
 * @param y Fractional Y coordinate (row).
 * @return Unit vector on the sphere.
 * @details Snaps to the integer LUT path when both coordinates are
 * near-integer; otherwise builds the vector analytically from spherical angles.
 */
template <int W, int H> Vector pixel_to_vector(float x, float y) {
  if (std::abs(x - floor(x)) < TOLERANCE &&
      std::abs(y - floor(y)) < TOLERANCE) {
    return pixel_to_vector<W, H>(static_cast<int>(x), static_cast<int>(y));
  }
  // y_to_phi<H> already accounts for H_OFFSET internally; pass H, not H_VIRT.
  return Vector(Spherical((x * 2 * PI_F) / W, y_to_phi<H>(y)));
}

/**
 * @brief Converts a 3D unit vector back to 2D pixel coordinates.
 * @tparam W The width.
 * @tparam H The height.
 * @param v The input unit vector.
 * @return The 2D PixelCoords.
 */
template <int W, int H> PixelCoords vector_to_pixel(const Vector &v) {
  constexpr int H_VIRT = H + hs::H_OFFSET;
  // Skip normalize: fast_atan2 is scale-invariant, and inputs from
  // the rasterizer are already unit vectors.
  float theta = fast_atan2(v.z, v.x);
  float phi = fast_acos(hs::clamp(v.y, -1.0f, 1.0f));
  PixelCoords p({wrap((theta * W) / (2 * PI_F), W), phi_to_y(phi, H_VIRT)});
  return p;
}

/**
 * @brief Converts Log-Polar coordinates (rho, theta) to a vector on the unit
 * sphere. Maps: Log-Polar -> Complex Plane -> Inverse Stereographic -> Sphere
 * @param rho The log-radius (natural logarithm of the radius on the complex
 * plane).
 * @param theta The angle in radians.
 * @return Unit vector on the sphere (unit by construction).
 */
inline Vector logPolarToVector(float rho, float theta) {
  const float R = expf(rho);
  const float R2 = R * R;
  const float y = (R2 - 1.0f) / (R2 + 1.0f);
  const float r_xz = sqrtf(std::max(0.0f, 1.0f - y * y));
  // Already unit-length: r_xz^2 + y^2 = (1 - y^2) + y^2 = 1 (the cos/sin fold the
  // xz components back to r_xz^2), so no normalize() is needed.
  return Vector(r_xz * cosf(theta), y, r_xz * sinf(theta));
}

/**
 * @brief Converts a vector on the unit sphere to Log-Polar coordinates.
 * Maps: Sphere -> Stereographic -> Complex Plane -> Log-Polar
 * @param v Normalized vector on the unit sphere.
 * @return Log-Polar coordinates.
 */
inline LogPolar vectorToLogPolar(const Vector &v) {
  // Both poles are stereographic singularities: at the north pole (v.y -> +1)
  // the planar radius R -> +inf, at the south pole (v.y -> -1) it goes R -> 0,
  // so rho = 0.5*log((1+y)/(1-y)) tends to +inf and -inf respectively. Clamp
  // each to a symmetric finite sentinel so neither pole leaks a non-finite rho
  // into downstream arithmetic.
  const float numer = 1.0f + v.y;
  const float denom = 1.0f - v.y;
  if (std::abs(denom) < 0.00001f) {
    return {10.0f, 0.0f}; // North pole sentinel (rho -> +inf)
  }
  if (std::abs(numer) < 0.00001f) {
    return {-10.0f, 0.0f}; // South pole sentinel (rho -> -inf)
  }
  const float rho = 0.5f * logf(numer / denom);
  const float theta = fast_atan2(v.z, v.x);
  return {rho, theta};
}

/**
 * @brief Calculates a point on the Fibonacci spiral on the unit sphere.
 * @param n The total number of points in the spiral.
 * @param eps The epsilon offset for the spiral.
 * @param i The index of the point to calculate.
 * @return The point on the unit sphere.
 */
inline Vector fib_spiral(int n, float eps, int i) {
  // Clamp before acosf: at i == n-1 (or eps near 1) float rounding can push the
  // argument just past -1, where acosf returns NaN; every sibling acos call in
  // the codebase clamps the same way.
  float phi = acosf(hs::clamp(1.0f - (2.0f * (static_cast<float>(i) + eps)) /
                                         static_cast<float>(n),
                              -1.0f, 1.0f));
  float theta = fmodf((2.0f * PI_F * static_cast<float>(i) * G), (2.0f * PI_F));
  // Y-up convention
  return Vector(sinf(phi) * cosf(theta), cosf(phi), sinf(phi) * sinf(theta))
      .normalized();
}

/**
 * @brief Class managing the current rotation state of an object, maintaining
 * history for interpolation.
 * @tparam CAP Maximum number of orientation frames retained in history.
 * @details Stores a list of Quaternions (`orientations`) generated during the
 * current frame step.
 */
template <int CAP = 4> class Orientation {
public:
  static constexpr int CAPACITY = CAP;
  /**
   * @brief Default constructor (identity rotation).
   */
  Orientation() : num_frames(0) { set(Quaternion()); }

  /**
   * @brief Constructs with a specific initial quaternion.
   * @param q The initial quaternion.
   */
  Orientation(const Quaternion &q) : num_frames(0) { set(q); }

  /**
   * @brief Gets the number of recorded orientation frames in the current step.
   * @return The number of frames.
   */
  int length() const { return num_frames; }

  /**
   * @brief Rotates a vector by the current (latest) orientation.
   * @param v The vector to orient.
   * @return The rotated vector.
   */
  Vector orient(const Vector &v) const {
    HS_CHECK(num_frames >= 1);
    return rotate(v, orientations[num_frames - 1]);
  }

  /**
   * @brief Rotates a vector by an orientation at a specific historical frame
   * index.
   * @param v The vector to orient.
   * @param i The frame index (0 being oldest, length-1 being current).
   * @return The rotated vector.
   */
  Vector orient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return rotate(v, orientations[i]);
  }

  /**
   * @brief Rotates a vector backward by the inverse of the current (latest)
   * orientation.
   * @param v The vector to unorient.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v) const {
    HS_CHECK(num_frames >= 1);
    return rotate(v, orientations[num_frames - 1].conjugate());
  }

  /**
   * @brief Rotates a vector backward by the inverse of a specific historical
   * orientation.
   * @param v The vector to unorient.
   * @param i The frame index.
   * @return The unrotated vector.
   */
  Vector unorient(const Vector &v, int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return rotate(v, orientations[i].conjugate());
  }

  /**
   * @brief Gets the current (latest) quaternion.
   * @return The Quaternion reference.
   */
  const Quaternion &get() const {
    HS_CHECK(num_frames >= 1);
    return orientations[num_frames - 1];
  }

  /**
   * @brief Gets the quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  const Quaternion &get(int i) const {
    HS_CHECK(i >= 0 && i < num_frames);
    return orientations[i];
  }

  /**
   * @brief Sets the orientation, clearing all history.
   * @param q The new orientation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation &set(const Quaternion &q) {
    orientations[0] = q;
    num_frames = 1;
    return *this;
  }

  /**
   * @brief Pushes a new quaternion onto the history, tracking a motion step.
   * @param q The new rotation quaternion.
   * @return Reference to the Orientation object.
   */
  Orientation &push(const Quaternion &q) {
    // Overflow means a caller emitted more motion substeps than the history
    // holds — an invariant violation, not a recoverable transient. Fail fast
    // rather than silently dropping the frame (a soft log is invisible
    // on-device and corrupts the motion-blur trail).
    HS_CHECK(num_frames < CAPACITY);
    orientations[num_frames++] = q;
    return *this;
  }

  /**
   * @brief Collapses the orientation history, retaining only the latest
   * quaternion.
   * @details Used after rendering a motion step to reset the motion blur
   * history.
   * @return Reference to the Orientation object.
   */
  Orientation &collapse() {
    if (num_frames > 1) {
      orientations[0] = orientations[num_frames - 1];
      num_frames = 1;
    }
    return *this;
  }

  /**
   * @brief Access a mutable quaternion at a specific historical frame index.
   * @param i The frame index.
   * @return The Quaternion reference.
   */
  Quaternion &at(int i) {
    HS_CHECK(i >= 0 && i < num_frames);
    return orientations[i];
  }

  /**
   * @brief Increases the resolution of the history to 'count' steps, preserving
   * shape via Slerp.
   * @param count The target number of steps in the history.
   * @note When `count > CAPACITY` the trail is upsampled to `CAPACITY` instead.
   * This is intentional graceful degradation, not a trap: a `count` past
   * capacity means the per-frame angular speed exceeds what the trail can
   * sub-sample within the `MAX_ANGLE` (one-column) smoothness threshold — a
   * transient artistic input (e.g. a cranked speed slider), not a structural
   * invariant violation. The clamp keeps the write provably in-bounds, the
   * current orientation is always exact (`collapse()` retains the latest), and
   * the only consequence is that the motion-blur smear is sampled more coarsely
   * at extreme speed. Fail-fast (`HS_CHECK`) is reserved for cold seams and
   * memory/structural invariants; this is a per-frame hot path and the bounded
   * soft-degrade is the correct response here (cf. the DMA overrun-drop). To
   * trade RAM + per-frame slerp work for a smoother fast smear, raise `CAP` on
   * the affected `Orientation` rather than reintroducing a trap.
   */
  void upsample(int count) {
    // A history must hold at least one frame; count < 1 is a caller bug. This
    // also forecloses the 0/0 in t = i/(count-1) below: num_frames is always
    // >= 1 (every constructor routes through set()), so the early return makes
    // the interpolation loop reachable only when count >= 2 (i.e. count-1 >= 1).
    HS_CHECK(count >= 1);
    if (num_frames >= count)
      return;
    if (count > CAPACITY) // soft-degrade past capacity — see @note above
      count = CAPACITY;

    std::array<Quaternion, CAPACITY> old_orientations;
    std::copy(orientations.begin(), orientations.begin() + num_frames,
              old_orientations.begin());

    int old_num_frames = num_frames;

    for (int i = 0; i < count; ++i) {
      float t = static_cast<float>(i) / (count - 1);
      float source_float_index = t * (old_num_frames - 1);
      int idx = static_cast<int>(source_float_index);
      float frac = source_float_index - idx;

      orientations[i] = slerp(
          old_orientations[idx],
          old_orientations[std::min((int)old_num_frames - 1, idx + 1)], frac);
    }
    num_frames = count;
  }

private:
  std::array<Quaternion, CAPACITY>
      orientations; /**< Storage for historical quaternions. */
  int num_frames;   /**< The current number of active frames in history. */
};

/**
 * @brief Generates a truly random 3D unit vector (direction) using Marsaglia's
 * method.
 * @return A normalized random Vector.
 */
inline Vector random_vector() {
  float v1, v2, s;
  do {
    v1 = 2.0f * hs::rand_f() - 1.0f;
    v2 = 2.0f * hs::rand_f() - 1.0f;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0f || s == 0.0f);

  float sqrt_s = sqrtf(1.0f - s);
  return Vector(2.0f * v1 * sqrt_s, 2.0f * v2 * sqrt_s, 1.0f - 2.0f * s);
}

/**
 * @brief Parameters defining a Lissajous curve.
 */
struct LissajousParams {
  float m1; /**< Frequency coefficient for the axial components (X and Z). */
  float m2; /**< Frequency coefficient for the orbital component (Y). */
  float a;  /**< Phase shift in radians (matches the daydream lissajous tool). */
  float domain; /**< The total duration (t) over which the curve is drawn. */
};

/**
 * @brief Calculates a 3D point on the unit sphere corresponding to a spherical
 * Lissajous curve.
 * @param m1 Frequency coefficient for XZ plane.
 * @param m2 Frequency coefficient for Y axis.
 * @param a Phase shift in radians, used as-is. Matches the daydream lissajous
 *          designer (tools/lissajous.html), whose radians-labelled slider and
 *          exported snippet feed straight into this function.
 * @param t Time variable (or position along the domain).
 * @return The calculated 3D point (unit vector).
 */
inline Vector lissajous(float m1, float m2, float a, float t) {
  Vector v(sinf(m2 * t) * cosf(m1 * t - a), cosf(m2 * t),
           sinf(m2 * t) * sinf(m1 * t - a));
  return v.normalized();
}

/**
 * @brief An orthonormal basis { u, v, w }.
 */
struct Basis {
  Vector u, v, w; /**< Orthonormal axes; v is the normal, u and w span the plane. */
};

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 * @param orientation The orientation quaternion.
 * @param normal The normal vector; after rotation it becomes the 'v' axis.
 * @return The constructed Basis.
 */
inline Basis make_basis(const Quaternion &orientation, const Vector &normal) {
  // rotate() is only a pure rotation for a unit quaternion; assert it rather
  // than silently scale/shear the basis (fail-fast doctrine).
  HS_CHECK(std::abs(orientation.squared_magnitude() - 1.0f) < 0.01f);
  Vector v = rotate(normal, orientation).normalized();
  // Pick the reference axis least parallel to v, testing the *rotated* vectors
  // that are actually crossed below.
  Vector ref = rotate(X_AXIS, orientation).normalized();
  if (std::abs(dot(v, ref)) > math::COS_AXIS_PARALLEL) {
    ref = rotate(Y_AXIS, orientation).normalized();
  }
  Vector u = cross(v, ref).normalized();
  Vector w = cross(v, u).normalized();
  return {u, v, w};
}

/**
 * @brief Adjusted basis and radius for drawing on the opposite side of the
 * sphere.
 * @param basis The current basis {u, v, w}.
 * @param radius Angular radius (0-2).
 * @return A pair containing the adjusted Basis and radius.
 */
inline std::pair<Basis, float> get_antipode(const Basis &basis, float radius) {
  if (radius > 1.0f) {
    Basis new_basis;
    new_basis.u = -basis.u; // Flip U to maintain chirality
    new_basis.v = -basis.v; // Flip V (Antipode)
    new_basis.w = basis.w;  // W stays (Rotation axis)
    return {new_basis, 2.0f - radius};
  }
  return {basis, radius};
}
