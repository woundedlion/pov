/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <type_traits>
#include <utility>
#include "render/sdf.h"
#include "render/shading.h"
#include "mesh/mesh_class_types.h"
#include "color/color.h"
#include "render/filter.h"
#include "containers/static_circular_buffer.h"
#include "render/canvas.h"
#include "platform/platform.h"
#include "render/scan/raster.h"

/**
 * @file volume.h
 * @brief Scan::Volume: the raymarcher, and the TransformedVolume placement wrapper.
 */

namespace Scan {

// Volume::draw has exactly one caller (effects/Raymarch.h), and
// TransformedVolume is instantiated only there, so no other effect pays ITCM
// for this region.
HS_O3_BEGIN
/**
 * @brief Generic wrapper that places an SDF in world space via a center point
 *        and a rotation quaternion. Satisfies the Volume::draw shape concept.
 * @tparam SDF Underlying signed-distance shape type.
 * @details The quaternion q maps local→world: world_p = center +
 * rotate(local_p, q). ray_to_local uses q.inverse() to map world→local.
 */
template <typename SDF> struct TransformedVolume {
  const SDF &sdf;   /**< Underlying SDF evaluated in local space. */
  Vector center;    /**< World-space origin of the local frame. */
  Quaternion q_inv; /**< Precomputed inverse rotation (world→local). */

  /**
   * @brief Constructs the transform from a center and a local→world rotation.
   * @param sdf Underlying SDF, stored by reference.
   * @param center World-space origin of the local frame.
   * @param q Local→world rotation; its inverse is precomputed.
   */
  TransformedVolume(const SDF &sdf, const Vector &center, const Quaternion &q)
      : sdf(sdf), center(center), q_inv(q.inverse()) {}

  /**
   * @brief Transforms a ray origin and direction from world to local space.
   * @param ro Ray origin in world space.
   * @param vd Ray direction in world space.
   * @return Pair of {local origin, local direction}.
   */
  std::pair<Vector, Vector> ray_to_local(const Vector &ro,
                                         const Vector &vd) const {
    return {rotate(ro - center, q_inv), rotate(vd, q_inv)};
  }

  /**
   * @brief Transforms only a ray origin from world to local space.
   * @param ro Ray origin in world space.
   * @return Local-space origin.
   * @details The local direction is constant across the draw, so the volume loop
   * precomputes it once and calls this per pixel to transform only the origin.
   */
  Vector origin_to_local(const Vector &ro) const {
    return rotate(ro - center, q_inv);
  }

  /**
   * @brief Evaluates the underlying SDF at a local-space point.
   * @param local_p Query point in local space.
   * @return Signed distance to the surface in local units.
   */
  float distance(const Vector &local_p) const { return sdf.distance(local_p); }
};

/**
 * @brief Raymarch volume renderer with orthographic projection.
 * @details The render loop is internal: callers provide a Shape with a
 * `float distance(const Vector&) const` method (evaluated per march step) and a
 * fragment shader (evaluated once per hit to populate Fragment registers for
 * shading).
 *
 * Coordinate-space contract:
 *   - `view_dir` is the normalized direction all rays travel (camera → scene).
 *   - Ray origins are computed via orthographic projection: each pixel's
 *     position is projected onto the plane perpendicular to `view_dir`,
 *     then offset backward along `view_dir`.
 *   - `bounds_center` and `view_dir` must both be in physical LED space.
 *   - Filter::World::Orient rotates the *output* position passed to the
 *     canvas, not the ray.
 */
struct Volume {
  /** Sphere-trace step overrelaxation; 1 is the conservative march. */
  static constexpr float OVERRELAX_OMEGA = 1.3f;
  /** Occluder-probe march steps. */
  static constexpr int PROBE_STEPS = 24;
  /** Probe steps taken at the near step floor before the far one applies. */
  static constexpr int PROBE_NEAR_STEPS = 6;
  /** Near probe step floor, as a fraction of the bounding radius. */
  static constexpr float PROBE_FLOOR_NEAR = 0.04f;
  /** Far probe step floor, as a fraction of the bounding radius. */
  static constexpr float PROBE_FLOOR_FAR = 0.12f;

  /**
   * @brief Sphere-traces a ray in local space, recording the closest approach
   *        to the first surface reached.
   * @param shape Volume shape providing distance().
   * @param local_ro Ray origin in local space.
   * @param local_vd Unit ray direction in local space.
   * @param bounds_radius Bounding sphere radius (past-the-back early-out).
   * @param max_steps Maximum sphere-tracing steps.
   * @param aa_width Anti-aliasing band half-width (deep-hit early-out).
   * @param closest_local Output: local-space point of closest approach.
   * @return Signed distance at the closest approach (FLT_MAX if never sampled).
   * @details Inside the AA band, stops at the first rising local minimum (the
   * silhouette graze owning the pixel's coverage); marching past it would let a
   * deeper occluded surface steal the closest approach.
   *
   * Steps are overrelaxed by OVERRELAX_OMEGA (Keinert et al., "Enhanced Sphere
   * Tracing"): each sample's unbounding sphere must overlap its predecessor's,
   * and a step that breaks that overlap is rewound to the predecessor's surface
   * and the ray finished conservatively, so the trace still cannot cross a
   * surface undetected.
   */
  template <typename Shape>
  static __attribute__((always_inline)) float
  trace_closest(const Shape &shape, const Vector &local_ro,
                const Vector &local_vd, float bounds_radius, int max_steps,
                float aa_width, Vector &closest_local) {
    HS_PROFILE_DEEP(vol_trace);
    Vector local_p = local_ro;
    closest_local = local_ro;
    // Sentinel for "no surface seen yet": any real signed distance the
    // trace reports is smaller, so the first sample always wins.
    float closest_d = FLT_MAX;
    // Drops to 1 for the rest of the ray once a step fails the overlap test.
    float omega = OVERRELAX_OMEGA;
    float prev_r = 0.0f;
    float step_len = 0.0f;

    for (int i = 0; i < max_steps; ++i) {
      // Early out: ray has exited the back of the bounding sphere. The
      // local-space dot is compared against the world-space bounds_radius, valid
      // because ray_to_local is length-preserving (unit local_vd) and the caller
      // passes the shape center as bounds_center. Both are HS_CHECKed once per
      // draw at the top.
      if (local_p.x * local_vd.x + local_p.y * local_vd.y +
              local_p.z * local_vd.z >
          bounds_radius)
        break;

      float d = shape.distance(local_p);

      float r = d < 0.0f ? -d : d;
      if (omega > 1.0f && r + prev_r < step_len) {
        // Overrelaxed step left the previous unbounding sphere, so the interval
        // it skipped is unverified: rewind to that sphere's surface and finish
        // the ray conservatively. The rejected sample updates nothing.
        float back = prev_r - step_len;
        local_p =
            Vector(local_p.x + local_vd.x * back, local_p.y + local_vd.y * back,
                   local_p.z + local_vd.z * back);
        omega = 1.0f;
        prev_r = 0.0f;
        step_len = 0.0f;
        continue;
      }
      prev_r = r;

      if (d < closest_d) {
        closest_d = d;
        closest_local = local_p;
        // A frontal hit converges on the surface from outside, so the
        // d < -aa_width break below never fires for it.
        if (closest_d <= aa_width * 0.02f)
          break;
      } else if (closest_d < aa_width) {
        // Rising past the first in-band local minimum: stop before a surface
        // behind the graze steals the closest approach.
        break;
      }

      if (d < -aa_width)
        break;

      // 1e-5 absolute stall-guard for the precision trace (fine steps near the
      // surface), bounded by max_steps and the early-out above. The probe loop
      // below instead uses a bounds_radius-relative floor for coarse punch-through.
      step_len = std::max(d * 0.9f * omega, 1e-5f);
      local_p = Vector(local_p.x + local_vd.x * step_len,
                       local_p.y + local_vd.y * step_len,
                       local_p.z + local_vd.z * step_len);
    }
    return closest_d;
  }

  /**
   * @brief Result of probing behind an AA halo for an occluded surface.
   */
  struct Occluder {
    bool solid; /**< A solid surface sits behind the halo (behind is valid). */
    Vector behind; /**< Local-space hit point when solid, else the grazed
                        background edge's closest approach when soft > 0. */
    float
        soft; /**< Coverage of a grazed background edge, for the corner fill. */
  };

  /**
   * @brief Marches behind an AA-halo pixel to find any surface the foreground edge
   *        occludes.
   * @param shape Volume shape providing distance().
   * @param closest_local Local-space closest approach (probe seed).
   * @param local_vd Unit ray direction in local space.
   * @param bounds_radius Bounding sphere radius (probe reach + step floor).
   * @param hit_threshold Solid-hit distance threshold.
   * @param aa_width Anti-aliasing band half-width (soft-occlusion falloff scale).
   * @return An Occluder: a solid hit point to antialias the edge over, or a
   *         grazed edge's closest approach and coverage for the corner where
   *         two edges meet.
   */
  template <typename Shape>
  static __attribute__((always_inline)) Occluder probe_occluder(
      const Shape &shape, const Vector &closest_local, const Vector &local_vd,
      float bounds_radius, float hit_threshold, float aa_width) {
    HS_PROFILE_DEEP(vol_probe);
    // March forward from the closest approach for a surface this halo occludes;
    // a solid hit is a self-occlusion edge (antialias over it). Step is floored
    // to punch past the stalled foreground; termination is the bounding sphere's
    // back face. With no solid hit, report a grazed background edge (local min of
    // pd) and its coverage for the corner fill.
    Vector probe = closest_local;
    float prev = FLT_MAX;  // previous step's distance
    bool climbing = false; // pd has risen off the foreground graze
    float min_behind = FLT_MAX;
    Vector min_pos = closest_local;
    // Bracket samples around the running minimum, as offsets along the ray from
    // min_pos, for the parabolic refinement below.
    float s = 0.0f, prev_s = 0.0f, min_s = 0.0f;
    float bef_s = 0.0f, bef_pd = FLT_MAX;
    float aft_s = 0.0f, aft_pd = FLT_MAX;
    bool need_aft = false;
    for (int i = 0; i < PROBE_STEPS; ++i) {
      // Stop at the back of the bounding sphere: nothing left to occlude this halo.
      if (probe.x * local_vd.x + probe.y * local_vd.y + probe.z * local_vd.z >
          bounds_radius)
        break;
      float pd = shape.distance(probe);
      if (pd < hit_threshold)
        return {true, probe, 0.0f}; // solid surface behind the edge
      if (need_aft) {
        aft_s = s;
        aft_pd = pd;
        need_aft = false;
      }
      if (pd > prev)
        climbing = true; // moving away from the foreground graze
      else if (climbing && pd < min_behind) {
        min_behind = pd; // descending toward a surface behind
        min_pos = probe;
        min_s = s;
        bef_s = prev_s;
        bef_pd = prev;
        need_aft = true;
        aft_pd = FLT_MAX;
      }
      prev = pd;
      prev_s = s;
      float floor = bounds_radius *
                    (i < PROBE_NEAR_STEPS ? PROBE_FLOOR_NEAR : PROBE_FLOOR_FAR);
      float step = std::max(pd * 0.9f, floor);
      probe = Vector(probe.x + local_vd.x * step, probe.y + local_vd.y * step,
                     probe.z + local_vd.z * step);
      s += step;
    }

    // The coarse floored stride quantizes the graze minimum, and the sampling
    // phase shifts every frame — corner coverage shimmers under motion. One
    // parabolic-interpolation step through the bracket tightens the minimum
    // (one extra distance eval, graze pixels only) and recovers a thin solid
    // chord the stride stepped over.
    if (min_behind < 2.0f * aa_width && bef_pd != FLT_MAX &&
        aft_pd != FLT_MAX) {
      float p = min_s - bef_s;
      float q = min_s - aft_s;
      float yb = bef_pd - min_behind;
      float ya = aft_pd - min_behind;
      float den = q * yb - p * ya;
      if (den < -1e-12f) {
        float ds = -0.5f * (q * q * yb - p * p * ya) / den;
        if (ds > -p && ds < -q) {
          Vector rp(min_pos.x + local_vd.x * ds, min_pos.y + local_vd.y * ds,
                    min_pos.z + local_vd.z * ds);
          float rpd = shape.distance(rp);
          if (rpd < min_behind) {
            min_behind = rpd;
            min_pos = rp;
          }
          if (min_behind < hit_threshold)
            return {true, min_pos, 0.0f};
        }
      }
    }

    float soft = (min_behind < aa_width)
                     ? quintic_kernel(1.0f - (min_behind - hit_threshold) /
                                                 (aa_width - hit_threshold))
                     : 0.0f;
    return {false, min_pos, soft};
  }

  /**
   * @brief Raymarches and shades a volume shape over its bounding sphere.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Shape Volume shape satisfying the concept below.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param bounds_center Bounding sphere center in physical LED space; must be
   *        a unit vector (on the canvas sphere).
   * @param bounds_radius Bounding sphere radius in world units.
   * @param view_dir Ray direction (camera → scene) in LED space; must point
   *        straight at bounds_center (a radial view). The scanned band is a cap
   *        around bounds_center, which is the orthographic footprint only under
   *        that view; a tilted one slides the footprint off the band and drops
   *        covered columns.
   * @param shape Volume shape providing ray_to_local() and distance().
   * @param frag_fn Fragment shader invoked once per hit.
   * @param max_steps Maximum sphere-tracing steps per ray.
   * @param aa_width Anti-aliasing band half-width in world units.
   * @details Shape concept:
   *   std::pair<Vector, Vector> ray_to_local(const Vector &ro, const Vector
   *   &vd) const; Vector origin_to_local(const Vector &ro) const; float
   *   distance(const Vector &local_point) const;
   *
   * Fragments are plotted at pixel centers by integer coordinates, so the
   * pipeline receives no sub-pixel positions from this draw.
   */
  template <int W, int H, typename Shape>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Vector &bounds_center,
       float bounds_radius, const Vector &view_dir, const Shape &shape,
       FragmentShaderFn frag_fn, int max_steps = 15, float aa_width = 0.01f) {
    check_canvas_dims<W, H>(canvas);
    check_fragment_shader(frag_fn);

    float vd_len = sqrtf(view_dir.x * view_dir.x + view_dir.y * view_dir.y +
                         view_dir.z * view_dir.z);
    float vd_inv = (vd_len > TOLERANCE) ? 1.0f / vd_len : 1.0f;
    Vector vd(view_dir.x * vd_inv, view_dir.y * vd_inv, view_dir.z * vd_inv);

    // Ray must start behind the farthest extent of the shape.
    float start_offset = 1.0f + bounds_radius;

    // bounds_center projected onto the view plane (⊥ vd).
    float bc_dot_vd = bounds_center.x * vd.x + bounds_center.y * vd.y +
                      bounds_center.z * vd.z;
    Vector bc_proj(bounds_center.x - bc_dot_vd * vd.x,
                   bounds_center.y - bc_dot_vd * vd.y,
                   bounds_center.z - bc_dot_vd * vd.z);
    float bounds_r2 = bounds_radius * bounds_radius;

    // Precompute the local-space view direction (shared across all pixels) and,
    // on the same cold transform, validate the volume preconditions. The per-step
    // early-out below compares a local-space dot against the world-space
    // bounds_radius, which holds only if ray_to_local is length-preserving
    // (|local_vd| == 1) and bounds_center maps to the shape's local origin (~0).
    // Trap a scaling shape or off-center bounds_center here, once per draw.
    auto [local_bc, local_vd] = shape.ray_to_local(bounds_center, vd);
    HS_CHECK(fabsf(local_vd.x * local_vd.x + local_vd.y * local_vd.y +
                   local_vd.z * local_vd.z - 1.0f) < TOLERANCE);
    HS_CHECK(local_bc.x * local_bc.x + local_bc.y * local_bc.y +
                 local_bc.z * local_bc.z <
             TOLERANCE);
    // The scan band below is a cap around bounds_center of angular radius
    // asin(bounds_radius), which equals the orthographic footprint only for a
    // radial view of a unit-length center: BoundingSphere reads center.y as
    // cos(phi), and a tilted view slides the footprint off the cap. Unit length
    // also backs the ray start offset above — farther out along the view axis a
    // ray can start in front of the shape.
    HS_CHECK(fabsf(dot(bounds_center, bounds_center) - 1.0f) < TOLERANCE);
    const Vector radial_err = cross(bounds_center, vd);
    HS_CHECK(bc_dot_vd < 0.0f && dot(radial_err, radial_err) < TOLERANCE);
    // aa_width > 0 is the contract: the slow-path AA divides by (aa_width -
    // hit_threshold) == 0.9*aa_width, so a zero band-width gives 0/0 -> NaN.
    HS_CHECK(aa_width > 0.0f);

    BoundingSphere<W, H> bounds(bounds_center, bounds_radius);

    // Tier 2: Clamp volume bounds to clip region
    const auto &cr = canvas.clip();
    const auto vol_xc = cr.x_clip();
    int vol_y_lo, vol_y_hi;
    if (!clamp_rows_to_clip(bounds, cr, vol_y_lo, vol_y_hi))
      return;

    scan_region<W, H>(
        vol_y_lo, vol_y_hi,
        [&](int y, auto &&out) { return bounds.get_intervals(y, out); },
        [&](int px, int py, const Vector &p, int max_run) {
          // Back-face cull
          float facing = p.x * vd.x + p.y * vd.y + p.z * vd.z;
          if (facing >= 0.0f)
            return 1;

          // Orthographic ray-sphere cull
          float pp_x = p.x - facing * vd.x;
          float pp_y = p.y - facing * vd.y;
          float pp_z = p.z - facing * vd.z;
          float dx = pp_x - bc_proj.x;
          float dy = pp_y - bc_proj.y;
          float dz = pp_z - bc_proj.z;
          if (dx * dx + dy * dy + dz * dz > bounds_r2)
            return 1;

          // Orthographic ray origin: outside the unit sphere
          Vector ro(pp_x - vd.x * start_offset, pp_y - vd.y * start_offset,
                    pp_z - vd.z * start_offset);

          // Transform the ray origin to local space once per pixel. The local
          // direction is constant across the draw (local_vd, computed above), so
          // only the origin is transformed here.
          Vector local_ro = shape.origin_to_local(ro);
          Vector closest_local;

          // --- Sphere tracing in local space ---
          float closest_d =
              trace_closest(shape, local_ro, local_vd, bounds_radius, max_steps,
                            aa_width, closest_local);

          if (closest_d >= aa_width) {
            // Shifting the ray origin one column moves every sampled point by
            // at most that chord, and the SDF is 1-Lipschitz, so clearance
            // beyond the offer's own width holds for the whole block. Traced
            // hits stay per-column: shading varies along the surface, so a
            // splat would band.
            if constexpr (pole_lod_blocks<decltype(shape)>) {
              if (max_run > 1) {
                const float block_slack =
                    pole_lod_block_slack<W>(max_run, p.y, shape);
                if (pole_lod_block_settles<decltype(shape)>(closest_d, aa_width,
                                                            block_slack))
                  return max_run;
              }
            }
            return 1;
          }

          // --- Fragment shading ---
          Fragment frag;
          frag.pos = closest_local;
          frag.size = closest_d;
          {
            HS_PROFILE_DEEP(vol_shade);
            frag_fn(closest_local, frag);
          }

          // One-sided AA with quintic kernel
          float hit_threshold = aa_width * 0.1f;
          float edge_alpha;

          if (closest_d <= hit_threshold) {
            // FAST PATH: Solid hit. No probe needed.
            edge_alpha = 1.0f;
          } else {
            // SLOW PATH: fuzzy AA border. Standard one-sided AA coverage...
            edge_alpha = quintic_kernel(1.0f - (closest_d - hit_threshold) /
                                                   (aa_width - hit_threshold));

            // ...then probe behind the halo for a surface this edge occludes.
            Occluder occ =
                probe_occluder(shape, closest_local, local_vd, bounds_radius,
                               hit_threshold, aa_width);
            if (occ.solid) {
              // Self-occlusion edge: antialias the foreground over the surface it
              // covers — lay the shaded background down, then blend the foreground
              // over it by the edge coverage. Smooth, vs. fading to black (fringe)
              // or snapping to opaque (jagged).
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              {
                HS_PROFILE_DEEP(vol_shade);
                frag_fn(occ.behind, bg);
              }
              if (bg.color.alpha > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, bg.color.color, 0.0f,
                              bg.color.alpha);
              }
              if (frag.color.alpha * edge_alpha > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, frag.color.color, 0.0f,
                              frag.color.alpha * edge_alpha);
              }
              return 1;
            }
            // No solid behind; a grazed background edge fills the corner,
            // shaded at its own point so the fill carries the background's
            // color, then the foreground blends over it.
            if (occ.soft > MIN_ALPHA) {
              Fragment bg;
              bg.pos = occ.behind;
              bg.size = 0.0f;
              {
                HS_PROFILE_DEEP(vol_shade);
                frag_fn(occ.behind, bg);
              }
              if (bg.color.alpha * occ.soft > MIN_ALPHA) {
                HS_PROFILE_DEEP(vol_plot);
                pipeline.plot(canvas, px, py, bg.color.color, 0.0f,
                              bg.color.alpha * occ.soft);
              }
            }
          }

          if (frag.color.alpha * edge_alpha > MIN_ALPHA) {
            HS_PROFILE_DEEP(vol_plot);
            pipeline.plot(canvas, px, py, frag.color.color, 0.0f,
                          frag.color.alpha * edge_alpha);
          }
          return 1;
        },
        vol_xc);
  }
};
HS_O3_END

} // namespace Scan
