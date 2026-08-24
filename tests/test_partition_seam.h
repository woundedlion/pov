/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Partition-op rasterizer seam calibration
 * (docs/specs/opchain_morph_spec.md, "Validation contract"). Measures, at the
 * shipping canvas size, the framebuffer delta a `kis` or `dual` swap produces
 * under a single flat fill colour, which isolates the
 * coverage discontinuity from the shading gradient.
 *
 * Doubles as the gated swap's seam regression: every swap is asserted against
 * the envelope this calibration fixed. Set HS_SEAM_DUMP=<dir> to write the
 * captures as PNG and to run the non-asserting gradient sweep.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "core/animation/opleg.h"
#include "core/mesh/conway.h"
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "core/render/canvas.h"
#include "core/render/scan.h"
#include "core/render/shading.h"
#include "tests/mesh_test_util.h"
#include "tests/pixel_test_util.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace partition_seam_tests {

constexpr int PS_W = 288; /**< IslamicStars canvas width. */
constexpr int PS_H = 144; /**< IslamicStars canvas height. */

inline uint8_t ps_geom_buf[512 * 1024];
inline uint8_t ps_temp_buf[512 * 1024];
inline uint8_t ps_aux_buf[512 * 1024];
inline uint8_t ps_scan_buf[256 * 1024];

/** Flat fill: full-scale grey, so a channel delta reads directly as a fraction
 * of the parent's coverage. */
constexpr uint16_t FILL = 65535;

/** Changed-pixel threshold, 1 % of full scale. */
constexpr int DELTA_THRESH = 655;

/** Deep-seam threshold, 20 % of full scale. */
constexpr int DEEP_THRESH = FILL / 5;

/** Radius in pixels of the parent-vertex neighbourhood. */
constexpr int VERT_RADIUS = 3;

/**
 * @brief Renders one compiled mesh through a fragment shader and captures the
 *        frame.
 * @param out Receives the PS_W x PS_H frame, row-major.
 * @param mesh Compiled mesh.
 * @param shader Fragment shader body, `void(const Vector &, Fragment &)`.
 */
template <typename ShaderFn>
inline void render(std::vector<Pixel> &out, const MeshState &mesh,
                   ShaderFn &&shader) {
  hs_test::StubEffect fx(PS_W, PS_H);
  {
    Canvas c(fx);
    Pipeline<PS_W, PS_H> pipe;
    Arena scan_scratch(ps_scan_buf, sizeof(ps_scan_buf));
    Scan::Mesh::draw<PS_W, PS_H>(pipe, c, mesh, shader, scan_scratch);
  }
  fx.advance_display();
  capture_frame<PS_W, PS_H>(fx, out);
}

/**
 * @brief Renders a mesh filled flat with FILL on every face.
 */
inline void render_flat(std::vector<Pixel> &out, const MeshState &mesh) {
  render(out, mesh, [](const Vector &, Fragment &f) {
    f.color = Color4(Pixel(FILL, FILL, FILL), 1.0f);
  });
}

/**
 * @brief Renders a mesh with the shipping inset gradient (grey ramp, no
 *        palette) at the given gain.
 */
inline void render_gradient(std::vector<Pixel> &out, const MeshState &mesh,
                            float gain) {
  render(out, mesh, [gain](const Vector &, Fragment &f) {
    const float t = hs::clamp(fragment_edge_dist(f) * gain, 0.0f, 1.0f);
    const uint16_t v = static_cast<uint16_t>(t * FILL);
    f.color = Color4(Pixel(v, v, v), 1.0f);
  });
}

/**
 * @brief Green channel of a capture, the measured channel under a grey fill.
 */
inline int lum(const Pixel &p) { return static_cast<int>(p.g); }

/**
 * @brief Statistics of one before/after capture pair.
 */
struct SeamStats {
  size_t lit_a = 0;          /**< Non-black pixels in the first capture. */
  size_t lit_b = 0;          /**< Non-black pixels in the second capture. */
  size_t full_a = 0;         /**< Pixels at full fill in the first capture. */
  size_t changed = 0;        /**< Pixels differing by more than DELTA_THRESH. */
  size_t changed_clean = 0;  /**< Changed, and full-fill before (new dark). */
  size_t changed_edge = 0;   /**< Changed, and already dark before. */
  size_t changed_near_v = 0; /**< Changed, within VERT_RADIUS of a parent
                                  vertex. */
  float max_dark = 0.0f;     /**< Largest darkening, fraction of full. */
  float max_bright = 0.0f;   /**< Largest brightening, fraction of full. */
  double mean_abs = 0.0;     /**< Mean |delta| over changed pixels. */
  double mean_near = 0.0;    /**< Mean |delta| near a parent vertex. */
  double mean_far = 0.0;     /**< Mean |delta| elsewhere. */
  double energy = 0.0;       /**< Sum(a-b) over the canvas, / (pixels*full). */
  double abs_energy = 0.0;   /**< Sum|a-b| over changed pixels only,
                                  / (pixels*full). */
  size_t deep = 0;           /**< Pixels differing by >= DEEP_THRESH. */
  double mean_clean = 0.0;   /**< Mean |delta| over changed clean-interior. */
  double mean_edge = 0.0;    /**< Mean |delta| over changed existing-edge. */
  float max_clean = 0.0f;    /**< Max |delta| over changed clean-interior. */
  float max_edge = 0.0f;     /**< Max |delta| over changed existing-edge. */
  double mean_band = 0.0;    /**< Mean min(row-run, column-run) of changed. */
  int max_band = 0;          /**< Max min(row-run, column-run) of changed. */
};

/**
 * @brief Screen pixel of a direction, by nearest LUT sample.
 * @param v Unit direction.
 * @param px Receives the column.
 * @param py Receives the row.
 */
inline void project(const Vector &v, int &px, int &py) {
  float best = -2.0f;
  px = py = 0;
  for (int y = 0; y < PS_H; ++y) {
    const float sp = TrigLUT<PS_W, PS_H>::sin_phi[y];
    const float cp = TrigLUT<PS_W, PS_H>::cos_phi[y];
    for (int x = 0; x < PS_W; ++x) {
      const Vector p(sp * TrigLUT<PS_W, PS_H>::cos_theta(x), cp,
                     sp * TrigLUT<PS_W, PS_H>::sin_theta[x]);
      const float d = dot(p, v);
      if (d > best) {
        best = d;
        px = x;
        py = y;
      }
    }
  }
}

/**
 * @brief Compares two captures and reports the seam statistics.
 * @param a Before capture.
 * @param b After capture.
 * @param verts Parent vertices, for the near-vertex split; may be empty.
 */
inline SeamStats compare(const std::vector<Pixel> &a,
                         const std::vector<Pixel> &b,
                         const std::vector<Vector> &verts) {
  SeamStats st;
  const size_t n = a.size();
  std::vector<uint8_t> changed(n, 0);

  double sum_abs = 0.0;
  long long sum_signed = 0;
  for (size_t i = 0; i < n; ++i) {
    if (lum(a[i]) > 0)
      ++st.lit_a;
    if (lum(b[i]) > 0)
      ++st.lit_b;
    const bool full = lum(a[i]) >= FILL - DELTA_THRESH;
    if (full)
      ++st.full_a;
    const int d = lum(a[i]) - lum(b[i]);
    sum_signed += d;
    if (std::abs(d) <= DELTA_THRESH)
      continue;
    changed[i] = 1;
    ++st.changed;
    const float ad = std::abs(d) / float(FILL);
    if (full) {
      ++st.changed_clean;
      st.mean_clean += ad;
      st.max_clean = std::max(st.max_clean, ad);
    } else {
      ++st.changed_edge;
      st.mean_edge += ad;
      st.max_edge = std::max(st.max_edge, ad);
    }
    sum_abs += std::abs(d);
    if (std::abs(d) >= DEEP_THRESH)
      ++st.deep;
    st.max_dark = std::max(st.max_dark, d / float(FILL));
    st.max_bright = std::max(st.max_bright, -d / float(FILL));
  }
  st.mean_abs = st.changed ? sum_abs / double(st.changed) / FILL : 0.0;
  st.energy = double(sum_signed) / double(n) / FILL;
  st.abs_energy = sum_abs / double(n) / FILL;
  if (st.changed_clean)
    st.mean_clean /= double(st.changed_clean);
  if (st.changed_edge)
    st.mean_edge /= double(st.changed_edge);

  // Band width: shorter of the changed pixel's row run (theta wraps) and its
  // column run.
  double band_sum = 0.0;
  for (int y = 0; y < PS_H; ++y) {
    for (int x = 0; x < PS_W; ++x) {
      const size_t i = size_t(y) * PS_W + x;
      if (!changed[i])
        continue;
      int hrun = 1;
      for (int k = 1; k < PS_W; ++k) {
        const int xx = (x + k) % PS_W;
        if (!changed[size_t(y) * PS_W + xx])
          break;
        ++hrun;
      }
      // Bound by the pixels the forward scan did not already take, so a fully
      // changed row reports PS_W rather than counting the wrap twice.
      const int back_limit = PS_W - hrun;
      for (int k = 1; k <= back_limit; ++k) {
        const int xx = (x - k + PS_W * 2) % PS_W;
        if (!changed[size_t(y) * PS_W + xx])
          break;
        ++hrun;
      }
      int vrun = 1;
      for (int k = 1; y + k < PS_H; ++k) {
        if (!changed[size_t(y + k) * PS_W + x])
          break;
        ++vrun;
      }
      for (int k = 1; y - k >= 0; ++k) {
        if (!changed[size_t(y - k) * PS_W + x])
          break;
        ++vrun;
      }
      const int band = std::min(hrun, vrun);
      band_sum += band;
      st.max_band = std::max(st.max_band, band);
    }
  }
  st.mean_band = st.changed ? band_sum / double(st.changed) : 0.0;

  if (verts.empty())
    return st;

  std::vector<uint8_t> near_v(n, 0);
  for (const Vector &v : verts) {
    int px = 0, py = 0;
    project(v, px, py);
    for (int dy = -VERT_RADIUS; dy <= VERT_RADIUS; ++dy) {
      const int y = py + dy;
      if (y < 0 || y >= PS_H)
        continue;
      for (int dx = -VERT_RADIUS; dx <= VERT_RADIUS; ++dx) {
        if (dx * dx + dy * dy > VERT_RADIUS * VERT_RADIUS)
          continue;
        const int x = (px + dx + PS_W * 2) % PS_W;
        near_v[size_t(y) * PS_W + x] = 1;
      }
    }
  }
  double near_sum = 0.0, far_sum = 0.0;
  size_t far_n = 0;
  for (size_t i = 0; i < n; ++i) {
    if (!changed[i])
      continue;
    const double d = std::abs(lum(a[i]) - lum(b[i])) / double(FILL);
    if (near_v[i]) {
      ++st.changed_near_v;
      near_sum += d;
    } else {
      ++far_n;
      far_sum += d;
    }
  }
  st.mean_near = st.changed_near_v ? near_sum / double(st.changed_near_v) : 0.0;
  st.mean_far = far_n ? far_sum / double(far_n) : 0.0;
  return st;
}

/** Gated-swap envelope. The spec's "<= 2 % of pixels" (section 9.5) is
 * unreachable at this canvas size: the swap's irreducible coverage delta is a
 * 2-4 px band along the child edges, 6.6-15.3 % of the canvas, so the bound is
 * on the changed fraction, the changed pixels' absolute energy and the deepest
 * pixel instead. Every bound is two-sided around its measured value: a widened
 * band and a collapsed seam — children that no longer partition the parent, and
 * so leave the capture untouched — are both out of envelope. Each swap's
 * changed fraction is bracketed by its own measured value plus or minus two
 * percentage points, enough for rounding drift either way. Fraction and energy
 * both count only pixels past DELTA_THRESH, so a uniform sub-threshold shift is
 * invisible to the gate. */
constexpr double CHANGED_FRAC_MARGIN = 0.02;
constexpr double MEASURED_CHANGED_FRAC_KIS_ICOSA = 0.1533;
constexpr double MEASURED_CHANGED_FRAC_KIS_CUBE = 0.1019;
constexpr double MEASURED_CHANGED_FRAC_KIS_DODECA = 0.1448;
constexpr double MEASURED_CHANGED_FRAC_DUAL_ICOSA = 0.1217;
constexpr double MEASURED_CHANGED_FRAC_DUAL_CUBE = 0.0662;
constexpr double MEASURED_CHANGED_FRAC_DUAL_DODECA = 0.1217;
constexpr double MAX_ABS_ENERGY = 0.02;

/** Half the smallest absolute energy measured over the six swaps (0.96 %, dual
 * cube). A seam that stopped moving pixels reads as zero energy. */
constexpr double MIN_ABS_ENERGY = 0.005;

/** Deepest pixel measured over the six calibration swaps, the same under
 * -ffast-math as under IEEE: a seam pixel that both children claim at half
 * coverage composites to 3/4 of the parent's fill. A child that loses its
 * share of such a pixel leaves half, so the margin stays well under 0.5. All
 * six swaps reach it exactly, so the darkening is bracketed on both sides. */
constexpr float MAX_MEASURED_PIXEL_DELTA = 0.25f;
constexpr float PIXEL_DELTA_MARGIN = 0.05f;
constexpr float MAX_PIXEL_DELTA = MAX_MEASURED_PIXEL_DELTA + PIXEL_DELTA_MARGIN;
constexpr float MIN_PIXEL_DELTA = MAX_MEASURED_PIXEL_DELTA - PIXEL_DELTA_MARGIN;

/**
 * @brief Asserts one swap's statistics against the gated-swap envelope.
 * @param st Statistics of the swap.
 * @param measured_changed_frac The swap's calibrated changed fraction; its
 *        bracket is that plus or minus CHANGED_FRAC_MARGIN.
 */
inline void expect_within_envelope(const SeamStats &st,
                                   double measured_changed_frac) {
  const double frac = st.changed / (double(PS_W) * PS_H);
  HS_EXPECT_LE(frac, measured_changed_frac + CHANGED_FRAC_MARGIN);
  HS_EXPECT_GE(frac, measured_changed_frac - CHANGED_FRAC_MARGIN);
  HS_EXPECT_LE(st.abs_energy, MAX_ABS_ENERGY);
  HS_EXPECT_GE(st.abs_energy, MIN_ABS_ENERGY);
  HS_EXPECT_LE(st.max_dark, MAX_PIXEL_DELTA);
  HS_EXPECT_GE(st.max_dark, MIN_PIXEL_DELTA);
  HS_EXPECT_LE(st.max_bright, MAX_PIXEL_DELTA);
}

/**
 * @brief Prints one statistics row.
 */
inline void report(const char *tag, const SeamStats &st) {
  const double n = double(PS_W) * PS_H;
  std::printf("  [%s] changed %zu px (%.2f%% of canvas, %.2f%% of lit); "
              "max dark %.1f%% bright %.1f%%; mean |d| %.1f%%; "
              "band mean %.2f max %d px\n",
              tag, st.changed, 100.0 * st.changed / n,
              st.lit_a ? 100.0 * st.changed / double(st.lit_a) : 0.0,
              100.0 * st.max_dark, 100.0 * st.max_bright, 100.0 * st.mean_abs,
              st.mean_band, st.max_band);
  std::printf("      lit before %zu / after %zu; full-fill before %zu; "
              "changed on clean interior %zu, on existing edge %zu; "
              ">=%.0f%% deep %zu px; signed energy %.4f%%, abs energy %.4f%%\n",
              st.lit_a, st.lit_b, st.full_a, st.changed_clean, st.changed_edge,
              100.0 * DEEP_THRESH / FILL, st.deep, 100.0 * st.energy,
              100.0 * st.abs_energy);
  std::printf("      clean interior: mean |d| %.1f%% max %.1f%%; existing "
              "edge: mean |d| %.1f%% max %.1f%%\n",
              100.0 * st.mean_clean, 100.0 * st.max_clean, 100.0 * st.mean_edge,
              100.0 * st.max_edge);
  if (st.changed_near_v || st.mean_far > 0.0)
    std::printf("      near parent vertex (<=%d px): %zu changed, mean |d| "
                "%.1f%% vs %.1f%% elsewhere\n",
                VERT_RADIUS, st.changed_near_v, 100.0 * st.mean_near,
                100.0 * st.mean_far);
}

// ---------------------------------------------------------------------------
// PNG dump (HS_SEAM_DUMP=<dir>), stored-deflate, 8-bit RGB.
// ---------------------------------------------------------------------------

/**
 * @brief CRC-32 of a byte range, PNG polynomial.
 */
inline uint32_t png_crc(const uint8_t *p, size_t n,
                        uint32_t crc = 0xFFFFFFFFu) {
  for (size_t i = 0; i < n; ++i) {
    crc ^= p[i];
    for (int k = 0; k < 8; ++k)
      crc = (crc >> 1) ^ (0xEDB88320u & (~(crc & 1) + 1));
  }
  return crc;
}

/**
 * @brief Appends a PNG chunk with its length, type, payload and CRC.
 */
inline void png_chunk(std::vector<uint8_t> &out, const char type[5],
                      const std::vector<uint8_t> &data) {
  const uint32_t n = static_cast<uint32_t>(data.size());
  for (int i = 3; i >= 0; --i)
    out.push_back(static_cast<uint8_t>(n >> (8 * i)));
  std::vector<uint8_t> body(type, type + 4);
  body.insert(body.end(), data.begin(), data.end());
  out.insert(out.end(), body.begin(), body.end());
  const uint32_t crc = png_crc(body.data(), body.size()) ^ 0xFFFFFFFFu;
  for (int i = 3; i >= 0; --i)
    out.push_back(static_cast<uint8_t>(crc >> (8 * i)));
}

/**
 * @brief Dump directory named by HS_SEAM_DUMP, or nullptr when unset.
 */
inline const char *seam_dump_dir() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  return std::getenv("HS_SEAM_DUMP");
#pragma clang diagnostic pop
}

/**
 * @brief Writes a capture to <HS_SEAM_DUMP>/<name>.png; no-op when unset.
 */
inline void dump_png(const char *name, const std::vector<Pixel> &px) {
  const char *dir = seam_dump_dir();
  if (!dir)
    return;

  std::vector<uint8_t> raw;
  raw.reserve(size_t(PS_H) * (1 + 3 * PS_W));
  for (int y = 0; y < PS_H; ++y) {
    raw.push_back(0);
    for (int x = 0; x < PS_W; ++x) {
      const Pixel &p = px[size_t(y) * PS_W + x];
      raw.push_back(static_cast<uint8_t>(p.r >> 8));
      raw.push_back(static_cast<uint8_t>(p.g >> 8));
      raw.push_back(static_cast<uint8_t>(p.b >> 8));
    }
  }

  // zlib stream with stored deflate blocks.
  std::vector<uint8_t> z{0x78, 0x01};
  for (size_t off = 0; off < raw.size();) {
    const size_t len = std::min<size_t>(65535, raw.size() - off);
    z.push_back(off + len >= raw.size() ? 1 : 0);
    z.push_back(static_cast<uint8_t>(len & 0xFF));
    z.push_back(static_cast<uint8_t>(len >> 8));
    z.push_back(static_cast<uint8_t>(~len & 0xFF));
    z.push_back(static_cast<uint8_t>((~len >> 8) & 0xFF));
    z.insert(z.end(), raw.begin() + off, raw.begin() + off + len);
    off += len;
  }
  uint32_t s1 = 1, s2 = 0;
  for (uint8_t b : raw) {
    s1 = (s1 + b) % 65521;
    s2 = (s2 + s1) % 65521;
  }
  const uint32_t adler = (s2 << 16) | s1;
  for (int i = 3; i >= 0; --i)
    z.push_back(static_cast<uint8_t>(adler >> (8 * i)));

  std::vector<uint8_t> ihdr;
  for (int i = 3; i >= 0; --i)
    ihdr.push_back(static_cast<uint8_t>(uint32_t(PS_W) >> (8 * i)));
  for (int i = 3; i >= 0; --i)
    ihdr.push_back(static_cast<uint8_t>(uint32_t(PS_H) >> (8 * i)));
  ihdr.insert(ihdr.end(), {8, 2, 0, 0, 0});

  std::vector<uint8_t> png{0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
  png_chunk(png, "IHDR", ihdr);
  png_chunk(png, "IDAT", z);
  png_chunk(png, "IEND", {});

  char path[512];
  std::snprintf(path, sizeof(path), "%s/%s.png", dir, name);
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  FILE *f = std::fopen(path, "wb");
#pragma clang diagnostic pop
  if (f) {
    std::fwrite(png.data(), 1, png.size(), f);
    std::fclose(f);
  }
}

// ---------------------------------------------------------------------------
// Gate 6 — kis fan seam, per seed face degree.
// ---------------------------------------------------------------------------

/**
 * @brief Amplified difference image: red for darkening, blue for brightening,
 *        both scaled 4x.
 */
inline std::vector<Pixel> diff_image(const std::vector<Pixel> &a,
                                     const std::vector<Pixel> &b) {
  std::vector<Pixel> d(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    const int v = std::min<int>(FILL, std::abs(lum(a[i]) - lum(b[i])) * 4);
    const uint16_t m = static_cast<uint16_t>(v);
    d[i] = lum(a[i]) >= lum(b[i]) ? Pixel(m, 0, 0) : Pixel(0, 0, m);
  }
  return d;
}

/**
 * @brief Copies a mesh's vertices out for the near-vertex split.
 */
inline std::vector<Vector> vertex_list(const PolyMesh &m) {
  std::vector<Vector> v;
  for (size_t i = 0; i < m.vertices.size(); ++i)
    v.push_back(m.vertices[i]);
  return v;
}

/**
 * @brief Measures the flat-shaded delta of `kis` on one seed.
 * @tparam Solid Seed solid descriptor.
 * @param name Tag for the report and PNG names.
 * @param measured_changed_frac This swap's calibrated changed fraction.
 */
template <typename Solid>
inline void measure_kis(const char *name, double measured_changed_frac) {
  Arena geom(ps_geom_buf, sizeof(ps_geom_buf));
  Arena temp(ps_temp_buf, sizeof(ps_temp_buf));
  Arena aux(ps_aux_buf, sizeof(ps_aux_buf));

  PolyMesh base;
  build_solid<Solid>(base, geom);
  PolyMesh kis = MeshOps::kis(base, geom, temp);

  MeshState base_ms, kis_ms;
  MeshOps::compile(base, base_ms, aux, temp);
  MeshOps::compile(kis, kis_ms, aux, temp);

  std::vector<Pixel> a, b;
  render_flat(a, base_ms);
  render_flat(b, kis_ms);

  char tag[128];
  std::snprintf(tag, sizeof(tag), "kis %s deg=%d F %zu->%zu", name,
                base.face_counts[0], base.face_counts.size(),
                kis.face_counts.size());
  const SeamStats st = compare(a, b, vertex_list(base));
  report(tag, st);

  char png[128];
  std::snprintf(png, sizeof(png), "%s_base_flat", name);
  dump_png(png, a);
  std::snprintf(png, sizeof(png), "%s_kis_flat", name);
  dump_png(png, b);
  std::snprintf(png, sizeof(png), "%s_kis_diff", name);
  dump_png(png, diff_image(a, b));

  // The fan partitions the parent's spherical patches, so total coverage is
  // unchanged; the delta is dominated by darkening on the new spokes, with a
  // smaller brightening where the children's own gnomonic frames narrow the
  // parent's existing edge bands.
  HS_EXPECT_EQ(st.lit_a, st.lit_b);
  HS_EXPECT_GT(st.energy, 0.0);
  HS_EXPECT_GT(st.max_dark, st.max_bright);
  expect_within_envelope(st, measured_changed_frac);
}

/**
 * @brief Measures the flat-shaded delta of `dual` on one seed.
 * @tparam Solid Seed solid descriptor.
 * @param name Tag for the report and PNG names.
 * @param measured_changed_frac This swap's calibrated changed fraction.
 */
template <typename Solid>
inline void measure_dual(const char *name, double measured_changed_frac) {
  Arena geom(ps_geom_buf, sizeof(ps_geom_buf));
  Arena temp(ps_temp_buf, sizeof(ps_temp_buf));
  Arena aux(ps_aux_buf, sizeof(ps_aux_buf));

  PolyMesh base;
  build_solid<Solid>(base, geom);
  PolyMesh dual = MeshOps::dual(base, geom, temp);

  MeshState base_ms, dual_ms;
  MeshOps::compile(base, base_ms, aux, temp);
  MeshOps::compile(dual, dual_ms, aux, temp);

  std::vector<Pixel> a, b;
  render_flat(a, base_ms);
  render_flat(b, dual_ms);

  char tag[128];
  std::snprintf(tag, sizeof(tag), "dual %s F %zu->%zu", name,
                base.face_counts.size(), dual.face_counts.size());
  const SeamStats st = compare(a, b, vertex_list(base));
  report(tag, st);

  char png[128];
  std::snprintf(png, sizeof(png), "%s_dual_flat", name);
  dump_png(png, b);
  std::snprintf(png, sizeof(png), "%s_dual_diff", name);
  dump_png(png, diff_image(a, b));

  // The dual's faces partition the sphere the same way the seed's do, so total
  // coverage is unchanged; unlike kis the delta has no fixed sign, since the
  // new edges neither contain nor are contained by the old ones.
  HS_EXPECT_EQ(st.lit_a, st.lit_b);
  expect_within_envelope(st, measured_changed_frac);
}

/**
 * @brief Measures the gradient's contribution on top of the coverage seam, at
 *        the gains the §3.3 gate would open and close between.
 * @details Reporting only, no assertions; run under HS_SEAM_DUMP.
 * @tparam Solid Seed solid descriptor.
 */
template <typename Solid> inline void measure_gradient(const char *name) {
  Arena geom(ps_geom_buf, sizeof(ps_geom_buf));
  Arena temp(ps_temp_buf, sizeof(ps_temp_buf));
  Arena aux(ps_aux_buf, sizeof(ps_aux_buf));

  PolyMesh base;
  build_solid<Solid>(base, geom);
  PolyMesh kis = MeshOps::kis(base, geom, temp);
  PolyMesh dual = MeshOps::dual(base, geom, temp);

  MeshState base_ms, kis_ms, dual_ms;
  MeshOps::compile(base, base_ms, aux, temp);
  MeshOps::compile(kis, kis_ms, aux, temp);
  MeshOps::compile(dual, dual_ms, aux, temp);

  const float GAINS[] = {1.0f, 0.5f, 0.25f, 0.15f, 0.1f, 0.05f, 0.02f};
  for (float g : GAINS) {
    std::vector<Pixel> a, b, c;
    render_gradient(a, base_ms, g);
    render_gradient(b, kis_ms, g);
    render_gradient(c, dual_ms, g);
    char tag[128];
    std::snprintf(tag, sizeof(tag), "kis %s gradient gain=%.2f", name, g);
    report(tag, compare(a, b, {}));
    std::snprintf(tag, sizeof(tag), "dual %s gradient gain=%.2f", name, g);
    report(tag, compare(a, c, {}));
    if (g == 1.0f) {
      char png[128];
      std::snprintf(png, sizeof(png), "%s_base_gain1", name);
      dump_png(png, a);
      std::snprintf(png, sizeof(png), "%s_kis_gain1", name);
      dump_png(png, b);
      std::snprintf(png, sizeof(png), "%s_dual_gain1", name);
      dump_png(png, c);
    }
  }
}

/**
 * @brief Gate 6: flat-shaded kis and dual seam calibration over three seeds
 *        spanning face degree (3, 4, 5).
 */
inline void test_partition_seam_calibration() {
  std::printf("  [gate6] canvas %dx%d, flat fill, threshold %.1f%% of full\n",
              PS_W, PS_H, 100.0 * DELTA_THRESH / FILL);
  measure_kis<Solids::Icosahedron>("icosa", MEASURED_CHANGED_FRAC_KIS_ICOSA);
  measure_kis<Solids::Cube>("cube", MEASURED_CHANGED_FRAC_KIS_CUBE);
  measure_kis<Solids::Dodecahedron>("dodeca", MEASURED_CHANGED_FRAC_KIS_DODECA);

  measure_dual<Solids::Icosahedron>("icosa", MEASURED_CHANGED_FRAC_DUAL_ICOSA);
  measure_dual<Solids::Cube>("cube", MEASURED_CHANGED_FRAC_DUAL_CUBE);
  measure_dual<Solids::Dodecahedron>("dodeca",
                                     MEASURED_CHANGED_FRAC_DUAL_DODECA);

  if (seam_dump_dir())
    measure_gradient<Solids::Dodecahedron>("dodeca");
}

/**
 * @brief Runs the partition-seam calibration.
 * @return The module's failure count.
 */
inline int run_partition_seam_tests() {
  hs_test::ModuleFixture fixture("partition_seam");
  test_partition_seam_calibration();
  return fixture.result();
}

} // namespace partition_seam_tests
} // namespace hs_test
