/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 *
 * Host driver for the mesh-scan AA coverage audit (HS_AA_AUDIT).
 *
 * For every Islamic-registry solid over a spread of orientations it
 *   - brute-forces each scanned row to find shaded pixels the emitted column
 *     runs never visit, and
 *   - renders the same frame twice (emitted runs vs. forced full-width scan)
 *     and diffs the framebuffers, so neighbour-face coverage is accounted for.
 * Optionally dumps both frames as raw RGB16 for PNG diffing. Not a CTest.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/engine/memory.h"
#include "core/mesh/mesh.h"
#include "core/mesh/solids.h"
#include "core/render/canvas.h"
#include "core/render/scan.h"
#include "tests/aa_audit.h"
#include "tests/test_fixture.h"

namespace {

constexpr int W = 288, H = 144;

alignas(32) uint8_t g_seed_a[3 * 1024 * 1024];
alignas(32) uint8_t g_seed_b[3 * 1024 * 1024];
alignas(32) uint8_t g_geom[3 * 1024 * 1024];
alignas(32) uint8_t g_scratch[3 * 1024 * 1024];

struct MeshFx : public Effect {
  MeshFx() : Effect(W, H) {}
  void draw_frame() override {}
};

bool g_face_color = false;
bool g_no_ref = false;

/** @brief Flat white, or a per-face hue cycle when g_face_color is set. */
void white(const Vector &, Fragment &f) {
  if (!g_face_color) {
    f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
    return;
  }
  int i = static_cast<int>(f.v2);
  static const uint16_t lut[6][3] = {{65535, 0, 0},     {0, 65535, 0},
                                     {0, 0, 65535},     {65535, 65535, 0},
                                     {0, 65535, 65535}, {65535, 0, 65535}};
  const uint16_t *c = lut[i % 6];
  f.color = Color4(Pixel(c[0], c[1], c[2]), 1.0f);
}

/** @brief Rotates v by yaw about Y then pitch about X. */
Vector rotate(const Vector &v, float yaw, float pitch) {
  float cy = cosf(yaw), sy = sinf(yaw);
  float x = v.x * cy - v.z * sy;
  float z = v.x * sy + v.z * cy;
  float cp = cosf(pitch), sp = sinf(pitch);
  float y = v.y * cp - z * sp;
  float z2 = v.y * sp + z * cp;
  return Vector(x, y, z2);
}

void render(MeshState &mesh, Arena &scratch, std::vector<uint16_t> &out) {
  MeshFx fx;
  {
    Canvas c(fx);
    Pipeline<W, H> pipe;
    scratch.reset();
    Scan::Mesh::draw<W, H>(pipe, c, mesh, white, scratch);
  }
  fx.advance_display();
  out.resize(static_cast<size_t>(W) * H * 3);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      Pixel p = fx.get_pixel(x, y);
      size_t i = (static_cast<size_t>(y) * W + x) * 3;
      out[i] = p.r;
      out[i + 1] = p.g;
      out[i + 2] = p.b;
    }
}

void dump(const char *path, const std::vector<uint16_t> &buf) {
  std::FILE *f = std::fopen(path, "wb");
  std::fwrite(buf.data(), 2, buf.size(), f);
  std::fclose(f);
}

} // namespace

int main(int argc, char **argv) {
  const char *only = nullptr;
  const char *dump_prefix = nullptr;
  int orientations = 8;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--solid") && i + 1 < argc)
      only = argv[++i];
    else if (!strcmp(argv[i], "--dump") && i + 1 < argc)
      dump_prefix = argv[++i];
    else if (!strcmp(argv[i], "--orientations") && i + 1 < argc)
      orientations = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--noref"))
      g_no_ref = true;
    else if (!strcmp(argv[i], "--legacy"))
      hs_aa::g_audit.legacy_pad = true;
    else if (!strcmp(argv[i], "--facecolor"))
      g_face_color = true;
  }

  Arena seed_a(g_seed_a, sizeof(g_seed_a));
  Arena seed_b(g_seed_b, sizeof(g_seed_b));
  Arena geom(g_geom, sizeof(g_geom));
  Arena scratch(g_scratch, sizeof(g_scratch));

  hs_test::reset_globals();

  long long tot_probes = 0, tot_ref_probes = 0, tot_painted = 0,
            tot_missed = 0, tot_diff = 0;
  int draws = 0;

  std::printf("%-58s %5s %10s %8s %8s %8s %7s\n", "solid", "F", "probes",
              "painted", "missed", "px/draw", "fbdiff");

  std::vector<uint16_t> a_buf, b_buf;

  for (const auto &entry : Solids::islamic_registry) {
    if (only && strcmp(entry.name, only))
      continue;
    seed_a.reset();
    seed_b.reset();
    geom.reset();
    PolyMesh poly = entry.generate(seed_a, seed_b);
    MeshState mesh;
    MeshOps::compile(poly, mesh, geom, scratch);

    hs_aa::g_audit.reset();

    long long fb_diff = 0, max_delta = 0, ref_probes = 0, holes = 0,
              big_delta = 0;

    for (int o = 0; o < orientations; ++o) {
      float yaw = 0.37f * static_cast<float>(o);
      float pitch = 0.21f * static_cast<float>(o) + 0.05f;
      std::vector<Vector> saved(mesh.vertices.size());
      for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        saved[i] = mesh.vertices[i];
        mesh.vertices[i] = rotate(saved[i], yaw, pitch);
      }

      hs_aa::g_audit.enabled = true;
      hs_aa::g_audit.full_scan = false;
      render(mesh, scratch, a_buf);
      hs_aa::g_audit.enabled = false;

      long long before = hs_aa::g_audit.probes;
      if (!g_no_ref) {
      hs_aa::g_audit.enabled = true;
      hs_aa::g_audit.full_scan = true;
      render(mesh, scratch, b_buf);
      hs_aa::g_audit.enabled = false;
      hs_aa::g_audit.full_scan = false;
      ref_probes += hs_aa::g_audit.probes - before;
      hs_aa::g_audit.probes = before;
      }

      for (size_t p = 0; !g_no_ref && p < a_buf.size(); p += 3) {
        long long d = 0;
        for (int c = 0; c < 3; ++c) {
          long long dc =
              std::abs(static_cast<long long>(a_buf[p + c]) - b_buf[p + c]);
          if (dc > d)
            d = dc;
        }
        if (!d)
          continue;
        ++fb_diff;
        if (d > max_delta)
          max_delta = d;
        if (d > 6553)
          ++big_delta;
        const bool ship_black =
            !a_buf[p] && !a_buf[p + 1] && !a_buf[p + 2];
        const bool ref_lit = b_buf[p] || b_buf[p + 1] || b_buf[p + 2];
        if (ship_black && ref_lit)
          ++holes;
      }

      if (dump_prefix && o == 0) {
        char p[512];
        std::snprintf(p, sizeof(p), "%s_ship.raw", dump_prefix);
        dump(p, a_buf);
        std::snprintf(p, sizeof(p), "%s_ref.raw", dump_prefix);
        dump(p, b_buf);
      }

      for (size_t i = 0; i < mesh.vertices.size(); ++i)
        mesh.vertices[i] = saved[i];
      ++draws;
    }

    const auto &a = hs_aa::g_audit;
    std::printf("%-58s %5zu %10lld %8lld %8lld %8.2f %7lld\n", entry.name,
                mesh.get_face_counts_size(), a.probes, a.painted, a.missed,
                static_cast<double>(a.missed) / orientations,
                fb_diff / orientations);
    std::printf("    ref probes %lld (%.2fx) | missed alpha max %.4f mean "
                "%.4f | max col gap %d | max fb delta %lld | >10%% px/draw "
                "%lld | HOLES %lld\n",
                ref_probes,
                a.probes ? static_cast<double>(ref_probes) / a.probes : 0.0,
                a.missed_alpha_max,
                a.missed ? a.missed_alpha_sum / a.missed : 0.0, a.max_gap_cols,
                max_delta, big_delta / orientations, holes);
    std::printf("    probe rows:");
    for (int y = 0; y < H; ++y)
      std::printf(" %d:%lld", y, a.probe_rows[y] / orientations);
    std::printf("\n");
    if (a.missed) {
      std::printf("    missed rows (per draw):");
      for (int y = 0; y < H; ++y)
        if (a.missed_rows[y] >= orientations)
          std::printf(" %d:%lld", y, a.missed_rows[y] / orientations);
      std::printf("\n    missed alpha hist (deciles):");
      for (int i = 0; i < 10; ++i)
        std::printf(" %lld", a.alpha_hist[i]);
      std::printf("\n");
    }
    tot_probes += a.probes;
    tot_ref_probes += ref_probes;
    tot_painted += a.painted;
    tot_missed += a.missed;
    tot_diff += fb_diff;
  }

  std::printf("\nTOTAL over %d draws: probes %lld (%.0f/draw) ref probes %lld "
              "(%.0f/draw, %.2fx) painted %lld missed %lld (%.2f/draw) "
              "fb-diff px %lld (%.2f/draw)\n",
              draws, tot_probes, static_cast<double>(tot_probes) / draws,
              tot_ref_probes, static_cast<double>(tot_ref_probes) / draws,
              static_cast<double>(tot_ref_probes) / tot_probes, tot_painted,
              tot_missed, static_cast<double>(tot_missed) / draws, tot_diff,
              static_cast<double>(tot_diff) / draws);
  return 0;
}
