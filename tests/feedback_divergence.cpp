/**
 * @file
 * @brief Divergence trace for perturbations of the feedback colour path.
 * @details Drives Filter::Pixel::Feedback::flush over N frames from a fixed
 * seed with sustained emission, then either dumps every frame or compares two
 * dumps frame by frame. The colour path re-enters as input each frame, so a
 * perturbation that looks negligible in one pass can compound; the trace shape
 * is the answer (flat or decaying = absorbed, rising = compounding).
 *
 * Two builds are needed because the paths being compared are compile-time
 * choices; keeping the switch out of shipping code costs nothing here:
 *
 *   (at the base revision)  feedback_divergence --dump base.bin
 *   (with the change)       feedback_divergence --dump alt.bin
 *                           feedback_divergence --compare base.bin alt.bin
 *
 * Emission and noise advance are pure functions of the frame index, so both
 * runs see identical input and the only variable is the code under test.
 */
// Host tool: the Windows CRT deprecates fopen in favour of fopen_s, which the
// other toolchains do not have.
#define _CRT_SECURE_NO_WARNINGS

#include "tests/test_filter.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

constexpr int W = 288, H = 144;
constexpr int CHANS = W * H * 3;

/** @brief One preset's flush driver: fixed seed, deterministic emission. */
struct Run {
  hs_test::filter_tests::PipeFx fx{W, H};
  Feedback::Style style;
  NoiseParams noise;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  explicit Run(const Feedback::Style &s) : style(s) {
    style.noise = &noise;
    noise.amplitude = style.amplitude;
    noise.frequency = style.frequency;
    noise.speed = style.speed;
    noise.scale = style.scale;
    noise.time = 0.0f;
    noise.sync();
    style.sync_hue();
  }

  // Saturated seed: divergences cluster on primaries driving the gamut-clip
  // path, so a washed-out seed understates the perturbation.
  void seed() {
    {
      Canvas c(fx);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          switch (((x / 24) + (y / 24)) % 6) {
          case 0: c(x, y) = Pixel(65535, 0, 0); break;
          case 1: c(x, y) = Pixel(0, 65535, 0); break;
          case 2: c(x, y) = Pixel(0, 0, 65535); break;
          case 3: c(x, y) = Pixel(65535, 65535, 0); break;
          case 4: c(x, y) = Pixel(65535, 0, 65535); break;
          default: c(x, y) = Pixel(60000, 20000, 5000); break;
          }
        }
    }
    fx.advance_display();
  }

  // Without fresh geometry the seed decays to black and a zero divergence
  // would mean the image died, not that the perturbation was absorbed.
  static void emit(Canvas &c, int frame) {
    const float t = frame * 0.05f;
    const int cx = (int)(W * 0.5f + W * 0.30f * fast_cosf(t));
    const int cy = (int)(H * 0.5f + H * 0.30f * fast_sinf(t * 0.7f));
    for (int dy = -10; dy <= 10; ++dy)
      for (int dx = -10; dx <= 10; ++dx) {
        if (dx * dx + dy * dy > 100) continue;
        const int x = ((cx + dx) % W + W) % W;
        const int y = cy + dy;
        if (y < 0 || y >= H) continue;
        switch ((frame / 7) % 4) {
        case 0: c(x, y) = Pixel(65535, 0, 0); break;
        case 1: c(x, y) = Pixel(0, 65535, 0); break;
        case 2: c(x, y) = Pixel(0, 0, 65535); break;
        default: c(x, y) = Pixel(65535, 65535, 0); break;
        }
      }
  }

  void step(int frame) {
    noise.time += 1.0f / 16.0f;
    auto trail = [](float, float, float) {
      return Color4(Pixel(0, 0, 0), 0.0f);
    };
    {
      Canvas c(fx);
      pipe.flush(c, ScreenTrailFn(trail), 1.0f);
      emit(c, frame);
    }
    fx.advance_display();
  }
};

const Feedback::Style presets[] = {Feedback::Style::Smoke(),
                                   Feedback::Style::Frozen()};
const char *preset_names[] = {"Smoke", "Frozen"};
constexpr int NPRESET = 2;

/** @brief Writes every frame of every preset to `path`. */
int dump(const char *path, int frames) {
  FILE *f = fopen(path, "wb");
  if (!f) {
    printf("cannot open %s for writing\n", path);
    return 1;
  }
  const int magic = 0x46424456; // "FBDV"
  fwrite(&magic, sizeof(magic), 1, f);
  fwrite(&frames, sizeof(frames), 1, f);
  const int np = NPRESET;
  fwrite(&np, sizeof(np), 1, f);

  std::vector<uint16_t> row(CHANS);
  for (int p = 0; p < NPRESET; ++p) {
    Run r(presets[p]);
    r.seed();
    for (int i = 0; i < frames; ++i) {
      r.step(i);
      uint16_t *dst = row.data();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel &px = r.fx.get_pixel(x, y);
          *dst++ = px.r;
          *dst++ = px.g;
          *dst++ = px.b;
        }
      fwrite(row.data(), sizeof(uint16_t), CHANS, f);
    }
    printf("dumped %s: %d frames\n", preset_names[p], frames);
  }
  fclose(f);
  return 0;
}

std::vector<uint16_t> load(const char *path, int &frames, int &np) {
  FILE *f = fopen(path, "rb");
  if (!f) {
    printf("cannot open %s\n", path);
    exit(1);
  }
  int magic = 0;
  if (fread(&magic, sizeof(magic), 1, f) != 1 || magic != 0x46424456) {
    printf("%s: not a divergence dump\n", path);
    exit(1);
  }
  size_t got = fread(&frames, sizeof(frames), 1, f);
  got += fread(&np, sizeof(np), 1, f);
  if (got != 2) {
    printf("%s: truncated header\n", path);
    exit(1);
  }
  std::vector<uint16_t> data((size_t)frames * np * CHANS);
  if (fread(data.data(), sizeof(uint16_t), data.size(), f) != data.size()) {
    printf("%s: truncated body\n", path);
    exit(1);
  }
  fclose(f);
  return data;
}

/** @brief Prints the per-frame divergence trace between two dumps. */
int compare(const char *pa, const char *pb) {
  int fa = 0, na = 0, fb = 0, nb = 0;
  std::vector<uint16_t> a = load(pa, fa, na);
  std::vector<uint16_t> b = load(pb, fb, nb);
  if (fa != fb || na != nb) {
    printf("dump mismatch: %d frames/%d presets vs %d/%d\n", fa, na, fb, nb);
    return 1;
  }

  for (int p = 0; p < na; ++p) {
    printf("\n=== %s: %d frames ===\n",
           p < NPRESET ? preset_names[p] : "preset", fa);
    printf("%6s %12s %10s %10s\n", "frame", "differing", "max|d|", "mean|d|");
    int worst = 0, worst_frame = 0;
    long worst_count = 0, spikes = 0, spike_chans = 0;
    long sp_first = 0, sp_second = 0;
    int sp_first_max = 0, sp_second_max = 0;
    double base_first = 0.0, base_second = 0.0;
    long nf_first = 0, nf_second = 0;

    for (int i = 0; i < fa; ++i) {
      const size_t off = ((size_t)p * fa + i) * CHANS;
      long differing = 0, over1 = 0;
      int maxd = 0;
      double sum = 0.0;
      for (int c = 0; c < CHANS; ++c) {
        const int d = std::abs((int)a[off + c] - (int)b[off + c]);
        if (!d) continue;
        ++differing;
        sum += d;
        if (d > 1) ++over1;
        if (d > maxd) maxd = d;
      }
      if (maxd > worst) { worst = maxd; worst_frame = i + 1; }
      if (differing > worst_count) worst_count = differing;
      const bool second = i >= fa / 2;
      if (maxd > 1) {
        ++spikes;
        spike_chans += over1;
        if (second) { ++sp_second; if (maxd > sp_second_max) sp_second_max = maxd; }
        else { ++sp_first; if (maxd > sp_first_max) sp_first_max = maxd; }
      }
      if (second) { base_second += differing; ++nf_second; }
      else { base_first += differing; ++nf_first; }
      if (i < 5 || (i + 1) % 100 == 0 || i == fa - 1)
        printf("%6d %12ld %10d %10.3f\n", i + 1, differing, maxd,
               differing ? sum / differing : 0.0);
    }
    printf("  peak: max|d| = %d LSB at frame %d; peak differing = %ld of %d "
           "(%.3f%%)\n", worst, worst_frame, worst_count, CHANS,
           100.0 * worst_count / CHANS);
    printf("  differing/frame: first half %.1f, second half %.1f\n",
           nf_first ? base_first / nf_first : 0.0,
           nf_second ? base_second / nf_second : 0.0);
    printf("  frames with max|d|>1: %ld of %d; mean %.1f channels each\n",
           spikes, fa, spikes ? (double)spike_chans / spikes : 0.0);
    printf("  trend: first half %ld (max %d LSB), second half %ld (max %d "
           "LSB)\n", sp_first, sp_first_max, sp_second, sp_second_max);
  }
  return 0;
}

void usage() {
  printf("usage: feedback_divergence --dump <file> [frames]\n"
         "       feedback_divergence --compare <base> <alt>\n");
}

} // namespace

int main(int argc, char **argv) {
  if (argc >= 3 && std::strcmp(argv[1], "--dump") == 0)
    return dump(argv[2], argc > 3 ? atoi(argv[3]) : 1000);
  if (argc >= 4 && std::strcmp(argv[1], "--compare") == 0)
    return compare(argv[2], argv[3]);
  usage();
  return 1;
}
