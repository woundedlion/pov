/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/vendor/FastNoiseLite.h"
#include <cstdio>

int main() {
  FastNoiseLite noise(1337);
  noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.SetFrequency(0.125f);
  std::puts("// OpenSimplex2 3D: 5x5x5, spacing 1.5, origin -3");
  std::puts("static constexpr float GOLDEN[] = {");
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k)
        std::printf(
            "    %.9ef,\n",
            noise.GetNoise(i * 1.5f - 3.0f, j * 1.5f - 3.0f, k * 1.5f - 3.0f));
  std::puts("};\n// OpenSimplex2 2D: 8x8, spacing 0.75, origin -3");
  std::puts("static constexpr float GOLDEN[] = {");
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j)
      std::printf("    %.9ef,\n",
                  noise.GetNoise(i * 0.75f - 3.0f, j * 0.75f - 3.0f));
  std::puts("};");
  if (std::fflush(stdout) != 0 || std::ferror(stdout)) {
    std::perror("noise_golden_gen: stdout");
    return 1;
  }
  return 0;
}
