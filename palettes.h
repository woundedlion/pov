/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "color.h"

namespace Palettes {

    // Procedural Palettes
    static constexpr ProceduralPalette darkRainbow({ 0.367f, 0.367f, 0.367f }, { 0.500f, 0.500f, 0.500f }, { 1.000f, 1.000f, 1.000f }, { 0.000f, 0.330f, 0.670f });
    static constexpr ProceduralPalette bloodStream({ 0.169f, 0.169f, 0.169f }, { 0.313f, 0.313f, 0.313f }, { 0.231f, 0.231f, 0.231f }, { 0.036f, 0.366f, 0.706f });
    static constexpr ProceduralPalette vintageSunset({ 0.256f, 0.256f, 0.256f }, { 0.500f, 0.080f, 0.500f }, { 0.277f, 0.277f, 0.277f }, { 0.000f, 0.330f, 0.670f });
    static constexpr ProceduralPalette richSunset({ 0.309f, 0.500f, 0.500f }, { 1.000f, 1.000f, 0.500f }, { 0.149f, 0.148f, 0.149f }, { 0.132f, 0.222f, 0.521f });
    static constexpr ProceduralPalette undersea({ 0.000f, 0.000f, 0.000f }, { 0.500f, 0.276f, 0.423f }, { 0.296f, 0.296f, 0.296f }, { 0.374f, 0.941f, 0.000f });
    static constexpr ProceduralPalette lateSunset({ 0.337f, 0.500f, 0.096f }, { 0.500f, 1.000f, 0.176f }, { 0.261f, 0.261f, 0.261f }, { 0.153f, 0.483f, 0.773f });
    static constexpr ProceduralPalette mangoPeel({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.080f, 0.500f }, { 0.431f, 0.431f, 0.431f }, { 0.566f, 0.896f, 0.236f });
    static constexpr ProceduralPalette iceMelt({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.500f, 0.500f }, { 0.083f, 0.147f, 0.082f }, { 0.579f, 0.353f, 0.244f });
    static constexpr ProceduralPalette lemonLime({ 0.455f, 0.455f, 0.455f }, { 0.571f, 0.151f, 0.571f }, { 0.320f, 0.320f, 0.320f }, { 0.087f, 0.979f, 0.319f });
    static constexpr ProceduralPalette algae({ 0.210f, 0.210f, 0.210f }, { 0.500f, 1.000f, 0.021f }, { 0.086f, 0.086f, 0.075f }, { 0.419f, 0.213f, 0.436f });
    static constexpr ProceduralPalette embers({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.500f, 0.500f }, { 0.265f, 0.285f, 0.198f }, { 0.577f, 0.440f, 0.358f });
    static constexpr ProceduralPalette fireGlow({ 0.000f, 0.000f, 0.000f }, { 0.560f, 0.560f, 0.560f }, { 0.216f, 0.346f, 0.174f }, { 0.756f, 0.542f, 0.279f });
    static constexpr ProceduralPalette darkPrimary({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.610f, 0.500f }, { 0.746f, 0.347f, 0.000f }, { 0.187f, 0.417f, 0.670f });
    static constexpr ProceduralPalette mauveFade({ 0.583f, 0.000f, 0.583f }, { 1.000f, 0.000f, 1.000f }, { 0.191f, 0.348f, 0.191f }, { 0.175f, 0.045f, 0.150f });
    static constexpr ProceduralPalette lavenderLake({ 0.473f, 0.473f, 0.473f }, { 0.500f, 0.500f, 0.500f }, { 0.364f, 0.124f, 0.528f }, { 0.142f, 0.378f, 0.876f });
    static constexpr ProceduralPalette desertRose({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.270f, 0.442f }, { 0.303f, 1.012f, 0.585f }, { 0.985f, 0.720f, 0.212f });
    static constexpr ProceduralPalette bruisedMoss({ 0.500f, 0.500f, 0.500f }, { 0.500f, 0.500f, 0.500f }, { 0.142f, 0.252f, 0.000f }, { 0.492f, 0.200f, 0.670f });
    static constexpr ProceduralPalette bruisedBanana({ 0.620f, 0.620f, 0.620f }, { 0.742f, 0.742f, 0.742f }, { 0.162f, 0.286f, 0.012f }, { 0.235f, 0.205f, 0.688f });
    static constexpr ProceduralPalette brightSunrise({ 0.620f, 0.620f, 0.620f }, { 0.742f, 0.742f, 0.742f }, { 0.162f, 0.286f, 0.012f }, { 0.090f, 0.205f, 0.688f });

    // Gradient Palettes
    static const Gradient rainbow({
      {0.0f / 16, CPixel(0xFF0000)},
      {1.0f / 16, CPixel(0xD52A00)},
      {2.0f / 16, CPixel(0xAB5500)},
      {3.0f / 16, CPixel(0xAB7F00)},
      {4.0f / 16, CPixel(0xABAB00)},
      {5.0f / 16, CPixel(0x56D500)},
      {6.0f / 16, CPixel(0x00FF00)},
      {7.0f / 16, CPixel(0x00D52A)},
      {8.0f / 16, CPixel(0x00AB55)},
      {9.0f / 16, CPixel(0x0056AA)},
      {10.0f / 16, CPixel(0x0000FF)},
      {11.0f / 16, CPixel(0x2A00D5)},
      {12.0f / 16, CPixel(0x5500AB)},
      {13.0f / 16, CPixel(0x7F0081)},
      {14.0f / 16, CPixel(0xAB0055)},
      {15.0f / 16, CPixel(0xD5002B)},
      {16.0f / 16, CPixel(0xD5002B)}
      });

    static const Gradient rainbowStripes({
      {0.0f / 16, CPixel(0xFF0000)},
      {1.0f / 16, CPixel(0x000000)},
      {2.0f / 16, CPixel(0xAB5500)},
      {3.0f / 16, CPixel(0x000000)},
      {4.0f / 16, CPixel(0xABAB00)},
      {5.0f / 16, CPixel(0x000000)},
      {6.0f / 16, CPixel(0x00FF00)},
      {7.0f / 16, CPixel(0x000000)},
      {8.0f / 16, CPixel(0x00AB55)},
      {9.0f / 16, CPixel(0x000000)},
      {10.0f / 16, CPixel(0x0000FF)},
      {11.0f / 16, CPixel(0x000000)},
      {12.0f / 16, CPixel(0x5500AB)},
      {13.0f / 16, CPixel(0x000000)},
      {14.0f / 16, CPixel(0xAB0055)},
      {15.0f / 16, CPixel(0x000000)},
      {16.0f / 16, CPixel(0xFF0000)}
      });

    static const Gradient rainbowThinStripes({
      {0.0f, CPixel(0xFF0000)},
      {1.0f / 32, CPixel(0x000000)},
      {3.0f / 32, CPixel(0x000000)},
      {4.0f / 32, CPixel(0xAB5500)},
      {5.0f / 32, CPixel(0x000000)},
      {7.0f / 32, CPixel(0x000000)},
      {8.0f / 32, CPixel(0xABAB00)},
      {9.0f / 32, CPixel(0x000000)},
      {11.0f / 32, CPixel(0x000000)},
      {12.0f / 32, CPixel(0x00FF00)},
      {13.0f / 32, CPixel(0x000000)},
      {15.0f / 32, CPixel(0x000000)},
      {16.0f / 32, CPixel(0x00AB55)},
      {17.0f / 32, CPixel(0x000000)},
      {19.0f / 32, CPixel(0x000000)},
      {20.0f / 32, CPixel(0x0000FF)},
      {21.0f / 32, CPixel(0x000000)},
      {23.0f / 32, CPixel(0x000000)},
      {24.0f / 32, CPixel(0x5500AB)},
      {25.0f / 32, CPixel(0x000000)},
      {27.0f / 32, CPixel(0x000000)},
      {28.0f / 32, CPixel(0xAB0055)},
      {29.0f / 32, CPixel(0x000000)},
      {32.0f / 32, CPixel(0x000000)}
      });

    static const Gradient emeraldForest({
      {0.0f, CPixel(0x004E64)},
      {0.2f, CPixel(0x0B6E4F)},
      {0.4f, CPixel(0x08A045)},
      {0.6f, CPixel(0x6BBF59)},
      {0.8f, CPixel(0x138086)},
      {1.0f, CPixel(0x000000)}
      });
}
