#pragma once

#include <tuple>
#include <utility>
#include <type_traits>
#include <cmath>
#include <vector>
#include "geometry.h" 
#include "color.h" 
#include "static_circular_buffer.h"

// Filter Traits
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
};
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
};

struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
};

template <int W, typename... Filters>
struct Pipeline;

/**
 * @brief Plots a color to the canvas, ensuring the y-coordinate is within bounds (0 to H-1).
 * * This function serves as a safety check for vertical bounds, plotting only to the
 * physical LED segments (y < H).
 * * @param canvas The canvas object (buffer).
 * @param x The horizontal coordinate (LED index within the buffer, modulo W is often handled upstream).
 * @param y The vertical coordinate (row index).
 * @param c The color (Pixel) to plot.
 */
void plot_virtual(Canvas& canvas, int x, int y, const Pixel& c) {
  if (y >= 0 && y < H) {
    canvas(XY(x, y)) = c;
  }
}

// Base Case: Canvas Sink
template <int W>
struct Pipeline<W> {
  static constexpr bool is_2d = true;

  // 2D Sink
  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha) {
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    xi = wrap(xi, W);
    auto p = blend_alpha(alpha)(cv(xi, yi), c);
    plot_virtual(cv, xi, yi, p);
  }

  // 3D Sink
  void plot(Canvas& cv, const Vector& v, const Pixel& c, float age, float alpha) {
    auto p = vector_to_pixel<W>(v);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  // History
  void trail(Canvas&, TrailFn auto, float) {}
};

// Recursive Case
template <int W, typename Head, typename... Tail>
struct Pipeline<W, Head, Tail...> : public Head {
  using Next = Pipeline<W, Tail...>;
  Next next;

  Pipeline(Head h, Tail... t) : Head(std::move(h)), next(std::move(t)...) {}
  Pipeline() = default;

  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha) {
    if constexpr (Head::is_2d) {
      // 2D -> 2D
      Head::plot(x, y, c, age, alpha,
        [&](float x, float y, const Pixel& c, float age, float alpha) {
          next.plot(cv, x, y, c, age, alpha);
        });
    }
    else {
      // 2D -> 3D Mismatch
      Vector v = pixel_to_vector<W>(x, y);
      plot(cv, v, c, age, alpha);
    }
  }

  void plot(Canvas& cv, const Vector& v, const Pixel& c, float age, float alpha) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha,
        [&](const Vector& v, const Pixel& c, float age, float alpha) {
          if constexpr (Next::is_2d) {
            // 3D -> 2D
            auto p = vector_to_pixel<W>(v);
            next.plot(cv, p.x, p.y, c, age, alpha);
          }
          else {
            // 3D -> 3D
            next.plot(cv, v, c, age, alpha);
          }
        });
    }
    else {
      // 3D -> 2D Mismatch
      auto p = vector_to_pixel<W>(v);
      plot(cv, p.x, p.y, c, age, alpha);
    }
  }

  void trail(Canvas& cv, TrailFn auto trailFn, float alpha) {
    if constexpr (Head::has_history) {
      Head::trail(trailFn, alpha,
        [&](float x, float y, const Pixel& c, float age, float alpha) {
          next.plot(cv, x, y, c, age, alpha);
        });
    }
    next.trail(cv, trailFn, alpha);
  }
};

///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Quintic Smootherstep Kernel function (Perlin's 5th order polynomial).
 * @details Guarantees zero velocity and acceleration at t=0 and t=1 for smooth anti-aliasing.
 * @param t The fractional coordinate (0.0 to 1.0).
 * @return The smoothed weight (0.0 to 1.0).
 */
inline float quintic_kernel(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  // W(t) = 6*t^5 - 15*t^4 + 10*t^3
  return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

/**
 * @brief A filter that applies the current global orientation to the pixel position.
 * @details Necessary for effects that rotate the entire scene. Converts the 2D pixel to a 3D vector,
 * rotates it by the Orientation quaternion, and re-projects it to 2D.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterOrient : public Is3D {
public:
  /**
   * @brief Initializes the filter with a reference to the global orientation state.
   * @param orientation Reference to the Orientation object.
   */
  FilterOrient(Orientation& orientation) :
    orientation(orientation)
  {
  }

  /**
   * @brief Rotates the pixel's 3D position based on the current orientation.
   */
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    pass(orientation.orient(v), color, age, alpha);
  }

private:

  Orientation& orientation; /**< Reference to the global Orientation state. */
};

/**
 * @brief Applies different orientations to points in n latitude bands defined by axis.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterOrientSlice : public Is3D {
public:
  /**
   * @brief Initializes the filter with a list of orientations and an axis.
   * @param orientations Vector of pointers to Orientation objects.
   * @param axis The axis defining the poles for slicing.
   */
  FilterOrientSlice(const std::vector<Orientation*>& orientations, const Vector& axis) :
    orientations(orientations),
    axis(axis)
  {
  }

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    float d = std::max(-1.0f, std::min(1.0f, dot(v, axis)));
    float t = 1.0f - std::acos(d) / PI_F;

    int idx = static_cast<int>(std::floor(t * orientations.size()));
    if (idx >= static_cast<int>(orientations.size())) idx = orientations.size() - 1;
    if (idx < 0) idx = 0;

    pass(orientations[idx]->orient(v), color, age, alpha);
  }

private:
  std::vector<Orientation*> orientations;
  Vector axis;
};

/**
 * @brief A filter that replicates the plotted pixel horizontally multiple times.
 * @details Used for visual effects like simulating multiple light sources or widening a sparse plot.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterReplicate : public Is3D {
public:

  /**
   * @brief Initializes the filter with the number of times to replicate.
   * @param count The number of copies to create horizontally. Clamped between 1 and W.
   */
  FilterReplicate(int count) :
    count(std::clamp(count, 1, W)),
    step(make_rotation(Y_AXIS, 2 * PI_F / count))
  {
  }

  /**
   * @brief Plots the fragment at its original position and at evenly spaced offsets.
   */
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    Vector r = v;
    pass(r, color, age, alpha);
    for (int i = 1; i < count; i++) {
      r = rotate(r, step);
      pass(r, color, age, alpha);
    }
  }

private:

  int count; /**< The number of times to replicate the pixel. */
  Quaternion step;
};

/**
 * @brief Applies an alpha falloff based on distance from an origin point on the sphere.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterHole : public Is3D {
public:
  /**
   * @brief Initializes the filter with an origin and radius.
   * @param origin The center point of the hole (normalized).
   * @param radius The radius (in radians) at which fading starts.
   */
  FilterHole(const Vector& origin, float radius) :
    origin(origin),
    radius(radius)
  {
  }

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) {
      pass(v, color, age, alpha);
    }
    else {
      float t = d / radius;
      float factor = quintic_kernel(t);
      Pixel c = color;
      c.r = static_cast<uint8_t>(c.r * factor);
      c.g = static_cast<uint8_t>(c.g * factor);
      c.b = static_cast<uint8_t>(c.b * factor);
      pass(v, c, age, alpha);
    }
  }

private:
  Vector origin;
  float radius;
};

/**
 * @brief Applies an alpha falloff based on distance from an origin point on the sphere.
 * Operates on a reference to the origin origin vector.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterHoleRef : public Is3D {
public:
  /**
   * @brief Initializes the filter with an origin and radius.
   * @param origin The center point of the hole (normalized).
   * @param radius The radius (in radians) at which fading starts.
   */
  FilterHoleRef(const Vector& origin, float radius) : origin(origin), radius(radius) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) {
      pass(v, color, age, alpha);
    }
    else {
      float t = d / radius;
      // Quintic kernel: 6t^5 - 15t^4 + 10t^3
      float factor = t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
      Pixel c = color;
      c.r = static_cast<uint8_t>(c.r * factor);
      c.g = static_cast<uint8_t>(c.g * factor);
      c.b = static_cast<uint8_t>(c.b * factor);
      pass(v, c, age, alpha);
    }
  }
private:
  const Vector& origin;
  float radius;
};


/**
 * @brief Applies a Mobius Transformation to 3D vectors.
 * @details Projects sphere -> complex plane -> transform -> sphere.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterMobius : public Is3D {
public:
  /**
   * @brief Initializes the filter with a reference to the Mobius parameters.
   * @param params Reference to the MobiusParams object.
   */
  FilterMobius(MobiusParams& params) :
    params(params)
  {
  }

  /**
   * @brief Applies the Mobius transformation to the vector.
   */
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    // 1. Stereo Projection -> 2. Mobius Transform -> 3. Inverse Stereo
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha);
  }

private:
  MobiusParams& params;
};

///////////////////////////////////////////////////////////////////////////////
// 2D Filters
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief A filter implementing anti-aliasing via Quintic Interpolation (Smootherstep).
 * @details Distributes the color/alpha of a sub-pixel position across its four nearest integer pixels
 * using a zero-acceleration weighting curve.
 * @tparam W The width of the effect.
 */
template <int W>
class FilterAntiAlias : public Is2D {
public:
  FilterAntiAlias() {}

  void plot(float x, float y, const Pixel& c, float age, float alpha, auto pass) {
    float x_i = 0;
    float x_m = std::modf(x, &x_i);
    float y_i = 0;
    float y_m = std::modf(y, &y_i);
    float xs = quintic_kernel(x_m);
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);  // Top-Left weight (xi, yi)
    float v10 = xs * (1 - ys);        // Top-Right weight (xi+1, yi)
    float v01 = (1 - xs) * ys;        // Bottom-Left weight (xi, yi+1)
    float v11 = xs * ys;              // Bottom-Right weight (xi+1, yi+1)

    if (v00 > 0.0001) {
      pass(x_i, y_i, c, age, alpha * v00);
    }
    if (v10 > 0.0001) {
      pass(wrap((x_i + 1), W), y_i, c, age, alpha * v10);
    }
    if (y_i < H - 1) {
      if (v01 > 0.0001) {
        pass(x_i, y_i + 1, c, age, alpha * v01);
      }
      if (v11 > 0.0001) {
        pass(wrap((x_i + 1), W), y_i + 1, c, age, alpha * v11);
      }
    }
  }
};

/**
 * @brief A filter that shifts the RGB channels slightly to create a chromatic aberration effect.
 * @tparam W The width of the effect.
 */
template<int W>
class FilterChromaticShift : public Is2D {
public:

  FilterChromaticShift()
  {
  }

  /**
   * @brief Plots the original color, then plots R, G, and B components at slightly offset positions.
   */
  void plot(float x, float y, const Pixel& color, float age, float alpha, auto pass) {
    CRGB r(color.r, 0, 0);
    CRGB g(0, color.g, 0);
    CRGB b(0, 0, color.b);
    pass(x, y, color, age, alpha);
    pass(wrap(x + 1.0f, W), y, r, age, alpha);
    pass(wrap(x + 2.0f, W), y, g, age, alpha);
    pass(wrap(x + 3.0f, W), y, b, age, alpha);
  }
};

/**
 * @brief A filter that manages time-to-live (TTL) for fragments, creating trails.
 * * @details Fragments are not passed immediately but tracked internally. On decay,
 * they are plotted as a trail using a color function based on their remaining TTL.
 * * @tparam W The width of the effect.
 * @tparam MAX_PIXELS The maximum number of decaying pixels to track.
 */
template <int W, int MAX_PIXELS>
class FilterDecay : public Is2DWithHistory {
public:

  /**
   * @brief Initializes the filter with a fragment lifetime.
   * @param lifetime The number of frames a trail fragment persists.
   */
  FilterDecay(int lifetime) :
    lifetime(lifetime),
    num_pixels(0)
  {
  }

  /**
   * @brief Plots a new fragment, optionally tracking it for decay.
   * @details If age >= 0, the fragment is tracked for trails. If age <= 0,
   * it is plotted immediately as a new head.
   */
  void plot(float x, float y, const Pixel& color, float age, float alpha, auto pass) {
    if (age >= 0) {
      if (num_pixels < MAX_PIXELS) {
        ttls[num_pixels++] = { static_cast<float>(x), static_cast<float>(y), static_cast<float>(lifetime - age) };
      }
      else {
        Serial.println("FilterDecay full!");
      }
    }
    if (age <= 0) {
      pass(x, y, color, age, alpha);
    }
  }

  /**
  * @brief Re-plots all existing trail fragments with reduced brightness/color based on age.
  * @param canvas The canvas to plot onto.
  * @param trailFn A function to calculate the color based on normalized age (0.0 to 1.0).
  * @param alpha The base opacity for the trails.
  */
  void trail(TrailFn auto trailFn, float alpha, auto pass) {
    for (int i = 0; i < num_pixels; ++i) {
      auto color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      pass(ttls[i].x, ttls[i].y, color, lifetime - ttls[i].ttl, alpha);
    }
    decay();
  }

  /**
   * @brief Decrements the TTL of all tracked fragments and removes expired ones.
   */
  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls[i].ttl < TOLERANCE) {
        num_pixels--;
        if (i < num_pixels) {
          ttls[i] = std::move(ttls[num_pixels]);
          i--;
        }
      }
    }
  }

private:

  /**
   * @brief Structure holding the state of a single decaying pixel fragment.
   */
  struct DecayPixel {
    float x; /**< X coordinate. */
    float y; /**< Y coordinate. */
    float ttl; /**< Time to Live (remaining frames). */
  };

  int lifetime; /**< The initial lifetime for a fragment. */
  std::array<DecayPixel, MAX_PIXELS> ttls; /**< Array of tracked decay pixels. */
  int num_pixels; /**< The current number of active decay pixels. */
};