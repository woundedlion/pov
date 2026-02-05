/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../scan.h"
#include "../animation.h"
#include "../filter.h"
#include "../palettes.h"
#include <vector>
#include <cmath>

// Math helpers
namespace SHMath {
    inline float factorial(int n) {
        if (n <= 1) return 1.0f;
        float result = 1.0f;
        for (int i = 2; i <= n; i++) result *= i;
        return result;
    }

    inline float associatedLegendre(int l, int m, float x) {
        float pmm = 1.0f;
        if (m > 0) {
            float somx2 = sqrtf(std::max(0.0f, (1.0f - x) * (1.0f + x)));
            float fact = 1.0f;
            for (int i = 1; i <= m; i++) {
                pmm *= -fact * somx2;
                fact += 2.0f;
            }
        }
        if (l == m) return pmm;

        float pmmp1 = x * (2.0f * m + 1.0f) * pmm;
        if (l == m + 1) return pmmp1;

        float pll = 0;
        for (int ll = m + 2; ll <= l; ll++) {
            pll = ((2.0f * ll - 1.0f) * x * pmmp1 - (ll + m - 1.0f) * pmm) / (ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }
        return pll;
    }

    inline float sphericalHarmonic(int l, int m, float theta, float phi) {
        int absM = std::abs(m);
        float N = sqrtf(((2.0f * l + 1.0f) / (4.0f * PI_F)) * (factorial(l - absM) / factorial(l + absM)));
        float P = associatedLegendre(l, absM, cosf(phi));

        if (m > 0) {
            return sqrtf(2.0f) * N * P * cosf(m * theta);
        } else if (m < 0) {
            return sqrtf(2.0f) * N * P * sinf(absM * theta);
        } else {
            return N * P;
        }
    }
}

template <int W>
class SphericalHarmonics : public Effect {
public:
    struct HarmonicBlob {
        int l, m;
        float amplitude;
        Quaternion orientation;
        
        HarmonicBlob(int l, int m, float amp, Quaternion q) : l(l), m(m), amplitude(amp), orientation(q) {}
        
        SDF::Bounds get_vertical_bounds() const {
             return { 0, H - 1 }; // Full Scan fallback
        }
        
        template<typename OutputIt>
        bool get_horizontal_intervals(int y, OutputIt out) const {
             return false; // Full Scan fallback
        }
        
        SDF::DistanceResult distance(const Vector& p) const {
             // Rotate p into local space
             // orientation is object->world. Inverse is world->object.
             Vector local = rotate(p, orientation.conjugate()); 
             
             // Spherical coords
             float r = 1.0f; // Unit sphere sampling
             float theta = atan2f(local.k, local.i); 
             if (theta < 0) theta += 2 * PI_F;
             
             // Phi is angle from Y axis? 
             // 3dmath/Scan usually uses Y as up. Phi 0 is +Y?
             // check y_to_phi: y=0 -> phi=0 (+Y). y=H -> phi=PI (-Y).
             float phi = acosf(std::clamp(local.j, -1.0f, 1.0f));
             float val = SHMath::sphericalHarmonic(l, m, theta, phi);
             

             
             // The shape surface is defined where radius R = 1 + amplitude * val?
             // Or R = amplitude * val?
             // JS: "out.rawDist contains the signed harmonic value"
             // JS Fragment Shader uses `val` to color. 
             // But what defines the geometry?
             // Usually Harmonic Blob is r = |Ylm|.
             // SDF distance is roughly (r - |Ylm|).
             // But Scan::rasterize expects signed distance to SURFACE.
             // If we just want to visualize the value on the sphere, we can use a sphere SDF and pass value in raw_dist.
             // But JS draws a 3D blob shape. 
             // "distance" should be (1.0 - (radius of blob at this angle)).
             // Blob radius R = abs(val) * amplitude.
             // Dist = 1.0 (screen radius) - R? No.
             // We are raymarching or rasterizing? 
             // Scan::rasterize is for Spherical intersection.
             // If we assume a unit sphere canvas, we just return negative dist if we want to hit it?
             // Actually, if we return `dist <= 0`, it draws.
             // If we want to color the WHOLE sphere with the harmonic pattern, we say dist = -1.0 (always inside).
             
             // JS Logic: new SDF.HarmonicBlob(...) passed to Scan.rasterize.
             // Scan.rasterize iterates pixels. If dist < 0, it draws fragment.
             // If HarmonicBlob returns distance to the harmonic surface...
             // But Holosphere is a SPHERICAL DISPLAY. We usually draw on the surface (r=1).
             // Visualization of harmonics often R = |Ylm|. This makes a 3D blob.
             // IF we project this blob onto the LED sphere (r=1), we are just seeing it from 0 or infinity?
             // Users usually want to see the SHAPE.
             // If we are mapping R to brightness, that's one thing.
             // If we are mapping it as a 3D object inside the sphere, we need ray intersection.
             // Since Scan::rasterize is simple, maybe it just maps values to colors on the sphere?
             // "Digital Twin Fragment Shader" uses `val` and `absVal` to color.
             // It does NOT seem to do raymarching.
             // It likely draws on the sphere surface based on the value at that angle.
             // So we return `dist = -1.0f` (Always intersect) and pass `val` as `raw_dist`.
             
             return { -1.0f, 0.0f, val }; // t is unused
        }
    };

    SphericalHarmonics() : Effect(W), filters() {
        // Spin
        Vector axis = Vector(0.5f, 1.0f, 0.2f).normalize();
        timeline.add(0, Rotation<W>(orientation, axis, 2 * PI_F * 100, 10000, ease_mid, true));
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        int idx = (int)mode;
        int l = (int)sqrtf((float)idx);
        int m = idx - l * l - l;
        
        HarmonicBlob blob(l, m, amplitude, orientation.get());
        
        auto shader = [&](const Vector& p, const Fragment& frag) {
             // float val = frag.v2; // Unused
             // HarmonicBlob::distance returns { -1.0f, 0.0f, val }.
             // Scan::rasterize maps: v0=t, v1=raw_dist.
             // Wait, `Scan::rasterize` (generic) populates `f` from `dist_res`.
             // `f.v0 = res.t`, `f.v1 = res.raw_dist`.
             // So `val` is in `frag.v1`.
             
             float abs_val = std::abs(frag.v1);
             
             Color4 base;
             if (frag.v1 >= 0) {
                 base = Palettes::richSunset.get(std::min(1.0f, abs_val * amplitude));
             } else {
                 Color4 p = Palettes::richSunset.get(std::min(1.0f, abs_val * amplitude));
                 base = Color4(Pixel(p.color.b, static_cast<uint8_t>(p.color.g * 0.8f), p.color.r), p.alpha);
             }
             
             // Ambient Occlusion
             float shadow = std::clamp((abs_val * amplitude - 0.0f) / 0.4f, 0.0f, 1.0f); // approx smoothstep
             float occlusion = 0.15f + 0.85f * shadow;
             base.color = base.color * occlusion;
             
             return base;
        };
        
        Scan::rasterize<W>(filters, canvas, blob, shader);
    }
    
    // Params
    float mode = 6.0f;
    float amplitude = 3.2f;

private:
    Orientation orientation;
    Timeline timeline;
    Pipeline<W> filters; // No filters needed? JS has pipeline.
};
