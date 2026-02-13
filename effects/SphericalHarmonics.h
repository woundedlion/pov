/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
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

template <int W, int H>
class SphericalHarmonics : public Effect {
public:
    struct HarmonicBlob {
        int l, m;
        float amplitude;
        Quaternion orientation;
        
        HarmonicBlob(int l, int m, float amp, Quaternion q) : l(l), m(m), amplitude(amp), orientation(q) {}
        
        template <int H_>
        SDF::Bounds get_vertical_bounds() const {
             return { 0, H_ - 1 }; // Full Scan fallback
        }
        
        template<int H_, typename OutputIt>
        bool get_horizontal_intervals(int y, OutputIt out) const {
             return false; // Full Scan fallback
        }
        
        SDF::DistanceResult distance(const Vector& p) const {
            return distance<true>(p);
        }

        template <bool ComputeUVs = true>
        SDF::DistanceResult distance(const Vector& p) const {
             Vector local = rotate(p, orientation.conjugate()); 
             
             float theta = atan2f(local.k, local.i); 
             if (theta < 0) theta += 2 * PI_F;
             float phi = acosf(std::clamp(local.j, -1.0f, 1.0f));
             float val = SHMath::sphericalHarmonic(l, m, theta, phi);
             
             return { -1.0f, 0.0f, val }; // t is unused
        }
    };

    SphericalHarmonics() : Effect(W, H), filters() {
        // Spin
        Vector axis = Vector(0.5f, 1.0f, 0.2f).normalize();
        timeline.add(0, Animation::Rotation<W>(orientation, axis, 2 * PI_F * 100, 10000, ease_mid, true));
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        int idx = (int)mode;
        int l = (int)sqrtf((float)idx);
        int m = idx - l * l - l;
        
        HarmonicBlob blob(l, m, amplitude, orientation.get());
        
        auto shader = [&](const Vector& p, const Fragment& frag) -> Fragment {
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
             
             Fragment out = frag;
             out.color = base;
             return out;
        };
        
        Scan::rasterize<W, H>(filters, canvas, blob, shader);
    }
    
    // Params
    float mode = 6.0f;
    float amplitude = 3.2f;

private:
    Orientation<W> orientation;
    Timeline<W> timeline;
    Pipeline<W, H> filters; // No filters needed? JS has pipeline.
};
