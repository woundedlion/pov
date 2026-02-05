/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../filter.h"
#include "../palettes.h"
#include <vector>
#include <cmath>

template <int W>
class MetaballEffect : public Effect {
public:
    struct Ball {
        Vector p;
        Vector v;
        float r;
    };

    MetaballEffect() : Effect(W), palette(Palettes::richSunset) {
        init_balls();
    }

    bool show_bg() const override { return false; } // We overwrite all pixels

    void draw_frame() override {
        Canvas canvas(*this); // Not used if we write directly to buffer? 
        // Effect usually exposes buffer access via operator() or plot functions.
        // Canvas wraps 'Effect&'.
        
        // 1. Physics
        for (auto& b : balls) {
            Vector force = b.p * (-gravity);
            b.v += force;
            b.p += b.v;
            // No bounds check? JS version doesn't seem to enforce sphere bounds strictly, 
            // relying on gravity to pull them back.
        }
        
        // 2. Render
        // We iterate coordinates. 
        // Optimization: Use LUT for unique vectors if W is small.
        // Or just compute.
        
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                Vector v = pixel_to_vector<W>(x, y);
                
                float sum = 0.0f;
                for (const auto& b : balls) {
                    float dist_sq = distance_squared(v, b.p);
                    // Avoid div by zero
                    if (dist_sq < 0.0001f) dist_sq = 0.0001f;
                    sum += (b.r * b.r) / dist_sq;
                }
                
                float t = sum / max_influence;
                if (t > 1.0f) t = 1.0f;
                
                // Color
                Color4 c = palette.get(t);
                
                // Write
                // canvas(x, y) = c;
                // Effect::operator() usually takes XY(x,y)
                // canvas(x, y) calls Effect::set/blend.
                canvas(x, y) = c; 
            }
        }
    }
    
    // Params
    float max_influence = 10.0f;
    float gravity = 0.005f;
    int num_balls = 16;
    float radius_scale = 1.0f;
    float velocity_scale = 1.0f;

private:
    std::vector<Ball> balls;
    const Palette& palette;

    void init_balls() {
        balls.clear();
        balls.reserve(num_balls);
        for (int i = 0; i < num_balls; ++i) {
            Vector p = random_vector() * 0.5f; // Start inside
            float r = (random(0.5f, 0.8f)) * radius_scale;
            Vector v = random_vector() * (0.05f * velocity_scale);
            balls.push_back({p, v, r});
        }
    }
    
    // Random helper
    float random(float min, float max) {
        return min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min)));
    }
    
    Vector random_vector() {
        // Uniform random unit vector
        // Use rejection or Marsaglia
        // Simple approximation for init:
        float x = random(-1, 1);
        float y = random(-1, 1);
        float z = random(-1, 1);
        Vector v(x, y, z);
        return v.normalize();
    }
    
    float distance_squared(const Vector& a, const Vector& b) {
        float dx = a.i - b.i;
        float dy = a.j - b.j;
        float dz = a.k - b.k;
        return dx*dx + dy*dy + dz*dz;
    }
};
