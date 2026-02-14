/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include <cmath>

template <int W, int H>
class MetaballEffect : public Effect {
public:
    struct Ball {
        Vector p;
        Vector v;
        float r;
    };

    MetaballEffect() : Effect(W, H), palette(Palettes::richSunset) {
        init_balls();
    }

    bool show_bg() const override { return false; } // We overwrite all pixels

    void draw_frame() override {
        Canvas canvas(*this); 

        // 1. Physics
        for (auto& b : balls) {
            Vector force = b.p * (-gravity);
            b.v += force;
            b.p += b.v;
        }
        
        // 2. Render
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                Vector v = pixel_to_vector<W, H>(x, y);
                
                float sum = 0.0f;
                for (const auto& b : balls) {
                    float dist_sq = distance_squared(v, b.p);
                    // Avoid div by zero
                    if (dist_sq < 0.0001f) dist_sq = 0.0001f;
                    sum += (b.r * b.r) / dist_sq;
                }
                
                float t = sum / max_influence;
                if (t > 1.0f) t = 1.0f;
                
                Color4 c = palette.get(t);
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
            float r = hs::rand_f(0.5f, 0.8f) * radius_scale;
            Vector v = random_vector() * (0.05f * velocity_scale);
            balls.push_back({p, v, r});
        }
    }
            
    float distance_squared(const Vector& a, const Vector& b) {
        float dx = a.i - b.i;
        float dy = a.j - b.j;
        float dz = a.k - b.k;
        return dx*dx + dy*dy + dz*dz;
    }
};
