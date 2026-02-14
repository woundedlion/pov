/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include <cmath>


template <int W, int H>
class HopfFibration : public Effect {
public:
    static constexpr int MAX_TRAILS = 50000; 

    HopfFibration() : Effect(W, H), 
        filters(Filter::World::Trails<W, MAX_TRAILS>(40), Filter::World::Orient<W>(orientation), Filter::Screen::AntiAlias<W, H>())
    {
        persist_pixels = false;
        init_fibers();
        timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 600 , ease_mid, true)); 
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        
        timeline.step(canvas);
        
        // Update Params
        flow_offset += 0.02f * flow_speed * 0.2f;
        tumble_angle_x += 0.003f * tumble_speed;
        tumble_angle_y += 0.005f * tumble_speed;
        
        float cx = cosf(tumble_angle_x);
        float sx = sinf(tumble_angle_x);
        float cy = cosf(tumble_angle_y);
        float sy = sinf(tumble_angle_y);
        float fold_base = sinf(tumble_angle_x * 0.5f) * 0.5f;

        for (size_t i = 0; i < fibers.size(); ++i) {
            const Vector& base = fibers[i];
            
            // Hopf Fiber Params (S2 base)
            float theta = acosf(base.j); // y is up
            float phi = atan2f(base.k, base.i);
            
            // Folding
            float eta = theta / 2.0f;
            float folding_val = sinf(phi * 2.0f + tumble_angle_y + fold_base) * 0.1f * tumble_speed * folding;
            eta += folding_val;
            
            // Twist
            phi += eta * twist;
            
            // Dot Generation
            float phase = i * (PI_F / fibers.size());
            float beta = flow_offset + phase;
            
            // 1. Construct point on S3
            // q = [cos(eta)cos(phi+beta), cos(eta)sin(phi+beta), sin(eta)cos(beta), sin(eta)sin(beta)]
            float q0 = cosf(eta) * cosf(phi + beta);
            float q1 = cosf(eta) * sinf(phi + beta);
            float q2 = sinf(eta) * cosf(beta);
            float q3 = sinf(eta) * sinf(beta);
            
            // 2. Apply Tumble (R_xw, R_yz)
            // R_xw
            float q0_r = q0 * cx - q3 * sx;
            float q3_r = q0 * sx + q3 * cx;
            q0 = q0_r; q3 = q3_r;

            // R_yz
            float q1_r = q1 * cy - q2 * sy;
            float q2_r = q1 * sy + q2 * cy;
            q1 = q1_r; q2 = q2_r;
            
            // 3. Stereographic Projection S3 -> R3
            float div = 1.001f - q3;
            float factor = 1.0f / div;
            Vector v(q0 * factor, q1 * factor, q2 * factor);
            v = v.normalize();
            
            Color4 c = Palettes::richSunset.get(0.0f);
            c.alpha = alpha; // parameter alpha
            
            if (i < prev_positions.size()) {
                const Vector& prev = prev_positions[i];
                // Draw line segment
                auto fragment_shader = [&](const Vector&, Fragment& f) {
                    f.color = c;
                };
                Plot::Line::draw<W, H>(filters, canvas, prev, v, fragment_shader);
                
                prev_positions[i] = v;
            } else {
                // First frame
                filters.plot(canvas, v, c.color, 0, c.alpha);
                prev_positions.push_back(v);
            }
        }
        
        // Render Trails
        filters.flush(canvas, [](const Vector& v, float t) {
            Color4 c = Palettes::richSunset.get(t);
            c.alpha *= (1.0f - t);
            return c;
        }, 1.0f);
    }
    
    // Params
    int num_fibers = 200;
    float flow_speed = 10.0f;
    float tumble_speed = 4.0f;
    float folding = 0.5f;
    float twist = 0.0f;
    float alpha = 0.4f;

private:
    float flow_offset = 0.0f;
    float tumble_angle_x = 0.0f;
    float tumble_angle_y = 0.0f;
    
    std::vector<Vector> fibers;
    std::vector<Vector> prev_positions;
    
    Orientation<W> orientation;
    Timeline<W> timeline;
    
    // Pipeline
    Pipeline<W, H, 
        Filter::World::Trails<W, MAX_TRAILS>, 
        Filter::World::Orient<W>, 
        Filter::Screen::AntiAlias<W, H>> filters;

    
    void init_fibers() {
        fibers.clear();
        prev_positions.clear();
        
        int side = static_cast<int>(ceilf(sqrtf(num_fibers)));
        int rings = side;
        int per_ring = static_cast<int>(ceilf(static_cast<float>(num_fibers) / rings));
        
        fibers.reserve(rings * per_ring);
        
        for (int i = 0; i < rings; ++i) {
            float theta = PI_F * (i + 0.5f) / rings;
            float y = cosf(theta);
            float r = sinf(theta);
            
            for (int j = 0; j < per_ring; ++j) {
                float phi = 2 * PI_F * j / per_ring; // 0..2PI
                
                // Y-up
                float x = r * cosf(phi);
                float z = r * sinf(phi);
                
                fibers.emplace_back(x, y, z);
            }
        }
    }
};
