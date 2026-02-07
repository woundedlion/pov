/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../plot.h"
#include "../animation.h"
#include "../filter.h"
#include "../solids.h"
#include "../palettes.h"
#include <vector>

template <int W>
class MindSplatter : public Effect {
public:
    MindSplatter() : Effect(W),
        filters(FilterOrient<W>(orientation), FilterAntiAlias<W>()),
        particle_system(2048)
    {
        rebuild();
        
        // Random Walk
        timeline.add(0, RandomWalk<W>(orientation, Y_AXIS)); 

        start_warp();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        particle_system.step(canvas);
        
        draw_particles(canvas, 1.0f);
    }
    
private:
    Orientation orientation;
    Timeline timeline;
    Pipeline<W, FilterOrient<W>, FilterAntiAlias<W>> filters;
    ParticleSystem particle_system;
    
    // Params
    float friction = 0.85f;
    float well_strength = 1.0f;
    float initial_speed = 0.025f;
    float angular_speed = 0.2f;
    
    // Warp params
    MobiusParams mobius;
    float warp_scale = 0.6f;

    void rebuild() {
        particle_system.reset(friction, 0.001f);
        
        using EmitSolid = Dodecahedron;
        using AttractSolid = Icosahedron; // Dual
        
        // Add Attractors
        for(const auto& v : AttractSolid::vertices) {
            particle_system.add_attractor(v, well_strength, 0.003f, 0.2f);
        }
        
        // Emitter Hues
        static std::vector<float> emitter_hues;
        if(emitter_hues.empty()) {
            for(size_t i=0; i<EmitSolid::NUM_VERTS; ++i) {
                emitter_hues.push_back(static_cast<float>(rand()) / RAND_MAX);
            }
        }
        
        // Add Emitters
        for(size_t i=0; i<EmitSolid::NUM_VERTS; ++i) {
            Vector axis = EmitSolid::vertices[i];
            
            particle_system.add_emitter([=, this](ParticleSystem& sys) mutable {
                static int Counter = 0;
                float angle = (Counter++) * angular_speed;
                
                // Basis
                Vector ref = (std::abs(dot(axis, X_AXIS)) > 0.99f) ? Y_AXIS : X_AXIS;
                Vector u = cross(axis, ref).normalize();
                Vector w = cross(axis, u).normalize();
                
                Vector vel = (u * cosf(angle) + w * sinf(angle)) * initial_speed;
                
                // Update Hue
                emitter_hues[i] = fmodf(emitter_hues[i] + G * 0.1f, 1.0f);
                
                GenerativePalette palette(
                    GradientShape::STRAIGHT,
                    HarmonyType::COMPLEMENTARY,
                    BrightnessProfile::DESCENDING,
                    SaturationProfile::MID,
                    static_cast<int>(emitter_hues[i] * 255.0f)
                );
                
                sys.spawn(axis, vel, palette, 160);
            });
        }
    }
    
    void draw_particles(Canvas& canvas, float opacity = 1.0f) {
        auto vertex_shader = [&](Fragment f) {            
            f.pos = mobius_transform(orientation.orient(f.pos), mobius);
            float holeAlpha = 1.0f;
            for(const auto& attr : particle_system.attractors) {
                 float d = angle_between(f.pos, attr.position);
                  if (d < attr.event_horizon) {
                        float t = d / attr.event_horizon;
                        holeAlpha *= quintic_kernel(t);
                  }
             }
             f.v0 = holeAlpha; // Overwrite progress
             return f;
        };

        auto fragment_shader = [&](const Vector& v, const Fragment& f) -> Fragment {
             float holeAlpha = f.v0; 
             int p_idx = static_cast<int>(f.v2 + 0.5f);
             
             Fragment f_out = f;
             if (p_idx < 0 || p_idx >= particle_system.active_count) {
                 f_out.color = Color4(CRGB::Black, 0.0f);
                 f_out.blend = BLEND_OVER;
                 return f_out;
             }

             const auto& p = particle_system.pool[p_idx];             
             float norm_id = static_cast<float>(p_idx) / particle_system.active_count;
             Color4 c = get_color(p.palette, norm_id);
             
             f_out.color = Color4(c.color, c.alpha * holeAlpha * opacity);
             f_out.blend = BLEND_ADD;
             
             return f_out;
        };
        
        Plot::ParticleSystem::draw<W>(filters, canvas, particle_system, fragment_shader, vertex_shader);
    }
    
    void start_warp() {
        schedule_warp();
    }
    
    void schedule_warp() {
        auto timer = RandomTimer(180, 300, [this](Canvas&) {
            perform_warp();
        });
        timeline.add(0, timer);
    }
    
    void perform_warp() {
        auto warp = MobiusWarp(mobius, warp_scale, 160, false);
        warp.then([this]() {
            schedule_warp();
        });
        timeline.add(0, warp);
    }
};
