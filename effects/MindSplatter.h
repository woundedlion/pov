/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include "../static_circular_buffer.h"

template <int W, int H>
class MindSplatter : public Effect {
public:
    MindSplatter() : Effect(W, H),
        filters(Filter::World::Orient<W>(orientation), Filter::Screen::AntiAlias<W, H>()),
        particle_system(friction, 0.001f, 25)
    {
        persist_pixels = false;
        timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS)); 
        rebuild();
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

    static const int NUM_PARTICLES = 2048;
    
    typedef Solids::Dodecahedron EmitSolid;
    typedef Solids::Icosahedron AttractSolid;
    
    typedef Animation::ParticleSystem<W, NUM_PARTICLES> ParticleSystem;
    typedef StaticCircularBuffer<GenerativePalette, NUM_PARTICLES> PaletteBuffer;
    
    // Params
    float friction = 0.85f;
    float well_strength = 1.0f;
    float initial_speed = 0.025f;
    float angular_speed = 0.2f;

    Orientation<W> orientation;
    Timeline<W> timeline;
    Pipeline<W, H, Filter::World::Orient<W>, Filter::Screen::AntiAlias<W, H>> filters;
    ParticleSystem particle_system;
    PaletteBuffer palette_buffer;
    std::array<float, EmitSolid::NUM_VERTS> emitter_hues;
    std::array<int, EmitSolid::NUM_VERTS> emit_counters;
    

    // Warp params
    MobiusParams mobius;
    float warp_scale = 0.6f;

    void rebuild() {
        particle_system.reset(friction, 0.001f);
        
        // Add Attractors
        for(const auto& v : AttractSolid::vertices) {
            particle_system.add_attractor(v, well_strength, 0.003f, 0.2f);
        }
        
        // Emitter Hues
        for(size_t i=0; i<EmitSolid::NUM_VERTS; ++i) {
            emitter_hues[i] = hs::rand_f();
        }
        
        palette_buffer.clear();
        emit_counters.fill(0);
        
        // Add Emitters
        for(size_t i=0; i<EmitSolid::NUM_VERTS; ++i) {
            Vector axis = EmitSolid::vertices[i];
            
            particle_system.add_emitter([this, i, axis](ParticleSystem& sys) mutable {
                float angle = (emit_counters[i]++) * angular_speed;
                
                // Basis
                auto basis = make_basis(Quaternion(), axis);
                Vector vel = (basis.u * cosf(angle) + basis.w * sinf(angle)) * initial_speed;
                
                // Update Hue
                emitter_hues[i] = fmodf(emitter_hues[i] + G * 0.1f, 1.0f);
                
                // Spawn with unique palette per particle
                if (particle_system.active_count < NUM_PARTICLES) {
                    auto& pal = palette_buffer.emplace_back(
                        GradientShape::STRAIGHT,
                        HarmonyType::COMPLEMENTARY,
                        BrightnessProfile::FLAT,
                        SaturationProfile::MID,
                        static_cast<uint8_t>(emitter_hues[i] * 255)
                    );
                    particle_system.spawn(axis, vel, pal, 160);
                }
            });
        }
    }
    
    void draw_particles(Canvas& canvas, float opacity = 1.0f) {
        auto vertex_shader = [&](Fragment& f) {
            Vector original_pos = f.pos;
            float holeAlpha = 1.0f;
            for(const auto& attr : particle_system.attractors) {
                 float d = angle_between(original_pos, attr.position);
                  if (d < attr.event_horizon) {
                        float t = d / attr.event_horizon;
                        holeAlpha *= quintic_kernel(t);
                  }
             }

             // Apply Transforms: Mobius THEN Orient
            f.pos = mobius_transform(f.pos, mobius);
            f.pos = orientation.orient(f.pos);
            f.v3 *= holeAlpha; 
        };

        auto fragment_shader = [&](const Vector& v, Fragment& f) {
             float alpha = std::min(f.v0, f.v3);
             int p_idx = static_cast<int>(f.v2 + 0.5f);
             
             if (p_idx < 0 || p_idx >= particle_system.active_count) {
                #ifdef DEBUG
                assert(false);
                #endif
                 f.color = Color4(CRGB(0, 0, 0), 0.0f);
                 return;
             }

             const auto& p = particle_system.pool[p_idx];             
             Color4 c = get_color(p.palette, f.v0);        
             c.alpha  = c.alpha * alpha * alpha * opacity;    
             f.color = c;
        };
        
        Plot::ParticleSystem::draw<W, H>(filters, canvas, particle_system, fragment_shader, vertex_shader);
    }
    
    void start_warp() {
        schedule_warp();
    }
    
    void schedule_warp() {
        auto timer = Animation::RandomTimer(180, 300, [this](Canvas&) {
            perform_warp();
        });
        timeline.add(0, timer);
    }
    
    void perform_warp() {
        auto warp = Animation::MobiusWarp(mobius, warp_scale, 160, false);
        warp.then([this]() {
            schedule_warp();
        });
        timeline.add(0, warp);
    }
};
