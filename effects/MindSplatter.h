/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
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
        particle_system(256) // Lower capacity for embedded?
    {
        // Initial build
        rebuild();
        
        // Random Walk
        timeline.add(0, Motion<W>(orientation, random_path, 2000, true)); // Placeholder for RandomWalk
        
        // Warp
        start_warp();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        particle_system.step(canvas);
        
        draw_particles(canvas);
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
    Transition* warp_anim = nullptr;
    
    // Placeholder random path for now (Animation.h needs RandomWalk port?)
    // Using simple procedural rotation.
    ProceduralPath<std::function<Vector(float)>> random_path{
        [](float t) {
             float x = cosf(t * 5) * sinf(t * 3);
             float y = sinf(t * 5) * sinf(t * 3);
             float z = cosf(t * 3);
             return Vector(x, y, z).normalize();
        }
    };

    void rebuild() {
        particle_system.reset(friction, 0.001f);
        
        // Dodecahedron Emitters
        // We can just iterate vertices of Dodecahedron
        // Need Dodecahedron static data.
        // It's in solids.h.
        
        // Dodecahedron Vertices
        // Assuming we can access them like load_mesh_data does
        using EmitSolid = Dodecahedron;
        using AttractSolid = Icosahedron; // Dual
        
        // Add Attractors
        for(const auto& v : AttractSolid::vertices) {
            particle_system.add_attractor(v, well_strength, 0.003f, 0.2f);
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
                
                sys.spawn(axis, vel, Color4(CRGB::White), 160);
            });
        }
    }
    
    void draw_particles(Canvas& canvas) {
        // Iterate active particles
        // Only active particles are in [0, active_count)
        for (int i=0; i < particle_system.active_count; ++i) {
            const auto& p = particle_system.pool[i];
            
            // Draw history
            int len = p.history_length();
            if (len < 2) continue;
            
            // Transform history frames to points
            // p.position is local anchor. 
            // Point = Frame * p.position.
            // Wait, logic in update_particle: `Vector pos = p.orientation.orient(p.position);`
            // `p.orientation` accumulates rotations.
            
            Vector prev_pos;

            
            for (int k=0; k < len; ++k) {
                // history[k] is a Quaternion. New API: get_history(k) where 0 is newest.
                Quaternion q = p.get_history(k);
                Vector pos = rotate(p.position, q);
                
                // Apply Mobius
                pos = mobius_transform(pos, mobius); 
                
                if (k > 0) {
                     // k=0 is newest (head), k=len-1 is oldest (tail)
                     // t should be 0 at head, 1 at tail? Or we want alpha.
                     // Alpha should be 1 at head, 0 at tail.
                     float t_norm = (float)k / (len - 1);
                     float alpha_val = 1.0f - t_norm;
                     
                     // Color: Fade out tail
                     Color4 c = p.color; 
                     c.alpha = alpha_val;
                     
                     auto fragment_shader = [&](const Vector&, const Fragment&){ return c; };
                     Plot::Line::draw<W>(filters, canvas, prev_pos, pos, fragment_shader);
                }
                prev_pos = pos;
            }
        }
    }
    
    void start_warp() {
        // Schedule next warp
         timeline.add(0, RandomTimer(180, 300, [this](Canvas&){
             perform_warp();
         }, false));
    }
    
    void perform_warp() {
        // Target new mobius params
        // Animate mobius.a.re/im etc.
        // Simplified: Just animate `mobius.zoom`? or a standard param.
        // JS MobiusWarp animates complex parameters.
        // For C++, let's just do a zoom bounce or similar.
        
        static float target_zoom = 1.0f;
        target_zoom = (target_zoom == 1.0f) ? 0.5f : 1.0f;
        
        // Add transition
        // timeline.add(0, Transition(mobius.zoom.re, target_zoom, 160, ease_in_out_sine)); 
        // MobiusParams members might likely be std::complex or floats.
        // Ensure 3dmath.h MobiusParams struct supports this.
        // If not, we skip warp animation or just set it.
        start_warp();
    }
};
