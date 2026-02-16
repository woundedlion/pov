/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <random>

template <int W, int H>
class IslamicStars : public Effect {
public:
    enum class SolidType {
        Icosahedron_Hk59_Bitruncate033,
        Octahedron_Hk17_Ambo_Hk72,
        Icosahedron_Kis_Gyro,
        TruncatedIcosidodecahedron_Truncate05_Ambo_Dual,
        Icosidodecahedron_Truncate05_Ambo_Dual,
        SnubDodecahedron_Truncate05_Ambo_Dual,
        Octahedron_Hk34_Ambo_Hk72,
        Rhombicuboctahedron_Hk63_Ambo_Hk63,
        TruncatedIcosahedron_Hk54_Ambo_Hk72,
        Dodecahedron_Hk54_Ambo_Hk72,
        Dodecahedron_Hk72_Ambo_Dual_Hk20,
        TruncatedIcosahedron_Truncate05_Ambo_Dual,
        Last // Sentinel
    };

    IslamicStars() : Effect(W, H), 
        filters(Filter::Screen::AntiAlias<W, H>())
    {
        registerParam("Duration", &params.duration, 48.0f, 192.0f);
        registerParam("Fade", &params.fade, 16.0f, 64.0f);
        registerParam("Ripp Amp", &ripple.amplitude, 0.0f, 0.5f);
        registerParam("Ripp Freq", &ripple.frequency, 1.0f, 100.0f);
        registerParam("Ripp Decay", &ripple.decay, 0.1f, 5.0f);

        persist_pixels = false;
        timeline.add(0, Animation::RandomWalk<W>(orientation, UP)); // Slow continuous spin
        
        // Init Ripple Defaults
        ripple.amplitude = 0.25f; 
        ripple.frequency = 8.0f; // Frequency of the source oscillations
        ripple.decay = 0.7f; // Spatial decay of the packet
        ripple.lifespawn = 0.0f; // Start invisible
        
        // Example: Trigger a ripple every 60 frames at a random location
        timeline.add(0, Animation::PeriodicTimer(120, [this](Canvas&) {
            Vector hit = random_vector();
            // Add the Ripple animation to the timeline
            // It will hijack 'this->ripple' params for 120 frames
            timeline.add(0, Animation::Ripple(this->ripple, hit, 0.4f, 120));
        }, true));
        
        spawn_shape();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
    }

private:
    Orientation<W> orientation;
    Timeline<W> timeline;
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> filters;
    RippleParams ripple;

    MeshState poly_to_state(const PolyMesh& src) const {
        MeshState dst;
        dst.vertices = src.vertices;
        
        dst.faces.clear();
        dst.face_counts.clear();
        dst.face_counts.reserve(src.faces.size());
        
        for(const auto& f : src.faces) {
            dst.face_counts.push_back((uint8_t)f.size());
            for(int idx : f) dst.faces.push_back(idx);
        }

        build_bvh(dst);        
        return dst;
    }

    
    int solid_idx = -1;

    void spawn_shape() {
        // Load Next Mesh
        solid_idx = (solid_idx + 1) % (int)SolidType::Last;
        PolyMesh mesh = generate_solid((SolidType)solid_idx);
        
        // Classify Topology
        auto faceIndices = MeshOps::classify_faces_by_topology(mesh);
        
        // Log Shape Name
        const char* names[] = {
            "Icosahedron_Hk59_Bitruncate033",
            "Octahedron_Hk17_Ambo_Hk72",
            "Icosahedron_Kis_Gyro",
            "TruncatedIcosidodecahedron_Truncate05_Ambo_Dual",
            "Icosidodecahedron_Truncate05_Ambo_Dual",
            "SnubDodecahedron_Truncate05_Ambo_Dual",
            "Octahedron_Hk34_Ambo_Hk72",
            "Rhombicuboctahedron_Hk63_Ambo_Hk63",
            "TruncatedIcosahedron_Hk54_Ambo_Hk72",
            "Dodecahedron_Hk54_Ambo_Hk72",
            "Dodecahedron_Hk72_Ambo_Dual_Hk20",
            "TruncatedIcosahedron_Truncate05_Ambo_Dual"
        };
        if(solid_idx >= 0 && solid_idx < (int)SolidType::Last) {
            hs::log("Spawning Shape: %s (V=%d, F=%d)", names[solid_idx], (int)mesh.vertices.size(), (int)mesh.faces.size());
        }

        // Flatten for Rendering
        MeshState mesh_state = poly_to_state(mesh);
        
        // Prepare Palettes
        std::vector<const Palette*> palettes = {
            &Palettes::embers,
            &Palettes::richSunset,
            &Palettes::brightSunrise,
            &Palettes::bruisedMoss,
            &Palettes::lavenderLake
        };
        static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
        std::shuffle(palettes.begin(), palettes.end(), g);
        
        // Create Sprite
        int duration = 96;
        int fade_in = 32;
        int fade_out = 32;
        
        auto draw_fn = [this, mesh_state, faceIndices, palettes](Canvas& canvas, float opacity) {
            // Rotate
            MeshState rotated_mesh = mesh_state; // Deep copy vectors
            Quaternion q = orientation.get();

            // Capture the current state of ripple params
            RippleParams current_ripple = this->ripple; 

            for(auto& v : rotated_mesh.vertices) {
                // 1. Standard Rotation
                v = rotate(v, q);
                
                // 2. Apply Ripple Transform
                // Only if ripple is active to save cycles
                if (current_ripple.lifespawn > 0.001f) {
                    v = ripple_transform(v, current_ripple);
                }
            }
            
            auto fragment_shader = [&](const Vector& p, Fragment& frag) {
                 int faceIdx = (int)std::round(frag.v2);
                 int topoIdx = 0;
                 if (faceIdx >= 0 && faceIdx < (int)faceIndices.size()) {
                     topoIdx = faceIndices[faceIdx];
                 }
                 const Palette* pal = palettes[topoIdx % palettes.size()];
                 
                 float distFromEdge = -frag.v1;
                 float size = frag.size;
                 float intensity = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
                 intensity = std::clamp(intensity, 0.0f, 1.0f);
                 
                 frag.color = pal->get(intensity);
                 frag.color.alpha = opacity;
            };
            
            Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, fragment_shader);
        };
        
        timeline.add(0, Animation::Sprite(draw_fn, duration, fade_in, ease_mid, fade_out, ease_mid));
        
        // Schedule Next
        // Overlap = fade. Next delay = duration - overlap.
        int next_delay = duration - fade_out;
        timeline.add(next_delay, Animation::PeriodicTimer(0, [this](Canvas&){ this->spawn_shape(); }, false));
    }
    
    PolyMesh generate_solid(SolidType type) {
        switch(type) {
            case SolidType::Icosahedron_Hk59_Bitruncate033: return Solids::IslamicStarPatterns::icosahedron_hk59_bitruncate033();
            case SolidType::Octahedron_Hk17_Ambo_Hk72: return Solids::IslamicStarPatterns::octahedron_hk17_ambo_hk72();
            case SolidType::Icosahedron_Kis_Gyro: return Solids::IslamicStarPatterns::icosahedron_kis_gyro();
            case SolidType::TruncatedIcosidodecahedron_Truncate05_Ambo_Dual: return Solids::IslamicStarPatterns::truncatedIcosidodecahedron_truncate05_ambo_dual();
            case SolidType::Icosidodecahedron_Truncate05_Ambo_Dual: return Solids::IslamicStarPatterns::icosidodecahedron_truncate05_ambo_dual();
            case SolidType::SnubDodecahedron_Truncate05_Ambo_Dual: return Solids::IslamicStarPatterns::snubDodecahedron_truncate05_ambo_dual();
            case SolidType::Octahedron_Hk34_Ambo_Hk72: return Solids::IslamicStarPatterns::octahedron_hk34_ambo_hk72();
            case SolidType::Rhombicuboctahedron_Hk63_Ambo_Hk63: return Solids::IslamicStarPatterns::rhombicuboctahedron_hk63_ambo_hk63();
            case SolidType::TruncatedIcosahedron_Hk54_Ambo_Hk72: return Solids::IslamicStarPatterns::truncatedIcosahedron_hk54_ambo_hk72();
            case SolidType::Dodecahedron_Hk54_Ambo_Hk72: return Solids::IslamicStarPatterns::dodecahedron_hk54_ambo_hk72();
            case SolidType::Dodecahedron_Hk72_Ambo_Dual_Hk20: return Solids::IslamicStarPatterns::dodecahedron_hk72_ambo_dual_hk20();
            case SolidType::TruncatedIcosahedron_Truncate05_Ambo_Dual: return Solids::IslamicStarPatterns::truncatedIcosahedron_truncate05_ambo_dual();
            default: return Solids::Platonic::dodecahedron();
        }
    }

    struct Params {
        float duration = 96.0f;
        float fade = 32.0f;
    } params;
};
