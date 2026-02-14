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
        persist_pixels = false;
        timeline.add(0, Animation::RandomWalk<W>(orientation, UP)); // Slow continuous spin        
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

    // Helper struct to own the data for a mesh (MeshState is just a view)
    struct OwnedMesh {
        std::vector<Vector> vertices;
        std::vector<uint8_t> face_counts;
        std::vector<int> faces;

        MeshState get_state() const {
             MeshState s;
             s.num_vertices = vertices.size();
             if (s.num_vertices > MeshState::MAX_VERTS) s.num_vertices = MeshState::MAX_VERTS;
             
             for(size_t i=0; i<s.num_vertices; ++i) s.vertices[i] = vertices[i];
             
             s.num_faces = face_counts.size();
             s.face_counts = const_cast<uint8_t*>(face_counts.data());
             s.faces = const_cast<int*>(faces.data());
             return s;
        }

        static OwnedMesh from_poly(const PolyMesh& src) {
            OwnedMesh dst;
            dst.vertices = src.vertices;
            for(const auto& f : src.faces) {
                dst.face_counts.push_back((uint8_t)f.size());
                for(int idx : f) dst.faces.push_back(idx);
            }
            return dst;
        }
    };

    
    int solid_idx = -1;

    void spawn_shape() {
        // 1. Load Next Mesh
        solid_idx = (solid_idx + 1) % (int)SolidType::Last;
        PolyMesh mesh = generate_solid((SolidType)solid_idx);
        
        // 2. Classify Topology
        auto topology = MeshOps::classify_faces_by_topology(mesh);
        
        // 3. Log Shape Name
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

        // 4. Flatten for Rendering
        OwnedMesh owned_mesh = OwnedMesh::from_poly(mesh);
        
        // 5. Prepare Palettes
        std::vector<const Palette*> palettes = {
            &Palettes::embers,
            &Palettes::richSunset,
            &Palettes::brightSunrise,
            &Palettes::bruisedMoss,
            &Palettes::lavenderLake
        };
        static std::mt19937 g(12345 + (int)timeline.t); // Simple seed variation
        std::shuffle(palettes.begin(), palettes.end(), g);
        
        // 5. Create Sprite
        // duration 96, fade 32 (matching JS)
        int duration = 96;
        int fade_in = 32;
        int fade_out = 32;
        
        auto draw_fn = [this, owned_mesh, topology, palettes](Canvas& canvas, float opacity) {
            MeshState mesh_state = owned_mesh.get_state();
            
            // Rotate
            MeshState rotated_mesh;
            rotated_mesh.num_vertices = mesh_state.num_vertices;
            rotated_mesh.num_faces = mesh_state.num_faces;
            rotated_mesh.face_counts = mesh_state.face_counts;
            rotated_mesh.faces = mesh_state.faces;
            
            Quaternion q = orientation.get();
            for(size_t i=0; i<mesh_state.num_vertices; ++i) {
                rotated_mesh.vertices[i] = rotate(mesh_state.vertices[i], q);
            }
            
            auto shader = [&](const Vector& p, Fragment& frag) {
                 int faceIdx = (int)std::round(frag.v2);
                 int topoIdx = 0;
                 if (faceIdx >= 0 && faceIdx < (int)topology.faceColorIndices.size()) {
                     topoIdx = topology.faceColorIndices[faceIdx];
                 }
                 const Palette* pal = palettes[topoIdx % palettes.size()];
                 
                 float distFromEdge = -frag.v1;
                 float size = frag.size;
                 float intensity = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
                 intensity = std::clamp(intensity, 0.0f, 1.0f);
                 
                 frag.color = pal->get(intensity);
                 frag.color.alpha = opacity;
            };
            
            Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader);
        };
        
        timeline.add(0, Animation::Sprite(draw_fn, duration, fade_in, ease_mid, fade_out, ease_mid));
        
        // 6. Schedule Next
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

};

