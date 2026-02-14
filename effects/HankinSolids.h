/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

#include <vector>
#include <string>

template <int W, int H>
class HankinSolids : public Effect {
public:
    enum class SolidMode {
        Tetrahedron,
        Cube,
        Octahedron,
        Dodecahedron,
        Icosahedron,
        TruncatedTetrahedron,
        Cuboctahedron,
        TruncatedCube,
        TruncatedOctahedron,
        Rhombicuboctahedron,
        TruncatedCuboctahedron,
        SnubCube,
        Icosidodecahedron,
        TruncatedDodecahedron,
        TruncatedIcosahedron,
        Rhombicosidodecahedron,
        TruncatedIcosidodecahedron,
        SnubDodecahedron,
        Last // Sentinel
    };

    // Buffer for morphing operations
    Animation::MorphBuffer morph_buffer;

    // Meshes necessary for cross-fade morphing
    MeshState primary_mesh;
    MeshState secondary_mesh; // Target mesh during morph

    HankinSolids() : Effect(W, H), 
        filters(Filter::Screen::AntiAlias<W, H>())
    {
        persist_pixels = false;
        
        // Continuous Random Walk
        timeline.add(0, Animation::RandomWalk<W>(orientation, Y_AXIS, Animation::RandomWalk<W>::Options::Languid()));
        
        // Init State
        solid_idx = (int)SolidMode::Dodecahedron;
        
        // Initial Geometry Generation
        PolyMesh base = Solids::get(solid_idx);
        if (enable_dual) base = MeshOps::dual(base);
        
        if (enable_hankin) {
            compiled_hankin = MeshOps::compile_hankin(base);
            base = MeshOps::update_hankin<PolyMesh>(compiled_hankin, hankin_angle);
        }
        load_poly_to_mesh(base, primary_mesh);
        
        // Start the Loop
        start_hankin_cycle();
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

    int solid_idx = 0;
    float hankin_angle = PI_F / 4.0f;
    bool enable_dual = false; 
    bool enable_hankin = true;
    float intensity = 1.2f;

    // Data Helpers
    std::vector<std::string> solids_list; // Not really used in C++ enum but logical parity
    CompiledHankin compiled_hankin; // Keep compiled state for the active solid
    
    void start_hankin_cycle() {
        constexpr int DURATION = 64;
        
        // 1. Mutation: Animate Hankin Angle
        // JS: this.timeline.add(0, new Animation.Mutation(this.params, 'hankinAngle', ...))
        // Here we animate member variable hankin_angle
        timeline.add(0, Animation::Mutation(hankin_angle, sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f), DURATION, ease_mid, false)
            .then([this]() {
                this->start_morph_cycle();
            }));
            
        // 2. Sprite: Draw the mesh
        timeline.add(0, Animation::Sprite([this](Canvas& c, float opacity) {
            if (enable_hankin) {
                // Update geometry from compiled hankin
                // JS: this.renderMesh = MeshOps.updateHankin(...)
                PolyMesh updated = MeshOps::update_hankin<PolyMesh>(compiled_hankin, hankin_angle);
                load_poly_to_mesh(updated, primary_mesh);
            }
            draw_mesh(c, primary_mesh, opacity);
        }, DURATION));
    }
    
    void start_morph_cycle() {
        constexpr int DURATION = 16;
        
        // Identify Next Solid
        int next_idx = (solid_idx + 1) % (int)SolidMode::Last;
        
        // Generate Target Mesh (Secondary)
        PolyMesh next_base = Solids::get(next_idx);
        if (enable_dual) next_base = MeshOps::dual(next_base);
        
        PolyMesh next_poly;
        if (enable_hankin) {
            // Must compile for target to get correct topology/vertices
            auto next_compiled = MeshOps::compile_hankin(next_base);
            // Use CURRENT hankin angle for smooth transition
            next_poly = MeshOps::update_hankin<PolyMesh>(next_compiled, hankin_angle); 
        } else {
            next_poly = next_base;
        }
        load_poly_to_mesh(next_poly, secondary_mesh);
        
        // 1. Morph Animation
        timeline.add(0, Animation::MeshMorph(&primary_mesh, &secondary_mesh, &morph_buffer, primary_mesh, secondary_mesh, DURATION, false, ease_in_out_sin)
            .then([this, next_idx]() {
                // Commit State
                this->solid_idx = next_idx;
                
                // Re-compile Hankin for the NEW solid
                PolyMesh new_base = Solids::get(solid_idx);
                if (enable_dual) new_base = MeshOps::dual(new_base);
                
                if (enable_hankin) {
                    compiled_hankin = MeshOps::compile_hankin(new_base);
                    PolyMesh initial = MeshOps::update_hankin<PolyMesh>(compiled_hankin, hankin_angle);
                } 
                
                // Loop
                this->start_hankin_cycle();
            }));
            
        // 2. Sprite: Outgoing
        timeline.add(0, Animation::Sprite([this](Canvas& c, float opacity) {
            draw_mesh(c, primary_mesh, opacity);
        }, DURATION, 0, ease_mid, DURATION, ease_mid));
        
        // 3. Sprite: Incoming
        timeline.add(0, Animation::Sprite([this](Canvas& c, float opacity) {
            draw_mesh(c, secondary_mesh, opacity);
        }, DURATION, DURATION, ease_mid, 0, ease_mid));
    }
    
    // Storage for Primary (0) and Secondary (1) topology
    // REMOVED: faces_store_0, faces_store_1 (MeshState owns data now)

    void load_poly_to_mesh(const PolyMesh& src, MeshState& dst) {
        // Copy Vertices
        dst.vertices = src.vertices;
        
        // Copy Topology
        dst.faces.clear();
        dst.face_counts.clear();
        dst.face_counts.reserve(src.faces.size());
        
        for(const auto& f : src.faces) {
            dst.face_counts.push_back((uint8_t)f.size());
            for(int idx : f) dst.faces.push_back(idx);
        }
        
        // Rebuild BVH (and offsets)
        build_bvh(dst);
    }

    void draw_mesh(Canvas& canvas, const MeshState& mesh, float opacity) {
        if (mesh.vertices.empty()) return;
        
        if (opacity < 0.01f) return;

        // Create rotated copy
        MeshState rotated_mesh = mesh; // Deep copy of vectors
        
        Quaternion q = orientation.get();
        for(auto& v : rotated_mesh.vertices) {
            v = rotate(v, q);
        }
        
        auto shader = [&](const Vector& p, Fragment& f) {
             int faceIdx = (int)std::round(f.v2);
             // Safety check for face index
             if (faceIdx < 0 || faceIdx >= rotated_mesh.face_counts.size()) return;
             
             int n = (int)rotated_mesh.face_counts[faceIdx];
             const Palette* pal = &Palettes::richSunset;
             
             switch (n) {
                 case 3:
                 case 4:
                 case 5:
                 case 8:
                 case 10:
                    pal = &Palettes::lavenderLake;
                    break;
                 case 6:
                    pal = &Palettes::emeraldForest;
                    break;
                 default:
                    pal = &Palettes::richSunset;
                    break;
             }
             
             // Edge Distance Intensity (v1 is -dist)
             float distFromEdge = -f.v1;             
             float size = f.size;
             float normalizedDist = (size > 0.0001f) ? (distFromEdge / size) : 0.0f;
             float t = std::clamp(normalizedDist * intensity, 0.0f, 1.0f);
             
             f.color = pal->get(t);
             f.color.alpha = opacity;
        };
        
        Scan::Mesh::draw<W, H>(filters, canvas, rotated_mesh, shader);
    }
};
