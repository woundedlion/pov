/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include "../effects_engine.h"

#include <vector>
#include <string>

template <int W>
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

    HankinSolids() : Effect(W), 
        filters(Filter::Screen::AntiAlias<W>())
    {
        // Continuous Rotation
        timeline.add(0, Animation::Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 8000, ease_mid, true));
        
        // Hankin Angle Mutation
        // sin_wave(from, to, freq, phase)
        timeline.add(0, Animation::Mutation(hankin_angle, sin_wave(0.0f, PI_F / 2.0f, 1.0f, 0.0f), 4000, ease_mid, true));
        
        // Shape Cycling
        timeline.add(0, Animation::PeriodicTimer(2500, [this](Canvas& c) { this->start_morph(); }, true));
        
        // Init
        solid_idx = (int)SolidMode::Dodecahedron;
        reload_solid();
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        // Update Procedural Mesh
        update_mesh();
        
        draw_mesh(canvas, cached_hankin_mesh, 1.0f);
    }

private:
    Orientation orientation;
    Timeline timeline;
    Pipeline<W, Filter::Screen::AntiAlias<W>> filters;

    int solid_idx = 0;
    float hankin_angle = PI_F / 4.0f;
    bool enable_dual = false; // Parameter from JS, defaulting to false
    bool enable_hankin = true;

    // Mesh Data
    PolyMesh base_mesh;
    CompiledHankin compiled_hankin;
    MeshState cached_hankin_mesh; // Renderable
    
    // Morphing State (Simplified: Instant switch for now, full morph requires dual mesh state management)
    // JS does cross-fade or mesh morph. For C++ MVP let's implement instant switch + angle animation first
    // to match the core visual.
    
    void start_morph() {
        solid_idx = (solid_idx + 1) % (int)SolidMode::Last;
        reload_solid();
    }

    void reload_solid() {
        // Regenerate Base Mesh
        switch ((SolidMode)solid_idx) {
            case SolidMode::Tetrahedron: base_mesh = Solids::Platonic::tetrahedron(); break;
            case SolidMode::Cube: base_mesh = Solids::Platonic::cube(); break;
            case SolidMode::Octahedron: base_mesh = Solids::Platonic::octahedron(); break;
            case SolidMode::Dodecahedron: base_mesh = Solids::Platonic::dodecahedron(); break;
            case SolidMode::Icosahedron: base_mesh = Solids::Platonic::icosahedron(); break;
            
            case SolidMode::TruncatedTetrahedron: base_mesh = Solids::Archimedean::truncatedTetrahedron(); break;
            case SolidMode::Cuboctahedron: base_mesh = Solids::Archimedean::cuboctahedron(); break;
            case SolidMode::TruncatedCube: base_mesh = Solids::Archimedean::truncatedCube(); break;
            case SolidMode::TruncatedOctahedron: base_mesh = Solids::Archimedean::truncatedOctahedron(); break;
            case SolidMode::Rhombicuboctahedron: base_mesh = Solids::Archimedean::rhombicuboctahedron(); break;
            case SolidMode::TruncatedCuboctahedron: base_mesh = Solids::Archimedean::truncatedCuboctahedron(); break;
            case SolidMode::SnubCube: base_mesh = Solids::Archimedean::snubCube(); break;
            case SolidMode::Icosidodecahedron: base_mesh = Solids::Archimedean::icosidodecahedron(); break;
            case SolidMode::TruncatedDodecahedron: base_mesh = Solids::Archimedean::truncatedDodecahedron(); break;
            case SolidMode::TruncatedIcosahedron: base_mesh = Solids::Archimedean::truncatedIcosahedron(); break;
            case SolidMode::Rhombicosidodecahedron: base_mesh = Solids::Archimedean::rhombicosidodecahedron(); break;
            case SolidMode::TruncatedIcosidodecahedron: base_mesh = Solids::Archimedean::truncatedIcosidodecahedron(); break;
            case SolidMode::SnubDodecahedron: base_mesh = Solids::Archimedean::snubDodecahedron(); break;
            default: base_mesh = Solids::Platonic::dodecahedron(); break;
        }

        if (enable_dual) {
            base_mesh = MeshOps::dual(base_mesh);
        }

        if (enable_hankin) {
            compiled_hankin = MeshOps::compile_hankin(base_mesh);
        }
    }

    void update_mesh() {
        if (enable_hankin) {
            // Update geometry based on angle
            PolyMesh updated = MeshOps::update_hankin<PolyMesh>(compiled_hankin, hankin_angle);
            convert_to_meshstate(updated, cached_hankin_mesh);
        } else {
            convert_to_meshstate(base_mesh, cached_hankin_mesh);
        }
    }

    void convert_to_meshstate(const PolyMesh& src, MeshState& dst) {
        if (src.vertices.size() > MeshState::MAX_VERTS) {
            dst.num_vertices = 0;
            return;
        }
        dst.num_vertices = src.vertices.size();
        for(size_t i=0; i<src.vertices.size(); ++i) {
            dst.vertices[i] = src.vertices[i];
        }
        dst.num_faces = src.faces.size();
        
        static std::vector<int> faces_buffer;
        static std::vector<uint8_t> counts_buffer;
        faces_buffer.clear();
        counts_buffer.clear();
        
        for(const auto& f : src.faces) {
            counts_buffer.push_back((uint8_t)f.size());
            for(int idx : f) {
                faces_buffer.push_back(idx);
            }
        }
        
        dst.face_counts = counts_buffer.data();
        dst.faces = faces_buffer.data();
    }

    void draw_mesh(Canvas& canvas, const MeshState& mesh, float opacity) {
        MeshState rotated_mesh;
        rotated_mesh.num_vertices = mesh.num_vertices;
        rotated_mesh.num_faces = mesh.num_faces;
        rotated_mesh.face_counts = mesh.face_counts;
        rotated_mesh.faces = mesh.faces;
        
        Quaternion q = orientation.get();
        for(size_t i=0; i<mesh.num_vertices; ++i) {
            rotated_mesh.vertices[i] = rotate(mesh.vertices[i], q);
        }
        
        auto shader = [&](const Vector& p, const Fragment& f) -> Fragment {
             int faceIdx = (int)std::round(f.v2);
             int n = (int)mesh.face_counts[faceIdx];
             
             const Palette* pal = &Palettes::lavenderLake; // Default
             if (n == 6) pal = &Palettes::emeraldForest;
             
             Fragment out = f;
             out.color = pal->get(0.5f);
             out.color.alpha = opacity;
             return out;
        };
        
        Scan::Mesh::draw<W>(filters, canvas, rotated_mesh, shader);
    }
};
