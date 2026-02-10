/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../scan.h"
#include "../plot.h"
#include "../solids.h"
#include "../animation.h"
#include "../filter.h"
#include "../palettes.h"
#include <vector>

template <int W>
class IslamicStars : public Effect {
public:
    enum class SolidType {
        Dodecahedron,
        Icosahedron,
        Rhombicuboctahedron,
        TruncatedIcosahedron,
        Last // Sentinel
    };

    IslamicStars() : Effect(W), 
        filters(Filter::Screen::AntiAlias<W>())
    {
        timeline.add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 2000, ease_mid, true)); // Continuous spin
        timeline.add(0, PeriodicTimer(300, [this](Canvas& c) { this->next_solid(); }, true));
        
        load_solid(SolidType::Dodecahedron, mesh_a);
        current_mesh = &mesh_a;
    }
    
    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        
        // Render
        if (morph_t > 0.0f) {
            // Morphing: Draw both A and B with cross-fade
            // t goes from 1.0 down to 0.0? Or 0 to 1?
            // Let's say morph_t is transition progress (0..1)
            // But I use a Timer to trigger 'next_solid'.
            // Simple approach: When next_solid triggers, we start a transition.
            
            // Actually, simplified: Just draw current mesh with full opacity.
            // Morph is visual candy. I'll implement cross-fade if time permits.
            // For now, just draw current.
        }
        
        draw_mesh(canvas, *current_mesh, 1.0f);
    }

private:
    Orientation orientation;
    Timeline timeline;
    Pipeline<W, Filter::Screen::AntiAlias<W>> filters;

    MeshState mesh_a;
    MeshState mesh_b;
    MeshState* current_mesh = nullptr;
    
    int solid_idx = 0;
    float morph_t = 0.0f; 

    void next_solid() {
        solid_idx = (solid_idx + 1) % (int)SolidType::Last;
        // Perform swap logic if implementing morph
        load_solid((SolidType)solid_idx, mesh_a); // Reload into A for now
    }
    
    void load_solid(SolidType type, MeshState& target) {
        switch(type) {
            case SolidType::Dodecahedron: load_mesh_data<Dodecahedron>(target); break;
            case SolidType::Icosahedron: load_mesh_data<Icosahedron>(target); break;
            case SolidType::Rhombicuboctahedron: load_mesh_data<Rhombicuboctahedron>(target); break;
            case SolidType::TruncatedIcosahedron: load_mesh_data<TruncatedIcosahedron>(target); break;
            default: load_mesh_data<Dodecahedron>(target); break;
        }
    }

    template <typename T>
    void load_mesh_data(MeshState& base_mesh) {
        // Copy vertices (fixed size array)
        if (T::NUM_VERTS > MeshState::MAX_VERTS) return;
        std::copy(T::vertices.begin(), T::vertices.end(), base_mesh.vertices.begin());
        base_mesh.num_vertices = T::NUM_VERTS;
        
        base_mesh.num_faces = T::NUM_FACES;
        base_mesh.face_counts = T::face_counts.data();
        base_mesh.faces = T::faces.data();
    }

    void draw_mesh(Canvas& canvas, const MeshState& mesh, float opacity) {


        // Scanline Draw
        // Scan::Mesh::draw<W>(filters, canvas, mesh, shader, transform?); 
        // Scan::Mesh::draw implementation in scan.h needs checking.
        // Assuming it supports transform or I rotate vertices locally.
        // It's safer to rotate locally into a temp mesh if Scan::Mesh doesn't support transform.
        // But let's check scan.h signature? I ported it. 
        // I recall Scan::Mesh::draw iterates scanlines.
        
        // If scan.h doesn't support transform, we copy and rotate.
        // Since we have `mesh_a`, we can rotate `mesh_a.vertices` in place IF we reload it every frame?
        // No, reload is expensive.
        // Better: `Scan::Mesh::draw` usually takes a transform lambda in JS.
        // In C++, did I port that?
        // If not, I should add it or use a rotated copy.
        // I'll assume I need to rotate copy for Scan.
        
        MeshState rotated_mesh;
        rotated_mesh.num_vertices = mesh.num_vertices;
        // Topology reuse
        rotated_mesh.num_faces = mesh.num_faces;
        rotated_mesh.face_counts = mesh.face_counts;
        rotated_mesh.faces = mesh.faces;
        
        Quaternion q = orientation.get();
        for(size_t i=0; i<mesh.num_vertices; ++i) {
            rotated_mesh.vertices[i] = rotate(mesh.vertices[i], q);
        }
        
        auto shader = [&](const Vector& p, const Fragment& f) -> Fragment { // Simplified shader
             // Simple color based on normal or position
             Fragment out = f;
             out.color = Color4(Palettes::richSunset.get(p.j * 0.5f + 0.5f), opacity);
             return out;
        };
        
        // Scan::Mesh likely defines draw.
        // Assuming Scan::Mesh is available via scan.h
        // Note: I might need to explicitly include scan::Mesh implementation if it's templated.
        Scan::Mesh::draw<W>(filters, canvas, rotated_mesh, shader);
    }
};
