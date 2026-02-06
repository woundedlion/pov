/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit permission.
 */
#pragma once

#include "../led.h"
#include "../geometry.h"
#include "../scan.h"
#include "../solids.h"
#include "../animation.h"
#include "../filter.h"
#include "../palettes.h"
#include <vector>
#include <map>

// Custom AlphaFalloffPalette
class AlphaFalloffPalette : public Palette {
public:
    AlphaFalloffPalette(std::function<float(float)> falloff_fn, const Palette& base_palette)
        : falloff_fn(falloff_fn), base_palette(base_palette) {}

    Color4 get(float t) const override {
        Color4 c = base_palette.get().get(t);
        c.alpha *= falloff_fn(t);
        return c;
    }

private:
    std::function<float(float)> falloff_fn;
    std::reference_wrapper<const Palette> base_palette;
};

template <int W>
class DreamBalls : public Effect {
public:
    enum class SolidType {
        Rhombicuboctahedron,
        Rhombicosidodecahedron,
        TruncatedCuboctahedron,
        Icosidodecahedron
    };

    struct Params {
        SolidType solid_type;
        int num_copies;
        float offset_radius;
        float offset_speed;
        float warp_scale;
        const Palette* palette;
        float alpha;
        bool enable_slice = false;
    };

    DreamBalls() : Effect(W), 
        filters(FilterOrientSlice<W>(orientations, Y_AXIS), FilterAntiAlias<W>()),
        slice_filter(filters) // Filters inherits Head (FilterOrientSlice)
    {
        // Initialize Orientations
        orientations.resize(2); // 2 slices

        // Initialize Presets
        setup_presets();
        
        // Initial load
        set_preset(0); // Load first preset
        
        timeline.add(0, PeriodicTimer(160, [this](Canvas& c) { this->spin_slices(); }, true));
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        
        slice_filter.get().enabled = params.enable_slice;        
        timeline.step(canvas);
        
        // Warp Animation Step
        warp_anim.scale = params.warp_scale;
        warp_anim.step(canvas);
        
        // Draw
        t += 0.01f;
        draw_scene(canvas, params, 1.0f);
    }
    
    Params params; 

private:
    // State
    float t = 0;
    
    // Mesh Data (Dynamic)
    MeshState base_mesh;
    MeshState displaced_mesh;
    struct Tangent { Vector u; Vector v; };
    std::vector<Tangent> tangents;
    
    // Animation Integ
    MobiusParams mobius_params;
    MobiusWarp warp_anim{mobius_params, 1.0f, 200, true};
    Timeline timeline;
    
    // Slices
    std::vector<Orientation> orientations;
    
    // Pipeline
    Pipeline<W, FilterOrientSlice<W>, FilterAntiAlias<W>> filters;
    std::reference_wrapper<FilterOrientSlice<W>> slice_filter;

    // Presets
    std::vector<Params> presets;
    int current_preset_idx = 0;
    
    AlphaFalloffPalette bloodStreamFalloff {
        [](float t){ return 1.0f - t; }, 
        Palettes::bloodStream
    };

    void setup_presets() {
        presets.push_back({
            SolidType::Rhombicuboctahedron, 18, 0.3f, 0.4f, 0.3f, 
            &bloodStreamFalloff, 
            0.7f
        });
        
        presets.push_back({
             SolidType::Rhombicosidodecahedron, 6, 0.05f, 1.0f, 1.8f,
             &bloodStreamFalloff,
             0.7f
        });

        presets.push_back({
             SolidType::TruncatedCuboctahedron, 6, 0.16f, 1.0f, 2.0f,
             &Palettes::richSunset,
             0.3f
        });
        
        presets.push_back({
             SolidType::Icosidodecahedron, 10, 0.16f, 1.0f, 0.5f,
             &Palettes::lavenderLake,
             0.3f
        });
    }

    void set_preset(int idx) {
        if (idx < 0 || static_cast<size_t>(idx) >= presets.size()) return;
        current_preset_idx = idx;
        params = presets[idx];
        load_solid(params.solid_type);
    }

    void load_solid(SolidType type) {
        switch(type) {
            case SolidType::Rhombicuboctahedron: load_mesh_data<Rhombicuboctahedron>(); break;
            case SolidType::Rhombicosidodecahedron: load_mesh_data<Rhombicosidodecahedron>(); break;
            case SolidType::TruncatedCuboctahedron: load_mesh_data<TruncatedCuboctahedron>(); break;
            case SolidType::Icosidodecahedron: load_mesh_data<Icosidodecahedron>(); break;
        }
    }

    template <typename T>
    void load_mesh_data() {
        // Copy vertices (mutable)
        if (T::NUM_VERTS > MeshState::MAX_VERTS) return;
        std::copy(T::vertices.begin(), T::vertices.end(), base_mesh.vertices.begin());
        base_mesh.num_vertices = T::NUM_VERTS;

        // Point to static topology (immutable)
        base_mesh.num_faces = T::NUM_FACES;
        base_mesh.face_counts = T::face_counts.data();
        base_mesh.faces = T::faces.data();
        
        // Prepare Displaced Mesh
        displaced_mesh.num_vertices = base_mesh.num_vertices;
        displaced_mesh.num_faces = base_mesh.num_faces;
        displaced_mesh.face_counts = base_mesh.face_counts;
        displaced_mesh.faces = base_mesh.faces;

        // Compute Tangents
        tangents.clear();
        tangents.reserve(base_mesh.num_vertices);
        for(const auto& p : base_mesh.vertices) {
            Vector axis = (std::abs(p.j) > 0.99f) ? X_AXIS : Y_AXIS;
            Vector u = cross(p, axis).normalize();
            Vector v = cross(p, u).normalize();
            tangents.push_back({u, v});
        }
    }

    void update_displaced_mesh(const Params& p, float angle_offset) {
        size_t count = base_mesh.num_vertices;
        float r = p.offset_radius;
        float speed = p.offset_speed;

        for(size_t i=0; i<count; ++i) {
            const Vector& v = base_mesh.vertices[i];
            const Tangent& tan = tangents[i];

            float phase = i * 0.1f;
            float angle = t * speed * 2 * PI_F + phase + angle_offset;

            float cosA = cosf(angle);
            float sinA = sinf(angle);

            Vector disp = v + (tan.u * cosA + tan.v * sinA) * r;
            displaced_mesh.vertices[i] = disp.normalize();
        }
    }

    void draw_mesh(Canvas& canvas, const MeshState& mesh, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader) {
        std::vector<std::pair<int, int>> edges;
        edges.reserve(mesh.num_faces * 4); // Avg 4 edges per face

        int offset = 0;
        for (size_t i = 0; i < mesh.num_faces; ++i) {
            int count = mesh.face_counts[i];
            for (int k = 0; k < count; ++k) {
                int u = mesh.faces[offset + k];
                int v = mesh.faces[offset + (k + 1) % count];
                if (u > v) std::swap(u, v);
                edges.push_back({u, v});
            }
            offset += count;
        }
        
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

        for (const auto& edge : edges) {
            Plot::Line::draw<W>(filters, canvas, 
                mesh.vertices[edge.first], mesh.vertices[edge.second], 
                fragment_shader, vertex_shader);
        }
    }

    void draw_scene(Canvas& canvas, const Params& p, float opacity) {
        
        auto vertex_shader = [&](Fragment f) {
            f.pos = mobius_transform(f.pos, mobius_params);
            return f;
        };

        auto fragment_shader = [&](const Vector& v, const Fragment& f) {
            Color4 c = p.palette->get(f.v0); 
            c.alpha *= p.alpha * opacity;
            return c;
        };

        for(int i=0; i<p.num_copies; ++i) {
            float offset = (static_cast<float>(i) / p.num_copies) * 2 * PI_F;
            update_displaced_mesh(p, offset);
            draw_mesh(canvas, displaced_mesh, fragment_shader, vertex_shader);
        }
    }

    void spin_slices() {
        Vector axis = random_vector();
        slice_filter.get().axis = axis; // Update axis
        
        for (int i = 0; i < orientations.size(); ++i) {
            float direction = (i % 2 == 0) ? 1.0f : -1.0f;
            timeline.add(0, Rotation<W>(orientations[i], axis, direction * 2 * PI_F, 80, ease_in_out_sin, false));
        }
    }
};
