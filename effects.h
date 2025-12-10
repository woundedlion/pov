#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "effects_engine.h"

///////////////////////////////////////////////////////////////////////////////

template <int W>
class RingShower : public Effect {
public:

    RingShower() :
        Effect(W)
    {
        persist_pixels = false;

        timeline.add(0,
            RandomTimer(1, 24,
                [this](auto&) { this->spawn_ring(); },
                true)
        );

    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
    }

private:

    struct Ring {
        Vector normal;
        float speed;
        float radius;
        float duration;
        GenerativePalette palette;

        Ring() :
            normal(random_vector()),
            speed(0.0),
            radius(0.0),
            duration(0.0),
            palette(GradientShape::CIRCULAR, HarmonyType::ANALOGOUS, BrightnessProfile::FLAT)
        {
        }
    };

    void spawn_ring() {
        for (size_t i = 0; i < MAX_RINGS; ++i) {
            if (rings[i].duration <= 0) {
                Ring& ring = rings[i];
                ring.normal = random_vector();
                ring.duration = hs::rand_int(16, 96);
                ring.radius = 0;
                ring.palette = GenerativePalette(
                    GradientShape::CIRCULAR,
                    HarmonyType::ANALOGOUS,
                    BrightnessProfile::FLAT);

                timeline.add(0,
                    Sprite(
                        [this, i](Canvas& canvas, float opacity) { this->draw_ring(canvas, opacity, i); },
                        ring.duration,
                        4, ease_mid,
                        0, ease_mid)
                );

                timeline.add(0,
                    Transition(ring.radius, 2, ring.duration, ease_mid, false, false)
                    .then([&ring]() { ring.duration = 0; })
                );

                return;
            }
        }
    }

    void draw_ring(Canvas& canvas, float opacity, size_t index) {
        Ring& ring = rings[index];
        dots.clear();
        ::draw_ring<W>(dots, orientation.get(), ring.normal, ring.radius,
            [&](auto& v, auto t) {
                return ring.palette.get(t);
            },
            0);
        plot_dots<W>(dots, filters, canvas, 0, opacity * alpha);
    }


    static constexpr size_t MAX_RINGS = 8;
    Ring rings[MAX_RINGS];
    Pipeline<W, FilterAntiAlias<W>> filters;
    Orientation orientation;
    Timeline timeline;
    Dots dots;
    static constexpr float alpha = 0.2;
};

///////////////////////////////////////////////////////////////////////////////

template <int W>
class Comets : public Effect {
public:
    Comets() :
        Effect(W),
        palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS, BrightnessProfile::ASCENDING),
        next_palette(GradientShape::STRAIGHT, HarmonyType::ANALOGOUS, BrightnessProfile::ASCENDING),
        trails(trail_length),
        filters(
            FilterOrient<W>(orientation),
            FilterChromaticShift<W>(),
            FilterAntiAlias<W>())
    {
        persist_pixels = false;
        update_path();

        for (int i = 0; i < NUM_NODES; ++i) {
            spawn_node(i);
        }

        timeline.add(0,
            PeriodicTimer(2 * cycle_duration, [this](auto&) {
                increment_path();
                update_palette();
                }, true)
        );

        timeline.add(0, RandomWalk<W>(orientation, random_vector()));
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);

        trails.render(canvas, filters,
            [this](const Vector& v, float t) {
                return palette.get(1.0f - t);
            }
        );
    }

private:

    struct Node {
        Orientation orientation;
        Vector v;
        const Path<W>& path;

        Node(const Path<W>& p) :
            v(Y_AXIS),
            path(p)
        {
        }
    };

    void increment_path() {
        cur_function = (cur_function + 1) % functions.size();
        update_path();
    }

    void update_path() {
        const auto& f = functions[cur_function];
        path.collapse();
        path.append_segment([f](float t) -> Vector {
            return lissajous(f.m1, f.m2, f.a, t);
            }, f.domain, 1024, ease_mid);
    }

    void update_palette() {
        next_palette = GenerativePalette(
            GradientShape::STRAIGHT,
            HarmonyType::ANALOGOUS,
            BrightnessProfile::ASCENDING
        );

        timeline.add(0,
            ColorWipe(palette, next_palette, wipe_duration, ease_mid)
        );
    }

    void spawn_node(int i) {
        nodes.emplace_back(this->path);

        timeline.add(0,
            Sprite(
                [this, i](Canvas& canvas, float opacity) { draw_node(canvas, opacity, i); },
                -1,
                16, ease_mid,
                0, ease_mid
            )
        );

        timeline.add(i * spacing,
            Motion<W>(nodes[i].orientation, nodes[i].path, cycle_duration, true)
        );
    }

    void draw_node(Canvas& canvas, float opacity, int i) {
        Node& node = nodes[i];
        tween(node.orientation, [this, &canvas, opacity, &node](auto& q, auto t) {
            dots.clear();
            auto v = rotate(node.v, node.orientation.get()).normalize();
            ::draw_vector<W>(dots, v,
                [this](const auto& v, auto t) {
                    return palette.get(1.0f - t);
                });

            // Record 3D dots into DecayBuffer instead of plotting immediately
            trails.record(dots, t, this->alpha * opacity);
            });
        node.orientation.collapse();
    }

    static constexpr int NUM_NODES = 1;
    float alpha = 1.0f;
    size_t cycle_duration = 80;
    size_t trail_length = 80;
    size_t wipe_duration = 48;
    size_t spacing = 48;

    const std::array<LissajousParams, 13> functions = { {
      {1.06f, 1.06f, 0.0f, 5.909f},
      {6.06, 1.0f, 0.0f, 2.0f * PI_F},
      {1.0f, 5.0f, 0.0f, 2.0f * PI_F},
      {6.02, 4.01, 0.0f, 3.132f},
      {46.62f, 62.16f, 0.0f, 0.404f},
      {46.26f, 69.39f, 0.0f, 0.272f},
      {19.44f, 9.72f, 0.0f, 0.646f},
      {8.51f, 17.01f, 0.0f, 0.739f},
      {7.66f, 6.38f, 0.0f, 4.924f},
      {8.75f, 5.0f, 0.0f, 5.027f},
      {11.67f, 14.58f, 0.0f, 2.154f},
      {11.67f, 8.75f, 0.0f, 2.154f},
      {10.94f, 8.75f, 0.0f, 2.872f}
    } };

    Path<W> path;
    size_t cur_function = 0;
    Orientation orientation;
    GenerativePalette palette;
    GenerativePalette next_palette;

    DecayBuffer<W, 3000> trails;

    Pipeline<W,
        FilterOrient<W>,
        FilterChromaticShift<W>,
        FilterAntiAlias<W>
    > filters;

    StaticCircularBuffer<Node, NUM_NODES> nodes;
    Timeline timeline;
    Dots dots;
};

///////////////////////////////////////////////////////////////////////////////

template <int W>
class RingSpin : public Effect {
public:

    struct Ring {
        Ring(const Vector& normal, const Palette& palette, uint8_t trail_length) :
            normal(normal),
            palette(palette),
            filters(
                FilterDecay<W, 6000>(trail_length),
                FilterAntiAlias<W>()
            )
        {
        }

        const Vector normal;
        const Palette& palette;
        Orientation orientation;

        // UPDATED: Use StaticCircularBuffer instead of std::vector
        Points base_points;

        Pipeline<W,
            FilterDecay<W, 6000>,
            FilterAntiAlias<W>
        > filters;
    };

    RingSpin() :
        Effect(W)
    {
        persist_pixels = false;
        rings.reserve(NUM_RINGS);
        for (int i = 0; i < NUM_RINGS; ++i) {
            spawn_ring(X_AXIS, *palettes[i]);
        }
    }

    bool show_bg() const { return false; }

    void spawn_ring(const Vector& normal, const Palette& palette) {
        auto ring_index = rings.size();
        rings.emplace_back(normal, palette, trail_length);
        sample_ring<W>(rings.back().base_points, Quaternion(), normal, 1.0f);

        timeline.add(0,
            Sprite(
                [this, ring_index](Canvas& canvas, float opacity) { draw_ring(canvas, opacity, ring_index); },
                -1,
                4, ease_mid,
                0, ease_mid
            ));

        timeline.add(0,
            RandomWalk<W>(rings[ring_index].orientation, rings[ring_index].normal));
    }

    void draw_ring(Canvas& canvas, float opacity, size_t ring_index) {
        auto& ring = rings[ring_index];
        tween(ring.orientation, [this, &canvas, opacity, &ring](auto& q, auto t) {
            dots.clear();
            Points points;
            for (const auto& p : ring.base_points) {
                points.push_back(rotate(p, q));
            }
            rasterize<W>(dots, points, [&](auto& v, auto t) { return vignette(ring.palette)(0); }, true);
            plot_dots<W>(dots, ring.filters, canvas, t, alpha * opacity);
            });
        ring.orientation.collapse();

        ring.filters.trail(canvas,
            [&](float x, float y, float t) { return vignette(ring.palette)(t); },
            alpha * opacity);
    }

    void draw_frame() {
        Canvas canvas(*this);
        timeline.step(canvas);
    }

private:

    static constexpr int NUM_RINGS = 4;
    std::vector<Ring> rings;
    static constexpr float alpha = 0.2;
    static constexpr float trail_length = 22;
    std::array<const Palette*, NUM_RINGS> palettes = { &richSunset, &mangoPeel, &undersea, &iceMelt };
    Timeline timeline;
    Dots dots;
};


template <int W>
class Dynamo : public Effect {
public:

    struct Node {
        Node() :
            x(0), y(0), v(0)
        {
        }

        int x;
        int y;
        int v;
    };

    Dynamo() :
        Effect(W),
        palettes({ GenerativePalette(
          GradientShape::VIGNETTE,
          HarmonyType::ANALOGOUS,
          BrightnessProfile::ASCENDING) }),
        palette_normal(Z_AXIS),
        speed(2),
        gap(3),
        trail_length(8),
        trails(trail_length),
        filters(
            FilterReplicate<W>(3),
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;

        for (size_t i = 0; i < NUM_NODES; ++i) {
            nodes[i].y = i;
        }

        timeline
            .add(0,
                RandomTimer(4, 64, [this](auto&) { reverse(); }, true))
            .add(0,
                RandomTimer(160, 160, [this](auto&) { color_wipe(); }, true))
            .add(0,
                RandomTimer(48, 160, [this](auto&) { rotate(); }, true));
    }

    bool show_bg() const { return false; }

    void reverse() {
        speed *= -1;
    }

    void rotate() {
        timeline.add(0,
            Rotation<W>(
                orientation, random_vector(), PI_F, 40, ease_in_out_sin, false));
    }

    void color_wipe() {
        if (palettes.is_full()) {
            Serial.println("Palettes full, dropping color wipe!");
            return;
        }
        palettes.push_front(GenerativePalette(
            GradientShape::VIGNETTE,
            HarmonyType::ANALOGOUS,
            BrightnessProfile::ASCENDING));
        palette_boundaries.push_front(0);
        timeline.add(0,
            Transition(palette_boundaries.front(), PI_F + 1, wipe_duration, ease_mid)
            .then([this]() {
                palette_boundaries.pop_back();
                palettes.pop_back();
                })
        );
    }

    Pixel color(const Vector& v, float t) {
        constexpr float blend_width = PI_F / 4;
        float a = angle_between(v, palette_normal);

        for (size_t i = 0; i < palette_boundaries.size(); ++i) {
            float boundary = palette_boundaries[i];
            auto lower_edge = boundary - blend_width;
            auto upper_edge = boundary + blend_width;

            if (a < lower_edge) {
                return std::visit([t](auto& p) { return p.get(1 - t); }, palettes[i]);
            }

            if (a >= lower_edge && a <= upper_edge) {
                auto blend_factor = (a - lower_edge) / (2 * blend_width);
                auto clamped_blend_factor = std::clamp(blend_factor, 0.0f, 1.0f);
                auto c1 = std::visit([t](auto& p) { return p.get(1 - t); }, palettes[i]);
                auto c2 = std::visit([t](auto& p) { return p.get(1 - t); }, palettes[i + 1]);
                return c1.lerp16(c2, to_short(clamped_blend_factor));
            }

            auto next_boundary_lower_edge = (i + 1 < palette_boundaries.size()
                ? palette_boundaries[i + 1] - blend_width
                : PI_F + 1);
            if (a > upper_edge && a < next_boundary_lower_edge) {
                return std::visit([t](auto& p) { return p.get(1 - t); }, palettes[i + 1]);
            }
        }

        return std::visit([t](auto& p) { return p.get(1 - t); }, palettes[0]);
    }

    void draw_frame() {
        Canvas canvas(*this);

        for (int i = std::abs(speed) - 1; i >= 0; --i) {
            pull(0);
            draw_nodes(canvas, static_cast<float>(i) / std::abs(speed));
        }

        trails.render(canvas, filters,
            [this](const Vector& v, float t) {
                return color(v, t);
            }
        );

        timeline.step(canvas);
    }

private:

    float node_y(const Node& node) const {
        return (static_cast<float>(node.y) / (nodes.size() - 1)) * (H_VIRT - 1);
    }

    void draw_nodes(Canvas& canvas, float age) {
        dots.clear();
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (i == 0) {
                auto from = pixel_to_vector<W>(nodes[i].x, node_y(nodes[i]));
                draw_vector<W>(dots, from, [this](auto& v, auto t) { return color(v, 0); });
            }
            else {
                auto from = pixel_to_vector<W>(nodes[i - 1].x, node_y(nodes[i - 1]));
                auto to = pixel_to_vector<W>(nodes[i].x, node_y(nodes[i]));
                draw_line<W>(dots, from, to, [this](auto& v, auto t) { return color(v, 0); }, 0, 1, false, true);
            }
        }
        trails.record(dots, age, alpha);
    }

    void pull(int leader) {
        nodes[leader].v = dir(speed);
        move(nodes[leader]);
        for (int i = leader - 1; i >= 0; --i) {
            drag(nodes[i + 1], nodes[i]);
        }
        for (size_t i = leader + 1; i < nodes.size(); ++i) {
            drag(nodes[i - 1], nodes[i]);
        }
    }

    void drag(Node& leader, Node& follower) {
        int dest = wrap(follower.x + follower.v, W);
        if (shortest_distance(dest, leader.x, W) > gap) {
            follower.v = leader.v;
            while (shortest_distance(follower.x, leader.x, W) > gap) {
                move(follower);
            }
        }
        else {
            move(follower);
        }
    }

    void move(Node& node) {
        node.x = wrap(node.x + node.v, W);
    }

    int dir(int speed) const {
        return speed < 0 ? -1 : 1;
    }

    Timeline timeline;

    static constexpr size_t MAX_PALETTES = 16;
    static constexpr size_t NUM_NODES = H_VIRT;
    StaticCircularBuffer <PaletteVariant, MAX_PALETTES> palettes;
    StaticCircularBuffer <float, MAX_PALETTES - 1> palette_boundaries;

    Vector palette_normal;
    std::array<Node, NUM_NODES> nodes;
    int speed;
    int gap;
    float alpha = 0.2;
    int wipe_duration = 160;
    uint32_t trail_length;
    Orientation orientation;
    Dots dots;

    DecayBuffer<W, 10000> trails;

    Pipeline<W,
        FilterReplicate<W>,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
};

///////////////////////////////////////////////////////////////////////////////

template <int W>
class MobiusGrid : public Effect {
public:
    MobiusGrid() :
        Effect(W),
        palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
        next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT),
        holeN(Z_AXIS),
        holeS(-Z_AXIS),
        filters(
            FilterHole<W>(holeN, 1.2f),
            FilterHole<W>(holeS, 1.2f),
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;

        timeline
            .add(0, MobiusWarp(params, num_rings, 160, true))
            .add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 400, ease_mid, true))
            .add(0, PeriodicTimer(120, [this](auto&) { wipe_palette(); }, true))
            .add(0, Mutation(num_rings, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true))
            .add(160, Mutation(num_lines, sin_wave(12.0f, 1.0f, 1.0f, 0.0f), 320, ease_mid, true));
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        dots.clear();

        float phase = fmodf(static_cast<float>(timeline.t), 120.0f) / 120.0f;

        draw_axis_rings(Z_AXIS, num_rings, phase);
        draw_longitudes(num_lines, phase);

        // Calculate Stabilizing Counter-Rotation
        Vector n_in = Z_AXIS;
        Vector n_trans = inv_stereo(mobius(stereo(n_in), params));
        Vector s_in = -Z_AXIS;
        Vector s_trans = inv_stereo(mobius(stereo(s_in), params));
        Vector mid = (n_trans + s_trans).normalize();
        Quaternion q = make_rotation(mid, Z_AXIS);

        // Apply counter-rotation
        for (auto& dot : dots) {
            dot.position = rotate(dot.position, q).normalize();
        }

        // Update hole origins to match the rotated geometry
        holeN = rotate(n_trans, q).normalize();
        holeS = rotate(s_trans, q).normalize();

        plot_dots<W>(dots, filters, canvas, 0, alpha);
    }

private:
    void wipe_palette() {
        next_palette = GenerativePalette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::FLAT);
        timeline.add(0, ColorWipe(palette, next_palette, 60, ease_mid));
    }

    void draw_axis_rings(const Vector& normal, float num, float phase) {
        const float log_min = -2.5f;
        const float log_max = 2.5f;
        const float range = log_max - log_min;
        int count = static_cast<int>(std::ceil(num));

        for (int i = 0; i < count; ++i) {
            float t = wrap((static_cast<float>(i) / num) + phase, 1.0f);
            float log_r = log_min + t * range;
            float r_val = expf(log_r);
            float radius = (4.0f / PI_F) * atanf(1.0f / r_val);

            Points points;
            sample_ring<W>(points, Quaternion(), normal, radius);

            Points transformed_points;
            for (const auto& p : points) {
                transformed_points.push_back(inv_stereo(mobius(stereo(p), params)));
            }

            float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
            rasterize<W>(dots, transformed_points, [&](const Vector&, float) {
                return palette.get(static_cast<float>(i) / num) * opacity;
                }, true);
        }
    }

    void draw_longitudes(float num, float phase) {
        int count = static_cast<int>(std::ceil(num));
        for (int i = 0; i < count; ++i) {
            float theta = (static_cast<float>(i) / num) * PI_F;
            Vector normal(cosf(theta), sinf(theta), 0.0f);

            Points points;
            sample_ring<W>(points, Quaternion(), normal, 1.0f);
            Points transformed_points;
            for (const auto& p : points) {
                transformed_points.push_back(inv_stereo(mobius(stereo(p), params)));
            }

            float opacity = std::clamp(num - static_cast<float>(i), 0.0f, 1.0f);
            rasterize<W>(dots, transformed_points, [&](const Vector&, float t_line) {
                // Approximate original Z to calculate log-gradient
                float original_idx = t_line * points.size();
                int idx1 = static_cast<int>(original_idx) % points.size();
                int idx2 = (idx1 + 1) % points.size();
                float f = original_idx - std::floor(original_idx);
                float z = points[idx1].k * (1.0f - f) + points[idx2].k * f;

                float R = sqrtf((1.0f + z) / (1.0f - z));
                float log_r = logf(R);
                const float log_min = -2.5f;
                const float log_max = 2.5f;
                float t = (log_r - log_min) / (log_max - log_min);

                return palette.get(wrap(t - phase, 1.0f)) * opacity;
                }, true);
        }
    }

    float alpha = 0.2f;
    float num_rings = 0;
    float num_lines = 0;
    GenerativePalette palette;
    GenerativePalette next_palette;
    MobiusParams params;
    Orientation orientation;
    Timeline timeline;
    Dots dots;

    Vector holeN;
    Vector holeS;

    Pipeline<W,
        FilterHoleRef<W>,
        FilterHoleRef<W>,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
};

////////////////////////////////////////////////////////////////////////////////

template <int W>
class Moire : public Effect {
public:
    Moire() :
        Effect(W),
        base_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
        base_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::BELL),
        int_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
        int_next_palette(GradientShape::CIRCULAR, HarmonyType::SPLIT_COMPLEMENTARY, BrightnessProfile::CUP),
        filters(
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;

        timeline
            .add(0, PeriodicTimer(80, [this](auto&) { color_wipe(); }))
            .add(0, Rotation<W>(orientation, Y_AXIS, 2 * PI_F, 300, ease_mid, true))
            .add(0, Transition(rotation, 2 * PI_F, 160, ease_mid, false, true)
                .then([this]() { rotation = 0.0f; }))
            .add(0, Mutation(amp, sin_wave(0.1f, 0.5f, 1.0f, 0.0f), 160, ease_mid, true)
            );
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        dots.clear();

        draw_layer([this](const Vector& p) { return inv_transform(p); }, base_palette);
        draw_layer([this](const Vector& p) { return transform(p); }, int_palette);

        plot_dots<W>(dots, filters, canvas, 0, alpha);
    }

private:
    void color_wipe() {
        base_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);
        int_next_palette = GenerativePalette(GradientShape::STRAIGHT, HarmonyType::TRIADIC, BrightnessProfile::ASCENDING);

        timeline.add(0, ColorWipe(base_palette, base_next_palette, 80, ease_mid));
        timeline.add(0, ColorWipe(int_palette, int_next_palette, 80, ease_mid));
    }

    Vector transform(const Vector& p) {
        Quaternion q_z = make_rotation(Z_AXIS, rotation);
        Quaternion q_x = make_rotation(X_AXIS, rotation);
        p = rotate(p, q_z);
        p = rotate(p, q_x);
        return p;
    }

    Vector inv_transform(const Vector& p) {
        Quaternion q_nx = make_rotation(-X_AXIS, rotation);
        Quaternion q_nz = make_rotation(-Z_AXIS, rotation);
        p = rotate(p, q_nx);
        p = rotate(p, q_nz);
        return p;
    }

    void draw_layer(TransformFn auto trans_fn, const GenerativePalette& pal) {
        int count = static_cast<int>(std::ceil(density));
        for (int i = 0; i <= count; ++i) {
            float t = static_cast<float>(i) / count;
            float r = t * 2.0f;

            Points points;
            sample_fn<W>(points, Quaternion(1, 0, 0, 0), Z_AXIS, r,
                sin_wave(-amp, amp, 4.0f, 0.0f));

            Points transformed_points;
            for (const auto& p : points) {
                transformed_points.push_back(trans_fn(p));
            }

            rasterize<W>(dots, transformed_points,
                [&](const Vector&, float) { return pal.get(t); }, true);
        }
    }

    float alpha = 0.1f;
    float density = 11.0f;
    float rotation = 0.0f;
    float amp = 0.0f;

    GenerativePalette base_palette;
    GenerativePalette base_next_palette;
    GenerativePalette int_palette;
    GenerativePalette int_next_palette;

    Orientation orientation;
    Timeline timeline;
    Dots dots;

    Pipeline<W,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
};

///////////////////////////////////////////////////////////////////////////////
// EXPERIMENTAL EFFECTS BELOW

template <int W>
class FlowField : public Effect {
public:
    FlowField() :
        Effect(W),
        palette(iceMelt),
        trails(k_trail_length),
        filters(
            FilterOrient<W>(orientation),
            FilterAntiAlias<W>()
        )
    {
        persist_pixels = false;
        noise_generator.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
        noise_generator.SetSeed(hs::rand_int(0, 65535));
        reset_particles();
    }

    bool show_bg() const override { return false; }

    void draw_frame() override {
        Canvas canvas(*this);
        timeline.step(canvas);
        dots.clear();

        for (Particle& p : particles) {
            // 1. Get acceleration from the 3D noise field (simulating 4D)
            Vector accel = get_noise_force(p.pos, timeline.t);

            // 2. Apply a gravity force pulling the particle to the center
            Vector gravity_force = p.pos * -k_gravity;
            accel = accel + gravity_force;

            // 3. Update velocity
            p.vel = p.vel + accel;

            // 4. Clamp velocity to max speed (manual implementation of clampLength)
            float speed = p.vel.length();
            if (speed > k_max_speed) {
                p.vel = (p.vel / speed) * k_max_speed;
            }

            // 5. Update position and re-normalize to keep it on the sphere
            p.pos = p.pos + p.vel;
            p.pos.normalize();

            // 6. Get color based on speed (ratio of current speed to max)
            float speed_ratio = speed / k_max_speed;
            Pixel color = palette.get(speed_ratio);

            // 7. Add the particle's "head" to the dot buffer
            dots.emplace_back(Dot(p.pos, color));
        }

        trails.record(dots, 0, k_alpha);

        trails.render(canvas, filters,
            [this](const Vector& v, float t) {
                return palette.get(1.0f - t) * (1.0f - t);
            }
        );
    }

private:

    struct Particle {
        Vector pos;
        Vector vel;

        Particle() : pos(random_vector()), vel(0, 0, 0) {}
    };

    static constexpr int k_num_particles = 100;
    static constexpr int k_trail_length = 8;
    static constexpr float k_noise_scale = 100;
    static constexpr float k_force_scale = 0.001;
    static constexpr float k_gravity = 0.001;
    static constexpr float k_max_speed = 0.1;
    static constexpr float k_alpha = 0.7;

    static constexpr float k_time_scale = 0.01;
    static constexpr int k_max_trail_dots = 1024;
    Timeline timeline;
    FastNoiseLite noise_generator;
    const ProceduralPalette& palette;
    std::array<Particle, k_num_particles> particles;

    Orientation orientation;
    DecayBuffer<W, k_max_trail_dots> trails;

    Pipeline<W,
        FilterOrient<W>,
        FilterAntiAlias<W>
    > filters;
    Dots dots;

    void reset_particles() {
        for (int i = 0; i < k_num_particles; ++i) {
            particles[i] = Particle();
        }
    }

    Vector get_noise_force(const Vector& pos, float t) {
        float t_scaled = t * k_time_scale;
        Vector n_pos = pos * k_noise_scale;

        float x_force = noise_generator.GetNoise(n_pos.i, n_pos.j, t_scaled);
        float y_force = noise_generator.GetNoise(n_pos.j, n_pos.k, t_scaled + 100.0f);
        float z_force = noise_generator.GetNoise(n_pos.k, n_pos.i, t_scaled + 200.0f);

        return Vector(x_force, y_force, z_force) * k_force_scale;
    }
};