#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

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
            rasterize<W>(dots, points, [&](auto& v, auto t) { return VignettePalette(ring.palette).get(0); }, true);
            plot_dots<W>(dots, ring.filters, canvas, t, alpha * opacity);
            });
        ring.orientation.collapse();

        ring.filters.trail(canvas,
            [&](float x, float y, float t) { return VignettePalette(ring.palette).get(t); },
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
