#pragma once
#include <functional>
#include <memory>
#include <vector>
#include <map>
#include "../effects_engine.h"

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
