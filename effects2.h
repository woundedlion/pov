#pragma once

#include "effects_infra.h"

DEFINE_GRADIENT_PALETTE(lucid_dream_p) {
  0, 6, 4, 47,  
  128, 162, 84, 84,
  255, 252, 114, 0 
};

template <int W>
class Test : public Effect {
public:
  Test() :
    Effect(W),
    decay(5, lucid_dream_p)
  {  
  antialias.chain(decay);
    sprites.emplace_back(Sprite<W, H>());
   for (int i = 0; i < 1; ++i) {
      sprites.back().build(i, 10, CHSV(0, 0, 255));
    }
  }

  bool show_bg() const { return false; }

  void draw_frame() {
    Canvas c(*this);
    decay.age(c);
    sprites.front()
      .rotate(Vector(0, 1, 0), 18, 1);
    paint(c, antialias, sprites.front());
  }

private:
  typedef std::vector<Sprite<W, H>> Sprites;

  FilterDecay<W, H> decay;
  FilterAntiAlias<W, H> antialias;
  Sprites sprites;
};