#pragma once

template <int W, int H, int N>
class ParticleSystem {
public:

	class Vector {
	public:

		Vector() :
			x(0),
			y(0)
		{}

		Vector(float x, float y) :
			x(x),
			y(y)
		{}

		float x;
		float y;
	};

	class Particle {
	public:

		Particle() :
			x(0),
			y(0)
		{}

		Particle(float x, float y, const Vector& v) :
			x(x),
			y(y),
			v(v)
		{}

		float x;
		float y;
		Vector v;
	};

	class ExtentMask
	{
	public:

		ExtentMask() {
			memset(mask, 0, sizeof(mask));
		}

		ExtentMask operator &(const ExtentMask& m) {
			ExtentMask r;
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				r[i] = mask[i] & m[i];
			}
			return r;
		}

		ExtentMask operator |(const ExtentMask& m) {
			ExtentMask r;
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				r[i] = mask[i] | m[i];
			}
			return r;
		}

		ExtentMask& operator &=(const ExtentMask& m) {
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				mask[i] &= m[i];
			}
			return *this;
		}

		ExtentMask& operator |=(const ExtentMask& m) {
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				mask[i] |= m[i];
			}
			return *this;
		}

		ExtentMask& mark(int i) {
			mask[i / sizeof(int)] |= 0x1 << (i % sizeof(int));
			return *this;
		}

		ExtentMask& unmark(int i) {
			mask[i / sizeof(int)] ^= 0x1 << (i % sizeof(int));
			return *this;
		}

		unsigned int mask[((N - 1) / sizeof(int)) + 1];
	};

	ParticleSystem() :
		num_particles(0)
	{}

	int size() const { return num_particles; }

	Particle& operator[](int i) { return particles[i]; }

	void add_particle(const Particle& p) {
		if (num_particles >= N) {
			return;
		}
		particles[num_particles++] = p;
		mark_extent(p);
	}

	void step(int dt) {
		for (int i = 0; i < num_particles; ++i) {
			gravity(particles[i], dt);
			restrict(particles[i]);
			move(particles[i], dt);
		}
	}

private:

	void mark_extent(const Particle& p) {
		int left = floor(p.x);
		int right = ceil(p.x);
		for (int i = left; i != right; i = (i + 1) % W) {
			x_mask[num_particles - 1].mark(i);
		}

		int top = floor(p.y);
		int bottom = ceil(p.y);
		for (int i = top; i <= bottom; ++i) {
			y_mask[num_particles - 1].mark(i);
		}
	}

	void unmark_extent(const Particle& p) {
		int left = floor(p.x);
		int right = ceil(p.x);
		for (int i = left; i != right; i = (i + 1) % W) {
			x_mask[num_particles - 1].unmark(i);
		}

		int top = floor(p.y);
		int bottom = ceil(p.y);
		for (int i = top; i <= bottom; ++i) {
			y_mask[num_particles - 1].unmark(i);
		}
	}

	void restrict(Particle& p){
		if (p.y >= H - 1) {
			p.v.y = 0;
		}
	}

	void gravity(Particle& p, int dt) {
		p.v.y += G * dt;
	}

	void move(Particle& p, int dt) {
		unmark_extent(p);
		p.x = wrap(p.x + (dt * p.v.x), W);
		p.y = min(max(0, p.y + p.v.y), H - 1);
		mark_extent(p);
	}

	constexpr static const float G = 1.0 / 16;

	Particle particles[N];
	int num_particles;
	ExtentMask x_mask[W];
	ExtentMask y_mask[H];
};

