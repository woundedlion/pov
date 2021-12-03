#pragma once

template<typename V>
V operator*(const V& v, float a) {
	return V(v.x * a, v.y * a);
}

template<typename V>
V operator*(float a, const V& v) {
	return v * a;
}

template<typename V>
V& operator*=(V& v, float a) {
	v.x *= a;
	v.y *= a;
	return v;
}

template<typename V>
V operator/(const V& v, float a) {
	return Vector(v.x / a, v.y / a);
}

template<typename V>
V operator/(float a, const V& v) {
	return v / a;
}

template<typename V>
V& operator/=(V& v, float a) {
	v.x /= a;
	v.y /= a;
	return v;
}

template<typename V>
V operator+(const V& v1, const V& v2) {
	return V(v1.x + v2.x, v1.y + v2.y);
}

template<typename V>
V& operator+=(V& v1, const V& v2) {
	v1.x += v2.x;
	v1.y += v2.y;
	return v1;
}

template<typename V>
V operator-(const V& v1, const V& v2) {
	return V(v1.x - v2.x, v1.y - v2.y);
}

template<typename V>
V& operator-=(V& v1, const V& v2) {
	v1.x -= v2.x;
	v1.y -= v2.y;
	return v1;
}

////////////////////////////////////////////////////////////////////////////////////////

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

		bool get(int i) const {
			return (mask[i / sizeof(int)] >> (i % sizeof(int))) & 0x1;
		}

		ExtentMask operator &(const ExtentMask& m) const {
			ExtentMask r;
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				r.mask[i] = mask[i] & m.mask[i];
			}
			return r;
		}

		ExtentMask operator |(const ExtentMask& m) const{
			ExtentMask r;
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				r.mask[i] = mask[i] | m.mask[i];
			}
			return r;
		}

		ExtentMask& operator &=(const ExtentMask& m) {
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				mask[i] &= m.mask[i];
			}
			return *this;
		}

		ExtentMask& operator |=(const ExtentMask& m) {
			for (size_t i = 0; i < sizeof(mask) / sizeof(int); ++i) {
				mask[i] |= m.mask[i];
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

		template <typename F>
		const ExtentMask& each(F&& f) const {
			for (int i = 0; i < N; ++i) {
				if (get(i)) {
					f(i);
				}
			}
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
		mark_extent(num_particles - 1);
	}

	void step(float dt) {
		float weights[N];
		ExtentMask neighbors[N];
		for (int i = 0; i < num_particles; ++i) {
			neighbors[i] = x_overlap(i) & y_overlap(i);
			neighbors[i].each([&](int j) {
				weights[i] = weight(i, j);
			});
		}
		for (int i = 0; i < num_particles; ++i) {
			gravity(i, dt);
			viscosity(i, neighbors[i], dt);
			pressure(i, neighbors[i], weights, dt);
			restrict(i);
			move(i, dt);
		}
	}

private:

	void mark_extent(int i) {
		Particle& p = particles[i];
		int left = floor(p.x);
		int right = ceil(p.x);
		for (int i = left; i != right; i = (i + 1) % W) {
			x_mask[i].mark(i);
		}

		int top = floor(p.y);
		int bottom = ceil(p.y);
		for (int i = top; i <= bottom; ++i) {
			y_mask[i].mark(i);
		}
	}

	void unmark_extent(int i) {
		Particle& p = particles[i];
		int left = floor(p.x);
		int right = ceil(p.x);
		for (int i = left; i != right; i = (i + 1) % W) {
			x_mask[i].unmark(i);
		}

		int top = floor(p.y);
		int bottom = ceil(p.y);
		for (int i = top; i <= bottom; ++i) {
			y_mask[i].unmark(i);
		}
	}

	void restrict(int i){
		Particle & p = particles[i];
		if (p.y >= H - 1) {
			p.v.y = 0;
		}
	}

	ExtentMask x_overlap(int i) {
		Particle& p = particles[i];
		ExtentMask overlap;
		for (int x = static_cast<int>(floor(p.x));
			x < static_cast<int>(ceil(p.x)); ++x)
		{
			overlap |= x_mask[x];
		}
		return overlap;
	}

	ExtentMask y_overlap(int i) {
		Particle& p = particles[i];
		ExtentMask overlap;
		for (int y = static_cast<int>(floor(p.y));
			y < static_cast<int>(ceil(p.y)); ++y)
		{
			overlap |= y_mask[y];
		}
		return overlap;
	}

	float distance(const Particle& p, const Particle& q) {
		return sqrt(pow(q.x - p.x, 2) + pow(q.y - p.y, 2));
	}

	float weight(int i, int j) {
		return max(0, 1.0 - distance(particles[i], particles[j]) / 1.0); // TODO: variable diameter, other smoothng kernels?
	}

	Vector normal(int i, int j) {
		Vector v = particles[j].v - particles[i].v;
		float length = sqrt(pow(v.x, 2) + pow(v.y, 2));
		v /= length;
		return v;
	}

	void pressure(int i, const ExtentMask& neighbors, float *weights, float dt) {
		neighbors.each([&](int j) {	
			particles[i].v += this->normal(i, j) * PRESSURE * dt * (weights[i] + weights[j]);
		});
	}

	void viscosity(int i, const ExtentMask& neighbors, float dt) {
		neighbors.each([&](int j) {
			particles[i].v += (particles[j].v - particles[i].v) * VISCOSITY * weight(i, j) * dt;
		});
	}

	void gravity(int i, float dt) {
		Particle& p = particles[i];
		p.v.y += G * dt;
	}

	void move(int i, float dt) {
		Particle& p = particles[i];
		unmark_extent(i);
		p.x = wrap(p.x + (dt * p.v.x), W);
		p.y = min(max(0, p.y + p.v.y), H - 1);
		mark_extent(i);
	}

	constexpr static const float G = 1.0 / 16;
	constexpr static const float PRESSURE = 0.2;
	constexpr static const float VISCOSITY = 0.1;

	Particle particles[N];
	int num_particles;
	ExtentMask x_mask[W];
	ExtentMask y_mask[H];
};

