
const float pi = 3.1415926535897f;
const float tau = 2 * pi;
const float radians = pi / 180;

inline float mod_tau(float n) {
	if (n > tau) return n - tau;
	else if (n < 0) return n + tau;
	return n;
}

inline float modf(float x, int m) {
	return (static_cast<int>(x) % m) + (x - static_cast<int>(x));
}

template <uint8_t W, uint8_t H>
class Projection
{
public:

	struct Point
	{
		Point(int x, int y) :
			x(x),
			y(y),
			lambda(x* tau / W - pi),
			phi((H - y)* pi / H - pi / 2)
		{}

		Point(uint8_t x, uint8_t y) :
			x(x),
			y(y),
			lambda(x* tau / W - pi),
			phi((H - y)* pi / H - pi / 2)
		{}

		Point(float x, float y) :
			x(x),
			y(y),
			lambda(x* tau / W - pi),
			phi((H - y)* pi / H - pi / 2)
		{}

		Point(const Point& p) :
			x(p.x),
			y(p.y),
			lambda(p.lambda),
			phi(p.phi)
		{}

		int xi() { return static_cast<int>(x + 0.5f) % W; }
		int yi() { return static_cast<int>(y + 0.5f); }

		float x;
		float y;
		float lambda;
		float phi;
	};

	Projection()
	{}

	Point project(uint8_t bx, uint8_t by) const {
		Point p(bx, by);
		return project(p);
	}

	Point project(const Point& src) const {
		Point p(src);

		// rotate lambda
		p.lambda += delta_lambda;
		if (p.lambda > pi) {
			p.lambda -= tau;
		}
		else if (p.lambda < -pi) {
			p.lambda += tau;
		}

		// convert to cartesian x, y, z
		float cos_p = cosf(p.phi);
		float x = cosf(p.lambda) * cos_p;
		float y = sinf(p.lambda) * cos_p;
		float z = sinf(p.phi);

		// rotate phi gamma
		float k = z * cos_dp + x * sin_dp;
		p.lambda = atan2f(y * cos_dg - k * sin_dg, x * cos_dp - z * sin_dp);
		p.phi = asinf(k * cos_dg + y * sin_dg);

		// convert to equirectangular x, y
		p.x = modf((p.lambda + pi) * W / tau, W);
		p.y = H - ((p.phi + pi / 2) * H / pi);

		return p;
	}

	Projection& rotate(uint16_t dl, uint16_t dp, uint16_t dg) {
		delta_lambda = mod_tau(delta_lambda + ((dl % 360) * radians));
		delta_phi = mod_tau(delta_phi + ((dp % 360) * radians));
		delta_gamma = mod_tau(delta_gamma + ((dg % 360) * radians));

		cos_dp = cosf(delta_phi);
		sin_dp = sinf(delta_phi);
		cos_dg = cosf(delta_gamma);
		sin_dg = sinf(delta_gamma);

		return *this;
	}

	Projection& reset() {
		delta_lambda = 0;
		delta_phi = 0;
		delta_gamma = 0;

		cos_dp = cosf(delta_phi);
		sin_dp = sinf(delta_phi);
		cos_dg = cosf(delta_gamma);
		sin_dg = sinf(delta_gamma);

		return *this;
	}

private:

	float delta_lambda = 0;
	float delta_phi = 0;
	float delta_gamma = 0;
	float cos_dp = cosf(delta_phi);
	float sin_dp = sinf(delta_phi);
	float cos_dg = cosf(delta_gamma);
	float sin_dg = sinf(delta_gamma);
};

