#pragma once

int wrap(int x, int m) {
	return (x >= 0 ?
		x % m :
		((x % m) + m) % m);
}

int wrap(float x, int m) {
	return x >= 0 ?
		fmod(x, m) :
		fmod(fmod(x, m) + m, m);
}

int wrap(float x, float m) {
	return x >= 0 ?
		fmod(x, m) :
		fmod(fmod(x, m) + m, m);
}

class Oscillator
{
public:

	Oscillator(int min, int max, int end_delay = 0) :
		t_(min),
		min_(min),
		max_(max),
		end_delay_(end_delay),
		end_count_(0),
		dir_(1)
	{}

	int get() {
		int r = t_;
		if ((t_ == max_ && dir_ == 1) || (t_ == min_ && dir_ == -1)) {
			if (end_count_++ >= end_delay_) {
				dir_ = dir_ * -1;
				end_count_ = 0;
			}
			else {
				return r;
			}
		}
		t_ += dir_;
		return r;
	}

	int dir() { return dir_; }
	void set_delay(int d) { end_delay_ = d; }

private:

	int t_;
	int min_;
	int max_;
	int end_delay_;
	int end_count_;
	int dir_;
};

