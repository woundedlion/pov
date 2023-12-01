//#include <fixmath.h>
#include <SPI.h>
#include <FastLED.h>
#include "led.h"
//#include "effects.h"

//#include "images.h"

namespace {
	POVDisplay<NUM_PIXELS, RPM> *pov;
}

void setup() {
  Serial.begin(9600);
  delay(1000);
	pov = new POVDisplay<NUM_PIXELS, RPM>();

}

void loop() {
	/*
	Serial.printf("%d\n", millis());
	Dot<40, 20> d(30, 10, CHSV(0, 0, 255));
	Serial.printf("d_x_y: (%f, %f)\n", d.x, d.y);
	Serial.printf("d_lambda_phi: (%f, %f)\n", d.lambda, d.phi);
	Serial.printf("d_cartesian: %f, %f, %f\n", d.x_cartesian(), d.y_cartesian(), d.z_cartesian());
	Quaternion q(1, 0, 0, 0);
	Serial.printf("q: (%f, %f, %f, %f)\n", q.r, q.v.i, q.v.j, q.v.k);
	Quaternion p(0, d.x_cartesian(), d.y_cartesian(), d.z_cartesian());
	Serial.printf("p: (%f, %f, %f, %f)\n", p.r, p.v.i, p.v.j, p.v.k);
	auto r = q.inverse() * p * q;
	Serial.printf("r: (%f, %f, %f, %f)\n", r.r, r.v.i, r.v.j, r.v.k);
	Dot<40, 20> d2(r.v.i, r.v.j, r.v.k, d.color);
	Serial.printf("d2_cartesian: %f, %f, %f\n", d2.x_cartesian(), d2.y_cartesian(), d2.z_cartesian());
	Serial.printf("d2_lambda_phi: (%f, %f)\n", d2.lambda, d2.phi);
	Serial.printf("d2_x_y: (%f, %f)\n", d2.x, d2.y);
	delay(1000);
	*/
    	pov->show<Test<96> >(9999);

//	pov->show<TheMatrix<40, 135> >(60);
/*	pov.show<ChainWiggle<96> >(150);
	pov.show<RingRotate<96> >(120);
	pov.show<RingShower<96> >(90);
	pov.show<Curves<96> >(120);
	pov.show<RingTwist<96> >(150);
	pov.show<Kaleidoscope<96> >(120);
	pov.show<StarsFade<40> >(60);
	pov.show<RingTrails<96> >(90);
	pov.show<DotTrails<96> >(120);
	pov.show<Burnout<40, 0, 5> >(52);
	pov.show<WaveTrails<96> >(60);
	pov.show<Spinner<48> >(50);

	pov.show<Spiral<4, 0> >(10);
	pov.show<Spiral<8, 0> >(10);
	pov.show<Spiral<16, 0> >(10);
	pov.show<Spiral<48, 0> >(30);
	pov.show<Spiral<48, 1> >(20);

	*/

	//	pov.show<Fire<40, 150, 120> >(45);


	/*
	pov.show<TheMatrix<40, 135> >(15);
	pov.show<DotTrails<96> >(30);
	pov.show<ChainWiggle<96> >(30);
	pov.show<RingRotate<96> >(30);
	pov.show<RingShower<96> >(15);
	pov.show<Curves<96> >(30);
	pov.show<RingTwist<96> >(30);
	pov.show<Kaleidoscope<96> >(30);
	pov.show<StarsFade<40> >(15);
	pov.show<RingTrails<96> >(15);
//	pov.show<Burnout<40, 0, 5> >(52);
	pov.show<WaveTrails<96> >(60);
	pov.show<Spinner<48> >(50);
	pov.show<Fire<40, 150, 120> >(15);

	pov.show<Spiral<4, 0> >(2);
	pov.show<Spiral<8, 0> >(2);
	pov.show<Spiral<16, 0> >(2);
	pov.show<Spiral<48, 0> >(5);
	pov.show<Spiral<48, 1> >(5);
	pov.show<Spiral<48, 2> >(2);
	*/
}


