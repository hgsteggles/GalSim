/*
 * Random.cpp
 *
 *  Created on: 11 Dec 2014
 *      Author: harry
 */

#include "Random.hpp"

std::default_random_engine & Random::getRandomEngine() {
	static std::default_random_engine u{};
	return u;
}

void Random::randomise() {
	static std::random_device rd{};
	getRandomEngine().seed( rd() );
}

int Random::randomInteger( int from, int upto ) {
	static std::uniform_int_distribution<> d{};
	using parm_t = decltype(d)::param_type;
	return d( getRandomEngine(), parm_t{from, upto} );
}

double Random::randomDouble( double from, double upto ) {
	static std::uniform_real_distribution<> d{};
	using parm_t = decltype(d)::param_type;
	return d( getRandomEngine(), parm_t{from, upto} );
}

double Random::randomNormal( double mean, double stdev ) {
	static std::normal_distribution<> d{};
	using parm_t = decltype(d)::param_type;
	return d(getRandomEngine(), parm_t{mean, stdev});
}

double Random::randomPowerlaw( double index, double min, double max ) {
	double p1 = index + 1.0;
	double lo = min;
	double hi = max;

	if (hi < lo)
		std::swap(hi, lo);

	double r = Random::randomDouble(0.0, 1.0);

	if (index != -1.0) {
		double norm = 1.0/(std::pow(hi, p1) - std::pow(lo, p1));
		double expo = std::log10(r/norm + std::pow(lo, p1)) / p1;
		return std::pow(10.0, expo);
	}
	else {
		double norm = 1.0/(std::log(hi) - std::log(lo));
		return std::exp((r/norm) + std::log(lo));
	}
}


