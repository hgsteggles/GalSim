/*
 * KroupaIMF.cpp
 *
 *  Created on: 11 Dec 2014
 *      Author: harry
 */

#include "KroupaIMF.hpp"
#include "Random.hpp"

#include <cmath>
#include <iostream>

KroupaIMF::KroupaIMF() {
	p[0] = -0.3;
	p[1] = -1.8;
	p[2] = -2.7;
	p[3] = -2.35;

	massRanges[0] = std::make_pair(0.01, 0.08);
	massRanges[1] = std::make_pair(0.08, 0.5);
	massRanges[2] = std::make_pair(0.5, 1.0);
	massRanges[3] = std::make_pair(1.0, 120.0);

	double sum = 0.0;

	std::array<double, 4> probs;
	for (int i = 0; i < 4; i++) {
		if (p[i] != -1.0)
			probs[i] = (std::pow(massRanges[i].second, p[i]+1) - std::pow(massRanges[i].first, p[i]+1))/(p[i]+1);
		else
			probs[i] = std::log(massRanges[i].second/massRanges[i].first);
		for (int j = i - 1; j >= 0; --j)
			probs[i] *= std::pow(massRanges[j].second, p[j] - p[j+1]);
		sum += probs[i];
	}

	for (int i = 0; i < 4; i++) {
		rangeCDF[i] = probs[i]/sum;
		if (i != 0)
			rangeCDF[i] += rangeCDF[i-1];
	}
}

double KroupaIMF::randomMass() const {
	double r = Random::randomDouble(0.0, 1.0);
	int seg = 3;

	for (int i = 0; i < 4; i++) {
		if (r < rangeCDF[i])
			return Random::randomPowerlaw(p[i], massRanges[i].first, massRanges[i].second);
	}

	throw std::runtime_error("KroupaIMF::randomMass: rangeCDF incorrectly set up.");
}


