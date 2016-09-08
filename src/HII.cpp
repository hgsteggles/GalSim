//
// Created by harry on 11/05/16.
//

#include <cmath>
#include "HII.hpp"
#include "Constants.hpp"

double HII::stromgrenRadius(double ne, double logq, double T) {
	double alphaH = 2.59e-13; // cm3.s-1
	alphaH *= std::pow(T / 10000.0, -0.7);
	return std::pow(3.0 * std::pow(10.0, logq) / (4.0 * CST::PI * ne * ne * alphaH), 1.0 / 3.0);
}

double HII::spitzerRadius(double stromRad, double soundSpeed, double time) {
	double ts = stromRad / soundSpeed;
	return stromRad * std::pow(1.0 + (7.0 / 4.0) * time / ts, 4.0 / 7.0);
}