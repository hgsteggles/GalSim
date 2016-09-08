/*
 * Random.hpp
 *
 *  Created on: 17 Nov 2014
 *      Author: harry
 */

#ifndef RANDOM_HPP_
#define RANDOM_HPP_

#include <random>

namespace Random {

std::default_random_engine & getRandomEngine();

void randomise();

int randomInteger( int from, int thru );

double randomDouble( double from, double upto );

double randomNormal( double mean, double stdev );

double randomPowerlaw( double index, double min, double max);

}



#endif /* RANDOM_HPP_ */
