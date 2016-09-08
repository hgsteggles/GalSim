/*
 * Statistics.hpp
 *
 *  Created on: 11 Dec 2014
 *      Author: harry
 */

#ifndef STATISTICS_HPP_
#define STATISTICS_HPP_

#include <cmath>
#include <algorithm>

namespace Statistics {

template<typename T>
double mean(const T& vec) {
	double sum = std::accumulate(std::begin(vec), std::end(vec), 0.0);
	return sum / (double)vec.size();
}

template<typename T>
double stddev(const T& vec) {
	double m = Statistics::mean(vec);

	double accum = 0.0;
	std::for_each (std::begin(vec), std::end(vec), [&](const double d) {
	    accum += (d - m) * (d - m);
	});

	double stdev = std::sqrt(accum / (double)(vec.size()-1));
}

}



#endif /* STATISTICS_HPP_ */
