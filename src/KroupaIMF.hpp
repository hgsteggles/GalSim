/*
 * KroupaIMF.hpp
 *
 *  Created on: 11 Dec 2014
 *      Author: harry
 */

#ifndef KROUPAIMF_HPP_
#define KROUPAIMF_HPP_

#include <array>

class KroupaIMF {
public:
	KroupaIMF();
	double randomMass() const;
private:
	std::array<double, 4> p;
	std::array<std::pair<double, double>, 4> massRanges;
	std::array<double, 4> rangeCDF;
};


#endif /* KROUPAIMF_HPP_ */
