#ifndef SPIRALARMS2_HPP
#define SPIRALARMS2_HPP

#include <vector>

#include "func.hpp"
#include "Constants.hpp"
#include "GalaxyParameters.hpp"

class SpiralArms {
	public:
		SpiralArms(const GalaxyParameters& gparams);
		void calcDensity();
		std::vector<std::vector<std::vector<double> > > density;
		int nocells[3];
	private:
		GalaxyParameters params;
};

#endif
