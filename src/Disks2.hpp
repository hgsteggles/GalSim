//
// Created by harry on 08/08/16.
//

#ifndef DISKS2_HPP
#define DISKS2_HPP

#include "GalaxyParameters.hpp"

#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <iomanip> //std::setprecision
#include <stdio.h>
#include <stdlib.h>
#include <vector> //std::std::vector

class ThinDisk {
public:
	ThinDisk(const GalaxyParameters& gparams);
	GalaxyParameters params;
	std::vector<std::vector<std::vector<double> > > density;
};

class ThickDisk {
public:
	ThickDisk(const GalaxyParameters& gparams);
	GalaxyParameters params;
	std::vector<std::vector<std::vector<double> > > density;
};

#endif //DISKS2_HPP
