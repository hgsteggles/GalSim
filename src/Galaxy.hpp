#ifndef GALAXY
#define GALAXY

#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <iomanip> //std::setprecision
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <random>
#include <algorithm>
#include <memory>

#include "GalaxyParameters.hpp"
#include "Disks2.hpp"
#include "SpiralArms2.hpp"
#include "Star.hpp"
#include "func.hpp"
#include "HosTracks.hpp"
#include "Recipes.hpp"

class Galaxy {
public:
	Galaxy(const GalaxyParameters& gal_params);
	~Galaxy();

	void readDensity(const std::string& filename);
	void calcDensity(const std::string& filename);

	void print(const std::string& filename);
	void print_xz(const std::string& filename);
	void print3D(const std::string& filename);

	void calcExtinction(const std::string& filename);
	void populate();
	void putStars();
	void putStars2();
	void currStarMass();
	void calcFluxes();

	double pix2kpc(double pixel, int dim) const;
	int kpc2pix(double x, int dim) const;

	GalaxyParameters params;

	std::default_random_engine rng;

	std::unique_ptr<ThinDisk> thinDisk;
	std::unique_ptr<ThickDisk> thickDisk;
	std::unique_ptr<SpiralArms> spiralArms;

	std::vector<double> gal_x;
	std::vector<double> gal_y;
	std::vector<double> gal_z;

	std::vector<std::vector<std::vector<double> > > density;
	std::vector<std::vector<std::vector<double> > > normdensity;
	//std::vector<Star> stars;
	StarPop stars;

	int nocells[3];
	double extinction_a, extinction_b;
	double logq_mms;
};

#endif
