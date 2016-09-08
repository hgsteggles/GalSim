/*
 * GalaxyParameters.hpp
 *
 *  Created on: 5 Dec 2014
 *      Author: harry
 */

#ifndef GALAXYPARAMETERS_HPP_
#define GALAXYPARAMETERS_HPP_

#include <string>
#include <iostream>
#include <array>

using ArmArray = std::array<double, 4>;

class GalaxyParameters {
public:
	//void printInfo();

	double sfr = 3.0;
	double total_time = 1.0e6;

	double radius = 0;
	double height = 0;
	double resolution = 0;
	std::array<double, 3> solar_position = std::array<double, 3>{0.0, 0.0, 0.0};
	double A_inner = 0;
	double Aa = 0;
	double A1 = 0;
	double A2 = 0;
	double H1 = 0;
	double H2 = 0;
	double h1 = 0;
	double n1 = 0;
	double n2 = 0;
	double rmin = 0;
	std::array<double, 4> wa_wj = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> na_fj = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> ha_hj = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> sp_joins = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> sp_rmin = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> sp_thmin = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> sp_asp = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::array<double, 4> sp_rmax = std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
	std::string starpopfile = "";
	std::string densityfile = "";
	std::string outputDirectory = "tmp";
};

/*
void GalaxyParameters::printInfo() {
	std::cout << "========= GalaxyParameters =========\n";
	std::cout << "radius = " << radius << '\n';
	std::cout << "height = " << radius << '\n';
	for (int i = 0; i < 3; ++i)
		std::cout << "solar_position[" << i << "] = " << solar_position[i] << '\n';
	std::cout << "A_inner = " << A_inner << '\n';
	std::cout << "A_a = " << Aa << '\n';
	for (int i = 0; i < 5; ++i)
		std::cout << "wa_wj[" << i << "] = " << wa_wj[i] << '\n';
	for (int i = 0; i < 5; ++i)
		std::cout << "na_fj[" << i << "] = " << na_fj[i] << '\n';
	for (int i = 0; i < 5; ++i)
		std::cout << "ha_hj[" << i << "] = " << ha_hj[i] << '\n';
	for (int i = 0; i < 4; ++i)
		std::cout << "arm_joinpoints[" << i << "] = " << sp_joins[i] << '\n';
	std::cout << "rmin_sparms = " << rmin << '\n';
	std::cout << "starpopfile = " << starpopfile << '\n';
	std::cout << "densityfile = " << densityfile << '\n';
}
*/


#endif /* GALAXYPARAMETERS_HPP_ */
