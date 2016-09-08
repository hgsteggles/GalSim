#ifndef HOS_TRACKS_HPP
#define HOS_TRACKS_HPP

#include <string>
#include <vector>
#include "func.hpp"

class HosData{
	public:
		HosData();
		std::vector<double> mass;
		std::vector<double> rstar;
		std::vector<double> rphot;
		std::vector<double> lstar;
		std::vector<double> ltot;
		std::vector<double> time;
};

class HosTrack{
	public:
		HosTrack();
		HosTrack(double md, const std::vector<double>& m, const std::vector<double>& r, const std::vector<double>& rph, const std::vector<double>& l, const std::vector<double>& lt, const std::vector<double>& t);
		double mdot;
		HosData data;
};

class HosTracks{
	public:
		HosTracks(std::string filename);
		std::vector<HosTrack> tracks;
};

#endif
