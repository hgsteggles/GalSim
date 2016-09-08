#ifndef FUNC_HPP
#define FUNC_HPP

#include <iostream> //std::cout
#include <fstream>
#include <sstream>
#include <iomanip> //std::setprecision
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <string>
#include <random>
#include <limits>
#include <algorithm>

#include "Recipes.hpp"

std::vector<double> deriv(const std::vector<double>& y);
double interpol(const std::vector<double>& f, const std::vector<double>& x, double subx);
double interpol2D(const std::vector<std::vector<double> >& f, int xc, int yc);
std::vector<double> cubeInterpol(const std::vector<std::vector<std::vector<double> > >& V, const std::vector<double>& xc, const std::vector<double>& yc, const std::vector<double>& zc);
std::vector<double> slice(const std::vector<double>& v, unsigned int imin, unsigned int imax);
std::vector<double> vinterpol(const std::vector<double>& f, const std::vector<double>& x, const std::vector<double>& subx);
double vmin(const std::vector<double>& f);
double vmax(const std::vector<double>& f);
double v3min(const std::vector<std::vector<std::vector<double> > >& vec);
double v3max(const std::vector<std::vector<std::vector<double> > >& vec);
std::vector<double> findgen(int n, double fac);
std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn);
double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2);
std::vector<double> vsplint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, const std::vector<double>& x2);
double int_tabulated(const std::vector<double>& x, const std::vector<double>& f);
double randomu(std::default_random_engine& engine);
double vtotal(const std::vector<double>& vec);
int round_int(double x);
double round_dp(double x, int dp);
std::vector<double> xhist(double binmin, double binmax, int nbins);
std::vector<int> yhist(const std::vector<double>& y, const std::vector<double>& x);
std::vector<int> random_loc(const std::vector<std::vector<int> >& vec2, const std::vector<std::vector<std::vector<double> > >& vec3,
 double val, std::default_random_engine& engine);
std::vector<int> random_loc2(const std::vector<std::vector<std::vector<double> > >& vec3, double val, std::default_random_engine& engine);
std::vector<int> rand_loc_gt(const std::vector<std::vector<std::vector<double> > >& vec3, double value, std::default_random_engine& engine);
void printScatter(const std::vector<double>& vec1, const std::vector<double>& vec2, std::string filename);
double planck(double wave, double temp);
double planck2(double freq, double temp);
double vmean(const std::vector<double>& vec);
double vstddev(const std::vector<double>& vec);
double random_kroupa(std::default_random_engine& engine);

//void readcols(std::string filename, std::string delimiter, std::vector<int>& v1, std::vector<string>& v2, std::vector<double>& v3, std::vector<double>& v4, std::vector<double>& v5, std::vector<double>& v6, std::vector<double>& v7, std::vector<double>& v8, std::vector<double>& v9);

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// Initialize original index locations.
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// Sort indexes based on comparing values in v.
	std::sort(idx.begin(), idx.end(),
			  [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

template <typename type1, typename UnaryPredicate>
std::vector<int> indicesWhereTrue(const std::vector<type1>& oldvec, UnaryPredicate&& f) {
	std::vector<int> vec;
	for (unsigned int i = 0; i < oldvec.size(); i++) {
		if ( f(oldvec[i]) )
			vec.push_back(i);
	}
	return vec;
}

template <class type1>
std::vector<type1> copyElements(const std::vector<type1>& oldvec, const std::vector<int>& coords) {
	std::vector<type1> newvec;
	for (unsigned int i = 0; i < oldvec.size(); i++) {
		for (unsigned int j = 0; j < coords.size(); j++) {
			if (coords[j] == i) {
				newvec.push_back(oldvec[i]);
			}
		}
	}
	return newvec;
}

template <class type1>
int nmatches(const std::vector<type1>& v1, type1 value) {
	int result = 0;
	int n = v1.size();
	for (int i = 0; i < n; i++) {
		if (v1[i] == value)
			result++;
	}
	return result;
}

template <class t1>
void readcol(std::string filename, std::string delimiter, int position, std::vector<t1>& variable) {
	variable.clear();
	size_t pos = 0;
	std::string token;
	std::ifstream ifile(filename.c_str());
	if (!ifile)
		std::cerr << "ERROR: unable to open " << filename << " for reading." << std::endl;
	std::string line;
	getline(ifile, line);
	while (getline(ifile, line)) {
		for (int i = 1; i <= position; i++) {
			//cerr << i << endl;
			pos = line.find(delimiter);
			token = line.substr(0, pos);
			if (i == position) {
				std::stringstream convert(token);
				t1 result;
				if (!(convert >> result))
					result = 0;
				variable.push_back(result);
				break;
			}
			line.erase(0, pos + delimiter.length());
		}
	}
	ifile.close();
}

template <class t1, class t2, class t3, class t4, class t5, class t6, class t7, class t8, class t9>
void readcols(std::string filename, std::string delimiter, std::vector<t1>& v1, std::vector<t2>& v2, std::vector<t3>& v3, std::vector<t4>& v4, std::vector<t5>& v5, std::vector<t6>& v6, std::vector<t7>& v7, std::vector<t8>& v8, std::vector<t9>& v9) {
	v1.clear();
	v2.clear();
	v3.clear();
	v4.clear();
	v5.clear();
	v6.clear();
	v7.clear();
	v8.clear();
	v9.clear();
	t1 result1;
	t2 result2;
	t3 result3;
	t4 result4;
	t5 result5;
	t6 result6;
	t7 result7;
	t8 result8;
	t9 result9;
	size_t pos = 0;
	std::string token;
	std::ifstream ifile(filename.c_str());
	if (!ifile)
		std::cerr << "ERROR: unable to open " << filename << " for reading." << std::endl;
	std::string line;
	while (getline(ifile, line)) {
		bool errorline = false;
		int i;
		for (i = 1; i <= 9; i++) {
			pos = line.find(delimiter);
			token = line.substr(0, pos);
			std::stringstream convert(token);
			if (i == 1)
				errorline = !(convert >> result1);
			else if (i == 2)
				errorline = !(convert >> result2);
			else if (i == 3)
				errorline = !(convert >> result3);
			else if (i == 4)
				errorline = !(convert >> result4);
			else if (i == 5)
				errorline = !(convert >> result5);
			else if (i == 6)
				errorline = !(convert >> result6);
			else if (i == 7)
				errorline = !(convert >> result7);
			else if (i == 8)
				errorline = !(convert >> result8);
			else if (i == 9)
				errorline = !(convert >> result9);
			line.erase(0, pos + delimiter.length());
			if (errorline)
				break;
		}
		if (errorline == false) {
			v1.push_back(result1);
			v2.push_back(result2);
			v3.push_back(result3);
			v4.push_back(result4);
			v5.push_back(result5);
			v6.push_back(result6);
			v7.push_back(result7);
			v8.push_back(result8);
			v9.push_back(result9);
		}
	}
	ifile.close();
}

#endif
