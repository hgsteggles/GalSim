#ifndef STAR_HPP
#define STAR_HPP

#include <vector>
#include <array>

class Star {
public:
	Star();
	Star(int i, double m, double cla);
	int id;
	double mcur, mfin, age, dsun, d_gc, col_d, ext_v, ext_k, ext_21, sb_rad, fl_rad, r_ss, ang_ss, v_lsr, lbol, t_ms, ps_sw, logq, j21, xyz[3], gcoord[2];
};

class StarPop {
public:
	StarPop();
	void read(const std::string& filename);
	void print(const std::string& filename);
	void printFinal(const std::string& filename);
	void printSelected(const std::vector<int>& elements, const std::string& filename);
	void addStar(int i, double m, double a);
	std::vector<int> id, msx, kh;
	std::vector<double> age, dsun, d_gc, col_d, ext_v, ext_k, ext_21;
	std::vector<double> rad_sb, rad_fl, r_ss, ang_ss, v_lsr, lbol, mfin, mcur, t_ms, ps_sw, logq, j21, r0_ss, ne;
	std::array<std::vector<double>, 3> xyz;
	std::array<std::vector<double>, 3> gcoord;
	int size;
};

#endif
