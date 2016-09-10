#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Star.hpp"
/*
mcur    -- Msun
mfin    -- Msun
clmass  -- Msun
age     -- kyr
dsun    -- kpc
d_gc    -- kpc
col_d   -- cm^-3 kpc
ext_v   -- ?
ext_k   -- ?
ext_21  -- ?
sb_rad  -- ?
fl_rad  -- ?
r_ss    -- ?
v_lsr   -- ?
lbol    -- lsun
t_ms    -- ?
ps_sw   -- ?
logq    -- ?
j21     -- Jy
*/

StarPop::StarPop()
: size(0)
{}

void StarPop::read(const std::string& filename) {
	std::ifstream ifile(filename);
	if (!ifile)
		throw std::runtime_error("StarPop::read: unable to open \'" + filename + "\' for reading.");
	bool errorline = false;

	std::string dummy;
	ifile >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

	while (!errorline) {
		double d_mfin, d_age, d_dsun, d_d_gc, d_col_d, d_ext_v, d_ext_k, d_ext_21;
		double d_x, d_y, d_z, d_l, d_b;
		errorline = !(ifile >> d_mfin >> d_age >> d_dsun >> d_d_gc >> d_col_d >> d_ext_v >> d_ext_k >> d_ext_21 >> d_x >> d_y >> d_z >> d_l >> d_b);
		if (!errorline) {
			id.push_back(size);
			mcur.push_back(0);
			mfin.push_back(d_mfin);
			age.push_back(d_age);
			dsun.push_back(d_dsun);
			d_gc.push_back(d_d_gc);
			col_d.push_back(d_col_d);
			ext_v.push_back(d_ext_v);
			ext_k.push_back(d_ext_k);
			ext_21.push_back(d_ext_21);
			rad_sb.push_back(0);
			rad_fl.push_back(0);
			r_ss.push_back(0);
			ang_ss.push_back(0);
			v_lsr.push_back(0);
			lbol.push_back(0);
			t_ms.push_back(0);
			ps_sw.push_back(0);
			logq.push_back(0);
			j21.push_back(0);
			r0_ss.push_back(0);
			ne.push_back(0);
			msx.push_back(0);
			kh.push_back(0);
			xyz[0].push_back(d_x);
			xyz[1].push_back(d_y);
			xyz[2].push_back(d_z);
			gcoord[0].push_back(d_l);
			gcoord[1].push_back(d_b);
			size++;
		}
	}
	ifile.close();
}

void StarPop::addStar(int i, double m, double cla) {
	id.push_back(i);
	mcur.push_back(0);
	mfin.push_back(m);
	age.push_back(cla);
	dsun.push_back(0);
	d_gc.push_back(0);
	col_d.push_back(0);
	ext_v.push_back(0);
	ext_k.push_back(0);
	ext_21.push_back(0);
	rad_sb.push_back(0);
	rad_fl.push_back(0);
	r_ss.push_back(0);
	ang_ss.push_back(0);
	v_lsr.push_back(0);
	lbol.push_back(0);
	t_ms.push_back(0);
	ps_sw.push_back(0);
	logq.push_back(0);
	j21.push_back(0);
	r0_ss.push_back(0);
	ne.push_back(0);
	msx.push_back(0);
	kh.push_back(0);
	xyz[0].push_back(0);
	xyz[1].push_back(0);
	xyz[2].push_back(0);
	gcoord[0].push_back(0);
	gcoord[1].push_back(0);
	size++;
}

void StarPop::print(const std::string& filename) {
	// creating filename
	std::ofstream ofile(filename);
	if (!ofile)
		std::cerr << "ERROR: unable to open starpop.txt for writing." << std::endl;

	// writing data to file
	ofile << "M_fin[msun]\tage[kyr]\td_sun[kpc]\td_gc[kpc]\tcol_d[cm-3.kpc]\text_v\text_k\text_21\tx[kpc]\ty[kpc]\tz[kpc]\tg_long[deg]\tg_lat[deg]\n";

	ofile << std::setprecision(6) << std::fixed;
	for (int i = 0; i < size; i++) {
		ofile << mfin[i] << '\t';
		ofile << age[i] << '\t';
		ofile << dsun[i] << '\t';
		ofile << d_gc[i] << '\t';
		ofile << col_d[i] << '\t';
		ofile << ext_v[i] << '\t';
		ofile << ext_k[i] << '\t';
		ofile << ext_21[i] << '\t';
		ofile << xyz[0][i] << '\t';
		ofile << xyz[1][i] << '\t';
		ofile << xyz[2][i] << '\t';
		ofile << gcoord[0][i] << '\t';
		ofile << gcoord[1][i] << '\n';
	}
	ofile.close();
}

void StarPop::printFinal(const std::string& filename) {
	// Creating filename.
	std::ofstream ofile(filename);
	if (!ofile)
		throw std::runtime_error("ERROR: unable to open " + filename + " for writing.");

	// Writing data to file.
	ofile << "M_cur[msun]\tM_fin[msun]\tt_ms[kyr]\tage[kyr]\td_sun[kpc]\td_gc[kpc]\tcol_d[cm-3.kpc]\text_v\text_k\text_21\tlbol[lsun]\tj21[Jy]\tx[kpc]\ty[kpc]\tz[kpc]\tg_long[deg]\tg_lat[deg]\n";

	ofile << std::setprecision(9) << std::fixed;
	for (int i = 0; i < size; i++) {
		ofile << mcur[i] << '\t';
		ofile << mfin[i] << '\t';
		ofile << t_ms[i] << '\t';
		ofile << age[i] << '\t';
		ofile << dsun[i] << '\t';
		ofile << d_gc[i] << '\t';
		ofile << col_d[i] << '\t';
		ofile << ext_v[i] << '\t';
		ofile << ext_k[i] << '\t';
		ofile << ext_21[i] << '\t';
		ofile << lbol[i] << '\t';
		ofile << j21[i] << '\t';
		ofile << xyz[0][i] << '\t';
		ofile << xyz[1][i] << '\t';
		ofile << xyz[2][i] << '\t';
		ofile << gcoord[0][i] << '\t';
		ofile << gcoord[1][i] << std::endl;
	}
	ofile.close();
}

Star::Star() : id(0), mcur(0), mfin(0), age(0), dsun(0), d_gc(0), col_d(0), ext_v(0), ext_k(0), ext_21(0), sb_rad(0), fl_rad(0), r_ss(0), ang_ss(0), v_lsr(0), lbol(0), t_ms(0), ps_sw(0), logq(0), j21(0) {
	xyz[0] = 0;
	xyz[1] = 0;
	xyz[2] = 0;
	gcoord[0] = 0;
	gcoord[1] = 0;
}

Star::Star(int i, double m, double a) : id(i), mcur(0), mfin(m), age(a), dsun(0), d_gc(0), col_d(0), ext_v(0), ext_k(0), ext_21(0), sb_rad(0), fl_rad(0), r_ss(0), ang_ss(0), v_lsr(0), lbol(0), t_ms(0), ps_sw(0), logq(0), j21(0) {
	xyz[0] = 0;
	xyz[1] = 0;
	xyz[2] = 0;
	gcoord[0] = 0;
	gcoord[0] = 0;
}
