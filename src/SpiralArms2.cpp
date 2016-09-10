#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <iomanip> //std::setprecision
#include <math.h>

#include "SpiralArms2.hpp"
#include "Timer.hpp"
#include "ProgressBar.hpp"
#include "Logger.hpp"

static double h(double x) {
	return pow(1.0 / std::cosh(x), 2.0);
}

static double Ga(double x, double y, double z, double A_a, double A_inner, const ArmArray& na_fj, const ArmArray& wa_wj, const ArmArray& ha_hj, std::vector<double> sj) {
	double result = 0;
	for (int m = 0; m < 4; m++) {
		double r = std::sqrt(x * x + y * y + z * z);

		double ga = std::exp(-1.0 * std::pow(sj[m] / wa_wj[m], 2.0));

		if (r > A_a)
			ga *= h(0.5 * (r - A_a));

		result += na_fj[m] * ga * h(z / ha_hj[m]);
	}

	return result;
}

SpiralArms::SpiralArms(const GalaxyParameters& gparams) : params(gparams) {
	int nocells[3];
	nocells[0] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[1] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[2] = (int)(2.0*params.height/params.resolution + 0.5);
	for (int i = 0; i < nocells[0]; i++) {
		std::vector<std::vector<double> > w;
		density.push_back(w);
		for (int j = 0; j < nocells[1]; j++) {
			std::vector<double> v;
			density[i].push_back(v);
			for (int k = 0; k < nocells[2]; k++)
				density[i][j].push_back(0);
		}
	}

}

void SpiralArms::calcDensity() {
	Logger::Instance().print<SeverityType::NOTICE>("Setting up SpiralArms.\n");

	int nocells[3];
	nocells[0] = density.size();
	nocells[1] = density[0].size();
	nocells[2] = density[0][0].size();
	std::vector<double> xaxis(nocells[0]);
	std::vector<double> yaxis(nocells[1]);
	std::vector<double> zaxis(nocells[2]);
	for (unsigned int i = 0; i < xaxis.size(); i++)
		xaxis[i] = -params.radius + (i * params.resolution);
	for (unsigned int i = 0; i < yaxis.size(); i++)
		yaxis[i] = -params.radius + (i * params.resolution);
	for (unsigned int i = 0; i < zaxis.size(); i++)
		zaxis[i] = -params.height + (i * params.resolution);

	double rmax = std::sqrt(2.0) * params.radius;
	double sp_rmax[4] = {rmax, rmax, rmax, rmax};
	int nfit = 1000;
	
	Timer timer;
	timer.start();

	Logger::Instance().print<SeverityType::NOTICE>("Computing SJ array...\n");
	ProgressBar progBar(100, 1000);

	std::vector<std::vector<std::vector<double> > > sj(nocells[0], std::vector<std::vector<double> >(nocells[1], std::vector<double>(4, 0.0)));
	for (int m = 0; m < 4; ++m) {
		std::vector<double> sp_r(nfit);
		for (unsigned int i = 0; i < sp_r.size(); i++)
			sp_r[i] = std::pow(i / (nfit - 1.0), 3.0) * (sp_rmax[m] - params.rmin) + params.rmin;

		std::vector<double> sp_th(nfit);
		for (unsigned int i = 0; i < sp_th.size(); i++)
			sp_th[i] = params.sp_asp[m] * std::log10(sp_r[i] / params.sp_rmin[m] ) + params.sp_thmin[m];
		/* FIND JOIN POINTS AND START VECTORS FROM THERE */

		int imin = 0;
		double thmin = std::fabs(sp_th[0] - params.sp_joins[m]);
		for (unsigned int i = 0; i < sp_th.size(); i++) {
			if (std::fabs(sp_th[i] - params.sp_joins[m]) < thmin && sp_r[i] < 6) {
				thmin = std::fabs(sp_th[i] - params.sp_joins[m]);
				imin = i;
			}
		}
		if (imin != 0) {
			sp_r.erase(sp_r.begin(), sp_r.begin()+imin-1);
			sp_th.erase(sp_th.begin(), sp_th.begin()+imin-1);
		}

		/*
		if (m == 2) {
			for (unsigned int i = 0; i < sp_th.size(); i++) {
				if(sp_th[i] >= 30.0 && sp_th[i] < 110.0)
					sp_r[i] = sp_r[i] - 0.7*std::exp(-1.0*std::pow((70.0-sp_th[i])/30.0, 2.0));
			}
		}
		*/

		std::vector<double> sp_x(sp_r.size());
		std::vector<double> sp_y(sp_r.size());
		for (unsigned int i = 0; i < sp_x.size(); i++) {
			sp_x[i] = sp_r[i] * std::cos(sp_th[i] * CST::PI/180.0);
			sp_y[i] = sp_r[i] * std::sin(sp_th[i] * CST::PI/180.0);
		}

		for (int ix = 0; ix < nocells[0]; ix++) {
			for (int iy = 0; iy < nocells[1]; iy++) {
				std::vector<double> d2arm(sp_x.size());
				int mini = 0;
				for (unsigned int i = 0; i < d2arm.size(); i++) {
					d2arm[i] = std::sqrt(std::pow(sp_x[i]-xaxis[ix], 2.0) + std::pow(sp_y[i]-yaxis[iy], 2.0));
					if(d2arm[i] < d2arm[mini])
						mini = i;
				}

				int int_wid  = 3;
				int minp = std::max(mini - int_wid, 0);
				int maxp = std::min(mini + int_wid, (int)d2arm.size()-1);
				int nsub = maxp - minp + 1;
				
				std::vector<double> subd = slice(d2arm, minp, maxp);
				std::vector<double> grad = deriv(subd);
				std::vector<double> xin = findgen(nsub, 1.0);
				int int_res  = 5;
				std::vector<double> xout(int_res*nsub, 1.0/int_res);
				// FIND LOCAL MINIMUM
				if(vmin(grad) < 0 && vmax(grad) > 0) {
					double d0 = interpol(xin, grad, 0);
					double r0 = interpol(subd, xin, d0);
					sj[ix][iy][m] = r0;
				}
				else {
					sj[ix][iy][m] = vmin(subd);
					if (sj[ix][iy][m] != sj[ix][iy][m])
						throw std::runtime_error("SpiralArms::calcDensity: NaN Sj, check arm parameters (rmin > 0).");
				}
			}

			if (progBar.timeToUpdate()) {
				progBar.update(ix * 25.0 / nocells[0] + 25.0 * m);
				Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
			}
		}
	}
	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');
	
	Logger::Instance().print<SeverityType::NOTICE>("Computing SpiralArm densities...\n");
	progBar = ProgressBar(nocells[0], 1000);

	for (int i = 0; i < nocells[0]; i++) {
		for (int j = 0; j < nocells[1]; j++) {
			for (int k = 0; k < nocells[2]; k++) {
				double x = i*params.resolution - params.radius;
				double y = j*params.resolution - params.radius;
				double z = k*params.resolution - params.height;
				density[i][j][k] = Ga(x, y, z, params.Aa, params.A_inner, params.na_fj, params.wa_wj, params.ha_hj, sj[i][j]);

				if (density[i][j][k] != density[i][j][k])
					throw std::runtime_error("SpiralArms::calcDensity: NaN density.\n");
			}
		}

		if (progBar.timeToUpdate()) {
			progBar.update(i + 1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Finished setting up SpiralArms.\n");
}
