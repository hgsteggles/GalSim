#include "Galaxy.hpp"
#include "Random.hpp"
#include "MonteCarloCumulative.hpp"
#include "KroupaIMF.hpp"
#include "Timer.hpp"
#include "ProgressBar.hpp"
#include "HII.hpp"
#include "Logger.hpp"

#include <algorithm>

Galaxy::Galaxy(const GalaxyParameters& gal_params) : params(gal_params), thinDisk(nullptr), thickDisk(nullptr), spiralArms(nullptr) {
	rng.seed(time(0));
	nocells[0] = 2.0 * params.radius / params.resolution;
	nocells[1] = 2.0 * params.radius / params.resolution;
	nocells[2] = 2.0 * params.height / params.resolution;

	for (int i = 0; i < nocells[0]; i++)
		gal_x.push_back((i * params.resolution) - params.radius);
	for (int i = 0; i < nocells[1]; i++)
		gal_y.push_back((i * params.resolution) - params.radius);
	for (int i = 0; i < nocells[2]; i++)
		gal_z.push_back((i * params.resolution) - params.height);

	for (int i = 0; i < nocells[0]; i++) {
		std::vector<std::vector<double>> w;
		density.push_back(w);
		normdensity.push_back(w);
		for (int j = 0; j < nocells[1]; j++) {
			std::vector<double> v;
			density[i].push_back(v);
			normdensity[i].push_back(v);
			for (int k = 0; k < nocells[2]; k++) {
				density[i][j].push_back(0);
				normdensity[i][j].push_back(0);
			}
		}
	}

	Random::randomise();

	if (params.densityfile.compare("") != 0)
		readDensity(params.densityfile);
	else
		calcDensity(params.outputDirectory + "density3D.txt");
	print_xy_slice(params.outputDirectory + "density-xy-slice.txt");
	print_xy(params.outputDirectory + "density-xy.txt");
	print_xz(params.outputDirectory + "density-xz.txt");

	calcExtinction();

	if (params.starpopfile.compare("") != 0) {
		stars.read(params.starpopfile);
	}
	else {
		populate();
		putStars();
		stars.print(params.outputDirectory + "starpop.txt");
	}

	currStarMass();
	calcFluxes();
	stars.printFinal(params.outputDirectory + "starpopfinal.txt");
}

Galaxy::~Galaxy() {
}

double Galaxy::pix2kpc(double pixel, int dim) const {
	if (dim == 2)
		return (pixel * params.resolution - params.height);
	else
		return (pixel * params.resolution - params.radius);
}

int Galaxy::kpc2pix(double x, int dim) const {
	return (int)((x/params.resolution) + (nocells[dim]/2.0) + 0.5);
}

void Galaxy::readDensity(const std::string& filename) {
	Logger::Instance().print<SeverityType::NOTICE>("Setting up galactic density distribution.\n");
	Timer timer;
	timer.start();

	std::ifstream file(filename);
	if (!file)
		throw std::runtime_error("Galaxy::setupDensity: unable to open \'" + filename + "\' for reading.\n");

	thinDisk = std::unique_ptr<ThinDisk>(new ThinDisk(params));
	thickDisk = std::unique_ptr<ThickDisk>(new ThickDisk(params));
	spiralArms = std::unique_ptr<SpiralArms>(new SpiralArms(params));

	Logger::Instance().print<SeverityType::NOTICE>("Reading galaxyfile...\n");
	ProgressBar progBar(nocells[0], 1000);

	int n0, n1, n2;
	file >> n0 >> n1 >> n2;
	if (n0 != nocells[0] || n1 != nocells[1] || n2 != nocells[2])
		throw std::runtime_error("Galaxy::readDensity: mismatch between paramfile and densityfile dims.\n");
	std::string dummy;
	file >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

	for (int i = 0; i < nocells[0]; i++) {
		for (int j = 0; j < nocells[1]; j++) {
			for (int k = 0; k < nocells[2]; k++) {
				double gx, gy, gz, darms, dthick, dthin, dtot;
				file >> gx >> gy >> gz >> dthin >> dthick >> darms >> dtot;
				thinDisk->density[i][j][k] = dthin;
				thickDisk->density[i][j][k] = dthick;
				spiralArms->density[i][j][k] = darms;
				density[i][j][k] = dtot;
				if (dtot < 0)
					throw std::runtime_error("Galaxy::setupDensity: reading negative density.\n");

			}
		}
		if (progBar.timeToUpdate()) {
			progBar.update(i+1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	file.close();

	Logger::Instance().print<SeverityType::NOTICE>("Finished reading density.\n");
}

void Galaxy::calcDensity(const std::string& filename) {
	Logger::Instance().print<SeverityType::NOTICE>("Setting up galactic density distribution.\n");
	Timer timer;
	timer.start();

	spiralArms = std::unique_ptr<SpiralArms>(new SpiralArms(params));
	spiralArms->calcDensity();
	thickDisk = std::unique_ptr<ThickDisk>(new ThickDisk(params));
	thinDisk = std::unique_ptr<ThinDisk>(new ThinDisk(params));

	Logger::Instance().print<SeverityType::NOTICE>("Calculating density...\n");
	ProgressBar progBar(nocells[0], 1000);

	for (int i = 0; i < nocells[0]; i++) {
		for (int j = 0; j < nocells[1]; j++) {
			for (int k = 0; k < nocells[2]; k++) {
				if (spiralArms != nullptr)
					density[i][j][k] = spiralArms->density[i][j][k];
				if (thickDisk != nullptr)
					density[i][j][k] += thickDisk->density[i][j][k];
				if (thinDisk != nullptr)
					density[i][j][k] += thinDisk->density[i][j][k];
			}
		}

		if (progBar.timeToUpdate()) {
			progBar.update(i+1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}
	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	print3D(filename);

	Logger::Instance().print<SeverityType::NOTICE>("Finished calculating density.\n");
}

void Galaxy::print_xy_slice(const std::string& filename) {
	// creating filename
	std::ofstream ofile(filename);
	if (!ofile)
		throw std::runtime_error("Galaxy::print_xy_slice: unable to open \'" + filename + "\' for outputting\n.");

	// writing data to file
	Logger::Instance().print<SeverityType::NOTICE>("Printing 2D (xy) density slice...\n");
	ProgressBar progBar(density.size() * density[0].size(), 1000);

	ofile << gal_x.size() << '\t' << gal_y.size() << '\n';
	ofile << "gal_x[kpc]\tgal_y[kpc]\tden_thin_disk[cm-3]\tden_thick_disk[cm-3]\tden_spiral_arms[cm-3]\tden_total[cm-3]\n";

	ofile << std::setprecision(6) << std::fixed;
	for (unsigned int i = 0; i < density.size(); i++) {
		for (unsigned int j = 0; j < density[i].size(); j++) {
			ofile << gal_x[i] << '\t';
			ofile << gal_y[j] << '\t';

			double thinD = 0;
			double thickD = 0;
			double spiralA = 0;

			int k = density[i][j].size() / 2;

			if (thinDisk != nullptr)
				thinD += thinDisk->density[i][j][k];
			if (thickDisk != nullptr)
				thickD += thickDisk->density[i][j][k];
			if (spiralArms != nullptr)
				spiralA += spiralArms->density[i][j][k];

			ofile << thinD << '\t';
			ofile << thickD << '\t';
			ofile << spiralA << '\t';
			ofile << density[i][j][density[i][j].size()/2] << std::endl;

			if (progBar.timeToUpdate()) {
				progBar.update((i+1)*(j+1));
				Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
			}
		}
	}
	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	ofile.close();
}

void Galaxy::print_xy(const std::string& filename) {
	// creating filename
	std::ofstream ofile(filename);
	if (!ofile)
		throw std::runtime_error("Galaxy::print: unable to open \'" + filename + "\' for outputting\n.");

	// writing data to file
	Logger::Instance().print<SeverityType::NOTICE>("Printing 2D (xy) integrated density slice...\n");
	ProgressBar progBar(density.size() * density[0].size(), 1000);

	ofile << gal_x.size() << '\t' << gal_y.size() << '\n';
	ofile << "gal_x[kpc]\tgal_y[kpc]\tden_thin_disk[cm-3]\tden_thick_disk[cm-3]\tden_spiral_arms[cm-3]\tden_total[cm-3]\n";

	ofile << std::setprecision(6) << std::fixed;
	for (unsigned int i = 0; i < density.size(); i++) {
		for (unsigned int j = 0; j < density[i].size(); j++) {
			ofile << gal_x[i] << '\t';
			ofile << gal_y[j] << '\t';

			double thinD = 0;
			double thickD = 0;
			double spiralA = 0;

			for (unsigned int k = 0; k < density[i][j].size(); k++) {
				if (thinDisk != nullptr)
					thinD += thinDisk->density[i][j][k];
				if (thickDisk != nullptr)
					thickD += thickDisk->density[i][j][k];
				if (spiralArms != nullptr)
					spiralA += spiralArms->density[i][j][k];
			}

			ofile << thinD << '\t';
			ofile << thickD << '\t';
			ofile << spiralA << '\t';
			ofile << density[i][j][density[i][j].size()/2] << std::endl;

			if (progBar.timeToUpdate()) {
				progBar.update((i+1)*(j+1));
				Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
			}
		}
	}
	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	ofile.close();
}

void Galaxy::print_xz(const std::string& filename) {
	// creating filename
	std::ofstream ofile(filename);
	if (!ofile)
		throw std::runtime_error("Galaxy::print: unable to open \'" + filename + "\' for outputting.\n");

	// writing data to file
	Logger::Instance().print<SeverityType::NOTICE>("Printing 2D (xz) integrated density slice...\n");
	ProgressBar progBar(density.size() * density[0].size(), 1000);

	ofile << gal_x.size() << '\t' << gal_z.size() << '\n';
	ofile << "gal_x[kpc]\tgal_z[kpc]\tden_thin_disk[cm-3]\tden_thick_disk[cm-3]\tden_spiral_arms[cm-3]\tden_total[cm-3]\n";

	ofile << std::setprecision(6) << std::fixed;
	for (unsigned int i = 0; i < density.size(); i++) {
		for (unsigned int k = 0; k < density[i][0].size(); k++) {
			ofile << gal_x[i] << '\t';
			ofile << gal_z[k] << '\t';

			double thinD = 0;
			double thickD = 0;
			double spiralA = 0;
			double total = 0;

			for (unsigned int j = 0; j < density[i].size(); j++) {
				if (thinDisk != nullptr)
					thinD += thinDisk->density[i][j][k];
				if (thickDisk != nullptr)
					thickD += thickDisk->density[i][j][k];
				if (spiralArms != nullptr)
					spiralA += spiralArms->density[i][j][k];

				total += density[i][j][k];
			}

			ofile << thinD << '\t';
			ofile << thickD << '\t';
			ofile << spiralA << '\t';
			ofile << total << std::endl;

			if (progBar.timeToUpdate()) {
				progBar.update((i+1)*(k+1));
				Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
			}
		}
	}
	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	ofile.close();
}

void Galaxy::print3D(const std::string& filename) {
	std::ofstream ofile(filename);
	if (!ofile)
		throw std::runtime_error("Galaxy::print3D: unable to open \'" + filename + "\' for outputting.\n");

	// writing data to file
	Logger::Instance().print<SeverityType::NOTICE>("Printing 3D density data...\n");
	ProgressBar progBar(density.size(), 1000);

	ofile << gal_x.size() << '\t' << gal_y.size() << '\t' << gal_z.size() << '\n';
	ofile << "gal_x[kpc]\tgal_y[kpc]\tgal_z[kpc]\tden_thin_disk[cm-3]\tden_thick_disk[cm-3]\tden_spiral_arms[cm-3]\tden_total[cm-3]\n";
	ofile << std::setprecision(9) << std::fixed;

	for (unsigned int i = 0; i < density.size(); i++) {
		for (unsigned int j = 0; j < density[i].size(); j++) {
			for (unsigned int k = 0; k < density[i][j].size(); k++) {
				ofile << gal_x[i] << '\t';
				ofile << gal_y[j] << '\t';
				ofile << gal_z[k] << '\t';
				if (thinDisk != nullptr)
					ofile << thinDisk->density[i][j][k] << '\t';
				else
					ofile << 0 << '\t';
				if (thickDisk != nullptr)
					ofile << thickDisk->density[i][j][k] << '\t';
				else
					ofile << 0 << '\t';
				if (spiralArms != nullptr)
					ofile << spiralArms->density[i][j][k] << '\t';
				else
					ofile << 0 << '\t';
				ofile << density[i][j][k] << std::endl;

				if (progBar.timeToUpdate()) {
					progBar.update(i + 1);
					Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
				}
			}
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	ofile.close();
}

void Galaxy::calcExtinction() {
	Logger::Instance().print<SeverityType::NOTICE>("Calculating extinction parameters from cepheids...\n");

	Timer timer;
	timer.start();

	int nsamp = 20;
	std::vector<std::vector<double> > clust_r(3, std::vector<double>()), clustAk(2, std::vector<double>());
	std::vector<std::string> id1, id2;
	std::vector<double> clust_L, clust_B, per, vint, bvint, ebv, vamp, mv, dist, ra3, de2, cl, cb;
	std::vector<int> ra1, ra2, de1;
	std::vector<std::string> starn, starn2, clust_id, ref;

	//** CLUSTERS **//
	readcols("config/required/clust_ak.txt", ",", clust_id, clust_L, clust_B, clust_r[0], clust_r[1], clust_r[2], clustAk[0], clustAk[1], ref);

	int nclust = clust_r[0].size();
	for (int iclust = 0; iclust < nclust; iclust++) {
		clust_r[1][iclust] = clust_r[0][iclust]+clust_r[1][iclust];
		clust_r[2][iclust] = clust_r[0][iclust]-clust_r[2][iclust];
	}

	// Calculate x,y,z coords & errors of clusters from r, r+dr, r-dr, l, l+dl, l-dl, b, b+db, b-db
	std::vector<std::vector<double> > cl_x(3, std::vector<double>(nclust, 0));
	std::vector<std::vector<double> > cl_y(3, std::vector<double>(nclust, 0));
	std::vector<std::vector<double> > cl_z(3, std::vector<double>(nclust, 0));
	for (int iclust = 0; iclust < nclust; iclust++) {
		for (int j = 0; j < 3; j++) {
			cl_x[j][iclust] = clust_r[j][iclust] * std::cos(clust_B[iclust] * CST::DEG2RAD) * std::sin(clust_L[iclust]*CST::DEG2RAD);
			cl_y[j][iclust] = -1.0 * clust_r[j][iclust] * std::cos(clust_B[iclust] * CST::DEG2RAD) * std::cos(clust_L[iclust] * CST::DEG2RAD) + params.solar_position[1];
			cl_z[j][iclust] = clust_r[j][iclust] * std::sin(clust_B[iclust] * CST::DEG2RAD);
		}
	}

	std::vector<std::vector<double> > clustColDens(3, std::vector<double>(nclust, 0.0));

	Logger::Instance().print<SeverityType::NOTICE>("Calculating cluster column densities...\n");
	ProgressBar progBar(nclust, 1000);

	// Calculate cluster column densities.
	for (int iclust = 0; iclust < nclust; iclust++) {
		for (int j = 0; j < 3; j++) {
			std::vector<double> xsamp_px(nsamp, 0.0), ysamp_px(nsamp, 0.0), zsamp_px(nsamp, 0.0);
			for (int isamp = 0; isamp < nsamp; isamp++) {
				xsamp_px[isamp] = kpc2pix(isamp * (cl_x[j][iclust] - params.solar_position[0]) / nsamp + params.solar_position[0], 0);
				ysamp_px[isamp] = kpc2pix(isamp * (cl_y[j][iclust] - params.solar_position[1]) / nsamp + params.solar_position[1], 1);
				zsamp_px[isamp] = kpc2pix(isamp * (cl_z[j][iclust] - params.solar_position[2]) / nsamp + params.solar_position[2], 2);
			}
			std::vector<double> interpolatedDensities = cubeInterpol(density, xsamp_px, ysamp_px, zsamp_px);
			std::vector<double> distances(nsamp, 0.0);
			for (int isamp = 0; isamp < nsamp; isamp++) {
				distances[isamp] = isamp * clust_r[j][iclust] / (double)(nsamp);
			}
			clustColDens[j][iclust] = int_tabulated(distances, interpolatedDensities);
		}
		clustColDens[1][iclust] = std::max(clustColDens[1][iclust], clustColDens[2][iclust]) - clustColDens[0][iclust];
		clustColDens[2][iclust] = clustColDens[0][iclust] - std::min(clustColDens[1][iclust], clustColDens[2][iclust]);

		if (progBar.timeToUpdate()) {
			progBar.update(iclust + 1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	//** CEPHEIDS **//
	readcols("config/required/cepheids.csv", ",", id1, starn, per, vint, bvint, ebv, vamp, mv, dist);
	readcols("config/required/cepheids_pos.csv", ",", id2, starn2, ra1, ra2, ra3, de1, de2, cl, cb);

	std::vector<std::string> ceph_id;
	std::vector<double> ceph_l, ceph_b, ceph_r, ceph_bv;
	int nceph = id2.size();

	for (int iceph = 0; iceph < nceph; iceph++) {
		std::string curr_id2 = id2[iceph];
		std::vector<int> coord = indicesWhereTrue(id1, [curr_id2](const std::string& str){return str == curr_id2;});
		//std::vector<int> coord = where_eq(id1, id2[iceph], nmatch);
		if (coord.size() > 0) {
			ceph_l.push_back(cl[iceph]);
			ceph_b.push_back(cb[iceph]);
			ceph_r.push_back(dist[coord[0]] / 1000.0);
			ceph_bv.push_back(ebv[coord[0]]);
			ceph_id.push_back(id2[iceph]);
		}
	}

	nceph = ceph_r.size();

	// Calculate x, y, z coords of cepheids from r, l, b
	std::vector<double> ceph_x(nceph, 0.0), ceph_y(nceph, 0.0), ceph_z(nceph, 0.0);

	for (int i = 0; i < nceph; i++) {
		ceph_x[i] = ceph_r[i] * std::cos(ceph_b[i] * CST::DEG2RAD) * std::sin(ceph_l[i] * CST::DEG2RAD);
		ceph_y[i] = -1.0 * (ceph_r[i] * std::cos(ceph_b[i] * CST::DEG2RAD) * std::cos(ceph_l[i] * CST::DEG2RAD)) + params.solar_position[1];
		ceph_z[i] = ceph_r[i] * std::sin(ceph_b[i] * CST::DEG2RAD);
	}

	Logger::Instance().print<SeverityType::NOTICE>("Calculating cepheid column densities...\n");
	progBar = ProgressBar(100, 1000);

	// Calculate cepheid column densities.
	std::vector<double> cephColDens(nceph, 0.0);

	for (int iceph = 0; iceph < nceph; iceph++) {
		std::vector<double> xsamp_px(nsamp, 0.0), ysamp_px(nsamp, 0.0), zsamp_px(nsamp, 0.0);
		for (int isamp = 0; isamp < nsamp; isamp++) {
			//Sampling indices in pixels.
			xsamp_px[isamp] = kpc2pix(isamp * (ceph_x[iceph] - params.solar_position[0]) / nsamp + params.solar_position[0], 0);
			ysamp_px[isamp] = kpc2pix(isamp * (ceph_y[iceph] - params.solar_position[1]) / nsamp + params.solar_position[1], 1);
			zsamp_px[isamp] = kpc2pix(isamp * (ceph_z[iceph] - params.solar_position[2]) / nsamp + params.solar_position[2], 2);
		}
		// Calculate density profile along LOS.
		std::vector<double> interpolatedDensities = cubeInterpol(density, xsamp_px, ysamp_px, zsamp_px);
		std::vector<double> distances(nsamp, 0.0);
		for (int isamp = 0; isamp < nsamp; isamp++)
			distances[isamp] = isamp * ceph_r[iceph] / (double)(nsamp);
		cephColDens[iceph] = int_tabulated(distances, interpolatedDensities);

		if (progBar.timeToUpdate()) {
			progBar.update(100.0 * iceph / (double)(nceph));
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	std::vector<double> fitx = clustColDens[0];
	std::vector<double> fity = clustAk[0];
	std::vector<double> xsig = clustColDens[2];
	std::vector<double> ysig = clustAk[1];
	fitx.push_back(vmean(cephColDens));
	fity.push_back(vmean(ceph_bv) * 3.1 * 0.112);
	xsig.push_back(vstddev(cephColDens));
	ysig.push_back(vstddev(ceph_bv) * 3.1 * 0.112);
	double asig, bsig;

	Recipes::fitexy(fitx, fity, xsig, ysig, extinction_a, extinction_b, asig, bsig);
	Logger::Instance().print<SeverityType::NOTICE>("extpar_a = ", extinction_a, '\t', "extpar_b = ", extinction_b, '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Finished calculating extinction parameters. ", timer.formatTime(timer.getTicks()), '\n');
}

void Galaxy::populate() {
	Logger::Instance().print<SeverityType::NOTICE>("Generating stars.\n");

	Timer timer;
	timer.start();

	double sfr_gal = params.sfr; //SFR of Galaxy in Msun/yr.
	double total_time = params.total_time; //total time of simulation in yrs.
	double min_outmass = 6; //minimum mass to catalogue (Msun).

	double total_stellar_mass = sfr_gal*total_time;
	/* Generating Stellar Populations */
	/* Single Star Mode */
	Logger::Instance().print<SeverityType::NOTICE>(
		"Generating stellar pops in single star mode...");
	Logger::Instance().print<SeverityType::NOTICE>("Total stellar mass in = ",
		total_stellar_mass, '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Calc. star masses & ages...\n");
	ProgressBar progBar(total_stellar_mass, 1000);

	KroupaIMF kroupa;

	double sum = 0;
	for (int istar = 0, i = 0; i < (int)(10.0 * total_stellar_mass) && sum < total_stellar_mass; i++) {
		double st_mfin = kroupa.randomMass();
		sum += st_mfin;
		if (st_mfin >= min_outmass) {
			istar++;
			double st_age = Random::randomDouble(0.0, 1.0) * total_time / 1.0e3;
			stars.addStar(istar, st_mfin, st_age);
		}
		if (sum > total_stellar_mass)
			break;

		if (progBar.timeToUpdate()) {
			progBar.update(sum);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Total stellar mass out = ", sum, '\n');
	Logger::Instance().print<SeverityType::NOTICE>("No. of stars = ", stars.id.size(), '\n');
	Logger::Instance().print<SeverityType::NOTICE>("Finished generating stars.\n");
}

static std::array<int, 3> coordsFromIndex(int index, int nx, int ny) {
	std::array<int, 3> ret;
	ret[2] = index / (nx * ny);
	index -= ret[2] * nx * ny;

	ret[1] = index / nx;
	index -= ret[1] * nx;

	ret[0] = index / 1;

	return ret;
}

void Galaxy::putStars() {
	Logger::Instance().print<SeverityType::NOTICE>("Putting stars into galaxy.\n");

	Timer timer;
	timer.start();

	int N = nocells[0]*nocells[1]*nocells[2];

	std::vector<std::vector<std::vector<double>>>& dencube = spiralArms->density;

	// Find max density for normalisation.
	double maxDensity = 0;
	for (int i = 0; i < N; i++) {
		std::array<int, 3> pos = coordsFromIndex(i, nocells[0], nocells[1]);
		maxDensity = std::max(maxDensity, dencube[pos[0]][pos[1]][pos[2]]);
	}

	double SF_POWER = 1.4;
	std::vector<double> probs = std::vector<double>(N, 0.0);

	Logger::Instance().print<SeverityType::NOTICE>("Calculating P(formation|position)...\n");
	ProgressBar progBar(N, 1000);

	for (int i = 0; i < N; i++) {
		std::array<int, 3> pos = coordsFromIndex(i, nocells[0], nocells[1]);
		probs[i] = std::pow(dencube[pos[0]][pos[1]][pos[2]] / maxDensity, SF_POWER);

		if (progBar.timeToUpdate()) {
			progBar.update(i + 1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	std::vector<std::size_t> prob_indexes(N);
	std::iota(prob_indexes.begin(), prob_indexes.end(), 0);

	// Remove zero prob indexes.
	prob_indexes.erase(std::remove_if(prob_indexes.begin(),
									  prob_indexes.end(),
									  [&probs](std::size_t i){return probs[i] < 1.0e-20;}),
					   prob_indexes.end());

	// Sort indexes by indexed density.
	std::sort(prob_indexes.begin(), prob_indexes.end(),
			  [&probs](std::size_t i1, std::size_t i2) {return probs[i1] < probs[i2];});

	// Compute Extinction Array
	Logger::Instance().print<SeverityType::NOTICE>("Calculating star positions & extinctions...\n");
	progBar = ProgressBar(stars.size, 1000);

	for (int istar = 0; istar < stars.size; istar++) {
		// Calculate Star Position
		double rnum = Random::randomDouble(0, 1);
		// Find lowest index of probs_indexes that has high enough prob.
		int ilow = 0;
		for (int i = 0; prob_indexes.size(); ++i) {
			if (probs[prob_indexes[i]] > rnum) {
				ilow = i;
				break;
			}
		}
		// Randomly select from indexes with high enough indexed prob.
		int index = prob_indexes[Random::randomInteger(ilow, prob_indexes.size()-1)];
		std::array<int, 3> intCoords = coordsFromIndex(index, nocells[0], nocells[1]);

		std::array<double, 3> pixCoords;
		pixCoords[0] = intCoords[0];
		pixCoords[1] = intCoords[1];
		pixCoords[2] = intCoords[2];

		for (int i = 0; i < 3; i++) {
			pixCoords[i] += 0.5 * (2.0 * Random::randomDouble(0.0, 1.0) - 1.0);
			stars.xyz[i][istar] = pix2kpc(pixCoords[i], i);
		}

		// Calculate Distance to Sun
		double dist2sun2 = 0;
		for (int i = 0; i < 3; i++)
			dist2sun2 += std::pow(stars.xyz[i][istar] - params.solar_position[i], 2.0);
		stars.dsun[istar] = std::sqrt(dist2sun2);

		// Integrate Gas Column Density Along LOS
		int nsamp = 20;

		std::vector<double> xsamp_px(nsamp, 0.0), ysamp_px(nsamp, 0.0), zsamp_px(nsamp, 0.0);
		for (int isamp = 0; isamp < nsamp; isamp++) {
			xsamp_px[isamp] = kpc2pix(isamp * (stars.xyz[0][istar] - params.solar_position[0]) / nsamp + params.solar_position[0], 0);
			ysamp_px[isamp] = kpc2pix(isamp * (stars.xyz[1][istar] - params.solar_position[1]) / nsamp + params.solar_position[1], 1);
			zsamp_px[isamp] = kpc2pix(isamp * (stars.xyz[2][istar] - params.solar_position[2]) / nsamp + params.solar_position[2], 2);
		}

		std::vector<double> interpolatedDensities = cubeInterpol(dencube, xsamp_px, ysamp_px, zsamp_px);
		std::vector<double> distances(nsamp, 0.0);
		for (int isamp = 0; isamp < nsamp; isamp++)
			distances[isamp] = isamp * stars.dsun[istar] / nsamp;

		stars.col_d[istar] = int_tabulated(distances, interpolatedDensities); //cm-3 kpc
		if (stars.col_d[istar] < 0) {
			throw std::runtime_error("Galaxy::putStars: negative column density encountered.\n"); 
		}

		// V-band Extinction
		stars.ext_v[istar] = (stars.col_d[istar]) * 100.0;

		// K-band Extinction
		double buff = pow(10.0, std::log10(stars.col_d[istar] * CST::PC2CM * 1.0e3)-21); //10^21 cm^-2
		double kext = extinction_a + buff*extinction_b;
		if (kext < 0)
			kext = 0;
		stars.ext_k[istar] = kext;

		// 21um Extinction
		stars.ext_21[istar] = kext * 0.5;

		// Galactic Coords
		double sunstar_x = stars.xyz[0][istar] - params.solar_position[0];
		double sunstar_y = stars.xyz[1][istar] - params.solar_position[1];
		double sunstar_z = stars.xyz[2][istar] - params.solar_position[2];
		double lgal = std::atan(sunstar_x / (-1.0 * sunstar_y)) / CST::DEG2RAD;
		double bgal = std::atan(sunstar_z / std::sqrt(std::pow(sunstar_x, 2.0) + std::pow(sunstar_y, 2.0))) / CST::DEG2RAD;

		if (sunstar_y > 0)
			lgal += 180;
		else if (sunstar_y < 0 && sunstar_x < 0)
			lgal += 360;

		stars.gcoord[0][istar] = lgal;
		stars.gcoord[1][istar] = bgal;
		stars.d_gc[istar] = std::sqrt(pow(stars.xyz[0][istar], 2.0) + pow(stars.xyz[1][istar], 2.0) + pow(stars.xyz[2][istar], 2.0));

		// Progress
		if (progBar.timeToUpdate()) {
			progBar.update(istar);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Finished placing stars.\n");
}

void Galaxy::currStarMass() {
	Logger::Instance().print<SeverityType::NOTICE>("Calculating current star masses.\n");

	std::vector<double> isomass, isotemp, isologl, isorad, isobc, isomv, isologq, isotkh;
	std::vector<double> mdot_ms_arr, mass_ms_arr;
	std::vector<std::string> dum;
	std::string filename = "config/required/lyman_101108_corr-highmass.dat";
	std::string filename2 = "config/required/hosokawa_tracks_interp/mass_ms_fine.dat";
	std::string delimiter = "\t";
	readcols(filename, delimiter, isomass, isotemp, isologl, isorad, isobc, isomv, isologq, isotkh, dum);
	readcols(filename2, delimiter, mdot_ms_arr, mass_ms_arr, dum, dum, dum, dum, dum, dum, dum);

	int nstars = stars.id.size();
	std::vector<double> mt_mcur(nstars, 0), mt_acc(nstars, 0), mt_tform(nstars, 0), t20(nstars, 0);
	std::vector<double> mt_lacc(nstars, 0);
	double MASS_MS = 20; //Msun
	double SIG_cl = 1.0; //g cm^-3

	Logger::Instance().print<SeverityType::NOTICE>("Computing MS masses...\n");
	ProgressBar progBar(nstars, 1000);

	for (int i = 0; i < nstars; i++) {
		//mt_mcur[i] = 0.18*pow(stars.mfin[i]/30.0, 0.5)*pow(SIG_cl, 1.5)*pow(stars.age[i]*1.0e3/1.0e4, 2.0);
		mt_mcur[i] = 5.4e-8 * std::pow(stars.mfin[i] / 30.0, 1.5) * std::pow(SIG_cl, 1.5) * pow(stars.age[i] * 1.0e3, 2.0) / stars.mfin[i];
		if (mt_mcur[i] > stars.mfin[i]) {
			mt_mcur[i] = stars.mfin[i];
			mt_acc[i] = 0.0;
		}
		else
			mt_acc[i] = 4.6e-4 * std::pow(stars.mfin[i] / 30.0, 0.75) * std::pow(SIG_cl, 0.75) * std::pow(mt_mcur[i] / stars.mfin[i], 0.5);

		mt_acc[i] = round(std::log10(mt_acc[i]) * 50.0) / -50.0;

		// Find MASS_MS of accreting objects.
		mt_tform[i] = 1.29e2 * std::pow(stars.mfin[i] / 30.0, 0.25) * std::pow(SIG_cl, -0.75); // kyrs since finished accreting.
		t20[i] = 1.1e2 * std::pow(MASS_MS / 20, 0.5) * std::pow(stars.mfin[i] / 30.0, -0.25) * std::pow(SIG_cl, -0.75); //kyrs since MS.

		if (mt_mcur[i] < stars.mfin[i]) {
			stars.mcur[i] = mt_mcur[i];
		}
		else {
			stars.mcur[i] = stars.mfin[i];
			stars.lbol[i] = std::pow(10.0, interpol(isologl, isomass, stars.mfin[i]));
		}

		stars.t_ms[i] = std::max(0.0, stars.age[i] - std::min(mt_tform[i], t20[i]));

		if (progBar.timeToUpdate()) {
			progBar.update(i + 1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	HosTracks htracks("config/required/hosokawa_tracks_interp/int_md.dat"); //hostracks.cpp
	int ntracks = htracks.tracks.size();
	std::vector<double> hos_mdot_arr(ntracks, 0.0);
	for (int i = 0; i < ntracks; i++)
		hos_mdot_arr[i] = htracks.tracks[i].mdot;
	double hos_mdot_res = hos_mdot_arr[1] - hos_mdot_arr[0];

	// Luminosities of accreting objects.
	Logger::Instance().print<SeverityType::NOTICE>("Computing luminosities...\n");
	progBar = ProgressBar(nstars, 1000);

	for (int i = 0; i < nstars; i++) {
		double acc_use;
		if (mt_mcur[i] < stars.mfin[i]) {
			if (mt_acc[i] < 3.0)
				acc_use = 3.0;
			else if (mt_acc[i] > 5.0)
				acc_use = 5.0;
			else
				acc_use = mt_acc[i];

			std::vector<double> mdot_acc = hos_mdot_arr;
			for (int j = 0; j < ntracks; j++)
				mdot_acc[j] = fabs(mdot_acc[j] - acc_use);

			std::vector<int> tracki = indicesWhereTrue(mdot_acc, [](double val){return val <= 1.0e-4;});
			std::vector<int> trackj = indicesWhereTrue(mdot_acc, [hos_mdot_res](double val){return val <= hos_mdot_res;});
			//std::vector<int> tracki = where_lteq(mdot_acc, 1.0e-4, ntracki);
			//std::vector<int> trackj = where_lteq(mdot_acc, hos_mdot_res, ntracki);
			std::vector<double> track_masses = htracks.tracks[tracki[0]].data.mass;
			std::vector<double> track_ltots = htracks.tracks[tracki[0]].data.ltot;
			std::vector<double> track_lstar = htracks.tracks[tracki[0]].data.lstar;
			stars.lbol[i] = interpol(track_ltots, track_masses, stars.mcur[i]);
			mt_lacc[i] = stars.lbol[i] - interpol(track_lstar, track_masses, stars.mcur[i]);
		}

		// Luminosities of accretors above MASS_MS
		double mt_lacc_const = 0.5 * (6.4e22 * 6.67e-11) * (2.0e30 / 4.0e26) / 7.0e8;
		double mt_temp = std::pow(10.0, interpol(isotemp, isomass, stars.mcur[i]));
		double mt_rad = std::sqrt(4.0e26 / (4.0 * CST::PI * 5.67e-8)) * std::sqrt(stars.lbol[i]/std::pow(mt_temp, 4.0)) / 7.0e8;
		if (mt_mcur[i] < stars.mfin[i] && mt_mcur[i] > MASS_MS) {
			mt_rad = interpol(isorad, isomass, stars.mcur[i]);
			mt_lacc[i] = mt_lacc_const*std::pow(10.0, -1.0 * mt_acc[i]) * stars.mcur[i] / mt_rad;
			stars.lbol[i] += mt_lacc[i];
		}

		if (progBar.timeToUpdate()) {
			progBar.update(i + 1);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');

	Logger::Instance().print<SeverityType::NOTICE>("Finished calculating current masses.\n");
}

void Galaxy::calcFluxes() {
	Logger::Instance().print<SeverityType::NOTICE>("Calculating fluxes.\n");

	double LYMANSCALE = 1.0;
	double F_21_SIG = 0.26;
	double F_21_SC = 1.40;
	double WID_21 = 0.42e13;
	double SOL_LUM = 4.0e26;
	double MASS_MS = 20.0;
	double KPC2M = 3.08567758e19;
	std::vector<double> isomass, isotemp, isologl, isorad, isobc, isomv, isologq, isotkh;
	std::vector<std::string> dum;
	std::string filename = "config/required/lyman_101108_corr-highmass.dat";
	std::string delimiter = "\t";
	readcols(filename, delimiter, isomass, isotemp, isologl, isorad, isobc, isomv, isologq, isotkh, dum);
	int ndata = isomass.size();
	int nstars = stars.mcur.size();

	//Find 21um flux level.
	for (int i = 0; i < nstars; i++) {
		double j21_sc_exp = Random::randomNormal(0.0, 1.0) * F_21_SIG + F_21_SC;
		double j21_scale = std::pow(10.0, j21_sc_exp);
		stars.j21[i] = stars.lbol[i] * SOL_LUM * std::pow(10.0, stars.ext_21[i] / -2.5) / (4.0 * CST::PI * std::pow(stars.dsun[i] * KPC2M, 2.0) * WID_21 * 1e-26 * j21_scale);
	}

	//Recompute Lyman flux.
	std::vector<double> ref_mstar = isomass;
	std::vector<double> ref_Q0 = isologq;
	for (int i = 0; i < ndata; i++)
		ref_Q0[i] = ref_Q0[i] + std::log10(LYMANSCALE); //Scaling factor of Lyman flux.

	logq_mms = interpol(ref_Q0, ref_mstar, MASS_MS);
	for (int i = 0; i < nstars; i++)
		stars.logq[i] = interpol(ref_Q0, ref_mstar, stars.mcur[i]);

	Logger::Instance().print<SeverityType::NOTICE>("Finished calculating fluxes.\n");
}








