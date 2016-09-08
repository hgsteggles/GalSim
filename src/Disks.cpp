#include "Disks.hpp"

static double G1(double r, double z, double n1, double A1, double Aa, double H1){
	double g1 = pow(cosh(Aa/A1)/cosh(r/A1), 2.0);
	double h = pow(1.0/cosh(z/H1), 2.0);
	return n1*g1*h;
}

static double G2(double r, double z, double n2, double A2, double H2){
	double g2 = exp(-1.0* pow((r-A2)/1.8, 2.0));
	double h = pow(1.0/cosh(z/H2), 2.0);
	return n2 * g2 * h;
}

ThinDisk::ThinDisk(const GalaxyParameters& gparams) : params(gparams) {
	std::cout << "|_ Setting up ThinDisk..." << std::endl;
	double nocells[3];
	nocells[0] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[1] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[2] = (int)(2.0*params.height/params.resolution + 0.5);
	for(int i = 0; i < nocells[0]; i++){
		std::vector < std::vector < double > > w;
		density.push_back(w);
		for(int j = 0; j < nocells[1]; j++){
			std::vector <double> v;
			density[i].push_back(v);
			for(int k = 0; k < nocells[2]; k++){
				double x = i*params.resolution - params.radius;
				double y = j*params.resolution - params.radius;
				double z = k*params.resolution - params.height;
				double r = sqrt(x*x + y*y + z*z);
				density[i][j].push_back(G2(r, z, params.n2, params.A2, params.H2));

				if (density[i][j][k] != density[i][j][k])
					throw std::runtime_error("ThinDisk::ThinDisk: NaN density.");
			}
		}
	}
	std::cout << "|_ Finished setting up ThinDisk." << std::endl;
}

ThickDisk::ThickDisk(const GalaxyParameters& gparams) : params(gparams) {
	std::cout << "|_ Setting up ThickDisk..." << std::endl;
	//double gal_rmax = 18.00;
	//double gal_hmax = 2.00;
	//double gal_res = 0.10;
	//double Aa = 11; //kpc
	//double A1 = 17; //kpc
	//double H1 = 0.95; //kpc
	//double n1h1 = 0.0165; //cm-3 kpc
	//double n1 = n1h1 / 0.88; //cm-3
	double nocells[3];
	nocells[0] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[1] = (int)(2.0*params.radius/params.resolution + 0.5);
	nocells[2] = (int)(2.0*params.height/params.resolution + 0.5);

	for(int i = 0; i < nocells[0]; i++){
		std::vector<std::vector<double > > w;
		density.push_back(w);
		for(int j = 0; j < nocells[1]; j++){
			std::vector <double> v;
			density[i].push_back(v);
			for(int k = 0; k < nocells[2]; k++){
				double x = i*params.resolution - params.radius;
				double y = j*params.resolution - params.radius;
				double z = k*params.resolution - params.height;
				double r = sqrt(x*x + y*y + z*z);
				density[i][j].push_back(G1(r, z, params.n1, params.A1, params.Aa, params.H1));

				if (density[i][j][k] != density[i][j][k])
					throw std::runtime_error("ThickDisk::ThickDisk: NaN density.");
			}
		}
	}
	std::cout << "|_ Finished setting up ThickDisk." << std::endl;
}
