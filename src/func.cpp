#include "func.hpp"
#include "Random.hpp"

std::vector<double> deriv(const std::vector<double>& y) {
	if (y.size() < 3)
		exit(EXIT_FAILURE);
	int n = y.size();
	std::vector<double> d(n, 0);
	for (int i = 1; i < n-1; i++)
		d[i] = (y[i+1] - y[i-1])/2.0;
	d[0] = (-3.0*y[0] + 4.0*y[1] - y[2])/2.0;
	d[n-1] = (3.0*y[n-1] - 4.0*y[n-2] + y[n-3])/2.0;
	return d;
}

double interpol(const std::vector<double>& f, const std::vector<double>& x, double subx) {
	int nx = x.size();
	int il=-1, ir=-1;
	for (int i = 0; i < nx; i++) {
		if (x[i] - subx > 0) {
			if (ir < 0)
				ir = i;
			else if (x[i] < x[ir])
				ir = i;
		}
		else if (x[i] - subx < 0) {
			if (il < 0)
				il = i;
			else if (x[i] > x[il])
				il = i;
		}
		else if (x[i] - subx == 0)
			return f[x[i]];
	}
	if (ir == -1 && il == -1) {
		std::cerr << "ERROR: ir & il not assigned values in interpol()." << std::endl;
		exit(EXIT_FAILURE);
	}
	else if (il == -1)
		return f[0];
	else if (ir == -1)
		return f[nx-1];
	else
		return f[il] + (subx - x[il]) * (f[ir] - f[il]) / (x[ir] - x[il]);
}

double interpol2D(const std::vector<std::vector<double> >& f, int xc, int yc) {
	double x = 1 - xc - (int)xc;
	double y = 1 - yc - (int)yc;
	double a = f[(int)xc][(int)yc];
	double b = f[(int)xc + 1][(int)yc];
	double c = f[(int)xc][(int)yc + 1];
	double d = f[(int)xc + 1][(int)yc + 1];
	return a * x * y + b * (1 - x) * y + c * x * (1 - y) + d * (1 - x) * (1 - y);
}

std::vector<double> slice(const std::vector<double>& v, unsigned int imin, unsigned int imax) {
	if (imax < imin) {
		std::cerr << "imax < imin in function slice()" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (imax > v.size()-1)
		imax = v.size()-1;
	if (imin < 0)
		imin = 0;
	std::vector<double> result(imax-imin+1, 0.0);
	for (unsigned int i = 0; i < result.size(); i++)
		result[i] = v[i+imin];
	return result;
}

std::vector<double> vinterpol(const std::vector<double>& f, const std::vector<double>& x, const std::vector<double>& subx) {
	std::vector<double> result(subx.size(), 0);
	for (unsigned int i = 0; i < result.size(); i++)
		result[i] = interpol(f, x, subx[i]);
	return result;
}

double v3min(const std::vector<std::vector<std::vector<double> > >& vec) {
	int imin = 0, jmin = 0, kmin = 0;
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec[i].size(); j++) {
			for (unsigned int k = 0; k < vec[i][j].size(); k++) {
				if (vec[i][j][k] < vec[imin][jmin][kmin]) {
					imin = i;
					jmin = j;
					kmin = k;
				}
			}
		}
	}
	return vec[imin][jmin][kmin];
}

double v3max(const std::vector<std::vector<std::vector<double> > >& vec) {
	int imax = 0, jmax = 0, kmax = 0;
	for (unsigned int i = 0; i < vec.size(); i++) {
		for (unsigned int j = 0; j < vec[i].size(); j++) {
			for (unsigned int k = 0; k < vec[i][j].size(); k++) {
				if (vec[i][j][k] > vec[imax][jmax][kmax]) {
					imax = i;
					jmax = j;
					kmax = k;
				}
			}
		}
	}
	return vec[imax][jmax][kmax];
}

double vmin(const std::vector<double>& f) {
	int imin = 0;
	for (unsigned int i = 0; i < f.size(); i++) {
		if (f[i] < f[imin])
			imin = i;
	}
	return f[imin];
}

double vmax(const std::vector<double>& f) {
	int imax = 0;
	for (unsigned int i = 0; i < f.size(); i++) {
		if (f[i] > f[imax])
			imax = i;
	}
	return f[imax];
}

double vtotal(const std::vector<double>& vec) {
	double sum = 0.0;
	for (unsigned int i = 0; i < vec.size(); i++)
		sum += vec[i];
	return sum;
}

std::vector<double> findgen(int n, double fac) {
	std::vector<double> result(n, 0);
	for (unsigned int i = 0; i < result.size(); i++)
		result[i] = i*fac;
	return result;
}

std::vector<double> cubeInterpol(const std::vector<std::vector<std::vector<double> > >& V, const std::vector<double>& xc, const std::vector<double>& yc, const std::vector<double>& zc) {
	std::vector<double> result;
	double a[2][2];
	double b[2];
	double c;
	for (unsigned int i = 0; i < xc.size(); i++) {
		int x0[3] = {(int)xc[i],(int)yc[i],(int)zc[i]};
		int x1[3] = {(int)xc[i] + 1,(int)yc[i] + 1,(int)zc[i] + 1};
		double xd[3] = {xc[i] - x0[0], yc[i] - x0[1], zc[i] - x0[2]};
		a[0][0] = V[x0[0]][x0[1]][x0[2]]*(1.0-xd[0]) + V[x1[0]][x0[1]][x0[2]]*xd[0];
		a[1][0] = V[x0[0]][x1[1]][x0[2]]*(1.0-xd[0]) + V[x1[0]][x1[1]][x0[2]]*xd[0];
		a[0][1] = V[x0[0]][x0[1]][x1[2]]*(1.0-xd[0]) + V[x1[0]][x0[1]][x1[2]]*xd[0];
		a[1][1] = V[x0[0]][x1[1]][x1[2]]*(1.0-xd[0]) + V[x1[0]][x1[1]][x1[2]]*xd[0];
		b[0] = a[0][0]*(1.0-xd[1]) + a[1][0]*xd[1];
		b[1] = a[0][1]*(1.0-xd[1]) + a[1][1]*xd[1];
		c = b[0]*(1.0-xd[2]) + b[1]*xd[2];
		if (c < 0 || c != c) {
			std::cout << std::endl;
			std::cout << xc[i] << '\t' << yc[i] << '\t' << zc[i] << std::endl;
			std::cout << x0[0] << '\t' << x0[1] << '\t' << x0[2] << std::endl;
			std::cout << x1[0] << '\t' << x1[1] << '\t' << x1[2] << std::endl;
			std::cout << xd[0] << '\t' << xd[1] << '\t' << xd[2] << std::endl;
			std::cout << a[0][0] << '\t' << a[1][0] << '\t' << a[0][1] << '\t' << a[1][1] << '\t' << b[0] << '\t' << b[1] << std::endl;
			exit(EXIT_FAILURE);
		}
		result.push_back(c);
	}
	return result;
}

std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn) {
	int n = x.size();
	std::vector<double> y2(n, 0.0);
	std::vector<double> u(n-1, 0.0);
	int i,k;
	double p,qn,sig,un;
	//The lower boundary condition is set either to be "natural"
	//or else to have a specified first derivative.
	if (yp1 > 0.99e30) {
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else{
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((f[1]-f[0])/(x[1]-x[0])-yp1);
	}
	//This is the decomposition loop of the tridiagonal algorithm. 
	//y2 and u are used for temporary storage of the decomposed factors.
	for (i = 1; i < n-1; i++) {
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (f[i+1]-f[i])/(x[i+1]-x[i]) - (f[i]-f[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	//The upper boundary condition is set either to be "natural"
	// or else to have a specified first derivative.
	if (ypn > 0.99e30) { 
		qn = 0.0;
		un = 0.0;
	}
	else{
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(f[n-1]-f[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	//This is the backsubstitution loop of the tridiagonal algorithm.
	for (k = n-2; k >= 0; k--)
		y2[k] = y2[k]*y2[k+1]+u[k];
	return y2;
}
double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2) {
	int klo, khi, k;
	double h,b,a;
	int n = x.size();
	//We will find the right place in the table by means of bisection. This is 
	//optimal if sequential calls to this routine are at random values of x. 
	//If sequential calls are in order, and closely spaced, one would do better
	//to store previous values of klo and khi and test if they remain appropriate 
	//on the next call.
	/*
	cerr << "x = {" << x[0];
	for (unsigned int i = 1; i < x.size(); i++)
		cerr << ", " << x[i];
	cerr << "}" << std::endl;
	cerr << "x2 = " << x2 << std::endl;
	 */
	klo = 0;
	khi = n-1;
	while(khi-klo > 1) {
		k = (khi+klo) >> 1;
		if (x[k] > x2)
			khi = k;
		else
			klo = k;
	}
	//cerr << "klo = " << klo << '\t' << "khi = " << khi << std::endl;
	//klo and khi now bracket the input value of x.
	h = x[khi]-x[klo];
	if (h == 0.0) {
		std::cerr << "ERROR: Bad x input to routine splint()" << std::endl; //The xa's must be distinct.
		for (int i = 0; i < n; i++)
			std::cerr << "x[" << i << "] = " << x[i];
		exit(EXIT_FAILURE);
	}
	a = (x[khi]-x2)/h;
	b = (x2-x[klo])/h; //Cubic spline polynomial is now evaluated.

	return a*f[klo]+b*f[khi]+((a*a*a-a)*f2[klo]+(b*b*b-b)*f2[khi])*(h*h)/6.0;
}

std::vector<double> vsplint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, const std::vector<double>& x2) {
	std::vector<double> result(x2.size(), 0.0);
	for (unsigned int i = 0; i < x2.size(); i++)
		result[i] = splint(x, f, f2, x2[i]);
	return result;
}

double int_tabulated(const std::vector<double>& x, const std::vector<double>& f) {
	double result = 0;
	int xsegments = x.size() - 1;
	while(xsegments%4 != 0)
		xsegments++;
	double xmin = vmin(x);
	double xmax = vmax(x);
	//uniform step size.
	double h = (xmax-xmin)/xsegments;
	//calculate interpolating cubic spline.
	std::vector<double> f2 = spline(x, f, 0, 0);
	//define x values at which interpolated y values are desired.
	std::vector<double> x2 = findgen(xsegments+1, h);
	for (int i = 0; i < xsegments+1; i++)
		x2[i] += xmin;
	//compute the interpolates.
	std::vector<double> z = vsplint(x, f, f2, x2);
	//compute integral using 5-point Newton-Cotes formula.
	for (int i = 4; i <= xsegments; i += 4) {
		result += (2.0*h/45.0)*(7.0*(z[i-4]+z[i]) + 32.0*(z[i-3]+z[i-1]) + 12.0*z[i-2]);
	}
	return result;
}

int round_int(double x) {
	return (int)(x+0.5);
}

double round_dp(double x, int dp) {
	x = x*pow(10, dp);
	x = round_int(x);
	return x/pow(10, dp);
}

std::vector<double> xhist(double binmin, double binmax, int nbins) {
	std::vector<double> xarr(nbins, 0.0);
	double binsize = (binmax-binmin)/nbins;
	for (int i = 0; i < nbins; i++)
		xarr[i] = binmin + (i*binsize) + 0.5*binsize;
	return xarr;
}

std::vector<int> yhist(const std::vector<double>& y, const std::vector<double>& x) {
	int nbins = x.size();
	int ndata = y.size();
	if (nbins >= 2) {
		std::vector<int> histarr(nbins, 0);
		double binsize = x[1]-x[0];
		std::vector<int> yarr(nbins, 0);
		for (int i = 0; i < nbins; i++) {
			for (int j = 0; j < ndata; j++) {
				if (y[j] >= x[i] - (binsize/2.0) && y[j] < x[i] + (binsize/2.0))
					histarr[i]++;
			}
		}
		return histarr;
	}
	else{
		std::cout << "ERROR: Invalid xhist array passed to function yhist()." << std::endl;
		std::vector<int> histarr(1, 0);
		return histarr;
	}
}
/*
std::vector<double>::iterator remove_if_lt (std::vector<double>::iterator first, std::vector<double>::iterator last, double value) {
	std::vector<double>::iterator result = first;
	while(first != last) {
		if (*first < value) {
 *result = *first;
			++result;
		}
		++first;
	}
	return result;
}
 */
void remove_if_lt(std::vector<double>& vec, double value) {
	std::vector<double>::iterator first=vec.begin(), last=vec.end(), result=vec.begin();
	while(first != last) {
		if (*first < value) {
			*result = *first;
			++result;
		}
		++first;
	}
	vec.erase(result, last);
}

int binarySearch(const std::vector<double>& vec, double val) {
	int lo = 0;
	int hi = vec.size();
	while (lo != hi) {
		int mid = lo + (hi-lo)/2;
		if (val == vec[mid])
			return mid;
		else if (val > vec[mid])
			lo = mid+1;
		else
			hi = mid;
	}
	return lo;
}

std::vector<int> random_loc(const std::vector<std::vector<int> >& vec2, const std::vector<std::vector<std::vector<double> > >& vec3, double val, std::default_random_engine& engine) {
	int num = vec2.size();
	int r = (int)(Random::randomDouble(0.0, 1.0)*(num-1) + 0.5);
	int minmax[4] = {0,0,0,0};
	double r2 = Random::randomDouble(0.0, 1.0);
	if (r2 > 0.5) {
		minmax[0] = 0;
		minmax[1] = r;
		minmax[3] = r;
		minmax[4] = num;
	}
	else{
		minmax[0] = r;
		minmax[1] = num;
		minmax[3] = 0;
		minmax[4] = r;
	}
	for (int i = minmax[0]; i < minmax[1]; i++) {
		if (vec3[vec2[i][0]][vec2[i][1]][vec2[i][2]] > val)
			return vec2[i];
	}
	for (int i = minmax[2]; i < minmax[3]; i++) {
		if (vec3[vec2[i][0]][vec2[i][1]][vec2[i][2]] > val)
			return vec2[i];
	}
	std::vector<int> no_result(3, -1);
	return no_result;
}
std::vector<int> random_loc2(const std::vector<std::vector<std::vector<double> > >& vec3, double val, std::default_random_engine& engine) {
	int ni = vec3.size(), nj = vec3[0].size(), nk = vec3[0][0].size();
	int num = ni*nj*nk;
	std::vector<int> randi(ni*nj*nk, 0);
	for (int i = 0; i < ni*nj*nk; i++)
		randi[i] = i;
	random_shuffle(randi.begin(), randi.end());
	for (int i = 0; i < num; i++) {
		int li = (int)(randi[i]/(nj*nk));
		int lj = (int)((randi[i] - (li*nj*nk))/nk);
		int lk = (int)(randi[i] - (lj*nk) - (li*nj*nk));
		if (vec3[li][lj][lk] > val) {
			std::vector<int> result(3, 0);
			result[0] = li;
			result[1] = lj;
			result[2] = lk;
			return result;
		}
	}
	std::vector<int> no_result(3, -1);
	return no_result;
}
std::vector<int> rand_loc_gt(const std::vector<std::vector<std::vector<double> > >& vec3, double value, std::default_random_engine& engine) {
	int xsize = vec3.size(), ysize = vec3[0].size(), zsize = vec3[0][0].size();
	std::vector<int> result(3, -1);
	int ranx = (int)((xsize-1)*Random::randomDouble(0.0, 1.0));
	int rany = (int)((ysize-1)*Random::randomDouble(0.0, 1.0));
	int ranz = (int)((zsize-1)*Random::randomDouble(0.0, 1.0));
	for (int i = ranx; i < xsize; i++) {
		for (int j = rany; j < ysize; j++) {
			for (int k = ranz; k < zsize; k++) {
				if (vec3[i][j][k] > value) {
					result[0] = i;
					result[1] = j;
					result[2] = k;
					return result;
				}
			}
		}
	}
	for (int i = std::max(ranx-1, 0); i >= 0; i--) {
		for (int j = std::max(rany-1, 0); j >= 0; j--) {
			for (int k = std::max(ranz-1, 0); k >= 0; k--) {
				if (vec3[i][j][k] > value) {
					result[0] = i;
					result[1] = j;
					result[2] = k;
					return result;
				}
			}
		}
	}
	return result;
}
void printScatter(const std::vector<double>& vec1, const std::vector<double>& vec2, std::string filename) {
	// creating filename
	std::ofstream ofile(filename);
	if (!ofile)
		std::cerr << "ERROR: unable to open " << filename << std::endl;
	// writing data to file
	ofile << std::setprecision(6) << std::fixed;
	for (unsigned int i = 0; i < vec1.size(); i++) {
		if (i < vec2.size())
			ofile << vec1[i] << '\t' << vec2[i] << std::endl;
		else
			ofile << vec1[i] << std::endl;
	}
	ofile.close();
}

double planck(double wave, double temp) {
	double w, c1, c2, val;
	w = wave/1.0e8; // A to cm
	c1 = 3.7417749e-5; // = 2*h*c*c
	c2 = 1.4387687; // = h*c/k
	val = c2/(w*temp);
	if (val < std::log(std::numeric_limits<double>::max()))
		return c1 / (pow(w, 5.0)*(exp(val)-1.0));
	else
		return 0.0;
}
double planck2(double freq, double temp) {
	double c1, c2, val;
	c1 = 1.474499155e-50; // = 2*h/c*c
	c2 = 4.799237548e-11; // = h/k
	val = c2*freq/temp;
	if (val < std::log(std::numeric_limits<double>::max()))
		return c1*pow(freq, 3.0) / (exp(val)-1.0);
	else
		return 0.0;
}


double vmean(const std::vector<double>& vec) {
	double result = 0;
	int n = vec.size();
	for (int i = 0; i < n; i++)
		result += vec[i];
	return result/(double)(n);
}

double vstddev(const std::vector<double>& vec) {
	double mu = vmean(vec);
	double result = 0;
	int n = vec.size();
	for (int i = 0; i < n; i++) {
		result += (vec[i]-mu)*(vec[i]-mu);
	}
	return sqrt(result/(double)(n));
}

double random_kroupa(std::default_random_engine& engine) {
	double p[4] = {-0.3, -1.8, -2.7, -2.35};
	double m[4][2] = {{0.01, 0.08}, {0.08, 0.5}, {0.5, 1.0}, {1.0, 120.0}};
	double sumtotal = 0;
	double norms[4];

	for (int i = 0; i < 4; i++) {
		if (p[i] != -1.0)
			norms[i] = (p[i]+1.0)/(std::pow(m[i][1], p[i]+1.0) - std::pow(m[i][0], p[i]+1.0));
		else
			norms[i] = 1.0/(std::log(m[i][1]) - std::log(m[i][0]));
		sumtotal += 1.0/norms[i];
	}

	double probs[4] = {0.0, 0.0, 0.0, 0.0};

	for (int i = 0; i < 4; i++)
		probs[i] = 1.0/(norms[i]*sumtotal);

	double ben_norm = (1.0/20.0) + (1.0/11.0) + (1.0/1.1) + (1.0/1.0);
	probs[3] = (1.0/20.0)/ben_norm;
	probs[2] = (1.0/11.0)/ben_norm;
	probs[1] = (1.0/1.1)/ben_norm;
	probs[0] = (1.0/1.0)/ben_norm;

	for (int i = 1; i < 4; i++)
		probs[i] += probs[i-1];

	double r = Random::randomDouble(0.0, 1.0);
	int seg = 3;

	for (int i = 0; i < 4; i++) {
		if (r < probs[i]) {
			seg = i;
			break;
		}
	}

	return Random::randomPowerlaw(p[seg], m[seg][0], m[seg][1]);
}
