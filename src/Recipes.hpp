#ifndef RECIPES_HPP
#define RECIPES_HPP

#include <vector> //std::vector
#include <string>

class Recipes {
public:
	/* avevar:
	Given array data[1..n], returns its mean as ave and its variance as var.
	(Dependencies: none)
	*/
	static void avevar(const std::vector<double>& data, double& ave, double& var);
	/* gammln:
	Returns the value ln[Γ(xx)] for xx > 0.
	(Dependencies: none)
	*/
	static double gammln(double a);
	/* gser:
	Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
	Also returns ln Γ(a) as gln.
	(Dependencies: gammln)
	*/
	static void gser(double& gamser, double a, double x, double& gln);
	/* gcf:
	Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf. Also returns ln Γ(a) as gln.
	(Dependencies: gammln)
	*/
	static void gcf(double& gammcf, double a, double x, double& gln);
	/* gammq:
	Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
	(Dependencies: gcf, gser)
	*/
	static double gammq(double a, double x);
	/* brent:
	Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
	between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
	the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
	the minimum is returned as xmin, and the minimum function value is returned as brent, the
	returned function value.
	(Dependencies: none)
	*/
	static double brent(double ax, double bx, double cx, double (*fptr)(double), double tol, double& xmin);
	/* zbrent:
	Using Brent’s method, ﬁnd the root of a function func known to lie between x1 and x2. The
	root, returned as zbrent, will be reﬁned until its accuracy is tol
	(Dependencies: none)
	*/
	static double zbrent(double (*fptr)(double), double x1, double x2, double tol);
	/* mnbrak:
	Given a function func, and given distinct initial points ax and bx, this routine searches in
	the downhill direction (deﬁned by the function as evaluated at the initial points) and returns
	new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
	values at the three points, fa, fb, and fc.
	(Dependencies: none)
	*/
	static void mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, double (*fptr)(double));
	/* fit:
	Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
	sig[1..ndata], ﬁt them to a straight line y = a + bx by minimizing χ2. Returned are
	a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and the
	goodness-of-ﬁt probability q (that the ﬁt would have χ2 this large or larger). If mwt=0 on
	input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
	the normalization of chi2 is to unit standard deviation on all points.
	(Dependencies: gammq)
	*/
	static void fit(const std::vector<double>& x, const std::vector<double>& y, int ndata, const std::vector<double>& sig, int mwt, double& a, double& b, double& siga, double& sigb, double& chi2, double& q);
	/* chixy:
	Captive function of fitexy, returns the value of (χ2 − offs) for the slope b=tan(bang).
	Scaled data and offs are communicated via the global variables.
	(Dependencies: none)
	*/
	static double chixy(double bang);
	/* fitexy:
	Straight-line ﬁt to input data x[1..ndat] and y[1..ndat] with errors in both x and y, the re-
	spective standard deviations being the input quantities sigx[1..ndat] and sigy[1..ndat].
	Output quantities are a and b such that y = a + bx minimizes χ2, whose value is returned
	as chi2. The χ2 probability is returned as q, a small value indicating a poor ﬁt (sometimes
	indicating underestimated errors). Standard errors on a and b are returned as siga and sigb.
	These are not meaningful if either (i) the ﬁt is poor, or (ii) b is so large that the data are
	consistent with a vertical (infinite b) line. If siga and sigb are returned as BIG, then the data
	are consistent with all values of b.
	(Dependencies: avevar, gammq, brent, zbrent, mnbrak, chixy)
	*/
	static void fitexy(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& sigx, const std::vector<double>& sigy, double& a, double& b, double& siga, double& sigb);
	/* randomp:
	Returns a random number distributed as a powerlaw.
	(Dependencies: randomu)
	*/
	static double randomp(double p, double min, double max);
	/* randomn:
	Returns a normally distributed deviate with zero mean and unit variance, using randomu(engine)
	as source of uniform deviates.
	(Dependencies: randomu)
	*/
	static double randomn();
private:
	static double BIG;
	//#define ITMAX 100
	//#define EPS 3.0e-8
	//#define EPS2 3.0e-7
	static double ZEPS;
	static double CGOLD;
	//#define FPMIN 1.0e-30
	static double ACC;
	static double GOLD;
	static double GLIMIT;
	static double TINY;

	static int nn;
	static std::vector<double> xx, yy, sx, sy, ww;
	static double aa, offs;
};

#endif
