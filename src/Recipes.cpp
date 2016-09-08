#include <math.h>
#include <random>
#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <iomanip> //std::setprecision

#include "Recipes.hpp"
#include "Random.hpp"

double Recipes::BIG = 1.0e30;
double Recipes::ZEPS = 1.0e-10;
double Recipes::CGOLD = 0.3819660;
double Recipes::ACC = 1.0e-4;
double Recipes::GOLD = 1.618034;
double Recipes::GLIMIT = 100.0;
double Recipes::TINY = 1.0e-20;

int Recipes::nn;
std::vector<double> Recipes::xx, Recipes::yy, Recipes::sx, Recipes::sy, Recipes::ww;
double Recipes::aa, Recipes::offs;

static double sign(double a, double b) {
	return b >= 0.0 ? std::fabs(a) : -std::fabs(a);
}

static void shft(double& a, double& b, double& c, double d) {
	a = b;
	b = c;
	c = d;
}

void Recipes::avevar(const std::vector<double>& data, double& ave, double& var) {
	int n = data.size();
	double s, ep;
	ave = 0.0;
	for (int j = 0; j < n; j++)
		ave += data[j];
	ave /= n;
	ep = 0.0;
	var = 0.0;
	for (int j = 0; j < n; j++) {
		s = data[j]-(ave);
		ep += s;
		var += s*s;
	}
	var = (var-(ep*ep/n))/(n-1); //Corrected two-pass formula (NR 14.1.8).
}
double Recipes::gammln(double a) {
	//Returns the value ln[GammaFunc(a)] fo a > 0.
	double x, y, tmp, ser;
	static double cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	int j;
	x = a;
	y = a;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j < 6; j++)
		ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
void Recipes::gser(double& gamser, double a, double x, double& gln) {
	const double EPS = std::numeric_limits<double>::epsilon();
	const int ITMAX = 100;
	int n;
	double sum, del, ap;
	gln = gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) 
			std::cout << "ERROR: x less than 0 in routine gser." << std::endl;
		gamser = 0.0;
		return;
	}
	else{
		ap = a;
		sum = 1.0 / a;
		del = 1.0/a;
		for (n = 1; n <= ITMAX; n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				gamser = sum*exp(-x+a*log(x)-(gln));
				return;
			}
		}
		std::cout << "ERROR: a too large, ITMAX too small in routine gser." << std::endl;
		return;
	}
}
void Recipes::gcf(double& gammcf, double a, double x, double& gln) {
	//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf. Also returns ln Γ(a) as gln.
	const int ITMAX = 100;
	const double EPS = std::numeric_limits<double>::epsilon();
	const double FPMIN = std::numeric_limits<double>::min()/EPS;
	int i;
	double an, b, c, d, del, h;
	gln = gammln(a);
	//Set up for evaluating continued fraction by modiﬁed Lentz’s method (NR §5.2) with b0 = 0.
	b = x+1.0-a;
	c = 1.0/FPMIN;
	d = 1.0/b;
	h = d;
	//Iterate to convergence.
	for (i = 1; i <= ITMAX; i++) {
		an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = b + (an/c);
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0/d;
		del = d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS)
			break;
	}
	if (i > ITMAX) 
		std::cout << "ERROR: a too large, ITMAX too small in gcf." << std::endl;
	gammcf = exp(-x+a*log(x)-(gln))*h; //Put factors in front.
}

double Recipes::gammq(double a, double x) {
	double gamser, gammcf, gln;
	if (x < 0.0 || a <= 0.0)
		std::cout << "ERROR: Invalid arguments in routine gammq." << std::endl;
	if (x < (a+1.0)) {
		gser(gamser, a, x, gln);
		return 1.0-gamser;
	}
	else{
		gcf(gammcf, a, x, gln);
		return gammcf;
	}
}

double Recipes::brent(double ax, double bx, double cx, double (*fptr)(double), double tol, double& xmin) {
	const int ITMAX = 100;
	double a, b, d=0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;
	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	v = bx;
	w = bx;
	x = bx;
	fv = (*fptr)(v);
	fw = (*fptr)(w);
	fx = (*fptr)(x);
	for (int iter = 0; iter < ITMAX; iter++) {
		xm = 0.5*(a+b);
		tol2 = 2.0*(tol1=tol*fabs(x)+ZEPS);
		//Test for done here.
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin=x;
			return fx;
		}
		//Construct a trial parabolic fit.
		if (fabs(e) > tol1) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q-(x-w)*r;
			q = 2.0*(q-r);
			if (q > 0.0)
				p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d = CGOLD*(e=(x >= xm ? a-x : b-x));
			//The above conditions determine the acceptability of the parabolic ﬁt. Here we
			//take the golden section step into the larger of the two segments.
			else{
				d = p/q; //Take the parabolic step.
				u = x+d;
				if (u-a < tol2 || b-u < tol2)
					d = sign(tol1, xm-x);
			}
		}
		else
			d = CGOLD*(e = (x >= xm ? a-x : b-x));
		u = ((fabs(d) >= tol1) ? (x+d) : (x+sign(tol1,d)));
		fu = (*fptr)(u);
		//This is the one function evaluation per iteration.
		if (fu <= fx) {
			//Now decide what to do with our function evaluation.
			if (u >= x)
				a=x;
			else
				b=x;
			shft(v,w,x,u); //Housekeeping follows:
			shft(fv,fw,fx,fu);
		}
		else{
			if (u < x)
				a = u; 
			else
				b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		} 
		//Done with housekeeping. Back for another iteration.
	}
	std::cout << "ERROR: Too many iterations in brent." << std::endl;
	xmin = x; //Never get here.
	return fx;
}

double Recipes::zbrent(double (*fptr)(double), double x1, double x2, double tol) {
	const double EPS = std::numeric_limits<double>::epsilon();
	const int ITMAX = 100;
	int iter;
	double a=x1, b=x2, c=x2, d, e=x2-x1, min1, min2;
	double fa=(*fptr)(a), fb=(*fptr)(b), fc, p, q, r, s, tol1, xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		std::cout << "ERROR: Root must be bracketed in zbrent." << std::endl;
	fc = fb;
	for (iter = 0; iter < ITMAX; iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			d = b-a;
			e = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0*EPS*fabs(b) + 0.5*tol; //Convergence check.
		xm = 0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0)
			return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb/fa;
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0 - s;
			}
			else{
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)
				q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q - fabs(tol1*q);
			min2 = fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p/q;
			}
			else{
				d = xm;
				e = d;
			}
		}
		else{
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += sign(tol1, xm);
		fb = (*fptr)(b);
	}
	std::cout << "ERROR: Maximum number of iterations exceeded in zbrent." << std::endl;
	return 0.0;
}

//mnbrak
//fit

void Recipes::mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, double (*fptr)(double)) {
	double ulim, u, r, q, fu, dum;
	fa = (*fptr)(ax);
	fb = (*fptr)(bx);
	//Switch roles of a and b so that we can go downhill in the direction from a to b.
	if (fb > fa) {
		shft(dum,ax,bx,dum);
		shft(dum,fb,fa,dum);
	}
	cx = bx + GOLD*(bx-ax); //First guess for c.
	fc = (*fptr)(cx);
	while(fb > fc) { //Keep returning here until we bracket.
		r = (bx-ax)*(fb-fc); 
		//Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any possible division by zero.
		q = (bx-cx)*(fb-fa);
		u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(fmax(fabs(q-r), TINY), q-r));
		ulim = bx+(GLIMIT*(cx-bx));
		if ((bx-u)*(u-cx) > 0.0) {
			fu = (*fptr)(u);
			if (fu < fc) {
				ax = bx;
				bx = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb) {
				cx = u;
				fc = fu;
				return;
			}
			u = cx+(GOLD*(cx-bx));
			fu = (*fptr)(u);
		}
		else if ((cx-u)*(u-ulim) > 0.0) {
			fu = (*fptr)(u);
			if (fu < fc) {
				shft(bx,cx,u,cx+GOLD*(cx-bx));
				shft(fb,fc,fu,(*fptr)(u));
			}
		}
		else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u = ulim;
			fu = (*fptr)(u);
		}
		else{
			u = cx+(GOLD*(cx-bx));
			fu = (*fptr)(u);
		}
		shft(ax,bx,cx,u);
		shft(fa,fb,fc,fu);
	}
}

void Recipes::fit(const std::vector<double>& x, const std::vector<double>& y, int ndata, const std::vector<double>& sig, int mwt, double& a, double& b, double& siga, double& sigb, double& chi2, double& q) {
	int i;
	double wt, t, sxoss, sx=0.0, sy=0.0, st2=0.0, ss, sigdat;
	b = 0.0;
	//Accumulate sums ...
	if (mwt) {
		ss=0.0;
		//...with weights
		for (i = 0; i < ndata; i++) {
			wt = 1.0/(sig[i]*sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	}
	else{
		//...or without weights.
		for (i = 0; i < ndata; i++) {
			sx += x[i];
			sy += y[i];
		}
		ss = ndata;
	}
	sxoss = sx/ss;
	if (mwt) {
		for (i = 0; i < ndata; i++) {
			t = (x[i]-sxoss)/sig[i];
			st2 += t*t;
			b += t*y[i]/sig[i];
		}
	}
	else{
		for (i = 0; i < ndata;i++) {
			t = x[i]-sxoss;
			st2 += t*t;
			b += t*y[i];
		}
	}
	b /= st2; //Solve for a, b, σ, and σb.
	a = (sy-sx*b)/ss;
	siga = sqrt((1.0+sx*sx/(ss*st2))/ss);
	sigb = sqrt(1.0/st2);
	chi2 = 0.0; //Calculate chi^2.
	q = 1.0;
	if (mwt == 0) {
		for (i = 0; i < ndata; i++)
			chi2 += pow(y[i]-a-b*x[i], 2.0);
		//For unweighted data evaluate typical sig using chi2, and adjust the standard deviations.
		sigdat = sqrt((chi2)/(ndata-2)); 
		siga *= sigdat;
		sigb *= sigdat;
	}
	else{
		for (i = 0; i < ndata; i++)
			chi2 += pow((y[i]-a-b*x[i])/sig[i], 2.0);
		if (ndata > 2)
			q = gammq(0.5*(ndata-2),0.5*(chi2)); //Equation (NR 15.2.12).
	}
}

double Recipes::chixy(double bang) {
	int j;
	double ans, avex=0.0, avey=0.0, sumw=0.0, b;
	b = tan(bang);
	for (j = 0; j < nn; j++) {
		ww[j] = (b*sx[j])*(b*sx[j])+(sy[j])*(sy[j]);
		sumw += (ww[j] = (ww[j] < 1.0/BIG ? BIG : 1.0/ww[j]));
		avex += ww[j]*xx[j];
		avey += ww[j]*yy[j];
	}
	avex /= sumw;
	avey /= sumw;
	aa = avey-b*avex;
	for (ans = -offs, j = 0; j < nn; j++)
		ans += ww[j]*(yy[j]-aa-b*xx[j])*(yy[j]-aa-b*xx[j]);
	return ans;
}

void Recipes::fitexy(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& sigx, const std::vector<double>& sigy, double& a, double& b, double& siga, double& sigb) {
	double PI = 3.14159;
	//vector<double> xx(n), yy(n), sx(n), sy(n), ww(n);
	//int n = x.size();
	nn = x.size();
	xx.resize(nn);
	yy.resize(nn);
	sx.resize(nn);
	sy.resize(nn);
	ww.resize(nn);
	double xm, ym, xvar, yvar, scale, ch[6], chi2, r2=0.0, bmx=BIG, bmn=BIG, d1, d2, amx, amn;
	avevar(x, xm, xvar);
	avevar(y, ym, yvar);
	scale = sqrt(xvar/yvar);
	for (int i = 0; i < nn; i++) {
		xx[i] = x[i];
		yy[i] = y[i]*scale;
		sx[i] = sigx[i];
		sy[i] = sigy[i]*scale;
		ww[i] = sqrt(sx[i]*sx[i] + sy[i]*sy[i]);
	}
	double dum1, dum2, dum3, dum4, dum5;
	fit(xx, yy, nn, ww, 1, dum1, b, dum2, dum3, dum4, dum5);
	double ang[6] = {0.0, atan(b), 0.0, 0.0, atan(b), 1.571};
	offs = 0.0;
	for (int i = 3; i < 6; i++)
		ch[i] = chixy(ang[i]);
	mnbrak(ang[0], ang[1], ang[2], ch[0], ch[1], ch[2], chixy);
	chi2 = brent(ang[0], ang[1], ang[2], chixy, ACC, b);
	chi2 = chixy(b);
	a = aa;
	//q = gammq(0.5*(nn-2), chi2*0.5);
	for (int i = 0; i < nn; i++)
		r2 += ww[i];
	r2 = 1.0/r2;
	offs = chi2+1.0;
	for (int i = 0; i < nn; i++) {
		if (ch[i] > offs) {
			d1 = fabs(ang[i]-b);
			while(d1 >= PI)
				d1 -= PI;
			d2 = PI - d1;
			if (ang[i] < b) {
				double swap = d1;
				d1 = d2;
				d2 = swap;
			}
			if (d1 < bmx)
				bmx = d1;
			if (d2 < bmn)
				bmn = d2;
		}
	}
	if (bmx < BIG) {
		bmx = zbrent(chixy, b, b+bmx, ACC) - b;
		amx = aa-a;
		bmn = zbrent(chixy, b, b-bmn, ACC) - b;
		amn = aa-a;
		sigb = sqrt(0.5*(bmx*bmx + bmn*bmn))/(scale*cos(b)*cos(b));
		siga = sqrt(0.5*(amx*amx + amn*amn) + r2)/scale;
	}
	else{
		sigb = BIG;
		siga = BIG;
	}
	a /= scale;
	b = tan(b)/scale;
	xx.clear();
	yy.clear();
	sx.clear();
	sy.clear();
	ww.clear();
}

double Recipes::randomp(double p, double min, double max) {
	double p1, lo, hi, temp, norm, r, expo, x;
	p1 = p + 1.0;
	lo = min;
	hi = max;
	if (hi < lo) {
		temp = lo;
		lo = hi;
		hi = temp;
	}
	r = Random::randomDouble(0.0, 1.0);
	if (p != -1.0) {
		norm = 1.0/(pow(hi, p1) - pow(lo, p1));
		expo = log10(r/norm + pow(lo, p1))/p1;
		x = pow(10.0, expo);
	}
	else {
		norm = 1.0/(log(hi) - log(lo));
		x = exp((r/norm) + log(lo));
	}
	return x;
}

double Recipes::randomn() {
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if (iset == 0) {
		do {
			//Pick 2 uniform random numbers in square extending -1 to +1 in each direction.
			v1 = 2.0*Random::randomDouble(0.0, 1.0) - 1.0;
			v2 = 2.0*Random::randomDouble(0.0, 1.0) - 1.0;
			rsq = (v1*v1) + (v2*v2); //See if in unit circle.
		} while(rsq >= 1.0 || rsq == 0); //If not, try again.
		fac = sqrt(-2.0*log(rsq)/rsq);
		//Now make the Box-Muller transformation to get 2 normal deviates.
		//Return one and save the other for next time.
		gset = v1*fac;
		iset = 1; //Set flag.
		return v2*fac;
	} 
	else {
		//We have an extra deviate, so unset flag and return it.
		iset = 0;
		return gset;
	}
}
