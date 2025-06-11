// VBMicrolensing v5.0 (2025)
//
// This code has been developed by Valerio Bozza (University of Salerno) and collaborators.
// Check the repository at https://github.com/valboz/VBMicrolensing
// for the newest version.
// Any use of the code for scientific publications should be acknowledged by a citation
// to the appropriate publication, as detailed in the repository page.
//
// The code relies on the root solving algorithm by Jan Skworon and Andy Gould
// described in Skowron & Gould arXiv:1203.1034.
// Please also cite this paper if specifically relevant in your scientific publication.
// The original Fortran code available on http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
// has been translated to C++ by Tyler M. Heintz and Ava R. Hoag (2017)
//
// The Multipoly method for the calculation of multiple-lens microlensing was developed
// with Vito Saggese (2024).
// 
// Optimizations for complex functions and for high-mag regime have been developed 
// by Jiyuan Zhang (2025).
//
// GNU Lesser General Public License applies to all parts of this code.
// Please read the separate LICENSE.txt file for more details.


#ifndef __multilens
#define __multilens
#define __unmanaged

#define _L_1 x1-((x1+a/2.0)/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*(x1-a/2.0)/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q) // Used in PlotCrits
#define _L_2 x2-(x2/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*x2/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q)
#define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test
#define _J1c coefs[21]/((zc-coefs[20])*(zc-coefs[20]))+coefs[22]/(zc*zc) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J2 -2.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z))
#define _J3 6.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z*z))
#define _skew(p1,p2,q1,q2) p1*q2-p2*q1
#define _NP 200.0
#define __rsize_ESPL 151
#define __zsize_ESPL 101

#define _sign(x) ((x>0)? +1 : -1)

#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <random>
class _sols_for_skiplist_curve;
class _skiplist_curve;
/*******************************************   end   *******************************************/

class _curve;
class _sols;
class _theta;
class complex;
struct annulus;


class complex {
public:
	double re;
	double im;
	complex(double a, double b) { re = a; im = b; }
	complex(double a) { re = a; im = 0; }
	complex(void) { re = 0; im = 0; }
};



class VBMicrolensing
{
	const int maxiter = 10000;
	const double pseudorandom[12] = { 0.5, 0.4, 0.3, 0.7, 0.6, 0.5, 0.4, 0.7, 0.5, 0.6, 0.3,
		0.3 };

	int* ndatasat;
	double** tsat, *** possat;
	double Mag0;
	double* dist_mp, * q;
	int nim0, n, n2, nnm1, nroots, nrootsmp, * nrootsmp_mp;
	complex* zr, * zcr, ** pmza, ** pyaza, ** ppmy, * pza, * pza2, ** pmza2, * pdum, * ppy, * a, * s_offset, * pert, y, yc, * s;
	complex* y_mp, *** pmza_mp, ** pza_mp, *** pyaza_mp, *** ppmy_mp, ** ppy_mp, ** zr_mp;
	complex* zaltc, * J1, * J1c, ** za, ** za2;
	complex* coefs, ** coefs_mp;
	complex** a_mp, * s_sort;

	double* prodevs, * errs, err, L0f, Jacf;
	complex* devs, * init, * centralimages, * newseeds, * grads, zf, S2f, * S2s, * S3s, * S4s;
	int lencentralimages, lennewseeds, ngoodold, ngood, iter, iter2;
	////
	double* good, * Jacs, rho, rho2, * m;
	double** m_mp, * q_sort;
	int* worst;
	double e, phi, phip, phi0, Om, inc, u0, tE_inv, t0, alpha, pai1, pai2, PosAng, dPosAng, thetaE, d3, v3, GM, flagits;
	int iastro;
	double Obj[3], rad[3], tang[3], t0old;
	double Eq2000[3], Quad2000[3], North2000[3];
	double Et0[2], vt0[2], Et[2], Ehel[2];
	double ESPLout[__rsize_ESPL][__zsize_ESPL], ESPLin[__rsize_ESPL][__zsize_ESPL], ESPLoutastro[__rsize_ESPL][__zsize_ESPL], ESPLinastro[__rsize_ESPL][__zsize_ESPL];
	bool ESPLoff, multidark;
	double* LDtab, * rCLDtab, * CLDtab;
	double scr2, sscr2;
	int npLD;
	annulus* annlist;
	_skiplist_curve** cprec, ** cpres, ** cfoll;
	double** A;

	void ComputeCentroids(double* pr, double t, double* c1s, double* c2s, double* c1l, double* c2l);
	void ComputeParallax(double, double);
	double LDprofile(double r);
	double rCLDprofile(double tc, annulus*, annulus*);
	void initroot();
	int froot(complex);
	bool checkroot(_theta*);

	void SetLensGeometry_spnp(int n, double* q, complex* s);
	void SetLensGeometry_multipoly(int n, double* q, complex* s);
	void initrootpoly();
	_curve* NewImages(complex, complex*, _theta*);
	_curve* NewImages(_theta*);
	_curve* NewImagespoly(_theta*);
	_curve* NewImagesmultipoly(_theta*);
	double BinaryMagSafe(double s, double q, double y1, double y2, double rho, _sols_for_skiplist_curve** images);
	double MultiMagSafe(double y1, double y2, double rho, _sols_for_skiplist_curve** images);
	void OrderImages(_sols_for_skiplist_curve*, _curve*);

	void OrderMultipleImages(_sols_for_skiplist_curve*, _curve*);
	void cmplx_laguerre(complex*, int, complex*, int&, bool&);
	void cmplx_newton_spec(complex*, int, complex*, int&, bool&);
	void cmplx_laguerre2newton(complex*, int, complex*, int&, bool&, int);
	void solve_quadratic_eq(complex&, complex&, complex*);
	void solve_cubic_eq(complex&, complex&, complex&, complex*);
	void polyproduct(complex* p1, int n1, complex* p2, int n2, complex* pdest);
	void copypol(complex* p1, int n1, complex* pdest);
	void change_n(int nn);
	void change_n_mp(int nn);
	void polycoefficients();
	void polycoefficients_multipoly();
	void polycritcoefficients(complex eiphi);

public:
	
	double rootaccuracy;
	double samplingfactor;
	bool squarecheck;
	bool astrometry;
	bool turn_off_secondary_source;
	bool turn_off_secondary_lens;

	static char ESPLtablefile[1024];
	static void SetESPLtablefile(char* instring) { strcpy(ESPLtablefile, instring); }
	double Tol, RelTol, a1, a2, t0_par,corrquad, corrquad2, safedist;
	double mass_radius_exponent, mass_luminosity_exponent, lens_mass_luminosity_exponent;
	int satellite, parallaxsystem, t0_par_fixed, nsat;
	int minannuli, maxannuli, nannuli, NPS, NPcrit;
	int newtonstep;
	double y_1, y_2, av, therr, astrox1, astrox2;
	double (*CumulativeFunction)(double r, double* LDpars);

	// Critical curves and caustics calculation
	_sols* PlotCrit();
	_sols* PlotCrit(double a, double q);
	// Initialization for parallax calculation
	void SetObjectCoordinates(char* Coordinates_file, char* Directory_for_satellite_tables);
	void SetObjectCoordinates(char* CoordinateString);
	// Skowron & Gould root calculation
	void cmplx_roots_gen(complex*, complex*, int, bool, bool);
	void cmplx_roots_multigen(complex*, complex**, int, bool, bool);
	// Bozza optimization
	int findimagepoly(int iroot);
	int findimagemultipoly(int iroot);

	// Set Lens Geometry
	void SetLensGeometry(int n, double* q, complex* s);
	void SetLensGeometry(int n, double* pr);

	// Magnification calculation functions.

	double BinaryMag0(double s, double q, double y1, double y2, _sols_for_skiplist_curve** Images);
	double BinaryMag0(double s, double q, double y1, double y2);

	double BinaryMag(double s, double q, double y1, double y2, double rho, double accuracy, _sols_for_skiplist_curve** Images);
	double BinaryMag(double s, double q, double y1, double y2, double rho, double accuracy);
	double BinaryMag2(double s, double q, double y1, double y2, double rho);
	double BinaryMagDark(double s, double q, double y1, double y2, double rho, double accuracy);
	void BinaryMagMultiDark(double s, double q, double y1, double y2, double rho, double* a1_list, int n_filters, double* mag_list, double accuracy);
	double MultiMag0(double y1, double y2, _sols_for_skiplist_curve** Images);
	double MultiMag0(double y1, double y2);
	double MultiMag(double y1, double y2, double rho, double accuracy, _sols_for_skiplist_curve** Images);
	double MultiMag(double y1, double y2, double rho, double accuracy);
	double MultiMag(double y1, double y2, double rho);
	double MultiMag2(double y1, double y2, double rho);
	double MultiMagDark(double y1, double y2, double rho, double accuracy);

	// Limb Darkening control
	enum LDprofiles { LDlinear, LDquadratic, LDsquareroot, LDlog, LDuser };
	void SetLDprofile(double(*UserLDprofile)(double), int tablesampling);
	void SetLDprofile(LDprofiles);

	// Method control
	enum class Method { Singlepoly, Multipoly, Nopoly };
	void SetMethod(Method);

	//ESPL functions
	void LoadESPLTable(char* tablefilename);
	double ESPLMag(double u, double rho);
	double ESPLMag2(double u, double rho);
	double ESPLMagDark(double u, double rho);
	double PSPLMag(double u);


	// New (v2) light curve functions, operating on arrays

	void PSPLLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void PSPLLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void ESPLLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void ESPLLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);

	void BinaryLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveW(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);
	void BinaryLightCurveKepler(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);

	void BinSourceLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);
	void BinSourceExtLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceExtLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceSingleLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceBinLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);

	void TripleLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void TripleLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void LightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np, int nl);

	// Old (v1) light curve functions, for a single calculation
	double PSPLLightCurve(double* parameters, double t);
	double PSPLLightCurveParallax(double* parameters, double t);
	double ESPLLightCurve(double* parameters, double t);
	double ESPLLightCurveParallax(double* parameters, double t);

	double BinaryLightCurve(double* parameters, double t);
	double BinaryLightCurveW(double* parameters, double t);
	double BinaryLightCurveParallax(double* parameters, double t);
	double BinaryLightCurveOrbital(double* parameters, double t);
	double BinaryLightCurveKepler(double* parameters, double t);

	double BinSourceLightCurve(double* parameters, double t);
	double BinSourceLightCurveParallax(double* parameters, double t);
	double BinSourceLightCurveXallarap(double* parameters, double t);
	double BinSourceExtLightCurve(double* parameters, double t);
	double BinSourceExtLightCurveXallarap(double* parameters, double t);
	double BinSourceBinLensXallarap(double* parameters, double t);
	double BinSourceSingleLensXallarap(double* parameters, double t);
	double BinSourceBinLensPOX(double* parameters, double t);

	double TripleLightCurve(double* parameters, double t);
	double TripleLightCurveParallax(double* parameters, double t);

	// Astrometric functions
	void CombineCentroids(double* mags, double* c1s, double* c2s, double* c1l, double* c2l, double* c1ltot, double* c2tot, double g, int np);
	void PSPLAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void ESPLAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void BinaryAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void BinaryAstroLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* seps_array, int np);
	void BinaryAstroLightCurveKepler(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* seps_array, int np);
	void BinSourceAstroLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void TripleAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);

	// Constructor and destructor

	VBMicrolensing();
	~VBMicrolensing();

private: // Must be declared here at the end
	LDprofiles curLDprofile;
	Method SelectedMethod;

};

//std::string VBMicrolensing::ESPLtablefile = "";

double VBDefaultCumulativeFunction(double r, double* a1);

struct annulus {
	double bin;
	double cum;
	double Mag;
	double err;
	double f;
	int nim;
	double LDastrox1, LDastrox2;
	annulus* prev, * next;
};



class _theta {
public:
	double th, maxerr, Mag, errworst, astrox1, astrox2;
	int imlength;
	_theta* prev, * next;

	_theta(double);

};

class _thetas {
public:
	_theta* first, * last;
	int length;

	_thetas(void);
	~_thetas(void);
	_theta* insert(double);
	_theta* insert_at_certain_position(_theta*, double);
	// method: this method can only be used when inserting an element in the middle of linked list
	// 		   i.e. *first's 'th' < current 'th' < *last's 'th'
	// 		   and the new element is forced to be inserted between itheta and itheta->next, 
	// 		   which means it's the programmer's responsibility to guarantee itheta->th < th < itheta->next->th holds
	//         (O(1) complexity)
	void remove(_theta*);

};

#define max_skiplist_level 2

class _point {
public:
	double x1;
	double x2;
	double parab, ds, dJ, Mag, err, parabastrox1;
	complex d;								  	  // d is z'(theta) at this point
	_theta* theta;							  // pointer
	_point* next, * prev;                       // pointers that point to _point variable
	_point* next_array[max_skiplist_level + 1];

	double parabastrox2;
	_point(double, double, _theta*);
	double operator-(_point);
};

class _curve {
public:
	int length;
	_point* first, * last;
	_curve* next, * prev;
	_curve* partneratstart, * partneratend;
	double parabstart, Magstart, errstart, parabastrox1, parabastrox2;

	_curve(_point*);
	_curve(void);
	~_curve(void);

	_curve* divide(_point*);
	void drop(_point*);
	void append(double, double);
	void append(_point*);
	void prepend(double, double);
	//	void prepend(_point *);
	_curve* join(_curve*);
	_curve* joinbefore(_curve*);
	_curve* reverse(void);
	double closest(_point*, _point**);
	double closest2(_point*, _point**);
	void complement(_point**, int, _point**, int);
};

class _sols {
public:
	int length;
	_curve* first, * last;

	_sols(void);
	~_sols(void);
	void drop(_curve*);
	void append(_curve*);
	void prepend(_curve*);
	void join(_sols*);
};


class _skiplist_curve {							 // a _skiplist_curve class variable is a skip list of _point variables
public:
	_point* first, * last;

	_point* head;
	_point* last_array[max_skiplist_level + 1];
	int Level;

	int length_notation;

	_skiplist_curve* next, * prev;
	_skiplist_curve* partneratstart, * partneratend;
	double parabstart, Magstart, errstart, parabastrox1, parabastrox2;

	_skiplist_curve(_point* p1, int new_Level);
	_skiplist_curve(void);
	~_skiplist_curve(void);
	_skiplist_curve* join(_skiplist_curve* new_curve);
	void append(_point* pp, int append_Level);
	void append(double x1, double x2, int append_Level);
	_skiplist_curve* find_prev_then_divide(double th);

};




class _sols_for_skiplist_curve {                 //       a  _sols_for_skiplist_curve class variable is a linked list of _skiplist_curve variables, 
	// while a _skiplist_curve class variable           is a skip list   of _point variables
public:
	int length;
	_skiplist_curve* first, * last;

	_sols_for_skiplist_curve(void);
	~_sols_for_skiplist_curve(void);
	void drop(_skiplist_curve* ref);
	void append(_skiplist_curve* cc);

};

#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 101

//////////////////////////////
//////////////////////////////
////////complex methods and operators
//////////////////////////////
//////////////////////////////

inline double abs2(complex z) {
	return (z.re * z.re + z.im * z.im);
}

inline double abs(complex z) {
	return sqrt(z.re * z.re + z.im * z.im);
}

inline complex conj(complex z) {
	return complex(z.re, -z.im);
}

inline complex sqrt(complex z) {
	double md = sqrt(z.re * z.re + z.im * z.im);
	return (md > 0) ? complex(sqrt((md + z.re) / 2), (sqrt((md - z.re) / 2) * ((z.im > 0) ? 1 : -1))) : 0.0;
}



inline double real(complex z) {
	return z.re;
}

inline double imag(complex z) {
	return z.im;
}

inline complex operator+(complex p1, complex p2) {
	return complex(p1.re + p2.re, p1.im + p2.im);
}

inline complex operator-(complex p1, complex p2) {
	return complex(p1.re - p2.re, p1.im - p2.im);
}

inline complex operator*(complex p1, complex p2) {
	return complex(p1.re * p2.re - p1.im * p2.im, p1.re * p2.im + p1.im * p2.re);
}

inline complex operator/(complex p1, complex p2) {
	double md = p2.re * p2.re + p2.im * p2.im;
	return complex((p1.re * p2.re + p1.im * p2.im) / md, (p1.im * p2.re - p1.re * p2.im) / md);
}

inline complex operator+(complex z, double a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, double a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, double a) {
	return complex(z.re * a, z.im * a);
}

inline complex operator/(complex z, double a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(double a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(double a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(double a, complex z) {
	return complex(a * z.re, a * z.im);
}

inline complex operator/(double a, complex z) {
	double md = z.re * z.re + z.im * z.im;
	return complex(a * z.re / md, -a * z.im / md);
}

inline complex operator+(complex z, int a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, int a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, int a) {
	return complex(z.re * a, z.im * a);
}

inline complex operator/(complex z, int a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(int a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(int a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(int a, complex z) {
	return complex(a * z.re, a * z.im);
}

inline complex operator/(int a, complex z) {
	double md = z.re * z.re + z.im * z.im;
	return complex(a * z.re / md, -a * z.im / md);
}

inline complex operator-(complex z) {
	return complex(-z.re, -z.im);
}

inline bool operator==(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

inline bool operator!=(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

inline complex expcmplx(complex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return complex(r * cos(theta), r * sin(theta));
}

inline complex cbrt(complex z) {
	complex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 0.333333333333);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	return 	complex(r_cube * cos(theta_cube), r_cube * sin(theta_cube));
}


inline _point::_point(double x, double y, _theta* theta1) {
	x1 = x;
	x2 = y;
	theta = theta1;
	next = 0;
	prev = 0;
	for (int i = 0; i < (max_skiplist_level + 1); i++)
	{
		next_array[i] = 0;
	}
}

inline double _point::operator-(_point p2) {
	static double dx1, dx2;
	dx1 = x1 - p2.x1;
	dx2 = x2 - p2.x2;
	return dx1 * dx1 + dx2 * dx2;
}


#endif
