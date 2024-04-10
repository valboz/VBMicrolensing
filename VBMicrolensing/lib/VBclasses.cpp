
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>


//////////////////////////////
//////////////////////////////
////////_point methods
//////////////////////////////
//////////////////////////////


_point::_point(double x, double y, _theta *theta1) {
	x1 = x;
	x2 = y;
	theta = theta1;
}

double _point::operator-(_point p2) {
	static double dx1, dx2;
	dx1 = x1 - p2.x1;
	dx2 = x2 - p2.x2;
	return dx1 * dx1 + dx2 * dx2;
}

//////////////////////////////
//////////////////////////////
////////_curve methods
//////////////////////////////
//////////////////////////////

_curve::_curve(void) {
	length = 0;
	first = last = 0;
	partneratstart = partneratend = this;
}

_curve::_curve(_point *p1) {
	length = 1;
	first = last = p1;
	p1->prev = p1->next = 0;
	partneratstart = partneratend = this;
}

_curve::~_curve(void) {
	_point *scan1, *scan2;
	scan1 = first;
	for (int i = 0; i<length; i++) {
		scan2 = scan1->next;
		delete scan1;
		scan1 = scan2;
	}
}

_curve *_curve::divide(_point *ref) {
	_point *scan;
	_curve *nc;
	int l1;

	l1 = 1;
	for (scan = first; scan != ref; scan = scan->next) l1++;
	nc = new _curve();
	nc->first = ref->next;
	nc->first->prev = 0;
	nc->last = last;
	nc->length = length - l1;
	nc->partneratend = partneratend;
	if (partneratend) partneratend->partneratend = nc;

	length = l1;
	last = ref;
	ref->next = 0;
	partneratend = 0;
	return nc;
}


void _curve::append(double x1, double x2) {
	_point *pp;
	pp = new _point(x1, x2, 0);
	if (length == 0) {
		first = pp;
		last = pp;
		pp->prev = 0;
	}
	else {
		last->next = pp;
		pp->prev = last;
		last = pp;
	}
	pp->next = 0;
	length++;
}

void _curve::append(_point *pp) {

	pp->next = last->next;
	pp->prev = last;
	last->next = pp;
	last = pp;
	length++;
}

void _curve::prepend(double x1, double x2) {
	_point *pp;
	pp = new _point(x1, x2, 0);
	if (length == 0) {
		first = pp;
		last = pp;
		pp->next = 0;
	}
	else {
		first->prev = pp;
		pp->next = first;
		first = pp;
	}
	pp->prev = 0;
	length++;
}

_curve *_curve::join(_curve *nc) {
	if (length>0) {
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length>0) {
		nc->first->prev = last;
		last = nc->last;
	}
	length += nc->length;
	partneratend = nc->partneratend;
	if (partneratend) partneratend->partneratend = this;
	nc->first = 0;
	nc->last = 0;
	nc->length = 0;
	delete nc;
	return this;
}

_curve *_curve::joinbefore(_curve *nc) {
	if (length>0) {
		first->prev = nc->last;
	}
	else {
		last = nc->last;
	};
	if (nc->length>0) {
		nc->last->next = first;
		first = nc->first;
	}
	length += nc->length;
	nc->first = 0;
	nc->last = 0;
	nc->length = 0;
	delete nc;
	return this;
}

_curve *_curve::reverse(void) {
	_point *scan1, *scan2, *scambio;
	if (length>1) {
		scan1 = first;
		while (scan1) {
			scan2 = scan1->next;
			scambio = scan1->next;
			scan1->next = scan1->prev;
			scan1->prev = scambio;
			scan1 = scan2;
		}
		scambio = first;
		first = last;
		last = scambio;
	}
	return this;
}

void _curve::drop(_point *ref) {
	_point *scan;
	if (length) {
		for (scan = last; scan && (scan != ref); scan = scan->prev);
		if (scan) {
			if (length == 1) {
				first = last = 0;
			}
			else {
				if (ref->prev) {
					ref->prev->next = ref->next;
					if (ref == last) {
						last = ref->prev;
					}
				}
				if (ref->next) {
					ref->next->prev = ref->prev;
					if (ref == first) {
						first = ref->next;
					}
				}
			}
			length--;
		}
	}
}

double _curve::closest2(_point *ref, _point **clos2) {
	double mi = 1.e100, mi2 = 1.e100, FP;
	_point *scan, *clos;
	if (length>1) {
		clos = *clos2 = first;
		for (scan = first; scan != 0; scan = scan->next) {
			FP = *scan - *ref;
			if (FP<mi) {
				mi2 = mi;
				mi = FP;
				*clos2 = clos;
				clos = scan;
			}
			else if (FP<mi2) {
				mi2 = FP;
				*clos2 = scan;
			}
		}
	}
	else {
		*clos2 = 0;
	}
	return (**clos2 - *ref);
}

double _curve::closest(_point *ref, _point **clos) {
	double mi = 1.e100, FP;
	_point *scan;
	for (scan = first; scan != 0; scan = scan->next) {
		FP = *scan - *ref;
		if (FP<mi) {
			mi = FP;
			*clos = scan;
		}
	}
	return mi;
}

void _curve::complement(_point **sott, int lensott, _point **res, int lenres) {
	int flag, i;
	_point *scan;
	i = 0;
	for (scan = first; scan != 0; scan = scan->next) {
		flag = 0;
		for (int j = 0; (j<lensott) && (!flag); j++) {
			if (scan == sott[j]) {
				flag = 1;
			}
		}
		if ((!flag) && (i<lenres)) {
			res[i] = scan;
			i++;
		}
	}
}

//////////////////////////////
//////////////////////////////
////////_sols methods
//////////////////////////////
//////////////////////////////


_sols::_sols(void) {
	length = 0;
	first = last = 0;
}

_sols::~_sols(void) {
	_curve *scan1, *scan2;
	scan1 = first;
	while (scan1) {
		scan2 = scan1->next;
		delete scan1;
		scan1 = scan2;
	}
}

void _sols::append(_curve *cc) {
	if (length == 0) {
		first = cc;
		last = cc;
		cc->prev = 0;
	}
	else {
		last->next = cc;
		cc->prev = last;
		last = cc;
	}
	cc->next = 0;
	length++;
}

void _sols::prepend(_curve *cc) {
	if (length == 0) {
		first = cc;
		last = cc;
		cc->next = 0;
	}
	else {
		first->prev = cc;
		cc->next = first;
		first = cc;
	}
	cc->prev = 0;
	length++;
}

void _sols::drop(_curve *ref) {
	_curve *scan;
	if (length) {
		for (scan = last; scan && (scan != ref); scan = scan->prev);
		if (scan) {
			if (length == 1) {
				first = last = 0;
			}
			else {
				if (ref->prev) {
					ref->prev->next = ref->next;
					if (ref == last) {
						last = ref->prev;
					}
				}
				if (ref->next) {
					ref->next->prev = ref->prev;
					if (ref == first) {
						first = ref->next;
					}
				}
			}
			length--;
		}
	}
}

void _sols::join(_sols *nc) {
	if (length>0) {
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length>0) {
		nc->first->prev = last;
		last = nc->last;
	}
	length += nc->length;
	nc->first = 0;
	nc->last = 0;
	nc->length = 0;
	delete nc;
}

//////////////////////////////
//////////////////////////////
////////_theta methods
//////////////////////////////
//////////////////////////////

_theta::_theta(double th1) {
	th = th1;
}
_thetas::_thetas(void) {
	length = 0;
}

_thetas::~_thetas(void) {
	_theta *scan, *scan2;
	scan = first;
	while (scan) {
		scan2 = scan->next;
		delete scan;
		scan = scan2;
	}
}

_theta *_thetas::insert(double th) {
	_theta *scan, *scan2;

	scan2 = new _theta(th);
	if (length) {
		if (th<first->th) {
			first->prev = scan2;
			scan2->next = first;
			scan2->prev = 0;
			first = scan2;
		}
		else {
			if (th>last->th) {
				last->next = scan2;
				scan2->prev = last;
				scan2->next = 0;
				last = scan2;
			}
			else {
				scan = first;
				while (scan->th<th) scan = scan->next;
				scan2->next = scan;
				scan2->prev = scan->prev;
				scan->prev->next = scan2;
				scan->prev = scan2;
			}
		}
	}
	else {
		first = scan2;
		last = scan2;
		scan2->next = 0;
		scan2->prev = 0;
	}
	length++;
	//	scan2->maxerr=0.;
	return scan2;
}

void _thetas::remove(_theta* stheta) {
	_theta* scan;
	scan = first;
	while (scan != 0) {
		if (scan == stheta) {
			if (scan != first) scan->prev->next = stheta->next;
			if (scan != last) scan->next->prev = stheta->prev;
			delete stheta;
			length--;
			break;
		}
		scan = scan->next;
	}
}

//////////////////////////////
//////////////////////////////
////////complex methods and operators
//////////////////////////////
//////////////////////////////


complex::complex(double a, double b) {
	re = a;
	im = b;
}

complex::complex(double a) {
	re = a;
	im = 0;
}

complex::complex(void) {
	re = 0;
	im = 0;
}

double abs2(complex z) {
	return (z.re*z.re + z.im*z.im);
}

double abs(complex z) {
	return sqrt(z.re*z.re + z.im*z.im);
}

complex conj(complex z) {
	return complex(z.re, -z.im);
}

complex sqrt(complex z) {
	double md = sqrt(z.re*z.re + z.im*z.im);
	return (md>0) ? complex((sqrt((md + z.re) / 2)*((z.im>0) ? 1 : -1)), sqrt((md - z.re) / 2)) : 0.0;
}

double real(complex z) {
	return z.re;
}

double imag(complex z) {
	return z.im;
}

complex operator+(complex p1, complex p2) {
	return complex(p1.re + p2.re, p1.im + p2.im);
}

complex operator-(complex p1, complex p2) {
	return complex(p1.re - p2.re, p1.im - p2.im);
}

complex operator*(complex p1, complex p2) {
	return complex(p1.re*p2.re - p1.im*p2.im, p1.re*p2.im + p1.im*p2.re);
}

complex operator/(complex p1, complex p2) {
	double md = p2.re*p2.re + p2.im*p2.im;
	return complex((p1.re*p2.re + p1.im*p2.im) / md, (p1.im*p2.re - p1.re*p2.im) / md);
}

complex operator+(complex z, double a) {
	return complex(z.re + a, z.im);
}

complex operator-(complex z, double a) {
	return complex(z.re - a, z.im);
}

complex operator*(complex z, double a) {
	return complex(z.re*a, z.im*a);
}

complex operator/(complex z, double a) {
	return complex(z.re / a, z.im / a);
}

complex operator+(double a, complex z) {
	return complex(z.re + a, z.im);
}

complex operator-(double a, complex z) {
	return complex(a - z.re, -z.im);
}

complex operator*(double a, complex z) {
	return complex(a*z.re, a*z.im);
}

complex operator/(double a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a * z.im / md);
}


complex operator+(complex z, int a) {
	return complex(z.re + a, z.im);
}

complex operator-(complex z, int a) {
	return complex(z.re - a, z.im);
}

complex operator*(complex z, int a) {
	return complex(z.re*a, z.im*a);
}

complex operator/(complex z, int a) {
	return complex(z.re / a, z.im / a);
}

complex operator+(int a, complex z) {
	return complex(z.re + a, z.im);
}

complex operator-(int a, complex z) {
	return complex(a - z.re, -z.im);
}

complex operator*(int a, complex z) {
	return complex(a*z.re, a*z.im);
}

complex operator/(int a, complex z) {
	double md = z.re*z.re + z.im*z.im;
	return complex(a*z.re / md, -a * z.im / md);
}

complex operator-(complex z) {
	return complex(-z.re, -z.im);
}

bool operator==(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

bool operator!=(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

complex expcmplx(complex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return complex(r*cos(theta), r*sin(theta));
}

complex cbrt(complex z) {
	complex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 0.333333333333);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	return 	complex(r_cube*cos(theta_cube), r_cube*sin(theta_cube));
}
