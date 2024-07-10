
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#define _Jac\
	S2 = 0;\
	for (int ik = 0; ik < n; ik++) {\
		pza[ik] = z - a[ik];\
		pmza[ik][ik]=m[ik]/pza[ik];\
		pmza2[ik][ik] = pmza[ik][ik] / pza[ik];\
		S2 = S2 + pmza2[ik][ik];\
	}\
	Jac=1-abs2(S2);

#define _L0\
	S1 = 0;\
	for (int ik = 0; ik < n; ik++) {\
		S1=S1 +pmza[ik][ik];\
	}\
	Lv=yc - conj(z) + S1;\
	Lnew = abs2(Lv);

#define _S3\
	S3=0;\
	for (int ik = 0; ik < n; ik++) {\
		pmza2[ik][ik] = pmza2[ik][ik] / pza[ik];\
		S3= S3 + pmza2[ik][ik];\
	}

#define _S4\
	S4=0;\
	for (int ik = 0; ik < n; ik++) {\
		S4= S4 + pmza2[ik][ik]/pza[ik];\
	}
//////////////////////////////
//////////////////////////////
////////Geometry functions
//////////////////////////////
//////////////////////////////

// SetLensGeometry
// initializes the arrays describing the geometric configuration of the lens system.
////// Input parameters
// nn: number of lenses; 
// q: array of size n-1 with the mass ratios assuming the first lens to have mass 1;
// s: array of size n-1 with the separations in complex coordinates from the first mass.
////////// What does it make?
// VBML::n is set to the number of lenses.
// VBML::m is an array of size n with the masses (normalized to 1 total mass)
// Note that the first mass is now the lowest mass in the system to improve numerical precision.
// VBML::a is an array of size n with the separations from the first mass;
// VBML::pza is Product(z-a[i]) (a polynomial of degree n)
// VBML::pmza is an array of n polynomials m[i] * Product(z-a[j]) (polynomials of degree n-1)


void VBMicrolensing::SetLensGeometry(int nn, double* pr) {
	double* q = (double*)malloc(sizeof(double) * (nn));
	complex* s = (complex*)malloc(sizeof(complex) * (nn));

	for (int i = 0, j = 2; i < nn; ++i, j += 3) {
		q[i] = pr[j];
	}
	for (int i = 0, j = 0; i < nn; ++i, j += 3) {
		s[i] = complex(pr[j], pr[j + 1]);
	}

	SetLensGeometry(nn, q, s);

	free(q);
	free(s);
}

void VBMicrolensing::SetLensGeometry(int nn, double* q, complex *s) {
	switch (SelectedMethod)
	{
	case Method::Singlepoly:
		SetLensGeometry_spnp(nn,q,s);
		break;
	case Method::Multipoly:
		SetLensGeometry_multipoly(nn, q, s);
		break;
	case Method::Nopoly:
		SetLensGeometry_spnp(nn, q, s);
		break;
	}
}

void VBMicrolensing::SetMethod(Method Met) {
	SelectedMethod = Met;
}

void VBMicrolensing::SetLensGeometry_spnp(int nn, double* q, complex* s) {
	static double sumq, qmin, Jac;
	static int iqmin, dg;
	static complex pbin[2], z, S2, fac;
	static int i, j;

	change_n(nn);

	qmin = sumq = q[0];
	iqmin = 0;
	for (i = 1; i < n; i++) {
		if (q[i] < qmin) {
			qmin = q[i];
			iqmin = i;
		}
		sumq += q[i];
	}
	//	sumq = 1; // Do not want normalized masses.

	m[0] = q[iqmin] / sumq;
	a[0] = 0;
	*s_offset = s[iqmin];
	for (i = 0; i < iqmin; i++) {
		m[i + 1] = q[i] / sumq;
		a[i + 1] = s[i] - *s_offset;
	}
	for (i = iqmin + 1; i < n; i++) {
		m[i] = q[i] / sumq;
		a[i] = s[i] - *s_offset;
	}

	// Central images for close-by lenses. Only used if J<0 as additional initial conditions.
	// Central images are calculated here because they do not depend on the source.
	lencentralimages = 0;
	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			if (abs(a[i] - a[j]) < m[i] + m[j])
				centralimages[lencentralimages] = (m[j] * a[i] + m[i] * a[j]) / (m[i] + m[j]);
			z = centralimages[lencentralimages];
			_Jac
				if (Jac < 0) lencentralimages++;
		}
	}

	for (int i = 0; i < n; i++) {
		pmza[i][0] = m[i];
	}
	pza[0] = pbin[1] = 1.0;
	for (int i = 0; i < n; i++) {
		pbin[0] = -a[i];
		for (int j = 0; j < n; j++) {
			if (i == j) {
				copypol(pza, i, pdum);
				polyproduct(pdum, i, pbin, 1, pza);
			}
			else {
				dg = (j > i) ? i : i - 1;
				copypol(pmza[j], dg, pdum);
				polyproduct(pdum, dg, pbin, 1, pmza[j]);
			}
		}
	}
}

void VBMicrolensing::SetLensGeometry_multipoly(int nn, double* q, complex* s) {
	static int j, i, x, k, p, dg;
	static double tempq, sumq;
	static complex temps, pbin[2];

	change_n_mp(nn);

	sumq = 0;

	for (i = 0; i < n; i++) {
		q_sort[i] = q[i];
		s_sort[i] = s[i];
		sumq += q[i];
	}

	x = nn;

	for (j = 0; j < x - 1; j++) {
		do {
			k = 0;
			for (i = 0; i < x - 1; i++) {
				if (q_sort[i] > q_sort[i + 1]) {
					tempq = q_sort[i];
					q_sort[i] = q_sort[i + 1];
					q_sort[i + 1] = tempq;

					temps = s_sort[i];
					s_sort[i] = s_sort[i + 1];
					s_sort[i + 1] = temps;

					k = 1;
					p = i + 1;
				}
			}
			x = p;
		} while (k == 1);
	}

	//normalize
	for (i = 0; i < n; i++) {
		q_sort[i] = q_sort[i] / sumq;
	}

	//set all reference systems, with the nth smaller mass in the first place
	for (i = 0; i < n; i++) {
		m_mp[0][i] = q_sort[i];
		a_mp[0][i] = s_sort[i] - s_sort[0];
	}
	for (j = 1; j < n; j++) {
		m_mp[j][0] = q_sort[j];
		a_mp[j][0] = complex(0, 0);
		for (i = 1; i < n; i++) {
			if (i == j) {
				m_mp[j][i] = q_sort[0];
				a_mp[j][i] = s_sort[0] - s_sort[j];
			}
			else {
				m_mp[j][i] = q_sort[i];
				a_mp[j][i] = s_sort[i] - s_sort[j];
			}
		}
	}

	*s_offset = s_sort[0];
	for (int i = 0; i < n; i++) {
		m[i] = m_mp[0][i];
		a[i] = a_mp[0][i];
	}

	// first steps for the creation of the polynomials

	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			pmza_mp[k][i][0] = m_mp[k][i];
		}
		pza_mp[k][0] = pbin[1] = 1.0;
		for (int i = 0; i < n; i++) {
			pbin[0] = -a_mp[k][i];
			for (int j = 0; j < n; j++) {
				if (i == j) {
					copypol(pza_mp[k], i, pdum);
					polyproduct(pdum, i, pbin, 1, pza_mp[k]);
				}
				else {
					dg = (j > i) ? i : i - 1;
					copypol(pmza_mp[k][j], dg, pdum);
					polyproduct(pdum, dg, pbin, 1, pmza_mp[k][j]);
				}
			}
		}
	}
}

//////////////////////////////
//////////////////////////////
////////Basic magnification functions
//////////////////////////////
//////////////////////////////

#define EXECUTE_METHOD(METHOD, THETA)                     \
    switch (METHOD) {                                      \
        case Method::Singlepoly:                           \
            polycoefficients();                            \
            Prov = NewImagespoly(THETA);                   \
            break;                                         \
        case Method::Multipoly:                            \
            for (int i = 0; i < n; i++) {                  \
                y_mp[i] = y + s_sort[0] - s_sort[i];       \
            }                                              \
            polycoefficients_multipoly();                  \
            Prov = NewImagesmultipoly(THETA);              \
            break;                                         \
        case Method::Nopoly:                               \
            Prov = NewImages(THETA);                       \
            break;                                         \
    }

double VBMicrolensing::MultiMag0(complex yi, _sols **Images) {
	static double Mag = -1.0;
	_theta *stheta;
	_curve *Prov, *Prov2;
	_point *scan1, *scan2;

	stheta = new _theta(-1.);
	
	y = yi - *s_offset; // Source position relative to first (lowest) mass
	rho = rho2 = 0;

	(*Images) = new _sols;
	corrquad = corrquad2 = 0; // to be implemented for v2.0
	safedist = 10;

	EXECUTE_METHOD(SelectedMethod, stheta)
	
	Mag = 0.;
	nim0 = 0;
	for (scan1 = Prov->first; scan1; scan1 = scan2) {
		scan2 = scan1->next;
		Prov2 = new _curve(scan1);
		(*Images)->append(Prov2);
		Mag += fabs(1 / scan1->dJ);
		nim0++;
	}
	Prov->length = 0;
	delete Prov;
	delete stheta;
	return Mag;

}

double VBMicrolensing::MultiMag0(complex y) {
	static _sols* images;
	static double mag;
	mag = MultiMag0(y, &images);
	delete images;
	return mag;
}

double VBMicrolensing::MultiMag0(double y1, double y2) {
	static _sols* images;
	static double mag;
	complex y = complex(y1, y2);
	mag = MultiMag0(y,&images);
	delete images;
	return mag;
}

double VBMicrolensing::MultiMag(complex yi, double RSv, double Tol, _sols **Images) {
	static complex y0;
	static double Mag = -1.0, th, thoff = 0.01020304,thoff2= 0.7956012033974483; //0.01020304
	static double errimage, maxerr, currerr, Magold,rhorad2,th2;
	static int NPSmax, flag, NPSold,isquare, flagfinal;
	static _thetas *Thetas;
	static _theta *stheta,*itheta,*jtheta;
	static _curve *Prov, *Prov2;
	static _point *scan1, *scan2;
	static int lsquares[4];

	

	try {
		y0 = yi - *s_offset; // Source position relative to first (lowest) mass
		rho = RSv;
		rho2 = RSv * RSv;
			

		NPS = 1;
		// Two channels: Tol>1 uses Tol points in the contour; Tol<1 uses Tol as accuracy goal and RelTol as precision goal
		if (Tol > 1.) {
			errimage = 0.;
			NPSmax = (int)(Tol);
		}
		else {
			errimage = Tol * M_PI*RSv*RSv;
			NPSmax = 2000;
		}

		// Calculation of the images

		(*Images) = new _sols;
		Thetas = new _thetas;
		th = thoff;
		stheta = Thetas->insert(th);
		stheta->maxerr = 1.e100;
		y = y0 + complex(RSv*cos(thoff), RSv*sin(thoff)); // first image
		

#ifdef _PRINT_TIMES
		tim0 = Environment::TickCount;
#endif
		
		EXECUTE_METHOD(SelectedMethod, stheta)

#ifdef _PRINT_TIMES
		tim1 = Environment::TickCount;
		GM += tim1 - tim0;
#endif
		stheta = Thetas->insert(2.0*M_PI + thoff);
		stheta->maxerr = 0.;
		stheta->Mag = 0.;
		stheta->errworst = Thetas->first->errworst;
		for (scan1 = Prov->first; scan1; scan1 = scan2) {
			scan2 = scan1->next;
			Prov2 = new _curve(scan1);
			Prov2->append(scan1->x1, scan1->x2);
			Prov2->last->theta = stheta;
			Prov2->last->d = Prov2->first->d;
			Prov2->last->dJ = Prov2->first->dJ;
			Prov2->last->ds = Prov2->first->ds;
			(*Images)->append(Prov2);
		}
		Prov->length = 0;
		delete Prov;

		th = thoff;
		for (int i = 0; i < 3; i++) {
			th += M_PI_2;
			stheta = Thetas->insert(th);
			y = y0 + complex(RSv*cos(th), RSv*sin(th));

			EXECUTE_METHOD(SelectedMethod, stheta)

			OrderMultipleImages((*Images), Prov);
		}
		NPS = 4;

		// Let us check the number if images in a square around the source
		if (squarecheck) {
			itheta = Thetas->first;
			stheta = new _theta(-1);
			while (itheta->next) {
				th = (itheta->th + itheta->next->th)*0.5;
				rhorad2 = RSv / cos((itheta->next->th - itheta->th)*0.5);
				y = y0 + complex(rhorad2*cos(th), rhorad2*sin(th));

				EXECUTE_METHOD(SelectedMethod, stheta)
				
				delete Prov;
				if (itheta->imlength != stheta->imlength && itheta->next->imlength != stheta->imlength) {
					jtheta = Thetas->insert(th);
					y = y0 + complex(RSv*cos(th), RSv*sin(th));

					EXECUTE_METHOD(SelectedMethod, stheta)

					OrderMultipleImages((*Images), Prov);
					NPS++;
				}
				else {
					itheta = itheta->next;
				}
			}
			delete stheta;
		}

		// Initial error calculation

		maxerr = currerr = Mag = 0.;
		stheta = Thetas->first;
		while (stheta->next) {
			Mag += stheta->Mag;
			if (stheta->next->th - stheta->th > 1.e-8) {
				currerr += stheta->maxerr;
#ifndef _uniform
				if (stheta->maxerr > maxerr) {
					maxerr = stheta->maxerr;
#else
				if (stheta->next->th*0.99999 - stheta->th > maxerr) {
					maxerr = stheta->next->th - stheta->th;
#endif
					itheta = stheta;
				}
			}
			stheta = stheta->next;
		}
		th = (itheta->th + itheta->next->th) *0.5;

		// Main cycle: sampling continues until total error is below Tol.
		flag = 0;
		Magold = -1.;
		NPSold = NPS + 1;

		while (((currerr > errimage) && (currerr > RelTol*Mag) && (NPS < NPSmax) && (flag < NPSold))) {
			stheta = Thetas->insert(th);
			y = y0 + complex(RSv*cos(th), RSv*sin(th));
#ifdef _PRINT_TIMES
			tim0 = Environment::TickCount;
#endif
			EXECUTE_METHOD(SelectedMethod, stheta)

#ifdef _PRINT_TIMES
			tim1 = Environment::TickCount;
			GM += tim1 - tim0;
#endif
#ifdef _PRINT_ERRORS2
			int lim = Prov->length;
#endif
			// Assign new images to correct curves
			OrderMultipleImages((*Images), Prov);

			maxerr = currerr = Mag = 0.;
			stheta = Thetas->first;
			//	printf("\n");
			while (stheta->next) {
				//	printf("\n%lf %le", stheta->th, stheta->maxerr);
				Mag += stheta->Mag;
				if (stheta->next->th - stheta->th > 1.e-8) {
					currerr += stheta->maxerr;
#ifndef _uniform
					if (stheta->maxerr > maxerr) {
						maxerr = stheta->maxerr;
#else
					if (stheta->next->th*0.99999 - stheta->th > maxerr) {
						maxerr = stheta->next->th - stheta->th;
#endif
						itheta = stheta;
					}
				}
				stheta = stheta->next;
			}
			th = (itheta->th + itheta->next->th) *0.5;
			NPS++;

#ifndef _uniform
			if (fabs(Magold - Mag) * 2 < errimage) {
				flag++;
			}
			else {
				flag = 0;
				Magold = Mag;
				NPSold = NPS + 1;
			}
#else
			currerr = 2 * errimage;
			if (NPS == 2 * NPSold) {
				if (fabs(Magold - Mag) * 2 < errimage) {
					flag = NPSold;
				}
				else {
					flag = 0;
					NPSold = NPS;
					Magold = Mag;
				}
			}
#endif
#ifdef _PRINT_ERRORS2
			printf("\nNPS= %d nim=%d Mag = %lf maxerr= %lg currerr =%lg th = %lf", NPS, lim, Mag / (M_PI*RSv*RSv), maxerr / (M_PI*RSv*RSv), currerr / (M_PI*RSv*RSv), th);
#endif

		}
		Mag /= (M_PI*RSv*RSv);
		therr = currerr / (M_PI*RSv*RSv);

		delete Thetas;

		return Mag;

	}
	catch (...) {
		FILE *f = fopen("Geom.txt","w");
		fprintf(f,"\n%d\n", n);
		for (int i = 0; i < n; i++) {
		fprintf(f,"%.16lf %.16lf %.16lf\n", m[i], a[i].re, a[i].im);
		}
		fprintf(f,"\n%.16lf %.16lf %.16lf\n", y.re, y.im, rho);
		for (int i = 0; i < ngood; i++) {
		fprintf(f, "%.16lf %.16lf %.16lg %.16lf %.16lg\n", zr[i].re, zr[i].im, errs[i],Jacs[i],good[i]);
		}
		fclose(f);
		return -1;
	} 

	
	
}

double VBMicrolensing::MultiMag(complex y, double RSv, double Tol) {
	static _sols *images;
	static double mag;
	mag = MultiMag(y, RSv, Tol, &images);
	delete images;
	return mag;
}

double VBMicrolensing::MultiMag(complex y, double RSv) {
	static _sols* images;
	static double mag;
	mag = MultiMag(y, RSv, Tol, &images);
	delete images;
	return mag;
}

double VBMicrolensing::MultiMag(double y1, double y2, double RSv) {
	static _sols* images;
	static double mag;
	mag = MultiMag(complex(y1, y2), RSv, Tol, &images);
	delete images;
	return mag;
}



///////////////////////////////////////////////
///////////////////////////////////////////////
////////// Internal private functions
///////////////////////////////////////////////
///////////////////////////////////////////////



#define _MJacobians1 \
	za[0][i]=z;\
	za2[0][i]=z*z;\
	zaltc[i] = yc+m[0]/z; \
	J1[i]=m[0]/za2[0][i];\
	for(int j=1;j<n;j++){\
		za[j][i]=z-a[j];\
		za2[j][i]=za[j][i]*za[j][i];\
		zaltc[i]=zaltc[i]+ m[j] / za[j][i];\
		J1[i]=J1[i]+ m[j] / za2[j][i];\
	}\
	J1c[i]=conj(J1[i]);\
	Jacs[i]=(1-J1[i]*J1c[i]).re;\
	LL = zaltc[i] - zc;

#define _MJacobians2 \
	J2=0;\
	for(int j=0;j<n;j++){\
		za2[j][worst[i]]=za2[j][worst[i]]*za[j][worst[i]];\
		J2=J2+m[j] / za2[j][worst[i]];\
	}\
	J2=-J2*2;\


#define _MJacobiansgood\
	dy = complex(-sin(theta->th)*rho, cos(theta->th)*rho);\
	dz = (dy - J1c[worst[i]]*conj(dy)) / Jacs[worst[i]];\
	Prov->last->d = dz;\
	Prov->last->J2 = J2;\
	Prov->last->ds = (imag(dy*dz*dz*J2) + rho2) / Jacs[worst[i]];\

#define _MJacobians0goodnopoly\
	J3=6*S4s[i];\
	dJ2=Jacs[i]*Jacs[i];\
	J1c2=J1c[i]*J1c[i];\
	J3=J3*J1c2;\
	ob2=(J2.re*J2.re+J2.im*J2.im)*(6-6*Jacs[i]+dJ2);\
	J2=J2*J2*J1c2*J1c[i];\
	cq= 0.5*(fabs(ob2-6.*J2.re-2.*J3.re*Jacs[i])+3*fabs(J2.im))/fabs(Jacs[i]*dJ2*dJ2); 


#define _MJacobians0good\
	J3=0;\
	for(int j=0;j<n;j++){\
		za2[j][worst[i]]=za2[j][worst[i]]*za[j][worst[i]];\
		J3=J3+m[j] / za2[j][worst[i]];\
	}\
	J3=J3*6;\
	dJ2=Jacs[worst[i]]*Jacs[worst[i]];\
	J1c2=J1c[worst[i]]*J1c[worst[i]];\
	J3=J3*J1c2;\
	ob2=(J2.re*J2.re+J2.im*J2.im)*(6-6*Jacs[worst[i]]+dJ2);\
	J2=J2*J2*J1c2*J1c[worst[i]];\
	cq= 0.5*(fabs(ob2-6.*J2.re-2.*J3.re*Jacs[worst[i]])+3*fabs(J2.im))/fabs(Jacs[worst[i]]*dJ2*dJ2); 

#define _MJacobians0bad\
	Jaltc=0;\
	for (int j = 0; j < n; j++) {\
		Jaltc=Jaltc+ m[j] / (zaltc[worst[i]]-conj(a[j]));\
	}\
	Jalt = conj(Jaltc);\
	JJalt2=(1-J1c[worst[i]]*Jalt);\
	J3=J2*J1c[worst[i]]*JJalt2;\
	J3=(J3-conj(J3)*Jalt)/(JJalt2*JJalt2*Jacs[worst[i]]);\
	cq=(J3.re*J3.re+J3.im*J3.im);

void VBMicrolensing::initroot() {
	static complex fac,fac2,z,S2;
	static double Jac;

	// Main image, should be positive and close to source
	init[n] = y;
	for (int i = 0; i < n; i++) {
		fac = y - a[i];
		prodevs[i] = m[i] / abs2(fac);
		devs[i] = (prodevs[i] > 1.e-2) ? (sqrt(0.25 + prodevs[i]) - 0.5) : prodevs[i];
		devs[i] = devs[i] * fac;
		init[n] = init[n]+ devs[i];
	}
	z = init[n];
	_Jac
	while (Jac <= 0) {	// Check if main image has positive parity. Otherwise displace
		z = z + ((abs2(y) > 1.e-6) ? y : complex(0, 1));
		_Jac
	}
	init[n] = z;

	// Secondary images, should be negative and close to lenses
	for (int i = 0; i < n; i++) {
		init[i] = a[i] - devs[i];
		z = init[i];
		_Jac
		while (Jac >= 0) {	// Check if secondary image has negative parity. Otherwise displace
			devs[i] = devs[i]*0.5;
			z = z + devs[i];
			_Jac
		}
		init[i] = z;
	}
	// Central images are calculated in SetLensGeometry
}

int VBMicrolensing::froot(complex zi) {
	static complex z, zc, S1,S2, S2v,S3,zo,zo2,epso,epsbase,epsn,epsl,gradL,zl,fac,fac2,dz,dzo,TJold,TJnew,Lv,den;
	static int iter3,iter4, ipseudo,flagmain,flag;
	static double Lnew, Lold, fad,Jac,Jacold,prefac, epsbo;

	Lold = 100.;
	epso = 1.e100;
	zo2 = 0;
	Jacold = 1.;
	iter = iter2 = 1;
	z = zi;
	_Jac
	_L0
	_S3
	fad = 1. - Jac;
	epsbase = epsn = ((Jac>-15)? fad*(sqrt(sqrt(fad)) - 1) : fad) / (conj(S2)*S3)*0.5;
	zo = z + epsbase;
	flagmain = 0; // Reports the corrections made during the step
	ipseudo = 0;
	// Main cycle. Here are the stopping conditions:
	// epsbase is the basic step calculated by Newton Method. We stop if it smaller than 0.3e-10
	// epsn is the actual step after corrections (see below). We stop if it smaller than 1.e-14
	// Lnew is the lens equation squared. We stop if it is smaller than 1.e-29
	// iter2 is the number of consecutive oscillations
	// iter is the number of iterations
	while ((abs2(epsbase) > rootaccuracy) && (abs2(epsn) > 1.e-28) && (Lnew > 1.e-29) && (iter2 < 9) && (iter < maxiter)) {
		zo2 = zo; // Stores z of two steps before.
		zo = z; // Stores z of previous step.
		Lold = Lnew; // Stores previous value of lens eq.
		prefac = fabs(Jac / Jacold); // prefac warns if we are leaving a critical curve
		Jacold = Jac; // Stores old Jac
		epso = epsn; // Stores old value of epsn
	
		epsbase = epsn = (conj(Lv) - Lv * conj(S2))/Jac;
		flag = 1; // Checks that the new point is acceptable
		iter3 = 0; // Counts the number of corrections in step
		while (flag > 0 && iter3 < 9) {
			if (flagmain == 3) {
				epsn = epsn * pseudorandom[ipseudo]; // In case of oscillation, shorten step by pseudorandom
				ipseudo++;
				if (ipseudo > 11) ipseudo = 0;
			}
			flagmain = 0;
			fad = abs2(epsn);
			while (abs2(zo2 - zo - epsn) < 1.69*fad && fad > 1.e-28) {
				fad = fad;
				epsn = epsn * pseudorandom[ipseudo]; // In case of doubleoscillation, shorten step by pseudorandom
				fad = abs2(epsn);
				ipseudo++;
				if (ipseudo > 11) ipseudo = 0;
				flagmain = 4;
			}

			// New step cannot be twice longer than previous step (avoid explosion)
			fad = sqrt(abs2(zo2 - zo) /fad);
			fad*=(prefac < 1) ? (1 + prefac) : 2;
			if (fad < 1) {
				epsn = epsn * fad;
				flagmain = 5;
			}
			if (abs2(epsn) < 1.e-28) {
				break;
			}

			// Update z with corrected step
			z = zo + epsn;
			_Jac
				// Check that Jacobian does not change sign
			if (signbit(Jac) != signbit(Jacold)) {
				// Correct by stepping through closest critical curve. Factor 1.5 keeps away enough from critical curve.
				// But maybe better try to keep the same Jacobian as before(Jac -	Jacold) / Jac
				_S3
				fad = 1. - Jac;
				fad *= (1 - Jacold / Jac)*(sqrt(sqrt(fad)) - 1);
				dz = dzo =  fad/ (conj(S2)*S3);
				// If correction is below machine precision, go back a little bit 
				if (abs2(dz) < 1.e-28) {
					dz = -epsn * 0.15;
					flagmain = 3;
				}
				// If correction is higher than step shorten correction and increase iter2. We want to break cycles that do not converge
				fad = abs(epsn / dz);
				if (abs2(epsn + dz) < 0.25*abs2(epsn) || fad < 2) {
					dz = complex(0.33, 0.01)*dz*fad;
					flagmain = 3;
					if(abs2(epsn + dz)<abs2(zo2-zo)) iter2++;
				}
				// If correction is not enough to step through critical curve, increase it 
				iter4 = 0;
				while (signbit(Jac) != signbit(Jacold) && iter4 < 6) {
					z = z + dz;
					_Jac
						iter4++;
				}
				if (iter4 < 6) {
					// Accept new point and proceed to next step
					_L0
					if (Lnew > Lold) {
						flagmain = 2;
					}
					flag = 0;
				}
				else {
					// Do not accept. Repeat cycle with shorter step 
					iter3++;
					epsn = epsn * 0.5;
					flag = 1;
				}
			}
			else {
				_L0
				if (Lnew > Lold) {
					dz = z;
					z = zo;
					_Jac
					_S3
					TJold = S2 * conj(S3);
					z = dz;
					_Jac
					_S3
					TJnew = S2 * conj(S3);
					if (TJold.re*TJnew.re + TJold.im*TJnew.im < 0) {
						// If new point is worse and in different Jacobian gradient region,
						// shorten step (avoid overshooting) 
						iter4 = 0;
						while (Lnew > Lold && iter4 < 10) {
							epsn = epsn * 0.5;
							z = zo + epsn;
							_Jac
							_L0
							iter4++;
						}
					}
					else {
						// If new point is worse but in same Jacobian gradient region, correct toward critical curve (when Jac > 0)
						if (Jac > 0) {
							gradL = -(conj(Lv) + conj(S2v)*Lv);
							den = (TJnew.re*gradL.re + TJnew.im*gradL.im);
							epsl = -0.25*Lnew / den * TJnew;
							zl = z;
							z = zl + epsl;
							_Jac
							while (signbit(Jac) != signbit(Jacold)) {
								epsl = epsl * 0.5;
								z = zl + epsl;
								_Jac
							}
							_L0
						}
					}
				}
				// If sign of the Jacobian is correct and no problems on L, accept step 
				flag = 0;
			}
		}
		//Prevent oscillations by checking that we are not going back and forth
		if (abs2(z - zo2)  < 0.25*abs2(z - zo)) {
			flagmain = 3;
			iter2++;
		}

		if(flagmain == 0) iter2 = 0;

		iter++;
	}
	newtonstep += iter;
	err = 3.163e-15 / Jac;
	err *= err;
	err+=abs2(epsbase);
	zf = z;
	Jacf = Jac;
	S2f = S2;
	L0f = Lnew;
	return iter2;
}

bool VBMicrolensing::checkroot(_theta *theta) {
	static double mn, fac;
	static int imn;
	static complex S3,z,S4;
	static double fad;
	if ((iter2 < 9 && iter < maxiter) || L0f < 1.e-29) {
		mn = 1.e100;
		imn = 0;
		for (int i = 0; i < ngood;i++) {
			if (signbit(Jacf) == signbit(Jacs[i])) {
				fac = abs2(zr[i] - zf) / (errs[i]+err + 4.e-20);
				if (fac < mn) {
					mn = fac;
					imn = i;
				}
			}	
		}
		if (mn > 10.) { // Good root, accept
			zr[ngood] = zf;
			errs[ngood] = err;
			Jacs[ngood] = Jacf;
			good[ngood] = L0f;
			S2s[ngood] = S2f;
			z = zf;
			_S3
			S3s[ngood] = S3;
			if (theta->th < 0) {
				_S4
				S4s[ngood] = S4;
			}
			fad = (1 - Jacf);
			grads[ngood] = fad*(sqrt(sqrt(fad)) - 1) / (conj(S2f)*S3);

			ngood++;
		}
		else { // Duplicate of old root
			if (err < errs[imn]) { // Replace old root
				zr[imn] = zf;
				errs[imn] = err;
				Jacs[imn] = Jacf;
				good[imn] = L0f;
			}
		}
		return true;
	} // Bad root, reject
	return false;
}


_curve *VBMicrolensing::NewImages(_theta *theta) {
	static _curve *Prov;
	static int nminus,nplus;
	static complex z, zc,dy,dz,J2,J3,Jalt,JJalt2,Jaltc,J1c2;
	static complex S2, S2c, S3, S3c,vec,newseed0;
	static double ob2,dJ2,cq,Jac,imul,phi;
	static complex newseedtrial[6] = {complex(1,0.5),complex(1,-0.5), 0.5, 0.75, 2., 4.,};
	static int imass, iphi,nsafe;

	yc = conj(y);
	initroot();

	ngood = ngoodold=0;
	for (int i = 0; i < n + 1; i++) {
		froot(init[i]);
		checkroot(theta);
	}
	for (int i = 0; i < lencentralimages; i++) {
		froot(centralimages[i]);
		checkroot(theta);
	}

	while (ngood > ngoodold) {
		lennewseeds = 0;
		for (int i = ngoodold; i < ngood; i++) {
			vec = complex(2, 1 + sqrt(errs[i]) / abs(grads[i]));
			newseeds[lennewseeds] = zr[i]+grads[i]*vec;
			z = newseeds[lennewseeds];
			_Jac
			iter = 0;
			newseed0 = newseeds[lennewseeds];
			while (signbit(Jac) == signbit(Jacs[i]) && iter<6) {
				newseeds[lennewseeds] = zr[i] + (newseed0 - zr[i])*newseedtrial[iter];
				z = newseeds[lennewseeds];
				_Jac
				iter++;
			}
			if(iter<6) lennewseeds++;
			newseeds[lennewseeds] = zr[i] + grads[i] * conj(vec);
			z = newseeds[lennewseeds];
			_Jac
			iter = 0;
			newseed0 = newseeds[lennewseeds];
			while (signbit(Jac) == signbit(Jacs[i]) && iter<6) {
				newseeds[lennewseeds] = zr[i] + (newseed0 - zr[i])*newseedtrial[iter];
				z = newseeds[lennewseeds];
				_Jac
				iter++;
			}
			if (iter<6) lennewseeds++;
		}
		ngoodold = ngood;
		for (int i = 0; i < lennewseeds; i++) {
			froot(newseeds[i]);
			checkroot(theta);
		}
	}

	nminus = nplus = 0;
	for (int i = 0; i < ngood; i++) {
		if (Jacs[i] > 0) {
			nplus++;
		}
		else {
			nminus++;
		}
	}

	imul = 1.;
	iphi = imass=nsafe=0;
	
	while (nminus - nplus + 1 != n) {
		
		phi = iphi * 2.61799; // 5*M_PI/6.
		z = imul*sqrt(m[imass])*complex(cos(phi), sin(phi))+a[imass];
		froot(z);
		checkroot(theta);
		if (ngood > ngoodold) {
			if (Jacs[ngoodold] > 0) {
				nplus++;
			}
			else {
				nminus++;
			}
		}
		imass++;
		if (imass == n) {
			imass = 0;
			iphi++;
			if (iphi == 12) {
				iphi = 0;
				if (imul > 1) {
					imul *= 1.1;
				}else{
					imul *= 0.9;
				}
				imul = 1. / imul;
			}
		}
		nsafe++;
	}
	
	Prov = new _curve;
	for (int i = 0; i < ngood; i++) {
		worst[i] = i;
		z = zr[i];
		zc = conj(z);
		Prov->append(z.re + s_offset->re, z.im + s_offset->im);
		Prov->last->dJ = Jacs[i];
		J2 = -2 * S3s[i];
		if (theta->th >= 0) {
			J1c[i] = conj(S2s[i]);
			_MJacobiansgood
		}
		else {
			_MJacobians0goodnopoly
			corrquad += cq;
		}
		Prov->last->theta = theta;
	}
	corrquad2 = 1.e200;
	theta->imlength = Prov->length;

	return Prov;
}

void VBMicrolensing::initrootpoly() {
	static double mrt;
	static complex dev, dev2, zplus, shear, alpha0;
	static int ir;
	zplus = y;
	ir = nroots - 1;
	for (int i = 0; i < n; i++) {
		J1[i] = sqrt(m[i]);
		dev2 = 0;
		for (int j = 0; j <n; j++) {
			if (i != j)
				dev2 = dev2 + m[j] / (a[i] - a[j]);
		}
		dev = (sqrt(0.25 + m[i] / abs2(y - a[i])) - 0.5)*(y - a[i]);
		if (abs2(dev2) < 1) {
			zr[ir] = a[i] - dev;
		}
		else {
			dev2 = dev2 + conj(y - a[i]);
			zr[ir] = a[i] - m[i] / dev2;
		}
		ir--;
		zplus = zplus + dev;
	}
	zr[ir] = zplus;
	ir--;

	for (int i = 0; i < n; i++) {
		dev = zplus - a[i];
		dev2 = dev / abs(dev)*J1[i] / 1.;
		for (int j = 0; j < 2 * (n - 1); j++) {
			dev = M_PI * (j / (n - 1.) + 0.5);
			zr[ir] = a[i] + dev2 * complex(cos(dev.re), sin(dev.re));
			ir--;
		}
	}
}

int VBMicrolensing::findimagepoly(int i) {
	static complex z, zc, yc, LL, zo, delta, lambda;
	static double dlmax = 1.0e-12, LLold, Jold, deltafac;
	static int iter, iter2, success;
	yc = conj(y);
	z = zr[i];
	zc = conj(z);
	_MJacobians1

	LLold = 1.e100;
	good[i] = abs2(LL);
	Jold = Jacs[i];
	iter = iter2 = 0;
	deltafac = 1.;

	while (good[i] > dlmax /*&& good[i]<LLold && signbit(Jold) == signbit(Jac[i])*/ && iter<0 && iter2<4) {
		LLold = good[i];
		Jold = Jacs[i];
		zo = z;
		lambda = 1.0;
		//for (int j = nroots-1; j > i; j--) {
		//	if(good[j]<dlmax)
		//		lambda = lambda + conj(LL) /(z - zr[j]);
		//}
		delta = (conj(LL*lambda) - LL * J1c[i]) / (Jacs[i] - 1 + abs2(lambda));
		iter2 = 0;
		deltafac = 1.;
		while (iter2<4 && (good[i] + dlmax > LLold /*|| signbit(Jold) != signbit(Jac[i])*/)) {
			z = zo + delta * deltafac;
			zc = conj(z);
			_MJacobians1
				good[i] = abs2(LL);
			deltafac *= 0.1;
			iter2++;
		}
		//		deltafac *= 4.;

		iter++;
	}
	success = (good[i] > dlmax) ? 0 : 1;
	zr[i] = z;
	return success;
}
_curve *VBMicrolensing::NewImagespoly(_theta *theta) {
	static complex  yc, z, zc, zo, delta, dy, dz, J2, J3, Jalt, Jaltc, JJalt2, LL, J1c2, dzita;
	static double dlmax = 1.0e-12, dzmax = 1.e-10, dJ2, ob2, cq, Jold, LLold;
	static int ngood, nplus, nminus, bad, isso, ncrit, igood, iter, iter2;
	static double mi, tst, isgood;
	static _curve *Prov;
	static _point *scan, *prin, *fifth, *left, *right, *center;

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	yc = conj(y);

#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif
	cmplx_roots_gen(zr, coefs, n2 + 1, false, false);

#ifdef _PRINT_TIMES
	tim1 = Environment::TickCount;
	inc += tim1 - tim0;
#endif
	
	for (int i = n2; i >= 0; i--) {
		findimagepoly(i);
		for (int j = n2; j > i; j--) {
			if (good[j] < 1.e9 && abs2(zr[j] - zr[i]) < dzmax && signbit(Jacs[j]) == signbit(Jacs[i])) {
				if (good[i] > good[j]) {
					good[i] = 1.e100;
					break;
				}
				else good[j] = 1.e100;
			}
		}
	}

	for (int i = 0; i < n2 + 1; i++) {
		int j = 0;
		while (j < i && good[worst[j]] < good[i]) j++; // worst contains root indices in increasing order of badness
		for (int k = i; k > j; k--) worst[k] = worst[k - 1];
		worst[j] = i;
	}
	for (int i = 0; i < n2 + 1; i++) {
	}
	// Determine good roots
	ngood = nplus = nminus = 0;
	bad = 0;
	isgood = good[worst[ngood]];
	while (ngood < n2 + 1 && (nminus != nplus + n - 1 || isgood < dlmax || nplus == 0)) {
		//		if (isgood < dlmax) {
		if (Jacs[worst[ngood]] > 0) {
			nplus++;
		}
		else nminus++;
		ngood++;
		if(ngood<n2+1)	isgood = good[worst[ngood]];
	}

	Prov = new _curve;
	for (int i = 0; i < ngood; i++) {
		z = zr[worst[i]];
		zc = conj(z);
		Prov->append(z.re + s_offset->re, z.im + s_offset->im);
		Prov->last->dJ = Jacs[worst[i]];
		_MJacobians2
			if (theta->th >= 0) {
				_MJacobiansgood
			}
			else {
				_MJacobians0good
					corrquad += cq;
			}
			Prov->last->theta = theta;
	}
	corrquad2 = 1.e200;
	theta->errworst = 1.e100;
	for (int i = ngood; i < n2 + 1; i++) {
		z = zr[worst[i]];
		if (theta->th < 0) {
			zc = conj(z);
			_MJacobians2
				_MJacobians0bad
				if (cq < corrquad2) corrquad2 = cq;
		}
		else {
			for (int j = ngood; j < i; j++) { // Find the two ghost roots with minimum distance.
				ob2 = abs2(zr[worst[j]] - z);
				if (ob2 < theta->errworst) theta->errworst = ob2;
			}
		}
	}
	theta->errworst = sqrt(theta->errworst);
	theta->imlength = Prov->length;

	return Prov;
}

_curve* VBMicrolensing::NewImagesmultipoly(_theta* theta) {
	static complex  yc, z, zc, zo, delta, dy, dz, J2, J3, Jalt, Jaltc, JJalt2, LL, J1c2, dzita;
	static double dlmax = 1.0e-12, dzmax = 1.e-10, dJ2, ob2, cq, Jold, LLold;
	static int ngood, nplus, nminus, bad, isso, ncrit, igood, iter, iter2;
	static double mi, tst, isgood;
	static _curve* Prov;
	static _point* scan, * prin, * fifth, * left, * right, * center;

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	yc = conj(y);

#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif

	cmplx_roots_multigen(zr, coefs_mp, n2 + 1, false, false);

	for (int i = n2; i >= 0; i--) {
		findimagepoly(i);
		for (int j = n2; j > i; j--) {
			if (good[j] < 1.e9 && abs2(zr[j] - zr[i]) < dzmax && signbit(Jacs[j]) == signbit(Jacs[i])) {
				if (good[i] > good[j]) {
					good[i] = 1.e100;
					break;
				}
				else good[j] = 1.e100;
			}
		}
	}

	for (int i = 0; i < n2 + 1; i++) {
		int j = 0;
		while (j < i && good[worst[j]] < good[i]) j++; // worst contains root indices in increasing order of badness
		for (int k = i; k > j; k--) worst[k] = worst[k - 1];
		worst[j] = i;
	}
	// Determine good roots
	ngood = nplus = nminus = 0;
	bad = 0;
	isgood = good[worst[ngood]];
	while (ngood < n2 + 1 && (nminus != nplus + n - 1 || isgood < dlmax || nplus == 0)) {
		//		if (isgood < dlmax) {
		if (Jacs[worst[ngood]] > 0) {
			nplus++;
		}
		else nminus++;
		ngood++;
		if (ngood < n2 + 1)	isgood = good[worst[ngood]];
	}	
	Prov = new _curve;
	for (int i = 0; i < ngood; i++) {
		z = zr[worst[i]];
		zc = conj(z);
		Prov->append(z.re + s_offset ->re, z.im + s_offset->im);
		Prov->last->dJ = Jacs[worst[i]];
		_MJacobians2
			if (theta->th >= 0) {
				_MJacobiansgood
			}
			else {
				_MJacobians0good
					corrquad += cq;
			}
		Prov->last->theta = theta;
	}
	corrquad2 = 1.e200;
	theta->errworst = 1.e100;
	for (int i = ngood; i < n2 + 1; i++) {
		z = zr[worst[i]];
		if (theta->th < 0) {
			zc = conj(z);
			_MJacobians2
				_MJacobians0bad
				if (cq < corrquad2) corrquad2 = cq;
		}
		else {
			for (int j = ngood; j < i; j++) { // Find the two ghost roots with minimum distance.
				ob2 = abs2(zr[worst[j]] - z);
				if (ob2 < theta->errworst) theta->errworst = ob2;
			}
		}
	}
	theta->errworst = sqrt(theta->errworst);
	theta->imlength = Prov->length;

	return Prov;
}

void VBMicrolensing::OrderMultipleImages(_sols *Sols, _curve *Newpts) {
	static _point *scan, *scan2, *scan3, *isso[2];
	static _curve *scurve, *scurve2;

	static _theta *theta;
	static double th, mi, cmp, cmp2, cmp_2,er3, parab1, parab2;
	static int nprec, npres, npres2, nfoll, issoc[2], ij;

	nprec = nfoll = 0;
	theta = Newpts->first->theta;
	th = theta->th;
	theta->Mag = theta->prev->Mag = theta->maxerr = theta->prev->maxerr = 0;

	// Calcolo dell'errore per le ghost images.
	switch (SelectedMethod)
	{
	case Method::Singlepoly:
		if (theta->next->imlength == theta->prev->imlength) {
			mi = theta->next->errworst - theta->errworst;
			if (mi > 2 * theta->errworst && mi < 1.e-5) {
				theta->prev->maxerr = theta->maxerr = 1.e-10 / theta->errworst;
			}
			else {
				mi = theta->prev->errworst - theta->errworst;
				if (mi > 2 * theta->errworst && mi < 1.e-5) {
					theta->prev->maxerr = theta->maxerr = 1.e-10 / theta->errworst;
				}
			}
		}
		break;
	case Method::Multipoly:
		if (theta->next->imlength == theta->prev->imlength) {
			mi = theta->next->errworst - theta->errworst;
			if (mi > 2 * theta->errworst && mi < 1.e-5) {
				theta->prev->maxerr = theta->maxerr = 1.e-10 / theta->errworst;
			}
			else {
				mi = theta->prev->errworst - theta->errworst;
				if (mi > 2 * theta->errworst && mi < 1.e-5) {
					theta->prev->maxerr = theta->maxerr = 1.e-10 / theta->errworst;
				}
			}
		}
		break;
	case Method::Nopoly:
		break;
	}
	// Per ciascuna immagine troviamo il punto in cui inserire i nuovi punti
	scurve = Sols->first;
	for (int i = 0; i < Sols->length; i++) {
		if (th < scurve->first->theta->th) {
			if (th > scurve->first->theta->prev->prev->th) { // sembra strano ma è giusto prev->prev->th (l'inserimento di theta c'è già stato)
				cfoll[nfoll] = scurve; // immagine coinvolta all'inizio
				nfoll++;
				scurve2 = scurve->next;
				Sols->drop(scurve);
				i--;
				scurve = scurve2;
			}
			else {
				scurve = scurve->next;
			}
		}
		else {
			if (th > scurve->last->theta->th) {
				if (th < scurve->last->theta->next->next->th) {
					cprec[nprec] = scurve; // immagine coinvolta alla fine
					nprec++;
				}
			}
			else {
				// immagine coinvolta al centro
				scan = scurve->last;
				while (scan->theta->th > th) {
					scan = scan->prev;
				}
				cfoll[nfoll] = scurve->divide(scan);
				nfoll++;
				cprec[nprec] = scurve;
				nprec++;
			}
			scurve = scurve->next;
		}
	}
	npres = npres2 = Newpts->length;


	// Costruzione matrice distanze con immagini precedenti
	mi = 1.e100;
	for (int i = 0; i < nprec; i++) {
		scan = Newpts->first;
		for (int j = 0; j < npres; j++) {
			if (signbit(scan->dJ) == signbit(cprec[i]->last->dJ)) {
				A[i][j] = *(cprec[i]->last) - *scan;
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
					isso[1] = scan;
				}
			}
			else A[i][j] = 1.e100;

			scan = scan->next;
		}
	}

	//  Associazione con le immagini che precedono
#ifdef _PRINT_ERRORS
	printf("\nPreceding ordinary");
#endif
	while (nprec && npres) {
		scan = cprec[issoc[0]]->last;
		scan2 = isso[1];
		cmp2 = mi / fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		er3 = mi*mi*mi*samplingfactor;
		cmp = (scan->theta->th - scan2->theta->th);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp *0.0416666666666667;
		parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
		// Nuova correzione parabolica
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds)*mi / 2, fabs(scan->parab)*(cmp2) / 10, fabs(scan->parab)*(1.5*fabs(cmp2 / (cmp*cmp) - 1)));
#endif

		mi = fabs((scan->ds - scan2->ds)*mi *0.5) + fabs(scan->parab*1.5*fabs(cmp2 / (cmp_2)-1))+er3;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		scan->Mag = ((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + scan->parab);
		scan->err = mi;
		theta->prev->Mag += scan->Mag;
		theta->prev->maxerr += scan->err;

		Newpts->drop(isso[1]);
		cprec[issoc[0]]->append(isso[1]);
		cprec[issoc[0]]->partneratend = 0;

		nprec--;
		npres--;
		cpres[npres] = cprec[issoc[0]];

		for (int i = issoc[0]; i < nprec; i++) {
			cprec[i] = cprec[i + 1];
			for (int j = 0; j < npres + 1; j++) {
				A[i][j] = A[i + 1][j];
			}
		}
		for (int j = issoc[1]; j < npres; j++) {
			for (int i = 0; i < nprec; i++) {
				A[i][j] = A[i][j + 1];
			}
		}
		mi = 1.e100;
		for (int i = 0; i < nprec; i++) {
			scan = Newpts->first;
			for (int j = 0; j < npres; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
					isso[1] = scan;
				}
				scan = scan->next;
			}
		}
	}

	// Costruzione matrice distanze tra immagini in creazione
	mi = 1.e100;
	scan = Newpts->first;
	for (int i = 0; i < npres - 1; i++) {
		scan2 = scan->next;
		for (int j = i + 1; j < npres; j++) {
			if (signbit(scan->dJ) != signbit(scan2->dJ)) {
				A[i][j] = *scan2 - *scan;
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
					isso[0] = scan;
					isso[1] = scan2;
				}
			}
			else A[i][j] = 1.e100;
			scan2 = scan2->next;
		}
		scan = scan->next;
	}

	// Caso di creazione nuove immagini 
#ifdef _PRINT_ERRORS
	printf("\nPreceding creation");
#endif
	while (npres > 1 && mi<1.e99) {
		Newpts->drop(isso[0]);
		Newpts->drop(isso[1]);
		scurve = new _curve(isso[0]);
		isso[0]->prev = isso[0]->next = 0;
		scurve2 = new _curve(isso[1]);
		isso[1]->prev = isso[1]->next = 0;
		scurve->partneratstart = scurve2;
		scurve2->partneratstart = scurve;
		Sols->append(scurve);
		Sols->append(scurve2);
		npres -= 2;
		cpres[npres] = scurve;
		cpres[npres + 1] = scurve2;
		scan = isso[0];
		scan2 = isso[1];

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		cmp_2 = mi / cmp2;
		er3 = sqrt(mi*mi*mi)*samplingfactor;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp*0.04166666667;
		parab1 = -(-scan->ds + scan2->ds) * mi;
		parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scurve->parabstart = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scurve->parabstart, (scan->ds + scan2->ds)*mi / 2, fabs(scurve->parabstart)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2 * cmp*cmp2)*cmp);
#endif

		mi = fabs((scan->ds + scan2->ds)*mi *0.5) + er3 + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2. * cmp*cmp2)*cmp;
#ifdef _noparab
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scurve->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		scurve->Magstart = -(((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + scurve->parabstart));
		scurve->errstart = mi;
		scurve2->parabstart = -scurve->parabstart;
		scurve2->Magstart = 0;
		scurve2->errstart = 0;
		theta->prev->Mag += scurve->Magstart;
		theta->prev->maxerr += scurve->errstart;

		// Aggiornamento matrice distanze

		ij = 0;
		for (int i = 0; i < npres + 1; i++) {
			if (i == issoc[1]) ij++;
			for (int j = (ij>0)? i +1 : issoc[1]; j < npres + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
		}
		ij = 0;
		for (int i = 0; i < npres; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[0]; j < npres; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
		}
		mi = 1.e100;
		scan = Newpts->first;
		for (int i = 0; i < npres - 1; i++) {
			scan2 = scan->next;
			for (int j = i + 1; j < npres; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
					isso[0] = scan;
					isso[1] = scan2;
				}
				scan2 = scan2->next;
			}
			scan = scan->next;
		}
	}

	// Immagini spaiate
	while (npres >0) {
		scan = Newpts->first;
		Newpts->drop(scan);
		scurve = new _curve(scan);
		scan->prev = scan->next = 0;
		scurve->partneratstart = 0;
		Sols->append(scurve);
		npres--;
		cpres[npres] = scurve;

		scurve->parabstart = 0;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, 0, 0, 0, 0, 0);
#endif

#ifdef _noparab
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scurve->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scurve->Magstart = 0;
		scurve->errstart = 0;
	}

	delete Newpts;

	// Costruzione matrice distanze tra immagini in distruzione
	mi = 1.e100;
	for (int i = 0; i < nprec - 1; i++) {
		for (int j = i + 1; j < nprec; j++) {
			if (signbit(cprec[i]->last->dJ) != signbit(cprec[j]->last->dJ)) {
				A[i][j] = *(cprec[i]->last) - *(cprec[j]->last);
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
			else A[i][j] = 1.e100;
		}
	}

	// Caso di distruzione immagini 
#ifdef _PRINT_ERRORS
	printf("\nPreceding destruction");
#endif
	while (nprec >1 && mi<1.e99) {

		cprec[issoc[0]]->partneratend = cprec[issoc[1]];
		cprec[issoc[1]]->partneratend = cprec[issoc[0]];

		scan = cprec[issoc[0]]->last;
		scan2 = cprec[issoc[1]]->last;

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		er3 = sqrt(mi*mi*mi)*samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp *0.04166666666667;
		parab1 = -(scan->ds - scan2->ds) * mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds)*mi / 2, fabs(scan->parab)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2 * cmp*cmp2)*cmp);
#endif

		mi = fabs((scan->ds + scan2->ds)*mi *0.5) + er3 + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2.0 * cmp*cmp2)*cmp;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		scan->Mag = ((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + scan->parab);
		scan->err = mi;
		scan2->Mag = 0;
		scan2->err = 0;
		scan2->parab = -scan->parab;
		theta->prev->Mag += scan->Mag;
		theta->prev->maxerr += scan->err;

		nprec -= 2;
		// Aggiornamento matrice distanze
		ij = 0;
		for (int i = 0; i < nprec + 1; i++) {
			if (i == issoc[1]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[1]; j < nprec + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cprec[i] = cprec[i + ij];
		}
		ij = 0;
		for (int i = 0; i < nprec; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[0]; j < nprec; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cprec[i] = cprec[i + ij];
		}
		mi = 1.e100;
		for (int i = 0; i < nprec - 1; i++) {
			for (int j = i + 1; j < nprec; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}
	}

	// Immagini spaiate
	while (nprec > 0) {
		cprec[0]->partneratend = 0;
		scan = cprec[0]->last;

		scan->parab = 0;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, 0, 0, 0, 0, 0);
#endif

#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = 0;
		scan->err = 0;

		nprec--;
		for (int i = 0; i < nprec; i++) {
			cprec[i] = cprec[i + 1];
		}
	}

	// immagini seguenti
	npres = npres2;
	// Costruzione matrice distanze con immagini seguenti
	mi = 1.e100;
	for (int i = 0; i < npres; i++) {
		for (int j = 0; j < nfoll; j++) {
			if (signbit(cpres[i]->last->dJ) == signbit(cfoll[j]->first->dJ)) {
				A[i][j] = *(cpres[i]->last) - *(cfoll[j]->first);
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
			else A[i][j] = 1.e100;
		}
	}

	//  Associazione con le immagini che seguono
#ifdef _PRINT_ERRORS
	printf("\nFollowing ordinary");
#endif
	while (nfoll && npres) {
		scan = cpres[issoc[0]]->last;
		scan2 = cfoll[issoc[1]]->first;
		cmp2 = mi / fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		er3 = mi*mi*mi*samplingfactor;
		cmp = (scan->theta->th - scan2->theta->th);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp *0.041666666667;
		parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
			// Nuova correzione parabolica
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds)*mi / 2, fabs(scan->parab)*(cmp2) / 10, fabs(scan->parab)*(1.5*fabs(cmp2 / (cmp*cmp) - 1)));
#endif

		mi = fabs((scan->ds - scan2->ds)*mi *0.5) + fabs(scan->parab*1.5*fabs(cmp2 / (cmp_2)-1))+er3;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
		scan->Mag = ((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + scan->parab);
		scan->err = mi;
		theta->Mag += scan->Mag;
		theta->maxerr += scan->err;

		cpres[issoc[0]]->join(cfoll[issoc[1]]);

		npres--;
		nfoll--;
		for (int i = issoc[0]; i < npres; i++) {
			cpres[i] = cpres[i + 1];
			for (int j = 0; j < nfoll + 1; j++) {
				A[i][j] = A[i + 1][j];
			}
		}
		for (int j = issoc[1]; j < nfoll; j++) {
			cfoll[j] = cfoll[j + 1];
			for (int i = 0; i < npres; i++) {
				A[i][j] = A[i][j + 1];
			}
		}
		mi = 1.e100;
		for (int i = 0; i < npres; i++) {
			for (int j = 0; j < nfoll; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}
	}

	// Costruzione matrice distanze tra immagini in creazione
	mi = 1.e100;
	for (int i = 0; i < nfoll - 1; i++) {
		for (int j = i + 1; j < nfoll; j++) {
			if (signbit(cfoll[i]->first->dJ) != signbit(cfoll[j]->first->dJ)) {
				A[i][j] = *(cfoll[i]->first) - *(cfoll[j]->first);
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
			else A[i][j] = 1.e100;
		}
	}

	// Caso di creazione nuove immagini
#ifdef _PRINT_ERRORS
	printf("\nFollowing creation");
#endif
	while (nfoll>1 && mi < 1.e99) {

		cfoll[issoc[0]]->partneratstart = cfoll[issoc[1]];
		cfoll[issoc[1]]->partneratstart = cfoll[issoc[0]];

		scan = cfoll[issoc[0]]->first;
		scan2 = cfoll[issoc[1]]->first;

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		er3 = sqrt(mi*mi*mi)*samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp *0.04166666666666667;
		parab1 = (scan->ds - scan2->ds) * mi;
		parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		cfoll[issoc[0]]->parabstart = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, cfoll[issoc[0]]->parabstart, (scan->ds + scan2->ds)*mi / 2, fabs(cfoll[issoc[0]]->parabstart)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2 * cmp*cmp2)*cmp);
#endif
		mi = fabs((scan->ds + scan2->ds)*mi *0.5) + er3 + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) - 2.0 * cmp*cmp2)*cmp;
#ifdef _noparab
		mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		cfoll[issoc[0]]->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			cfoll[issoc[0]]->Magstart = -(((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + cfoll[issoc[0]]->parabstart));
		cfoll[issoc[0]]->errstart = mi;
		cfoll[issoc[1]]->parabstart = -cfoll[issoc[0]]->parabstart;
		cfoll[issoc[1]]->Magstart = 0;
		cfoll[issoc[1]]->errstart = 0;
		theta->Mag += cfoll[issoc[0]]->Magstart;
		theta->maxerr += cfoll[issoc[0]]->errstart;

		Sols->append(cfoll[issoc[0]]);
		Sols->append(cfoll[issoc[1]]);
		nfoll -= 2;

		// Aggiornamento matrice distanze
		ij = 0;
		for (int i = 0; i < nfoll + 1; i++) {
			if (i == issoc[1]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[1]; j < nfoll + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cfoll[i] = cfoll[i + ij];
		}
		ij = 0;
		for (int i = 0; i < nfoll; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[0]; j < nfoll; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cfoll[i] = cfoll[i + ij];
		}
		mi = 1.e100;
		for (int i = 0; i < nfoll - 1; i++) {
			for (int j = i + 1; j < nfoll; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}
	}

	// Immagini spaiate
	while (nfoll > 0) {

		cfoll[0]->partneratstart = 0;
		scan = cfoll[0]->first;

		cfoll[0]->parabstart = 0;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, 0, 0, 0, 0, 0);
#endif
#ifdef _noparab
		mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		cfoll[issoc[0]]->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			cfoll[0]->Magstart = 0;
		cfoll[0]->errstart = 0;

		Sols->append(cfoll[0]);
		nfoll--;
		for (int i = 0; i < nfoll; i++) {
			cfoll[i] = cfoll[i + 1];
		}
	}

	// Costruzione matrice distanze tra immagini in distruzione
	mi = 1.e100;
	for (int i = 0; i < npres - 1; i++) {
		for (int j = i + 1; j < npres; j++) {
			if (signbit(cpres[i]->last->dJ) != signbit(cpres[j]->last->dJ)) {
				A[i][j] = *(cpres[i]->last) - *(cpres[j]->last);
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
			else A[i][j] = 1.e100;
		}
	}


	// Caso di distruzione immagini
#ifdef _PRINT_ERRORS
	printf("\nFollowing destruction");
#endif
	while (npres>1 && mi<1.e99) {
		cpres[issoc[0]]->partneratend = cpres[issoc[1]];
		cpres[issoc[1]]->partneratend = cpres[issoc[0]];

		scan = cpres[issoc[0]]->last;
		scan2 = cpres[issoc[1]]->last;

		cmp2 = fabs(scan->d.re*scan2->d.re + scan->d.im*scan2->d.im);
		er3 = sqrt(mi*mi*mi)*samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp *0.0416666666667;
		parab1 = -(scan->ds - scan2->ds) * mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds)*mi / 2, fabs(scan->parab)*(cmp*cmp) / 10, 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2 * cmp*cmp2)*cmp);
#endif

		mi = fabs((scan->ds + scan2->ds)*mi *0.5) + er3 + 1.5*fabs(((scan->d.re - scan2->d.re)*(scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im)*(scan->x2 - scan2->x2)) + 2.0 * cmp*cmp2)*cmp;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = ((scan->dJ > 0) ? -1 : 1)*((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) *0.5 + scan->parab);
		scan->err = mi;
		scan2->parab = -scan->parab;
		scan2->Mag = 0;
		scan2->err = 0;
		theta->Mag += scan->Mag;
		theta->maxerr += scan->err;

		npres -= 2;
		// Aggiornamento matrice distanze
		ij = 0;
		for (int i = 0; i < npres + 1; i++) {
			if (i == issoc[1]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[1]; j < npres + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cpres[i] = cpres[i + ij];
		}
		ij = 0;
		for (int i = 0; i < npres; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij>0) ? i + 1 : issoc[0]; j < npres; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cpres[i] = cpres[i + ij];
		}
		mi = 1.e100;
		for (int i = 0; i < npres - 1; i++) {
			for (int j = i + 1; j < npres; j++) {
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}

	}

	// Immagini spaiate
	while (npres > 0) {
		cpres[0]->partneratend = 0;

		scan = cpres[0]->last;

		scan->parab = 0;

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, 0, 0, 0, 0, 0);
#endif

#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2)*(scan2->x1 - scan->x1) / 2 * cmp*cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = 0;
		scan->err = 0;

		npres--;
		for (int i = 0; i < npres; i++) {
			cpres[i] = cpres[i + 1];
		}
	}

}
