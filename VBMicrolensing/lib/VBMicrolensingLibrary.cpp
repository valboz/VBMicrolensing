// VBMicorlensing v0.0 (2018)
//

#define _CRT_SECURE_NO_WARNINGS


#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _PRINT_ERRORS2
#define _PRINT_ERRORS

#ifndef __unmanaged
using namespace VBMicrolensingLibrary;
#endif

//////////////////////////////
//////////////////////////////
////////Constructor and destructor
//////////////////////////////
//////////////////////////////

VBMicrolensing::VBMicrolensing() {
	Obj[0] = -0.0397317;
	Obj[1] = 0.998164;
	Obj[2] = -0.045714;
	// reference is ecliptic with x-axis toward the equinox.
	// axial tilt at J2000 is 23:26:21.406 from JPL fundamental ephemeris
	Eq2000[0] = 1;
	Eq2000[1] = Eq2000[2] = Quad2000[0] = North2000[0] = 0;
	Quad2000[1] = 0.9174820003578725;
	Quad2000[2] = -0.3977772982704228;
	North2000[1] = 0.3977772982704228;
	North2000[2] = 0.9174820003578725;
	t0old = 0.;
	Tol = 1.e-2;
	RelTol = 0;
	tsat = 0;
	possat = 0;
	nsat = 0;
	ndatasat = 0;
	satellite = 0;
	parallaxsystem = 0;
	t0_par_fixed = -1;
	t0_par = 7000;
	minannuli = 1;
	curLDprofile = LDlinear;
	a1 = 0;
	npLD = 0;
	LDtab = rCLDtab = CLDtab = 0;
	n=0;
	zr = zcr= pza = pdum = ppy =a=coefs=0;
	s = s_sort = 0;
	m_mp = 0;
	zr_mp= coefs_mp=a_mp=ppy_mp=pza_mp = 0;
	good = Jacs = m=0;
	pmza = pyaza = ppmy =  0;
	pmza_mp = pyaza_mp = ppmy_mp = 0;
	dist_mp = 0;
	nrootsmp_mp = 0;
	y_mp = 0;
	init = 0;
	centralimages = 0;
	errs = 0;
	newseeds = 0;
	grads = 0;
	S2s= S3s= S4s = 0;
	s_offset = new complex;
	q = q_sort = 0;
	A = 0;
	cprec = cpres = cfoll = 0;
	worst=0;
	pert = 0;
	Mag0 = 0;
	NPcrit = 200;
	ESPLoff = true;
	multidark = false;
	astrometry = false;
	mass_luminosity_exponent = 4.0;
	mass_radius_exponent = 0.9;
	rootaccuracy = 9.e-22;
	samplingfactor = 0.125;
	squarecheck = false;
	CumulativeFunction = &VBDefaultCumulativeFunction;
	SelectedMethod = Method::Nopoly;
}

VBMicrolensing::~VBMicrolensing() {
	if (nsat) {
		for (int i = 0; i<nsat; i++) {
			for (int j = 0; j<ndatasat[i]; j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}

	if (m) {
		free(m);
		free(a);
	}
	if (zr) {
		free(zr);
		free(zcr);
		free(good);
		free(worst);
		free(pert);
		free(Jacs);
		free(zaltc);
		free(J1);
		free(J1c);
		free(prodevs);
		free(devs);
		free(init);
		free(centralimages);
		free(errs);
		free(newseeds);
		free(grads);
		free(S2s);
		free(S3s);
		free(S4s);
	}
	if (pmza) {
		for (int i = 0; i < n; i++) {
			free(pmza[i]);
			free(pmza2[i]);
			free(pyaza[i]);
			free(ppmy[i]);
			free(za[i]);
			free(za2[i]);
		}
		free(pmza);
		free(pmza2);
		free(pyaza);
		free(ppmy);
		free(pza);
		free(pza2);
		free(pdum);
		free(ppy);
		free(za);
		free(za2);
	}
	if (A) {
		for (int i = 0; i < nroots; i++) {
			free(A[i]);
		}
		free(A);
		free(cprec);
		free(cpres);
		free(cfoll);
	}
	if (npLD > 0) {
		free(LDtab);
		free(rCLDtab);
	}
	//multipoly
	if (coefs_mp) {
		for (int i = 0; i < n; i++) {
			free(coefs_mp[i]);
		}
		free(coefs_mp);
	}
	if (m_mp) {
		for (int i = 0; i < n; i++) {
			free(m_mp[i]);
			free(a_mp[i]);
		}
		free(m_mp);
		free(a_mp);
		free(q_sort);
		free(s_sort);
		free(y_mp);
	}
	if (pmza_mp) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				free(pmza_mp[j][i]);
				free(pyaza_mp[j][i]);
				free(ppmy_mp[j][i]);
			}
		}
		for (int i = 0; i < n; i++) {
			free(pmza_mp[i]);
			free(pyaza_mp[i]);
			free(ppmy_mp[i]);

			free(ppy_mp[i]);
			free(pza_mp[i]);

		}
		free(pmza_mp);
		free(pyaza_mp);
		free(ppmy_mp);
		free(pza_mp);
		free(ppy_mp);
	}
	if (zr_mp) {
		for (int j = 0; j < n; j++) {
			free(zr_mp[j]);
		}
		free(zr_mp);
		free(nrootsmp_mp);
	}

	//delete s_offset;

}

