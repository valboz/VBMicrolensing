// VBMicrolensing v4.1.1 (2024)
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
// GNU Lesser General Public License applies to all parts of this code.
// Please read the separate LICENSE.txt file for more details.

#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
char systemslash = '\\';
#else
char systemslash = '/';
#endif


#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define _PRINT_ERRORS2
//#define _PRINT_ERRORS
//#define _PRINT_ERRORS_DARK

#pragma region Constructor/destructor

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
	n = 0;
	zr = zcr = pza = pdum = ppy = a = coefs = 0;
	s = s_sort = 0;
	m_mp = 0;
	zr_mp = coefs_mp = a_mp = ppy_mp = pza_mp = 0;
	good = Jacs = m = 0;
	pmza = pyaza = ppmy = 0;
	pmza_mp = pyaza_mp = ppmy_mp = 0;
	dist_mp = 0;
	nrootsmp_mp = 0;
	y_mp = 0;
	init = 0;
	centralimages = 0;
	errs = 0;
	newseeds = 0;
	grads = 0;
	S2s = S3s = S4s = 0;
	s_offset = new complex;
	q = q_sort = 0;
	A = 0;
	cprec = cpres = cfoll = 0;
	worst = 0;
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
		for (int i = 0; i < nsat; i++) {
			for (int j = 0; j < ndatasat[i]; j++) free(possat[i][j]);
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

#pragma endregion

#pragma region single-source-mag


void VBMicrolensing::LoadESPLTable(char* filename) {
	FILE* f;

	if ((f = fopen(filename, "rb")) != 0) {
		fread(ESPLin, sizeof(double), __rsize * __zsize, f);
		fread(ESPLout, sizeof(double), __rsize * __zsize, f);
		fread(ESPLinastro, sizeof(double), __rsize * __zsize, f);
		fread(ESPLoutastro, sizeof(double), __rsize * __zsize, f);
		fclose(f);
		ESPLoff = false;
	}
	else {
		printf("\nESPL table not found !");
	}
}

double VBMicrolensing::PSPLMag(double u) {
	static double u2, u22;
	u2 = u * u;
	u22 = u2 + 2;
	if (astrometry) {
		astrox1 = u + u / u22;
	}
	return  u22 / sqrt(u2 * (u2 + 4));
}

double VBMicrolensing::ESPLMag(double u, double RSv) {
	double mag, z, fr, cz, cr, u2;
	int iz, ir;

	if (ESPLoff) {
		printf("\nLoad ESPL table first!");
		return 0;
	}

	fr = -10.857362047581296 * log(0.01 * RSv);
	if (fr > __rsize - 1) fr = __rsize - 1.000001;
	if (fr < 0) printf("Source too large!");
	ir = (int)floor(fr);
	fr -= ir;
	cr = 1 - fr;

	z = u / RSv;

	if (z < 1) {
		z *= __zsize - 1;
		iz = (int)floor(z);
		z -= iz;
		cz = 1 - z;
		mag = sqrt(1 + 4. / (RSv * RSv));
		mag *= ESPLin[ir][iz] * cr * cz + ESPLin[ir + 1][iz] * fr * cz + ESPLin[ir][iz + 1] * cr * z + ESPLin[ir + 1][iz + 1] * fr * z;
		if (astrometry) {
			astrox1 = (1 - 1. / (4 + RSv * RSv)) * u;
			astrox1 *= ESPLinastro[ir][iz] * cr * cz + ESPLinastro[ir + 1][iz] * fr * cz + ESPLinastro[ir][iz + 1] * cr * z + ESPLinastro[ir + 1][iz + 1] * fr * z;
		}
	}
	else {
		z = 0.99999999999999 / z;
		z *= __zsize - 1;
		iz = (int)floor(z);
		z -= iz;
		cz = 1 - z;

		u2 = u * u;
		mag = (u2 + 2) / sqrt(u2 * (u2 + 4));
		mag *= ESPLout[ir][iz] * cr * cz + ESPLout[ir + 1][iz] * fr * cz + ESPLout[ir][iz + 1] * cr * z + ESPLout[ir + 1][iz + 1] * fr * z;
		if (astrometry) {
			astrox1 = u * (u2 + 3) / (u2 + 2);
			astrox1 *= ESPLoutastro[ir][iz] * cr * cz + ESPLoutastro[ir + 1][iz] * fr * cz + ESPLoutastro[ir][iz + 1] * cr * z + ESPLoutastro[ir + 1][iz + 1] * fr * z;
		}
	}

	return mag;
}

double VBMicrolensing::ESPLMag2(double u, double rho) {
	double Mag, u2, u6, rho2Tol;
	int c = 0;

	//Tol1_4 = sqrt(2 / Tol);
	//u2 = u*u;
	//u3Tol = u2*u*Tol;

	//if (u2 < Tol1_4) {
	//	rho2 = rho*rho;
	//	if (u3Tol > rho2*(1 + Tol1_4*rho)) {

	u2 = u * u;
	rho2Tol = rho * rho / Tol;
	u6 = u2 * u2 * u2;

	if (u6 * (1 + 0.003 * rho2Tol) > 0.027680640625 * rho2Tol * rho2Tol) {
		Mag = (u2 + 2) / (u * sqrt(u2 + 4));
		if (astrometry) {
			astrox1 = u * (1 + 1 / (u2 + 2));
		}
	}
	else {
		Mag = ESPLMagDark(u, rho);
	}
	Mag0 = 0;
	return Mag;
}

double VBMicrolensing::ESPLMagDark(double u, double RSv) {
	double Mag = -1.0, Magold = 0., Tolv = Tol;
	double tc, rb, lc, rc, cb, u2;
	int c = 0, flag;
	double currerr, maxerr;
	annulus* first, * scan, * scan2;
	int nannold, totNPS = 1;
	double LDastrox1 = 0.0;

	while ((Mag < 0.9) && (c < 3)) {

		first = new annulus;
		first->bin = 0.;
		first->cum = 0.;

		u2 = u * u;
		first->Mag = Mag0 = (u2 + 2) / (u * sqrt(u2 + 4));
		first->nim = 2;
		if (astrometry) {
			astrox1 = u * (u2 + 3) / (u2 + 2);
			first->LDastrox1 = astrox1 * first->Mag;
		}

		scr2 = sscr2 = 0;
		first->f = LDprofile(0);
		first->err = 0;
		first->prev = 0;


		first->next = new annulus;
		scan = first->next;
		scan->prev = first;
		scan->next = 0;
		scan->bin = 1.;
		scan->cum = 1.;
		scan->Mag = ESPLMag(u, RSv);//ESPLMag(u, RSv, Tolv, &Images);
		if (astrometry) {
			scan->LDastrox1 = astrox1 * scan->Mag;
		}
		scan->nim = 2;
		scr2 = sscr2 = 1;
		scan->f = LDprofile(0.9999999);
		scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) / 4);

		Magold = Mag = scan->Mag;
		if (astrometry) {
			LDastrox1 = scan->LDastrox1;
		}
		//	
		//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
		currerr = scan->err;
		flag = 0;
		nannuli = nannold = 1;
		while (((flag < nannold + 5) && (currerr > Tolv) && (currerr > RelTol * Mag)) || (nannuli < minannuli)) {
			maxerr = 0;
			for (scan2 = first->next; scan2; scan2 = scan2->next) {
#ifdef _PRINT_ERRORS_DARK
				printf("\n%d %lf %le | %lf %le", nannuli, scan2->Mag, scan2->err, Mag, currerr);
#endif
				if (scan2->err > maxerr) {
					maxerr = scan2->err;
					scan = scan2;
				}
			}

			nannuli++;
			Magold = Mag;
			Mag -= (scan->Mag * scan->bin * scan->bin - scan->prev->Mag * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			if (astrometry) {
				LDastrox1 -= (scan->LDastrox1 * scan->bin * scan->bin - scan->prev->LDastrox1 * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);

			}
			currerr -= scan->err;
			rc = scan->cum;
			lc = scan->prev->cum;
			tc = (lc + rc) / 2;
			cb = rCLDprofile(tc, scan->prev, scan);

			scan->prev->next = new annulus;
			scan->prev->next->prev = scan->prev;
			scan->prev = scan->prev->next;
			scan->prev->next = scan;
			scan->prev->bin = cb;
			scan->prev->cum = tc;
			scan->prev->f = LDprofile(cb);
			scan->prev->Mag = ESPLMag(u, RSv * cb);
			if (astrometry) {
				scan->prev->LDastrox1 = astrox1 * scan->prev->Mag;

			}
			scan->prev->nim = 2;
			scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin) / 4);
			scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin) / 4);
			rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
			scan->prev->err += fabs(rb * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin));
			scan->err += fabs(rb * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin));

			Mag += (scan->bin * scan->bin * scan->Mag - cb * cb * scan->prev->Mag) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			Mag += (cb * cb * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
			if (astrometry) {
				LDastrox1 += (scan->bin * scan->bin * scan->LDastrox1 - cb * cb * scan->prev->LDastrox1) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
				LDastrox1 += (cb * cb * scan->prev->LDastrox1 - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->LDastrox1) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);

			}
			currerr += scan->err + scan->prev->err;

			if (fabs(Magold - Mag) * 2 < Tolv) {
				flag++;
			}
			else {
				flag = 0;
				nannold = nannuli;
			}

		}

		while (first) {
			scan = first->next;
			delete first;
			first = scan;
		}

		Tolv /= 10;
		c++;
	}
	therr = currerr;
	if (astrometry) {
		LDastrox1 /= Mag;
		astrox1 = LDastrox1;

	}
	return Mag;
}

#pragma endregion

#pragma region binary-mag


double VBMicrolensing::BinaryMag0(double a1, double q1, double y1v, double y2v, _sols** Images) {
	static complex a, q, m1, m2, y;
	static double av = -1.0, qv = -1.0;
	static complex  coefs[24], d1, d2, dy, dJ, dz;
	static double Mag, Ai;

	static _theta* stheta;
	static _curve* Prov, * Prov2;
	static _point* scan1, * scan2;

	Mag = Ai = -1.0;
	stheta = new _theta(-1.);
	if ((a1 != av) || (q1 != qv)) {
		av = a1;
		qv = q1;
		if (q1 < 1) {
			a = complex(-a1, 0);
			q = complex(q1, 0);
		}
		else {
			a = complex(a1, 0);
			q = complex(1 / q1, 0);
		}
		m1 = 1.0 / (1.0 + q);
		m2 = q * m1;

		coefs[20] = a;
		coefs[21] = m1;
		coefs[22] = m2;
		coefs[6] = a * a;
		coefs[7] = coefs[6] * a;
		coefs[8] = m2 * m2;
		coefs[9] = coefs[6] * coefs[8];
		coefs[10] = a * m2;
		coefs[11] = a * m1;
		coefs[23] = 0;

	}
	y = complex(y1v, y2v);
	(*Images) = new _sols;
	corrquad = corrquad2 = 0;
	safedist = 10;
	Prov = NewImages(y, coefs, stheta);
	if (Prov->length == 0) {
		delete Prov;
		delete stheta;
		return -1;
	}
	if (q.re < 0.01) {
		safedist = y1v + coefs[11].re - 1 / a.re;
		safedist *= safedist;
		safedist += y2v * y2v - 36 * q1 / (a1 * a1);
	}
	Mag = 0.;
	astrox1 = 0.;
	astrox2 = 0.;
	nim0 = 0;
	for (scan1 = Prov->first; scan1; scan1 = scan2) {
		scan2 = scan1->next;
		Prov2 = new _curve(scan1);
		//Prov2->append(scan1->x1, scan1->x2);
		//Prov2->last->theta = stheta;
		//Prov2->last->d = Prov2->first->d;
		//Prov2->last->dJ = Prov2->first->dJ;
		//Prov2->last->J2 = Prov2->first->J2;
		//Prov2->last->ds = Prov2->first->ds;
		(*Images)->append(Prov2);
		Ai = fabs(1 / scan1->dJ);
		Mag += Ai;
		if (astrometry) {
			astrox1 += scan1->x1 * Ai;
			astrox2 += (scan1->x2) * Ai;
		}
		nim0++;
	}
	Prov->length = 0;
	delete Prov;
	delete stheta;
	if (astrometry) {
		astrox1 /= (Mag);
		astrox1 -= coefs[11].re;
		astrox2 /= (Mag);
	}
	NPS = 1;
	return Mag;

}

double VBMicrolensing::BinaryMag0(double a1, double q1, double y1v, double y2v) {
	static _sols* images;
	static double mag;
	mag = BinaryMag0(a1, q1, y1v, y2v, &images);
	delete images;
	return mag;
}

double VBMicrolensing::BinaryMagSafe(double s, double q, double y1v, double y2v, double RS, _sols** images) {
	static double Mag, mag1, mag2, RSi, RSo, delta1, delta2;
	static int NPSsafe;
	Mag = BinaryMag(s, q, y1v, y2v, RS, Tol, images);
	RSi = RS;
	RSo = RS;
	NPSsafe = NPS;
	if (Mag < 0 || Mag * RS > 3) {
		mag1 = -1;
		delta1 = 3.33333333e-8;
		while ((mag1 < 0.1 || mag1 * RSi > 3) && RSi >= 0) {
			delete* images;
			delta1 *= 3.;
			RSi = RS - delta1;
			mag1 = (RSi > 0) ? BinaryMag(s, q, y1v, y2v, RSi, Tol, images) : BinaryMag0(s, q, y1v, y2v, images);
			//			printf("\n-safe1 %lf %lf %d", RSi, mag1, NPS);
			NPSsafe += NPS;
		}
		if (mag1 < 0) mag1 = 1.0;
		mag2 = -1;
		delta2 = 3.33333333e-8;
		while (mag2 < 0.1 || mag2 * RSo > 3) {
			delta2 *= 3.;
			RSo = RS + delta2;
			delete* images;
			mag2 = BinaryMag(s, q, y1v, y2v, RSo, Tol, images);
			//			printf("\n-safe2 %lf %lf %d", RSo,mag2,NPS);
			NPSsafe += NPS;
		}
		Mag = (mag1 * delta2 + mag2 * delta1) / (delta1 + delta2);
	}
	NPS = NPSsafe;

	return Mag;
}

double VBMicrolensing::BinaryMag(double a1, double q1, double y1v, double y2v, double RSv, double Tol, _sols** Images) {
	static complex a, q, m1, m2, y0, y, yc, z, zc;
	static double av = -1.0, qv = -1.0;
	static complex coefs[24], d1, d2, dy, dJ, dz;
	static double thoff = 0.01020304, errbuff;
	static double Mag, th;
	////////////////////////////  
	static double errimage, maxerr, currerr, Magold;
	static int NPSmax, flag, NPSold, flagbad, flagbadmax = 3;
	static _curve* Prov, * Prov2;
	static _point* scan1, * scan2;
	static _thetas* Thetas;
	static _theta* stheta, * itheta;

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	// Initialization of the equation coefficients

	if ((a1 != av) || (q1 != qv)) {
		av = a1;
		qv = q1;
		if (q1 < 1) {
			a = complex(-a1, 0);
			q = complex(q1, 0);
		}
		else {
			a = complex(a1, 0);
			q = complex(1 / q1, 0);
		}
		m1 = 1.0 / (1.0 + q);
		m2 = q * m1;

		coefs[20] = a;
		coefs[21] = m1;
		coefs[22] = m2;
		coefs[6] = a * a;
		coefs[7] = coefs[6] * a;
		coefs[8] = m2 * m2;
		coefs[9] = coefs[6] * coefs[8];
		coefs[10] = a * m2;
		coefs[11] = a * m1;

	}
	coefs[23] = RSv;

	y0 = complex(y1v, y2v);
	NPS = 1;
	if (Tol > 1.) {
		errimage = 0.;
		NPSmax = (int)(Tol);
	}
	else {
		errimage = Tol * M_PI * RSv * RSv;
		NPSmax = 10000; // era 32000
	}
	errbuff = 0;

	// Calculation of the images

	(*Images) = new _sols;
	Thetas = new _thetas;
	th = thoff;
	stheta = Thetas->insert(th);
	stheta->maxerr = 1.e100;
	y = y0 + complex(RSv * cos(thoff), RSv * sin(thoff));


#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif
	flag = 0;
	flagbad = 0;
	while (flag == 0) {
		Prov = NewImages(y, coefs, stheta);
		if (Prov->length > 0) {
			flag = 1;
		}
		else {
			delete Prov;
			stheta->th += 0.01;
			if (stheta->th > 2.0 * M_PI) {
				delete Thetas;
				return -1;
			}
			y = y0 + complex(RSv * cos(stheta->th), RSv * sin(stheta->th));
		}
	}
#ifdef _PRINT_TIMES
	tim1 = Environment::TickCount;
	GM += tim1 - tim0;
#endif
	stheta = Thetas->insert(2.0 * M_PI + Thetas->first->th);
	stheta->maxerr = 0.;
	stheta->Mag = 0.;
	stheta->astrox1 = 0.;
	stheta->astrox2 = 0.;
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

	th = M_PI + Thetas->first->th;
	flag = 0;
	Magold = -1.;
	NPSold = 2;
	currerr = 1.e100;
	do {
		stheta = Thetas->insert(th);
		y = y0 + complex(RSv * cos(th), RSv * sin(th));
#ifdef _PRINT_TIMES
		tim0 = Environment::TickCount;
#endif
		//if (NPS == 422) {
		//	NPS = NPS;
		//}

		Prov = NewImages(y, coefs, stheta);
#ifdef _PRINT_TIMES
		tim1 = Environment::TickCount;
		GM += tim1 - tim0;
#endif
		if (Prov->length > 0) {
			flagbad = 0;
			OrderImages((*Images), Prov);
			if ((stheta->th - stheta->prev->th) * RSv < 1.e-11/* || stheta->maxerr > jumperrfactor * currerr || stheta->prev->maxerr > jumperrfactor * currerr*/) {
				errbuff += stheta->maxerr + stheta->prev->maxerr;
				stheta->maxerr = 0;
				stheta->prev->maxerr = 0;
			}
		}
		else {
			delete Prov;
			flagbad++;
			if (flagbad == flagbadmax) {
				if (NPS < 16) {
					delete Thetas;
					return -1;
				}
				errbuff += stheta->prev->maxerr;
				stheta->prev->maxerr = 0;
				NPS--;
				NPSmax--;
			}
			else {
				th = (th - stheta->prev->th >= stheta->next->th - th) ? (th + flagbad * stheta->prev->th) / (1 + flagbad) : (th + flagbad * stheta->next->th) / (1 + flagbad);
			}
			Thetas->remove(stheta);
		}

		if (flagbad == 0 || flagbad == flagbadmax) {
			maxerr = currerr = Mag = 0.;

			astrox1 = astrox2 = 0.;
			stheta = Thetas->first;

			while (stheta->next) {
				currerr += stheta->maxerr;
				Mag += stheta->Mag;

				if (astrometry) {
					astrox1 += stheta->astrox1;
					astrox2 += stheta->astrox2;
				}
#ifndef _uniform
				if (stheta->maxerr > maxerr) {
					maxerr = stheta->maxerr;
#else
				if (stheta->next->th * 0.99999 - stheta->th > maxerr) {
					maxerr = stheta->next->th - stheta->th;
#endif
					itheta = stheta;
				}
				stheta = stheta->next;
#ifdef _selectimage
				if ((NPS == NPSmax - 1) && (fabs(floor(stheta->th / M_PI * _npoints / 2 + 0.5) - stheta->th / M_PI * _npoints / 2) < 1.e-8)) {
					printf("%d %.15le\n", (int)floor(stheta->th / M_PI * _npoints / 2 + 0.5), Mag);
				}
#endif
				}
			th = (itheta->th + itheta->next->th) / 2;
			NPS++;
#ifndef _uniform
			if (fabs(Magold - Mag) * 2 < errimage) {
				flag++;
			}
			else {
				flag = 0;
				Magold = Mag;
				NPSold = NPS + 8;
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
			printf("\nNPS= %d Mag = %lf maxerr= %lg currerr =%lg th = %lf", NPS, Mag / (M_PI * RSv * RSv), maxerr / (M_PI * RSv * RSv), currerr / (M_PI * RSv * RSv), th);
#endif
			}
		} while ((currerr > errimage) && (currerr > RelTol * Mag) && (NPS < NPSmax) && ((flag < NPSold)/* || NPS<8 ||(currerr>10*errimage)*/)/*&&(flagits)*/);
		if (astrometry) {
			astrox1 /= (Mag);
			astrox2 /= (Mag);
		}
		Mag /= (M_PI * RSv * RSv);
		therr = (currerr + errbuff) / (M_PI * RSv * RSv);

		delete Thetas;
		//	if (NPS == NPSmax) return 1.e100*Tol; // Only for testing
		return Mag;

	}


double VBMicrolensing::BinaryMag(double a1, double q1, double y1v, double y2v, double RSv, double Tol) {
	static _sols* images;
	static double mag;
	mag = BinaryMag(a1, q1, y1v, y2v, RSv, Tol, &images);
	delete images;
	return mag;
}

double VBMicrolensing::BinaryMag2(double s, double q, double y1v, double y2v, double rho) {
	static double Mag, rho2, y2a;//, sms , dy1, dy2;
	static int c;
	static _sols* Images;

	c = 0;

	y2a = fabs(y2v);

	Mag0 = BinaryMag0(s, q, y1v, y2a, &Images);
	delete Images;
	rho2 = rho * rho;
	corrquad *= 6 * (rho2 + 1.e-4 * Tol);
	corrquad2 *= (rho + 1.e-3);
	if (corrquad < Tol && corrquad2 < 1 && (/*rho2 * s * s<q || */ safedist > 4 * rho2)) {
		Mag = Mag0;
	}
	else {
		Mag = BinaryMagDark(s, q, y1v, y2a, rho, Tol);
	}
	Mag0 = 0;

	if (y2v < 0) {
		y_2 = y2v;
		astrox2 = -astrox2;
	}
	return Mag;
}

double VBDefaultCumulativeFunction(double cb, double* a1) {
	static double r2, cr2, scr2, cc;
	r2 = cb * cb;
	cr2 = 1 - r2;
	scr2 = sqrt(cr2);
	cc = (3 * r2 * (1 - *a1) - 2 * (*a1) * (scr2 * cr2 - 1)) / (3 - *a1);
	return cc;
}

double VBMicrolensing::BinaryMagDark(double a, double q, double y1, double y2, double RSv, double Tolnew) {
	static double Mag, Magold, Tolv;
	static double LDastrox1, LDastrox2;
	static double tc, lc, rc, cb, rb;
	static int c, flag;
	static double currerr, maxerr;
	static annulus* first, * scan, * scan2;
	static int nannold, totNPS;
	static _sols* Images;

	Mag = -1.0;
	Magold = 0.;
	Tolv = Tol;
	LDastrox1 = LDastrox2 = 0.0;
	c = 0;
	totNPS = 1;

	Tol = Tolnew;
	y_1 = y1;
	y_2 = y2;
	while ((Mag < 0.9) && (c < 3)) {

		first = new annulus;
		first->bin = 0.;
		first->cum = 0.;
		if (Mag0 > 0.5) {
			first->Mag = Mag0;
			first->nim = nim0;
		}
		else {
			first->Mag = BinaryMag0(a, q, y_1, y_2, &Images);
			first->nim = Images->length;
			delete Images;
		}
		if (astrometry) {
			first->LDastrox1 = astrox1 * first->Mag;
			first->LDastrox2 = astrox2 * first->Mag;
		}
		scr2 = sscr2 = 0;
		first->f = LDprofile(0);
		first->err = 0;
		first->prev = 0;


		first->next = new annulus;
		scan = first->next;
		scan->prev = first;
		scan->next = 0;
		scan->bin = 1.;
		scan->cum = 1.;
		scan->Mag = BinaryMagSafe(a, q, y_1, y_2, RSv, &Images);
		if (astrometry) {
			scan->LDastrox1 = astrox1 * scan->Mag;
			scan->LDastrox2 = astrox2 * scan->Mag;
		}
		totNPS += NPS;
		scan->nim = Images->length;
		delete Images;
		scr2 = sscr2 = 1;
		scan->f = LDprofile(0.9999999);
		if (scan->nim == scan->prev->nim) {
			scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
		}
		else {
			scan->err = fabs((scan->Mag) * (scan->prev->f - scan->f) / 4);
		}

		Magold = Mag = scan->Mag;
		if (astrometry) {
			LDastrox1 = scan->LDastrox1;
			LDastrox2 = scan->LDastrox2;
		}
		//			scan->err+=scan->Mag*Tolv*0.25; //Impose calculation of intermediate annulus at mag>4. Why?
		currerr = scan->err;
		flag = 0;
		nannuli = nannold = 1;
		while (((flag < nannold + 5) && (currerr > Tolv) && (currerr > RelTol * Mag)) || (nannuli < minannuli)) {
			maxerr = 0;
			for (scan2 = first->next; scan2; scan2 = scan2->next) {
#ifdef _PRINT_ERRORS_DARK
				printf("\n%d %lf %le | %lf %le", nannuli, scan2->Mag, scan2->err, Mag, currerr);
#endif
				if (scan2->err > maxerr) {
					maxerr = scan2->err;
					scan = scan2;
				}
			}

			nannuli++;
			Magold = Mag;
			Mag -= (scan->Mag * scan->bin * scan->bin - scan->prev->Mag * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			if (astrometry) {
				LDastrox1 -= (scan->LDastrox1 * scan->bin * scan->bin - scan->prev->LDastrox1 * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
				LDastrox2 -= (scan->LDastrox2 * scan->bin * scan->bin - scan->prev->LDastrox2 * scan->prev->bin * scan->prev->bin) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			}
			currerr -= scan->err;
			lc = scan->prev->cum;
			rc = scan->cum;
			tc = (lc + rc) * 0.5;
			cb = rCLDprofile(tc, scan->prev, scan);
			scan->prev->next = new annulus;
			scan->prev->next->prev = scan->prev;
			scan->prev = scan->prev->next;
			scan->prev->next = scan;
			scan->prev->bin = cb;
			scan->prev->cum = tc;
			scan->prev->f = LDprofile(cb);
			scan->prev->Mag = BinaryMagSafe(a, q, y_1, y_2, RSv * cb, &Images);
			if (astrometry) {
				scan->prev->LDastrox1 = astrox1 * scan->prev->Mag;
				scan->prev->LDastrox2 = astrox2 * scan->prev->Mag;
			}
			totNPS += NPS;
			scan->prev->nim = Images->length;
			if (scan->prev->prev->nim == scan->prev->nim) {
				scan->prev->err = fabs((scan->prev->Mag - scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin) / 4);
			}
			else {
				scan->prev->err = fabs((scan->prev->bin * scan->prev->bin * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->prev->f - scan->prev->f) / 4);
			}
			if (scan->nim == scan->prev->nim) {
				scan->err = fabs((scan->Mag - scan->prev->Mag) * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin) / 4);
			}
			else {
				scan->err = fabs((scan->bin * scan->bin * scan->Mag - scan->prev->bin * scan->prev->bin * scan->prev->Mag) * (scan->prev->f - scan->f) / 4);
			}
			rb = (scan->Mag + scan->prev->prev->Mag - 2 * scan->prev->Mag);
			scan->prev->err += fabs(rb * (scan->prev->prev->f - scan->prev->f) * (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin));
			scan->err += fabs(rb * (scan->prev->f - scan->f) * (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin));
#ifdef _PRINT_ERRORS_DARK
			printf("\n%d", Images->length);
#endif
			delete Images;

			Mag += (scan->bin * scan->bin * scan->Mag - cb * cb * scan->prev->Mag) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			Mag += (cb * cb * scan->prev->Mag - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->Mag) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
			currerr += scan->err + scan->prev->err;
			if (astrometry) {
				LDastrox1 += (scan->bin * scan->bin * scan->LDastrox1 - cb * cb * scan->prev->LDastrox1) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
				LDastrox1 += (cb * cb * scan->prev->LDastrox1 - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->LDastrox1) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
				LDastrox2 += (scan->bin * scan->bin * scan->LDastrox2 - cb * cb * scan->prev->LDastrox2) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
				LDastrox2 += (cb * cb * scan->prev->LDastrox2 - scan->prev->prev->bin * scan->prev->prev->bin * scan->prev->prev->LDastrox2) * (scan->prev->cum - scan->prev->prev->cum) / (scan->prev->bin * scan->prev->bin - scan->prev->prev->bin * scan->prev->prev->bin);
			}


			if (fabs(Magold - Mag) * 2 < Tolv) {
				flag++;
			}
			else {
				flag = 0;
				nannold = nannuli;
			}

		}

		if (multidark) {
			annlist = first;
		}
		else {
			while (first) {
				scan = first->next;
				delete first;
				first = scan;
			}
		}

		Tolv /= 10;
		c++;
	}
	NPS = totNPS;
	therr = currerr;
	if (astrometry) {
		LDastrox1 /= Mag;
		LDastrox2 /= Mag;
		astrox1 = LDastrox1;
		astrox2 = LDastrox2;
	}
	return Mag;
}

void VBMicrolensing::BinaryMagMultiDark(double a, double q, double y1, double y2, double RSv, double* a1_list, int nfil, double* mag_list, double Tol) {
	annulus* scan;
	int imax = 0;
	double Mag, r2, cr2, scr2, a1;

	multidark = true;

	for (int i = 1; i < nfil; i++) {
		if (a1_list[i] > a1_list[imax]) imax = i;
	}
	a1 = a1_list[imax];
	mag_list[imax] = BinaryMagDark(a, q, y1, y2, RSv, Tol);

	for (int i = 0; i < nfil; i++) {
		if (i != imax) {
			Mag = 0;
			a1 = a1_list[i];
			for (scan = annlist->next; scan; scan = scan->next) {
				r2 = scan->bin * scan->bin;
				cr2 = 1 - r2;
				scr2 = sqrt(cr2);
				scan->cum = (3 * r2 * (1 - a1) - 2 * a1 * (scr2 * cr2 - 1)) / (3 - a1);
				Mag += (scan->bin * scan->bin * scan->Mag - scan->prev->bin * scan->prev->bin * scan->prev->Mag) * (scan->cum - scan->prev->cum) / (scan->bin * scan->bin - scan->prev->bin * scan->prev->bin);
			}
			mag_list[i] = Mag;
		}
	}

	while (annlist) {
		scan = annlist->next;
		delete annlist;
		annlist = scan;
	}

	multidark = false;
}


#define _Jacobians1 \
	z=zr[i];\
	dza=z-coefs[20];\
	za2 = dza*dza;\
	zb2=z*z;\
	J1= coefs[21]/za2+coefs[22]/zb2;\
	J1c=conj(J1);\
	dJ=1-J1*J1c;\
	J2=-2.*(coefs[21]/(za2*dza)+coefs[22]/(zb2*z));

#define _Jacobians2\
	dy = complex(-sin(theta->th), cos(theta->th))*coefs[23];\
	dz = (dy - J1c*conj(dy)) / dJ.re;\
	Prov->last->x1 -= coefs[11].re;\
	Prov->last->dJ = dJ.re;\
	Prov->last->d = dz;\
	Prov->last->J2 = J2;\
	Prov->last->ds = (imag(dy*dz*dz*J2) + coefs[23].re*coefs[23].re) / dJ.re;

#define _Jacobians3\
	Prov->last->dJ = dJ.re;\
	J3=6.*(coefs[21]/(za2*za2)+coefs[22]/(zb2*zb2));\
	dJ2=dJ.re*dJ.re;\
	za2=J1c*J1c;\
	J3=J3*za2;\
	ob2=(J2.re*J2.re+J2.im*J2.im)*(6-6*dJ.re+dJ2);\
	J2=J2*J2*za2*J1c;\
	cq= 0.5*(fabs(ob2-6.*J2.re-2.*J3.re*dJ.re)+3*fabs(J2.im))/fabs(dJ.re*dJ2*dJ2); 

#define _Jacobians4\
	zaltc=yc+coefs[21]/dza+coefs[22]/z;\
	za2=(zaltc-coefs[20]);\
	Jaltc=coefs[21]/(za2*za2)+coefs[22]/(zaltc*zaltc);\
	Jalt = conj(Jaltc);\
	JJalt2=(1-J1c*Jalt);\
	J3=J2*J1c*JJalt2;\
	J3=(J3-conj(J3)*Jalt)/(JJalt2*JJalt2*dJ.re);\
	cq=(dJ.re<-100)? 0.0 : (J3.re*J3.re+J3.im*J3.im);


_curve* VBMicrolensing::NewImages(complex yi, complex * coefs, _theta * theta) {
	static complex  y, yc, z, zc, J1, J1c, dy, dz, dJ, J2, J3, dza, za2, zb2, zaltc, Jalt, Jaltc, JJalt2;
	static complex zr[5] = { 0.,0.,0.,0.,0. };
	static double dlmin = 1.0e-4, dlmax = 1.0e-3, good[5], dJ2, ob2, cq;
	static int worst1, worst2, worst3, bad, f1, checkJac;
	static double av = 0.0, m1v = 0.0, disim, disisso;
	static _curve* Prov;
	static _point* scan, * prin, * fifth, * left, * right, * center;

#ifdef _PRINT_TIMES
	static double tim0, tim1;
#endif

	y = yi + coefs[11];
	yc = conj(y);

	// coefs[6]=a*a; coefs[7]=a*a*a; coefs[8]=m2*m2; coefs[9]=a*a*m2*m2; coefs[10]=a*m2; coefs[11]=a*m1; coefs[20]=a; coefs[21]=m1; coefs[22]=m2;

	coefs[0] = coefs[9] * y;
	coefs[1] = coefs[10] * (coefs[20] * (coefs[21] + y * (2 * yc - coefs[20])) - 2 * y);
	coefs[2] = y * (1 - coefs[7] * yc) - coefs[20] * (coefs[21] + 2 * y * yc * (1 + coefs[22])) + coefs[6] * (yc * (coefs[21] - coefs[22]) + y * (1 + coefs[22] + yc * yc));
	coefs[3] = 2 * y * yc + coefs[7] * yc + coefs[6] * (yc * (2 * y - yc) - coefs[21]) - coefs[20] * (y + 2 * yc * (yc * y - coefs[22]));
	coefs[4] = yc * (2 * coefs[20] + y);
	coefs[4] = yc * (coefs[4] - 1) - coefs[20] * (coefs[4] - coefs[21]);
	coefs[5] = yc * (coefs[20] - yc);

	bad = 1;
	disim = -1.;
	f1 = 0;

#ifdef _PRINT_TIMES
	tim0 = Environment::TickCount;
#endif
	cmplx_roots_gen(zr, coefs, 5, true, true);

#ifdef _PRINT_TIMES
	tim1 = Environment::TickCount;
	inc += tim1 - tim0;
#endif
	// apply lens equation to check if it is really solved
	for (int i = 0; i < 5; i++) {
		z = zr[i];
		zc = conj(z);
		good[i] = abs(_LL); // Lens equation check
		switch (i) {
		case 0:
			worst1 = i;
			break;
		case 1:
			if (good[i] > good[worst1]) {
				worst2 = worst1;
				worst1 = i;
			}
			else worst2 = i;
			break;
		case 2:
			if (good[i] > good[worst1]) {
				worst3 = worst2;
				worst2 = worst1;
				worst1 = i;
			}
			else if (good[i] > good[worst2]) {
				worst3 = worst2;
				worst2 = i;
			}
			else worst3 = i;
			break;
		default:
			if (good[i] > good[worst1]) {
				worst3 = worst2;
				worst2 = worst1;
				worst1 = i;
			}
			else if (good[i] > good[worst2]) {
				worst3 = worst2;
				worst2 = i;
			}
			else if (good[i] > good[worst3]) {
				worst3 = i;
			}
		}
	}
	Prov = new _curve;
	checkJac = 0;
	//	if (!((good[worst3] < dlmin) && ((good[worst1] < dlmin) || (good[worst2] > dlmax)))) {  // old check for unacceptable roots

	// 3 good roots
	if (good[worst2] * dlmin > good[worst3] + 1.e-12) {
		for (int i = 0; i < 5; i++) {
			if ((i != worst1) && (i != worst2)) {
				//if((i==worst3)&&(good[i]>dlmax)&&(good[worst2]>1.e2*good[worst3])){
				//	zr[i]=(coefs[21].re<coefs[22].re)? 0.5*coefs[20]+coefs[21]/(0.5*coefs[20]-yc-coefs[22]/coefs[20]) : -0.5*coefs[20]+coefs[22]/(-0.5*coefs[20]-yc+coefs[21]/coefs[20]);
				//}
				Prov->append(zr[i].re, zr[i].im);

				_Jacobians1
					if (theta->th >= 0) {
						_Jacobians2
					}
					else {
						_Jacobians3
							corrquad += cq;
					}
				checkJac += (fabs(Prov->last->dJ) > 1.e-7) ? _sign(Prov->last->dJ) : 10;
				Prov->last->theta = theta;

			}
		}
		if (theta->th < 0) {
			dz = zr[worst2] - zr[worst1];

			int i = worst1;
			_Jacobians1
				_Jacobians4
				corrquad2 = cq;

			i = worst2;
			_Jacobians1
				_Jacobians4
				if (cq > corrquad2) corrquad2 = cq;
			//_Jacobians3
			//corrquad2 +=  1/cq;

		}
		else {
			theta->errworst = abs(zr[worst1] - zr[worst2]);
		}

	}
	else {
		if (good[worst2] * dlmax > good[worst3] + 1.e-12 && theta->th >= 0) { // Dubious cases. Better exclude them
			return Prov;
		}
		else {		// 5 good roots
			f1 = 0;
			for (int i = 0; i < 5; i++) {
				Prov->append(zr[i].re, zr[i].im);

				_Jacobians1
					if (theta->th >= 0) {
						_Jacobians2
					}
					else {
						_Jacobians3
							corrquad += cq;
					}
				checkJac += (fabs(Prov->last->dJ) > 1.e-7) ? _sign(Prov->last->dJ) : 10;
				Prov->last->theta = theta;

				if (fabs(dJ.re) < 1.e-5) f1 = 1;
			}
			theta->errworst = -1.e100;
			// check Jacobians in ambiguous cases
			if (f1) {
				left = right = center = fifth = 0;
				dJ.re = 0;
				for (scan = Prov->first; scan; scan = scan->next) {
					if (_sign(scan->x2) == _sign(y.im)) {
						prin = scan;
					}
					else {
						dz.re = fabs(scan->dJ);
						if (dz.re > dJ.re) {
							fifth = scan;
							dJ.re = dz.re;
						}
					}
				}
				for (scan = Prov->first; scan; scan = scan->next) {
					if ((scan != prin) && (scan != fifth)) {
						if (left) {
							if (scan->x1 < left->x1) {
								if (left != right) {
									center = left;
								}
								left = scan;
							}
							else {
								if (scan->x1 > right->x1) {
									if (left != right) {
										center = right;
									}
									right = scan;
								}
								else {
									center = scan;
								}
							}
						}
						else {
							left = right = center = scan;
						}
					}
				}
				if (left->dJ > 0) left->dJ = -left->dJ;
				if (center->dJ < 0) center->dJ = -center->dJ;
				if (right->dJ > 0) right->dJ = -right->dJ;
			}
		}
	}
	if (checkJac != -1) {
		//		printf("\ncheckJac!");
		if (theta->th < 0) {
			dJ = 0;
			for (scan = Prov->first; scan; scan = scan->next) {
				dJ = dJ + 1 / fabs(scan->dJ);
			}
			if (fabs(dJ.re - 1) < Tol) {
				checkJac = -1;
				corrquad = 0;
			}
		}
		if (checkJac != -1) {
			_point* scan2;
			for (scan = Prov->first; scan; scan = scan2) {
				scan2 = scan->next;
				Prov->drop(scan);
				delete scan;
			}
		}
	}
	return Prov;
}


void VBMicrolensing::OrderImages(_sols * Sols, _curve * Newpts) {
	static double A[5][5];
	static _curve* cprec[5];
	static _curve* cpres[5];
	static _curve* cfoll[5];
	static _point* scan, * scan2, * scan3, * isso[2];
	static _curve* scurve, * scurve2;

	_theta* theta;
	static double th, mi, cmp, cmp2, cmp_2, dx2, avgx2, avgx1, avg2x1, pref, d2x2, dx1, d2x1, avgwedgex1, avgwedgex2, parab1, parab2;

	int nprec = 0, npres, nfoll = 0, issoc[2], ij;

	theta = Newpts->first->theta;
	th = theta->th;
	theta->Mag = theta->prev->Mag = theta->maxerr = theta->prev->maxerr = 0;
	theta->astrox1 = theta->prev->astrox1 = theta->astrox2 = theta->prev->astrox2 = 0;
	if (Newpts->length == 3) {
		mi = theta->next->errworst - theta->errworst;
		if ((mi > theta->errworst) && (theta->prev->errworst > 0.)) {
			theta->prev->maxerr = mi * mi;
		}
		mi = theta->prev->errworst - theta->errworst;
		if ((mi > theta->errworst) && (theta->next->errworst > 0.)) {
			theta->maxerr = mi * mi;
		}
	}

	// Per ciascuna immagine troviamo il punto in cui inserire i nuovi punti 
	scurve = Sols->first;
	for (int i = 0; i < Sols->length; i++) {
		if (th < scurve->first->theta->th) {
			if (th > scurve->first->theta->prev->prev->th) {
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
	npres = Newpts->length;

	//if((theta->th>4.7116917419)&&(theta->th<4.711691759)){
	//	theta->th=theta->th;
	//}
	//if((theta->prev->th>4.7116917419)&&(theta->prev->th<4.711691759)){
	//	theta->th=theta->th;
	//}


	// Caso di creazione nuove immagini// 

	if (nprec < npres) {
		mi = 1.e100;
		scan = Newpts->first;
		for (int i = 0; i < Newpts->length - 1; i++) {
			scan2 = scan->next;
			for (int j = i + 1; j < Newpts->length; j++) {
				cmp = (*scan2) - (*scan);
				if (cmp < mi) {
					mi = cmp;
					isso[0] = scan;
					isso[1] = scan2;
				}
				scan2 = scan2->next;
			}
			scan = scan->next;
		}
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
		cpres[3] = scurve;
		cpres[4] = scurve2;
		scan = isso[0];
		scan2 = isso[1];

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		cmp = sqrt(mi / cmp2);   // Delta theta tilde

		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp * 0.04166666667;
		parab1 = -(-scan->ds + scan2->ds) * mi;
		parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scurve->parabstart = 0.5 * (parab1 + parab2);
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: created image:
		if (astrometry) {
			avgwedgex1 = -(-scan->x1 * scan->ds + scan2->x1 * scan2->ds) * mi;
			avgwedgex2 = -(-scan->x2 * scan->ds + scan2->x2 * scan2->ds) * mi;
			dx2 = -(-scan->d.im + scan2->d.im);
			d2x2 = dx2 * dx2;
			dx1 = -(-scan->d.re + scan2->d.re);
			d2x1 = dx1 * dx1;
			scurve->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
			scurve->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;
		}
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scurve->parabstart, (scan->ds + scan2->ds) * mi / 2, fabs(scurve->parabstart) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2 * cmp * cmp2) * cmp);
#endif

		mi = fabs((parab1 - parab2) * 0.5) + fabs(scurve->parabstart) * (cmp_2 * 0.1) + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2 * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scurve->parabstart = 0.;
#endif

#ifdef _selectimage
		if (_selectionimage)
#endif
			pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
		theta->prev->Mag -= ((scan->dJ > 0) ? -1 : 1) * (pref + scurve->parabstart);
		theta->prev->maxerr += mi;

#ifdef _ERRORS_ANALYTIC
		char filnam[32];
		sprintf(filnam, "%02dprev.txt", NPS);
		FILE* f = fopen(filnam, "a+");
		fprintf(f, "%.15le %.15le %.15le\n", -((scan->dJ > 0) ? -1 : 1) * (pref), -((scan->dJ > 0) ? -1 : 1) * (scurve->parabstart), mi);
		fclose(f);
#endif
		scurve2->parabstart = -scurve->parabstart;

		if (astrometry) {
			dx2 = scan2->x2 - scan->x2;
			avgx1 = scan->x1 + scan2->x1;
			avg2x1 = avgx1 * avgx1;
			avgx2 = scan->x2 + scan2->x2;
			theta->prev->astrox1 += ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + scurve->parabastrox1);
			theta->prev->astrox2 -= ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + scurve->parabastrox2);
			scurve2->parabastrox2 = -scurve->parabastrox2;
			scurve2->parabastrox1 = -scurve->parabastrox1;
		}



	}

	// Caso di distruzione immagini// 
	if (nprec > npres) {
		mi = 1.e100;
		for (int i = 0; i < nprec - 1; i++) {
			for (int j = i + 1; j < nprec; j++) {
				cmp = *(cprec[i]->last) - *(cprec[j]->last);
				if (cmp < mi) {
					mi = cmp;
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}
		cprec[issoc[0]]->partneratend = cprec[issoc[1]];
		cprec[issoc[1]]->partneratend = cprec[issoc[0]];

		scan = cprec[issoc[0]]->last;
		scan2 = cprec[issoc[1]]->last;

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		cmp = sqrt(mi / cmp2);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp * 0.04166666666667;
		parab1 = -(scan->ds - scan2->ds) * mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: destructed image:
		if (astrometry) {
			avgwedgex1 = -(scan->x1 * scan->ds - scan2->x1 * scan2->ds) * mi;
			avgwedgex2 = -(scan->x2 * scan->ds - scan2->x2 * scan2->ds) * mi;
			dx2 = -(scan->d.im - scan2->d.im);
			d2x2 = dx2 * dx2;
			dx1 = -(scan->d.re - scan2->d.re);
			d2x1 = dx1 * dx1;
			scan->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
			scan->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;

		}
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds) * mi / 2, fabs(scan->parab) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2 * cmp * cmp2) * cmp);
#endif

		mi = fabs((parab1 - parab2) * 0.5) + fabs(scan->parab) * (cmp * cmp * 0.1) + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
		theta->prev->Mag += ((scan->dJ > 0) ? -1 : 1) * (pref + scan->parab);
#ifdef _ERRORS_ANALYTIC
		char filnam[32];
		sprintf(filnam, "%02dprev.txt", NPS);
		FILE* f = fopen(filnam, "a+");
		fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1) * (pref), ((scan->dJ > 0) ? -1 : 1) * (scan->parab), mi);
		fclose(f);
#endif
		if (astrometry) {
			dx2 = scan2->x2 - scan->x2;
			avgx1 = scan->x1 + scan2->x1;
			avg2x1 = avgx1 * avgx1;
			avgx2 = scan->x2 + scan2->x2;
			theta->prev->astrox1 -= ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + scan->parabastrox1);
			theta->prev->astrox2 += ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + scan->parabastrox2);
		}
		theta->prev->maxerr += mi;
		scan2->parab = -scan->parab;
		if (astrometry) {
			scan2->parabastrox2 = -scan->parabastrox2;
			scan2->parabastrox1 = -scan->parabastrox1;
		}


		nprec -= 2;
		ij = 0;
		for (int i = 0; i < nprec; i++) {
			if (i == issoc[0]) ij++;
			if (i == issoc[1] - 1) ij++;
			cprec[i] = cprec[i + ij];
		}
	}

	// Costruzione matrice distanze con immagini precedenti// 
	mi = 1.e100;
	for (int i = 0; i < nprec; i++) {
		cpres[i] = cprec[i];
		scan = Newpts->first;
		for (int j = 0; j < nprec; j++) {
			A[i][j] = (signbit(cprec[i]->last->dJ) == signbit(scan->dJ)) ? *(cprec[i]->last) - *scan : 100;
			if (A[i][j] < mi) {
				mi = A[i][j];
				issoc[0] = i;
				issoc[1] = j;
				isso[1] = scan;
			}
			scan = scan->next;
		}
	}

	//  Associazione con le immagini che precedono// 
	while (nprec) {
		scan = cprec[issoc[0]]->last;
		scan2 = isso[1];

		cmp2 = mi / fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		cmp = (scan->theta->th - scan2->theta->th);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp * 0.0416666666666667;  ////// (1/24 cube(delta Teta))
		parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
		// Nuova correzione parabolica
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: ordinary image:
		if (astrometry) {
			avgwedgex1 = (scan->x1 * scan->ds + scan2->x1 * scan2->ds) * mi;
			avgwedgex2 = (scan->x2 * scan->ds + scan2->x2 * scan2->ds) * mi;
			dx2 = scan->d.im + scan2->d.im;
			d2x2 = dx2 * dx2;
			dx1 = scan->d.re + scan2->d.re;
			d2x1 = dx1 * dx1;
			scan->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
			scan->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;
		}
#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds) * mi / 2, fabs(scan->parab) * (cmp2) / 10, fabs(scan->parab) * (1.5 * fabs(cmp2 / (cmp * cmp) - 1)));
#endif

		mi = fabs((parab1 - parab2) * 0.5) + fabs(scan->parab * (cmp2 * 0.1 + 1.5 * fabs(cmp2 / (cmp_2)-1)));
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
		theta->prev->Mag += ((scan->dJ > 0) ? -1 : 1) * (pref + scan->parab);
#ifdef _ERRORS_ANALYTIC
		char filnam[32];
		sprintf(filnam, "%02dprev.txt", NPS);
		FILE* f = fopen(filnam, "a+");
		fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1) * (pref), ((scan->dJ > 0) ? -1 : 1) * (scan->parab), mi);
		fclose(f);
#endif
		if (astrometry) {
			dx2 = scan2->x2 - scan->x2;
			avgx1 = scan->x1 + scan2->x1;
			avg2x1 = avgx1 * avgx1;
			avgx2 = scan->x2 + scan2->x2;
			theta->prev->astrox1 -= ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + scan->parabastrox1);
			theta->prev->astrox2 += ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + scan->parabastrox2);
		}
		theta->prev->maxerr += mi;



		Newpts->drop(isso[1]);
		cprec[issoc[0]]->append(isso[1]);
		cprec[issoc[0]]->partneratend = 0;

		nprec--;
		for (int i = issoc[0]; i < nprec; i++) {
			cprec[i] = cprec[i + 1];
			for (int j = 0; j < nprec + 1; j++) {
				A[i][j] = A[i + 1][j];
			}
		}
		for (int j = issoc[1]; j < nprec; j++) {
			for (int i = 0; i < nprec; i++) {
				A[i][j] = A[i][j + 1];
			}
		}
		mi = 1.e100;
		for (int i = 0; i < nprec; i++) {
			scan = Newpts->first;
			for (int j = 0; j < nprec; j++) {
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
	delete Newpts;

#ifdef _PRINT_ERRORS
	printf("\nN");
#endif

	// immagini seguenti// 
	if (nfoll) {
		// Caso di creazione nuove immagini

		if (npres < nfoll) {
			mi = 1.e100;
			for (int i = 0; i < nfoll - 1; i++) {
				for (int j = i + 1; j < nfoll; j++) {
					cmp = *(cfoll[i]->first) - *(cfoll[j]->first);
					if (cmp < mi) {
						mi = cmp;
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}
			cfoll[issoc[0]]->partneratstart = cfoll[issoc[1]];
			cfoll[issoc[1]]->partneratstart = cfoll[issoc[0]];

			scan = cfoll[issoc[0]]->first;
			scan2 = cfoll[issoc[1]]->first;

			cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
			cmp = sqrt(mi / cmp2);
			cmp_2 = cmp * cmp;
			mi = cmp_2 * cmp * 0.04166666666666667;
			parab1 = (scan->ds - scan2->ds) * mi;
			parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
			cfoll[issoc[0]]->parabstart = 0.5 * (parab1 + parab2);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: created image:
			if (astrometry) {
				avgwedgex1 = (scan->x1 * scan->ds - scan2->x1 * scan2->ds) * mi;
				avgwedgex2 = (scan->x2 * scan->ds - scan2->x2 * scan2->ds) * mi;
				dx2 = -(-scan->d.im + scan2->d.im);
				d2x2 = dx2 * dx2;
				dx1 = -(-scan->d.re + scan2->d.re);
				d2x1 = dx1 * dx1;
				cfoll[issoc[0]]->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
				cfoll[issoc[0]]->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;
			}
#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, cfoll[issoc[0]]->parabstart, (scan->ds + scan2->ds) * mi / 2, fabs(cfoll[issoc[0]]->parabstart) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2 * cmp * cmp2) * cmp);
#endif
			mi = fabs((parab1 - parab2) * 0.5) + fabs(cfoll[issoc[0]]->parabstart) * (cmp * cmp * 0.1) + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
			mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
			cfoll[issoc[0]]->parabstart = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif
				pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
			theta->Mag -= ((scan->dJ > 0) ? -1 : 1) * (pref + cfoll[issoc[0]]->parabstart);
#ifdef _ERRORS_ANALYTIC
			char filnam[32];
			sprintf(filnam, "%02dfoll.txt", NPS);
			FILE* f = fopen(filnam, "a+");
			fprintf(f, "%.15le %.15le %.15le\n", -((scan->dJ > 0) ? -1 : 1) * (pref), -((scan->dJ > 0) ? -1 : 1) * (cfoll[issoc[0]]->parabstart), mi);
			fclose(f);
#endif
			if (astrometry) {
				dx2 = scan2->x2 - scan->x2;
				avgx1 = scan->x1 + scan2->x1;
				avg2x1 = avgx1 * avgx1;
				avgx2 = scan->x2 + scan2->x2;
				theta->astrox1 += ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + cfoll[issoc[0]]->parabastrox1);
				theta->astrox2 -= ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + cfoll[issoc[0]]->parabastrox2);
			}
			theta->maxerr += mi;

			cfoll[issoc[1]]->parabstart = -cfoll[issoc[0]]->parabstart;
			if (astrometry) {
				cfoll[issoc[1]]->parabastrox2 = -cfoll[issoc[0]]->parabastrox2;
				cfoll[issoc[1]]->parabastrox1 = -cfoll[issoc[0]]->parabastrox1;
			}
			Sols->append(cfoll[issoc[0]]);
			Sols->append(cfoll[issoc[1]]);
			nfoll -= 2;
			ij = 0;
			for (int i = 0; i < nfoll; i++) {
				if (i == issoc[0]) ij++;
				if (i == issoc[1] - 1) ij++;
				cfoll[i] = cfoll[i + ij];
			}
		}


		// Caso di distruzione immagini
		if (npres > nfoll) {
			mi = 1.e100;
			for (int i = 0; i < npres - 1; i++) {
				for (int j = i + 1; j < npres; j++) {
					cmp = *(cpres[i]->last) - *(cpres[j]->last);
					if (cmp < mi) {
						mi = cmp;
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}
			cpres[issoc[0]]->partneratend = cpres[issoc[1]];
			cpres[issoc[1]]->partneratend = cpres[issoc[0]];

			scan = cpres[issoc[0]]->last;
			scan2 = cpres[issoc[1]]->last;

			cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
			cmp = sqrt(mi / cmp2);
			cmp_2 = cmp * cmp;
			mi = cmp_2 * cmp * 0.0416666666667;
			parab1 = -(scan->ds - scan2->ds) * mi;
			parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
			scan->parab = 0.5 * (parab1 + parab2);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: destructed image:
			if (astrometry) {
				avgwedgex1 = -(scan->x1 * scan->ds - scan2->x1 * scan2->ds) * mi;
				avgwedgex2 = -(scan->x2 * scan->ds - scan2->x2 * scan2->ds) * mi;
				dx2 = -(scan->d.im - scan2->d.im);
				d2x2 = dx2 * dx2;
				dx1 = -(scan->d.re - scan2->d.re);
				d2x1 = dx1 * dx1;
				scan->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
				scan->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;
			}

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds) * mi / 2, fabs(scan->parab) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2 * cmp * cmp2) * cmp);
#endif

			mi = fabs((parab1 - parab2) * 0.5) + fabs(scan->parab) * (cmp * cmp * 0.1) + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
			mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
			scan->parab = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif
				pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
			theta->Mag += ((scan->dJ > 0) ? -1 : 1) * (pref + scan->parab);
#ifdef _ERRORS_ANALYTIC
			char filnam[32];
			sprintf(filnam, "%02dfoll.txt", NPS);
			FILE* f = fopen(filnam, "a+");
			fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1) * (pref), ((scan->dJ > 0) ? -1 : 1) * (scan->parab), mi);
			fclose(f);
#endif
			if (astrometry) {
				dx2 = scan2->x2 - scan->x2;
				avgx1 = scan->x1 + scan2->x1;
				avg2x1 = avgx1 * avgx1;
				avgx2 = scan->x2 + scan2->x2;
				theta->astrox1 -= ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + scan->parabastrox1);
				theta->astrox2 += ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + scan->parabastrox2);
			}
			theta->maxerr += mi;
			scan2->parab = -scan->parab;
			if (astrometry) {
				scan2->parabastrox2 = -scan->parabastrox2;
				scan2->parabastrox1 = -scan->parabastrox1;
			}
			npres -= 2;
			ij = 0;
			for (int i = 0; i < npres; i++) {
				if (i == issoc[0]) ij++;
				if (i == issoc[1] - 1) ij++;
				cpres[i] = cpres[i + ij];
			}
		}

		// Costruzione matrice distanze con immagini seguenti

		mi = 1.e100;
		for (int i = 0; i < npres; i++) {
			for (int j = 0; j < npres; j++) {
				A[i][j] = signbit(cpres[i]->last->dJ) == signbit(cfoll[j]->first->dJ) ? *(cpres[i]->last) - *(cfoll[j]->first) : 100;
				if (A[i][j] < mi) {
					mi = A[i][j];
					issoc[0] = i;
					issoc[1] = j;
				}
			}
		}


		//  Associazione con le immagini che seguono// 
		while (npres) {
			scan = cpres[issoc[0]]->last;
			scan2 = cfoll[issoc[1]]->first;
			cmp2 = mi / fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
			cmp = (scan->theta->th - scan2->theta->th);
			cmp_2 = cmp * cmp;
			mi = cmp_2 * cmp * 0.041666666667;
			parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
			// Nuova correzione parabolica
			parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
			scan->parab = 0.5 * (parab1 + parab2);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////astro: ordinary image:
			if (astrometry) {
				avgwedgex1 = (scan->x1 * scan->ds + scan2->x1 * scan2->ds) * mi;
				avgwedgex2 = (scan->x2 * scan->ds + scan2->x2 * scan2->ds) * mi;
				dx2 = scan->d.im + scan2->d.im;
				d2x2 = dx2 * dx2;
				dx1 = scan->d.re + scan2->d.re;
				d2x1 = dx1 * dx1;
				scan->parabastrox1 = -0.125 * d2x1 * dx2 * mi - avgwedgex1;
				scan->parabastrox2 = -0.125 * d2x2 * dx1 * mi + avgwedgex2;
			}

#ifdef _PRINT_ERRORS
			printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds) * mi / 2, fabs(scan->parab) * (cmp2) / 10, fabs(scan->parab) * (1.5 * fabs(cmp2 / (cmp * cmp) - 1)));
#endif

			mi = fabs((parab1 - parab2) * 0.5) + fabs(scan->parab * (cmp2 * 0.1 + 1.5 * fabs(cmp2 / (cmp_2)-1)));
#ifdef _noparab
			mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
			scan->parab = 0.;
#endif
#ifdef _selectimage
			if (_selectionimage)
#endif
				pref = (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5;
			theta->Mag += ((scan->dJ > 0) ? -1 : 1) * (pref + scan->parab);
			if (astrometry) {
				dx2 = scan2->x2 - scan->x2;
				avgx1 = scan->x1 + scan2->x1;
				avg2x1 = avgx1 * avgx1;
				avgx2 = scan->x2 + scan2->x2;
				theta->astrox1 -= ((scan->dJ > 0) ? -1 : 1) * (avg2x1 * dx2 * 0.125 + scan->parabastrox1);
				theta->astrox2 += ((scan->dJ > 0) ? -1 : 1) * (pref * avgx2 * 0.25 + scan->parabastrox2);
			}
			theta->maxerr += mi;
#ifdef _ERRORS_ANALYTIC
			char filnam[32];
			sprintf(filnam, "%02dfoll.txt", NPS);
			FILE* f = fopen(filnam, "a+");
			fprintf(f, "%.15le %.15le %.15le\n", ((scan->dJ > 0) ? -1 : 1) * (pref), ((scan->dJ > 0) ? -1 : 1) * (scan->parab), mi);
			fclose(f);
#endif
			cpres[issoc[0]]->join(cfoll[issoc[1]]);

			npres--;
			for (int i = issoc[0]; i < npres; i++) {
				cpres[i] = cpres[i + 1];
				for (int j = 0; j < npres + 1; j++) {
					A[i][j] = A[i + 1][j];
				}
			}
			for (int j = issoc[1]; j < npres; j++) {
				cfoll[j] = cfoll[j + 1];
				for (int i = 0; i < npres; i++) {
					A[i][j] = A[i][j + 1];
				}
			}
			mi = 1.e100;
			for (int i = 0; i < npres; i++) {
				for (int j = 0; j < npres; j++) {
					if (A[i][j] < mi) {
						mi = A[i][j];
						issoc[0] = i;
						issoc[1] = j;
					}
				}
			}
		}
	}

}


#pragma endregion

#pragma region multi-mag


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

void VBMicrolensing::SetLensGeometry(int nn, double* q, complex * s) {
	switch (SelectedMethod)
	{
	case Method::Singlepoly:
		SetLensGeometry_spnp(nn, q, s);
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

void VBMicrolensing::SetLensGeometry_spnp(int nn, double* q, complex * s) {
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

void VBMicrolensing::SetLensGeometry_multipoly(int nn, double* q, complex * s) {
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

double VBMicrolensing::MultiMag0(complex yi, _sols * *Images) {
	static double Mag = -1.0;
	_theta* stheta;
	_curve* Prov = 0, * Prov2;
	_point* scan1, * scan2;

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
	mag = MultiMag0(y, &images);
	delete images;
	return mag;
}

double VBMicrolensing::MultiMag(complex yi, double RSv, double Tol, _sols * *Images) {
	static complex y0;
	static double Mag = -1.0, th, thoff = 0.01020304, thoff2 = 0.7956012033974483; //0.01020304
	static double errimage, maxerr, currerr, Magold, rhorad2, th2;
	static int NPSmax, flag, NPSold, isquare, flagfinal;
	static _thetas* Thetas;
	static _theta* stheta, * itheta, * jtheta;
	static _curve* Prov, * Prov2;
	static _point* scan1, * scan2;
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
			errimage = Tol * M_PI * RSv * RSv;
			NPSmax = 2000;
		}

		// Calculation of the images

		(*Images) = new _sols;
		Thetas = new _thetas;
		th = thoff;
		stheta = Thetas->insert(th);
		stheta->maxerr = 1.e100;
		y = y0 + complex(RSv * cos(thoff), RSv * sin(thoff)); // first image


#ifdef _PRINT_TIMES
		tim0 = Environment::TickCount;
#endif

		EXECUTE_METHOD(SelectedMethod, stheta)

#ifdef _PRINT_TIMES
			tim1 = Environment::TickCount;
		GM += tim1 - tim0;
#endif
		stheta = Thetas->insert(2.0 * M_PI + thoff);
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
			y = y0 + complex(RSv * cos(th), RSv * sin(th));

			EXECUTE_METHOD(SelectedMethod, stheta)

				OrderMultipleImages((*Images), Prov);
		}
		NPS = 4;

		// Let us check the number if images in a square around the source
		if (squarecheck) {
			itheta = Thetas->first;
			stheta = new _theta(-1);
			while (itheta->next) {
				th = (itheta->th + itheta->next->th) * 0.5;
				rhorad2 = RSv / cos((itheta->next->th - itheta->th) * 0.5);
				y = y0 + complex(rhorad2 * cos(th), rhorad2 * sin(th));

				EXECUTE_METHOD(SelectedMethod, stheta)

					delete Prov;
				if (itheta->imlength != stheta->imlength && itheta->next->imlength != stheta->imlength) {
					jtheta = Thetas->insert(th);
					y = y0 + complex(RSv * cos(th), RSv * sin(th));

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
				if (stheta->next->th * 0.99999 - stheta->th > maxerr) {
					maxerr = stheta->next->th - stheta->th;
#endif
					itheta = stheta;
				}
				}
			stheta = stheta->next;
			}
		th = (itheta->th + itheta->next->th) * 0.5;

		// Main cycle: sampling continues until total error is below Tol.
		flag = 0;
		Magold = -1.;
		NPSold = NPS + 1;

		while (((currerr > errimage) && (currerr > RelTol * Mag) && (NPS < NPSmax) && (flag < NPSold))) {
			stheta = Thetas->insert(th);
			y = y0 + complex(RSv * cos(th), RSv * sin(th));
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
					if (stheta->next->th * 0.99999 - stheta->th > maxerr) {
						maxerr = stheta->next->th - stheta->th;
#endif
						itheta = stheta;
					}
					}
				stheta = stheta->next;
				}
			th = (itheta->th + itheta->next->th) * 0.5;
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
			printf("\nNPS= %d nim=%d Mag = %lf maxerr= %lg currerr =%lg th = %lf", NPS, lim, Mag / (M_PI * RSv * RSv), maxerr / (M_PI * RSv * RSv), currerr / (M_PI * RSv * RSv), th);
#endif

			}
		Mag /= (M_PI * RSv * RSv);
		therr = currerr / (M_PI * RSv * RSv);

		delete Thetas;

		return Mag;

		}
	catch (...) {
		FILE* f = fopen("Geom.txt", "w");
		fprintf(f, "\n%d\n", n);
		for (int i = 0; i < n; i++) {
			fprintf(f, "%.16lf %.16lf %.16lf\n", m[i], a[i].re, a[i].im);
		}
		fprintf(f, "\n%.16lf %.16lf %.16lf\n", y.re, y.im, rho);
		for (int i = 0; i < ngood; i++) {
			fprintf(f, "%.16lf %.16lf %.16lg %.16lf %.16lg\n", zr[i].re, zr[i].im, errs[i], Jacs[i], good[i]);
		}
		fclose(f);
		return -1;
	}



		}

double VBMicrolensing::MultiMag(complex y, double RSv, double Tol) {
	static _sols* images;
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
	static complex fac, fac2, z, S2;
	static double Jac;

	// Main image, should be positive and close to source
	init[n] = y;
	for (int i = 0; i < n; i++) {
		fac = y - a[i];
		prodevs[i] = m[i] / abs2(fac);
		devs[i] = (prodevs[i] > 1.e-2) ? (sqrt(0.25 + prodevs[i]) - 0.5) : prodevs[i];
		devs[i] = devs[i] * fac;
		init[n] = init[n] + devs[i];
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
				devs[i] = devs[i] * 0.5;
				z = z + devs[i];
				_Jac
			}
		init[i] = z;
	}
	// Central images are calculated in SetLensGeometry
}

int VBMicrolensing::froot(complex zi) {
	static complex z, zc, S1, S2, S2v, S3, zo, zo2, epso, epsbase, epsn, epsl, gradL, zl, fac, fac2, dz, dzo, TJold, TJnew, Lv, den;
	static int iter3, iter4, ipseudo, flagmain, flag;
	static double Lnew, Lold, fad, Jac, Jacold, prefac, epsbo;

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
	epsbase = epsn = ((Jac > -15) ? fad * (sqrt(sqrt(fad)) - 1) : fad) / (conj(S2) * S3) * 0.5;
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

		epsbase = epsn = (conj(Lv) - Lv * conj(S2)) / Jac;
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
			while (abs2(zo2 - zo - epsn) < 1.69 * fad && fad > 1.e-28) {
				fad = fad;
				epsn = epsn * pseudorandom[ipseudo]; // In case of doubleoscillation, shorten step by pseudorandom
				fad = abs2(epsn);
				ipseudo++;
				if (ipseudo > 11) ipseudo = 0;
				flagmain = 4;
			}

			// New step cannot be twice longer than previous step (avoid explosion)
			fad = sqrt(abs2(zo2 - zo) / fad);
			fad *= (prefac < 1) ? (1 + prefac) : 2;
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
					fad *= (1 - Jacold / Jac) * (sqrt(sqrt(fad)) - 1);
					dz = dzo = fad / (conj(S2) * S3);
					// If correction is below machine precision, go back a little bit 
					if (abs2(dz) < 1.e-28) {
						dz = -epsn * 0.15;
						flagmain = 3;
					}
					// If correction is higher than step shorten correction and increase iter2. We want to break cycles that do not converge
					fad = abs(epsn / dz);
					if (abs2(epsn + dz) < 0.25 * abs2(epsn) || fad < 2) {
						dz = complex(0.33, 0.01) * dz * fad;
						flagmain = 3;
						if (abs2(epsn + dz) < abs2(zo2 - zo)) iter2++;
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
							if (TJold.re * TJnew.re + TJold.im * TJnew.im < 0) {
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
									gradL = -(conj(Lv) + conj(S2v) * Lv);
									den = (TJnew.re * gradL.re + TJnew.im * gradL.im);
									epsl = -0.25 * Lnew / den * TJnew;
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
		if (abs2(z - zo2) < 0.25 * abs2(z - zo)) {
			flagmain = 3;
			iter2++;
		}

		if (flagmain == 0) iter2 = 0;

		iter++;
	}
	newtonstep += iter;
	err = 3.163e-15 / Jac;
	err *= err;
	err += abs2(epsbase);
	zf = z;
	Jacf = Jac;
	S2f = S2;
	L0f = Lnew;
	return iter2;
}

bool VBMicrolensing::checkroot(_theta * theta) {
	static double mn, fac;
	static int imn;
	static complex S3, z, S4;
	static double fad;
	if ((iter2 < 9 && iter < maxiter) || L0f < 1.e-29) {
		mn = 1.e100;
		imn = 0;
		for (int i = 0; i < ngood; i++) {
			if (signbit(Jacf) == signbit(Jacs[i])) {
				fac = abs2(zr[i] - zf) / (errs[i] + err + 4.e-20);
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
			grads[ngood] = fad * (sqrt(sqrt(fad)) - 1) / (conj(S2f) * S3);

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


_curve* VBMicrolensing::NewImages(_theta * theta) {
	static _curve* Prov;
	static int nminus, nplus;
	static complex z, zc, dy, dz, J2, J3, Jalt, JJalt2, Jaltc, J1c2;
	static complex S2, S2c, S3, S3c, vec, newseed0;
	static double ob2, dJ2, cq, Jac, imul, phi;
	static complex newseedtrial[6] = { complex(1,0.5),complex(1,-0.5), 0.5, 0.75, 2., 4., };
	static int imass, iphi, nsafe;

	yc = conj(y);
	initroot();

	ngood = ngoodold = 0;
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
			newseeds[lennewseeds] = zr[i] + grads[i] * vec;
			z = newseeds[lennewseeds];
			_Jac
				iter = 0;
			newseed0 = newseeds[lennewseeds];
			while (signbit(Jac) == signbit(Jacs[i]) && iter < 6) {
				newseeds[lennewseeds] = zr[i] + (newseed0 - zr[i]) * newseedtrial[iter];
				z = newseeds[lennewseeds];
				_Jac
					iter++;
			}
			if (iter < 6) lennewseeds++;
			newseeds[lennewseeds] = zr[i] + grads[i] * conj(vec);
			z = newseeds[lennewseeds];
			_Jac
				iter = 0;
			newseed0 = newseeds[lennewseeds];
			while (signbit(Jac) == signbit(Jacs[i]) && iter < 6) {
				newseeds[lennewseeds] = zr[i] + (newseed0 - zr[i]) * newseedtrial[iter];
				z = newseeds[lennewseeds];
				_Jac
					iter++;
			}
			if (iter < 6) lennewseeds++;
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
	iphi = imass = nsafe = 0;

	while (nminus - nplus + 1 != n) {

		phi = iphi * 2.61799; // 5*M_PI/6.
		z = imul * sqrt(m[imass]) * complex(cos(phi), sin(phi)) + a[imass];
		froot(z);
		checkroot(theta);
		if (ngood > ngoodold) {
			if (Jacs[ngoodold] > 0) {
				nplus++;
			}
			else {
				nminus++;
			}
			ngoodold = ngood;
		}
		imass++;
		if (imass == n) {
			imass = 0;
			iphi++;
			if (iphi == 12) {
				iphi = 0;
				if (imul > 1) {
					imul *= 1.1;
				}
				else {
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
		for (int j = 0; j < n; j++) {
			if (i != j)
				dev2 = dev2 + m[j] / (a[i] - a[j]);
		}
		dev = (sqrt(0.25 + m[i] / abs2(y - a[i])) - 0.5) * (y - a[i]);
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
		dev2 = dev / abs(dev) * J1[i] / 1.;
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

	while (good[i] > dlmax /*&& good[i]<LLold && signbit(Jold) == signbit(Jac[i])*/ && iter < 0 && iter2 < 4) {
		LLold = good[i];
		Jold = Jacs[i];
		zo = z;
		lambda = 1.0;
		//for (int j = nroots-1; j > i; j--) {
		//	if(good[j]<dlmax)
		//		lambda = lambda + conj(LL) /(z - zr[j]);
		//}
		delta = (conj(LL * lambda) - LL * J1c[i]) / (Jacs[i] - 1 + abs2(lambda));
		iter2 = 0;
		deltafac = 1.;
		while (iter2 < 4 && (good[i] + dlmax > LLold /*|| signbit(Jold) != signbit(Jac[i])*/)) {
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
_curve* VBMicrolensing::NewImagespoly(_theta * theta) {
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
		if (ngood < n2 + 1)	isgood = good[worst[ngood]];
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

_curve* VBMicrolensing::NewImagesmultipoly(_theta * theta) {
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

void VBMicrolensing::OrderMultipleImages(_sols * Sols, _curve * Newpts) {
	static _point* scan, * scan2, * scan3, * isso[2];
	static _curve* scurve, * scurve2;

	static _theta* theta;
	static double th, mi, cmp, cmp2, cmp_2, er3, parab1, parab2;
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
			if (th > scurve->first->theta->prev->prev->th) { // sembra strano ma h giusto prev->prev->th (l'inserimento di theta c'h gi` stato)
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
		cmp2 = mi / fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		er3 = mi * mi * mi * samplingfactor;
		cmp = (scan->theta->th - scan2->theta->th);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp * 0.0416666666666667;
		parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
		// Nuova correzione parabolica
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds) * mi / 2, fabs(scan->parab) * (cmp2) / 10, fabs(scan->parab) * (1.5 * fabs(cmp2 / (cmp * cmp) - 1)));
#endif

		mi = fabs((scan->ds - scan2->ds) * mi * 0.5) + fabs(scan->parab * 1.5 * fabs(cmp2 / (cmp_2)-1)) + er3;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = ((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + scan->parab);
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
	while (npres > 1 && mi < 1.e99) {
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

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		cmp_2 = mi / cmp2;
		er3 = sqrt(mi * mi * mi) * samplingfactor;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp * 0.04166666667;
		parab1 = -(-scan->ds + scan2->ds) * mi;
		parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scurve->parabstart = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scurve->parabstart, (scan->ds + scan2->ds) * mi / 2, fabs(scurve->parabstart) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2 * cmp * cmp2) * cmp);
#endif

		mi = fabs((scan->ds + scan2->ds) * mi * 0.5) + er3 + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2. * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scurve->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scurve->Magstart = -(((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + scurve->parabstart));
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
			for (int j = (ij > 0) ? i + 1 : issoc[1]; j < npres + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
		}
		ij = 0;
		for (int i = 0; i < npres; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij > 0) ? i + 1 : issoc[0]; j < npres; j++) {
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
	while (npres > 0) {
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
		mi = fabs(scurve->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
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
	while (nprec > 1 && mi < 1.e99) {

		cprec[issoc[0]]->partneratend = cprec[issoc[1]];
		cprec[issoc[1]]->partneratend = cprec[issoc[0]];

		scan = cprec[issoc[0]]->last;
		scan2 = cprec[issoc[1]]->last;

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		er3 = sqrt(mi * mi * mi) * samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp * 0.04166666666667;
		parab1 = -(scan->ds - scan2->ds) * mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds) * mi / 2, fabs(scan->parab) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2 * cmp * cmp2) * cmp);
#endif

		mi = fabs((scan->ds + scan2->ds) * mi * 0.5) + er3 + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = ((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + scan->parab);
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
			for (int j = (ij > 0) ? i + 1 : issoc[1]; j < nprec + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cprec[i] = cprec[i + ij];
		}
		ij = 0;
		for (int i = 0; i < nprec; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij > 0) ? i + 1 : issoc[0]; j < nprec; j++) {
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
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
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
		cmp2 = mi / fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		er3 = mi * mi * mi * samplingfactor;
		cmp = (scan->theta->th - scan2->theta->th);
		cmp_2 = cmp * cmp;
		mi = cmp_2 * cmp * 0.041666666667;
		parab1 = (scan->ds + scan2->ds) * mi; // Vecchia Correzione parabolica
		// Nuova correzione parabolica
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im - scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re - scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds - scan2->ds) * mi / 2, fabs(scan->parab) * (cmp2) / 10, fabs(scan->parab) * (1.5 * fabs(cmp2 / (cmp * cmp) - 1)));
#endif

		mi = fabs((scan->ds - scan2->ds) * mi * 0.5) + fabs(scan->parab * 1.5 * fabs(cmp2 / (cmp_2)-1)) + er3;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = ((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + scan->parab);
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
	while (nfoll > 1 && mi < 1.e99) {

		cfoll[issoc[0]]->partneratstart = cfoll[issoc[1]];
		cfoll[issoc[1]]->partneratstart = cfoll[issoc[0]];

		scan = cfoll[issoc[0]]->first;
		scan2 = cfoll[issoc[1]]->first;

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		er3 = sqrt(mi * mi * mi) * samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp * 0.04166666666666667;
		parab1 = (scan->ds - scan2->ds) * mi;
		parab2 = -0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		cfoll[issoc[0]]->parabstart = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, cfoll[issoc[0]]->parabstart, (scan->ds + scan2->ds) * mi / 2, fabs(cfoll[issoc[0]]->parabstart) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2 * cmp * cmp2) * cmp);
#endif
		mi = fabs((scan->ds + scan2->ds) * mi * 0.5) + er3 + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) - 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		cfoll[issoc[0]]->parabstart = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			cfoll[issoc[0]]->Magstart = -(((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + cfoll[issoc[0]]->parabstart));
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
			for (int j = (ij > 0) ? i + 1 : issoc[1]; j < nfoll + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cfoll[i] = cfoll[i + ij];
		}
		ij = 0;
		for (int i = 0; i < nfoll; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij > 0) ? i + 1 : issoc[0]; j < nfoll; j++) {
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
		mi = fabs(cfoll[issoc[0]]->parabstart) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
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
	while (npres > 1 && mi < 1.e99) {
		cpres[issoc[0]]->partneratend = cpres[issoc[1]];
		cpres[issoc[1]]->partneratend = cpres[issoc[0]];

		scan = cpres[issoc[0]]->last;
		scan2 = cpres[issoc[1]]->last;

		cmp2 = fabs(scan->d.re * scan2->d.re + scan->d.im * scan2->d.im);
		er3 = sqrt(mi * mi * mi) * samplingfactor;
		cmp_2 = mi / cmp2;
		cmp = sqrt(cmp_2);
		mi = cmp_2 * cmp * 0.0416666666667;
		parab1 = -(scan->ds - scan2->ds) * mi;
		parab2 = 0.0833333333 * ((scan2->x1 - scan->x1) * (scan2->d.im + scan->d.im) - (scan2->x2 - scan->x2) * (scan2->d.re + scan->d.re)) * cmp;
		scan->parab = 0.5 * (parab1 + parab2);

#ifdef _PRINT_ERRORS
		printf("\n%le %le %le %le %le %le %le %le", scan->x1, scan->x2, scan->dJ, (scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2, scan->parab, (scan->ds + scan2->ds) * mi / 2, fabs(scan->parab) * (cmp * cmp) / 10, 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2 * cmp * cmp2) * cmp);
#endif

		mi = fabs((scan->ds + scan2->ds) * mi * 0.5) + er3 + 1.5 * fabs(((scan->d.re - scan2->d.re) * (scan->x1 - scan2->x1) + (scan->d.im - scan2->d.im) * (scan->x2 - scan2->x2)) + 2.0 * cmp * cmp2) * cmp;
#ifdef _noparab
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
		scan->parab = 0.;
#endif
#ifdef _selectimage
		if (_selectionimage)
#endif
			scan->Mag = ((scan->dJ > 0) ? -1 : 1) * ((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) * 0.5 + scan->parab);
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
			for (int j = (ij > 0) ? i + 1 : issoc[1]; j < npres + 1; j++) {
				A[i][j] = A[i + ij][j + 1];
			}
			cpres[i] = cpres[i + ij];
		}
		ij = 0;
		for (int i = 0; i < npres; i++) {
			if (i == issoc[0]) ij++;
			for (int j = (ij > 0) ? i + 1 : issoc[0]; j < npres; j++) {
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
		mi = fabs(scan->parab) * 2 + 0 * fabs((scan->x2 + scan2->x2) * (scan2->x1 - scan->x1) / 2 * cmp * cmp / 6);
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

#pragma endregion

#pragma region lightcurves


//////////////////////////////
//////////////////////////////
////////New (v2) light curve functions
//////////////////////////////
//////////////////////////////

void VBMicrolensing::PSPLLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) * tE_inv;
		u = tn * tn + u0 * u0;

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i] = (u + 2) / sqrt(u * (u + 4));

	}
}

void VBMicrolensing::PSPLLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, pai1 = pr[3], pai2 = pr[4];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u1 = u0 + pai1 * Et[1] - pai2 * Et[0];
		u = tn * tn + u1 * u1;

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = (u + 2) / sqrt(u * (u + 4));
	}

}


void VBMicrolensing::ESPLLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, rho = exp(pr[3]);

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) * tE_inv;
		u = sqrt(tn * tn + u0 * u0);

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i] = ESPLMag2(u, rho);

	}
}

void VBMicrolensing::ESPLLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, rho = exp(pr[3]), pai1 = pr[4], pai2 = pr[5];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u1 = u0 + pai1 * Et[1] - pai2 * Et[0];
		u = sqrt(tn * tn + u1 * u1);

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = ESPLMag2(u, rho);
	}

}


void VBMicrolensing::BinaryLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]);
	double salpha = sin(pr[3]), calpha = cos(pr[3]);

	//	_sols *Images; double Mag; // For debugging

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - pr[6]) * tE_inv;
		y1s[i] = pr[2] * salpha - tn * calpha;
		y2s[i] = -pr[2] * calpha - tn * salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);

		//Mag=BinaryMag(s, q, y1s[i], y2s[i], rho, Tol,&Images); // For debugging
		//delete Images;
		//mags[i] -= Mag;
		//if (fabs(mags[i]) > Tol) {
		//	printf("\n%lf %lf %lf", y1s[i], y2s[i], mags[i]);
		//}
	}
}



void VBMicrolensing::TripleLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), di, mindi;
	double q[3] = { 1, exp(pr[1]),exp(pr[8]) };
	complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);

	s[0] = exp(pr[0]) / (q[0] + q[1]);
	s[1] = s[0] * q[0];
	s[0] = -q[1] * s[0];
	s[2] = exp(pr[7]) * complex(cbeta, sbeta) + s[0];
	//	_sols *Images; double Mag; // For debugging

	SetLensGeometry(3, q, s);

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - pr[6]) * tE_inv;
		y1s[i] = pr[2] * salpha - tn * calpha;
		y2s[i] = -pr[2] * calpha - tn * salpha;
		mindi = 1.e100;
		for (int i = 0; i < n; i++) {
			di = fabs(y1s[i] - s[i].re) + fabs(y2s[i] - s[i].im);
			di /= sqrt(q[i]);
			if (di < mindi) mindi = di;
		}
		if (mindi >= 10.) {

			mags[i] = 1.;
		}
		else {
			mags[i] = MultiMag(complex(y1s[i], y2s[i]), rho, Tol);
		}
	}
}

void VBMicrolensing::TripleLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), di, mindi, u, u0 = pr[2], t0 = pr[6], pai1 = pr[10], pai2 = pr[11];
	double q[3] = { 1, exp(pr[1]),exp(pr[8]) };
	complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);
	double Et[2];

	s[0] = exp(pr[0]) / (q[0] + q[1]);
	s[1] = s[0] * q[0];
	s[0] = -q[1] * s[0];
	s[2] = exp(pr[7]) * complex(cbeta, sbeta) + s[0];
	//	_sols *Images; double Mag; // For debugging

	SetLensGeometry(3, q, s);

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u = u0 + pai1 * Et[1] - pai2 * Et[0];
		y1s[i] = u * salpha - tn * calpha;
		y2s[i] = -u * calpha - tn * salpha;
		mindi = 1.e100;
		for (int i = 0; i < n; i++) {
			di = fabs(y1s[i] - s[i].re) + fabs(y2s[i] - s[i].im);
			di /= sqrt(q[i]);
			if (di < mindi) mindi = di;
		}
		if (mindi >= 10.) {

			mags[i] = 1.;
		}
		else {
			mags[i] = MultiMag(complex(y1s[i], y2s[i]), rho, Tol);
		}
	}
}

void VBMicrolensing::LightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np, int nl) {
	double rho = exp(pr[2]), tn, tE_inv = exp(-pr[1]), di, mindi;

	double* q = (double*)malloc(sizeof(double) * (nl));
	complex* s = (complex*)malloc(sizeof(complex) * (nl));

	q[0] = 1.;
	for (int i = 1, j = 1; i < nl; ++i, j += 3) {
		q[i] = pr[j + 5];
	}

	s[0] = complex(0, pr[3]);
	for (int i = 1, j = 1; i < nl; ++i, j += 3) {
		s[i] = complex(pr[j + 3], pr[j + 4]);
	}

	SetLensGeometry(nl, q, s);


	free(q);
	free(s);

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - pr[0]) * tE_inv;
		y1s[i] = -tn;
		y2s[i] = 0.;
		mindi = 1.e100;
		for (int i = 0; i < n; i++) {
			di = fabs(y1s[i] - a[i].re) + fabs(y2s[i] - a[i].im);
			di /= sqrt(m[i]);
			if (di < mindi) mindi = di;
		}
		if (mindi >= 10.) {

			mags[i] = 1.;
		}
		else {
			mags[i] = MultiMag(complex(y1s[i], y2s[i]), rho, Tol);
		}
	}

}


void VBMicrolensing::BinaryLightCurveW(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0, u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xc;

	xc = (s - 1 / s) / (1 + q);
	if (xc < 0) xc = 0.;
	t0 = pr[6] + xc * calpha / tE_inv;
	u0 = pr[2] + xc * salpha;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t0) * tE_inv;
		y1s[i] = u0 * salpha - tn * calpha;
		y2s[i] = -u0 * calpha - tn * salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);
	}
}


void VBMicrolensing::BinaryLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, u, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u = u0 + pai1 * Et[1] - pai2 * Et[0];
		y1s[i] = u * salpha - tn * calpha;
		y2s[i] = -u * calpha - tn * salpha;
		mags[i] = BinaryMag2(s, q, y1s[i], y2s[i], rho);
	}
}


void VBMicrolensing::BinaryLightCurveOrbital(double* pr, double* ts, double* mags, double* y1s, double* y2s, double* seps, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm, s_true;
	double w13, w123, den, den0, u;
	t0old = 0;

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		inc = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0; // orbital radius
	COm = (Cphi0 * calpha + Cinc * salpha * Sphi0) / den0;
	SOm = (Cphi0 * salpha - Cinc * calpha * Sphi0) / den0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		phi = (ts[i] - t0_par) * w + phi0;
		Cphi = cos(phi);
		Sphi = sin(phi);
		den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
		seps[i] = s_true * den; // projected separation at time ts[i]

		u = u0 + pai1 * Et[1] - pai2 * Et[0];
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		y1s[i] = (Cphi * (u * SOm - tn * COm) + Cinc * Sphi * (u * COm + tn * SOm)) / den;
		y2s[i] = (-Cphi * (u * COm + tn * SOm) - Cinc * Sphi * (tn * COm - u * SOm)) / den;
		mags[i] = BinaryMag2(seps[i], q, y1s[i], y2s[i], rho);
	}
}

void VBMicrolensing::BinaryLightCurveKepler(double* pr, double* ts, double* mags, double* y1s, double* y2s, double* seps, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], alpha = pr[3], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11], szs = pr[12], ar = pr[13] + 1.e-8;
	double Et[2];
	double u, w22, w11, w33, w12, w23, szs2, ar2, EE, dE;
	double wt2, smix, sqsmix, e, h, snu, co1EE0, co2EE0, cosE, sinE, co1tperi, tperi, EE0, M, a, St, psi, dM, conu, n;
	double arm1, arm2;
	double X[3], Y[3], Z[3], r[2], x[2];
	t0old = 0;

	smix = 1 + szs * szs;
	sqsmix = sqrt(smix);
	w22 = w2 * w2;
	w11 = w1 * w1;
	w33 = w3 * w3;
	w12 = w11 + w22;
	w23 = w22 + w33;
	wt2 = w12 + w33;

	szs2 = szs * szs;
	ar2 = ar * ar;
	arm1 = ar - 1;
	arm2 = 2 * ar - 1;
	//	n = sqrt(wt2) / (ar*sqrt(-1 + 2 * ar)*sqrt(smix));
	n = sqrt(wt2 / arm2 / smix) / ar;
	Z[0] = -szs * w2;
	Z[1] = szs * w1 - w3;
	Z[2] = w2;
	h = sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
	for (int i = 0; i < 3; i++) Z[i] /= h;
	X[0] = -ar * w11 + arm1 * w22 - arm2 * szs * w1 * w3 + arm1 * w33;
	X[1] = -arm2 * w2 * (w1 + szs * w3);
	X[2] = arm1 * szs * w12 - arm2 * w1 * w3 - ar * szs * w33;
	e = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	for (int i = 0; i < 3; i++) X[i] /= e;
	e /= ar * sqsmix * wt2;
	Y[0] = Z[1] * X[2] - Z[2] * X[1];
	Y[1] = Z[2] * X[0] - Z[0] * X[2];
	Y[2] = Z[0] * X[1] - Z[1] * X[0];

	//	h = sqrt((smix)*w22 + (szs*w1 - w3)*(szs*w1 - w3));
	//	co1e = (1 - ar)*(1 - ar) + ar2 * szs2 + (-1 + 2 * ar)*(w11*(1 - szs2) - szs2 * w22 + 2 * szs*w1*w3) / wt2;
	//	co1nu = ar2 * szs2 + arm2 * (w11*(1 - szs2) - szs2 * w22 + 2 * szs*w1*w3) / wt2;
	//	co1e = arm1*arm1 + co1nu;
	//	coe2 = ar2 * smix;
	//	e = sqrt(co1e/coe2);
	//	co1nu = (-1 + 2 * ar)*sqrt(smix)*(w1 + szs * w3)*(szs2*(w12)-2 * szs*w1*w3 + w23);
	//	co2nu = ar * e*h*(smix)*sqrt((smix))*wt2;
	conu = (X[0] + X[2] * szs) / sqsmix;
	co1EE0 = conu + e;
	co2EE0 = 1 + e * conu;
	cosE = co1EE0 / co2EE0;
	EE0 = acos(cosE);
	snu = (Y[0] + Y[2] * szs);
	EE0 *= (snu > 0) ? 1 : -1;
	sinE = sqrt(1 - cosE * cosE) * ((snu > 0) ? 1 : -1);
	co1tperi = e * sinE;
	tperi = t0_par - (EE0 - co1tperi) / n;
	//	coX = ar * e*sqrt(smix)*wt2;
		//coX1 = -ar * w11 + (-1 + ar)*w22 + (1 - 2 * ar)*szs*w1*w3 + (-1 + ar)*w33;
		//coX2 = (-1 + 2 * ar)*w1*w23 + szs2 * w1*((-1 + ar)*w12 - ar * w33) + szs * w3*((2 - 3 * ar)*w11 + ar * w23);
		//coY1 = -(-1 + 2 * ar)*w2*(w1 + szs * w3);
		//coY2 = w2 * (-szs2 * w12 + 2 * szs*w1*w3 - w23 + ar * (-4 * szs*w1*w3 + szs2 * (w12 - w33) + (-w11 + w23)));
	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);
		while (M > M_PI) M -= 2 * M_PI;
		while (M < -M_PI) M += 2 * M_PI;
		EE = M + e * sin(M);
		dE = 1;
		while (fabs(dE) > 1.e-8) {
			dM = M - (EE - e * sin(EE));
			dE = dM / (1 - e * cos(EE));
			EE += dE;
			if (EE > M_PI) EE = M_PI;
			if (EE < -M_PI) EE = -M_PI;
		}
		a = ar * s * sqrt(smix);

		r[0] = a * (cos(EE) - e);
		r[1] = a * sqrt(1 - e * e) * sin(EE);
		x[0] = r[0] * X[0] + r[1] * Y[0];  // (coX1*x[1] + coX2 * y[1] / h) / coX;
		x[1] = r[0] * X[1] + r[1] * Y[1];   //(coY1*x[1] + y[1] * coY2 / h) / coX;
		St = sqrt(x[0] * x[0] + x[1] * x[1]);
		psi = atan2(x[1], x[0]);// +((ar > 1) ? 0 : M_PI);
		u = u0 + pai1 * Et[1] - pai2 * Et[0];
		tn = (ts[i] - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		y1s[i] = -tn * cos(alpha + psi) + u * sin(alpha + psi);
		y2s[i] = -u * cos(alpha + psi) - tn * sin(alpha + psi);
		seps[i] = St;

		mags[i] = BinaryMag2(seps[i], q, y1s[i], y2s[i], rho);

	}
}


void VBMicrolensing::BinSourceLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t01) * tE_inv;
		u = tn * tn + u1 * u1;

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = (u + 2) / sqrt(u * (u + 4));

		tn = (ts[i] - t02) * tE_inv;
		u = tn * tn + u2 * u2;

		mags[i] += FR * (u + 2) / sqrt(u * (u + 4));
		mags[i] /= (1 + FR);

	}

}


void VBMicrolensing::BinSourceLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], w1 = pr[8], w2 = pr[9], w3 = pr[10];
	double Et[2];
	t0old = 0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		tn = (ts[i] - t01) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u0 = u1 + pai1 * Et[1] - pai2 * Et[0];
		u = tn * tn + u0 * u0;

		y1s[i] = -tn;
		y2s[i] = -u0;
		mags[i] = (u + 2) / sqrt(u * (u + 4));

		tn = (ts[i] - t02) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
		u0 = u2 + pai1 * Et[1] - pai2 * Et[0];
		u = tn * tn + u0 * u0;

		mags[i] += FR * (u + 2) / sqrt(u * (u + 4));
		mags[i] /= (1 + FR);
	}
}


void VBMicrolensing::BinSourceLightCurveXallarap(double* pr, double* ts, double* mags, double* y1s, double* y2s, double* seps, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], q = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double th, Cth, Sth;
	double Et[2];
	double s, s_true, w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm;
	double w13, w123, den, den0, du0, dt0;
	t0old = 0;

	s = sqrt((u1 - u2) * (u1 - u2) + (t01 - t02) * (t01 - t02) * (tE_inv * tE_inv));
	th = atan2((u1 - u2), (tE_inv * (t01 - t02)));
	Cth = cos(th);
	Sth = sin(th);
	u0 = (u1 + u2 * q) / (1 + q);
	t0 = (t01 + t02 * q) / (1 + q);

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		inc = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0;
	COm = (Cphi0 * Cth + Cinc * Sth * Sphi0) / den0;
	SOm = (Cphi0 * Sth - Cinc * Cth * Sphi0) / den0;


	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		phi = (ts[i] - t0_par) * w + phi0;
		Cphi = cos(phi);
		Sphi = sin(phi);
		den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
		seps[i] = s_true * den;

		dt0 = s_true * (COm * Cphi - Cinc * SOm * Sphi) / (1 + q) * q;  //Position of the primary component with respect to center of mass
		du0 = s_true * (SOm * Cphi + Cinc * COm * Sphi) / (1 + q) * q;

		tn = -((ts[i] - t0_par) * tE_inv + dt0 + pai1 * Et[0] + pai2 * Et[1]);
		u = -(u0 + du0 + pai1 * Et[1] - pai2 * Et[0]);
		y1s[i] = tn;
		y2s[i] = u;
		u = tn * tn + u * u;

		mags[i] = (u + 2) / sqrt(u * (u + 4));

		tn = -((ts[i] - t0_par) * tE_inv - dt0 / q + pai1 * Et[0] + pai2 * Et[1]); // Position of the secondary component
		u = -(u0 - du0 / q + pai1 * Et[1] - pai2 * Et[0]);
		u = tn * tn + u * u;

		mags[i] += FR * (u + 2) / sqrt(u * (u + 4));
		mags[i] /= (1 + FR);
	}
}

void VBMicrolensing::BinSourceExtLightCurve(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), rho = exp(pr[6]), rho2, tn, u;

	for (int i = 0; i < np; i++) {
		tn = (ts[i] - t01) * tE_inv;
		u = tn * tn + u1 * u1;

		y1s[i] = -tn;
		y2s[i] = -u1;
		mags[i] = ESPLMag2(sqrt(u), rho);

		tn = (ts[i] - t02) * tE_inv;
		u = tn * tn + u2 * u2;
		rho2 = rho * pow(FR, mass_radius_exponent / mass_luminosity_exponent);
		mags[i] += FR * ESPLMag2(sqrt(u), rho2);
		mags[i] /= (1 + FR);

	}

}

void VBMicrolensing::BinSourceBinLensXallarap(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xi1 = pr[7], xi2 = pr[8], omega = pr[9], inc = pr[10], phi = pr[11], qs = exp(pr[12]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, qs4;



	if (t0_par_fixed == 0) t0_par = pr[6];


	for (int i = 0; i < np; i++) {

		phit = omega * (ts[i] - t0_par);

		disp[0] = sin(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

		disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

		Xal[0] = xi1 * disp[0] + xi2 * disp[1];
		Xal[1] = xi2 * disp[0] - xi1 * disp[1];
		tn = (ts[i] - pr[6]) * tE_inv + Xal[0];
		u0 = pr[2] + Xal[1];
		y1s[i] = u0 * salpha - tn * calpha;
		y2s[i] = -u0 * calpha - tn * salpha;
		Mag = BinaryMag2(s, q, y1s[i], y2s[i], rho);

		disp2[0] = -sin(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

		disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

		Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
		Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
		tn2 = (ts[i] - pr[6]) * tE_inv + Xal2[0];
		u02 = pr[2] + Xal2[1];
		y1s2 = u02 * salpha - tn2 * calpha;
		y2s2 = -u02 * calpha - tn2 * salpha;
		rho2 = rho * pow(qs, mass_radius_exponent);
		Mag2 = BinaryMag2(s, q, y1s2, y2s2, rho2);
		qs4 = pow(qs, mass_luminosity_exponent);
		mags[i] = (Mag + qs4 * Mag2) / (1 + qs4);
	}
}

void VBMicrolensing::BinSourceSingleLensXallarap(double* pr, double* ts, double* mags, double* y1s, double* y2s, double* y1s2, double* y2s2, int np) {
	double  t0 = pr[1], rho = exp(pr[3]), tn, tE_inv = exp(-pr[2]), u0;
	double  xi1 = pr[4], xi2 = pr[5], omega = pr[6], inc = pr[7], phi = pr[8], qs = exp(pr[9]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, qs4, u, u2;



	t0_par = pr[1];


	for (int i = 0; i < np; i++) {

		phit = omega * (ts[i] - t0_par);

		disp[0] = cos(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

		disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

		Xal[0] = xi1 * disp[0] + xi2 * disp[1];
		Xal[1] = xi2 * disp[0] - xi1 * disp[1];
		tn = (ts[i] - pr[1]) * tE_inv + Xal[0];
		u0 = pr[0] + Xal[1];
		u = sqrt(tn * tn + u0 * u0);

		y1s[i] = -tn;
		y2s[i] = -u0;
		Mag = ESPLMag2(u, rho);  /*If you want only the second source put =0, otherwise replace ESPLMag2(u, rho);*/


		disp2[0] = -cos(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

		disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

		Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
		Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
		tn2 = (ts[i] - pr[1]) * tE_inv + Xal2[0];
		u02 = pr[0] + Xal2[1];
		u2 = sqrt(tn2 * tn2 + u02 * u02);
		y1s2[i] = -tn2;
		y2s2[i] = -u02;
		rho2 = rho * pow(qs, mass_radius_exponent);
		Mag2 = ESPLMag2(u2, rho2);  /*If you want only the second source put =0, otherwise replace ESPLMag2(u2, rho2);*/
		qs4 = pow(qs, mass_luminosity_exponent);
		mags[i] = (Mag + qs4 * Mag2) / (1 + qs4);
	}
}



//////////////////////////////
//////////////////////////////
////////Old (v1) light curve functions
//////////////////////////////
//////////////////////////////


double VBMicrolensing::PSPLLightCurve(double* pr, double t) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u;

	tn = (t - t0) * tE_inv;
	u = tn * tn + u0 * u0;

	y_1 = -tn;
	y_2 = -u0;
	return (u + 2) / sqrt(u * (u + 4));

}


double VBMicrolensing::PSPLLightCurveParallax(double* pr, double t) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, pai1 = pr[3], pai2 = pr[4];
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u1 = u0 + pai1 * Et[1] - pai2 * Et[0];
	u = tn * tn + u1 * u1;

	y_1 = -tn;
	y_2 = -u1;
	return (u + 2) / sqrt(u * (u + 4));

}


double VBMicrolensing::ESPLLightCurve(double* pr, double t) {
	double u0 = exp(pr[0]), t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, rho = exp(pr[3]);

	tn = (t - t0) * tE_inv;
	u = sqrt(tn * tn + u0 * u0);

	y_1 = -tn;
	y_2 = -u0;
	return ESPLMag2(u, rho);

}

double VBMicrolensing::ESPLLightCurveParallax(double* pr, double t) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, rho = exp(pr[3]), pai1 = pr[4], pai2 = pr[5];
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u1 = u0 + pai1 * Et[1] - pai2 * Et[0];
	u = sqrt(tn * tn + u1 * u1);

	y_1 = -tn;
	y_2 = -u1;
	return ESPLMag2(u, rho);

}


double VBMicrolensing::BinaryLightCurve(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]);
	double salpha = sin(pr[3]), calpha = cos(pr[3]);

	tn = (t - pr[6]) * tE_inv;
	y_1 = pr[2] * salpha - tn * calpha;
	y_2 = -pr[2] * calpha - tn * salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}

double VBMicrolensing::TripleLightCurve(double* pr, double t) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), di, mindi;
	static double q[3];
	static double prold[] = { 0,0,0,0,0 };
	static Method oldmethod = Method::Nopoly;
	static complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);
	int inew = 0;
	bool changed = false;
	for (int i = 0; i < 5; i++) {
		if (pr[inew] != prold[i]) {
			changed = true;
			prold[i] = pr[inew];
		}
		inew += (i == 1) ? 6 : 1;
	}
	if (SelectedMethod != oldmethod) {
		changed = true;
		oldmethod = SelectedMethod;
	}
	if (changed) {
		q[0] = 1;
		q[1] = exp(pr[1]);
		q[2] = exp(pr[8]);
		s[0] = exp(pr[0]) / (q[0] + q[1]);
		s[1] = s[0] * q[0];
		s[0] = -q[1] * s[0];
		s[2] = exp(pr[7]) * complex(cbeta, sbeta) + s[0];
		SetLensGeometry(3, q, s);
	}

	tn = (t - pr[6]) * tE_inv;
	y_1 = pr[2] * salpha - tn * calpha;
	y_2 = -pr[2] * calpha - tn * salpha;

	mindi = 1.e100;
	for (int i = 0; i < n; i++) {
		di = fabs(y_1 - s[i].re) + fabs(y_2 - s[i].im);
		di /= sqrt(q[i]);
		if (di < mindi) mindi = di;
	}
	if (mindi >= 10.) {

		return 1.;
	}
	else {
		return MultiMag(complex(y_1, y_2), rho, Tol);
	}
}

double VBMicrolensing::TripleLightCurveParallax(double* pr, double t) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), di, mindi, u, u0 = pr[2], t0 = pr[6], pai1 = pr[10], pai2 = pr[11];
	double q[3] = { 1, exp(pr[1]), exp(pr[8]) };
	static complex s[3];
	static double prold[] = { 0,0,0,0,0 };
	static Method oldmethod = Method::Nopoly;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);
	static double Et[2];
	int inew = 0;
	bool changed = false;
	for (int i = 0; i < 5; i++) {
		if (pr[inew] != prold[i]) {
			changed = true;
			prold[i] = pr[inew];
		}
		inew += (i == 1) ? 6 : 1;
	}
	if (SelectedMethod != oldmethod) {
		changed = true;
		oldmethod = SelectedMethod;
	}
	if (changed) {
		q[0] = 1;
		q[1] = exp(pr[1]);
		q[2] = exp(pr[8]);
		s[0] = exp(pr[0]) / (q[0] + q[1]);
		s[1] = s[0] * q[0];
		s[0] = -q[1] * s[0];
		s[2] = exp(pr[7]) * complex(cbeta, sbeta) + s[0];
		SetLensGeometry(3, q, s);
	}

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	y_1 = u * salpha - tn * calpha;
	y_2 = -u * calpha - tn * salpha;
	mindi = 1.e100;
	for (int i = 0; i < n; i++) {
		di = fabs(y_1 - s[i].re) + fabs(y_2 - s[i].im);
		di /= sqrt(q[i]);
		if (di < mindi) mindi = di;
	}
	if (mindi >= 10.) {

		return 1.;
	}
	else {
		return MultiMag(complex(y_1, y_2), rho, Tol);
	}
}

double VBMicrolensing::BinaryLightCurveW(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0, u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xc;

	xc = (s - 1 / s) / (1 + q);
	if (xc < 0) xc = 0.;
	t0 = pr[6] + xc * calpha / tE_inv;
	u0 = pr[2] + xc * salpha;

	tn = (t - t0) * tE_inv;
	y_1 = u0 * salpha - tn * calpha;
	y_2 = -u0 * calpha - tn * salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}


double VBMicrolensing::BinaryLightCurveParallax(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, u, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	y_1 = u * salpha - tn * calpha;
	y_2 = -u * calpha - tn * salpha;
	return BinaryMag2(s, q, y_1, y_2, rho);

}


double VBMicrolensing::BinaryLightCurveOrbital(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm, s_true;
	double w13, w123, den, den0, u;

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		inc = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}

	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0;
	COm = (Cphi0 * calpha + Cinc * salpha * Sphi0) / den0;
	SOm = (Cphi0 * salpha - Cinc * calpha * Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phi = (t - t0_par) * w + phi0;
	Cphi = cos(phi);
	Sphi = sin(phi);
	den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
	av = s_true * den;

	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	y_1 = (Cphi * (u * SOm - tn * COm) + Cinc * Sphi * (u * COm + tn * SOm)) / den;
	y_2 = (-Cphi * (u * COm + tn * SOm) - Cinc * Sphi * (tn * COm - u * SOm)) / den;
	return BinaryMag2(av, q, y_1, y_2, rho);
}

double VBMicrolensing::BinaryLightCurveKepler(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], alpha = pr[3], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11], szs = pr[12], ar = pr[13] + 1.e-8;
	double Et[2];
	double u, w22, w11, w33, w12, w23, szs2, ar2, EE, dE;
	double wt2, smix, sqsmix, e, h, snu, co1EE0, co2EE0, cosE, sinE, co1tperi, tperi, EE0, M, a, St, psi, dM, conu, n;
	double arm1, arm2;
	double X[3], Y[3], Z[3], r[2], x[2];
	t0old = 0;

	smix = 1 + szs * szs;
	sqsmix = sqrt(smix);
	w22 = w2 * w2;
	w11 = w1 * w1;
	w33 = w3 * w3;
	w12 = w11 + w22;
	w23 = w22 + w33;
	wt2 = w12 + w33;

	szs2 = szs * szs;
	ar2 = ar * ar;
	arm1 = ar - 1;
	arm2 = 2 * ar - 1;
	n = sqrt(wt2 / arm2 / smix) / ar;
	Z[0] = -szs * w2;
	Z[1] = szs * w1 - w3;
	Z[2] = w2;
	h = sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
	for (int i = 0; i < 3; i++) Z[i] /= h;
	X[0] = -ar * w11 + arm1 * w22 - arm2 * szs * w1 * w3 + arm1 * w33;
	X[1] = -arm2 * w2 * (w1 + szs * w3);
	X[2] = arm1 * szs * w12 - arm2 * w1 * w3 - ar * szs * w33;
	e = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	for (int i = 0; i < 3; i++) X[i] /= e;
	e /= ar * sqsmix * wt2;
	Y[0] = Z[1] * X[2] - Z[2] * X[1];
	Y[1] = Z[2] * X[0] - Z[0] * X[2];
	Y[2] = Z[0] * X[1] - Z[1] * X[0];

	conu = (X[0] + X[2] * szs) / sqsmix;
	co1EE0 = conu + e;
	co2EE0 = 1 + e * conu;
	cosE = co1EE0 / co2EE0;
	EE0 = acos(cosE);
	snu = (Y[0] + Y[2] * szs);
	EE0 *= (snu > 0) ? 1 : -1;
	sinE = sqrt(1 - cosE * cosE) * ((snu > 0) ? 1 : -1);
	co1tperi = e * sinE;
	tperi = t0_par - (EE0 - co1tperi) / n;

	ComputeParallax(t, t0, Et);
	M = n * (t - tperi);
	while (M > M_PI) M -= 2 * M_PI;
	while (M < -M_PI) M += 2 * M_PI;
	EE = M + e * sin(M);
	dE = 1;
	while (fabs(dE) > 1.e-8) {
		dM = M - (EE - e * sin(EE));
		dE = dM / (1 - e * cos(EE));
		EE += dE;
		if (EE > M_PI) EE = M_PI;
		if (EE < -M_PI) EE = -M_PI;
	}
	a = ar * s * sqrt(smix);

	r[0] = a * (cos(EE) - e);
	r[1] = a * sqrt(1 - e * e) * sin(EE);
	x[0] = r[0] * X[0] + r[1] * Y[0];  // (coX1*x[1] + coX2 * y[1] / h) / coX;
	x[1] = r[0] * X[1] + r[1] * Y[1];   //(coY1*x[1] + y[1] * coY2 / h) / coX;
	St = sqrt(x[0] * x[0] + x[1] * x[1]);
	psi = atan2(x[1], x[0]);// +((ar > 1) ? 0 : M_PI);

	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	y_1 = -tn * cos(alpha + psi) + u * sin(alpha + psi);
	y_2 = -u * cos(alpha + psi) - tn * sin(alpha + psi);

	return BinaryMag2(St, q, y_1, y_2, rho);

}


double VBMicrolensing::BinSourceLightCurve(double* pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, mag;

	tn = (t - t01) * tE_inv;
	u = tn * tn + u1 * u1;

	y_1 = -tn;
	y_2 = -u1;
	mag = (u + 2) / sqrt(u * (u + 4));

	tn = (t - t02) * tE_inv;
	u = tn * tn + u2 * u2;

	mag += FR * (u + 2) / sqrt(u * (u + 4));
	mag /= (1 + FR);

	return mag;

}


double VBMicrolensing::BinSourceLightCurveParallax(double* pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], w1 = pr[8], w2 = pr[9], w3 = pr[10];
	double Et[2], mag;

	ComputeParallax(t, t0, Et);

	tn = (t - t01) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u0 = u1 + pai1 * Et[1] - pai2 * Et[0];
	u = tn * tn + u0 * u0;

	y_1 = -tn;
	y_2 = -u0;
	mag = (u + 2) / sqrt(u * (u + 4));

	tn = (t - t02) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u0 = u2 + pai1 * Et[1] - pai2 * Et[0];
	u = tn * tn + u0 * u0;

	mag += FR * (u + 2) / sqrt(u * (u + 4));
	mag /= (1 + FR);

	return mag;
}


double VBMicrolensing::BinSourceLightCurveXallarap(double* pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), tn, u, u0, pai1 = pr[6], pai2 = pr[7], q = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double th, Cth, Sth;
	double Et[2], mag;
	double s, s_true, w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm;
	double w13, w123, den, den0, du0, dt0;

	s = sqrt((u1 - u2) * (u1 - u2) + (t01 - t02) * (t01 - t02) * (tE_inv * tE_inv));
	th = atan2((u1 - u2), (tE_inv * (t01 - t02)));
	Cth = cos(th);
	Sth = sin(th);
	u0 = (u1 + u2 * q) / (1 + q);
	t0 = (t01 + t02 * q) / (1 + q);

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		inc = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		inc = 0.;
		phi0 = 0.;
	}
	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(inc);
	Sinc = sin(inc);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0;
	COm = (Cphi0 * Cth + Cinc * Sth * Sphi0) / den0;
	SOm = (Cphi0 * Sth - Cinc * Cth * Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phi = (t - t0_par) * w + phi0;
	Cphi = cos(phi);
	Sphi = sin(phi);
	den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
	av = s_true * den;

	dt0 = s_true * (COm * Cphi - Cinc * SOm * Sphi) / (1 + q) * q;  //Position of the primary component with respect to center of mass
	du0 = s_true * (SOm * Cphi + Cinc * COm * Sphi) / (1 + q) * q;

	tn = -((t - t0_par) * tE_inv - dt0 + pai1 * Et[0] + pai2 * Et[1]);
	u = -(u0 + du0 + pai1 * Et[1] - pai2 * Et[0]);
	y_1 = tn;
	y_2 = u;
	u = tn * tn + u * u;

	mag = (u + 2) / sqrt(u * (u + 4));

	tn = -((t - t0_par) * tE_inv + dt0 / q + pai1 * Et[0] + pai2 * Et[1]); // Position of the secondary component
	u = -(u0 - du0 / q + pai1 * Et[1] - pai2 * Et[0]);
	u = tn * tn + u * u;

	mag += FR * (u + 2) / sqrt(u * (u + 4));
	mag /= (1 + FR);

	return mag;

}

double VBMicrolensing::BinSourceExtLightCurve(double* pr, double t) {
	double u1 = pr[2], u2 = pr[3], t01 = pr[4], t02 = pr[5], tE_inv = exp(-pr[0]), FR = exp(pr[1]), rho = exp(pr[6]), rho2, tn, u, mag;

	tn = (t - t01) * tE_inv;
	u = tn * tn + u1 * u1;

	y_1 = -tn;
	y_2 = -u1;
	mag = ESPLMag2(sqrt(u), rho);

	tn = (t - t02) * tE_inv;
	u = tn * tn + u2 * u2;
	rho2 = rho * pow(FR, mass_radius_exponent / mass_luminosity_exponent);

	mag += FR * ESPLMag2(sqrt(u), rho2);
	mag /= (1 + FR);

	return mag;

}


double VBMicrolensing::BinSourceBinLensXallarap(double* pr, double t) {

	double s = exp(pr[0]), q = exp(pr[1]), rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), u0;
	double salpha = sin(pr[3]), calpha = cos(pr[3]), xi1 = pr[7], xi2 = pr[8], omega = pr[9], inc = pr[10], phi = pr[11], qs = exp(pr[12]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, y1s, y2s, mags, qs4;



	if (t0_par_fixed == 0) t0_par = pr[6];




	phit = omega * (t - t0_par);

	disp[0] = cos(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	tn = (t - pr[6]) * tE_inv + Xal[0];
	u0 = pr[2] + Xal[1];
	y1s = u0 * salpha - tn * calpha;
	y2s = -u0 * calpha - tn * salpha;
	Mag = BinaryMag2(s, q, y1s, y2s, rho);

	disp2[0] = -cos(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
	tn2 = (t - pr[6]) * tE_inv + Xal2[0];
	u02 = pr[2] + Xal2[1];
	y1s2 = u02 * salpha - tn2 * calpha;
	y2s2 = -u02 * calpha - tn2 * salpha;
	rho2 = rho * pow(qs, mass_radius_exponent);
	Mag2 = BinaryMag2(s, q, y1s2, y2s2, rho2);
	qs4 = pow(qs, mass_luminosity_exponent);
	mags = (Mag + qs4 * Mag2) / (1 + qs4);

	return mags;


}

double VBMicrolensing::BinSourceSingleLensXallarap(double* pr, double t) {

	double  t0 = pr[1], rho = exp(pr[3]), tn, tE_inv = exp(-pr[2]), u0;
	double  xi1 = pr[4], xi2 = pr[5], omega = pr[6], inc = pr[7], phi = pr[8], qs = exp(pr[9]);

	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u02, rho2, tn2, y1s2, y2s2, qs4, u, y1s, y2s, mags, u2;



	t0_par = pr[1];

	phit = omega * (t - t0_par);

	disp[0] = cos(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));

	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	tn = (t - pr[1]) * tE_inv + Xal[0];
	u0 = pr[0] + Xal[1];
	u = sqrt(tn * tn + u0 * u0);

	y1s = -tn;
	y2s = -u0;
	Mag = ESPLMag2(u, rho); /*If you want only the second source put =0, otherwise replace ESPLMag2(u, rho);*/


	disp2[0] = -cos(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));

	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];
	tn2 = (t - pr[1]) * tE_inv + Xal2[0];
	u02 = pr[0] + Xal2[1];
	u2 = sqrt(tn2 * tn2 + u02 * u02);
	y1s2 = -tn2;
	y2s2 = -u02;
	rho2 = rho * pow(qs, mass_radius_exponent);
	Mag2 = ESPLMag2(u2, rho2); /*If you want only the second source put =0, otherwise replace ESPLMag2(u2, rho2);*/
	qs4 = pow(qs, mass_luminosity_exponent);
	mags = (Mag + qs4 * Mag2) / (1 + qs4);
	return mags;
}

double VBMicrolensing::BinSourceBinLensPOX(double* pr, double t) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tE_inv = exp(-pr[5]), t0 = pr[6];
	double pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double tn, w, phi0, phil, incl, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm, s_true;
	double w13, w123, den, den0, u;

	double xi1 = pr[12], xi2 = pr[13], omega = pr[14], inc = pr[15], phi = pr[16], qs = exp(pr[17]);
	double Xal[2], phit, disp[2], Xal2[2], disp2[2];
	double Mag, Mag2, u01, u02, rho2, tn1, tn2, mags, qs4;

	w13 = w1 * w1 + w3 * w3;
	w123 = sqrt(w13 + w2 * w2);
	w13 = sqrt(w13);
	if (w13 > 1.e-8) {
		w3 = (w3 > 1.e-8) ? w3 : 1.e-8;
		w = w3 * w123 / w13;
		incl = acos(w2 * w3 / w13 / w123);
		phi0 = atan2(-w1 * w123, w3 * w13);
	}
	else {
		w = w2;
		incl = 0.;
		phi0 = 0.;
	}

	Cphi0 = cos(phi0);
	Sphi0 = sin(phi0);
	Cinc = cos(incl);
	Sinc = sin(incl);
	den0 = sqrt(Cphi0 * Cphi0 + Cinc * Cinc * Sphi0 * Sphi0);
	s_true = s / den0;
	COm = (Cphi0 * calpha + Cinc * salpha * Sphi0) / den0;
	SOm = (Cphi0 * salpha - Cinc * calpha * Sphi0) / den0;

	ComputeParallax(t, t0, Et);

	phil = (t - t0_par) * w + phi0;
	Cphi = cos(phil);
	Sphi = sin(phil);
	den = sqrt(Cphi * Cphi + Cinc * Cinc * Sphi * Sphi);
	av = s_true * den;
	u = u0 + pai1 * Et[1] - pai2 * Et[0];
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];

	phit = omega * (t - t0_par);

	disp[0] = sin(inc) * (-cos(phi) + cos(phi + phit) + phit * sin(phi));
	disp[1] = -phit * cos(phi) - sin(phi) + sin(phi + phit);
	disp2[0] = -sin(inc) * (cos(phi) + cos(phi + phit) / qs - phit * sin(phi));
	disp2[1] = phit * cos(phi) + sin(phi) + sin(phi + phit) / qs;

	Xal[0] = xi1 * disp[0] + xi2 * disp[1];
	Xal[1] = xi2 * disp[0] - xi1 * disp[1];
	Xal2[0] = xi1 * disp2[0] - xi2 * disp2[1];
	Xal2[1] = xi2 * disp2[0] + xi1 * disp2[1];

	tn1 = tn + Xal[0];
	u01 = u + Xal[1];
	tn2 = tn + Xal2[0];
	u02 = u + Xal2[1];
	rho2 = rho * pow(qs, mass_radius_exponent);
	qs4 = pow(qs, mass_luminosity_exponent);


	/*	y1s = u01 * salpha - tn1 * calpha;
		y2s = -u01 * calpha - tn1 * salpha;
		y1s2 = u02 * salpha - tn2 * calpha;
		y2s2 = -u02 * calpha - tn2 * salpha;*/

	y_1 = (Cphi * (u02 * SOm - tn2 * COm) + Cinc * Sphi * (u02 * COm + tn2 * SOm)) / den;
	y_2 = (-Cphi * (u02 * COm + tn2 * SOm) - Cinc * Sphi * (tn2 * COm - u02 * SOm)) / den;
	Mag2 = BinaryMag2(av, q, y_1, y_2, rho2);

	y_1 = (Cphi * (u01 * SOm - tn1 * COm) + Cinc * Sphi * (u01 * COm + tn1 * SOm)) / den;
	y_2 = (-Cphi * (u01 * COm + tn1 * SOm) - Cinc * Sphi * (tn1 * COm - u01 * SOm)) / den;
	Mag = BinaryMag2(av, q, y_1, y_2, rho);

	mags = (Mag + qs4 * Mag2) / (1 + qs4);

	return mags;
}



#pragma endregion

#pragma region limbdarkening


double VBMicrolensing::LDprofile(double r) {
	static int ir;
	static double rr, ret;
	switch (curLDprofile) {
	case LDuser:
		rr = r * npLD;
		ir = (int)rr;
		rr -= ir;
		ret = LDtab[ir] * (1 - rr) + LDtab[ir + 1] * rr;
		break;
	case LDlinear:
		ret = 3 / (3 - a1) * (1 - a1 * scr2);
		break;
	case LDsquareroot:
		ret = 3 / (3 - a1 - 0.6 * a2) * (1 - a1 * scr2 - a2 * sscr2);
	case LDquadratic:
		ret = 3 / (3 - a1 - 0.5 * a2) * (1 - a1 * scr2 - a2 * sscr2);
		break;
	case LDlog:
		ret = 3 / (3 - a1 + 0.666666666666 * a2) * (1 - a1 * scr2 - a2 * sscr2);
		break;
	}
	return ret;
}

double VBMicrolensing::rCLDprofile(double tc, annulus * left, annulus * right) {
	static int ic;
	static double rc, cb, lc, r2, cr2, cc, lb, rb;

	switch (curLDprofile) {
	case LDuser:
		rc = tc * npLD;
		ic = (int)rc;
		rc -= ic;
		cb = rCLDtab[ic] * (1 - rc) + rCLDtab[ic + 1] * rc;
		break;
	case LDlinear:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc) * (rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = 1 - sqrt(cr2);
			cc = (3 * r2 - a1 * (r2 - 2 * scr2 * cr2)) / (3 - a1);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDsquareroot:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc) * (rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = sqrt(cr2);
			sscr2 = 1 - sqrt(scr2);
			scr2 = 1 - scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2 * cr2) - 0.6 * a2 * (r2 - 4 * sscr2 * cr2)) / (3 - a1 - 0.6 * a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDquadratic:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc) * (rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = 1 - sqrt(cr2);
			sscr2 = scr2 * scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2 * cr2) + a2 * (4 * scr2 - (2 + 4 * scr2) * r2 + 1.5 * r2 * r2)) / (3 - a1 - 0.5 * a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	case LDlog:
		lb = left->bin;
		rb = right->bin;
		lc = left->cum;
		rc = right->cum;
		do {
			cb = rb + (tc - rc) * (rb - lb) / (rc - lc);
			r2 = cb * cb;
			cr2 = 1 - r2;
			scr2 = sqrt(cr2);
			sscr2 = scr2 * log(scr2);
			scr2 = 1 - scr2;
			cc = (3 * r2 - a1 * (r2 - 2 * scr2 * cr2) + 2 * a2 * (scr2 * (1 + scr2 * (scr2 / 3 - 1)) + sscr2 * cr2)) / (3 - a1 + 0.6666666666666666 * a2);
			if (cc > tc) {
				rb = cb;
				rc = cc;
			}
			else {
				lb = cb;
				lc = cc;
			}
		} while (fabs(cc - tc) > 1.e-5);
		break;
	}

	return cb;
}

void VBMicrolensing::SetLDprofile(double (*UserLDprofile)(double), int newnpLD) {
	int ic, ir;
	if (npLD > 0) {
		free(LDtab);
		free(rCLDtab);
	}
	if (newnpLD > 0) {
		npLD = newnpLD;
		double npLD2 = npLD * npLD;
		LDtab = (double*)malloc(sizeof(double) * (npLD + 1));
		CLDtab = (double*)malloc(sizeof(double) * (npLD + 1));
		rCLDtab = (double*)malloc(sizeof(double) * (npLD + 1));

		LDtab[0] = UserLDprofile(0.);
		CLDtab[0] = 0.;
		for (int i = 1; i <= npLD; i++) {
			LDtab[i] = UserLDprofile(((double)i) / npLD);
			CLDtab[i] = CLDtab[i - 1] + (LDtab[i] * i + LDtab[i - 1] * (i - 1));
		}
		for (int i = 0; i <= npLD; i++) {
			LDtab[i] *= npLD2 / CLDtab[npLD];
			CLDtab[i] /= CLDtab[npLD];
		}
		ic = 1;
		rCLDtab[0] = 0;
		ir = 1;
		while (ic < npLD) {
			while (CLDtab[ir] * npLD < ic && ir < npLD) ir++;
			rCLDtab[ic] = ((CLDtab[ir] - ((double)ic) / npLD) * (ir - 1) + (((double)ic) / npLD - CLDtab[ir - 1]) * ir) / (CLDtab[ir] - CLDtab[ir - 1]) / npLD;
			ic++;
		}
		rCLDtab[npLD] = 1;


		//printf("\n\n--------------------");
		//annulus left, right;
		//left.cum = left.bin=0;
		//right.cum = right.bin=1;
		//for (int i = 0; i <= npLD; i++) {
		//	double rl, fl;
		//	rl = ((double)i) / npLD;
		//	scr2 = 1 - sqrt(1 - rl * rl);
		//	fl = LDprofile(rl);
		//	rl = rCLDprofile(((double) i) / npLD,&left,&right);
		//	printf("\n%lf %lf %lf %lf", fl, LDtab[i], rl, rCLDtab[i]);
		//}
		//printf("\n--------------------\n\n");

		free(CLDtab);
		curLDprofile = LDuser;
	}
	else {
		npLD = 0;
		curLDprofile = LDlinear;
	}
}

void VBMicrolensing::SetLDprofile(LDprofiles LDval) {
	if (npLD > 0) {
		npLD = 0;
		free(LDtab);
		free(rCLDtab);
	}
	curLDprofile = LDval;
}

#pragma endregion

#pragma region parallax


//////////////////////////////
//////////////////////////////
////////Parallax settings and computation
//////////////////////////////
//////////////////////////////


void VBMicrolensing::SetObjectCoordinates(char* modelfile, char* sateltabledir) {
	double RA, Dec, dis, phiprec;
	FILE* f;
	char CoordinateString[512];
	char filename[256];
	int ic;

	f = fopen(modelfile, "r");
	if (f != 0) {
		fscanf(f, "%[^\n]s", CoordinateString);
		fclose(f);
		SetObjectCoordinates(CoordinateString);

		// Looking for satellite table files in the specified directory
		sprintf(filename, "%s%csatellite*.txt", sateltabledir, systemslash);
		nsat = 0;
		for (unsigned char c = 32; c < 255; c++) {
			filename[strlen(filename) - 5] = c;
			f = fopen(filename, "r");
			if (f != 0) {
				nsat++;
				fclose(f);
			}
		}


		tsat = (double**)malloc(sizeof(double*) * nsat);
		possat = (double***)malloc(sizeof(double**) * nsat);
		ndatasat = (int*)malloc(sizeof(int) * nsat);

		// Reading satellite table files
		ic = 0;
		for (unsigned char c = 32; c < 255; c++) {
			filename[strlen(filename) - 5] = c;
			f = fopen(filename, "r");
			if (f != 0) {
				int flag2 = 0;
				long startpos = 0;
				char teststring[1000];
				ndatasat[ic] = 1;

				// Finding start of data
				while (!feof(f)) {
					fscanf(f, "%s", teststring);
					if (!feof(f)) {
						fseek(f, 1, SEEK_CUR);
						teststring[5] = 0;
						if (strcmp(teststring, "$$SOE") == 0) {
							flag2 = 1;
							break;
						}
					}
				}
				// Finding end of data
				if (flag2) {
					flag2 = 0;
					startpos = ftell(f);
					while (!feof(f)) {
						fscanf(f, "%[^\n]s", teststring);
						if (!feof(f)) {
							fseek(f, 1, SEEK_CUR);
							teststring[5] = 0;
							if (strcmp(teststring, "$$EOE") == 0) {
								flag2 = 1;
								break;
							}
							else {
								ndatasat[ic]++;
							}
						}
					}
				}

				// Allocating memory according to the length of the table
				tsat[ic] = (double*)malloc(sizeof(double) * ndatasat[ic]);
				possat[ic] = (double**)malloc(sizeof(double*) * ndatasat[ic]);
				for (int j = 0; j < ndatasat[ic]; j++) {
					possat[ic][j] = (double*)malloc(sizeof(double) * 3);
				}
				ndatasat[ic]--;

				// Reading data
				if (f) {
					fseek(f, startpos, SEEK_SET);
					for (int id = 0; id < ndatasat[ic]; id++) {

						if (fscanf(f, "%lf %lf %lf %lf %lf", &(tsat[ic][id]), &RA, &Dec, &dis, &phiprec) == 5) {
							tsat[ic][id] -= 2450000;
							RA *= M_PI / 180;
							Dec *= M_PI / 180;
							for (int i = 0; i < 3; i++) {
								possat[ic][id][i] = dis * (cos(RA) * cos(Dec) * Eq2000[i] + sin(RA) * cos(Dec) * Quad2000[i] + sin(Dec) * North2000[i]);
							}
						}
						else {
							ndatasat[ic] = id;
							break;
						}
					}
					fclose(f);
				}

				ic++;
			}
		}
	}
	else {
		printf("\nFile not found!\n");
	}
}

void VBMicrolensing::SetObjectCoordinates(char* CoordinateString) {
	double RA, Dec, hr, mn, sc, deg, pr, ssc;

	if (nsat) {
		for (int i = 0; i < nsat; i++) {
			for (int j = 0; j < ndatasat[i]; j++) free(possat[i][j]);
			free(tsat[i]);
			free(possat[i]);
		}
		free(tsat);
		free(possat);
		free(ndatasat);
	}
	sscanf(CoordinateString, "%lf:%lf:%lf %lf:%lf:%lf", &hr, &mn, &sc, &deg, &pr, &ssc);
	RA = (hr + mn / 60 + sc / 3600) * M_PI / 12,
		Dec = (fabs(deg) + pr / 60 + ssc / 3600) * M_PI / 180;
	if (deg < 0) Dec = -Dec;

	for (int i = 0; i < 3; i++) {
		Obj[i] = (cos(RA) * cos(Dec) * Eq2000[i] + sin(RA) * cos(Dec) * Quad2000[i] + sin(Dec) * North2000[i]);
		rad[i] = Eq2000[i];
		tang[i] = North2000[i];
	}

	if (t0_par_fixed == -1) t0_par_fixed = 0;

}

void VBMicrolensing::ComputeParallax(double t, double t0, double* Et) {
	static double a0 = 1.00000261, adot = 0.00000562; // Ephemeris from JPL website 
	static double e0 = 0.01671123, edot = -0.00004392;
	static double inc0 = -0.00001531, incdot = -0.01294668;
	static double L0 = 100.46457166, Ldot = 35999.37244981;
	static double om0 = 102.93768193, omdot = 0.32327364;
	static double deg = M_PI / 180;
	static double a, e, inc, L, om, M, EE, dE, dM;
	static double x1, y1, vx, vy, Ear[3], vEar[3];
	static double Et0[2], vt0[2], r, sp, ty, Spit;
	int c = 0, ic;

	if (t0_par_fixed == 0) t0_par = t0;
	if (t0_par_fixed == -1) {
		printf("\nUse SetObjectCoordinates to input target coordinates");
	}
	else {

		if (t0_par != t0old) {
			t0old = t0_par;
			ty = (t0_par - 1545) / 36525.0;

			a = a0 + adot * ty;
			e = e0 + edot * ty;
			inc = (inc0 + incdot * ty) * deg;
			L = (L0 + Ldot * ty) * deg;
			om = (om0 + omdot * ty) * deg;

			M = L - om;
			M -= floor((M + M_PI) / (2 * M_PI)) * 2 * M_PI;

			EE = M + e * sin(M);
			dE = 1;
			while (fabs(dE) > 1.e-8) {
				dM = M - (EE - e * sin(EE));
				dE = dM / (1 - e * cos(EE));
				EE += dE;
			}
			x1 = a * (cos(EE) - e);
			y1 = a * sqrt(1 - e * e) * sin(EE);
			//		r=a*(1-e*cos(EE));
			vx = -a / (1 - e * cos(EE)) * sin(EE) * Ldot * deg / 36525;
			vy = a / (1 - e * cos(EE)) * cos(EE) * sqrt(1 - e * e) * Ldot * deg / 36525;

			Ear[0] = x1 * cos(om) - y1 * sin(om);
			Ear[1] = x1 * sin(om) * cos(inc) + y1 * cos(om) * cos(inc);
			Ear[2] = x1 * sin(om) * sin(inc) + y1 * cos(om) * sin(inc);
			vEar[0] = vx * cos(om) - vy * sin(om);
			vEar[1] = vx * sin(om) * cos(inc) + vy * cos(om) * cos(inc);
			vEar[2] = vx * sin(om) * sin(inc) + vy * cos(om) * sin(inc);

			sp = 0;
			switch (parallaxsystem) {
			case 1:
				for (int i = 0; i < 3; i++) sp += North2000[i] * Obj[i];
				for (int i = 0; i < 3; i++) rad[i] = -North2000[i] + sp * Obj[i];
				break;
			default:
				for (int i = 0; i < 3; i++) sp += Ear[i] * Obj[i];
				for (int i = 0; i < 3; i++) rad[i] = Ear[i] - sp * Obj[i];
				break;
			}

			r = sqrt(rad[0] * rad[0] + rad[1] * rad[1] + rad[2] * rad[2]);
			rad[0] /= r;
			rad[1] /= r;
			rad[2] /= r;
			tang[0] = rad[1] * Obj[2] - rad[2] * Obj[1];
			tang[1] = rad[2] * Obj[0] - rad[0] * Obj[2];
			tang[2] = rad[0] * Obj[1] - rad[1] * Obj[0];

			Et0[0] = Et0[1] = vt0[0] = vt0[1] = 0;
			for (int i = 0; i < 3; i++) {
				Et0[0] += Ear[i] * rad[i];
				Et0[1] += Ear[i] * tang[i];
				vt0[0] += vEar[i] * rad[i];
				vt0[1] += vEar[i] * tang[i];
			}
		}

		ty = (t - 1545) / 36525.0;

		a = a0 + adot * ty;
		e = e0 + edot * ty;
		inc = (inc0 + incdot * ty) * deg;
		L = (L0 + Ldot * ty) * deg;
		om = (om0 + omdot * ty) * deg;

		M = L - om;
		M -= floor((M + M_PI) / (2 * M_PI)) * 2 * M_PI;

		EE = M + e * sin(M);
		dE = 1;
		while (dE > 1.e-8) {
			dM = M - (EE - e * sin(EE));
			dE = dM / (1 - e * cos(EE));
			EE += dE;
		}
		x1 = a * (cos(EE) - e);
		y1 = a * sqrt(1 - e * e) * sin(EE);
		//	r=a*(1-e*cos(EE));

		Ear[0] = x1 * cos(om) - y1 * sin(om);
		Ear[1] = x1 * sin(om) * cos(inc) + y1 * cos(om) * cos(inc);
		Ear[2] = x1 * sin(om) * sin(inc) + y1 * cos(om) * sin(inc);
		Et[0] = Et[1] = 0;
		for (int i = 0; i < 3; i++) {
			Et[0] += Ear[i] * rad[i];
			Et[1] += Ear[i] * tang[i];
		}
		Et[0] += -Et0[0] - vt0[0] * (t - t0_par);
		Et[1] += -Et0[1] - vt0[1] * (t - t0_par);

		if (satellite > 0 && satellite <= nsat) {
			if (ndatasat[satellite - 1] > 2) {
				int left, right;
				if (t < tsat[satellite - 1][0]) {
					ic = 0;
				}
				else {
					if (t > tsat[satellite - 1][ndatasat[satellite - 1] - 1]) {
						ic = ndatasat[satellite - 1] - 2;
					}
					else {
						left = 0;
						right = ndatasat[satellite - 1] - 1;
						while (right - left > 1) {
							ic = (right + left) / 2;
							if (tsat[satellite - 1][ic] > t) {
								right = ic;
							}
							else {
								left = ic;
							}
						}
						ic = left;
					}
				}
				ty = t - tsat[satellite - 1][ic];
				for (int i = 0; i < 3; i++) {
					Spit = possat[satellite - 1][ic][i] * (1 - ty) + possat[satellite - 1][ic + 1][i] * ty;
					Et[0] += Spit * rad[i];
					Et[1] += Spit * tang[i];
				}
			}
		}
	}
}


#pragma endregion

#pragma region criticalcurves-caustics


//////////////////////////////
//////////////////////////////
////////Critical curves and caustics
//////////////////////////////
//////////////////////////////


_sols* VBMicrolensing::PlotCrit() {
	complex ej, y, z;
	int NPS = 200;
	_sols* CriticalCurves;
	_curve* Prov, * Prov2, * isso;
	_point* pisso;
	double SD, MD, CD;

	CriticalCurves = new _sols;
	for (int i = 0; i < 2 * n; i++) {
		Prov = new _curve;
		CriticalCurves->append(Prov);
	}

	for (int j = 0; j < NPcrit; j++) {
		ej = complex(cos(2 * j * M_PI / NPcrit), -sin(2 * j * M_PI / NPcrit));
		polycritcoefficients(ej);
		cmplx_roots_gen(zcr, coefs, 2 * n, true, true);
		if (j > 0) {
			Prov2 = new _curve();
			for (int i = 0; i < 2 * n; i++) {
				Prov2->append(zcr[i].re + s_offset->re, zcr[i].im + s_offset->im);
			}
			for (Prov = CriticalCurves->first; Prov; Prov = Prov->next) {
				Prov2->closest(Prov->last, &pisso);
				Prov2->drop(pisso);
				Prov->append(pisso);
			}
		}
		else {
			Prov = CriticalCurves->first;
			for (int i = 0; i < 2 * n; i++) {
				Prov->append(zcr[i].re + s_offset->re, zcr[i].im + s_offset->im);
				Prov = Prov->next;
			}
		}
	}

	Prov = CriticalCurves->first;
	while (Prov->next) {
		SD = *(Prov->first) - *(Prov->last);
		MD = 1.e100;
		isso = 0;
		for (Prov2 = Prov->next; Prov2; Prov2 = Prov2->next) {
			CD = *(Prov2->first) - *(Prov->last);
			if (CD < MD) {
				MD = CD;
				isso = Prov2;
			}
		}
		if (MD < SD) {
			CriticalCurves->drop(isso);
			Prov->join(isso);
		}
		else {
			Prov = Prov->next;
		}
	}

	// Caustics

	for (Prov = CriticalCurves->last; Prov; Prov = Prov->prev) {
		Prov2 = new _curve;
		for (_point* scanpoint = Prov->first; scanpoint; scanpoint = scanpoint->next) {
			y = z = complex(scanpoint->x1 - s_offset->re, scanpoint->x2 - s_offset->im);
			for (int i = 0; i < n; i++) {
				y = y - m[i] / conj(z - a[i]);
			}
			Prov2->append(y.re + s_offset->re, y.im + s_offset->im);
		}
		CriticalCurves->append(Prov2);
	}
	return CriticalCurves;
}

_sols* VBMicrolensing::PlotCrit(double a1, double q1) {
	complex  a, q, ej, zr[4], x1, x2;
	int NPS = 200;
	_sols* CriticalCurves;
	_curve* Prov, * Prov2, * isso;
	_point* pisso;
	double SD, MD, CD, centeroffset;

	a = complex(a1, 0.0);
	q = complex(q1, 0.0);
	centeroffset = a1 / 2.0 * (1.0 - q1) / (1.0 + q1);

	CriticalCurves = new _sols;
	for (int i = 0; i < 4; i++) {
		Prov = new _curve;
		CriticalCurves->append(Prov);
	}

	for (int j = 0; j < NPcrit; j++) {
		ej = complex(cos(2 * j * M_PI / NPcrit), -sin(2 * j * M_PI / NPcrit));
		complex  coefs[5] = { a * a / 16.0 * (4.0 - a * a * ej) * (1.0 + q),a * (q - 1.0),(q + 1.0) * (1.0 + a * a * ej / 2.0),0.0,-(1.0 + q) * ej };
		cmplx_roots_gen(zr, coefs, 4, true, true);
		if (j > 0) {
			Prov2 = new _curve();
			for (int i = 0; i < 4; i++) {
				Prov2->append(zr[i].re + centeroffset, zr[i].im);
			}
			for (Prov = CriticalCurves->first; Prov; Prov = Prov->next) {
				Prov2->closest(Prov->last, &pisso);
				Prov2->drop(pisso);
				Prov->append(pisso);
			}
		}
		else {
			Prov = CriticalCurves->first;
			for (int i = 0; i < 4; i++) {
				Prov->append(zr[i].re + centeroffset, zr[i].im);
				Prov = Prov->next;
			}
		}
	}

	Prov = CriticalCurves->first;
	while (Prov->next) {
		SD = *(Prov->first) - *(Prov->last);
		MD = 1.e100;
		isso = 0;
		for (Prov2 = Prov->next; Prov2; Prov2 = Prov2->next) {
			CD = *(Prov2->first) - *(Prov->last);
			if (CD < MD) {
				MD = CD;
				isso = Prov2;
			}
		}
		if (MD < SD) {
			CriticalCurves->drop(isso);
			Prov->join(isso);
		}
		else {
			Prov = Prov->next;
		}
	}

	// Caustics

	for (Prov = CriticalCurves->last; Prov; Prov = Prov->prev) {
		Prov2 = new _curve;
		for (_point* scanpoint = Prov->first; scanpoint; scanpoint = scanpoint->next) {
			x1 = complex(scanpoint->x1 - centeroffset, 0.0);
			x2 = complex(scanpoint->x2, 0.0);
			Prov2->append(real(_L1) + centeroffset, real(_L2));
		}
		CriticalCurves->append(Prov2);
	}
	return CriticalCurves;
}

#pragma endregion

#pragma region polynomials


//////////////////////////////
//////////////////////////////
// Polynomial methods
//////////////////////////////
//////////////////////////////


void VBMicrolensing::polyproduct(complex * __restrict p1, int n1, complex * __restrict p2, int n2, complex * __restrict pdest) {
	// polynomials are with increasing degree: p1[0]+p1[1]*z+p1[2]*z^2 + ...
	for (int i = 0; i <= n1 + n2; i++) pdest[i] = 0;
	for (int i = 0; i <= n1; i++) {
		for (int j = 0; j <= n2; j++) {
			pdest[i + j] = pdest[i + j] + p1[i] * p2[j];
		}
	}
}

void VBMicrolensing::copypol(complex * __restrict p1, int n1, complex * __restrict pdest) {
	for (int i = 0; i <= n1; i++) {

		pdest[i] = p1[i];
	}
}

void VBMicrolensing::change_n(int nn) {
	if (coefs) free(coefs);
	if (m) {
		free(m);
		free(a);
	}
	if (zr) {
		free(zr);
		free(zcr);
		free(good);
		free(Jacs);
		free(worst);
		free(pert);
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
	if (coefs_mp) {
		for (int i = 0; i < n; i++) {
			free(coefs_mp[i]);
			coefs_mp[i] = NULL;
		}
		free(coefs_mp);
		coefs_mp = NULL;
	}

	if (m_mp) {
		for (int i = 0; i < n; i++) {
			free(m_mp[i]);
			m_mp[i] = NULL;
			free(a_mp[i]);
			a_mp[i] = NULL;
		}
		free(m_mp);
		m_mp = NULL;
		free(a_mp);
		a_mp = NULL;
		free(q_sort);
		q_sort = NULL;
		free(s_sort);
		s_sort = NULL;
		free(y_mp);
		y_mp = NULL;
	}

	if (pmza_mp) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				free(pmza_mp[j][i]);
				pmza_mp[j][i] = NULL;
				free(pyaza_mp[j][i]);
				pyaza_mp[j][i] = NULL;
				free(ppmy_mp[j][i]);
				ppmy_mp[j][i] = NULL;
			}
		}
		for (int i = 0; i < n; i++) {
			free(pmza_mp[i]);
			pmza_mp[i] = NULL;
			free(pyaza_mp[i]);
			pyaza_mp[i] = NULL;
			free(ppmy_mp[i]);
			ppmy_mp[i] = NULL;
			free(ppy_mp[i]);
			ppy_mp[i] = NULL;
			free(pza_mp[i]);
			pza_mp[i] = NULL;
		}
		free(pmza_mp);
		pmza_mp = NULL;
		free(pyaza_mp);
		pyaza_mp = NULL;
		free(ppmy_mp);
		ppmy_mp = NULL;
		free(pza_mp);
		pza_mp = NULL;
		free(ppy_mp);
		ppy_mp = NULL;
	}

	if (zr_mp) {
		for (int j = 0; j < n; j++) {
			free(zr_mp[j]);
			zr_mp[j] = NULL;
		}
		free(zr_mp);
		zr_mp = NULL;
		free(nrootsmp_mp);
		nrootsmp_mp = NULL;
	}

	n = nn;

	n2 = n * n;
	nnm1 = n2 - n;
	nroots = 2 * n2 + 1;
	coefs = (complex*)malloc(sizeof(complex) * (nroots + 1));
	pmza = (complex**)malloc(sizeof(complex*) * n);
	pmza2 = (complex**)malloc(sizeof(complex*) * n);
	pyaza = (complex**)malloc(sizeof(complex*) * n);
	ppmy = (complex**)malloc(sizeof(complex*) * n);
	za = (complex**)malloc(sizeof(complex*) * n);
	za2 = (complex**)malloc(sizeof(complex*) * n);
	for (int i = 0; i < n; i++) {
		pmza[i] = (complex*)malloc(sizeof(complex) * n);
		pmza2[i] = (complex*)malloc(sizeof(complex) * (2 * n - 1));
		pyaza[i] = (complex*)malloc(sizeof(complex) * (n + 1));
		ppmy[i] = (complex*)malloc(sizeof(complex) * (nnm1 + 1));
		za[i] = (complex*)malloc(sizeof(complex) * (nroots));
		za2[i] = (complex*)malloc(sizeof(complex) * (nroots));
	}

	pza = (complex*)malloc(sizeof(complex) * (n + 1));
	pza2 = (complex*)malloc(sizeof(complex) * (2 * n + 1));
	pdum = (complex*)malloc(sizeof(complex) * (nroots + 1));
	ppy = (complex*)malloc(sizeof(complex) * (nroots + 1));

	zr = (complex*)malloc(sizeof(complex) * (nroots));
	zcr = (complex*)malloc(sizeof(complex) * (2 * n));
	for (int i = 0; i < nroots; i++) {
		zr[i] = 0;
	}
	for (int i = 0; i < 2 * n; i++) {
		zcr[i] = 0;
	}


	good = (double*)malloc(sizeof(double) * (nroots));
	Jacs = (double*)malloc(sizeof(double) * (nroots));
	J1 = (complex*)malloc(sizeof(complex) * (nroots));
	J1c = (complex*)malloc(sizeof(complex) * (nroots));
	zaltc = (complex*)malloc(sizeof(complex) * (nroots));
	worst = (int*)malloc(sizeof(int) * (nroots));
	pert = (complex*)malloc(sizeof(complex) * n);

	m = (double*)malloc(sizeof(double) * n);
	a = (complex*)malloc(sizeof(complex) * n);
	prodevs = (double*)malloc(sizeof(double) * n);
	devs = (complex*)malloc(sizeof(complex) * n);
	init = (complex*)malloc(sizeof(complex) * (n + 1));
	centralimages = (complex*)malloc(sizeof(complex) * (nnm1) / 2);
	errs = (double*)malloc(sizeof(double) * nroots);
	newseeds = (complex*)malloc(sizeof(complex) * (2 * nroots));
	grads = (complex*)malloc(sizeof(complex) * (nroots));
	S2s = (complex*)malloc(sizeof(complex) * (nroots));
	S3s = (complex*)malloc(sizeof(complex) * (nroots));
	S4s = (complex*)malloc(sizeof(complex) * (nroots));

	cprec = (_curve**)malloc(sizeof(_curve*) * nroots);
	cpres = (_curve**)malloc(sizeof(_curve*) * (nroots));
	cfoll = (_curve**)malloc(sizeof(_curve*) * (nroots));
	A = (double**)malloc(sizeof(double*) * (nroots));
	for (int i = 0; i < nroots; i++) {
		A[i] = (double*)malloc(sizeof(double) * (nroots));
	}

}

void VBMicrolensing::change_n_mp(int nn) {
	if (coefs) free(coefs);
	if (m) {
		free(m);
		free(a);
	}
	if (zr) {
		free(zr);
		free(zcr);
		free(good);
		free(Jacs);
		free(worst);
		free(pert);
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

	n = nn;

	n2 = n * n;
	nnm1 = n2 - n;
	nroots = 2 * n2 + 1;
	nrootsmp_mp = (int*)malloc(sizeof(int) * n);

	q_sort = (double*)malloc(sizeof(double) * n);
	s_sort = (complex*)malloc(sizeof(complex) * n);
	y_mp = (complex*)malloc(sizeof(complex) * n);

	m_mp = (double**)malloc(sizeof(double*) * n);
	a_mp = (complex**)malloc(sizeof(complex*) * n);
	for (int i = 0; i < n; i++) {
		m_mp[i] = (double*)malloc(sizeof(double) * n);
		a_mp[i] = (complex*)malloc(sizeof(complex) * n);
	}

	coefs_mp = (complex**)malloc(sizeof(complex*) * n);
	pza_mp = (complex**)malloc(sizeof(complex*) * n);
	ppy_mp = (complex**)malloc(sizeof(complex*) * n);
	for (int j = 0; j < n; j++) {
		coefs_mp[j] = (complex*)malloc(sizeof(complex) * (nroots + 1));
		pza_mp[j] = (complex*)malloc(sizeof(complex) * (n + 1));
		ppy_mp[j] = (complex*)malloc(sizeof(complex) * (nroots + 1));
	}

	pmza_mp = (complex***)malloc(sizeof(complex**) * n);
	pyaza_mp = (complex***)malloc(sizeof(complex**) * n);
	ppmy_mp = (complex***)malloc(sizeof(complex**) * n);
	for (int j = 0; j < n; j++) {
		pmza_mp[j] = (complex**)malloc(sizeof(complex*) * n);
		pyaza_mp[j] = (complex**)malloc(sizeof(complex*) * n);
		ppmy_mp[j] = (complex**)malloc(sizeof(complex*) * n);
		for (int i = 0; i < n; i++) {
			pmza_mp[j][i] = (complex*)malloc(sizeof(complex) * n);
			pyaza_mp[j][i] = (complex*)malloc(sizeof(complex) * (n + 1));
			ppmy_mp[j][i] = (complex*)malloc(sizeof(complex) * (nnm1 + 1));
		}
	}
	dist_mp = (double*)malloc(sizeof(double) * n);

	zr_mp = (complex**)malloc(sizeof(complex*) * n);
	for (int j = 0; j < n; j++) {
		zr_mp[j] = (complex*)malloc(sizeof(complex) * nroots);
		for (int i = 0; i < n; i++) {
			zr_mp[j][i] = 0;
		}
	}

	coefs = (complex*)malloc(sizeof(complex) * (nroots + 1));
	pmza = (complex**)malloc(sizeof(complex*) * n);
	pmza2 = (complex**)malloc(sizeof(complex*) * n);
	pyaza = (complex**)malloc(sizeof(complex*) * n);
	ppmy = (complex**)malloc(sizeof(complex*) * n);
	za = (complex**)malloc(sizeof(complex*) * n);
	za2 = (complex**)malloc(sizeof(complex*) * n);
	for (int i = 0; i < n; i++) {
		pmza[i] = (complex*)malloc(sizeof(complex) * n);
		pmza2[i] = (complex*)malloc(sizeof(complex) * (2 * n - 1));
		pyaza[i] = (complex*)malloc(sizeof(complex) * (n + 1));
		ppmy[i] = (complex*)malloc(sizeof(complex) * (nnm1 + 1));
		za[i] = (complex*)malloc(sizeof(complex) * (nroots));
		za2[i] = (complex*)malloc(sizeof(complex) * (nroots));
	}

	pza = (complex*)malloc(sizeof(complex) * (n + 1));
	pza2 = (complex*)malloc(sizeof(complex) * (2 * n + 1));
	pdum = (complex*)malloc(sizeof(complex) * (nroots + 1));
	ppy = (complex*)malloc(sizeof(complex) * (nroots + 1));

	zr = (complex*)malloc(sizeof(complex) * (nroots));
	zcr = (complex*)malloc(sizeof(complex) * (2 * n));
	for (int i = 0; i < nroots; i++) {
		zr[i] = 0;
	}
	for (int i = 0; i < 2 * n; i++) {
		zcr[i] = 0;
	}


	good = (double*)malloc(sizeof(double) * (nroots));
	Jacs = (double*)malloc(sizeof(double) * (nroots));
	J1 = (complex*)malloc(sizeof(complex) * (nroots));
	J1c = (complex*)malloc(sizeof(complex) * (nroots));
	zaltc = (complex*)malloc(sizeof(complex) * (nroots));
	worst = (int*)malloc(sizeof(int) * (nroots));
	pert = (complex*)malloc(sizeof(complex) * n);

	m = (double*)malloc(sizeof(double) * n);
	a = (complex*)malloc(sizeof(complex) * n);
	prodevs = (double*)malloc(sizeof(double) * n);
	devs = (complex*)malloc(sizeof(complex) * n);
	init = (complex*)malloc(sizeof(complex) * (n + 1));
	centralimages = (complex*)malloc(sizeof(complex) * (nnm1) / 2);
	errs = (double*)malloc(sizeof(double) * nroots);
	newseeds = (complex*)malloc(sizeof(complex) * (2 * nroots));
	grads = (complex*)malloc(sizeof(complex) * (nroots));
	S2s = (complex*)malloc(sizeof(complex) * (nroots));
	S3s = (complex*)malloc(sizeof(complex) * (nroots));
	S4s = (complex*)malloc(sizeof(complex) * (nroots));

	cprec = (_curve**)malloc(sizeof(_curve*) * nroots);
	cpres = (_curve**)malloc(sizeof(_curve*) * (nroots));
	cfoll = (_curve**)malloc(sizeof(_curve*) * (nroots));
	A = (double**)malloc(sizeof(double*) * (nroots));
	for (int i = 0; i < nroots; i++) {
		A[i] = (double*)malloc(sizeof(double) * (nroots));
	}
}

void VBMicrolensing::polycritcoefficients(complex eiphi) {
	static int dg;
	static complex pbin[3];

	for (int i = 0; i < n; i++) {
		pmza2[i][0] = m[i];
	}
	pza2[0] = pbin[2] = 1.0;
	for (int i = 0; i < n; i++) {
		pbin[0] = a[i] * a[i];
		pbin[1] = -2 * a[i];
		for (int j = 0; j < n; j++) {
			if (i == j) {
				copypol(pza2, i * 2, pdum);
				polyproduct(pdum, i * 2, pbin, 2, pza2);
			}
			else {
				dg = (j > i) ? i * 2 : (i - 1) * 2;
				copypol(pmza2[j], dg, pdum);
				polyproduct(pdum, dg, pbin, 2, pmza2[j]);
			}
		}
	}
	copypol(pza2, 2 * n, coefs);
	for (int j = 0; j <= 2 * n - 2; j++) {
		coefs[j] = coefs[j] * eiphi;
		for (int i = 0; i < n; i++) {
			coefs[j] = coefs[j] - pmza2[i][j];
		}
	}
	coefs[2 * n - 1] = coefs[2 * n - 1] * eiphi;
	coefs[2 * n] = coefs[2 * n] * eiphi;
}

// polycoefficients
// Calculates the polynomial coefficients necessary for the solution of the lens equation
// VBML::coefs polynomial of coefficients of degree n^2+1.

void VBMicrolensing::polycoefficients() {

	static int dg;
	static complex pbin[2], lam;


	for (int k = 0; k < n; k++) {
		pyaza[0][k] = 0;
		for (int j = 0; j < n; j++) {
			pyaza[0][k] = pyaza[0][k] + pmza[j][k];
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		lam = conj(y - a[i]);
		for (int k = 0; k < n; k++) {
			pyaza[i][k] = lam * pza[k] + pyaza[0][k];
		}
		pyaza[i][n] = lam * pza[n];
	}

	for (int i = 0; i < n; i++) {
		ppmy[i][0] = m[i];
	}
	ppy[0] = 1.0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				dg = i * n;
				copypol(ppy, dg, pdum);
				polyproduct(pdum, dg, pyaza[i], n, ppy);
			}
			else {
				dg = (j > i) ? i * n : (i - 1) * n;
				copypol(ppmy[j], dg, pdum);
				polyproduct(pdum, dg, pyaza[i], n, ppmy[j]);
			}
		}
	}

	for (int i = 0; i < nnm1 + 1; i++) {
		pdum[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < nnm1 + 1; j++) {
			pdum[j] = pdum[j] + ppmy[i][j];
		}
	}
	polyproduct(pdum, nnm1, pza, n, coefs);

	pbin[0] = y;
	pbin[1] = -1.0;
	copypol(ppy, n2, pdum);
	polyproduct(pdum, n2, pbin, 1, ppy);

	for (int i = 0; i < n2 + 1; i++) {
		coefs[i] = coefs[i] + ppy[i];
	}
	coefs[n2 + 1] = ppy[n2 + 1];
}

void VBMicrolensing::polycoefficients_multipoly() {

	static int dg;
	static complex pbin[2], lam;

	for (int l = 0; l < n; l++) {
		for (int i = 0; i < nnm1 + 1; i++) {
			pdum[i] = 0;
		}
		for (int k = 0; k < n; k++) {
			pyaza_mp[l][0][k] = 0;
			for (int j = 0; j < n; j++) {
				pyaza_mp[l][0][k] = pyaza_mp[l][0][k] + pmza_mp[l][j][k];
			}
		}
		for (int i = n - 1; i >= 0; i--) {
			lam = conj(y_mp[l] - a_mp[l][i]);
			for (int k = 0; k < n; k++) {
				pyaza_mp[l][i][k] = lam * pza_mp[l][k] + pyaza_mp[l][0][k];
			}
			pyaza_mp[l][i][n] = lam * pza_mp[l][n];
		}

		for (int i = 0; i < n; i++) {
			ppmy_mp[l][i][0] = m_mp[l][i];
		}
		ppy_mp[l][0] = 1.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					dg = i * n;
					copypol(ppy_mp[l], dg, pdum);
					polyproduct(pdum, dg, pyaza_mp[l][i], n, ppy_mp[l]);
				}
				else {
					dg = (j > i) ? i * n : (i - 1) * n;
					copypol(ppmy_mp[l][j], dg, pdum);
					polyproduct(pdum, dg, pyaza_mp[l][i], n, ppmy_mp[l][j]);
				}
			}
		}

		for (int i = 0; i < nnm1 + 1; i++) {
			pdum[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < nnm1 + 1; j++) {
				pdum[j] = pdum[j] + ppmy_mp[l][i][j];
			}
		}
		polyproduct(pdum, nnm1, pza_mp[l], n, coefs_mp[l]);

		pbin[0] = y_mp[l];
		pbin[1] = -1.0;
		copypol(ppy_mp[l], n2, pdum);
		polyproduct(pdum, n2, pbin, 1, ppy_mp[l]);

		for (int i = 0; i < n2 + 1; i++) {
			coefs_mp[l][i] = coefs_mp[l][i] + ppy_mp[l][i];
		}
		coefs_mp[l][n2 + 1] = ppy_mp[l][n2 + 1];
	}

}


#pragma endregion

#pragma region Skowron-Gould


//////////////////////////////
//////////////////////////////
////////_Skowron & Gould functions, translated by Tyler M. Heintz and Ava R. Hoag
//////////////////////////////
//////////////////////////////
// See copyright notice for these functions


void VBMicrolensing::cmplx_roots_gen(complex * roots, complex * poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points) {
	//roots - array which will hold all roots that had been found.
	//If the flag 'use_roots_as_starting_points' is set to
	//.true., then instead of point(0, 0) we use value from
	//this array as starting point for cmplx_laguerre

	//poly - is an array of polynomial cooefs, length = degree + 1,
	//poly[0] x ^ 0 + poly[1] x ^ 1 + poly[2] x ^ 2 + ...

	//degree - degree of the polynomial and size of 'roots' array

	//polish_roots_after - after all roots have been found by dividing
	//original polynomial by each root found,
	//you can opt in to polish all roots using full
	//polynomial

	//use_roots_as_starting_points - usually we start Laguerre's 
	//method from point(0, 0), but you can decide to use the
	//values of 'roots' array as starting point for each new
	//root that is searched for.This is useful if you have
	//very rough idea where some of the roots can be.
	//

	static complex poly2[MAXM];
	static int i, j, n, iter;
	static bool success;
	static complex coef, prev;

	if (!use_roots_as_starting_points) {
		for (int jj = 0; jj < degree; jj++) {
			roots[jj] = complex(0, 0);
		}
	}

	for (j = 0; j <= degree; j++) poly2[j] = poly[j];

	// Don't do Laguerre's for small degree polynomials
	if (degree <= 1) {
		if (degree == 1) roots[0] = -poly[0] / poly[1];
		return;
	}

	for (n = degree; n >= 3; n--) {
		cmplx_laguerre2newton(poly2, n, &roots[n - 1], iter, success, 2);
		if (!success) {
			roots[n - 1] = complex(0, 0);
			cmplx_laguerre(poly2, n, &roots[n - 1], iter, success);
		}

		// Divide by root
		coef = poly2[n];
		for (i = n - 1; i >= 0; i--) {
			prev = poly2[i];
			poly2[i] = coef;
			coef = prev + roots[n - 1] * coef;
		}
	}

	//Find the to last 2 roots
	solve_quadratic_eq(roots[1], roots[0], poly2);
	//cmplx_laguerre2newton(poly2, 2, &roots[1], iter, success, 2);
	//if (!success) {
	//	solve_quadratic_eq(roots[1], roots[0], poly2);
	//}
	//else {
	//	roots[0] = -(roots[1] + poly2[1] / poly2[2]); // Viete's Formula for the last root
	//}



	if (polish_roots_after) {
		for (n = 0; n < degree; n++) {
			cmplx_newton_spec(poly, degree, &roots[n], iter, success); // Polish roots with full polynomial
		}
	}

	return;
}

void VBMicrolensing::cmplx_roots_multigen(complex * roots, complex * *poly, int degree, bool polish_roots_after, bool use_roots_as_starting_points) {

	static complex poly2[MAXM];
	static int l, j, i, k, nl, ind, degreenew, croots, n;
	static double dif0, br;
	static bool success;
	static complex coef, prev, przr;


	nl = sqrt(degree - 1);
	for (l = 0; l < nl; l++) nrootsmp_mp[l] = 0;
	for (l = 0; l < nl; l++) {
		for (i = 0; i < degree; i++) {
			zr_mp[l][i] = complex(0., 0.);
		}
	}
	//Cycle reference systems
	for (l = 0; l < nl; l++) {

		br = false;
		//copy poly coefs
		for (j = 0; j <= degree; j++) poly2[j] = poly[l][j];
		//Don't do Lagierre's for small degree polybnomials
		if (l != nl - 1) {
			if (degree <= 1) {
				if (degree == 1) zr_mp[l][0] = -poly[l][0] / poly[l][1];
				nrootsmp_mp[l] = 1;
				break;
			}
			//Do Laguerre for degree >=3
			for (n = degree; n >= 3; n--) {
				cmplx_laguerre2newton(poly2, n, &zr_mp[l][n - 1], iter, success, 2);
				if (!success) {
					zr_mp[l][n - 1] = complex(0, 0);
					cmplx_laguerre(poly2, n, &zr_mp[l][n - 1], iter, success);
				}
				nrootsmp_mp[l]++;
				//distance check
				dif0 = abs2(zr_mp[l][n - 1]);
				for (i = 1; i < nl; i++) {
					if (abs2(zr_mp[l][n - 1] - a_mp[l][i]) < dif0) {
						dist_mp[l] = abs2(zr_mp[l][n - 1] - a_mp[l][i]);
						zr_mp[l][n - 1] = complex(0, 0);
						br = true;
						nrootsmp_mp[l]--;
						break;
					}
				}
				if (br) break;
				//Divide by root
				coef = poly2[n];
				for (i = n - 1; i >= 0; i--) {
					prev = poly2[i];
					poly2[i] = coef;
					coef = prev + zr_mp[l][n - 1] * coef;
				}
			}
			if (br) continue;
			//find the last 2 roots
			solve_quadratic_eq(zr_mp[l][1], zr_mp[l][0], poly2);
			nrootsmp_mp[l] += 2;
			for (i = 1; i < nl; i++) {
				if (abs2(zr_mp[l][1] - a_mp[l][i]) < abs2(zr_mp[l][1])) {
					zr_mp[l][1] = zr_mp[l][0];
					zr_mp[l][0] = complex(0, 0);
					dist_mp[l] = abs2(zr_mp[l][1] - a_mp[l][i]);
					nrootsmp_mp[l]--;
					break;
				}
			}
			k = degree - nrootsmp_mp[l];
			for (i = 1; i < nl; i++) {
				if (abs2(zr_mp[l][k] - a_mp[l][i]) < abs2(zr_mp[l][k])) {
					zr_mp[l][k] = complex(0, 0);
					nrootsmp_mp[l]--;
					break;
				}
			}
		}

		//LAST lens
		if (l == nl - 1) {
			//Set previous roots
			ind = 0;
			for (int ll = 0; ll < nl - 1; ll++) {
				for (int i = ind; i < ind + nrootsmp_mp[ll]; i++) {
					zr_mp[l][degree - i - 1] = zr_mp[ll][degree - 1 - i + ind] + s_sort[ll] - s_sort[l];
				}
				ind += nrootsmp_mp[ll];
			}


			//divide by previous roots

			degreenew = degree;
			for (int i = 0; i < nl - 1; i++) {
				degreenew -= nrootsmp_mp[i];
			}

			for (int n = degree; n > degreenew; n--) {
				coef = poly2[n];
				for (i = n - 1; i >= 0; i--) {
					prev = poly2[i];
					poly2[i] = coef;
					coef = prev + zr_mp[l][n - 1] * coef;
				}
			}

			if (degreenew <= 1) {
				if (degreenew == 1) zr_mp[l][0] = -poly2[0] / poly2[1];
				nrootsmp_mp[l] = 1;

				break;
			}

			for (n = degreenew; n >= 3; n--) {
				cmplx_laguerre2newton(poly2, n, &zr_mp[l][n - 1], iter, success, 2);
				if (!success) {
					zr_mp[l][n - 1] = complex(0, 0);
					cmplx_laguerre(poly2, n, &zr_mp[l][n - 1], iter, success);
				}
				nrootsmp_mp[l] += 1;

				// Divide by root
				coef = poly2[n];
				for (i = n - 1; i >= 0; i--) {
					prev = poly2[i];
					poly2[i] = coef;
					coef = prev + zr_mp[l][n - 1] * coef;
				}
			}
			solve_quadratic_eq(zr_mp[l][1], zr_mp[l][0], poly2);
			nrootsmp_mp[l] += 2;

		}
	}

	ind = degree - 1;
	for (l = 0; l < nl - 1; l++) {
		for (i = 0; i < nrootsmp_mp[l]; i++) {
			roots[ind] = zr_mp[l][degree - 1 - i] + s_sort[l] - s_sort[0];
			ind--;
		}
	}
	for (i = 0; i < nrootsmp_mp[nl - 1]; i++) {
		roots[ind] = zr_mp[nl - 1][i] + s_sort[nl - 1] - s_sort[0];
		ind--;
	}

	return;
}

void VBMicrolensing::solve_quadratic_eq(complex & x0, complex & x1, complex * poly) {
	static complex a, b, c, b2, delta;
	a = poly[2];
	b = poly[1];
	c = poly[0];
	b2 = b * b;
	delta = sqrt(b2 - 4 * a * c);
	if (real(conj(b) * delta) >= 0) {
		x0 = -0.5 * (b + delta);
	}
	else {
		x0 = -0.5 * (b - delta);
	}
	if (x0 == complex(0., 0.)) {
		x1 = complex(0., 0.);
	}
	else { //Viete's formula
		x1 = c / x0;
		x0 = x0 / a;
	}
	return;

}

void VBMicrolensing::solve_cubic_eq(complex & x0, complex & x1, complex & x2, complex * poly) {
	//Cubic equation solver for comples polynomial (degree=3)
	//http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
	// poly is an array of polynomial cooefs, length = degree+1, poly[0] is constant
	//	0				1				2			3
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + poly[3] x^3
	static complex zeta = complex(-0.5, 0.8660254037844386);
	static complex zeta2 = complex(-0.5, -0.8660254037844386);
	static double third = 0.3333333333333333;
	static complex s0, s1, s2;
	static complex E1; //x0+x1+x2
	static complex E2; //x0*x1+x1*x2+x2*x0
	static complex E3; //x0*x1*x2
	static complex A, B, a_1, E12, delta, A2;

	static complex val, x;
	a_1 = 1 / poly[3];
	E1 = -poly[2] * a_1;
	E2 = poly[1] * a_1;
	E3 = -poly[0] * a_1;

	s0 = E1;
	E12 = E1 * E1;
	A = 2.0 * E1 * E12 - 9.0 * E1 * E2 + 27.0 * E3;
	B = E12 - 3.0 * E2;
	//quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
	A2 = A * A;
	delta = sqrt(A2 - 4.0 * (B * B * B));
	if (real(conj(A) * delta) >= 0.0) { // scalar product to decide the sign yielding bigger magnitude
		s1 = cbrt(0.5 * (A + delta));
	}
	else
	{
		s1 = cbrt(0.5 * (A - delta));
	}
	if (s1.re == 0.0 && s1.im == 0.0) {
		s2 = complex(0, 0);
	}
	else {
		s2 = B / s1;
	}

	x0 = third * (s0 + s1 + s2);
	x1 = third * (s0 + s1 * zeta2 + s2 * zeta);
	x2 = third * (s0 + s1 * zeta + s2 * zeta2);

	return;

}

void VBMicrolensing::cmplx_laguerre(complex * poly, int degree, complex * root, int& iter, bool& success) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method. In every loop it calculates simplified 
	//Adams' stopping criterion for the value of the polynomial.
	//
	//Uses 'root' value as a starting point(!!!!!)
	//Remember to initialize 'root' to some initial guess or to
	//point(0, 0) if you have no prior knowledge.
	//
	//poly - is an array of polynomial cooefs
	//
	//length = degree + 1, poly(1) is constant
	//	1              2				3
	//poly(1) x ^ 0 + poly(2) x ^ 1 + poly(3) x ^ 2 + ...
	//
	//degree - a degree of the polynomial
	//
	//root - input: guess for the value of a root
	//output : a root of the polynomial
	//iter - number of iterations performed(the number of polynomial
	//evaluations and stopping criterion evaluation)
	//
	//success - is false if routine reaches maximum number of iterations
	//
	//For a summary of the method go to :
	//http://en.wikipedia.org/wiki/Laguerre's_method
	//
	static int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	static double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297,
		0.91577881, 0.25921289, 0.50487203,
		0.08177045, 0.13653241, 0.306162,
		0.37794326, 0.04618805, 0.75132137 }; // some random numbers

	static double faq; //jump length
	static double FRAC_ERR = 2.0e-15; //Fractional Error for double precision
	static complex p, dp, d2p_half; //value of polynomial, 1st derivative, and 2nd derivative
	static int i, j, k;
	static bool good_to_go;
	static complex denom, denom_sqrt, dx, newroot;
	static double ek, absroot, abs2p;
	static complex fac_newton, fac_extra, F_half, c_one_nth;
	static double one_nth, n_1_nth, two_n_div_n_1;
	static complex c_one = complex(1, 0);
	static complex zero = complex(0, 0);
	static double stopping_crit2;

	//--------------------------------------------------------------------------------------------

	//EXTREME FAILSAFE! not usually needed but kept here just to be on the safe side. Takes care of first coefficient being 0
	if (false) {
		if (degree < 0) {
			printf("Error: cmplx_laguerre: degree<0");
			return;
		}
		if (poly[degree] == complex(0, 0)) {
			if (degree == 0) return;
			cmplx_laguerre(poly, degree - 1, root, iter, success);
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false; // we just checked if poly[0] is zero and it isnt
				printf("Warning: cmplx_laguerre: degree = 0 and poly[0] does not equal zero, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	} // End of EXTREME failsafe

	good_to_go = false;
	one_nth = 1.0 / degree;
	n_1_nth = (degree - 1.0) * one_nth;
	two_n_div_n_1 = 2.0 / n_1_nth;
	c_one_nth = complex(one_nth, 0.0);
	for (i = 1; i <= MAXIT; i++) {
		ek = abs(poly[degree]); // Preparing stopping criterion
		absroot = abs(*root);
		// Calculate the values of polynomial and its first and second derivatives
		p = poly[degree];
		dp = zero;
		d2p_half = zero;
		for (k = degree - 1; k >= 0; k--) {
			d2p_half = dp + d2p_half * (*root);
			dp = p + dp * *root;
			p = poly[k] + p * (*root); // b_k
			//Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
			//Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
			//ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
			//Eq 8.
			ek = absroot * ek + abs(p);
		}
		iter += 1;

		abs2p = real(conj(p) * p);
		if (abs2p == 0) return;
		stopping_crit2 = pow(FRAC_ERR * ek, 2.0);
		if (abs2p < stopping_crit2) {
			//(simplified a little Eq. 10 of Adams 1967)
			//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2) {
				return; // we are at a good place!
			}
			else {
				good_to_go = true;
			}
		}
		else {
			good_to_go = false;
		}

		faq = 1.0;
		denom = zero;
		if (dp != zero) {
			fac_newton = p / dp;
			fac_extra = d2p_half / dp;
			F_half = fac_newton * fac_extra;
			denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half);

			//NEXT LINE PROBABLY CAN BE COMMENTED OUT. Check if compiler outputs positive real
			if (real(denom_sqrt) >= 0.0) {
				denom = c_one_nth + n_1_nth * denom_sqrt;
			}
			else {
				denom = c_one_nth - n_1_nth * denom_sqrt;
			}
		}

		if (denom == 0) {
			dx = (absroot + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = fac_newton / denom;
		}


		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes so return
		if (good_to_go) {
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { //decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq * dx; // do jump of semi-random length
		}
		*root = newroot;
	}
	success = false; // too many iterations here
	return;
}

void VBMicrolensing::cmplx_newton_spec(complex * poly, int degree, complex * root, int& iter, bool& success) {
	//Subroutine finds one root of a complex polynomial
	//Newton's method. It calculates simplified Adams' stopping 
	//criterion for the value of the polynomial once per 10 iterations (!),
	//after initial iteration. This is done to speed up calculations
	//when polishing roots that are known preety well, and stopping
	// criterion does significantly change in their neighborhood.

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess.
	//Do not initilize 'root' to point (0,0) if the polynomial 
	//coefficients are strictly real, because it will make going 
	//to imaginary roots impossible.

	// poly - is an array of polynomial cooefs
	//	length = degree+1, poly(1) is constant 
	//0					1				2
	//poly[0] x^0 + poly[1] x^1 + poly[2] x^2 + ...
	//degree - a degree of the polynomial
	// root - input: guess for the value of a root
	//		  output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial evaluations)
	//success - is false if routine reaches maximum number of iterations

	//For a summary of the method go to: 
	//http://en.wikipedia.org/wiki/Newton's_method

	static int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	static double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers
	static double faq; //jump length
	static double FRAC_ERR = 2e-15;
	static complex p; //value of polynomial
	static complex dp; //value of 1st derivative
	static int i, k;
	static bool good_to_go;
	static complex dx, newroot;
	static double ek, absroot, abs2p;
	static complex zero = complex(0, 0);
	static double stopping_crit2;

	iter = 0;
	success = true;

	//the next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version
	if (false) { //change false to true if you would like to use caustion about haveing first coefficient == 0
		if (degree < 0) {
			printf("Error: cmplx_newton_spec: degree<0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_newton_spec(poly, degree, root, iter, success);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {
				success = false;
				printf("Warning: cmplx_newton_spec: degree=0 and poly[0]!=0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME Failsafe
	good_to_go = false;

	stopping_crit2 = 0.0; //value not important, will be initialized anyway on the first loop
	for (i = 1; i <= MAXIT; i++) {
		faq = 1.0;
		//prepare stopping criterion
		//calculate value of polynomial and its first two derivatives
		p = poly[degree];
		dp = zero;
		if (i % 10 == 1) { //calculate stopping criterion every tenth iteration
			ek = abs(poly[degree]);
			absroot = abs(*root);
			for (k = degree - 1; k >= 0; k--) {
				dp = p + dp * (*root);
				p = poly[k] + p * (*root); //b_k
				//Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
				//Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
				//ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
				//Eq. 8
				ek = absroot * ek + abs(p);
			}
			stopping_crit2 = pow(FRAC_ERR * ek, 2);
		}
		else { // calculate just the value and derivative
			for (k = degree - 1; k >= 0; k--) { //Horner Scheme, see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
				dp = p + dp * (*root);
				p = poly[k] + p * (*root);
			}
		}

		iter = iter + 1;

		abs2p = real(conj(p) * p);
		if (abs2p == 0.0) return;
		if (abs2p < stopping_crit2) { //simplified a little Eq. 10 of Adams 1967
			if (dp == zero) return; //if we have problem with zero, but we are close to the root, just accept
			//do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2) return; //return immediatley because we are at very good place
			else {
				good_to_go = true; //do one iteration more
			}
		}

		else {
			good_to_go = false; //reset if we are outside the zone of the root
		}
		if (dp == zero) {
			//problem with zero
			dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI));
		}
		else {
			dx = p / dp; // Newton method, see http://en.wikipedia.org/wiki/Newton's_method
		}
		newroot = *root - dx;
		if (newroot == *root) return; //nothing changes -> return
		if (good_to_go) {//this was jump already after stopping criterion was met
			*root = newroot;
			return;
		}
		if (i % FRAC_JUMP_EVERY == 0) { // decide whether to do a jump of modified length (to break cycles)
			faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
			newroot = *root - faq * dx;
		}
		*root = newroot;
	}
	success = false;
	return;
	//too many iterations here
}

void VBMicrolensing::cmplx_laguerre2newton(complex * poly, int degree, complex * root, int& iter, bool& success, int starting_mode) {
	//Subroutine finds one root of a complex polynomial using
	//Laguerre's method, Second-order General method and Newton's
	//method - depending on the value of function F, which is a 
	//combination of second derivative, first derivative and
	//value of polynomial [F=-(p"*p)/(p'p')].

	//Subroutine has 3 modes of operation. It starts with mode=2
	//which is the Laguerre's method, and continues until F
	//becames F<0.50, at which point, it switches to mode=1,
	//i.e., SG method (see paper). While in the first two
	//modes, routine calculates stopping criterion once per every
	//iteration. Switch to the last mode, Newton's method, (mode=0)
	//happens when becomes F<0.05. In this mode, routine calculates
	//stopping criterion only once, at the beginning, under an
	//assumption that we are already very close to the root.
	//If there are more than 10 iterations in Newton's mode,
	//it means that in fact we were far from the root, and
	//routine goes back to Laguerre's method (mode=2).

	//Uses 'root' value as a starting point (!!!!!)
	//Remember to initialize 'root' to some initial guess or to 
	//point (0,0) if you have no prior knowledge.

	//poly - is an array of polynomial cooefs
	//	0					1				2
	//	poly[0] x^0 + poly[1] x^1 + poly[2] x^2
	//degree - a degree of the polynomial
	//root - input: guess for the value of a root
	//		output: a root of the polynomial
	//iter - number of iterations performed (the number of polynomial
	//		 evaluations and stopping criterion evaluation)
	//success - is false if routine reaches maximum number of iterations
	//starting_mode - this should be by default = 2. However if you  
	//				  choose to start with SG method put 1 instead.
	//				  Zero will cause the routine to
	//				  start with Newton for first 10 iterations, and
	//				  then go back to mode 2.

	//For a summary of the method see the paper: Skowron & Gould (2012)

	static int FRAC_JUMP_EVERY = 10;
	const int FRAC_JUMP_LEN = 10;
	static double FRAC_JUMPS[FRAC_JUMP_LEN] = { 0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045, 0.13653241, 0.306162, 0.37794326, 0.04618805, 0.75132137 }; //some random numbers

	static double faq; //jump length
	static double FRAC_ERR = 2.0e-15;

	static complex p; //value of polynomial
	static complex dp; //value of 1st derivative
	static complex d2p_half; //value of 2nd derivative
	static int i, j, k;
	static bool good_to_go;
	//complex G, H, G2;
	static complex denom, denom_sqrt, dx, newroot;
	static double ek, absroot, abs2p, abs2_F_half;
	static complex fac_netwon, fac_extra, F_half, c_one_nth;
	static double one_nth, n_1_nth, two_n_div_n_1;
	static int mode;
	static complex c_one = complex(1, 0);
	static complex zero = complex(0, 0);
	static double stopping_crit2;

	iter = 0;
	success = true;
	stopping_crit2 = 0; //value not important, will be initialized anyway on the first loop

	//next if block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
	if (false) {//change false to true if you would like to use caution about having first coefficent == 0
		if (degree < 0) {
			printf("Error: cmplx_laguerre2newton: degree < 0");
			return;
		}
		if (poly[degree] == zero) {
			if (degree == 0) return;
			cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode);
			return;
		}
		if (degree <= 1) {
			if (degree == 0) {//// we know from previous check that poly[0] not equal zero
				success = false;
				printf("Warning: cmplx_laguerre2newton: degree = 0 and poly[0] = 0, no roots");
				return;
			}
			else {
				*root = -poly[0] / poly[1];
				return;
			}
		}
	}
	//end EXTREME failsafe

	j = 1;
	good_to_go = false;

	mode = starting_mode; // mode = 2 full laguerre, mode = 1 SG, mode = 0 newton

	for (;;) { //infinite loop, just to be able to come back from newton, if more than 10 iteration there

		////////////
		///mode 2///
		////////////

		if (mode >= 2) {//Laguerre's method
			one_nth = 1.0 / (degree); ///
			n_1_nth = (degree - 1) * one_nth; ////
			two_n_div_n_1 = 2.0 / n_1_nth;
			c_one_nth = complex(one_nth, 0.0);

			for (i = 1; i <= MAXIT; i++) {
				faq = 1.0;

				//prepare stoping criterion
				ek = abs(poly[degree]);
				absroot = abs(*root);
				//calculate value of polynomial and its first two derivative
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
					d2p_half = dp + d2p_half * (*root);
					dp = p + dp * (*root);
					p = poly[k - 1] + p * (*root); // b_k
					//Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
					//Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
					//ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
					//Eq 8.
					ek = absroot * ek + abs(p);
				}
				abs2p = real(conj(p) * p); // abs(p)
				iter = iter + 1;
				if (abs2p == 0) return;

				stopping_crit2 = pow(FRAC_ERR * ek, 2);
				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
					//do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01 * stopping_crit2) return; // ten times better than stopping criterion
					//return immediately, because we are at very good place
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}

				denom = zero;
				if (dp != zero) {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.0625) {//F<0.50, F/2<0.25
						//go to SG method
						if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.02
							mode = 0; //go to Newton's
						}
						else {
							mode = 1; //go to SG
						}
					}

					denom_sqrt = sqrt(c_one - two_n_div_n_1 * F_half);

					//NEXT LINE PROBABLY CAN BE COMMENTED OUT 
					if (real(denom_sqrt) > 0.0) {
						//real part of a square root is positive for probably all compilers. You can y
						//test this on your compiler and if so, you can omit this check
						denom = c_one_nth + n_1_nth * denom_sqrt;
					}
					else {
						denom = c_one_nth - n_1_nth * denom_sqrt;
					}
				}
				if (denom == zero) {//test if demoninators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) + expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					dx = fac_netwon / denom;
				}
				newroot = *root - dx;
				if (newroot == *root) return; // nothing changes -> return
				if (good_to_go) {//this was jump already after stopping criterion was met
					*root = newroot;
					return;
				}
				if (mode != 2) {
					*root = newroot;
					j = i + 1; //remember iteration index
					break; //go to Newton's or SG
				}
				if ((i % FRAC_JUMP_EVERY) == 0) {//decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[((i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN)];
					newroot = *root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
				}
				*root = newroot;
			} //do mode 2

			if (i >= MAXIT) {
				success = false;
				return;
			}
		}

		////////////
		///mode 1///
		////////////

		if (mode == 1) {//SECOND-ORDER GENERAL METHOD (SG)

			for (i = j; i <= MAXIT; i++) {
				faq = 1.0;
				//calculate value of polynomial and its first two derivatives
				p = poly[degree];
				dp = zero;
				d2p_half = zero;
				if ((i - j) % 10 == 0) {
					//prepare stopping criterion
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
						//Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
						//Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
						//ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
						//Eq 8.
						ek = absroot * ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR * ek, 2);
				}
				else {
					for (k = degree; k >= 1; k--) {//Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
						d2p_half = dp + d2p_half * (*root);
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root); //b_k
					}
				}
				abs2p = real(conj(p) * p); //abs(p)**2
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {//(simplified a little Eq. 10 of Adams 1967)
					if (dp == zero) return;
					//do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01 * stopping_crit2) return; //ten times better than stopping criterion
					//ten times better than stopping criterion
					else {
						good_to_go = true; //do one iteration more
					}
				}
				else {
					good_to_go = false; //reset if we are outside the zone of the root
				}
				if (dp == zero) {//test if denominators are > 0.0 not to divide by zero
					dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, FRAC_JUMPS[i % FRAC_JUMP_LEN] * 2 * M_PI)); //make some random jump
				}
				else {
					fac_netwon = p / dp;
					fac_extra = d2p_half / dp;
					F_half = fac_netwon * fac_extra;

					abs2_F_half = real(conj(F_half) * F_half);
					if (abs2_F_half <= 0.000625) {//F<0.05, F/2<0.025
						mode = 0; //set Newton's, go there after jump
					}
					dx = fac_netwon * (c_one + F_half); //SG
				}
				newroot = *root - dx;
				if (newroot == *root) return; //nothing changes -> return
				if (good_to_go) {
					*root = newroot; //this was jump already after stopping criterion was met
					return;
				}
				if (mode != 1) {
					*root = newroot;
					j = i + 1; //remember iteration number
					break; //go to Newton's
				}
				if ((i % FRAC_JUMP_EVERY) == 0) {// decide whether to do a jump of modified length (to break cycles)
					faq = FRAC_JUMPS[(i / FRAC_JUMP_EVERY - 1) % FRAC_JUMP_LEN];
					newroot = *root - faq * dx; //do jump of some semi random lenth (0 < faq < 1)		
				}
				*root = newroot;
			}
			if (i >= MAXIT) {
				success = false;
				return;
			}

		}
		//------------------------------------------------------------------------------- mode 0
		if (mode == 0) { // Newton's Method

			for (i = j; i <= j + 10; i++) { // Do only 10 iterations the most then go back to Laguerre
				faq = 1.0;

				//calc polynomial and first two derivatives
				p = poly[degree];
				dp = zero;
				if (i == j) { // Calculating stopping criterion only at the beginning
					ek = abs(poly[degree]);
					absroot = abs(*root);
					for (k = degree; k >= 1; k--) {
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root);
						ek = absroot * ek + abs(p);
					}
					stopping_crit2 = pow(FRAC_ERR * ek, 2.0);
				}
				else {
					for (k = degree; k >= 1; k--) {
						dp = p + dp * (*root);
						p = poly[k - 1] + p * (*root);
					}
				}
				abs2p = real(conj(p) * p);
				iter = iter + 1;
				if (abs2p == 0.0) return;

				if (abs2p < stopping_crit2) {
					if (dp == zero) return;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01 * stopping_crit2) {
						return; // return immediately since we are at a good place
					}
					else {
						good_to_go = true; // do one more iteration
					}
				}
				else {
					good_to_go = false;
				}

				if (dp == zero) {
					dx = (abs(*root) + 1.0) * expcmplx(complex(0.0, 2 * M_PI * FRAC_JUMPS[i % FRAC_JUMP_LEN])); // make a random jump
				}
				else {
					dx = p / dp;
				}

				newroot = *root - dx;
				if (newroot == *root) return;
				if (good_to_go) {
					*root = newroot;
					return;
				}
				*root = newroot;
			}
			if (iter >= MAXIT) {
				//too many iterations
				success = false;
				return;
			}
			mode = 2; //go back to Laguerre's. Happens when could not converge with 10 steps of Newton
		}

	}/// end of infinite loop
}

#pragma endregion

#pragma region classes


_point::_point(double x, double y, _theta * theta1) {
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

_curve::_curve(_point * p1) {
	length = 1;
	first = last = p1;
	p1->prev = p1->next = 0;
	partneratstart = partneratend = this;
}

_curve::~_curve(void) {
	_point* scan1, * scan2;
	scan1 = first;
	for (int i = 0; i < length; i++) {
		scan2 = scan1->next;
		delete scan1;
		scan1 = scan2;
	}
}

_curve* _curve::divide(_point * ref) {
	_point* scan;
	_curve* nc;
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
	_point* pp;
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

void _curve::append(_point * pp) {

	pp->next = last->next;
	pp->prev = last;
	last->next = pp;
	last = pp;
	length++;
}

void _curve::prepend(double x1, double x2) {
	_point* pp;
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

_curve* _curve::join(_curve * nc) {
	if (length > 0) {
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length > 0) {
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

_curve* _curve::joinbefore(_curve * nc) {
	if (length > 0) {
		first->prev = nc->last;
	}
	else {
		last = nc->last;
	};
	if (nc->length > 0) {
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

_curve* _curve::reverse(void) {
	_point* scan1, * scan2, * scambio;
	if (length > 1) {
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

void _curve::drop(_point * ref) {
	_point* scan;
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

double _curve::closest2(_point * ref, _point * *clos2) {
	double mi = 1.e100, mi2 = 1.e100, FP;
	_point* scan, * clos;
	if (length > 1) {
		clos = *clos2 = first;
		for (scan = first; scan != 0; scan = scan->next) {
			FP = *scan - *ref;
			if (FP < mi) {
				mi2 = mi;
				mi = FP;
				*clos2 = clos;
				clos = scan;
			}
			else if (FP < mi2) {
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

double _curve::closest(_point * ref, _point * *clos) {
	double mi = 1.e100, FP;
	_point* scan;
	for (scan = first; scan != 0; scan = scan->next) {
		FP = *scan - *ref;
		if (FP < mi) {
			mi = FP;
			*clos = scan;
		}
	}
	return mi;
}

void _curve::complement(_point * *sott, int lensott, _point * *res, int lenres) {
	int flag, i;
	_point* scan;
	i = 0;
	for (scan = first; scan != 0; scan = scan->next) {
		flag = 0;
		for (int j = 0; (j < lensott) && (!flag); j++) {
			if (scan == sott[j]) {
				flag = 1;
			}
		}
		if ((!flag) && (i < lenres)) {
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
	_curve* scan1, * scan2;
	scan1 = first;
	while (scan1) {
		scan2 = scan1->next;
		delete scan1;
		scan1 = scan2;
	}
}

void _sols::append(_curve * cc) {
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

void _sols::prepend(_curve * cc) {
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

void _sols::drop(_curve * ref) {
	_curve* scan;
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

void _sols::join(_sols * nc) {
	if (length > 0) {
		last->next = nc->first;
	}
	else {
		first = nc->first;
	};
	if (nc->length > 0) {
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
	_theta* scan, * scan2;
	scan = first;
	while (scan) {
		scan2 = scan->next;
		delete scan;
		scan = scan2;
	}
}

_theta* _thetas::insert(double th) {
	_theta* scan, * scan2;

	scan2 = new _theta(th);
	if (length) {
		if (th < first->th) {
			first->prev = scan2;
			scan2->next = first;
			scan2->prev = 0;
			first = scan2;
		}
		else {
			if (th > last->th) {
				last->next = scan2;
				scan2->prev = last;
				scan2->next = 0;
				last = scan2;
			}
			else {
				scan = first;
				while (scan->th < th) scan = scan->next;
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

void _thetas::remove(_theta * stheta) {
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



#pragma endregion
