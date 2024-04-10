
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
#ifdef _PRINT_ERRORS_DARK
			printf("\n%d", Images->length);
#endif

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