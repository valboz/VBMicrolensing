
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


//////////////////////////////
//////////////////////////////
////////Critical curves and caustics
//////////////////////////////
//////////////////////////////


_sols *VBMicrolensing::PlotCrit() {
	complex ej, y, z;
	int NPS = 200;
	_sols *CriticalCurves;
	_curve *Prov, *Prov2, *isso;
	_point *pisso;
	double SD, MD, CD;

	CriticalCurves = new _sols;
	for (int i = 0; i<2 * n; i++) {
		Prov = new _curve;
		CriticalCurves->append(Prov);
	}

	for (int j = 0; j<NPcrit; j++) {
		ej = complex(cos(2 * j*M_PI / NPcrit), -sin(2 * j*M_PI / NPcrit));
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
		for (Prov2 = Prov->next; Prov2; Prov2 = Prov2->next) {
			CD = *(Prov2->first) - *(Prov->last);
			if (CD<MD) {
				MD = CD;
				isso = Prov2;
			}
		}
		if (MD<SD) {
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
		for (_point *scanpoint = Prov->first; scanpoint; scanpoint = scanpoint->next) {
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

_sols *VBMicrolensing::PlotCrit(double a1, double q1) {
	complex  a, q, ej, zr[4], x1, x2;
	int NPS = 200;
	_sols *CriticalCurves;
	_curve *Prov, *Prov2, *isso;
	_point *pisso;
	double SD, MD, CD, centeroffset;

	a = complex(a1, 0.0);
	q = complex(q1, 0.0);
	centeroffset = a1 / 2.0*(1.0 - q1) / (1.0 + q1);

	CriticalCurves = new _sols;
	for (int i = 0; i<4; i++) {
		Prov = new _curve;
		CriticalCurves->append(Prov);
	}

	for (int j = 0; j<NPcrit; j++) {
		ej = complex(cos(2 * j*M_PI / NPcrit), -sin(2 * j*M_PI / NPcrit));
		complex  coefs[5] = { a*a / 16.0*(4.0 - a * a*ej)*(1.0 + q),a*(q - 1.0),(q + 1.0)*(1.0 + a * a*ej / 2.0),0.0,-(1.0 + q)*ej };
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
		for (Prov2 = Prov->next; Prov2; Prov2 = Prov2->next) {
			CD = *(Prov2->first) - *(Prov->last);
			if (CD<MD) {
				MD = CD;
				isso = Prov2;
			}
		}
		if (MD<SD) {
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
		for (_point *scanpoint = Prov->first; scanpoint; scanpoint = scanpoint->next) {
			x1 = complex(scanpoint->x1 - centeroffset, 0.0);
			x2 = complex(scanpoint->x2, 0.0);
			Prov2->append(real(_L1) + centeroffset, real(_L2));
		}
		CriticalCurves->append(Prov2);
	}
	return CriticalCurves;
}


//////////////////////////////
//////////////////////////////
////////Parallax settings and computation
//////////////////////////////
//////////////////////////////


void VBMicrolensing::SetObjectCoordinates(char *modelfile, char *sateltabledir) {
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

void VBMicrolensing::SetObjectCoordinates(char *CoordinateString) {
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

void VBMicrolensing::ComputeParallax(double t, double t0, double *Et) {
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
			inc = (inc0 + incdot * ty)*deg;
			L = (L0 + Ldot * ty)*deg;
			om = (om0 + omdot * ty)*deg;

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
			y1 = a * sqrt(1 - e * e)*sin(EE);
			//		r=a*(1-e*cos(EE));
			vx = -a / (1 - e * cos(EE))*sin(EE)*Ldot*deg / 36525;
			vy = a / (1 - e * cos(EE))*cos(EE)*sqrt(1 - e * e)*Ldot*deg / 36525;

			Ear[0] = x1 * cos(om) - y1 * sin(om);
			Ear[1] = x1 * sin(om)*cos(inc) + y1 * cos(om)*cos(inc);
			Ear[2] = x1 * sin(om)*sin(inc) + y1 * cos(om)*sin(inc);
			vEar[0] = vx * cos(om) - vy * sin(om);
			vEar[1] = vx * sin(om)*cos(inc) + vy * cos(om)*cos(inc);
			vEar[2] = vx * sin(om)*sin(inc) + vy * cos(om)*sin(inc);

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
		inc = (inc0 + incdot * ty)*deg;
		L = (L0 + Ldot * ty)*deg;
		om = (om0 + omdot * ty)*deg;

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
		y1 = a * sqrt(1 - e * e)*sin(EE);
		//	r=a*(1-e*cos(EE));

		Ear[0] = x1 * cos(om) - y1 * sin(om);
		Ear[1] = x1 * sin(om)*cos(inc) + y1 * cos(om)*cos(inc);
		Ear[2] = x1 * sin(om)*sin(inc) + y1 * cos(om)*sin(inc);
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
