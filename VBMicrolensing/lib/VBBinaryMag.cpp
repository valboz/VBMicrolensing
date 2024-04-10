
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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
	if (Mag < 0) {
		mag1 = -1;
		delta1 = 3.33333333e-8;
		while (mag1 < 0.1 && RSi >= 0) {
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
		while (mag2 < 0.1) {
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

double VBDefaultCumulativeFunction(double cb, double *a1) {
	static double r2, cr2, scr2, cc;
	r2 = cb * cb;
	cr2 = 1 - r2;
	scr2 = sqrt(cr2);
	cc = (3 * r2*(1 - *a1) - 2 * (*a1)*(scr2*cr2 - 1)) / (3 - *a1);
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


_curve* VBMicrolensing::NewImages(complex yi, complex* coefs, _theta* theta) {
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


void VBMicrolensing::OrderImages(_sols* Sols, _curve* Newpts) {
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
