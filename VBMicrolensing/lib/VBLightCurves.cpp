
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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



void VBMicrolensing::TripleLightCurve(double *pr, double *ts, double *mags, double *y1s, double *y2s, int np) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]),di,mindi;
	double q[3] = { 1, exp(pr[1]),exp(pr[8]) };
	complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta=sin(pr[9]), cbeta=cos(pr[9]);

	s[0] = exp(pr[0]) / (q[0] + q[1]);
	s[1] = s[0] * q[0];
	s[0] = -q[1]* s[0];
	s[2] = exp(pr[7])*complex(cbeta, sbeta)+s[0];
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

			mags[i]= 1.;
		}
		else {
			mags[i] = MultiMag(complex(y1s[i], y2s[i]), rho, Tol);
		}
	}
}

void VBMicrolensing::TripleLightCurveParallax(double* pr, double* ts, double* mags, double* y1s, double* y2s, int np) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), di, mindi, u, u0=pr[2], t0=pr[6], pai1 = pr[10], pai2 = pr[11];
	double q[3] = { 1, exp(pr[1]),exp(pr[8]) };
	complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);
	double Et[2];

	s[0] = exp(pr[0]) / (q[0] + q[1]);
	s[1] = s[0] * q[0];
	s[0] = -q[1] * s[0];
	s[2] = exp(pr[7]) * complex(cbeta, sbeta)+s[0];
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

	double* q= (double*)malloc(sizeof(double) * (nl));
	complex* s= (complex*)malloc(sizeof(complex) * (nl));

	q[0] = 1.;
	for (int i = 1, j =1; i < nl; ++i, j+=3) {
		q[i] = pr[j + 5];
	}

	s[0] = complex(0, pr[3]);
	for (int i = 1, j=1; i < nl; ++i, j+=3) {
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


void VBMicrolensing::BinaryLightCurveOrbital(double *pr, double *ts, double *mags, double *y1s, double *y2s, double *seps, int np) {
	double s = exp(pr[0]), q = exp(pr[1]), u0 = pr[2], rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]), t0 = pr[6], pai1 = pr[7], pai2 = pr[8], w1 = pr[9], w2 = pr[10], w3 = pr[11];
	double salpha = sin(pr[3]), calpha = cos(pr[3]);
	double Et[2];
	double w, phi0, inc, phi, Cinc, Sinc, Cphi, Sphi, Cphi0, Sphi0, COm, SOm,s_true;
	double w13, w123, den, den0, u;
	t0old = 0;

	w13 = w1*w1 + w3*w3;
	w123 = sqrt(w13 + w2*w2);
	w13 = sqrt(w13);
	if (w13>1.e-8) {
		w3 = (w3>1.e-8) ? w3 : 1.e-8;
		w = w3*w123 / w13;
		inc = acos(w2*w3 / w13 / w123);
		phi0 = atan2(-w1*w123, w3*w13);
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
	den0 = sqrt(Cphi0*Cphi0 + Cinc*Cinc*Sphi0*Sphi0);
	s_true = s / den0; // orbital radius
	COm = (Cphi0*calpha + Cinc*salpha*Sphi0) / den0;
	SOm = (Cphi0*salpha - Cinc*calpha*Sphi0) / den0;

	for (int i = 0; i < np; i++) {
		ComputeParallax(ts[i], t0, Et);

		phi = (ts[i] - t0_par)*w + phi0;
		Cphi = cos(phi);
		Sphi = sin(phi);
		den = sqrt(Cphi*Cphi + Cinc*Cinc*Sphi*Sphi);
		seps[i] = s_true*den; // projected separation at time ts[i]

		u = u0 + pai1*Et[1] - pai2*Et[0];
		tn = (ts[i] - t0) * tE_inv + pai1*Et[0] + pai2*Et[1];
		y1s[i] = (Cphi*(u*SOm - tn*COm) + Cinc*Sphi*(u*COm + tn*SOm)) / den;
		y2s[i] = (-Cphi*(u*COm + tn*SOm) - Cinc*Sphi*(tn*COm - u*SOm)) / den;
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
		M = n * (ts[i] - tperi);
		EE = M + e * sin(M);
		dE = 1;
		while (fabs(dE) > 1.e-8) {
			dM = M - (EE - e * sin(EE));
			dE = dM / (1 - e * cos(EE));
			EE += dE;
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

double VBMicrolensing::ESPLLightCurveParallax(double *pr, double t) {
	double u0 = pr[0], t0 = pr[2], tE_inv = exp(-pr[1]), tn, u, u1, rho = exp(pr[3]), pai1 = pr[4], pai2 = pr[5];
	double Et[2];

	ComputeParallax(t, t0, Et);
	tn = (t - t0) * tE_inv + pai1 * Et[0] + pai2 * Et[1];
	u1 = u0 + pai1 * Et[1] - pai2 * Et[0];
	u = sqrt(tn*tn + u1 * u1);

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

double VBMicrolensing::TripleLightCurve(double *pr, double t) {
	double rho = exp(pr[4]), tn, tE_inv = exp(-pr[5]),di,mindi;
	static double q[3];
	static double prold[] = { 0,0,0,0,0 };
	static Method oldmethod = Method::Nopoly;
	static complex s[3];
	double salpha = sin(pr[3]), calpha = cos(pr[3]), sbeta = sin(pr[9]), cbeta = cos(pr[9]);
	int inew=0;
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
		s[2] = exp(pr[7])*complex(cbeta, sbeta)+s[0];
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
	EE = M + e * sin(M);
	dE = 1;
	while (fabs(dE) > 1.e-8) {
		dM = M - (EE - e * sin(EE));
		dE = dM / (1 - e * cos(EE));
		EE += dE;
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

