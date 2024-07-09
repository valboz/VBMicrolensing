
#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//////////////////////////////
//////////////////////////////
// Polynomial methods
//////////////////////////////
//////////////////////////////


void VBMicrolensing::polyproduct(complex *p1, int n1, complex *p2, int n2, complex *pdest) {
	// polynomials are with increasing degree: p1[0]+p1[1]*z+p1[2]*z^2 + ...
	for (int i = 0; i <= n1 + n2; i++) pdest[i] = 0;
	for (int i = 0; i <= n1; i++) {
		for (int j = 0; j <= n2; j++) {
			pdest[i + j] = pdest[i + j] + p1[i] * p2[j];
		}
	}
}

void VBMicrolensing::copypol(complex *p1, int n1, complex *pdest) {
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
				dg = (j > i) ? i * n : (i - 1)* n;
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
