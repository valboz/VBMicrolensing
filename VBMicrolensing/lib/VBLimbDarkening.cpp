#include "VBMicrolensingLibrary.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

double VBMicrolensing::rCLDprofile(double tc, annulus* left, annulus* right) {
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