// MultipleLens.cpp : main project file.

#include "VBMicrolensingLibrary.h"
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace System;

void PrintCau();
void PrintImages();
void PrintConfiguration();

// If you wnat critical curves, caustics, save imagese or a particular configuration
//you can use this functions Print Cau(), PrintImages() and PrintCofiguration() in the main.

struct Point {
	double x;
	double y;
};

//Here you can insert your specific configuration.

const int n0 = 10; // Number of lenses inserted in the lists
int n1 = 4;			// The calculations will consider the first "n1" lenses.

//q denotes the list of lens masses.
//s represents the list of lens positions on the complex plane.

double q[n0] = { 1.,1.e-1,1.1e-4,1.1e-6,1.5e-1,0.0015,0.001,0.001,0.001,0.0002 };
complex s[n0] = { complex(0.,0.),complex(1.,-0.7),complex(2,0.7),complex(0.6,-.6), complex(0.3,0.),complex(-0.3,0.1),0.,0.,0.,0. };

_sols* images; //images
complex y = complex(0.1, 0.04); // source position
double rho = .0011; //source radius


// Declaration of an instance to VBMicroLensing class
VBMicrolensing VBML;

int main()
{
	//////////////////////////////////
	//// 
	//// Microlensing magnification for a specific source position
	////
	/////////////////////////////////
	// 
	//double Mag;
	//double Toll = 0.001;
	//// Choose the algorithm you want to use (Nopoly, Multipoly, Singlepoly).
	//VBML.SetMethod(VBMicrolensing::Method::Multipoly);
	//
	//// Call the function setLensGeometry, which initializes the arrays 
	//// describing the geometric configuration of the lens system based 
	//// on the chosen algorithm and configuration.
	//VBML.SetLensGeometry(n1,q,s);

	//// Accuracy control
	//VBML.Tol = 1.e-2;
	//VBML.RelTol = 1.e-3;

	//// Call to the MultiMag function with these parameters
	//Mag = VBML.MultiMag(y, rho, Toll, &images);
	//printf("Microlensing Magnification = %lf\n", Mag);

	//getchar();
	//return 0;

	//////////////////////////////////
	//// 
	//// Draw lightcurve from standard parameters, TRIPLE LENS
	////
	/////////////////////////////////
	//const int np = 1000;
	//double tmin = -50., tmax = 50., t;
	//double pr[10] = { log(0.765),log(0.00066),0.0060,3.212,log(0.00567),log(50.13),0,log(1.5),log(0.000001),-1.5 };
	//double ts[np], y1s[np], y2s[np], mags[np];
	//char buffer[1024];
	//FILE* ff;

	//VBML.Tol = 1.e-2;
	//VBML.RelTol = 1.e-3;
	////Method
	//VBML.SetMethod(VBMicrolensing::Method::Nopoly);

	//for (int i = 0; i < np; i++) {
	//	t = pr[6] + tmin + i * (tmax - tmin) / np;
	//	ts[i] = t;
	//}
	//VBML.TripleLightCurve(pr, ts, mags, y1s, y2s, np);
	//FILE* f = fopen("LightCurve.txt", "w");
	//printf("Done\n");
	//for (int i = 0; i < np; i++) {
	//	fprintf(f, "%lf %lf\n", ts[i], mags[i]);
	//}
	//fclose(f);
	//PrintCau();
	//getchar();
	//return 0;


	////////////////////////////////
	// 
	// Draw lightcurve from standard parameters, MULTIPLE LENS
	//
	///////////////////////////////
	const int np = 1000, nl = 4;
	double tmin = -50., tmax = 50., t;
	double pr[3 * nl + 4];
	double ts[np], y1s[np], y2s[np], mags[np];
	char buffer[1024];
	FILE* ff;

	double q[nl-1] = { 1.e-1,1.1e-4,1.1e-6};//lens masses
	double s1[nl-1] = { 1.0, 2.0, 0.6 }; //real part lens positions without first lens
	double s2[nl-1] = { -0.7, 0.7, -0.6};//im part lens positions

	pr[0] = 0.; // t0
	pr[1] = 0.; // log(tE)
	pr[2] = -2.52287874528; // log(rho)
	pr[3] = 0.; // im position first lens 

	for (int i = 4, j = 0; i < 1 + 3*nl; i += 3, ++j) {
		pr[i] = s1[j];
		pr[i + 1] = s2[j];
		pr[i + 2] = q[j];
	}

	VBML.Tol = 1.e-3;
	VBML.RelTol = 1.e-3;
	//Method
	VBML.SetMethod(VBMicrolensing::Method::Nopoly);

	for (int i = 0; i < np; i++) {
		t = pr[0] + tmin + i * (tmax - tmin) / np;
		ts[i] = t;
	}
	VBML.LightCurve(pr, ts, mags, y1s, y2s, np, nl);
	FILE* f = fopen("LightCurve.txt", "w");
	printf("Done\n");
	for (int i = 0; i < np; i++) {
		fprintf(f, "%lf %lf\n", ts[i], mags[i]);
	}
	fclose(f);
	PrintCau();
	getchar();
	return 0;
}




void PrintCau() {
	_sols *CriticalCurves;
	_curve *scancurve;
	_point *scanpoint;
	FILE *f;
	int ncc;

	CriticalCurves = VBML.PlotCrit();
	f = fopen("outcurves.causticdata", "w");
	ncc = CriticalCurves->length;

	scancurve = CriticalCurves->first;
	for (int i = 0; i < ncc; i++) {
		fprintf(f, "Curve: %d\n", i + 1);
		for (scanpoint = scancurve->first; scanpoint; scanpoint = scanpoint->next) {
			fprintf(f, "%lf %lf\n", scanpoint->x1, scanpoint->x2);
		}
		scancurve = scancurve->next;
	}
	fclose(f);
	delete CriticalCurves;
}

void PrintImages() {
	FILE *f;

	f=fopen("outcurves.txt","w");
	_curve *c = images->first;
	for(int i=0;i<images->length;i++){
		fprintf(f, "Curve: %d\n", i + 1);
		for(_point *p=c->first;p;p=p->next){
			fprintf(f,"%.16lf %.16lf\n",p->x1,p->x2);
		}
		c = c->next;
	}
	fclose(f);
}

void PrintConfiguration() {
	FILE* f;

	f = fopen("geometry.txt", "w");
	fprintf(f, "%d\n", n1);
	for (int i = 0; i < n1; i++) {
		fprintf(f, "%.16lf %.16lf %.16lf\n", q[i], s[i].re, s[i].im);
	}
	fprintf(f, "%.16lf %.16lf %.16lf\n", y.re, y.im, rho);
	fclose(f);
}