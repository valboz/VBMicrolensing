#include <stdio.h>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>


#include"VBMicrolensingLibrary.h"

void PrintCau(double,double);
void PrintMultiCau();
void PrintImages(_sols *);

double y_1, y_2, s, q, rho, mag,mag2;
VBMicrolensing VBM;
const int np = 3000;

//int main(){
//	VBMicrolensing VBM;
//
//	double pr[7]; // Array of parameters
//	double u0, t0, tE, rho, alpha, s, q;
//	const int np = 100; // Number of points in the light curve
//	double times[np]; // Array of times 
//	double mags[np]; // Array of magnifications where to store the output
//	double y1s[np]; // Array of source coordinates y1 for each time
//	double y2s[np]; // Array of source coordinates y2 for each time
//
//
//	u0 = -0.01; // Impact parameter
//	t0 = 7550.4; // Time of closest approach to the center of mass
//	tE = 100.3; // Einstein time
//	rho = 0.01; // Source radius
//	s = 0.8; // Separation between the two lenses
//	q = 0.1; // Mass ratio
//	alpha = 0.53; // Angle of the source trajectory
//
//	pr[0] = log(s);
//	pr[1] = log(q);
//	pr[2] = u0;
//	pr[3] = alpha;
//	pr[4] = log(rho);
//	pr[5] = log(tE);
//	pr[6] = t0;
//
//	for (int i = 0; i < np; i++) {
//		times[i] = t0 - tE + 2 * tE / np * i;
//	}
//	VBM.astrometry = true;
//	VBM.BinaryLightCurve(pr,times, mags, y1s,y2s,np); // Calculates the Binary Lens magnification at time t with parameters in pr
//	
//	for (int i = 0; i < np; i++) {
//		printf("\n%lf %lf %lf %lf", times[i], mags[i], y1s[i], y2s[i]);
//	}
//	getchar();
//}

int main() {
	VBMicrolensing VBM;

	int nn = 3; // Number of lenses

	double pr[] = {     //parameters
		0., 0.0, 1,    // First lens: x1_1, x1_2, m1
		0.8, 0., 0.1,  // Second lens: x2_2, x2_2, m2
		0, 0.8, 1e-2,   // Third lens: x3_re, x3_im, m3
	};

	VBM.SetLensGeometry(nn, pr); //Initialize the lens configuration

	double x_cm = (0.0 * 1 + 0.8 * 0.1) / (1 + 0.1);
	double y1 = 0.01 + x_cm; // add center of mass x to y1
	double y2 = 0.01;
	double rho = 0.01; //source radius

	VBM.astrometry = true;

	double Mag = VBM.MultiMag(y1, y2, rho);
	printf("Multiple Lens Magnification = %f", Mag); 
	printf("\nCentroid shift = (%lf,%lf)\n", VBM.astrox1 - y1, VBM.astrox2 - y2);  // Output should be (-0.164...,-0.074...)
	getchar();

	

}

//
//
//int main() {
//	VBMicrolensing VBM;
//
//	double Mag, s, q, y1, y2, rho;
//
//	s = 0.8; //separation between the two lenses
//	q = 0.1; // mass ratio
//	y1 = 0.01; // y1 is the source coordinate along the axis parallel to the line joining the two lenses 
//	y2 = 0.01; // y2 is the source coordinate orthogonal to the first one
//	rho = 0.01; // Source radius in Einstein radii
//
//	VBM.astrometry = true; // We want astrometry
//
//	Mag = VBM.BinaryMag0(s, q, y1, y2); // Call to the BinaryMag2 function with these parameters
//	printf("Binary lens Magnification = %lf\n", Mag); // Output should be 18.28....
//	printf("\nCentroid shift = (%lf,%lf)\n", VBM.astrox1 - y1, VBM.astrox2 - y2);  // Output should be (-0.164...,-0.074...)
//	printf("(y1,y2): (%ef,%ef)", y1, y2);
//	getchar();
//
//	//Binarymag0: -0.163203, -0.073453 mag 18.185448
//
//}

//int main()
//{
//
//	VBMicrolensing VBM;
//
//	VBM.SetMethod(VBMicrolensing::Method::Multipoly); //Choose the method: Nopoly, Multipoly, Singlepoly
//
//
//	double pr[10]; // Array of parameters
//	double u0, t0, tE, rho, alpha, s12, q2, s13, q3, psi;
//	double t, Mag;
//	const int np = 100; // Number of points in the light curve
//	double times[np]; // Array of times 
//	double mags[np]; // Array of magnifications where to store the output
//	double y1s[np]; // Array of source coordinates y1 for each time
//	double y2s[np]; // Array of source coordinates y2 for each time
//
//	u0 = 0.0060; // Impact parameter with respect to center of mass of the first two lenses
//	t0 = 0.0; // Time of closest approach to center of mass of the first two lenses
//	tE = 50.13; // Einstein radius crossing time
//	rho = 0.0567; // Source radius 
//	s12 = 0.765; // Separation between the second and first lens
//	q2 = 0.00066; // Mass ratio of the second lens to the primary
//	alpha = 3.212; // Angle of the source trajectory
//	s13 = 1.5; // Separation between the third lens and the primary
//	q3 = 0.000001; // Mass ratio of the third lens to the primary
//	psi = -1.5; // Angle between second and third lens as shown in figure
//
//	pr[0] = log(s12);
//	pr[1] = log(q2);
//	pr[2] = u0;
//	pr[3] = alpha;
//	pr[4] = log(rho);
//	pr[5] = log(tE);
//	pr[6] = t0;
//	pr[7] = log(s13);
//	pr[8] = log(q3);
//	pr[9] = psi;
//
//	// Suppose we want to simulate the light curve with a fixed time sampling
//	for (int i = 0; i < np; i++) {
//		times[i] = t0 - tE + 2 * tE / np * i;
//	}
//
//	VBM.TripleLightCurve(pr, times, mags, y1s, y2s, np); // Calculates the Triple Lens magnification at time t with parameters in pr
//	
//	// Let's print the results
//	for (int i = 0; i < np; i++) {
//		printf("\n%lf %lf %lf %lf", times[i], mags[i], y1s[i], y2s[i]);
//	}
//
//	// Scrittura su file con fprintf
//	FILE* fout = fopen("lightcurve.txt", "w");
//	if (fout == NULL) {
//		perror("Errore nell'apertura del file lightcurve.txt");
//		return 1;
//	}
//
//	for (int i = 0; i < np; i++) {
//		fprintf(fout, "%lf\t%lf\t%lf\t%lf\n", times[i], mags[i], y1s[i], y2s[i]);
//	}
//	fclose(fout);
//	getchar();
//	PrintMultiCau();
//	getchar();
//
//}


void PrintCau(double s, double q) {
	_sols* CriticalCurves;
	_curve* scancurve;
	_point* scanpoint;
	FILE* f;
	int ncc;

	CriticalCurves = VBM.PlotCrit(s,q);
	f = fopen("outcurves.causticdata", "w");
	ncc = CriticalCurves->length;

	scancurve = CriticalCurves->first;
	for (int i = 0; i < ncc; i++) {
		fprintf(f, "Curve: %d\n", i + 1);
		for (scanpoint = scancurve->first; scanpoint; scanpoint = scanpoint->next) {
			fprintf(f, "%.16lf %.16lf\n", scanpoint->x1, scanpoint->x2);
		}
		scancurve = scancurve->next;
	}
	fclose(f);
	delete CriticalCurves;
}


void PrintMultiCau() {
	_sols* CriticalCurves;
	_curve* scancurve;
	_point* scanpoint;
	FILE* f;
	int ncc;

	CriticalCurves = VBM.PlotCrit();
	f = fopen("outcurves.causticdata", "w");
	ncc = CriticalCurves->length;

	scancurve = CriticalCurves->first;
	for (int i = 0; i < ncc; i++) {
		fprintf(f, "Curve: %d\n", i + 1);
		for (scanpoint = scancurve->first; scanpoint; scanpoint = scanpoint->next) {
			fprintf(f, "%.16lf %.16lf\n", scanpoint->x1, scanpoint->x2);
		}
		scancurve = scancurve->next;
	}
	fclose(f);
	delete CriticalCurves;
}

void PrintImages(_sols * images) {
	FILE* f;

	f = fopen("outcurves.txt", "w");
	_curve* c = images->first;
	for (int i = 0; i < images->length; i++) {
		fprintf(f, "Curve: %d\n", i + 1);
		for (_point* p = c->first; p; p = p->next) {
			fprintf(f, "%.16lf %.16lf\n", p->x1, p->x2);
		}
		c = c->next;
	}
	fclose(f);

	f = fopen("geometry.txt", "w");
	fprintf(f, "%d\n", 2);
	fprintf(f, "%.16lf %.16lf %.16lf\n", 1 / (1 + q), -q*s / (1 + q), 0);
	fprintf(f, "%.16lf %.16lf %.16lf\n", q/(1+q), s/(1+q), 0);
	fprintf(f, "%.16lf %.16lf %.16lf\n", y_1, y_2, rho);
	fclose(f);
	//f = fopen("geometry.txt", "w");
	//fprintf(f, "%d\n", n1);
	//for (int i = 0; i < n1; i++) {
	//	fprintf(f, "%.16lf %.16lf %.16lf\n", q[i], s[i].re, s[i].im);
	//}
	//fprintf(f, "%.16lf %.16lf %.16lf\n", y.re, y.im, rho);
	//fclose(f);

}