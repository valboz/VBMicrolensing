[Back to **Binary lenses**](BinaryLenses.md)


# Binary lenses with shear and convergence

The binary lens model can be extended by considering external mass sheet, which adds shear (float) and convergence (complex) parameters. 

```
#include "VBMicrolensingLibrary.h"

VBMicrolensing VBM;

s = 0.8; // separation between the two lenses
q = 0.1; // mass ratio
y1 = 0.01; // y1 is the source coordinate along the axis parallel to the line joining the two lenses 
y2 = 0.01; // y2 is the source coordinate orthogonal to the first one

kappa = 0.05; // external mass sheet convergence
gamma_real = 0.1; // external mass sheet shear - real part
gamma_imag = -0.06; // external mass sheet shear - imaginary part

Mag = VBM.BinaryMag0_shear(s, q, y1, y2, kappa, gamma_real, gamma_imag);
printf("Binary lens Magnification = %lf\n", Mag)  # Output should be 9.84...
```

The resolution of the lens equation is obtained by recasting the lens equation as a nineth order complex polynomial, whose roots are found by the [Skowron & Gould algorithm](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/), as discussed in the introduction. 

You may check that in the limit of zero shear and zero convergence, the magnification returned is the same as in standard points-source calculations:

```
Mag_1 = VBM.BinaryMag0_shear(s, q, y1, y2, 0., 0., 0.);
Mag_2 = VBM.BinaryMag0(s, q, y1, y2); 
printf("Test: %lf - %lf = %lf\n", Mag_1, Mag_2, Mag_1-Mag_2); // Output should be 18.18... - 18.18... = -1.2...e-12
```

