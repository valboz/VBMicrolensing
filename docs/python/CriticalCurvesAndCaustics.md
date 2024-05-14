[Back to **Binary lenses**](BinaryLenses.md)

# Critical curves and caustics

The gravitational lensing phenomenology is strictly related to the existence of critical curves and caustics (see the recommended [reviews](reviews.md) for a full treatment).

VBMicroLensing offers the calculation of critical curves and caustics with an arbitrary number of points through the functions ```Caustics``` and ```CriticalCurves```.


The result is a list object, which is a collection of pairs of lists for each curve. Each pair of lists represents the x-coordinate and the y-coordinate of the point.

The use of these objects is very intuitive, as illustrated by this examples:

BinaryLens

```
VBML = VBMicrolensing.VBMicrolensing() # Instance to VBMicroLensing

# Parameters of our binary lens
s=2.5;  # separation between the two lenses
q=0.1;  # mass ratio

Caustics = VBM.CausticsBin(s,q)
CriticalCurves=VBM.CriticalCurvesBin(s,q)

```

Critical curves and caustics are calculated through the resolution of a fourth order complex polynomial (see [reviews](reviews.md)) by the [Skowron & Gould algorithm](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/). 

The **number of points** calculated for the critical curves is controlled by ```VBML.NPcrit```, which can be changed by the user according to the desired sampling. The default value is 200.

[Go to: **Limb Darkening**](LimbDarkening.md)
