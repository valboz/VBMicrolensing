[Back to **Critical curves and caustics**](CriticalCurvesAndCaustics.md)

# Limb darkening

Real stars are not uniform disks. The intensity declines from the center to the edge of the source, as described by limb darkening laws. VBMicrolensing has several built-in limb darkening laws and allows the user to specify an arbitrary law.

## Linear limb darkening

The default limb darkening profile is linear, according to

$$ I(\nu) = I(0) \left[1 - a1 (1 - \nu)\right] $$

with $\nu=\sqrt{1-r^2/\rho^2}$.

The coefficient `VBM.a1` specifies the linear limb darkening coefficient. We note that there are other popular expressions for linear limb darkening that can be reduced to the one adopted in VBMicrolensing through simple transformations (e.g. for $\Gamma_1$ see [An et al. ApJ 572:521 (2002)](https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract), Eq. (11)).

In order to use linear limb darkening in all calculations with extended sources in VBMicrolensing, it is sufficient to give a non-vanishing value to `VBM.a1`:

```
VBM.a1 = 0.51  # Linear limb darkening coefficient.

s = 0.8 #separation between the two lenses
q = 0.1 # mass ratio
y1 = 0.01 # y1 is the source coordinate along the axis parallel to the line joining the two lenses 
y2 = 0.01 # y2 is the source coordinate orthogonal to the first one
rho = 0.01 # Source radius in Einstein radii

Mag = VBM.BinaryMag2(s, q, y1, y2, rho)  # Call to the BinaryMag2 
print(f"Magnification with limb darkened source = {Mag}")  # Output should be 18.27.....

```

In order to go back to uniform sources, just set `VBM.a1 = 0`. In general, a calculation with limb darkening is slower than a calculation with a uniform source, because it is repeated on several concentric annuli.

There are three more limb darkening laws that are already available in VBMicrolensing. You may switch from one to another using the function `VBM.SetLDprofile`.

If you want to go back to linear limb darkening,  you may use ```VBM.SetLDprofile(VBM.LDlinear);```

## Square root limb darkening

$$I(\nu) = I(0) \left[1 - a1 (1 - \nu)- a2 (1-\sqrt{\nu}) \right]$$

This law has two parameters that are given through `VBM.a1` and `VBM.a2`:

```
VBM.SetLDprofile(VBM.LDprofiles.LDsquareroot)
VBM.a1 = 0.51
VBM.a2 = 0.3
Mag = VBM.BinaryMag2(s, q, y1, y2, rho)
print(f"Magnification with square root limb darkened source = {Mag}")  # Output should be 18.2730.....
```

## Quadratic limb darkening

$$I(\nu) = I(0) \left[1 - a1 (1 - \nu)- a2 (1-\nu)^2 \right]$$

```
VBM.SetLDprofile(VBM.LDquadratic)
VBM.a1 = 0.51
VBM.a2 = 0.3
Mag = VBM.BinaryMag2(s, q, y1, y2, rho)
print(f"Magnification with quadratic limb darkened source = {Mag}") # Output should be 18.2728.....
```

## Logarithmic limb darkening

$$I(\nu) = I(0) \left[ 1 - a1 (1 - \nu)- a2 ~ \nu \ln{\nu} \right]$$

```
VBM.SetLDprofile(VBM.LDlog)
VBM.a1 = 0.51
VBM.a2 = 0.3
Mag = VBM.BinaryMag2(s, q, y1, y2, rho)
print(f"Magnification with logarithmic limb darkened source = {Mag}") # Output should be 18.2788.....
```

## User-defined limb darkening

The user-defined limb darkening is available only in the C++ version of VBMicrolensing, see the [corresponding documentation](/docs/C%2B%2B/LimbDarkening.md).

## Multi-band observations

Sometimes, simultaneous observations in multi-band are available. In this case, it is possible to repeat the calculation of the magnification for each band changing `VBM.a1` as appropriate. 

However, there is the possibility to save some time by `BinaryMagMultiDark`. This function exploits the same sequence of concentric annuli to calculate the magnification with different limb darkening coefficients simultaneously. The computational time is then comparable to a single calculation. Here is an example:

```
import numpy as np

nbands = 4 # number of bands
a1_list = np.array([0.2, 0.3, 0.51, 0.6]) # list of limb darkening coefficients
mag_list = np.empty(nbands) # array where to store the output
Tol = 1e-3 # Accuracy (see section about accuracy control)

VBM.BinaryMagMultiDark(s, q, y1, y2, rho, a1_list, nbands, mag_list, Tol);

for i in range(4):
    print("BinaryMagMultiDark at band {}: {}".format(i, mag_list[i]))
```

`BinaryMagMultiDark` currently works with linear limb darkening only. It is a particular function that has been introduced after a specific scientific request, but we believe it can be useful to general users with easy customization if necessary.

[Go to: **Accuracy control**](AccuracyControl.md)
