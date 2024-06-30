[Back to **Limb Darkening**](LimbDarkening.md)

# Accuracy control

The accuracy goal of the calculation can be controlled by the user through the property ```VBM.Tol```. In fact, ```ESPLMag2```, ```BinaryMag2``` and ```MultiMag``` refine their calculations until they match the required accuracy goal. The result will be ```Mag``` $\pm$ ```VBM.Tol``` (absolute accuracy).

Keep in mind that the computational time typically scales as ```VBM.Tol^(-1/2)```. By default ```VBM.Tol``` is set to ```1.e-2```.

Here are some examples:

```
VBM = VBMicrolensing.VBMicrolensing()
Mag, s, q, y1, y2, Rs = 0, 0.8, 0.1, 0.01, 0.01, 0.01

VBM.Tol = 1e-3
Mag = VBM.BinaryMag2(s, q, y1, y2, Rs)
print("Magnification (accuracy at 1.e-3) =", Mag)  # Output should be 18.283....

VBM.Tol = 1e-4
Mag = VBM.BinaryMag2(s, q, y1, y2, Rs)
print("Magnification (accuracy at 1.e-4) =", Mag)  # Output should be 18.2833....
```

## Precision control

In general, the photometric precision of ground observatories is not better than 0.001. Therefore, it makes sense to set a relative precision goal, instead of an asbolute accuracy goal. This is set by ```VBM.RelTol```. The result will be ```Mag```$\pm$ ```Mag*VBM.RelTol``` (relative precision).

For example, let us set a poor 10% precision, just to see the difference:

```
VBM.RelTol = 1e-1
Mag = VBM.BinaryMag2(s, q, y1, y2, Rs)
print("Magnification (relative precision at 1.e-1) =", Mag)  # Output should be 18.24....
```

If you do not want to use relative precision, just set ```VBM.RelTol = 0;``` which is the default value.

In general, the calculation stops when the first of the two goals is reached, either absolute accuracy or relative precision.

## Accuracy in astrometry

VBMicrolensing allows full direct control on the accuracy of the magnification calculation. The accuracy on the astrometric centroid calculation scales similarly, but is also affected by the source size. As discussed in [V. Bozza, E. Khalouei and E. Bachelet, MNRAS 505 (2021) 126](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505..126B/abstract), a useful formula to track the astrometric accuracy $\delta_{ast}$ as a function of the magnification accuracy $\delta_{phot}$ is

$$\delta_{ast} \simeq 50 \rho ~ \delta_{phot}$$ 

[Go to **Light Curve Functions**](LightCurves.md)
