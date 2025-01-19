[Back to **Summary**](readme.md)

# Single lenses


The lens equation for a single lens is

$$ u=x-\frac{1}{x}$$

For a given lens-source angular separation $u$ there are two images $x$ solving this equation. We note that all angular distances are in units of the Einstein angle $\theta_E$. 

## Point-Source-Point-Lens

If the source is point-like, the magnification is given by the famous Paczynski formula

$$ \mu = \frac{u^2+2}{u\sqrt{u^2+4}}$$

In VBMicrolensing, this formula is obtained through the function ```PSPLMag``` as follows:

```
VBM = VBMicrolensing.VBMicrolensing()

u=0.1;  # Source-lens separation in Einstein radii
Mag = VBM.PSPLMag(u)
print(f"PSPL Magnification = {Mag}") # Output should be 10.037...
```

## Extended-Source-Point-Lens

For extended sources, the magnification depends on $\rho$, the **source radius** normalized to the Einstein angle, and can be calculated through elliptic integrals. In order to make VBMicrolensing as fast as possible, we provide **pre-calculated tables** in the file "ESPL.tbl". This file is automatically loaded upon the first calculation involving Extended-Source-Point-Lenses (ESPL).

```
VBM = VBMicrolensing.VBMicrolensing()

u = 0.1 # Source-lens separation in Einstein radii
rho = 0.01 # Source radius in units of the Einstein angle

Mag = VBM.ESPLMag2(u, rho) # Call to the ESPLMag2 function with these parameters
print(f"\nMagnification of Extended-source-point-lens = {Mag}\n")  # Output should be 10.050.....

```

The current range for the pre-calculated table is $10^{-4} \leq \rho \leq 10^{+2}$. Sources smaller than the minimum are considered equal to the minimum. Sources larger than the maximum generate an error message. 

By default, VBMicrolensing works with **uniform sources**. We will introduce **Limb Darkening** in a [later section](LimbDarkening.md): arbitrary Limb Darkening laws can be implemented in VBMicrolensing.

## Astrometry

For a Point-Source, in the reference frame in which the **lens is in the origin**, the centroid of the images is simply

$$ \bar x = \frac{u}{u^2+2} + u$$

If you need astrometry calculations together with magnification, you have to turn astrometry on by ```VBM.astrometry = True``` and read the results in ```VBM.astrox1```. This works in the same way for ```PSPLMag``` and ```ESPLMag2```.

```
VBM = VBMicrolensing.VBMicrolensing()

u = 0.1 # Source-lens separation in Einstein radii
rho = 0.01 # Source radius in units of the Einstein angle

VBM.astrometry = True # We want astrometry

Mag = VBM.ESPLMag2(u, rho) # Call to the ESPLMag2 function with these parameters
print(f"Magnification of Extended-source-point-lens = {Mag}\n")  # Output should be 10.050.....
print(f"Centroid shift = {VBM.astrox1 - u}\n")  # Output should be 0.0493.....

```

Note that ```VBM.astrox1``` reports the **centroid position with respect to the lens**. The **centroid position with respect to the source** is ```VBM.astrox1 - u```.

[Go to: **Binary lenses**](BinaryLenses.md)
