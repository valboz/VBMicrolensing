[Back to **Light Curve Functions**](LightCurves.md)

# Parallax

In the basic light curve functions discussed in the [Light Curve Functions](LightCurves.md) section, we assumed that observer, lens and source move rectilinearly. Therefore, the relative lens-source angular motion is rectilinear and described by the two parameters $u_0$ and $\alpha$.

In reality, the observer motion is curvilinear as the Earth orbits the Sun. For satellite observations, the situation is similar, since all spacecraft will describe curved orbits. In this section we will introduce a set of functions similar to the "static" ones discussed before, but including the observer motion.

First we will treat geocentric observations and then we will move to observations from space.

## Target coordinates

We need to specify J2000.0 equatorial coordinates for our microlensing target. This is done by function `SetObjectCoordinates`:

```
VBM.SetObjectCoordinates("17:59:02.3 -29:04:15.2") # Assign RA and Dec to our microlensing event
```
The string given as argument here obviously contains Right Ascension and Declination in standard notations.

## Parallax system

Then we should decide in which coordinates we want to express the parallax vector $\vec \pi_E$. In the literature, there are two popular choices: North-East system $(\pi_{E_,N},\pi_{E,E})$ and parallel/orthogonal to the Earth acceleration direction $(\pi_{E,\parallel},\pi_{E,\perp})$. In VBMicrolensing you have both possibilities by setting `VBM.parallaxsystem` to 1 or 0 respectively. The default value is 0, corresponding to the parallel/orthogonal system.

## Reference time for parallax $t_{0,par}$

Finally, we have to decide the reference time for the parallax effect $t_{0,par}$. Whatever the values of the parallax components, the source position at $t=t_{0,par}$ remains fixed. By default, VBMicrolensing uses $t_{0,par}=t_0$, so that the light curve is unchanged at the time of closest approach to the center of mass of the lens. However, if you want to keep the source position at another time fixed, you can set `VBM.t0_par_fixed = 1;` and choose your reference time via `VBM.t0_par`.

## Light curve functions with parallax

All light curve functions defined in the [Light Curves](LightCurves.md) section have their corresponding counterpart including the parallax effect:

```
PSPLLightCurveParallax
ESPLLightCurveParallax
BinaryLightCurveParallax
TripleLightCurveParallax
```

The only difference is that the array of parameters must include two more entries for the components of the parallax vector. Here is a full example demonstrating the use of `BinaryLightCurveParallax`:

```
import VBMicrolensing
import math
import numpy as np
import matplotlib.pyplot as plt

VBM = VBMicrolensing.VBMicrolensing()

s = 0.9       # Separation between the lenses
q = 0.1       # Mass ratio
u0 = 0.0       # Impact parameter with respect to center of mass
alpha = 1.0       # Angle of the source trajectory
rho = 0.01       # Source radius
tE = 30.0      # Einstein time in days
t0 = 7500      # Time of closest approach to center of mass
paiN = 0.3     # North component of the parallax vector
paiE = -0.2     # East component of the parallax vector

# Array of parameters. Note that s, q, rho and tE are in log-scale
pr = [math.log(s), math.log(q), u0, alpha, math.log(rho), math.log(tE), t0, paiN, paiE]

t = np.linspace(t0-tE, t0+tE, 300) # Array of times

VBM.parallaxsystem = 1       # Set parallax system to North-East
VBM.SetObjectCoordinates("17:59:02.3 -29:04:15.2") # Assign RA and Dec to our microlensing event

magnifications, y1, y2 = VBM.BinaryLightCurve(pr,t)      # Calculation of static binary-lens light curve
magnificationspar, y1par, y2par = VBM.BinaryLightCurveParallax(pr,t)      # Calculation of light curve with parallax

plt.plot(t,magnifications,"g")
plt.plot(t,magnificationspar,"m")
```

<img src="BinaryLens_lightcurve_parallax.png" width = 400>

The light curve including parallax is in magenta in this plot. 

And here is the corresponding source trajectory
```
caustics = VBM.Caustics(s,q)
for cau in caustics:
    plt.plot(cau[0],cau[1])
plt.plot(y1,y2,"g")
plt.plot(y1par,y2par,"m")
```

<img src="BinaryLens_lightcurve_parallax_caustics.png" width = 400>

In this example we have not set `VBM.t0_par`, which means that $t_{0,par}=t_0$ here.

## Satellite Parallax

VBMicrolensing can calculate the magnification as seen from a spacecraft. In order to do that, it is necessary to have the ephemerides of the satellite in the format given by the [NASA Horizons system](http://ssd.jpl.nasa.gov/horizons.cgi).

In particular, we assume five columns:
- JD
- RA (degrees)
- Dec (degrees)
- Distance from Earth (AU)
- Distance rate change (not really needed but included by default in Horizons).

Examples of valid satellite ephemerid tables are in [https://github.com/valboz/VBMicrolensing/tree/master/VBMicrolensing/data](https://github.com/valboz/VBMicrolensing/tree/master/VBMicrolensing/data).

The satellite table(s) should be named "satellite*.txt" (with * replaced by a single character) and placed in a single directory. In order to inform VBMicrolensing of these tables, there is an alternative version of the `VBM.SetObjectCoordinates` function:
```
coordinatefile = "OB151212coordinates.txt"
satellitedir = "/usr/data/satellitetables"
VBM.SetObjectCoordinates(coordinatefile, satellitedir)
```
The first argument is the name of a file containing the coordinates of the event, i.e. a string with RA and Dec as before (this file is for the event [OGLE-2015-BLG-1212](https://ui.adsabs.harvard.edu/abs/2016ApJ...820...79B/abstract)). The second argument is the path to the directory containing the satellite tables. When the `VBM.SetObjectCoordinates` is executed, the satellite tables are pre-loaded so that they are ready for use in any calculations.

If you want the magnification as seen from satellite 1, then just set VBM.satellite to 1 before the parallax calculation.
```
magnificationspar, y1par, y2par = VBM.BinaryLightCurveParallax(pr,t)      # Calculation of light curve with parallax
VBM.satellite = 1
magnificationssat, y1sat, y2sat = VBM.BinaryLightCurveParallax(pr,t)      # Calculation of light curve seen from Spitzer

plt.plot(t,magnificationspar,"m")
plt.plot(t,magnificationssat,"r")
```
<img src="BinaryLens_lightcurve_satellite.png" width = 400>

The light curve as seen from Spitzer (satellite 1) is in red. And here we have the source trajectories

```
caustics = VBM.Caustics(s,q)
for cau in caustics:
    plt.plot(cau[0],cau[1])
plt.plot(y1par,y2par,"m")
plt.plot(y1sat,y2sat,"r")
```
<img src="BinaryLens_lightcurve_satellite_caustics.png" width = 400>

If you want to return to the ground do not forget to set VBM.satellite back to 0!

[Go to **Orbital motion**](OrbitalMotion.md)
