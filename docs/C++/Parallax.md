[Back to **Light Curve Functions**](LightCurves.md)

# Parallax

In the basic light curve functions discussed in the [Light Curve Functions](LightCurves.md) section, we assumed that observer, lens and source move rectilinearly. Therefore, the relative lens-source angular motion is rectilinear and described by the two parameters $u_0$ and $\alpha$.

In reality, the observer motion is curvilinear as the Earth orbits the Sun. For satellite observations, the situation is similar, since all spacecraft will describe curved orbits. In this section we will introduce a set of functions similar to the "static" ones discussed before, but including the observer motion.

First we will treat geocentric observations and then we will move to observations from space.

## Target coordinates

We need to specify J2000.0 equatorial coordinates for our microlensing target. This can be done by preparing a text file containing R.A. and Dec., similar to the following example:

```
17:57:05 -30:22:59
```

Then we inform VBMicrolensing of these coordinates by the function `SetObjectCoordinates`

```
VBMicrolensing VBM;

VBM.SetObjectCoordinates("OB151212coords.txt",".");
```

The first argument is the name of the file we have just prepared (these are the coordinates of the event [OGLE-2015-BLG-1212](https://ui.adsabs.harvard.edu/abs/2016ApJ...820...79B/abstract)). The second argument will only be used in case of observations from space (see [below](Parallax.md#satellite-parallax)).

## Ephemeris table

In order to calculate the parallax effect, we need to track the Earth position around the Sun. By default, VBMicrolensing uses an ephemeris table (provided in the [/VBMicrolensing/data](/VBMicrolensing/data)) that should be loaded before the first parallax computation. 

```
VBM.LoadSunTable("SunEphemeris.txt");
```

At this point, we are ready for our first calculation.

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
VBMicrolensing VBM; // Declare instance to VBMicrolensing

double pr[9]; // Array of parameters
double s, q, u0, alpha, rho, tE, t0, paiN, paiE, t;

VBM.SetObjectCoordinates("OB151212coords.txt", ".");  // Read target coordinates in file
VBM.LoadSunTable("SunEphemeris.txt"); // Load ephemeris table for parallax calculation

u0 = -0.01; // Impact parameter
t0 = 7550.4; // Time of closest approach to the center of mass
tE = 100.3; // Einstein time
rho = 0.01; // Source radius
s = 0.8; // Separation between the two lenses
q = 0.1; // Mass ratio
alpha = 0.53; // Angle between a vector pointing to the left and the source velocity
paiN = 0.3; // Parallax component in the North direction
paiE = 0.13; // Parallax component in the East direction

// Let us put all parameters in our array
pr[0] = log(s);
pr[1] = log(q);
pr[2] = u0;
pr[3] = alpha;
pr[4] = log(rho);
pr[5] = log(tE);
pr[6] = t0;
pr[7] = paiN;
pr[8] = paiE;

t = 7551.6; // Time at which we want to calculate the magnification

Mag = VBM.BinaryLightCurveParallax(pr, t); // Calculates the Binary Lens magnification at time t with parameters in pr
printf("Binary Light Curve with Parallax at time t: %lf", Mag); // Output should be 31.01...

```

Finally, we add that all light curve functions are available in two versions as explained in [Light Curves](LightCurves.md): the version performing a single calculation of the magnification at time t (as in the example above) and the version calculating the full light curve with one single call (see [Light Curves](LightCurves.md) for details).

## Satellite Parallax

VBMicrolensing can calculate the magnification as seen from a spacecraft. In order to do that, it is necessary to have the ephemerides of the satellite in the format given by the [NASA Horizons system](http://ssd.jpl.nasa.gov/horizons.cgi).

In particular, we assume five columns:
- JD
- RA (degrees)
- Dec (degrees)
- Distance from Earth (AU)
- Distance rate change (not really needed but included by default in Horizons).

Examples of valid satellite ephemerid tables are in [https://github.com/valboz/VBMicrolensing/tree/master/VBMicrolensing/data](https://github.com/valboz/VBMicrolensing/tree/master/VBMicrolensing/data).

The satellite table(s) should be named "satellite*.txt" (with * replaced by a single character). The satellite table files should be in the directory specified as second argument in the `VBM.SetObjectCoordinates` function:

```
VBM.SetObjectCoordinates("OB151212coords.txt","/satellitedir");
```

When the `VBM.SetObjectCoordinates` is executed, the satellite tables are pre-loaded so that they are ready for use in any calculations.

If you want the magnification as seen from satellite 1, then just set VBM.satellite to 1 before the parallax calculation.

```
VBM.SetObjectCoordinates("OB151212coords.txt","/satellitedir"); // Specify coordinates and satellite table directory
VBM.satellite = 1; // All following calculations will be performed as seen from satellite 1 (Spitzer in this example)
Mag = VBM.BinaryLightCurveParallax(pr, t);
printf("Magnification as seen from satellite 1: %lf", Mag); // Output should be 3.88...
```

If you want to return to the ground do not forget to set VBM.satellite back to 0!

If the input time is outside the range of the satellite ephemeris table, `VBM.parallaxextrapolation` is set to 2. This can be used to check that the calculation was performed correctly.

## Parallax system

By default, the parallax components are expressed in the North-East system $(\pi_{E_,N},\pi_{E,E})$. An alternative possibility is to express the parallax vector in the parallel/orthogonal components to the Earth acceleration direction $(\pi_{E,\parallel},\pi_{E,\perp})$. In VBMicrolensing you have both possibilities by setting `VBM.parallaxsystem` to 1 or 0 respectively. The default value is 1, corresponding to the North-East system.

## Reference time for parallax $t_{0,par}$

The parallax effect is introduced as a deviation of the observer from a frame centered on the Earth at a specific reference time $t_{0,par}$, in such a way that the position and the velocity of the source at time $t=t_{0,par}$ remains fixed as seen from the observer. By default, VBMicrolensing uses $t_{0,par}=t_0$, so that the light curve is unchanged at the time of closest approach to the center of mass of the lens. However, if you want to keep the source position at another time fixed, you can set `VBM.t0_par_fixed = 1;` and choose your reference time via `VBM.t0_par`.

## JD vs HJD

When we include parallax, it is important to clarify whether the input time specifications are in JD or HJD. By default, VBMicrolensing assumes that times are given in $HJD' = HJD - 2450000$. However, if you want to calculate a light curve with JD' on your horizontal axis, you should just set `VBM.t_in_HJD = false;` before the execution of the light curve function. All conversions are made by VBMicrolensing internally.

## Implementation of parallax calculations

In order to calculate the parallax effect, we need to track the Earth position around the Sun. We have two possibilities in VBMicrolensing: using an ephemeris table from [Horizons](https://ssd.jpl.nasa.gov/horizons/app.html), or calculate the Earth position solving the Kepler equation from [orbital elements and their secular changes](https://ssd.jpl.nasa.gov/planets/approx_pos.html).

### Ephemeris table

The ephemeris lookup table requires fewer calculations than the Kepler equation and is more accurate. The default ephemeris runs from 1990 to 2050 in steps of one day. If the input time is outside the range of the Sun ephemeris table, `VBM.parallaxextrapolation` is set to 1. This can be used to check that the calculation was performed correctly.
If the user needs a different time window and smaller steps, it is possible to change the ephemeris table to a different file by

```
VBM.LoadSunTable("mySunEphemeris.txt");
```

The file "mySunEphemeris.txt" should be the output of the [Horizons](https://ssd.jpl.nasa.gov/horizons/app.html) system for a geocentric ephemeris of the Sun. The output columns should be JD, RA, Dec, range, rangedot. The [default table file](/VBMicrolensing/data/SunEphemeris.txt) can be taken as reference.

### Kepler's equation calculation

In alternative, VBMicrolensing can calculate ephemeris dynamically for any time. There is no need to pre-load the ephemeris table if you choose this alternative. To follow this route, just add the line

```
VBM.parallaxephemeris = false;
```

The resolution of Kepler equation is slower and retrieves the Earth-Moon barycenter rather than the Earth center. So, it is also less accurate, but provides an alternative reference to check for consistency of parallax calculations. 

## Terrestrial parallax

Different observers on the Earth surface see slightly different microlensing light curves due to the difference in the observation points. In order to keep this difference into account, you should generate different ephemeris tables for each observers. This is possible in [Horizons](https://ssd.jpl.nasa.gov/horizons/app.html) by changing the observer location. Then, you should separately load the right table before each light curve calculation:

```
VBM.LoadSunTable("SunEphemerisfromAfrica.txt");
magAfr = VBM.BinaryLightCurveParallax(pr,t);      // Calculation of light curve with parallax seen from Africa
VBM.LoadSunTable("SunEphemerisfromChile.txt");
magChi = VBM.BinaryLightCurveParallax(pr,t);      // Calculation of light curve with parallax seen from Chile
```

In order to appreciate the differences, you also need to choose very dense sampling in the table generation.


[Go to **Orbital motion**](OrbitalMotion.md)
