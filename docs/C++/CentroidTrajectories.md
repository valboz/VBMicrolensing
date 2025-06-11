[Back to **Binary Sources**](BinarySources.md)

# Centroid Trajectories

All basic magnification functions illustrated in [Single Lenses](SingleLenses.md) and [Binary Lenses](BinaryLenses.md) allow the calculation of the centroid of the images in a frame centered on the lens in Einstein units.

However, in order to exploit this information and calculate the centroid positions in the sky during a microlensing event, we need some additional information:

1 - The source heliocentric proper motion in declination and right ascension in mas/yr.

2 - The source parallax in mas.

3 - The Einstein angle in mas.

The lens parallax and proper motion components are obtained from the source parallax and proper motion components once the microlensing parallax components $\pi_N,\pi_E$, Einstein angle and Einstein time are given. Therefore, they do not represent independent parameters.

## Astro-photometric functions

`VBMicrolensing` contains the following astro-photometric functions:

- **PSPLAstroLightCurve(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, np)** <br/>extending PSPLLightCurveParallax

- **ESPLAstroLightCurve(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, np)** <br/>extending ESPLLightCurveParallax

- **BinaryAstroLightCurve(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, np)** <br/>extending BinaryLightCurveParallax

- **BinaryAstroLightCurveOrbital(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, separations_list, np)** <br/>extending BinaryLightCurveOrbital

- **BinaryAstroLightCurveKepler(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, separations_list, np)** <br/>extending BinaryLightCurveKepler

- **BinSourceAstroLightCurveXallarap(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_s1_list, y2_s1_list, y1_s2_list, y2_s2_list, np)** <br/>extending BinSourceExtLightCurveXallarap

- **TripleAstroLightCurve(parameters, times, magnifications, source_centroids_dec, source_centroids_ra, lens_centroids_dec, lens_centroids_ra, y1_list, y2_list, np)** <br/>extending TripleLightCurveParallax

Similarly to their corresponding original functions, these new functions take arrays for the parameters list and the observation times as arguments. The results of the calculations are stored in the arrays provided in the argument list. These include the coordinates for the centroid of the magnified source, the coordinates for the lens and the source position in the reference frame centered in the lens, the separations between orbiting binary lenses. The number of points in the arrays is given as the last parameter np.

Here is a full example with the `PSPLAstroLightCurve`. For the other functions we just have to change the standard parameters accordingly, as explained in the corresponding sections [Light Curves](LightCurves.md), [Parallax](Parallax.md), [Orbital Motion](OrbitalMotion.md), [Binary Sources](BinarySources.md).

```
VBMicrolensing VBM;
const int np = 10000;

# Standard parameters for PSPL with parallax
double t0 = 5034.0;
double tE = 27.0;
double u0 = 0.1;
double paiN = -0.1;
double paiE = +0.2;

# Additional parameters required for centroid trajectory
double muS_Dec = -3.597; // Source proper motion (Dec) in mas/yr
double muS_RA = -2.263; // Source proper motion (RA) in mas/yr
double paiS = 0.12; // Source parallax in mas
double thetaE = 5.15; // Einstein angle in mas

// Here we fill the array of parameters
double pr[] = {u0,log(tE),t0, paiN,paiE,     // Standard light curve parameters for PSPL including parallax
     muS_Dec,muS_RA, paiS, thetaE};      // Additional parameters required for centroid trajectory

// Here we declare all needed arrays
double t[np],mags[np],y1[np],y2[np],c1s[np],c2s[np], c1l[np], c2l[np];

// Let us fill the array of observation epochs
for (int i = 0; i < np; i++) {
	t[i] = -3*365.25 + t0 + 6 * 365.25 / np * i;
}


VBM.SetObjectCoordinates("17:51:40.2082 -29:53:26.502");  // Coordinates of the microlensing event
VBM.astrometry = true; // Do not forget the turn astrometry on!

VBM.PSPLAstroLightCurve(pr, t, mags, c1s, c2s, c1l, c2l, y1, y2, np); // Astro-photometric function call
// Now each array is filled with the corresponding results.
// The centroid for the source is in c1s and c2s, the centroid for the lens in c1l, c2l.
```

## Blending

We may combine lens and source centroid if we know the blending ratio $g = F_L/F_S$ using the function 'CombineCentroids'.

```
double g = 0.4; // Blending
double c1comb[np], c2comb[np]; // Arrays where to store the combined centroid coordinates

VBM.CombineCentroids(mags, c1s, c2s, c1l, c2l, c1comb, c2comb, g, np);
```

## Binary lenses and binary sources

The lens centroid for binary lenses is not the center of mass but the center of light, which depends on the flux ratio, if both components are luminous. The flux ratio between the two lenses is determined using the mass ratio parameter $q$ and applying a mass-luminosity relation: $FR = q^p$, where $p$ can be set by the property `VBM.lens_mass_luminosity_exponent`. The default value is 4. This works also for planetary systems, since the flux ratio would be negligible. If we want to force the secondary to be dark, then we may set `VBM.turn_off_secondary_lens = true`. In this way, all the lens flux will come from the primary.

For binary sources, the centroid is determined using the flux ratio as already explained in [Binary Sources](BinarySources.md).

[Go to **Advanced control**](AdvancedControl.md)
