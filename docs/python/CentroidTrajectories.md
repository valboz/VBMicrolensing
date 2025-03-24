[Back to **Binary Sources**](BinarySources.md)

# Centroid Trajectories

We have seen that all basic magnification functions illustrated in [Single Lenses](SingleLenses.md) and [Binary Lenses](BinaryLenses.md) allow the calculation of the centroid of the images in a frame centered on the lens in Einstein units.

However, in order to exploit this information and calculate the centroid positions in the sky all along a microlensing event, we need some additional information:
1 - The lens offset in the sky at time t0 in declination and right ascension from the reference coordinates (these can be set to zero, but are useful if you need to fit an astrometric series).
2 - The source heliocentric proper motion in declination and right ascension in mas/yr.
3 - The source parallax in mas.
4 - The Einstein angle in mas.

The lens parallax and proper motion components are obtained from the source parallax and proper motion components once the microlensing parallax components $\pi_N,\pi_E$ and the Einstein angle are given. Therefore, they do not represent independent parameters.

`VBMicrolensing` contains the following astro-photometric functions:

```
PSPLAstroLightCurve
ESPLAstroLightCurve
BinaryAstroLightCurve
BinaryAstroLightCurveOrbital
BinaryAstroLightCurveKepler
BinSourceAstroLightCurveXallarap
```



[Go to **Advanced control**](AdvancedControl.md)
