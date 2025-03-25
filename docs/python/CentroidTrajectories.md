[Back to **Binary Sources**](BinarySources.md)

# Centroid Trajectories

All basic magnification functions illustrated in [Single Lenses](SingleLenses.md) and [Binary Lenses](BinaryLenses.md) allow the calculation of the centroid of the images in a frame centered on the lens in Einstein units.

However, in order to exploit this information and calculate the centroid positions in the sky during a microlensing event, we need some additional information:

1 - The lens offset in the sky at time t0 in declination and right ascension from the reference coordinates (these can be set to zero, but are useful if you need to fit an astrometric series).

2 - The source heliocentric proper motion in declination and right ascension in mas/yr.

3 - The source parallax in mas.

4 - The Einstein angle in mas.

The lens parallax and proper motion components are obtained from the source parallax and proper motion components once the microlensing parallax components $\pi_N,\pi_E$, Einstein angle and Einstein time are given. Therefore, they do not represent independent parameters.

`VBMicrolensing` contains the following astro-photometric functions:

- **PSPLAstroLightCurve(parameters, times)**, extending PSPLLightCurveParallax <br/>
returns `[magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_list, y2_list]`

`ESPLAstroLightCurve(parameters, times)              # extending ESPLLightCurveParallax
- returns [magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_list, y2_list]

BinaryAstroLightCurve(parameters, times)            # extending BinaryLightCurveParallax
- returns [magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_list, y2_list]

BinaryAstroLightCurveOrbital(parameters, times)     # extending BinaryLightCurveOrbital
- returns [magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_list, y2_list, separations_list]

BinaryAstroLightCurveKepler(parameters, times)      # extending BinaryLightCurveKepler
- returns [magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_list, y2_list, separations_list]

BinSourceAstroLightCurveXallarap(parameters, times)   # extending BinSourceExtLightCurveXallarap
- returns [magnifications, source_centroid_dec, source_centroid_ra, lens_centroid_dec, lens_centroid_ra, y1_s1_list, y2_s1_list, y1_s2_list, y2_s2_list]
```

Similarly to their corresponding original functions, these new functions take a parameters list and a list of observation times as arguments. The output contains a list of magnifications calculated at the epochs in `times`, centroid positions in (dec,ra) for source and lens, source positions in the lens reference frame.

Here is a full example with the `PSPLAstroLightCurve`:






[Go to **Advanced control**](AdvancedControl.md)
