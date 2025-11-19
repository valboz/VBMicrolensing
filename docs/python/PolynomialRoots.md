# Polynomial Roots

Most of the calculations in `VBMicrolensing` require the resolution of complex polynomials by [Skowron & Gould algorithm](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/). For the user convenience, we also offer a function for calculation of roots of a generic polynomial by this algorithm:

```
VBM = VBMicrolensing.VBMicrolensing()

coefficients = [(1.0,0.0), (-2.0,0.0), (1.0,0.0)]

roots = VBM.cmplx_roots_gen(coefficients) 
print(roots) # Output should be [[1.0,0.0],[1.0,0.0]]
```

The list of coefficients is an array of t-uples representing real and imaginary part of a complex number. The polynomial to be solved is then
```
coefficients[0] z^0 + coefficients[1] z^1 + coefficients[2] z^2 + ...
```

The roots are returned in a similar form, as a list of t-uples with real and imaginary part of each root.
