# Polynomial Roots

Most of the calculations in `VBMicrolensing` require the resolution of complex polynomials by 

A function for calculation of roots of a generic polynomial by this algorithm is also offered:

```
VBM = VBMicrolensing.VBMicrolensing()

s = 0.8 # separation between the two lenses
q = 0.1 # mass ratio
y1 = 0.01 # y1 is the source coordinate along the axis parallel to the line joining the two lenses 
y2 = 0.01 # y2 is the source coordinate orthogonal to the first one

Mag = VBM.BinaryMag0(s, q, y1, y2) # Call to the BinaryMag0 function with these parameters
print(f"Magnification of a point-source = {Mag}\n") # Output should be 18.18.....

```
