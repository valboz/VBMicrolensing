# <span style="color:red">VBMicrolensing</span>

[Back to **Binary lenses**](BinaryLenses.md)

# Critical curves and caustics

The gravitational lensing phenomenology is strictly related to the existence of critical curves and caustics (see the recommended [reviews](reviews.md) for a full treatment).

VBMicroLensing offers the calculation of critical curves and caustics with an arbitrary number of points through the functions ```Caustics``` and ```CriticalCurves```.

The result is a list object, which is a collection of pairs of lists for each curve. Each pair of lists represents the x-coordinate and the y-coordinate of the point.

The use of these objects is very intuitive, as illustrated by this examples:

## Binary Lens

```
import VBMicrolensing
from matplotlib import pyplot as plt

VBML = VBMicrolensing.VBMicrolensing() # Instance to VBMicroLensing

# Parameters of our binary lens
s=2.5;  # separation between the two lenses
q=0.1;  # mass ratio

caustics = VBML.Caustics(s,q)
criticalcurves=VBML.Criticalcurves(s,q)

#plot
fig = plt.figure(figsize=(5, 5))
for cau in caustics:
        plt.plot(cau[0], cau[1], 'k-', markersize=0.1)
plt.xlim(-0.4, 0.1)
plt.ylim(-0.2, 0.2)
plt.title('Caustics')
plt.xlabel('X')
plt.ylabel('Y')
plt.minorticks_on()
plt.tick_params(axis='both', which='major', width=2, length=10, direction='in', bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', which='minor', width=1, length=5, direction='in', bottom=True, top=True, left=True, right=True)

```
<img src="Caustics_binary.png" width = 600>


## Multi Lens

```
import VBMicrolensing
from matplotlib import pyplot as plt

VBML = VBMicrolensing.VBMicrolensing() # Instance to VBMicroLensing

#parameters of the system
#s = [(real position, imaginary position, lens mass), (...), ...]

s=[(0,0,1),(-1.2,0.5,0.5),(-1,0.4,1.1e-1),(0.6,0,1.1e-2)]

s_real = [item[0] for item in s]
s_im = [item[1] for item in s]
s_q = [item[2] for item in s]

VBML.SetLensGeometry(s_q,s_real,s_im)

caustics = VBML.Multicaustics()
criticalcurves=VBML.Multicriticalcurves()

#plot
fig = plt.figure(figsize=(5, 5))
for cau in caustics:
        plt.plot(cau[0], cau[1], 'k-', markersize=0.1)
plt.xlim(-1.5, 1)
plt.ylim(-1, 1.5)
plt.title('Caustics')
plt.xlabel('X')
plt.ylabel('Y')
plt.minorticks_on()
plt.tick_params(axis='both', which='major', width=2, length=10, direction='in', bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='both', which='minor', width=1, length=5, direction='in', bottom=True, top=True, left=True, right=True)

```
<img src="Caustics_multi.png" width = 600>


Critical curves and caustics are calculated through the resolution of a fourth order complex polynomial (see [reviews](reviews.md)) by the [Skowron & Gould algorithm](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/). 

The **number of points** calculated for the critical curves is controlled by ```VBML.NPcrit```, which can be changed by the user according to the desired sampling. The default value is 200.

[Go to: **Limb Darkening**](LimbDarkening.md)
