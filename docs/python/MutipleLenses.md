# <span style="color:red">VBMicrolensing</span>

[Back to **Binary lenses**](SingleLenses.md)


# Multiple Lenses

The multiple lens equation reads

$$z_s=z-\sum_{k=1}^{N} \frac{m_k}{\overline{z} - \overline{z_k}}$$

Here $N$ is the number of lenses, $z_k$ is the position of the k-th lens in the complex plane, $z_s$ is the source position and $m_k$ is the mass fracrion of the k-th lens. All angular coordinates are in units of the total Einstein radius, defined in terms of the total mass of the system.

In VBMicrolensing, when we want to use functions that model multiple lenses, the first step is to initialize the lens configuration using the `SetLensGeometry` function and chose the method with `SetMethod`.

The `SetLensGeometry` function takes three different arguments: the masses of the lenses, the real positions of the lenses, and the imaginary positions of the lenses. 

In VBMicrolensign we have three different methods: `Singlepoly`, `SetLensGeometry` and `SetLensGeometry`.
Singlepoly solve the lens equation with the classical associate polynomia, it suffers from numerical errors, making it inaccurate even with configurations involving three lenses.
Multipoly use classical associate polynomia, to overcome the problem of the singlepoly, re-center the polynomial on each lens to find nearby roots.
Nopoly method uses an expansion of the lens equation, a Newton-Raphson method. The advantage of Nopoly is its speed compared to Multipoly. However, there is a risk that it may fail to find some of the images. We recommend using Nopoly and, in case of doubtful results, switching to Multipoly.
For full details, refer to the paper (currently in preparation).

Here is an example of how to initialize the lens configuration:

```
import VBMicrolensing

VBML = VBMicrolensing.VBMicrolensing()

s = [(0,0,1),(1,-0.7,1.e-4),(2,0.7,1.e-4),(0.6,-.6,1.e-6)] 

s_1 = [item[0] for item in s]
s_2 = [item[1] for item in s]
q = [item[2] for item in s]

VBML.SetMethod(VBML.Method.Multipoly) #Singlepoly, Multipoly, Nopoly

VBML.SetLensGeometry(q,s_1,s_2)

```

## Multiple lensing with point sources

For point sources, we can get the magnification with the `MultiMag0` function. This depends on  the source position $y_1$, $y_2$. Here is an example:

```
import VBMicrolensing

VBML = VBMicrolensing.VBMicrolensing()

s = [(0,0,1),(1,-0.7,1.e-4),(2,0.7,1.e-4),(0.6,-.6,1.e-6)] 

s_1 = [item[0] for item in s]
s_2 = [item[1] for item in s]
q = [item[2] for item in s]

VBML.SetMethod(VBML.Method.Multipoly)

VBML.SetLensGeometry(q,s_1,s_2)

y_1=1
y_2=-0.7

Mag=VBML.MultiMag0(y_1,y_2)
```

## Multiple lensing with extended sources

For extended sources, the function is `MultiMag`. This function also depends on $\rho$, the source radius in units of the total Einstein angle:

```
import VBMicrolensing

VBML = VBMicrolensing.VBMicrolensing()

s = [(0,0,1),(1,-0.7,1.e-4),(2,0.7,1.e-4),(0.6,-.6,1.e-6)] 

s_1 = [item[0] for item in s]
s_2 = [item[1] for item in s]
q = [item[2] for item in s]
rho=0.01

VBML.SetMethod(VBML.Method.Nopoly)

VBML.SetLensGeometry(q,s_1,s_2)

y_1=1.7
y_2=-0.7

Mag=VBML.MultiMag(y_1,y_2,rho)
print(Mag)

```

[Go to **Critical curves and caustics**](CriticalCurvesAndCaustics.md)