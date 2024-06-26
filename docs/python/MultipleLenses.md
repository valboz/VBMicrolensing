[Back to **Binary lenses**](BinaryLenses.md)


# Multiple Lenses

The multiple lens equation reads

$$ \vec{y} = \vec{x} - \sum_{k=1}^N m_k \frac{\vec{x}- \vec{x}_k}{|\vec{x}- \vec{x}_k|^2}

$$z_s=z-\sum_{k=1}^{N} \frac{m_k}{\overline{z} - \overline{z_k}}$$

Here $N$ is the number of lenses, $z_k$ is the position of the k-th lens in the complex plane, $z_s$ is the source position and $m_k$ is the mass fracrion of the k-th lens. All angular coordinates are in units of the total Einstein radius, defined in terms of the total mass of the system.

In VBMicrolensing, when we want to use functions that model multiple lenses, the first step is to initialize the lens configuration using the `SetLensGeometry` function and chose the method with `SetMethod`.

The `SetLensGeometry` function takes three different arguments: the masses of the lenses, the real positions of the lenses, and the imaginary positions of the lenses. 

In VBMicrolensign when selecting the method with`SetMethod`, we have three different alternatives: `Singlepoly`, `SetLensGeometry` and `SetLensGeometry`.
Singlepoly solve the lens equation with the classical associate polynomia, it suffers from numerical errors, making it inaccurate even with configurations involving three lenses.
Multipoly use classical associate polynomia, to overcome the problem of the singlepoly, re-center the polynomial on each lens to find nearby roots.
Nopoly method uses an expansion of the lens equation, a Newton-Raphson method. The advantage of Nopoly is its speed compared to Multipoly. However, there is a risk that it may fail to find some of the images.
We recommend using Nopoly and, in case of doubtful results, switching to Multipoly.
For full details, refer to the paper (V.Bozza, V.Saggese et al. ,currently in preparation).

Here is an example of how to initialize the lens configuration:

```
VBM = VBMicrolensing.VBMicrolensing()

parameters = [0,0,1,1,-0.7,1.e-4,2,0.7,1.e-4,0.6,-.6,1.e-6] 
#Within the parameter list are the lens details, as follows: [z1_re, z1_im, q1, z2_re, z2_im, q2,....,zn_re,zn_im,qn]

VBM.SetMethod(VBM.Method.Multipoly) #Singlepoly, Multipoly, Nopoly

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

```

## Multiple lensing with point sources

For point sources, we can get the magnification with the `MultiMag0` function. This depends on  the source position $y_1$, $y_2$. Here is an example:

```
VBM = VBMicrolensing.VBMicrolensing()

parameters = [0,0,1,1,-0.7,1.e-4,2,0.7,1.e-4,0.6,-.6,1.e-6] 
#Within the parameter list are the lens details, as follows: [z1_re, z1_im, q1, z2_re, z2_im, q2,....,zn_re,zn_im,qn]

VBM.SetMethod(VBM.Method.Nopoly) #Singlepoly, Multipoly, Nopoly

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

y=(0.1,-0.5) #Source position 

Mag=VBM.MultiMag0(y)  
print(f"Magnification of a point-source = {Mag}\n") #Output should be 2.14....
```

## Multiple lensing with extended sources

For extended sources, the function is `MultiMag`. This function also takes $\rho$ as an argument, the source radius in units of the Einstein angle:

```
VBM = VBMicrolensing.VBMicrolensing()

s = [(0,0,1),(1,-0.7,1.e-4),(2,0.7,1.e-4),(0.6,-.6,1.e-6)] 

parameters = [0,0,1,1,-0.7,1.e-4,2,0.7,1.e-4,0.6,-.6,1.e-6] 
#Within the parameter list are the lens details, as follows: [z1_re, z1_im, q1, z2_re, z2_im, q2,....,zn_re,zn_im,qn]

VBM.SetMethod(VBM.Method.Nopoly) #Singlepoly, Multipoly, Nopoly

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

y=(0.,-0.2) #Source position 
rho=0.01 #source radius

Mag=VBM.MultiMag(y,rho) 
print(f"Multiple Lens Magnification = {Mag}\n") #Output should be 5.07...

```

[Go to **Critical curves and caustics**](CriticalCurvesAndCaustics.md)
