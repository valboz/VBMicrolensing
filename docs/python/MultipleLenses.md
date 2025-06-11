[Back to **Binary lenses**](BinaryLenses.md)


# Multiple Lenses

The multiple lens equation reads

$$ \vec{y} = \vec{x} - \sum_{k=1}^N m_k \frac{\vec{x}- \vec{x}_k}{|\vec{x}- \vec{x}_k|^2} $$

Here $N$ is the number of lenses, $\vec{x}_k$ is the position of the k-th lens in the lens plane, $\vec{y}$ is the source position and $m_k$ is the mass of the k-th lens. All angular coordinates are in units of the Einstein radius for a unitary mass. 

## Setting the configuration of the lenses

In VBMicrolensing, before using any functions for multiple lenses, the first step is to initialize the lens configuration using the `SetLensGeometry` function. The `SetLensGeometry` function takes just a list of $3N$ values. For each of the $N$ lenses, we need to specify its position in the lens plane and its mass.

```
VBM = VBMicrolensing.VBMicrolensing()

parameters = [0,0,1,            # First lens: x1_1, x1_2, m1
              1,-0.7,1.e-4,     # Second lens: x2_1, x2_2, m2
              2,0.7,1.e-4,      # Third lens: x3_1, x3_2, m3
              0.6,-.6,1.e-6]    # Fourth lens: x4_1, x4_2, m4

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

```

Once the configuration is specified, we can make all calculations we want with the same configuration. If we want to change the lenses configuration, we need to call `SetLensGeometry` again.

Note that we do not impose that the total mass be 1. We remind that the coordinates are in units of the Einstein radius for a unitary mass.

## Multiple lensing with point sources

For point sources, we can calculate the magnification with the `MultiMag0` function. This depends on the source position $y_1$, $y_2$. Here is an example:

```
VBM = VBMicrolensing.VBMicrolensing()

parameters = [0,0,1,            # First lens: x1_1, x1_2, m1
              1,-0.7,1.e-4,     # Second lens: x2_1, x2_2, m2
              2,0.7,1.e-4,      # Third lens: x3_1, x3_2, m3
              0.6,-.6,1.e-6]    # Fourth lens: x4_1, x4_2, m4

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

y1=0.1 #Source position 
y2=-0.5

Mag=VBM.MultiMag0(y1,y2)   
print(f"Magnification of a point-source = {Mag}\n") #Output should be 2.14....
```

## Multiple lensing with extended sources

For extended sources, our recommended function is `MultiMag2`. This function also takes $\rho$ as an argument, the source radius in units of the Einstein angle:

```
VBM = VBMicrolensing.VBMicrolensing()

parameters = [0,0,1,            # First lens: x1_1, x1_2, m1
              1,-0.7,1.e-4,     # Second lens: x2_1, x2_2, m2
              2,0.7,1.e-4,      # Third lens: x3_1, x3_2, m3
              0.6,-.6,1.e-6]    # Fourth lens: x4_1, x4_2, m4

VBM.SetLensGeometry(parameters) #Initialize the lens configuration

y1=0.  #Source position
y2=-0.2 
rho=0.01 #Source radius

Mag=VBM.MultiMag2(y1,y2,rho) 
print(f"Multiple Lens Magnification = {Mag}\n") #Output should be 5.07...

```

Astrometry for multiple lenses is available with the same notation as for [Binary Lenses](BinaryLenses.md). We just note that  ```VBM.astrox1``` and ```VBM.astrox2``` express the centroid position in the origin as defined by the ```SetLensGeometry``` function.


## Three different algorithms

VBMicrolensing offers three different algorithms for multiple lenses calculations. The preferred algorithm can be selected with the `SetMethod` function:
```
VBM.SetMethod(VBM.Method.Nopoly)
```

The three different alternative methods are available as `Method.Singlepoly`, `Method.Multipoly` and `Method.Nopoly`.

`Method.Singlepoly` solves the lens equation with the classical associated complex polynomial. It suffers from numerical errors, making it inaccurate even with configurations involving three lenses when two lenses are small. It is offered here just as a reference, but it is not intended to be used in any scientific calculations.

`Method.Multipoly` still uses associated complex polynomials but, in order to avoid numerical problems, the reference frame is re-centered on each of the lenses for the calculation of the corresponding images. The computational time is longer, but this is a much safer algorithm. Problems come back only with 5 lenses or more.

`Method.Nopoly` uses a Newton-Raphson method on the lens equation without any manipulations. Nopoly is slower than Multipoly for 3 lenses but is faster with 4 or more lenses, where it proves to be the most reliable and accurate method. This is the default method if the user makes no choice.

The full details of all algorithms are described in [V. Bozza, v. Saggese, G. Covone, P. Rota & J. Zhang, A&A 694 (2025) 219](https://ui.adsabs.harvard.edu/abs/2025A%26A...694A.219B/abstract).

[Go to **Critical curves and caustics**](CriticalCurvesAndCaustics.md)
