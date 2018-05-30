# BEER_curve
A very small package to model the BEaming, Ellipsoidal variations, and Reflected/emitted light from low-mass companions
with thanks to [Faigler and Mazeh (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.415.3921F) for coming up with a great name.

### Installing
```
pip install BEER_curve
```

BEER_curve uses
[PyAstronomy](https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/index.html) to calculate transits and eclipses. In order to install the transit routines, you'll need also to compile associated transit files occultquad and occultnl. The following procedure works for my computer, running Mac OS 10.13.1 (as of 2018 May 30):

1. Make sure you have a working fortran compiler, such as [gcc](https://stackoverflow.com/questions/9353444/how-to-use-install-gcc-on-mac-os-x-10-8-xcode-4-4).
2. Navigate to PyAstronomy's forTrans directory - 
    1. To figure out where PyAstronomy lives on your computer, in a [terminal window](http://blog.teamtreehouse.com/introduction-to-the-mac-os-x-command-line), type `pip show PyAstronomy` (which assumes you have a working version of [pip](https://conda.io/docs/user-guide/install/index.html)).
    2. At the bottom of the PyAstronomy description that comes up, you should
       see something like:
       ```
       Location: /Users/bjackson/anaconda2/lib/python2.7/site-packages
       ```
       Next type, `cd
       /Users/bjackson/anaconda2/lib/python2.7/site-packages/PyAstronomy/modelSuite/XTran/forTrans`.
       You've made it to the forTrans directory!
3. Run [f2py](https://docs.scipy.org/doc/numpy-1.14.0/f2py/index.html) to
   generate the required .so files:
   * f2py -c occultquad.pyf occultquad.f
   * f2py -c occultnl.pyf occultnl.f
### Example
```
import matplotlib.pyplot as plt
import numpy as np
from BEER_curve import BEER_curve

# HAT-P-7 b parameters from Jackson et al. (2012)
params = {
    "per": 2.204733,
    "i": 83.1,
    "a": 4.15,
    "T0": 0.,
    "p": 1./12.85,
    "linLimb": 0.314709,
    "quadLimb": 0.312125,
    "b": 0.499,
    "Aellip": 37.e-6,
    "Abeam": 5.e-6,
    "F0": 0.,
    "Aplanet": 60.e-6,
    "phase_shift": 0.
    }

t = np.linspace(0, 2*params['per'], 1000)

BC = BEER_curve(t, params)
plt.scatter(t % params['per'], BC.all_signals())
plt.show() # not required if you're running the code in a jupyter notebook

```
