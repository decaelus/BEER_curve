# BEER_curve
A very small package to model the BEaming, Ellipsoidal variations, and Reflected/emitted light from low-mass companions
with thanks to [Faigler and Mazeh (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.415.3921F) for coming up with a great name.

### Installing
```
pip install BEER_curve
```
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
    "baseline": 0.,
    "Aellip": 37.e-6,
    "Abeam": 5.e-6,
    "F0": 0.,
    "Aplanet": 60.e-6,
    "eclipse_depth": 60.e-6,
    "phase_shift": 0.
    }

t = np.linspace(0, 2*params['per'], 1000)

BC = BEER_curve(t, params)
plt.scatter(t % params['per'], BC.all_signals())
plt.show() # not required if you're running the code in a jupyter notebook

```
