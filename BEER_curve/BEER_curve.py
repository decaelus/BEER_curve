import numpy as np
from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC

__all__ = ['BEER_curve']

class BEER_curve(object):
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted components

    :param Aellip: amplitude (unitless) of the ellipsoidal variations
    :type u: float

    :param Abeam: amplitude (unitless) of the beaming, RV signal
    :type u: float

    :param Aplanet: amplitude (unitless) of planet's reflected/emitted signal
    :type u: float

    """

    def __init__(self, t, params):

        self.t = t
        self.params = params

        #Just circular orbits and quadratic LD for now
        self.ma = MandelAgolLC(orbit="circular", ld="quad")

    def _calc_phi(self):
        """
        Calculates orbital phase
        """

        return ((self.t - self.params['T0']) % self.params['per'])/self.params['per']

    def planetary_curve(self):
        """
        Calculates planet's reflected/emitted component, i.e. R in BEER
        """
        Aplanet = self.params['Aplanet']
        phi = self._calc_phi()

        return -(Aplanet*np.cos(2.*np.pi*phi))

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    #HAT-P-7 b parameters from Jackson et al. (2013)
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
            "Aplanet": 60.e-6 
            }

    t = np.linspace(0, 2*params['per'])

    BC = BEER_curve(t, params)

    plt.plot(t, BC.planetary_curve())
    plt.show(block=True)
