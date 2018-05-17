import numpy as np
import copy
from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC

__all__ = ['BEER_curve']


class BEER_curve(object):
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted
    components (as well as transit and eclipse signals)

    :param Aellip: amplitude (unitless) of the ellipsoidal variations
    :type u: float

    :param Abeam: amplitude (unitless) of the beaming, RV signal
    :type u: float

    :param Aplanet: amplitude (unitless) of planet's reflected/emitted signal
    :type u: float

    """

    def __init__(self, time, params):

        self.time = time
        self.params = params

        self.phi = self.__calc_phi__()

        # Just circular orbits and quadratic LD for now
        self.ma = MandelAgolLC(orbit="circular", ld="quad")

        self.ma["per"] = params["per"]
        self.ma["i"] = params["i"]
        self.ma["a"] = params["a"]
        self.ma["T0"] = params["T0"]
        self.ma["p"] = params["p"]
        self.ma["linLimb"] = params["linLimb"]
        self.ma["quadLimb"] = params["quadLimb"]
        self.ma["b"] = params["b"]

    def __calc_phi__(self):
        """
        Calculates orbital phase
        """
        time = self.time
        T0 = self.params['T0']
        per = self.params['per']

        return ((time - T0) % per)/per

    def __calc_eclipse_time__(self):
        """
        Calculates mid-eclipse time -
        I've included this function here in anticipation of using eccentric
        orbits in the near future.
        """

        T0 = self.params['T0']
        per = self.params['per']

        return T0 + 0.5*per

    def reflected_emitted_curve(self):
        """
        Calculates planet's reflected/emitted component, i.e. R in BEER
        """
        Aplanet = self.params['Aplanet']
        phi = self.phi

        return -(Aplanet*np.cos(2.*np.pi*phi))

    def beaming_curve(self):
        """
        Calculates the beaming effect curve
        """
        Abeam = self.params['Abeam']
        phi = self.phi

        return Abeam*np.sin(2.*np.pi*phi)

    def ellipsoidal_curve(self):
        """
        Calculates the ellipsoidal variation curve
        """
        Aellip = self.params['Aellip']
        phi = self.phi

        return -Aellip*np.cos(2.*2.*np.pi*phi)

    def all_BEER_curves(self):
        """
        Calculates all BEER curves
        """
        return self.reflected_emitted_curve() +\
            self.beaming_curve() +\
            self.ellipsoidal_curve()

    def transit(self):
        """
        Uses PyAstronomy's quadratic limb-darkening routine to calculate
        the transit light curve
        """

        time = self.time
        ma = self.ma

        return ma.evaluate(time) - 1.0

    def eclipse(self):
        """
        Uses PyAstronomy's transit light curve routine with uniform
        limb to calculate eclipse
        """

        time = self.time
        ma = self.ma
        TE = self.__calc_eclipse_time__()
        eclipse_depth = self.params['eclipse_depth']

        # Make a copy of ma but set limb-darkening parameters to zero for
        # uniform disk
        cp = copy.deepcopy(ma)
        cp['linLimb'] = 0.
        cp['quadLimb'] = 0.
        cp['T0'] = TE
        cp['p'] = np.sqrt(eclipse_depth)

        # Shift so baseline out of transit and eclipse is zero
        return cp.evaluate(time) - 1.0

    def transit_and_eclipse(self):
        """
        Calculates both transit and eclipse signals
        """

        return self.transit() + self.eclipse()

    def all_signals(self):
        """
        Calculates BEER curves, as well as transit and eclipse signals
        """

        return self.transit_and_eclipse() + self.all_BEER_curves()


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # HAT-P-7 b parameters from Jackson et al. (2013)
    params = {
            "per": 2.204733,
            "i": 83.1,
            "a": 4.15,
            "T0": 0.,
            "eclipse_depth": 60.e-6,
            "p": 1./12.85,
            "linLimb": 0.314709,
            "quadLimb": 0.312125,
            "b": 0.499,
            "Aellip": 37.e-6,
            "Abeam": 5.e-6,
            "Aplanet": 60.e-6
            }

    t = np.linspace(0, 2*params['per'], 1000)

    BC = BEER_curve(t, params)

    plt.scatter(t % params['per'], BC.all_signals())
    plt.show(block=True)
