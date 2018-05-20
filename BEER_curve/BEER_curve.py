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

        # Orbital phase
        self.phi = self._calc_phi()

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

    def _calc_phi(self):
        """
        Calculates orbital phase
        """
        time = self.time
        T0 = self.params['T0']
        per = self.params['per']

        return ((time - T0) % per)/per

    def _calc_eclipse_time(self):
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
        
        F0 = self.params['F0']
        Aplanet = self.params['Aplanet']
        phase_shift = self.params['phase_shift']

        phi = self.phi

        return F0 - Aplanet*np.cos(2.*np.pi*(phi - phase_shift))

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

        return ma.evaluate(time)

    def eclipse(self):
        """
        Uses PyAstronomy's transit light curve routine with uniform
        limb to calculate eclipse
        """

        time = self.time
        ma = self.ma
        TE = self._calc_eclipse_time()
        eclipse_depth = self.params['F0'] + self.params['Aplanet']

        # Make a copy of ma but set limb-darkening parameters to zero for
        # uniform disk
        cp = copy.deepcopy(ma)
        cp['linLimb'] = 0.
        cp['quadLimb'] = 0.
        cp['T0'] = TE
        cp['p'] = np.sqrt(eclipse_depth)

        return cp.evaluate(time)

    def all_signals(self):
        """
        Calculates BEER curves, as well as transit and eclipse signals
        """

        transit = self.transit() - 1.
        eclipse = self.eclipse()

        # Shift so middle of eclipse sits at zero and reflected curve
        # disappers in eclipse
        scaled_eclipse = eclipse
        scaled_eclipse[eclipse == np.min(eclipse)] = 0.

        Be = self.beaming_curve()
        
        E = self.ellipsoidal_curve()
        E -= np.min(E)

        R = self.reflected_emitted_curve()

        full_signal = transit + Be + E + R*scaled_eclipse

        return full_signal


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # HAT-P-7 b parameters from Jackson et al. (2013)
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
    plt.show(block=True)
