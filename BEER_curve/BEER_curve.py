import numpy as np
import copy
from PyAstronomy.modelSuite.XTran.forTrans import MandelAgolLC
from PyAstronomy.pyasl import isInTransit

__all__ = ['BEER_curve']

class BEER_curve(object):
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted
    components (as well as transit and eclipse signals)
    """

    def __init__(self, time, params, data=None, zero_eclipse_method='mean'):
        """
        Parameters
        ----------
        time : numpy array
            observational time (same units at orbital period)

        params : dict of floats/numpy arrays
            params["per"] - orbital period (any units)
            params["i"] - orbital inclination (degrees); 
                if not given, impact parameter must be
            params["a"] - semi-major axis (stellar radius)
            params["T0"] - mid-transit time (same units as period)
            params["p"] - planet's radius (stellar radius)
            params["LDC"] - numpy array of limb-darkening coefficients
        params["b"] - impact parameter (stellar radius);
            if not given, inclination angle must be
        params["Aellip"] - amplitude of the ellipsoidal variations
        params["Abeam"] - amplitude of the beaming, RV signal
        params["Aplanet"] - amplitude of planet's reflected/emitted signal

        data : numpy array
            observational data (same units as BEER amplitudes)

        zero_eclipse_method : str
            which method to use to shift the data to zero the eclipse
        """

        self.time = time
        self.params = params

        # Orbital phase
        self.phi = self._calc_phi()

        self.ma = MandelAgolLC(orbit="circular", ld="quad")
        # If quadratic limb-darkening
        if(len(params["LDC"]) == 2):
            self.ma["linLimb"] = params["LDC"][0]
            self.ma["quadLimb"] = params["LDC"][1]

        # If non-linear limb-darkening
        elif(len(params["LDC"]) == 4):
            self.ma = MandelAgolLC(orbit="circular", ld="nl")

            self.ma["a1"] = params["LDC"][0]
            self.ma["a2"] = params["LDC"][1]
            self.ma["a3"] = params["LDC"][2]
            self.ma["a4"] = params["LDC"][3]

        self.ma["per"] = params["per"]
        self.ma["i"] = params["i"]
        self.ma["a"] = params["a"]
        self.ma["T0"] = params["T0"]
        self.ma["p"] = params["p"]
        self.ma["b"] = params["b"]

        if(data is not None):
            self.data = data
        if((zero_eclipse_method is not None) and (data is not None)):
            self.data -= self.fit_eclipse_bottom(
                    zero_eclipse_method=zero_eclipse_method)

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

        # Zero out all the limb-darkening coefficients
        if(len(self.params['LDC']) == 2):
            cp['linLimb'] = 0.
            cp['quadLimb'] = 0.
        if(len(self.params['LDC']) == 4):
            cp['a1'] = 0.
            cp['a2'] = 0.
            cp['a3'] = 0.
            cp['a4'] = 0.

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

        self.model_signal = full_signal

        return full_signal

    def transit_duration(self, which_duration="full"):
        """
        Calculates transit duration

        Parameters
        ----------
        which_duration : str 
            "full" - time from first to fourth contact
            "center" - time from contact to contact between planet's center and
                stellar limb
            "short" - time from second to third contact
    
        Returns
        -------
        transit_duration : float
            transit duration in same units as period
        """

        period = self.params['per']
        rp = self.params['p']
        b = self.params['b']
        sma = self.params['a']
    
        if(which_duration == "full"):
            return period/np.pi*np.arcsin(np.sqrt((1. + rp)**2 - b**2)/sma)
        elif(which_duration == "center"):
            return period/np.pi*np.arcsin(np.sqrt(1. - b**2)/sma)
        elif(which_duration == "short"):
            return period/np.pi*np.arcsin(np.sqrt((1. - rp)**2 - b**2)/sma)
        else:
            raise \
                ValueError("which_duration must be 'full', 'center', 'short'!")

    def fit_eclipse_bottom(self, zero_eclipse_method="mean"):
        """
        Calculates the eclipse bottom to set the zero-point in the data

        Parameters
        ----------
        zero_eclipse_method : str
            Which method used to set zero-point - 
                "mean" - Use in-eclipse average value
                "median" - Use in-eclipse median value
        """

        if(zero_eclipse_method == "mean"):
            calc_method = np.nanmean
        elif(zero_eclipse_method == "median"):
            calc_method = np.nanmedian
        else:
            raise ValueError("which_method should be mean or median!")

        # Find in-eclipse points
        TE = self._calc_eclipse_time()
        dur = self.transit_duration(which_duration="short")
        ind = isInTransit(self.time, TE, self.ma["per"], 0.5*dur,\
                boolOutput=False)

        eclipse_bottom = None
        if(ind.size > 0):
            eclipse_bottom = calc_method(self.data[ind])

        return eclipse_bottom

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
