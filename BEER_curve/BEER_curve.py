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

    def __init__(self, time, params, data=None, third_harmonic=False,
            zero_eclipse_method='mean', supersample_factor=1, exp_time=0):
        """
        Parameters
        ----------
        time : numpy array
            observational time (same units at orbital period)

        params : dict of floats/numpy arrays
            params["per"] - orbital period (any units)
            params["a"] - semi-major axis (stellar radius)
            params["T0"] - mid-transit time (same units as period)
            params["p"] - planet's radius (stellar radius)
            if quadratic limb-darkening,
                params["linLimb"] - linear limb-darkening coefficient
                params["quadLimb"] - quadratic limb-darkening coefficient
            if non-linear limb-darkening,
                params["a*"] - limb-darkening coefficients
            params["b"] - impact parameter (stellar radius);
                if not given, inclination angle must be
            params["Aellip"] - amplitude of the ellipsoidal variations
            params["Abeam"] - amplitude of the beaming, RV signal
            params["Aplanet"] - amplitude of planet's reflected/emitted signal

            if(third_harmonic):
            params["A3"] - amplitude of third harmonic
            params["theta3"] - phase offset for third harmonic

        data : numpy array
            observational data (same units as BEER amplitudes)

        third_harmonic : bool
            Whether or not to include a third harmonic in phase curve

        zero_eclipse_method : str
            which method to use to shift the data to zero the eclipse

        supersample_factor : int
            Number of points subdividing exposure

        exp_time : float
            Exposure time (in same units as `t`)
        """

        # Based on Kreidberg's batman approach
        self.supersample_factor = supersample_factor
        self.time = time
        self.exp_time = exp_time
        if self.supersample_factor > 1:
            t_offsets = np.linspace(-self.exp_time/2., 
                    self.exp_time/2., 
                    self.supersample_factor)
            self.time_supersample = (t_offsets +\
                    self.time.reshape(self.time.size, 1)).flatten()

        else: self.time_supersample = self.time

        self.params = params

        # Orbital phase
        self.phi = self._calc_phi()

        self.ma = MandelAgolLC(orbit="circular", ld="quad")
        # If quadratic limb-darkening
        if("linLimb" in params):
            self.ma["linLimb"] = params["linLimb"]
            self.ma["quadLimb"] = params["quadLimb"]

        # If non-linear limb-darkening
        elif("a1" in params):
            self.ma = MandelAgolLC(orbit="circular", ld="nl")

            self.ma["a1"] = params["a1"]
            self.ma["a2"] = params["a2"]
            self.ma["a3"] = params["a3"]
            self.ma["a4"] = params["a4"]

        self.ma["per"] = params["per"]
        # Set using the impact parameter
        self.ma["i"] = 180./np.pi*np.arccos(params["b"]/params["a"])
        self.ma["a"] = params["a"]
        self.ma["T0"] = params["T0"]
        self.ma["p"] = params["p"]

        self.third_harmonic = third_harmonic

        if(data is not None):
            self.data = data
        if((zero_eclipse_method is not None) and (data is not None)):
            self.data -= self.fit_eclipse_bottom(
                    zero_eclipse_method=zero_eclipse_method
                    )

    def _calc_phi(self):
        """
        Calculates orbital phase
        """
        time = self.time_supersample
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

    def _reflected_emitted_curve(self):
        """
        Calculates planet's reflected/emitted component, i.e. R in BEER
        """
        
        F0 = self.params['F0']
        Aplanet = self.params['Aplanet']
        phase_shift = self.params['phase_shift']

        phi = self.phi

        return F0 - Aplanet*np.cos(2.*np.pi*(phi - phase_shift))

    def _beaming_curve(self):
        """
        Calculates the beaming effect curve
        """
        Abeam = self.params['Abeam']
        phi = self.phi

        return Abeam*np.sin(2.*np.pi*phi)

    def _ellipsoidal_curve(self):
        """
        Calculates the ellipsoidal variation curve
        """
        Aellip = self.params['Aellip']
        phi = self.phi

        return -Aellip*np.cos(2.*2.*np.pi*phi)

    def _third_harmonic(self):
        """
        Returns third harmonic relevant to HAT-P-7 and KOI-13
        """
        A3 = self.params['A3']
        theta3 = self.params['theta3']
        phi = self.phi

        return A3*np.cos(3.*2.*np.pi*(phi - theta3))

    def _transit(self):
        """
        Uses PyAstronomy's quadratic limb-darkening routine to calculate
        the transit light curve
        """

        time = self.time_supersample
        ma = self.ma

        return ma.evaluate(time)

    def _eclipse(self):
        """
        Uses PyAstronomy's transit light curve routine with uniform
        limb to calculate eclipse
        """

        time_supersample = self.time_supersample
        ma = self.ma
        TE = self._calc_eclipse_time()
        eclipse_depth = self.params["F0"] + self.params["Aplanet"]

        # Make a copy of ma but set limb-darkening parameters to zero for
        # uniform disk
        cp = copy.deepcopy(ma)

        # Zero out all the limb-darkening coefficients
        if("linLimb" in self.params):
            cp["linLimb"] = 0.
            cp["quadLimb"] = 0.
        if("a1" in self.params):
            cp["a1"] = 0.
            cp["a2"] = 0.
            cp["a3"] = 0.
            cp["a4"] = 0.

        cp["T0"] = TE
        cp["p"] = np.sqrt(eclipse_depth)

        return cp.evaluate(time_supersample)

    def all_signals(self):
        """
        Calculates BEER curves, as well as transit and eclipse signals
        """

        time_supersample = self.time_supersample
        time = self.time

        transit = self._transit() - 1.
        eclipse = self._eclipse()

        # Shift so middle of eclipse sits at zero and reflected curve
        # disappers in eclipse
        scaled_eclipse = eclipse
        scaled_eclipse[eclipse == np.min(eclipse)] = 0.

        Be = self._beaming_curve()
        
        E = self._ellipsoidal_curve()
        E -= np.min(E)

        R = self._reflected_emitted_curve()

        full_signal = transit + Be + E + R*scaled_eclipse
        if(self.third_harmonic):
            full_signal += self._third_harmonic()

        self.model_signal = full_signal

        if(self.supersample_factor > 1): 
            self.model_signal =\
                    np.mean(full_signal.reshape(-1, self.supersample_factor),\
                    axis=1)
            full_signal =\
                    np.mean(full_signal.reshape(-1, self.supersample_factor),\
                    axis=1)

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
        time = self.time
        period = self.ma["per"]
        TE = self._calc_eclipse_time()
        dur = self.transit_duration(which_duration="short")
        ind = isInTransit(time, TE, period, 0.5*dur, boolOutput=True)

        eclipse_bottom = 0.
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
