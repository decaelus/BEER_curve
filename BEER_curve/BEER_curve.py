import numpy as np
import copy
from PyAstronomy.pyasl import isInTransit

__all__ = ['BEER_curve']


class BEER_curve(object):
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted
    components (as well as transit and eclipse signals)
    """

    def __init__(self, time, params, data=None, third_harmonic=False,
            supersample_factor=1, exp_time=0):
        """
        Parameters
        ----------
        time : numpy array
            observational time (same units at orbital period)

        params : dict of floats/numpy arrays
            params["per"] - orbital period (any units)
            params["T0"] - mid-transit time
            params["baseline"] - photometric baseline
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

        self.third_harmonic = third_harmonic

        if(data is not None):
            self.data = data

    def _calc_phi(self):
        """
        Calculates orbital phase
        """
        time = self.time_supersample
        T0 = self.params['T0']
        per = self.params['per']

        return ((time - T0) % per)/per

    def _reflected_emitted_curve(self):
        """
        Calculates planet's reflected/emitted component, i.e. R in BEER
        """
        
        Aplanet = self.params['Aplanet']
        phase_shift = self.params['phase_shift']

        phi = self.phi

        return -Aplanet*np.cos(2.*np.pi*(phi - phase_shift))

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

    def all_signals(self):
        """
        Calculates BEER curves
        """

        time_supersample = self.time_supersample
        time = self.time

        baseline = self.params["baseline"]
        Be = self._beaming_curve()
        E = self._ellipsoidal_curve()
        R = self._reflected_emitted_curve()

        full_signal = baseline + Be + E + R
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

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # HAT-P-7 b parameters from Jackson et al. (2013)
    params = {
            "per": 2.204733,
            "T0": 0.,
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
