import numpy as np
import copy
import batman

__all__ = ['BEER_curve', 'all_signals_func']


class BEER_curve(object):
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted
    components (as well as transit and eclipse signals)
    """

    def __init__(self, time, params, supersample_factor=1, exp_time=0):
        """
        Parameters
        ----------
        time : numpy array
            observational time (same units at orbital period)

        params : dict of floats/numpy arrays
            params["per"] - orbital period (any units)
            params["a"] - semi-major axis (in units of stellar radius)
            params["b"] - impact parameter (in units of stellar radius)
            params["p"] - planetary radius (in stellar radii)
            params["T0"] - mid-transit time
            params["baseline"] - photometric baseline
            params["Aellip"] - amplitude of the ellipsoidal variations
            params["eclipse_depth"] - depth of eclipse

            params["Abeam"] - amplitude of the beaming, RV signal
            params["Aplanet"] - amplitude of planet's reflected/emitted signal
            params["phase_shift"] - phase shift of planet's signal

            OR

            params["Asin"] - amplitude of sine term
            params["Acos"] - amplitude of cosine term

            params["A3"] - amplitude of third harmonic
            params["theta3"] - phase offset for third harmonic

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

    def _sine_term(self):
        """
        Returns Asin*sin(2*pi*phi)
        """

        Asin = self.params['Asin']
        phi = self.phi

        return Asin*np.sin(2.*np.pi*phi)

    def _cosine_term(self):
        """
        Returns -Acos*cos(2*pi*phi)
        """

        Acos = self.params['Acos']
        phi = self.phi

        return -Acos*np.cos(2.*np.pi*phi)

    def _third_harmonic(self):
        """
        Returns third harmonic
        """
        A3 = self.params['A3']
        theta3 = self.params['theta3']
        phi = self.phi

        return A3*np.cos(3.*2.*np.pi*(phi - theta3))

    def _eclipse(self):
        """
        Calculates eclipse signal
        """
        ma = batman.TransitParams()

        ma.t0 = self._calc_eclipse_time()
        ma.per = self.params['per']
        ma.rp = np.sqrt(np.abs(self.params['eclipse_depth']))
        ma.a = self.params['a']
        ma.inc = np.arccos(self.params['b']/self.params['a'])*180./np.pi
        ma.ecc = 0. # Assume zero eccentricity for now
        ma.w = 90. #long of periastron (in degrees); irrelevant for e = 0
        ma.u = [0., 0.] # uniform disk for eclipse
        ma.limb_dark = "quadratic"       #limb darkening model

        m = batman.TransitModel(ma, self.time_supersample) # initializes model
        flux = m.light_curve(ma) # calculates light curve

        return (flux - 1.)*np.sign(self.params['eclipse_depth'])

    def _calc_eclipse_time(self):
        """
        Returns eclipse time --
          For now (2018 Oct 31), assumes circular orbit
        """

        return self.params['T0'] + 0.5*self.params['per']

    def all_signals(self):
        """
        Calculates BEER curves
        """

        time_supersample = self.time_supersample
        time = self.time

        baseline = self.params["baseline"]
        E = self._ellipsoidal_curve()
        eclipse = self._eclipse()

        if(('Aplanet' in self.params.keys()) & ('Abeam' in self.params.keys())):
            Be = self._beaming_curve()
            R = self._reflected_emitted_curve()
        else:
            Be = self._sine_term()
            R = self._cosine_term()

        full_signal = baseline + Be + E + R + eclipse
        if('A3' in self.params):
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

def all_signals_func(time, per, a, b, p, T0, baseline, Aellip, Abeam, Aplanet,
        phase_shift, eclipse_depth, supersample_factor=1., exp_time=0.):
    """
    A standalone function version of the BEER_curve all_signals method

    Parameters
    ----------
    time : numpy array
        observational time (same units at orbital period)

    per - orbital period (any units)
    a - semi-major axis (in units of stellar radius)
    b - impact parameter (in units of stellar radius)
    p - planetary radius (in stellar radii)
    T0 - mid-transit time
    baseline - photometric baseline
    Aellip - amplitude of the ellipsoidal variations
    Abeam - amplitude of the beaming, RV signal
    Aplanet - amplitude of planet's reflected/emitted signal
    phase_shift - phase shift of planet's signal
    eclipse_depth - depth of eclipse

    supersample_factor : (optional) int
        Number of points subdividing exposure

    exp_time : (optional) float
        Exposure time (in same units as `t`)

    Returns
    -------
    The BEER curve signals, including the eclipse but NOT the transit

    """

    params = {'per': per, 'a': a, 'b': b, 'p': p, 'T0': T0,
            'baseline': baseline, 'Aellip': Aellip, 'Abeam': Abeam, 
            'Aplanet': Aplanet, 'phase_shift': phase_shift, 
            'eclipse_depth': eclipse_depth}

    BC = BEER_curve(time, params, 
            supersample_factor=supersample_factor, exp_time=exp_time)

    return BC.all_signals()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    # HAT-P-7 b parameters from Jackson et al. (2012)
    params = {
            "per": 2.204733,
            "a": 4.15,
            "b": 4.15*np.cos(83.1/180.*np.pi),
            "p": 1./12.85,
            "T0": 0.,
            "baseline": 0.,
            "Aellip": 37.e-6,
            "Abeam": 5.e-6,
            "Aplanet": 60.e-6,
            "phase_shift": 0.01,
            "eclipse_depth": 60.e-6
            }

    t = np.linspace(0, 2*params['per'], 1000)

    BC = BEER_curve(t, params)

    plt.plot(t % params['per'], BC.all_signals(), 'bo')

    # And also the sine and cosine versions
    params = {
            "per": 2.204733,
            "a": 4.15,
            "b": 4.15*np.cos(83.1/180.*np.pi),
            "p": 1./12.85,
            "T0": 0.,
            "baseline": 0.,
            "Aellip": 37.e-6,
            "Asin": 5.e-6 - 60.e-6*np.sin(2.*np.pi*0.01), 
            "Acos": 60.e-6*np.cos(2.*np.pi*0.01),
            "eclipse_depth": 60.e-6
            }

    BC = BEER_curve(t, params)

    plt.plot(t % params['per'], BC.all_signals(), 'r.')
    plt.show(block=True)

