import numpy as np
from batman import *

__all__ = ['BEER_params', 'BEER_curve']

class BEER_params(TransitParams):
    """
    Object in which to store the physical parameters of the system.

    :param Aellip: amplitude (unitless) of the ellipsoidal variations
    :type u: float 

    :param Abeam: amplitude (unitless) of the beaming, RV signal
    :type u: float 

    :param Aplanet: amplitude (unitless) of planet's reflected/emitted signal
    :type u: float 

    """
    def __init__(self, **kwargs):
        super(TransitParams, self).__init__(**kwargs)
        
        #Setting my own parameters
        self.Aellip = None
        self.Abeam = None
        self.Aplanet = None

class BEER_curve(object):
    """
        Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted
    """

    def planetary_curve(phi, A_refl):
        y = - (A_refl*cos(2.*pi*phi))
        y = y + abs(min(y))  
        return y
