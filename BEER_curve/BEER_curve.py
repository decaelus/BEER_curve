import numpy as np
from PyAstronomy.modelSuite import forTrans as ft

__all__ = ['BEER_params', 'BEER_curve']

class BEER_curve():
    """
    Calculates the BEaming, Ellipsoidal variation, and Reflected/emitted components

    :param Aellip: amplitude (unitless) of the ellipsoidal variations
    :type u: float

    :param Abeam: amplitude (unitless) of the beaming, RV signal
    :type u: float

    :param Aplanet: amplitude (unitless) of planet's reflected/emitted signal
    :type u: float

    """

    def __init__(self, params, t):
        #Setting my own parameters
        self.Aellip = params.Aellip
        self.Abeam = params.Abeam
        self.Aplanet = params.Aplanet

#   def planetary_curve(phi, A_refl):
#       y = -(A_refl*cos(2.*pi*phi))
#       y = y + abs(min(y))  
#       return y

if __name__ == "__main__":
    BC = BEER_curve()

    print(BC.params)
