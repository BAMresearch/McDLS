# -*- coding: utf-8 -*-
# models/dlssphere.py

import numpy
from numpy import pi, exp, sqrt, sin, cos
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import DLSModel
from utils.units import NM

class DLSSphere(DLSModel):
    """Sphere Form factor for DLS"""
    shortName = "Sphere (DLS)"
    parameters = (FitParameter("radius",
                    NM.toSi(100.), unit = NM,
                    displayName = "Hydrodynamic Radius",
                    valueRange = (0., numpy.inf),
                    activeRange = NM.toSi((1., 1000.)),
                    generator = RandomUniform,
                    decimals = 1), 
                  Parameter("withFF", True,
                    displayName = "Include the SAS sphere form factor?"),
                  )

    def __init__(self):
        super(DLSSphere, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        return (pi*4./3.) * self.radius()**3

    def _ffSphere(self):
        if not self.withFF():
            return 1.
        qr = self.scatteringVector * self.radius()
        return 3. * (sin(qr) - qr * cos(qr)) / (qr**3.) # usual sphere ff

    # how to get adjustable size sensitivity by compensationEXponent?
    # -> the SAXS way is by default not possible due to the normalized model == 1
    # idea: using a weight w/o compExp for normalization but with compExp in the model
    #   -> violates the specs but may lead to desired result ...
    def weight(self):
        # a compExp. != 1 just changes the volume of a scatterer
        return self._ffSphere() * self.volume()**self.compensationExponent

    def formfactor(self, data):
        return exp( data.tauGamma.sanitized / self.radius() )

DLSSphere.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
