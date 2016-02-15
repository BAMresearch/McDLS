# -*- coding: utf-8 -*-
# models/dlssphere.py

import numpy
from numpy import pi, exp, sqrt, sin, cos
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter
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
                  )

    def __init__(self):
        super(DLSSphere, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        # square root of the volume to be compatible with SAS-style volume
        # normalization (individual volumes are squared before summing up)
        return sqrt((pi*4./3.) * self.radius()**(3. * self.compensationExponent))

    def formfactor(self, data):
        qr = data.angles * self.radius()
        ff = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.) # usual sphere ff
        res = (exp( data.tauGamma.sanitized / self.radius() ))
        res = data.tauGamma.unflatten(res) * ff
        res = data.tauGamma.flatten(res)
        return res

DLSSphere.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
