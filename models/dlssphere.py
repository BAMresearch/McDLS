# -*- coding: utf-8 -*-
# models/dlssphere.py

import numpy
from numpy import pi, exp, sqrt
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
        return sqrt((pi*4./3.) * self.radius()**(3. * self.compensationExponent))

    def formfactor(self, data):
        return (exp( data.tauGamma.sanitized / self.radius() ))

DLSSphere.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
