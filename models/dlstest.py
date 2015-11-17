# -*- coding: utf-8 -*-
# models/dlstest.py

import numpy
from numpy import pi, exp, sqrt
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter
from scatteringmodel import DLSModel
from utils.units import NM

class DLSTest(DLSModel):
    """Test Form factor for DLS"""
    shortName = "DLS-Test"
    parameters = (FitParameter("radius",
                    NM.toSi(100.), unit = NM,
                    displayName = "Hydrodynamic Radius",
                    valueRange = (0., numpy.inf),
                    activeRange = NM.toSi((1., 1000.)),
                    generator = RandomUniform,
                    decimals = 1), 
                  )

    def __init__(self):
        super(DLSTest, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        return sqrt((pi*4./3.) * self.radius()**3.)

    def formfactor(self, data):
        return (exp( data.tauGamma.value / self.radius() ))

DLSTest.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
