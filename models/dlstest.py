# -*- coding: utf-8 -*-
# models/dlstest.py

import numpy
from numpy import pi, sin, cos, exp
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from utils.units import Length, NM, K, VIS, NoUnit

# Boltzmann constant in m²·kg·s⁻²·K⁻¹
KB = 1.38064852 * 1e23

class DLSTest(ScatteringModel):
    """Test Form factor for DLS"""
    shortName = "DLS-Test"
    parameters = (FitParameter("radius",
                    NM.toSi(100.), unit = NM,
                    displayName = "Hydrodynamic Radius",
                    valueRange = (0., numpy.inf),
                    activeRange = NM.toSi((1., 1000.)),
                    generator = RandomUniform,
                    decimals = 1), 
                  Parameter("temp", K.toSi(300.), unit = K,
                    displayName = K.name(),
                    valueRange = (0., numpy.inf),
                    decimals = 3),
                  Parameter("vis", VIS.toSi(1.), unit = VIS,
                    displayName = VIS.name(),
                    valueRange = (0., numpy.inf),
                    decimals = 3),
                  Parameter("refIdx", 1., unit = NoUnit,
                    displayName = "Refractive Index",
                    valueRange = (0., numpy.inf),
                    decimals = 3),
                  Parameter("wavelength", NM.toSi(600.), unit = NM,
                    displayName = "Wavelength",
                    valueRange = (0., numpy.inf),
                    decimals = 1),
                  )

    def __init__(self):
        super(DLSTest, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        v = (pi*4./3.) * self.radius()**(3 * .5) # gets squared in parent scope
        return v

    def formfactor(self, dataset):
        q = 4. * pi * self.refIdx() * sin(theta * .5) / self.wavelength()
        gamma = - q * q * self.temp() * KB / (6. * pi * self.vis() * self.radius())
        return exp(dataset.q * gamma * .5)

Sphere.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
