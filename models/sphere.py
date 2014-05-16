# -*- coding: utf-8 -*-
# models/sphere.py

import numpy
from numpy import pi, sin, cos
from cutesnake.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel

class Sphere(ScatteringModel):
    """Form factor of a sphere"""
    shortName = "Sphere"
    parameters = (FitParameter("radius", 1.0,
                    displayName = "Sphere radius",
                    valueRange = (0., numpy.inf),
                    generator = RandomUniform,
                    suffix = "nm", decimals = 1), )
    parameters[0].setActive(True)

    def __init__(self):
        super(Sphere, self).__init__()
        # this only works for people
        # defining lengths in angstrom or nm, not m.
        self.radius.setValueRange((1.0, 1e4))

    def volume(self):
        result = (pi*4./3.) * self.radius()**(3. * self.compensationExponent)
        return result

    def formfactor(self, dataset):
        qr = dataset.q * self.radius()
        result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
        return result

Sphere.factory()

# see GaussianChain for some notes on this
def test():
    Sphere.testRelErr = 1e-4
    for fn in ("sasfit_sphere-2-1.dat",
               "sasfit_sphere-10-1.dat",
               "sasfit_sphere-20-1.dat",
               "sasfit_sphere-50-1.dat",
               "sasfit_sphere-100-1.dat"):
        yield Sphere.test, fn

# vim: set ts=4 sts=4 sw=4 tw=0:
