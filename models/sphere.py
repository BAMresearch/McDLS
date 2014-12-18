# -*- coding: utf-8 -*-
# models/sphere.py

import numpy
from numpy import pi, sin, cos
from cutesnake.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from sasunit import Length, SLD

class Sphere(ScatteringModel):
    """Form factor of a sphere"""
    shortName = "Sphere"
    parameters = (FitParameter("radius", 1.0,
                    displayName = "Sphere radius",
                    valueRange = (0., numpy.inf),
                    generator = RandomUniform,
                    decimals = 1), 
                  Parameter("sld", 1e14,
                    displayName = "scattering length density difference",
                    valueRange = (0., numpy.inf),
                    decimals = 1), )
    parameters[0].setActive(True)
    # set units
    parameters[0].unit = Length(u'nm')
    parameters[1].unit = SLD(u'Å⁻²')
    #set suffix (normally set in above FitParameter definition) identical
    #to displayname (temporary). Eventually, GUI should use unit metadata
    #now done automatically through definition in ParameterFloat: 
    # parameters[0].setSuffix( parameters[0].unit.displayMagnitudeName )

    def __init__(self):
        super(Sphere, self).__init__()
        #stored in SI units. GUI input must convert upon ingestion
        self.radius.setValueRange((1.0e-10, 1e-6))

    def volume(self):
        result = (pi*4./3.) * self.radius()**(3. * self.compensationExponent)
        return result

    def absVolume(self):
        """
        Volume calculation taking the scattering length density into account
        """
        return self.volume() * self.sld()**2

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
