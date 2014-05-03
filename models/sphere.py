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
        ScatteringModel.__init__(self)
        # this only works for people
        # defining lengths in angstrom or nm, not m.
        self.radius.setValueRange((1.0, 1e4))

    def volume(self, paramValues):
        r = numpy.array((self.radius(),))
        if self.radius.isActive():
            r = paramValues[:, 0]
        result = (pi*4./3.) * r**(3. * self.compensationExponent)
        return result

    def formfactor(self, dataset, paramValues):
        r = numpy.array((self.radius(),))
        if self.radius.isActive():
            r = paramValues[:, 0]
        q = dataset.q
        qr = numpy.outer(q, r)
        result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
        return result

Sphere.factory()

# vim: set ts=4 sts=4 sw=4 tw=0:
