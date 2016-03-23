# -*- coding: utf-8 -*-
# models/sphere.py

import numpy
from numpy import pi, sin, cos
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import SASModel
from dataobj import SASData
from utils.units import Length, NM, SLD

class Sphere(SASModel):
    """Form factor of a sphere"""
    shortName = "Sphere"
    parameters = (FitParameter("radius",
                    NM.toSi(10.), unit = NM,
                    displayName = "Sphere radius",
                    valueRange = (0., numpy.inf),
                    activeRange = NM.toSi((1., 1000.)),
                    generator = RandomUniform,
                    decimals = 1), 
                  Parameter("sld", SLD(u'Å⁻²').toSi(1e-6), unit = SLD(u'Å⁻²'),
                    displayName = "scattering length density difference",
                    valueRange = (0., numpy.inf),
                    decimals = 1), )

    def __init__(self):
        super(Sphere, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        r"""Calculates the volume of a sphere defined by:

        :math:`v(r) = {4\pi \over 3} r^3`
        """
        result = (pi*4./3.) * self.radius()**3
        return result

    def absVolume(self):
        r"""Calculates the volume of a sphere taking the scattering length
        density difference :math:`\Delta\rho` into account:

        :math:`v_{abs}(r, \Delta\rho) = v_{sph}(r) \cdot \Delta\rho^2`
        """
        return self.volume() * self.sld()**2

    def formfactor(self, dataset):
        r"""Calculates the form factor of a sphere defined by:

        :math:`F(q, r) = { 3 ~ sin(qr) - qr \cdot cos(qr) \over (qr)^3 }`
        """
        if isinstance(dataset, SASData):
            q = dataset.q
        else:
            q = dataset
        qr = q * self.radius() 
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
