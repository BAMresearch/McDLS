# -*- coding: utf-8 -*-
# models/SphericalCoreShell.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from models.scatteringmodel import SASModel
from bases.algorithm import RandomExponential, RandomUniform
from utils.units import Length, SLD

class SphericalCoreShell(SASModel):
    r"""Form factor for a spherical core shell structure
    as defined in the SASfit manual (par. 3.1.4, Spherical Shell III).
    One modification is the ability to specify SLD for core, shell and
    solvent, identical to the notation used in the Core-shell ellipsoid.
    Compared wiht a SASfit-generated model (both with and without distribution)
    """
    shortName = "Core-Shell Sphere"
    parameters = (
            FitParameter("radius", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Core Radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf),
                    activeRange = Length(u'nm').toSi((0.1, 1e3))), # preset
            FitParameter("t", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Thickness of Shell",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf),
                    activeRange = Length(u'nm').toSi((0.1, 1e3))), # preset
            Parameter("eta_c", SLD(u'Å⁻²').toSi(3.16e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Core SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
            Parameter("eta_s", SLD(u'Å⁻²').toSi(2.53e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Shell SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
            Parameter("eta_sol", 0., unit = SLD(u'Å⁻²'),
                    displayName = "Solvent SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
    )

    def __init__(self):
        super(SphericalCoreShell, self).__init__()
        # some presets of parameters to fit
        self.radius.setActive(True)

    def formfactor(self, dataset):
        def k(q, r, dEta):
            # modified K, taken out the volume scaling
            qr = numpy.outer(q, r)
            k = dEta * 3 * (
                    sin(qr) - qr * cos(qr)
                    ) / (qr)**3
            return k

        # dToR = pi / 180. #degrees to radian

        vc = 4./3 * pi *  self.radius() **3
        vt = 4./3 * pi * (self.radius() + self.t()) ** 3
        vRatio = vc / vt

        ks = k(dataset.q, self.radius() + self.t(),
                          self.eta_s() - self.eta_sol())
        kc = k(dataset.q, self.radius(),
                          self.eta_s() - self.eta_c())
        return (ks - vRatio * kc).flatten()

    def volume(self):
        v = 4./3 * pi * (self.radius() + self.t())**3
        return v

    def absVolume(self):
        return self.volume() #TODO: check how to do this.

SphericalCoreShell.factory()

# vim: set ts=4 sts=4 sw=4 tw=0:
