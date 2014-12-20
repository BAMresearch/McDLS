# -*- coding: utf-8 -*-
# models/EllipsoidalCoreShell.py

import scipy, scipy.special
import numpy as np
from numpy import pi, sin, cos, sqrt
from utils.parameter import FitParameter, Parameter
from bases.algorithm import RandomUniform, RandomExponential
from scatteringmodel import ScatteringModel
from utils.units import Length, NoUnit, SLD

# parameters must not be inf

class EllipsoidsIsotropic(ScatteringModel):
    r"""Form factor for a spheroidal structure with semi-axes a = b, c.
    c can be set to be an aspect ratio with respect to a
    tested with Ellipsoid II from SASfit 20140626
    """
    shortName = "Isotropic Ellipsoids"
    parameters = (
            FitParameter("a", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Radius of semi-axes a, b",
                    generator = RandomExponential,
                    valueRange = Length(u'nm').toSi((0.1, 1e10)),
                    activeRange = Length(u'nm').toSi((0.1, 1e3))), # preset
            Parameter("useAspect", True,
                    displayName = "Use aspect ratio (checked) or length to define c-axis"),
            FitParameter("c", Length(u'nm').toSi(10.), unit = Length(u'nm'),
                    displayName = "Radius of semi-axes c",
                    generator = RandomExponential,
                    valueRange = Length(u'nm').toSi((0.1, 1e10)),
                    activeRange = Length(u'nm').toSi((1., 1e4))), # preset
            FitParameter("aspect", 10.0,
                    displayName = "aspect ratio of semi-axes c to a, b",
                    generator = RandomExponential,
                    valueRange = (1e-3, 1e3)),
            Parameter("intDiv", 100,
                    displayName = "Orientation Integration Divisions",
                    valueRange = (0, 1e4)),
            Parameter("sld", SLD(u'Å⁻²').toSi(1e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Scattering length density difference",
                    valueRange = (0, SLD(u'Å⁻²').toSi(1e-16))),
    )

    def __init__(self):
        super(EllipsoidsIsotropic, self).__init__()
        # some presets of parameters to fit
        self.a.setActive(True)

    def formfactor(self, dataset):
        #From Pedersen, adv. colloid interf. sci. 70 (1997), 171--210

        def rPlugin(Ra, Rc, alpha):
            """calculates replacement R to plug into the Rayleigh sc. func"""
            return np.sqrt(Ra**2 * sin(alpha)**2 + Rc**2 * cos(alpha)**2)

        Ra = self.a()
        q = dataset.q
        if self.useAspect():
            Rc = self.a() * self.aspect()
        else:
            Rc = self.c()

        intVal = np.linspace(0., pi / 2., self.intDiv())
        
        qrP = np.outer(q, rPlugin(Ra, Rc, intVal))
        fsplit = 3.* ( sin(qrP) - qrP * cos (qrP) ) / (qrP**3.) 
        
        # integrate over orientation
        return np.sqrt(np.mean(fsplit**2 * sin(intVal), axis=1)) # should be length q

    def volume(self):
        Ra = self.a()
        if self.useAspect():
            Rc = self.a() * self.aspect()
        else:
            Rc = self.c()

        v = 4./3. * pi * Ra**2. * Rc
        return v**self.compensationExponent

    def absVolume(self):
        return self.volume() * self.sld()**2

EllipsoidsIsotropic.factory()

