# -*- coding: utf-8 -*-
# models/EllipsoidalCoreShell.py

import scipy, scipy.special
import numpy as np
from numpy import pi, sin, cos, sqrt
from utils.parameter import FitParameter, Parameter
from cutesnake.algorithm import RandomUniform, RandomExponential
from scatteringmodel import ScatteringModel

# parameters must not be inf

class EllipsoidsIsotropic(ScatteringModel):
    r"""Form factor for a spheroidal structure with semi-axes a = b, c.
    c can be set to be an aspect ratio with respect to a
    tested with Ellipsoid II from SASfit 20140626
    """
    shortName = "Isotropic Ellipsoids"
    parameters = (
            FitParameter("a", 1.0,
                    displayName = "Radius of semi-axes a, b",
                    generator = RandomExponential,
                    valueRange = (0., np.inf), suffix = "nm"),
            FitParameter("c", 10.0,
                    displayName = "Radius of semi-axes c",
                    generator = RandomExponential,
                    valueRange = (0., np.inf), suffix = "nm"),
            FitParameter("intDiv", 100,
                    displayName = "Orientation Integration Divisions",
                    generator = RandomUniform,
                    valueRange = (0, 1e4), suffix = "-"),
            FitParameter("cIsAspect", True,
                    displayName = "Radius c value denotes aspect ratio"),
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets, are these still necessary? defined above..
        self.a.setValueRange((0.1, 1e3))
        self.c.setValueRange((1., 1e4))

    def formfactor(self, dataset):
        #From Pedersen, adv. colloid interf. sci. 70 (1997), 171--210

        def rPlugin(Ra, Rc, alpha):
            """calculates replacement R to plug into the Rayleigh sc. func"""
            return np.sqrt(Ra**2 * sin(alpha)**2 + Rc**2 * cos(alpha)**2)

        Ra = self.a()
        q = dataset.q
        if self.cIsAspect():
            Rc = self.a() * self.c()
        else:
            Rc = self.c()

        intVal = np.linspace(0., pi / 2., self.intDiv())
        
        qrP = np.outer(q, rPlugin(Ra, Rc, intVal))
        fsplit = 3.* ( sin(qrP) - qrP * cos (qrP) ) / (qrP**3.) 
        
        # integrate over orientation
        return np.sqrt(np.mean(fsplit**2 * sin(intVal), axis=1)) # should be length q

    def volume(self):
        v = 4./3. * pi * self.a()**2. * self.c()
        return v**self.compensationExponent

EllipsoidsIsotropic.factory()

