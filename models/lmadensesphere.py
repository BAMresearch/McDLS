# -*- coding: utf-8 -*-
# models/lmadensesphere.py

import numpy as np
from numpy import pi, sin, cos
import logging
from cutesnake.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from sasunit import Length, Fraction, NoUnit, SLD

class LMADenseSphere(ScatteringModel):
    """Form factor of a sphere convoluted with a structure factor,
    equations 15-17 from Pedersen, J. Appl. Cryst. 27 (1994), 595--608. 
    Correct eqn given in Kinning and Thomas, Macromolecules 17 (1984) 1712.
    Internally set parameters are volume fraction of the hard spheres,
    and the multiplication factor /mf/ for an additional stand-off distance
    between the hard spheres: Rh=mf*R where Rh is the hard-sphere radius
    ("interaction radius") used in the structure factor, R is the radius of
    the sphere, and mf is the multiplication factor.
    """

    shortName = "LMADenseSphere"
    parameters = (
            FitParameter("radius", 1.0e-9, unit = Length(u"nm"),
                    displayName = "Sphere radius",
                    valueRange = (1e-10, 1e1),
                    generator = RandomUniform,
                    decimals = 1),
            FitParameter("volFrac", 0.1, unit = Fraction(u"%"),
                    displayName = "Volume fraction of spheres",
                    valueRange = (0, 1.),
                    generator = RandomUniform,
                    decimals = 1),
            FitParameter("mf", -1., # auto
                    displayName = "standoff multiplier (-1 = auto)",
                    valueRange = (-1., 400.),
                    generator = RandomUniform,
                    decimals = 1,
                    displayValues = {-1.: "auto"}),
            Parameter("sld", 1e14, unit = SLD(u'Å⁻²'),
                    displayName = "Scattering length density difference",
                    valueRange = (0, 1e15))
            )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)

    def volume(self):
        result = (pi*4./3.) * self.radius()**(3. * self.compensationExponent)
        return result
    
    def absVolume(self):
        return self.volume() * self.sld()**2

    def formfactor(self, dataset):
        
        SFmu = self.volFrac()
        SFmf = self.mf()
        if SFmf == -1:
            SFmf = (0.634 / SFmu) **(1. / 3)

        def SFG(A, SFmu):
            alpha = ( 1 + 2 * SFmu )**2 / ( 1 - SFmu )**4
            beta = -6 * SFmu * ( 1 + SFmu / 2 )**2 / ( 1 - SFmu )**4
            gamma = SFmu * alpha / 2
            G = (
                    alpha * ( sin(A) - A * cos(A) ) / A**2 +
                    beta *(2 * A * sin(A) + (2 - A**2) * cos(A) - 2) / A**3 +
                    gamma*( -1 * A**4 * cos(A) + 4 * ((3 * A**2 - 6) * cos(A)
                        + (A**3 - 6 * A) * sin(A) + 6 ) ) / A**5
                    )
            return G


        qr = dataset.q * self.radius()
        result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
        #now we introduce the structure factor
        rhsq = 2. * dataset.q * (SFmf * self.radius()) 
        G = SFG(rhsq, SFmu)
        S = (( 1. + 24. * SFmu * G / rhsq ))**(-1)
        #print (S < 0).sum()
        # the above structure factor needs to be the square root as it is
        # taken into the form factor. Eventually, we want to calculate the
        # scattering as FF**2 * S, which we are now achieving as
        # (FF * S**0.5)**2. The minus sign in the exponent of the above
        # equation comes from the original equation in the literature.
        result = np.sqrt(result**2 * S)
        return result

LMADenseSphere.factory()

# vim: set ts=4 sts=4 sw=4 tw=0:
