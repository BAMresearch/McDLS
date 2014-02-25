# -*- coding: utf-8 -*-
# models/lmadensesphere.py

import numpy
from numpy import pi, sin, cos
import logging
from cutesnake.algorithm import RandomUniform
from utils.parameter import Parameter
from scatteringmodel import ScatteringModel

class LMADenseSphere(ScatteringModel):
    """Form factor of a sphere convoluted with a structure factor, 
    equations 15-17 from Pedersen, J. Appl. Cryst. 27 (1994), 595--608
    Internally set parameters are volume fraction of the hard spheres,
    and the multiplication factor /mf/ for an additional stand-off distance
    between the hard spheres: Rh=mf*R where Rh is the hard-sphere radius 
    ("interaction radius") used in the structure factor, R is the radius of 
    the sphere, and mf is the multiplication factor. 
    """

    shortName = "LMADenseSphere"
    parameters = (Parameter("radius", 1.0,
                    displayName = "Sphere radius",
                    valueRange = (0., numpy.inf),
                    generator = RandomUniform,
                    suffix = "nm", decimals = 1), )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        self.radius.setValueRange((1.0, 1e4)) #this only works for people
        #defining lengths in angstrom or nm, not m.

    def updateParamBounds(self, bounds):
        bounds = ScatteringModel.updateParamBounds(self, bounds)
        if len(bounds) < 1:
            return
        if len(bounds) == 1:
            logging.warning("Only one bound provided, "
                            "assuming it denotes the maximum.")
            bounds.insert(0, self.radius.valueRange(0))
        elif len(bounds) > 2:
            bounds = bounds[0:2]
        logging.info("Updating lower and upper contribution parameter bounds "
                     "to: ({0}, {1}).".format(min(bounds), max(bounds))) 
        #logging.info changed from bounds[0] and bounds[1] to reflect better 
        #what is done below:
        self.radius.setValueRange((min(bounds), max(bounds)))

    def vol(self, paramValues, compensationExponent = None):
        assert ScatteringModel.vol(self, paramValues)
        if compensationExponent is None:
            compensationExponent = self.compensationExponent
        result = (pi*4./3.) * paramValues**(3. * compensationExponent)
        return result

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)
        
        structureFactorVolumeFraction = 0.2
        structureFactorStandoff = (0.634/structureFactorVolumeFraction)**(1./3)
        #no idea if this will work
        SFmu=structureFactorVolumeFraction
        SFmf=structureFactorStandoff
        def SFG(A,SFmu):
            alpha = ( 1 + 2 * SFmu )**2 / ( 1 - SFmu )**4
            beta = -6 * SFmu * ( 1 + SFmu / 2 )**2 / ( 1 - SFmu )**2
            gamma = SFmu * alpha / 2
            G = ( 
                    alpha * ( sin(A) - A * cos(A) ) / A**2 +
                    beta *(2 * A * sin(A) + (2 - A**2) * cos(A) - 2) / A**3 +
                    gamma*( -1 * A**4 * cos(A) + 4 * ((3 * A**2 - 6) * cos(A)
                        + (A**3 - 6 * A) * sin(A) + 6 ) ) / A**5 
                    )
            return G


        r = paramValues.flatten()
        q = dataset[:, 0]
        qr = numpy.outer(q, r)
        result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
        #now we introduce the structure factor
        rhsq = numpy.outer(q,(SFmf * r))
        G = SFG(rhsq,SFmu)
        Ssqrt = (abs( 1. + 24. * SFmu * G / rhsq ))**(-0.5)
        # the above structure factor needs to be the square root as it is
        # taken into the form factor. Eventually, we want to calculate the
        # scattering as FF**2 * S, which we are now achieving as
        # (FF * S**0.5)**2. The minus sign in the exponent of the above
        # equation comes from the original equation in the literature.
        result = result * Ssqrt
        return result

LMADenseSphere.factory()

# vim: set ts=4 sts=4 sw=4 tw=0:
