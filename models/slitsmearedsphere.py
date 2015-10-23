# -*- coding: utf-8 -*-
# models/sphere.py

import numpy
import numpy as np
from numpy import pi, sin, cos
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter, Parameter
from scatteringmodel import SASModel
from utils.units import Length, NM, SLD, NoUnit, ScatteringVector

class slitSmearedSphere(SASModel):
    """Form factor of a sphere, smeared for slit-collimated instruments. 
    TEST USE ONLY"""
    shortName = "Slit-smeared Sphere"
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
                    decimals = 1), 
                  Parameter("slitWidth", ScatteringVector(u'nm⁻¹').toSi(-1.), 
                    unit = ScatteringVector(u'nm⁻¹'),
                    displayName = "slit width in Q [-1 is no smearing]",
                    valueRange = (-1., numpy.inf),
                    displayValues = {-1.e-9: "no smearing"},
                    decimals = 1), 
                  Parameter("smearingSteps", 25,
                    displayName = "number of smearing steps",
                    valueRange = (0, 1000),
                    decimals = 0), 
                  )

    def __init__(self):
        super(slitSmearedSphere, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        result = (pi*4./3.) * self.radius()**(3. * self.compensationExponent)
        return result

    def absVolume(self):
        """
        Volume calculation taking the scattering length density into account
        """
        return self.volume() * self.sld()**2

    def formfactor(self, dataset):
        def baseCalc(q, radius = None):
            assert (radius is not None)
            qr = q * radius
            result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
            return result

        def smear(q, fhandle, fparams = None, slitWidth = None, nIntSteps = 50):
            """ clean smearing program for obtaining slit-smeared scattering patterns. 
            usage: 
            *q*: A one-dimensional scattering vector
            *fhandle*: A function handle to a scattering pattern calculator. Must accept q
            *fparams*: A dictionary of keyword-parameter sets passed on to the calculator
            *slitwidth*: The width of the slit in units of q
            """
            assert isinstance(q, numpy.ndarray)
            assert (q.ndim == 1)
            # nIntSteps = q.size
            if slitWidth > q.max(): #TODO: use clip function
                slitWidth = q.max()
            dU = np.logspace(np.log10(q.min() / 10.), 
                    np.log10(slitWidth / 2.), num = nIntSteps) 
            dU = np.concatenate(([0,], dU)) [np.newaxis, :] 
            # assuming a square weighting function. 
            weightFunc = np.ones((dU.size,)) / slitWidth
            # calculate the intensities at sqrt(q**2 + dU **2)
            locs = np.sqrt((q[:,np.newaxis] + 0 * dU)**2 + (0 * q[:, np.newaxis] + dU)**2)
            if fparams is None:
                fVal = fhandle(locs) * (0 * q[:, np.newaxis] + weightFunc)
            else:
                fVal = fhandle(locs, **fparams) * (0 * q[:, np.newaxis] + weightFunc)
            sI = 2 * np.trapz(fVal, x = dU, axis = 1) # * 2.
            return sI

        if self.slitWidth() > 0.:
            return smear(dataset.q, baseCalc, 
                    fparams = {"radius": self.radius()}, 
                    slitWidth = self.slitWidth(), 
                    nIntSteps = self.smearingSteps())
        else:
            return baseCalc(dataset.q, self.radius())
        
            

slitSmearedSphere.factory()

# see GaussianChain for some notes on this
def test():
    slitSmearedSphere.testRelErr = 1e-4
    for fn in ("sasfit_sphere-2-1.dat",
               "sasfit_sphere-10-1.dat",
               "sasfit_sphere-20-1.dat",
               "sasfit_sphere-50-1.dat",
               "sasfit_sphere-100-1.dat"):
        yield Sphere.test, fn

# vim: set ts=4 sts=4 sw=4 tw=0:
