# -*- coding: utf-8 -*-
# models/cylinders.py

"""This model is a special case of the radially isotropic cylinders, where
the cylinders are also slightly tilted out of the plane parallel to the
detector. This out of plane tilt is described using a Gaussian. The
integration over this tilt angle is done over several segments of the Gaussian
PDF, with each segment occupying an equal cumulative probability. The centroid
value used for the integration is the mass-weighted centre.
"""

import numpy, scipy, scipy.special, scipy.stats
from numpy import pi, zeros, sin, cos, linspace, diff, sinc
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class CylindersRadiallyIsotropicTilted(ScatteringModel):
    r"""Form factor of cylinders *UNFINISHED*
    which are radially isotropic (so not spherically isotropic!)
    """
    shortName = "Cylinders defined by aspect ratio"
    parameters = (
            FitParameter("radius", 1.0,
                    displayName = "Cylinder radius",
                    valueRange = (0.1, numpy.inf),
                    activeRange = (0.1, 1e3),
                    suffix = "nm"),
            FitParameter("aspect", 10.0,
                    displayName = "Aspect ratio L/(2R) of the cylinder",
                    valueRange = (0.1, numpy.inf),
                    activeRange = (1, 20),
                    suffix = "-"),
            FitParameter("psiAngle", 0.1,
                    displayName = "in-plane cylinder rotation",
                    valueRange = (0.1, 180.1), suffix = "deg."),
            FitParameter("psiAngleDivisions", 303.,
                    displayName = "in-plane angle divisions",
                    valueRange = (1, numpy.inf), suffix = "-"),
            FitParameter("phiDistWidth", 10.0,
                    displayName = "out-of-plane axis distribution width",
                    #with 90 degrees fully out-of-plane (parallel to beam):
                    valueRange = (0.1, 90.1), suffix = "deg."),
            FitParameter("phiDistDivisions", 9.,
                    displayName = "out of plane integration divisions",
                    valueRange = (1, numpy.inf), suffix = "-"),
    )
    parameters[0].setActive(True)
    parameters[1].setActive(False) # not expected to vary
    parameters[2].setActive(False) # keep out for now.
    parameters[3].setActive(False) # not expected to vary
    parameters[4].setActive(False) # gaussian width
    parameters[5].setActive(False) # integration steps

    def __init__(self):
        super(CylindersRadiallyIsotropicTilted, self).__init__()
        # some presets

    def formfactor(self, dataset):
        #the remaining values are never active fitting parameters
        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #also shown in Figure 3.98 in the SASfit manual
        #used in the equation for a cylinder from Pedersen, 1997

        dToR = pi/180. #degrees to radian
        psiRange = self.psiAngle.valueRange()
        psi = numpy.linspace(psiRange[0], psiRange[1], self.psiAngleDivisions())
        # determine the divisions to integrate phi over: equal probability
        x = linspace(0, 0.99, self.phiDistDivisions() + 1.) # zero added
        phiAngles = scipy.stats.norm.interval(x)[1] # only get the positive values
        # now calculate the centroid value for each integration bin
        phiCtr = scipy.stats.norm.interval(x[:-1] + diff(x)/2.)[1]
        # each of these integration bins is equally probable. Can integrate
        # using mean function.

        # rotation can be used to get slightly better results, but
        # ONLY FOR RADIAL SYMMETRY, NOT SPHERICAL.
        fcyl = 0.
        for pIdx in range(len(phiCtr)):
            # TODO: Implementing from equations 3.263 in SASfit manual
            # leave the cylinder axis arbitrary psi rotation out of it for now.
            # calculate cos(gamma)
            # since we do not want to describe theta, we use an angle
            # omega describing cylinder axis projection in x-z plane (=psi),
            # combined with
            # phi describing cylinder axis projection onto x-y plane
            # theta=omega/cos(phi) (probably)
            # cosGammaP=sin(psi*dToR)*cos(psi*dToR)*cos(phiCtr[pidX]*dToR)\
            #         + cos(psi*dToR)*sin(psi*dToR)
            # cosGammaM=
            qRsina = numpy.outer(dataset.q, self.radius() * sin(psi * dToR))
            # approximation for small tilts, adjusts the length of the cylinder only!
            qLcosa = numpy.outer(dataset.q, self.radius() * self.aspect()
                                 * cos(phiCtr[pIdx] * dToR) * cos(psi * dToR))

            fsplit = (2. * scipy.special.j1(qRsina)/qRsina * sinc(qLcosa / pi))
            # integrate over orientation
            fcyl += numpy.sqrt(numpy.mean(fsplit**2, axis=1)) / len(phiCtr) # should be length q

        return fcyl

    def volume(self):
        v = pi * self.radius()**2 * (2. * self.radius() * self.aspect())
        return v**self.compensationExponent

CylindersRadiallyIsotropicTilted.factory()

#if __name__ == "__main__":
#    from bases.datafile import PDHFile, AsciiFile
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
#    model = CylindersRadiallyIsotropicTilted()
#    model.radius.setValue(1.)
#    model.radius.setActive(False)
#    model.aspect.setValue(100.)
#    model.aspect.setActive(False)
#    model.psiAngle.setValue(1.)
#    model.psiAngle.setActive(False)
#    model.psiAngleDivisions.setValue(303)
#    model.psiAngleDivisions.setActive(False)
#    intensity = model.formfactor(pf.data, None).reshape(-1)
#    q = pf.data[:, 0]
#    oldInt = pf.data[:, 1]
#    delta = abs(oldInt - intensity)
#    result = numpy.dstack((q, intensity, delta))[0]
#    AsciiFile.writeFile("CylindersRadiallyIsotropic.dat", result)
#    # call it like this:
#    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
