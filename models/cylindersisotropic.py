# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from cutesnake.algorithm import RandomExponential, RandomUniform

# parameters must not be inf

class CylindersIsotropic(ScatteringModel):
    r"""Form factor of cylinders
    previous version (length-fixed) checked against SASfit
    """
    shortName = "Isotropic Cylinders"
    parameters = (
            FitParameter("radius", 1.0,
                    displayName = "Cylinder Radius",
                    generator = RandomExponential,
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            FitParameter("length", 10.0,
                    displayName = "Length L of the Cylinder",
                    generator = RandomExponential,
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            FitParameter("psiAngle", 0.0,
                    displayName = "Internal Parameter -- ignore",
                    generator = RandomUniform,
                    valueRange = (0.01, 180.01), suffix = "deg."),
            FitParameter("psiAngleDivisions", 303.,
                    displayName = "Orientation Integration Divisions",
                    valueRange = (1, numpy.inf), suffix = "-"),
            FitParameter("lengthIsAspect", True,
                    displayName = "length value indicates aspect ratio"),
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((0.1, 1e3))
        self.length.setValueRange((1, numpy.inf))

    def formfactor(self, dataset):
        # psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        # used in the equation for a cylinder from Pedersen, 1997

        dToR = pi/180. #degrees to radian
        psiRange = self.psiAngle.valueRange()
        psi = numpy.linspace(psiRange[0], psiRange[1], self.psiAngleDivisions())

        if self.lengthIsAspect():
            halfLength = self.radius() * self.length()
        else:
            halfLength = 0.5 * self.length()
        qRsina = numpy.outer(dataset.q, self.radius() * sin((psi) * dToR))
        qLcosa = numpy.outer(dataset.q, halfLength * cos((psi) * dToR))
        fsplit = ((2.*scipy.special.j1(qRsina)/qRsina * sinc(qLcosa/pi))
                  * sqrt(abs(sin((psi) * dToR))[newaxis,:] + 0. * qRsina))
        #integrate over orientation
        return numpy.sqrt(numpy.mean(fsplit**2, axis=1)) # should be length q

    def volume(self):
        if self.lengthIsAspect():
            halfLength = self.radius() * self.length()
        else:
            halfLength = self.length() / 2.
        v = pi * self.radius()**2 * (halfLength * 2.)
        return v**self.compensationExponent

CylindersIsotropic.factory()

#if __name__ == "__main__":
#    from cutesnake.datafile import PDHFile, AsciiFile
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
#    model = CylindersIsotropic()
#    model.radius.setValue(1.)
#    model.radius.setActive(False)
#    model.length.setValue(100.)
#    model.length.setActive(False)
#    model.psiAngle.setValue(1.)
#    model.psiAngle.setActive(False)
#    model.psiAngleDivisions.setValue(303)
#    model.psiAngleDivisions.setActive(False)
#    intensity = (model.formfactor(pf.data, None).reshape(-1))**2
#    q = pf.data[:, 0]
#    oldInt = pf.data[:, 1]
#    delta = abs(oldInt - intensity)
#    result = numpy.dstack((q, intensity, delta))[0]
#    AsciiFile.writeFile("CylindersIsotropic.dat", result)
#    # call it like this:
#    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
