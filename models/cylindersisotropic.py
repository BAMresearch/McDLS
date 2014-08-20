# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from cutesnake.algorithm import RandomExponential, RandomUniform
from sasunit import SASUnit

# parameters must not be inf

class CylindersIsotropic(ScatteringModel):
    r"""Form factor of cylinders
    previous version (length-fixed) checked against SASfit
    """
    shortName = "Isotropic Cylinders"
    parameters = (
            FitParameter("radius", 1.0e-9,
                    displayName = "Cylinder Radius",
                    generator = RandomExponential,
                    valueRange = (0.1e-10, numpy.inf)),
            Parameter("useAspect", True,
                    displayName = "Use aspect ratio (checked) or length "),
            FitParameter("length", 10.0e-9,
                    displayName = "Length L of the Cylinder",
                    generator = RandomExponential,
                    valueRange = (1e-10, 1e1 )),
            FitParameter("aspect", 10.0,
                    displayName = "Aspect ratio of the Cylinder",
                    generator = RandomExponential,
                    valueRange = (1e-3, 1e3)),
            FitParameter("psiAngle", 0.1,
                    displayName = "Internal Parameter, not user adjustable",
                    generator = RandomUniform,
                    valueRange = (0.01, 2 * pi + 0.01)),
            Parameter("psiAngleDivisions", 303.,
                    displayName = "Orientation Integration Divisions",
                    valueRange = (1, numpy.inf)),
            Parameter("sld", 1e14,
                    displayName = "scattering length density",
                    valueRange = (0, numpy.inf))
    )
    parameters[0].setActive(True)
    parameters[0].unit = SASUnit(magnitudedict = "length",
            simagnitudename = u"m",
            displaymagnitudename = u'nm')
    parameters[2].unit = SASUnit(magnitudedict = "length",
            simagnitudename = u"m",
            displaymagnitudename = u'nm')
    parameters[3].unit = SASUnit(magnitudedict = "none",
            simagnitudename = u"-",
            displaymagnitudename = u'-')
    parameters[4].unit = SASUnit(magnitudedict = "angle",
            simagnitudename = u"rad",
            displaymagnitudename = u'˚')
    parameters[5].unit = SASUnit(magnitudedict = "none",
            simagnitudename = u"-",
            displaymagnitudename = u'-')
    parameters[6].unit = SASUnit(magnitudedict = "SLD",
            simagnitudename = u'm⁻²',
            displaymagnitudename = u'Å⁻²')

    def __init__(self):
        ScatteringModel.__init__(self)

    def formfactor(self, dataset):
        # psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        # used in the equation for a cylinder from Pedersen, 1997

        psiRange = self.psiAngle.valueRange()
        psi = numpy.linspace(psiRange[0], psiRange[1], self.psiAngleDivisions())

        if self.useAspect():
            halfLength = self.radius() * self.aspect()
        else:
            halfLength = 0.5 * self.length()
        qRsina = numpy.outer(dataset.q, self.radius() * sin((psi) ))
        qLcosa = numpy.outer(dataset.q, halfLength * cos((psi) ))
        fsplit = ((2.*scipy.special.j1(qRsina)/qRsina * sinc(qLcosa/pi))
                  * sqrt(abs(sin((psi) ))[newaxis,:] + 0. * qRsina))
        #integrate over orientation
        return numpy.sqrt(numpy.mean(fsplit**2, axis=1)) # should be length q

    def volume(self):
        if self.useAspect():
            halfLength = self.radius() * self.aspect()
        else:
            halfLength = self.length() / 2.
        v = pi * self.radius()**2 * (halfLength * 2.)
        return v**self.compensationExponent

    def absVolume(self):
        return self.volume() * self.sld()**2

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
