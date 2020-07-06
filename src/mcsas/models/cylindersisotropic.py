# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
import numpy as np
import logging
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from bases.model import SASModel
from bases.algorithm import RandomExponential, RandomUniform
from utils.units import Length, NoUnit, Angle, SLD

# parameters must not be inf

class CylindersIsotropic(SASModel):
    r"""Form factor of cylinders
    previous version (length-fixed) checked against SASfit
    """
    shortName = "SASfit Isotropic Cylinders"
    parameters = (
            FitParameter("radius", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Cylinder Radius",
                    generator = RandomExponential,
                    valueRange = (Length(u'nm').toSi(0.1), numpy.inf)),
            Parameter("useAspect", True,
                    displayName = "Use aspect ratio (checked) or length "),
            FitParameter("length", Length(u'nm').toSi(10.), unit = Length(u'nm'),
                    displayName = "Length L of the Cylinder",
                    generator = RandomExponential,
                    valueRange = (Length(u'nm').toSi(0.1),
                                  Length(u'nm').toSi(1e10))),
            FitParameter("aspect", 10.0,
                    displayName = "Aspect ratio of the Cylinder",
                    generator = RandomExponential,
                    valueRange = (1e-3, 1e3)),
            Parameter("intDiv", 100.,
                    displayName = "Orientation Integration Divisions",
                    valueRange = (1, 1e4)),
            Parameter("sld", SLD(u'Å⁻²').toSi(1e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Scattering length density difference",
                    valueRange = (0, numpy.inf))
    )

    def __init__(self):
        super(CylindersIsotropic, self).__init__()
        # some presets of parameters to fit
        self.radius.setActive(True)

    def formfactor(self, dataset):
        # psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        # used in the equation for a cylinder from Pedersen, 1997
        # reworked on 20171010, changed to SASfit function (eq.3.215, sasfit doc 0.94.6)

        # not sure we need the step size, since the integration uses np.average
        x, step = numpy.linspace(0., 1., self.intDiv(), endpoint = True, retstep = True)

        # replace x=0 and x=1 with more easy-to-calculate values, will be replaced below
        # avoiding infinities and nan's may save some time:
        x[0] = 0.5
        x[-1] = 0.5
        # x = x[1:] # clip zero value
        # logging.info("Shape cylinder int. x: {}, q: {}".format(x.shape, dataset.q.shape))

        if self.useAspect():
            halfLength = self.radius() * self.aspect()
        else:
            halfLength = 0.5 * self.length()

        QRsqrtx = numpy.outer(dataset.q, self.radius() * np.sqrt(1. - x**2.) )
        QLx     = numpy.outer(dataset.q, 2. * halfLength * x )
        # shape of the above two: [q, numInt]

        numerator   = scipy.special.j1(QRsqrtx) * np.sin(QLx / 2.)
        denominator = QRsqrtx * QLx
        fsplit = numerator / denominator

        # limit for the function where x -> 0: np.sin(QLx / 2.) / QLx = 0.5
        fsplit[:, 0] = 0.5 * (scipy.special.j1(dataset.q * self.radius()) / 
                (dataset.q * self.radius()))
        # not quite sure, but this might be the limit for x -> 1:
        fsplit[:, -1] = np.sin(dataset.q * halfLength) / (dataset.q * halfLength)

        # shape of fsplit: [q, numInt]

        # fsplit = ((scipy.special.j1(qRsina)/qRsina * sinc(qLcosa/pi))
        #           * sqrt(abs(sin((psi) ))[newaxis,:] + 0. * qRsina))

        # return 16 * np.sqrt(step * np.sum(fsplit**2, axis = 1))
        return np.sqrt(16 * np.trapz(fsplit**2, dx = step, axis = 1))

    def volume(self):
        if self.useAspect():
            halfLength = self.radius() * self.aspect()
        else:
            halfLength = self.length() / 2.
        v = pi * self.radius()**2 * (halfLength * 2.)
        return v

    def absVolume(self):
        return self.volume() * self.sld()**2

CylindersIsotropic.factory()

#if __name__ == "__main__":
#    from bases.datafile import PDHFile, AsciiFile
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
