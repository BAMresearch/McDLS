# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos

from mcsas.utils.parameter import FitParameter, Parameter
from mcsas.bases.model import SASModel
from mcsas.bases.algorithm import RandomUniform, RandomExponential
from mcsas.utils.units import Length, NoUnit, Angle, SLD

# parameters must not be inf

class CylindersRadiallyIsotropic(SASModel):
    r"""Form factor of cylinders
    which are radially isotropic (so not spherically isotropic!)
    !!!completed but not verified!!!
    """
    shortName = "Radially (in-plane) isotropic cylinders"
    parameters = (
            FitParameter("radius", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Cylinder radius",
                    generator = RandomExponential,
                    valueRange = Length(u'nm').toSi((0.1, numpy.inf)),
                    activeRange = Length(u'nm').toSi((0.1, 1e3))), # preset
            FitParameter("aspect", 10.0,
                    displayName = "Aspect ratio L/(2R) of the cylinder",
                    generator = RandomUniform,
                    valueRange = (0.1, numpy.inf),
                    activeRange = (1.0, 20.)), # preset
            FitParameter("psiAngle", 0.17, unit = Angle(u'°'),
                    displayName = "in-plane cylinder rotation",
                    generator = RandomUniform,
                    valueRange = (0.01, 2 * pi + 0.01)),
            Parameter("psiAngleDivisions", 303., #setting to int gives OverflowError
                    displayName = "in-plane angle divisions",
                    valueRange = (1, numpy.inf)),
            Parameter("sld", SLD(u'Å⁻²').toSi(1e-6), unit = SLD(u'Å⁻²'),
                    displayName = "scattering length density difference",
                    valueRange = (0., numpy.inf))
    )

    def __init__(self):
        super(CylindersRadiallyIsotropic, self).__init__()
        # some presets of parameters to fit
        self.radius.setActive(True)
        self.aspect.setActive(False)  # not expected to vary
        self.psiAngle.setActive(True) # better when random

    def formfactor(self, dataset):
        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #used in the equation for a cylinder from Pedersen, 1997

        # dToR = pi/180. #degrees to radian not necessary since unit conversion
        psiRange = self.psiAngle.valueRange()
        psi = numpy.linspace(psiRange[0], psiRange[1], self.psiAngleDivisions())

        ##replicate so we cover all possible combinations of psi, phi and psi
        #psiLong=psi[ numpy.sort( numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(psi) ) ] #indexed to 111222333444 etc
        #qLong=q[ numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(q) ] #indexed to 1234123412341234 etc

        #rotation can be used to get slightly better results, but
        #ONLY FOR RADIAL SYMMETRY, NOT SPHERICAL.
        qRsina = numpy.outer(dataset.q, self.radius() * sin(((psi - self.psiAngle()) )))
        qLcosa = numpy.outer(dataset.q, self.radius() * self.aspect() * cos(((psi - self.psiAngle()))))
        #leave the rotation out of it for now.
        #qRsina=numpy.outer(q,radi*sin(((psi)*dToR)))
        #qLcosa=numpy.outer(q,radi*asp*cos(((psi)*dToR)))
        fsplit = (2. * scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa)
        #integrate over orientation
        return numpy.sqrt(numpy.mean(fsplit**2, axis=1)) # should be length q

    def volume(self):
        v = pi * self.radius()**2 * (2. * self.radius() * self.aspect())
        return v

    def absVolume(self):
        return self.volume() * self.sld()**2

CylindersRadiallyIsotropic.factory()

#if __name__ == "__main__":
#    from bases.datafile import PDHFile, AsciiFile
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
#    model = CylindersRadiallyIsotropic()
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
