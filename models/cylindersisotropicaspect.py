# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis
from utils.parameter import FitParameter, Parameter
from models.scatteringmodel import SASModel
from utils.units import Length, Angle

# parameters must not be inf

class CylindersIsotropic(SASModel):
    r"""Form factor of cylinders
    which are radially isotropic (so not spherically isotropic!)
    !!!completed but not verified!!!
    """
    shortName = "Cylinders defined by aspect ratio"
    parameters = (
            FitParameter("radius", Length("nm").toSi(1.), unit = Length("nm"),
                    displayName = "Cylinder radius",
                    valueRange = (0., numpy.inf),
                    activeRange = Length("nm").toSi((0.1, 1e3)),
                    suffix = "nm"),
            FitParameter("aspect", 10.0,
                    displayName = "Aspect ratio L/(2R) of the cylinder",
                    valueRange = (0., numpy.inf),
                    activeRange = (1.0, 20),
                    suffix = "-"),
            FitParameter("psiAngle", Angle(u"°").toSi(10.), unit = Angle(u"°"),
                    displayName = "in-plane cylinder rotation",
                    valueRange = (0., Angle(u"°").toSi(180.))),
            FitParameter("psiAngleDivisions", 303.,
                    displayName = "in-plane angle divisions",
                    valueRange = (0, numpy.inf), suffix = "-"),
    )

    def __init__(self):
        super(CylindersIsotropic, self).__init__()
        # some presets of parameters to fit
        self.radius.setActive(True)
        self.aspect.setActive(False) # not expected to vary
        self.psiAngle.setActive(True)  # better when random
        self.psiAngleDivisions.setActive(False) # not expected to vary

    def formfactor(self, dataset):
        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #used in the equation for a cylinder from Pedersen, 1997

        dToR = pi/180. #degrees to radian
        psiRange = self.psiAngle.valueRange()
        psi = numpy.linspace(psiRange[0], psiRange[1], self.psiAngleDivisions())

        ##replicate so we cover all possible combinations of psi, phi and psi
        #psiLong=psi[ numpy.sort( numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(psi) ) ] #indexed to 111222333444 etc
        #qLong=q[ numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(q) ] #indexed to 1234123412341234 etc

        #this is wrong for spherical symmetry:
        #qRsina=numpy.outer(q,radi*sin(((psi-psiA)*dToR)%180))
        #qLcosa=numpy.outer(q,radi*asp*cos(((psi-psiA)*dToR)%180))
        #fsplit=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )*sqrt((sin((psi-psiA)*dToR))[newaxis,:]%180+0*qRsina)
        qRsina = numpy.outer(dataset.q, self.radius() * sin((psi * dToR) % 180.))
        qLcosa = numpy.outer(dataset.q, self.radius() * self.aspect() * cos((psi * dToR) % 180.))
        fsplit = ((2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa)
                  * sqrt((sin(psi * dToR))[newaxis,:]%180 + 0. * qRsina))
        #integrate over orientation
        return numpy.sqrt(numpy.mean(fsplit**2, axis=1)) #should be length q

    def volume(self):
        v = pi * self.radius()**2 * (2. * self.radius() * self.aspect())
        return v

CylindersIsotropic.factory()

#if __name__ == "__main__":
#    from bases.datafile import PDHFile, AsciiFile
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
#    model = CylindersIsotropic()
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
#    AsciiFile.writeFile("CylindersIsotropic.dat", result)
#    # call it like this:
#    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
