# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis
from utils.parameter import Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class CylindersIsotropic(ScatteringModel):
    r"""Form factor of cylinders
    which are radially isotropic (so not spherically isotropic!)
    !!!completed but not verified!!!
    """
    shortName = "Cylinders defined by aspect ratio"
    parameters = (
            Parameter("radius", 1.0,
                    displayName = "Cylinder radius",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("aspect", 10.0,
                    displayName = "Aspect ratio L/(2R) of the cylinder",
                    valueRange = (0., numpy.inf), suffix = "-"),
            Parameter("psiAngle", 10.0,
                    displayName = "in-plane cylinder rotation",
                    valueRange = (0., 180.), suffix = "deg."),
            Parameter("psiAngleDivisions", 303.,
                    displayName = "in-plane angle divisions",
                    valueRange = (0, numpy.inf), suffix = "-"),
    )
    parameters[0].setActive(True)
    parameters[1].setActive(False) # not expected to vary
    parameters[2].setActive(True)  # better when random
    parameters[3].setActive(False) # not expected to vary

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((0.1, 1e3))
        self.aspect.setValueRange((1, 20))

    def formfactor(self, dataset, paramValues):

        # vectorized data and arguments
        q = dataset.q
        radius = numpy.array((self.radius(),))
        aspect = numpy.array((self.aspect(),))
        psiAngle = numpy.array((self.psiAngle(),))
        psiAngleDivisions = numpy.array((self.psiAngleDivisions(),))

        if self.radius.isActive():
            radius = paramValues[:, 0]
        if self.aspect.isActive():
            idx = int(self.radius.isActive())
            aspect = paramValues[:, idx]
        if self.psiAngle.isActive():
            #Question: can we nog simply do idx+=1 here?
            idx = numpy.sum((
                    int(self.radius.isActive()),
                    int(self.aspect.isActive())))
            psiAngle = paramValues[:, idx]
        if self.psiAngleDivisions.isActive():
            #not expected to be variable.
            idx = numpy.sum((
                    int(self.radius.isActive()),
                    int(self.aspect.isActive()),
                    int(self.psiAngle.isActive())))
            psiAngleDivisions = paramValues[:, idx]

        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #used in the equation for a cylinder from Pedersen, 1997

        dToR=pi/180. #degrees to radian
        psiRange=self.psiAngle.valueRange()
        psi=numpy.linspace(
                psiRange[0],psiRange[1],psiAngleDivisions)

        ##replicate so we cover all possible combinations of psi, phi and psi
        #psiLong=psi[ numpy.sort( numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(psi) ) ] #indexed to 111222333444 etc
        #qLong=q[ numpy.array( range(
        #    (len(psi)*len(q))
        #    ) ) %len(q) ] #indexed to 1234123412341234 etc

        fsplit=zeros((len(q),len(psi)))
        Fcyl=zeros((len(q),len(radius)))
        for ri in range(len(radius)):
            psiA=psiAngle[ri%len(psiAngle)]
            asp=aspect[ri%len(aspect)]
            radi=radius[ri%len(radius)]
            #this is wrong for spherical symmetry:
            #qRsina=numpy.outer(q,radi*sin(((psi-psiA)*dToR)%180))
            #qLcosa=numpy.outer(q,radi*asp*cos(((psi-psiA)*dToR)%180))
            #fsplit=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )*sqrt((sin((psi-psiA)*dToR))[newaxis,:]%180+0*qRsina)
            qRsina=numpy.outer(q,radi*sin(((psi)*dToR)%180))
            qLcosa=numpy.outer(q,radi*asp*cos(((psi)*dToR)%180))
            fsplit=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )*sqrt((sin((psi)*dToR))[newaxis,:]%180+0*qRsina)
            #integrate over orientation
            Fcyl[:,ri]=numpy.sqrt(numpy.mean(fsplit**2,axis=1)) #should be length q

        return Fcyl

    def volume(self, paramValues, compensationExponent = None):
        if compensationExponent is None:
            compensationExponent = self.compensationExponent
        idx = 0
        radius=paramValues[:,0]
        if self.aspect.isActive():
            idx+=1
            aspect =paramValues[:,idx]
        else:
            aspect=self.aspect()

        v = pi*radius**2*(2*radius*aspect)
        return v**compensationExponent

CylindersIsotropic.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    # FIXME: use SASData.load() instead
    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
    model = CylindersIsotropic()
    model.radius.setValue(1.)
    model.radius.setActive(False)
    model.aspect.setValue(100.)
    model.aspect.setActive(False)
    model.psiAngle.setValue(1.)
    model.psiAngle.setActive(False)
    model.psiAngleDivisions.setValue(303)
    model.psiAngleDivisions.setActive(False)
    intensity = model.formfactor(pf.data, None).reshape(-1)
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("CylindersIsotropic.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
