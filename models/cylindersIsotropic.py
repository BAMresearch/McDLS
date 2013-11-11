# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from cutesnake.algorithm import Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class CylindersIsotropic(ScatteringModel):
    r"""Form factor of cylinders
    previous version (length-fixed) checked against SASfit
    """
    shortName = "Isotropic Cylinders"
    parameters = (
            Parameter("radius", 1.0,
                    displayName = "Cylinder radius",
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            Parameter("length", 10.0,
                    displayName = "length L of the cylinder",
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            Parameter("psiAngle", 0.0,
                    displayName = "internal parameter -- ignore",
                    valueRange = (0.01, 180.01), suffix = "deg."),
            Parameter("psiAngleDivisions", 303.,
                    displayName = "orientation integration divisions",
                    valueRange = (1, numpy.inf), suffix = "-"),
            Parameter("lengthIsAspect", True,
                    displayName = "length value indicates aspect ratio"),
    )
    parameters[0].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((0.1, 1e3))
        self.length.setValueRange((1, numpy.inf))

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)

        # vectorized data and arguments
        q = dataset[:, 0]
        radius = numpy.array((self.radius.value(),))
        length = numpy.array((self.length.value(),))
        psiAngleDivisions = numpy.array((self.psiAngleDivisions.value(),))
        #unused:
        psiAngle = numpy.array((self.psiAngle.value(),))

        if self.radius.isActive:
            radius = paramValues[:, 0]
        if self.length.isActive:
            idx = int(self.radius.isActive)
            length = paramValues[:, idx]
        #remaining parameters are never active fitting parameters    

        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #used in the equation for a cylinder from Pedersen, 1997

        dToR=pi/180. #degrees to radian
        psiRange=self.psiAngle.valueRange()
        psi=numpy.linspace(
                psiRange[0],psiRange[1],psiAngleDivisions)
        
        fsplit=zeros((len(q),len(psi)))
        Fcyl=zeros((len(q),len(radius)))
        for ri in range(len(radius)):
            radi=radius[ri%len(radius)]
            if self.lengthIsAspect.value():
                halfLength=radi*length[ri%len(length)]
            else:
                halfLength=0.5*length[ri%len(length)]
            qRsina=numpy.outer(q,radi*sin((psi)*dToR))
            qLcosa=numpy.outer(q,halfLength*cos((psi)*dToR))
            fsplit=( 2.*scipy.special.j1(qRsina)/qRsina * sinc(qLcosa/pi) )*sqrt(abs(sin((psi)*dToR))[newaxis,:]+0*qRsina)
            #integrate over orientation
            Fcyl[:,ri]=numpy.sqrt(numpy.mean(fsplit**2,axis=1)) #should be length q

        return Fcyl

    def vol(self, paramValues, compensationExponent = None):                   
        assert ScatteringModel.vol(self, paramValues)                          
        if compensationExponent is None:                                       
            compensationExponent = self.compensationExponent                   
        idx = 0
        radius=paramValues[:,0]                                                
        if self.length.isActive:                                               
            idx+=1                                                             
            length =paramValues[:,idx]                                         
        else:                                                                  
            length=self.length.value()                                         
        if self.lengthIsAspect.value():
            halfLength=radius*length
        else:
            halfLength=length/2.

        v = pi*radius**2*(halfLength*2)                                     
        return v**compensationExponent          

CylindersIsotropic.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
    model = CylindersIsotropic()
    model.radius.setValue(1.)
    model.radius.isActive = False
    model.length.setValue(100.)
    model.length.isActive = False
    model.psiAngle.setValue(1.)
    model.psiAngle.isActive = False
    model.psiAngleDivisions.setValue(303)
    model.psiAngleDivisions.isActive = False
    intensity = model.ff(pf.data, None).reshape(-1)
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("CylindersIsotropic.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
