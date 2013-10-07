# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis
from cutesnake.algorithm import Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class CylindersIsotropicLength(ScatteringModel):
    r"""Form factor of cylinders
    which are radially isotropic (so not spherically isotropic!)
    !!!completed but not verified!!!
    """
    shortName = "Cylinders defined by length"
    parameters = (
            Parameter("radius", 1.0,
                    displayName = "Cylinder radius",
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            Parameter("length", 10.0,
                    displayName = "length L of the cylinder",
                    valueRange = (0.1, numpy.inf), suffix = "nm"),
            Parameter("psiAngle", 10.0,
                    displayName = "in-plane cylinder rotation",
                    valueRange = (0.1, 360.1), suffix = "deg."),
            Parameter("psiAngleDivisions", 303.,
                    displayName = "in-plane angle divisions",
                    valueRange = (1, numpy.inf), suffix = "-"),
    )
    parameters[0].isActive = True
    parameters[1].isActive = False#not expected to vary
    parameters[2].isActive = True #better when random 
    parameters[3].isActive = False #not expected to vary

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
        psiAngle = numpy.array((self.psiAngle.value(),))
        psiAngleDivisions = numpy.array((self.psiAngleDivisions.value(),))

        if self.radius.isActive:
            radius = paramValues[:, 0]
        if self.length.isActive:
            idx = int(self.radius.isActive)
            length = paramValues[:, idx]
        if self.psiAngle.isActive:
            #Question: can we nog simply do idx+=1 here?
            idx = numpy.sum((
                    int(self.radius.isActive),
                    int(self.length.isActive)))
            psiAngle = paramValues[:, idx]
        if self.psiAngleDivisions.isActive:
            #not expected to be variable.
            idx = numpy.sum((
                    int(self.radius.isActive),
                    int(self.length.isActive),
                    int(self.psiAngle.isActive)))
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
            lengtH=length[ri%len(length)]
            radi=radius[ri%len(radius)]
            qRsina=numpy.outer(q,radi*sin((psi)*dToR))
            qLcosa=numpy.outer(q,lengtH/2.*cos((psi)*dToR))
            #this is wrong:
            #qRsina=numpy.outer(q,radi*sin((psi-psiA)*dToR))
            #qLcosa=numpy.outer(q,lengtH/2.*cos((psi-psiA)*dToR))
            #fsplit=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )*sqrt(abs(sin((psi-psiA)*dToR))[newaxis,:]+0*qRsina)
            fsplit=( 2*scipy.special.j1(qRsina)/qRsina * sin(qLcosa)/qLcosa )*sqrt(abs(sin((psi)*dToR))[newaxis,:]+0*qRsina)
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

        v = pi*radius**2*(length)                                     
        return v**compensationExponent          

CylindersIsotropicLength.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
    model = CylindersIsotropicLength()
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
    AsciiFile.writeFile("CylindersIsotropicLength.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
