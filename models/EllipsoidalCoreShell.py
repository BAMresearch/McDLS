# -*- coding: utf-8 -*-
# models/cylinders.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from cutesnake.algorithm import Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class EllipsoidalCoreShell(ScatteringModel):
    r"""Form factor for an ellipsoidal core shell structure
    as defined in the sasFIT manual (par. 3.2.3)
    """
    shortName = "Core-shell ellipsoid"
    parameters = (
            Parameter("a", 1.0,
                    displayName = "Principal core radius",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("b", 10.0,
                    displayName = "Equatorial core radius",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("t", 1.0,
                    displayName = "Thickness of shell",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("eta_c", 303.,
                    displayName = "core SLD",
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("eta_s", 303.,
                    displayName = "shell SLD",
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("eta_sol", 303.,
                    displayName = "solvent SLD",
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("intDiv", 100,
                    displayName = "orientation integration divisions"),
                    valueRange = (0, 1e4), suffix = "-"),
    )
    parameters[0].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.a.setValueRange((0.1, 1e3))
        self.b.setValueRange((1., 1e4))
        self.t.setValueRange((0.1, 1e3))

    def ff(self, dataset, paramValues):
        def j1(x):
            return ( sin(x) - x * cos(x) ) / (x**2)

        def xc(Q, a, b, mu):
            sfact = sqrt( 
                    a**2 * mu**2 + b**2 * (1 - mu**2)
                    )
            return numpy.outer(Q, sfact)

        def xt(Q, a, b, t, mu):
            sfact = sqrt( 
                    (a + t)**2 * mu**2 + (b + t)**2 * (1 - mu**2)
                    )
            return numpy.outer(Q, sfact)

        assert ScatteringModel.ff(self, dataset, paramValues)

        # vectorized data and arguments
        q = dataset[:, 0]
        a = numpy.array((self.a.value(),))
        b = numpy.array((self.b.value(),))
        t = numpy.array((self.t.value(),))
        eta_c = numpy.array((self.eta_c.value(),))
        eta_s = numpy.array((self.eta_s.value(),))
        eta_sol = numpy.array((self.eta_sol.value(),))
        intDiv = numpy.array((self.intDiv.value(),))
        #unused:

        if self.a.isActive:
            a = paramValues[:, 0]
        if self.b.isActive:
            idx = int(self.a.isActive)
            b = paramValues[:, idx]
        if self.t.isActive:
            idx = int(self.a.isActive) + int(self.b.isActive)
            t = paramValues[:, idx]
        #remaining parameters are never active fitting parameters    

        #psi and phi defined in fig. 1, Pauw et al, J. Appl. Cryst. 2010
        #used in the equation for a cylinder from Pedersen, 1997

        dToR = pi / 180. #degrees to radian
        intVal = numpy.linspace(0., 1., intDiv)
        
        Vc = 4./3 * pi * a * b **2
        Vt = 4./3 * pi * (a + t) * (b + t) ** 2
        VRatio = Vc / Vt 
        if not isinstance(VRatio, np.ndarray):
            VRatio = np.array(VRatio)
        
        Fell=zeros((len(q), len(radius)))
        for ri in range(len(a)):
            arad = a[ ri % len(a) ]
            brad = b[ ri % len(b) ]
            ti = t[ ri % len(t) ]
            qRsina=numpy.outer(q,radi*sin((psi)*dToR))
            qLcosa=numpy.outer(q,halfLength*cos((psi)*dToR))
            Xc = xc(q, a, b, intVal)
            Xt = xt(q, a, b, t, intVal)
            fsplit = ( 
                    (eta_c - eta_s) * VRatio[ri % len(VRatio)] * 
                    ( 3 * j1( Xc ) / Xc ) + 
                    (eta_s - eta_sol) * 1. * 
                    ( 3 * j1( Xt ) / Xt )
                    )
            #integrate over orientation
            Fcyl[:,ri]=numpy.sqrt(numpy.mean(fsplit**2, axis=1)) #should be length q

        return Fcyl

    def vol(self, paramValues, compensationExponent = None):                   
        assert ScatteringModel.vol(self, paramValues)                          
        if compensationExponent is None:                                       
            compensationExponent = self.compensationExponent                   

        a = numpy.array((self.a.value(),))
        b = numpy.array((self.b.value(),))
        t = numpy.array((self.t.value(),))

        if self.a.isActive:
            a = paramValues[:, 0]
        if self.b.isActive:
            idx = int(self.a.isActive)
            b = paramValues[:, idx]
        if self.t.isActive:
            idx = int(self.a.isActive) + int(self.b.isActive)
            t = paramValues[:, idx]

        v = 4./3 * pi * (a + t) * (b + t)**2                                     
        return v**compensationExponent          

EllipsoidalCoreShell.factory()

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
    intensity = (model.ff(pf.data, None).reshape(-1))**2
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("CylindersIsotropic.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
