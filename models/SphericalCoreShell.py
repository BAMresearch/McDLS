# -*- coding: utf-8 -*-
# models/SphericalCoreShell.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from cutesnake.algorithm import Parameter
from scatteringmodel import ScatteringModel

# parameters must not be inf

class SphericalCoreShell(ScatteringModel):
    r"""Form factor for a spherical core shell structure
    as defined in the SASfit manual (par. 3.1.4, Spherical Shell III).
    One modification is the ability to specify SLD for core, shell and 
    solvent, identical to the notation used in the Core-shell ellipsoid.
    Compared wiht a SASfit-generated model (no distribution)
    """
    shortName = "Core-shell ellipsoid"
    parameters = (
            Parameter("radius", 1.0,
                    displayName = "Core radius",
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
    )
    parameters[0].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((0.1, 1e3))
        self.t.setValueRange((0.1, 1e3))

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)

        def K(Q, R, DEta):
            #modified K, taken out the volume scaling
            QR = numpy.outer(Q, R)
            k = DEta * 3 * (
                    sin(QR) - QR * cos(QR)
                    ) / (QR)**3
            return k


        # vectorized data and arguments
        q = dataset[:, 0]
        radius = numpy.array((self.radius.value(),))
        t = numpy.array((self.t.value(),))
        eta_c = numpy.array((self.eta_c.value(),))
        eta_s = numpy.array((self.eta_s.value(),))
        eta_sol = numpy.array((self.eta_sol.value(),))
        #unused:

        idx = 0
        if self.radius.isActive:
            radius = paramValues[:, idx]
            idx += 1
        if self.t.isActive:
            t = paramValues[:, idx]
            idx += 1
        if self.eta_c.isActive:
            eta_c = paramValues[:, idx]
            idx += 1
        if self.eta_s.isActive:
            eta_s = paramValues[:, idx]
            idx += 1
        if self.eta_sol.isActive:
            eta_sol = paramValues[:, idx]
            idx += 1
        #remaining parameters are never active fitting parameters    

        dToR = pi / 180. #degrees to radian
        
        Vc = 4./3 * pi * radius **3
        Vt = 4./3 * pi * (radius + t) ** 3
        VRatio = Vc / Vt 
        if not isinstance(VRatio, numpy.ndarray):
            VRatio = numpy.array(VRatio)
        
        if paramValues is None:
            Fell=zeros((len(q), 1 ))
        else:
            Fell=zeros((len(q), len(paramValues[:,(idx-1)])))

        for ri in range(len(radius)):
            rad = radius[ ri % len(radius) ]
            ti = t[ ri % len(t) ]
            VRati = VRatio[ ri % len(VRatio) ] 
            Ks = K(q, (rad + ti), (eta_s - eta_sol))
            Kc = K(q, rad, (eta_s - eta_c))
            print('shape Ks: {}, shape Kc: {}, shape VRatio: {}'.format(
                numpy.shape(Ks), numpy.shape(Kc), numpy.shape(VRatio))) 
            Fell[:,ri] = ( Ks - VRati * Kc ).flatten()
            #integrate over orientation
            #Fell[:,ri]=numpy.sqrt(numpy.mean(fsplit**2, axis=1)) #should be length q

        return Fell

    def vol(self, paramValues, compensationExponent = None):                   
        assert ScatteringModel.vol(self, paramValues)                          
        if compensationExponent is None:                                       
            compensationExponent = self.compensationExponent                   

        radius = numpy.array((self.radius.value(),))
        t = numpy.array((self.t.value(),))

        if self.radius.isActive:
            radius = paramValues[:, 0]
        if self.t.isActive:
            idx = int(self.radius.isActive)
            t = paramValues[:, idx]

        v = 4./3 * pi * (radius + t)**3                                     
        return v**compensationExponent          

SphericalCoreShell.factory()

if __name__ == "__main__":
    import sys
    sys.path.append('..')
    sys.path.append('.')
    sys.path.append('../utils')
    sys.path.append('../cutesnake')
    from cutesnake.datafile import PDHFile, AsciiFile
    from models.SphericalCoreShell import SphericalCoreShell
    pf = PDHFile("testData/SphCoreShell_R100_dR150_c3p16_s2p53.csv")
    model = SphericalCoreShell()
    model.radius.setValue(100.)
    model.radius.isActive = False
    model.t.setValue(150.)
    model.t.isActive = False
    model.eta_c.setValue(3.16)
    model.eta_c.isActive = False
    model.eta_s.setValue(2.53)
    model.eta_s.isActive = False
    model.eta_sol.setValue(0.)
    model.eta_sol.isActive = False
    intensity = (model.ff(pf.data, None).reshape(-1))**2
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    #normalize
    intensity /= intensity.max()
    oldInt /= oldInt.max()
    delta = abs(oldInt - intensity)
    print('mean Delta: {}'.format(delta.mean()))
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("SphericalCoreShell.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python SphericalCoreShell.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
