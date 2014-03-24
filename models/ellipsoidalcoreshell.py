# -*- coding: utf-8 -*-
# models/EllipsoidalCoreShell.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import Parameter
from cutesnake.algorithm import RandomUniform, RandomExponential
from scatteringmodel import ScatteringModel

# parameters must not be inf

class EllipsoidalCoreShell(ScatteringModel):
    r"""Form factor for an ellipsoidal core shell structure
    as defined in the SASfit manual (par. 3.2.3)
    Tested 2014-01-21 against SASfit function with good agreement.
    """
    shortName = "Core-shell ellipsoid"
    parameters = (
            Parameter("a", 1.0,
                    displayName = "Principal core radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("b", 10.0,
                    displayName = "Equatorial core radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("t", 1.0,
                    displayName = "Thickness of shell",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("eta_c", 3.15,
                    displayName = "core SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("eta_s", 2.53,
                    displayName = "shell SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("eta_sol", 0.,
                    displayName = "solvent SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
            Parameter("intDiv", 100,
                    displayName = "orientation integration divisions",
                    valueRange = (0, 1e4), suffix = "-"),
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.a.setValueRange((0.1, 1e3))
        self.b.setValueRange((1., 1e4))
        self.t.setValueRange((0.1, 1e3))

    def formfactor(self, dataset, paramValues):
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


        # vectorized data and arguments
        q = dataset.q
        a = numpy.array((self.a(),))
        b = numpy.array((self.b(),))
        t = numpy.array((self.t(),))
        eta_c = numpy.array((self.eta_c(),))
        eta_s = numpy.array((self.eta_s(),))
        eta_sol = numpy.array((self.eta_sol(),))
        intDiv = numpy.array((self.intDiv(),))
        #unused:

        idx = 0
        if self.a.isActive():
            a = paramValues[:, idx]
            idx += 1
        if self.b.isActive():
            b = paramValues[:, idx]
            idx += 1
        if self.t.isActive():
            t = paramValues[:, idx]
            idx += 1
        if self.eta_c.isActive():
            eta_c = paramValues[:, idx]
            idx += 1
        if self.eta_s.isActive():
            eta_s = paramValues[:, idx]
            idx += 1
        if self.eta_sol.isActive():
            eta_sol = paramValues[:, idx]
            idx += 1
        #remaining parameters are never active fitting parameters

        dToR = pi / 180. #degrees to radian
        intVal = numpy.linspace(0., 1., intDiv)

        Vc = 4./3 * pi * a * b **2
        Vt = 4./3 * pi * (a + t) * (b + t) ** 2
        VRatio = Vc / Vt
        if not isinstance(VRatio, numpy.ndarray):
            VRatio = numpy.array(VRatio)

        if paramValues is None:
            Fell=zeros((len(q), 1 ))
        else:
            Fell=zeros((len(q), len(paramValues[:,(idx-1)])))

        for ri in range(len(a)):
            arad = a[ ri % len(a) ]
            brad = b[ ri % len(b) ]
            ti = t[ ri % len(t) ]
            etac = eta_c[ ri % len(eta_c) ]
            etas = eta_s[ ri % len(eta_s) ]
            etasol = eta_sol[ ri % len(eta_sol) ]
            Xc = xc(q, arad, brad, intVal)
            Xt = xt(q, arad, brad, ti, intVal)
            fsplit = (
                    (etac - etas) * VRatio[ri % len(VRatio)] *
                    ( 3 * j1( Xc ) / Xc ) +
                    (etas - etasol) * 1. *
                    ( 3 * j1( Xt ) / Xt )
                    )
            #integrate over orientation
            Fell[:,ri]=numpy.sqrt(numpy.mean(fsplit**2, axis=1)) #should be length q

        return Fell

    def volume(self, paramValues):
        a = numpy.array((self.a(),))
        b = numpy.array((self.b(),))
        t = numpy.array((self.t(),))

        if self.a.isActive():
            a = paramValues[:, 0]
        if self.b.isActive():
            idx = int(self.a.isActive())
            b = paramValues[:, idx]
        if self.t.isActive():
            idx = int(self.a.isActive()) + int(self.b.isActive())
            t = paramValues[:, idx]

        v = 4./3 * pi * (a + t) * (b + t)**2
        return v**self.compensationExponent

EllipsoidalCoreShell.factory()

if __name__ == "__main__":
    import sys
    sys.path.append('..')
    sys.path.append('.')
    sys.path.append('../utils')
    sys.path.append('../cutesnake')
    from cutesnake.datafile import PDHFile, AsciiFile
    from models.EllipsoidalCoreShell import EllipsoidalCoreShell
    # FIXME: use SASData.load() instead
    pf = PDHFile("testData/EllCoreShell_a100_b150_t500_c3p16_s2p53_sol0.csv")
    model = EllipsoidalCoreShell()
    model.a.setValue(100.)
    model.a.setActive(False)
    model.b.setValue(150.)
    model.b.setActive(False)
    model.t.setValue(500.)
    model.t.setActive(False)
    model.eta_c.setValue(3.16)
    model.eta_c.setActive(False)
    model.eta_s.setValue(2.53)
    model.eta_s.setActive(False)
    model.eta_sol.setValue(0.)
    model.eta_sol.setActive(False)
    model.intDiv.setValue(100)
    model.intDiv.setActive(False)
    intensity = (model.formfactor(pf.data, None).reshape(-1))**2
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    #normalize
    intensity /= intensity.max()
    oldInt /= oldInt.max()
    delta = abs(oldInt - intensity)
    print('mean Delta: {}'.format(delta.mean()))
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("EllipsoidalCoreShell.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python EllipsoidalCoreShell.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
