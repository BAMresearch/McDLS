# -*- coding: utf-8 -*-
# models/EllipsoidalCoreShell.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from bases.algorithm import RandomUniform, RandomExponential
from scatteringmodel import ScatteringModel
from utils.units import Length, NoUnit, SLD

# parameters must not be inf

class EllipsoidalCoreShell(ScatteringModel):
    r"""Form factor for an ellipsoidal core shell structure
    as defined in the SASfit manual (par. 3.2.3)
    Tested 2014-01-21 against SASfit function with good agreement.
    """
    shortName = "Core-Shell Ellipsoid"
    parameters = (
            FitParameter("a", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Principal Core Radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf)),
            FitParameter("b", Length(u'nm').toSi(10.), unit = Length(u'nm'),
                    displayName = "Equatorial Core Radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf)),
            FitParameter("t", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Thickness of Shell",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf)),
            Parameter("eta_c", SLD(u'Å⁻²').toSi(3.15e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Core SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
            Parameter("eta_s", SLD(u'Å⁻²').toSi(2.53e-6), unit = SLD(u'Å⁻²'),
                    displayName = "Shell SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
            Parameter("eta_sol", 0., unit = SLD(u'Å⁻²'),
                    displayName = "Solvent SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf)),
            Parameter("intDiv", 100,
                    displayName = "Orientation Integration Divisions",
                    generator = RandomUniform,
                    valueRange = (0, 1e4)),
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.a.setDisplayActiveRange((0.1, 1e3))
        self.b.setDisplayActiveRange((1.0, 1e4))
        self.t.setDisplayActiveRange((0.1, 1e3))

    def formfactor(self, dataset):
        def j1(x):
            return ( sin(x) - x * cos(x) ) / (x**2)

        def calcXc(q, a, b, mu):
            sfact = sqrt(
                    a**2 * mu**2 + b**2 * (1 - mu**2)
                    )
            return numpy.outer(q, sfact)

        def calcXt(q, a, b, t, mu):
            sfact = sqrt(
                    (a + t)**2 * mu**2 + (b + t)**2 * (1 - mu**2)
                    )
            return numpy.outer(q, sfact)

        intVal = numpy.linspace(0., 1., self.intDiv())

        vc = 4./3. * pi *  self.a() * self.b() **2.
        vt = 4./3. * pi * (self.a() + self.t()) * (self.b() + self.t()) **2.
        vRatio = vc / vt

        xc = calcXc(dataset.q, self.a(), self.b(), intVal)
        xt = calcXt(dataset.q, self.a(), self.b(), self.t(), intVal)
        fsplit = (
                (self.eta_c() - self.eta_s()) * vRatio *
                ( 3 * j1( xc ) / xc ) +
                (self.eta_s() - self.eta_sol()) * 1. *
                ( 3 * j1( xt ) / xt )
                )
        # integrate over orientation
        return numpy.sqrt(numpy.mean(fsplit**2, axis=1)) # should be length q

    def volume(self):
        v = 4./3 * pi * (self.a() + self.t()) * (self.b() + self.t())**2
        return v**self.compensationExponent

    def absVolume(self):
        return self.volume() #TODO: check how to do this.

EllipsoidalCoreShell.factory()

#if __name__ == "__main__":
#    import sys
#    sys.path.append('..')
#    sys.path.append('.')
#    sys.path.append('../utils')
#    from bases.datafile import PDHFile, AsciiFile
#    from models.EllipsoidalCoreShell import EllipsoidalCoreShell
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("testData/EllCoreShell_a100_b150_t500_c3p16_s2p53_sol0.csv")
#    model = EllipsoidalCoreShell()
#    model.a.setValue(100.)
#    model.a.setActive(False)
#    model.b.setValue(150.)
#    model.b.setActive(False)
#    model.t.setValue(500.)
#    model.t.setActive(False)
#    model.eta_c.setValue(3.16)
#    model.eta_c.setActive(False)
#    model.eta_s.setValue(2.53)
#    model.eta_s.setActive(False)
#    model.eta_sol.setValue(0.)
#    model.eta_sol.setActive(False)
#    model.intDiv.setValue(100)
#    model.intDiv.setActive(False)
#    intensity = (model.formfactor(pf.data, None).reshape(-1))**2
#    q = pf.data[:, 0]
#    oldInt = pf.data[:, 1]
#    #normalize
#    intensity /= intensity.max()
#    oldInt /= oldInt.max()
#    delta = abs(oldInt - intensity)
#    print('mean Delta: {}'.format(delta.mean()))
#    result = numpy.dstack((q, intensity, delta))[0]
#    AsciiFile.writeFile("EllipsoidalCoreShell.dat", result)
#    # call it like this:
#    # PYTHONPATH=..:../mcsas/ python EllipsoidalCoreShell.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
