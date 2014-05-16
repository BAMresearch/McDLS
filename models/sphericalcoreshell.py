# -*- coding: utf-8 -*-
# models/SphericalCoreShell.py

import numpy, scipy, scipy.special
from numpy import pi, zeros, sin, cos, sqrt, newaxis, sinc
from utils.parameter import FitParameter, Parameter
from scatteringmodel import ScatteringModel
from cutesnake.algorithm import RandomExponential, RandomUniform

class SphericalCoreShell(ScatteringModel):
    r"""Form factor for a spherical core shell structure
    as defined in the SASfit manual (par. 3.1.4, Spherical Shell III).
    One modification is the ability to specify SLD for core, shell and
    solvent, identical to the notation used in the Core-shell ellipsoid.
    Compared wiht a SASfit-generated model (both with and without distribution)
    """
    shortName = "Core-Shell Sphere"
    parameters = (
            FitParameter("radius", 1.0,
                    displayName = "Core Radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            FitParameter("t", 1.0,
                    displayName = "Thickness of Shell",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            FitParameter("eta_c", 3.16,
                    displayName = "Core SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
            FitParameter("eta_s", 2.53,
                    displayName = "Shell SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
            FitParameter("eta_sol", 0.,
                    displayName = "Solvent SLD",
                    generator = RandomUniform,
                    valueRange = (0, numpy.inf), suffix = "-"),
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.radius.setValueRange((0.1, 1e3))
        self.t.setValueRange((0.1, 1e3))

    def formfactor(self, dataset):
        def k(q, r, dEta):
            # modified K, taken out the volume scaling
            qr = numpy.outer(q, r)
            k = dEta * 3 * (
                    sin(qr) - qr * cos(qr)
                    ) / (qr)**3
            return k

        # dToR = pi / 180. #degrees to radian

        vc = 4./3 * pi *  self.radius() **3
        vt = 4./3 * pi * (self.radius() + self.t()) ** 3
        vRatio = vc / vt

        ks = k(dataset.q, self.radius() + self.t(),
                          self.eta_s() - self.eta_sol())
        kc = k(dataset.q, self.radius(),
                          self.eta_s() - self.eta_c())
        return ks - vRatio * kc
        #integrate over orientation
        #Fell[:,ri]=numpy.sqrt(numpy.mean(fsplit**2, axis=1)) #should be length q

    def volume(self):
        v = 4./3 * pi * (self.radius() + self.t())**3
        return v**self.compensationExponent

SphericalCoreShell.factory()

#if __name__ == "__main__":
#    import sys
#    sys.path.append('..')
#    sys.path.append('.')
#    sys.path.append('../utils')
#    sys.path.append('../cutesnake')
#    from cutesnake.datafile import PDHFile, AsciiFile
#    from models.SphericalCoreShell import SphericalCoreShell
#    # FIXME: use SASData.load() instead
#    pf = PDHFile("testData/SphCoreShell_R100_dR150_c3p16_s2p53.csv")
#    model = SphericalCoreShell()
#    model.radius.setValue(100.)
#    model.radius.setActive(False)
#    model.t.setValue(150.)
#    model.t.setActive(False)
#    model.eta_c.setValue(3.16)
#    model.eta_c.setActive(False)
#    model.eta_s.setValue(2.53)
#    model.eta_s.setActive(False)
#    model.eta_sol.setValue(0.)
#    model.eta_sol.setActive(False)
#    intensity = (model.formfactor(pf.data, None).reshape(-1))**2
#    q = pf.data[:, 0]
#    oldInt = pf.data[:, 1]
#    #normalize
#    intensity /= intensity.max()
#    oldInt /= oldInt.max()
#    delta = abs(oldInt - intensity)
#    print('mean Delta: {}'.format(delta.mean()))
#    result = numpy.dstack((q, intensity, delta))[0]
#    AsciiFile.writeFile("SphericalCoreShell.dat", result)
#    # call it like this:
#    # PYTHONPATH=..:../mcsas/ python SphericalCoreShell.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
