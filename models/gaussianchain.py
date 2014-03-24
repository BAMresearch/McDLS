# -*- coding: utf-8 -*-
# models/gaussianchain.py

import numpy
from utils.parameter import Parameter
from scatteringmodel import ScatteringModel
from cutesnake.algorithm import RandomUniform, RandomExponential

# parameters must not be inf

class GaussianChain(ScatteringModel):
    r"""Form factor of flexible polymer chains which are not selfavoiding
    and obey Gaussian statistics after [Debye47]_

    See also: http://sasfit.sf.net/manual/Gaussian_Chain#Gauss_2

    .. [Debye47] `P. Debye, Mollecular-weight determination by light
        scattering, Journal of Physical and Colloid Chemistry, 51:18--32,
        1947. <http://dx.doi.org/10.1021/j150451a002>`_

    I_0 = (bp - (k * Rg^2) * eta_s)^2 with k = 1 nm.
    k * Rg^2 = volume approximation
    """
    shortName = "Gaussian Chain"
    parameters = (
            Parameter("rg", 1.0,
                    displayName = "radius of gyration, Rg",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("bp", 100.0,
                    displayName = "scattering length of the polymer",
                    generator = RandomUniform,
                    valueRange = (0., numpy.inf), suffix = "cm"),
            Parameter("etas", 1.0,
                    displayName = "scattering length density of the solvent",
                    generator = RandomUniform,
                    valueRange = (0., numpy.inf), suffix = "cm<sup>-1</sup>"),
            Parameter("k", 1.0,
                    displayName = "volumetric scaling factor of Rg",
                    generator = RandomUniform,
                    valueRange = (0., numpy.inf), suffix = "nm")
    )
    parameters[0].setActive(True)

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.rg.setValueRange((1, 1e2))
        self.bp.setValueRange((0.1, 1e3))
        self.etas.setValueRange((0.1, 10.))
        self.k.setValueRange((0.1, 10.))

    def formfactor(self, dataset, paramValues = None):
        # vectorized data and arguments
        q = dataset.q
        rg = numpy.array((self.rg(),))
        if self.rg.isActive():
            rg = paramValues[:, 0]
        bp = numpy.array((self.bp(),))
        if self.bp.isActive():
            idx = int(self.rg.isActive())
            bp = paramValues[:, idx]
        etas = numpy.array((self.etas(),))
        if self.etas.isActive():
            idx = int(self.rg.isActive()) + int(self.bp.isActive())
            etas = paramValues[:, idx]
        k = numpy.array((self.k(),))
        if self.k.isActive():
            idx = int(self.rg.isActive()) + int(self.bp.isActive()) + int(self.etas.isActive())
            k = paramValues[:, idx]

        beta = bp - (k * rg**2) * etas
        u = numpy.outer(q, rg)**2 # a matrix usually
        result = numpy.sqrt(2.) * numpy.sqrt(numpy.expm1(-u) + u) / u
        for i, res in zip(range(0, result.shape[0]), result):
            if q[i] <= 0.0:
                result[i] = beta
            result[i] = res * beta
        return result

    def volume(self, paramValues):
        rg = numpy.array((self.rg(),))
        if self.rg.isActive():
            rg = paramValues[:, 0]
        k = numpy.array((self.k(),))
        if self.k.isActive():
            idx = sum([int(p.isActive()) for p in self.rg, self.bp, self.etas])
            k = paramValues[:, idx]
        v = k * rg**2
        if len(v) <= 1 and len(v) < len(paramValues):
            # duplicates of single value for dimension compatibility
            v = numpy.ones(len(paramValues)) * v
        return v**self.compensationExponent

GaussianChain.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    # FIXME: use SASData.load() instead
    pf = PDHFile("../brianpauw/sasfit_gauss2-5-1.5-2-1.dat")
    model = GaussianChain()
    model.rg.setValue(5.)
    model.rg.setActive(False)
    model.bp.setValue(1.5)
    model.bp.setActive(False)
    model.etas.setValue(1.)
    model.etas.setActive(False)
    model.k.setValue(0.08)
    model.k.setActive(False)
    intensity = model.formfactor(pf.data, None).reshape(-1)**2.
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("gauss.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
