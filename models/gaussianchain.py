# -*- coding: utf-8 -*-
# models/gaussianchain.py

import numpy
from numpy import pi
from mcsas.McSAS import ScatteringModel
from cutesnake.algorithm import Parameter, RandomUniform
from scatteringmodel import ScatteringModel

# parameters must not be inf

class GaussianChain(ScatteringModel):
    r"""Form factor of flexible polymer chains which are not selfavoiding
    and obey Gaussian statistics after [Debye47]_

    See also: http://sasfit.sf.net/manual/Gaussian_Chain#Gauss_2

    .. [Debye47] `P. Debye, Mollecular-weight determination by light
        scattering, Journal of Physical and Colloid Chemistry, 51:18--32,
        1947. <http://dx.doi.org/10.1021/j150451a002>`_
    """
    shortName = "Gaussian Chain"
    parameters = (
            Parameter("rg", 1.0,
                    displayName = "radius of gyration",
                    valueRange = (0., numpy.inf), suffix = "nm"),
            Parameter("bp", 100.0,
                    displayName = "scattering length of the polymer",
                    valueRange = (0., numpy.inf), suffix = "cm"),
            Parameter("etas", 1.0,
                    displayName = "scattering length density of the solvent",
                    valueRange = (0., numpy.inf), suffix = "cm<sup>-1</sup>"),
            Parameter("volume", 1.0,
                    displayName = "molecular volume of a single polymer molecule",
                    valueRange = (0., numpy.inf), suffix = "cm<sup>3</sup>")
    )
    parameters[0].isActive = True
    parameters[3].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        # some presets
        self.rg.setValueRange((0.1, 1e3))
        self.bp.setValueRange((0.1, 1e3))
        self.etas.setValueRange((0.1, 10.))
        self.volume.setValueRange((1.0, 1e3))

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)

        # vectorized data and arguments
        q = dataset[:, 0]
        rg = numpy.array((self.rg.value(),))
        if self.rg.isActive:
            rg = paramValues[:, 0]
        bp = numpy.array((self.bp.value(),))
        if self.bp.isActive:
            idx = int(self.rg.isActive)
            bp = paramValues[:, idx]
        etas = numpy.array((self.etas.value(),))
        if self.etas.isActive:
            idx = int(self.rg.isActive) + int(self.bp.isActive)
            etas = paramValues[:, idx]
        volume = numpy.array((self.volume.value(),))
        if self.volume.isActive:
            idx = int(self.rg.isActive) + int(self.bp.isActive) + int(self.etas.isActive)
            volume = paramValues[:, idx]

        u = numpy.outer(q, rg)**2 # a matrix usually
        beta = bp - volume * etas
        beta = beta*beta * 2.0
        result = (numpy.expm1(-u) + u) / (u*u)
        for i, res in zip(range(0, result.shape[0]), result):
            if q[i] <= 0.0:
                result[i] = beta
            result[i] = res * beta
        return result

    def vol(self, paramValues, compensationExponent = None):
        assert ScatteringModel.vol(self, paramValues)
        if compensationExponent is None:
            compensationExponent = self.compensationExponent
        idx = 0
        for p in self.rg, self.bp, self.etas:
            idx += int(p.isActive)
        v = numpy.array((self.volume.value(),))
        if self.volume.isActive:
            v = paramValues[:, idx]
        if len(v) <= 1 and len(v) < len(paramValues):
            # duplicates of single value for dimension compatibility
            v = numpy.ones(len(paramValues)) * v
        return v**compensationExponent

GaussianChain.factory()

if __name__ == "__main__":
    from cutesnake.datafile import PDHFile, AsciiFile
    pf = PDHFile("sasfit_gauss2-1-100-1-1.dat")
    model = GaussianChain()
    model.rg.setValue(1.)
    model.rg.isActive = False
    model.bp.setValue(100.)
    model.bp.isActive = False
    model.etas.setValue(1.)
    model.etas.isActive = False
    model.volume.setValue(1.)
    model.volume.isActive = False
    intensity = model.ff(pf.data, None).reshape(-1)
    q = pf.data[:, 0]
    oldInt = pf.data[:, 1]
    delta = abs(oldInt - intensity)
    result = numpy.dstack((q, intensity, delta))[0]
    AsciiFile.writeFile("gauss.dat", result)
    # call it like this:
    # PYTHONPATH=..:../mcsas/ python brianpauwgui/gaussianchain.py && gnuplot -p -e 'set logscale xy; plot "gauss.dat" using 1:2:3 with errorbars'

# vim: set ts=4 sts=4 sw=4 tw=0:
