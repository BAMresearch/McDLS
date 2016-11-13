# -*- coding: utf-8 -*-
# models/kholodenko.py

from __future__ import absolute_import
import logging
import numpy
from scipy.special import j1 as bessel_j1
from scipy.integrate import quad
from utils.parameter import FitParameter
from .scatteringmodel import SASModel
from bases.algorithm import RandomUniform, RandomExponential
from utils.units import Length

LASTMSG = set()

def core(z, qValue, kuhnLength, x):
    if z <= 0.0 or x <= 0.0:
        return 1.0
    ratio = 3.0 / kuhnLength
    if qValue < ratio:
        e = numpy.sqrt(1.0 - qValue*qValue*kuhnLength*kuhnLength/9.)
        fz = numpy.sinh(e*z) / (e*numpy.sinh(z))
    elif qValue > ratio:
        f = numpy.sqrt(qValue*qValue*kuhnLength*kuhnLength/9. - 1.0)
        fz = numpy.sin(f*z)  / (f*numpy.sinh(z))
    else: # qValue == ratio
        fz = z / numpy.sinh(z)
    res = fz * (2./x) * (1.0 - z/x)
    return res

def coreIntegral(qValue, kuhnLength, x):
    res = quad(core, 0, x,
               args = (qValue, kuhnLength, x),
               limit = 10000, full_output = 1, epsabs = 0.0, epsrel = 1e-10)
    if len(res) > 3:
        LASTMSG.add(res[-1])
    return numpy.sqrt(res[0])
vectorizedCoreIntegral = numpy.vectorize(coreIntegral)

def calcPcs(u):
    if u <= 0.0:
        return 1.0
    res = 2. * bessel_j1(u) / u
    return res
vectorizedPcs = numpy.vectorize(calcPcs)

class Kholodenko(SASModel):
    r"""Form factor of a worm-like structure after [Kholodenko93]_

    .. [Kholodenko93] `A. L. Kholodenko. Analytical calculation of the
        scattering function for polymers of arbitrary flexibility using the
        dirac propagator. Macromolecules, 26:4179–4183, 1993.
        <http://dx.doi.org/10.1021/ma00068a017>`_
    """
    shortName = "Kholodenko Worm"
    parameters = (
            FitParameter("radius", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "Radius",
                    generator = RandomExponential,
                    valueRange = (0., numpy.inf),
                    activeRange = Length(u'nm').toSi((1, 5))), # preset
            FitParameter("lenKuhn", Length(u'nm').toSi(1.), unit = Length(u'nm'),
                    displayName = "kuhn length",
                    generator = RandomUniform,
                    valueRange = (0., numpy.inf),
                    activeRange = Length(u'nm').toSi((10, 50))), # preset
            FitParameter("lenContour", Length(u'nm').toSi(2.), unit = Length(u'nm'),
                    displayName = "contour length",
                    generator = RandomUniform,
                    valueRange = (0., numpy.inf),
                    activeRange = Length(u'nm').toSi((100, 1000))), # preset
    )

    def __init__(self):
        super(Kholodenko, self).__init__()
        # some presets of parameters to fit
        self.radius.setActive(True)
        self.lenKuhn.setActive(True)
        self.lenContour.setActive(True)

    def formfactor(self, dataset):
        # vectorized data and arguments
        qr = dataset.q * self.radius() # a vector
        pcs = vectorizedPcs(qr)

        x = 3. * self.lenContour() / self.lenKuhn()
        p0 = vectorizedCoreIntegral(dataset.q, self.lenKuhn(), x)
        if len(LASTMSG):
            logging.warning("\n".join(["numpy.quad integration messages: "] + list(LASTMSG)))
        return p0 * pcs # non-squared as opposed to SASfit

    def volume(self):
        volume = numpy.pi * self.lenContour() * self.radius()**2
        return volume

Kholodenko.factory()

def test():
    # volume is already included in the model, how to exclude it?
    Kholodenko.testVolExp = 0.0
    for fn in ("sasfit_kho-1-10-1000.dat",
               ):
        yield Kholodenko.test, fn

# vim: set ts=4 sts=4 sw=4 tw=0:
