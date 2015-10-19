# -*- coding: utf-8 -*-
# mcsas/backgroundscalingfit.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from __future__ import absolute_import # PEP328

from scipy import optimize
import numpy
import warnings
warnings.simplefilter("ignore")
from lmfit import Parameter as lmfit_Parameter
from lmfit import minimize as lmfit_minimize
import sys
from models import SASModel

class BackgroundScalingFit(object):
    """
    Chi-squared convergence calculation happens here.
    Optimizes the scaling and background factor to match *intCalc* closest
    to intObs. 
    Returns an array with scaling factors. Input initial guess *sc* has 
    to be a two-element array with the scaling and background.

    **Input arguments:**

    :arg intObs: An array of *measured*, respectively *observed*
                 intensities
    :arg intCalc: An array of intensities which should be scaled to match
                  *intObs*
    :arg intError: An array of uncertainties to match *intObs*
    :arg sc: A 2-element array of initial guesses for scaling
             factor and background
    :arg ver: *(optional)* Can be set to 1 for old version, more robust
              but slow, default 2 for new version,
              10x faster than version 1, requires decent starting values
    :arg outputIntensity: *(optional)* Return the scaled intensity as
                          third output argument, default: False
    :arg background: *(optional)* Enables a flat background contribution,
                     default: True

    :returns: (*sc*, *conval*): A tuple of an array containing the
              intensity scaling factor and background and the reduced
              chi-squared value.
    """
    findBackground = None # True: find optimal background as well
    _model = None

    def __init__(self, findBackground, model):
        self.findBackground = findBackground
        self._model = model

    @staticmethod
    def chi(sc, dataMeas, dataErr, dataCalc):
        """Chi calculation, difference of measured and calculated signal.
        """
        return (dataMeas - sc[0] * dataCalc - sc[1]) / dataErr
    
    @staticmethod
    def chiNoBg(sc, dataMeas, dataErr, dataCalc):
        """Chi calculation, difference of measured and calculated signal,
        scaling only, no backgrund.
        """
        return (dataMeas - sc[0] * dataCalc) / dataErr

    @staticmethod
    def chiSqr(dataMeas, dataErr, dataCalc):
        """Reduced Chi-squared calculation, size of parameter-space not taken
        into account; for data with known intError.
        """
        return sum(((dataMeas - dataCalc)/dataErr)**2) / len(dataMeas)

    def dataScaled(self, data, sc):
        """Returns the input data scaled by the provided factor and background
        level applied if requested."""
        if self.findBackground:
            return (data * sc[0]) + sc[1] # apply background on request
        # else:
        return (data * sc[0])

# we need following external interface:
# - be able to call v1 & v2
# - specify volume separately in order to ignore it for DLS purposes
# - be able to call w/ or w/o volume
# -> 4 function: v1, v1Vol, v2, v2Vol
# error is associated to measured data, should be given next to it
# TODO: constraints!

    def fitLM(self, dataMeas, dataErr, dataCalc, sc):
        func = self.chiNoBg
        if self.findBackground:
            func = self.chi
        sc, success = optimize.leastsq(func, sc, args = (dataMeas, dataErr, dataCalc), full_output = False)
        return sc

    def fitSimplex(self, dataMeas, dataErr, dataCalc, sc):
        def residual(xsc):
            return self.chiSqr(dataMeas, dataErr, self.dataScaled(dataCalc, xsc))
        sc = optimize.fmin(residual, sc, full_output = False, disp = 0)
        return sc

    @staticmethod
    def params2sc(params):
        return numpy.array([params['scaling'].value,
                            params['background'].value])

    def fitLMConstrained(self, dataMeas, dataErr, dataCalc, sc):

        def residual(params, dMeas, dErr, dCalc):
            return self.chi(self.params2sc(params), dMeas, dErr, dCalc)

        # provide parameter as tuple improves performance in lmfit
        # (avoids costly deepcopy of Parameters)
        scaling = lmfit_Parameter('scaling', value = sc[0], min = 0.0, max = 1e4)
        background = lmfit_Parameter('background', value = sc[1], min = 0.0, max = 0.2)
        if not self.findBackground:
            background.set(value = 0.0, vary = False)

        out = lmfit_minimize(residual, (scaling, background),
                             args = (dataMeas, dataErr, dataCalc))
        return self.params2sc(out.params)

    # slower compared to the lmfit implementation, 2x for SAXS data
    def fitBFGSConstrained(self, dataMeas, dataErr, dataCalc, sc):
        def residual(xsc, dMeas, dErr, dCalc):
            return self.chiSqr(dMeas, dErr, self.dataScaled(dCalc, xsc))
        out = optimize.fmin_l_bfgs_b(residual, sc,
                                     args = (dataMeas, dataErr, dataCalc),
                                     bounds = [(0.0, 1e4), (0.0, 0.2)],
                                     approx_grad = True, disp = 0)
        return out[0]

    def calc(self, dataMeas, dataErr, dataCalc, sc, vol = None, ver = 2):
        dataMeas = dataMeas.flatten()
        dataErr  = dataErr.flatten()
        dataCalc = dataCalc.flatten()

        if vol is not None:
            if isinstance(self._model, SASModel):
#                print >>sys.__stderr__, "SASvol", vol
                dataCalc /= vol

        if True:
            # different data fit approaches: speed vs. stability (?)
            if ver == 2:
                sc = self.fitLM(dataMeas, dataErr, dataCalc, sc)
            else:
                sc = self.fitSimplex(dataMeas, dataErr, dataCalc, sc)
        else:
            #sc = self.fitLMConstrained(dataMeas, dataErr, dataCalc, sc)
            sc = self.fitBFGSConstrained(dataMeas, dataErr, dataCalc, sc)

        if not self.findBackground:
            sc[1] = 0.0
        # calculate convergence value
        conval = self.chiSqr(dataMeas, dataErr, self.dataScaled(dataCalc, sc))
        return sc, conval, dataCalc

# vim: set ts=4 sts=4 sw=4 tw=0:
