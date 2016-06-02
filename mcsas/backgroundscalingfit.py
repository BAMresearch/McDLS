# -*- coding: utf-8 -*-
# mcsas/backgroundscalingfit.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from __future__ import absolute_import # PEP328
from scipy import optimize
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
    _dataCalcSquared = False

    def __init__(self, findBackground, model):
        self.findBackground = findBackground
        self._dataCalcSquared = not isinstance(model, SASModel)

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

    def calc(self, data, modelData, sc, ver = 2):
        dataMeas = data.f.binnedData.flatten()
        dataErr  = data.f.binnedDataU.flatten()
        dataCalc = modelData.chisqrInt

        # different data fit approaches: speed vs. stability (?)
        if ver == 2:
            sc = self.fitLM(dataMeas, dataErr, dataCalc, sc)
        else:
            sc = self.fitSimplex(dataMeas, dataErr, dataCalc, sc)

        if not self.findBackground:
            sc[1] = 0.0
        # calculate convergence value
        conval = self.chiSqr(dataMeas, dataErr, self.dataScaled(dataCalc, sc))
        return sc, conval, dataCalc

# vim: set ts=4 sts=4 sw=4 tw=0:
