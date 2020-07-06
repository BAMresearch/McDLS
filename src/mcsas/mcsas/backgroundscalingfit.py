# -*- coding: utf-8 -*-
# mcsas/backgroundscalingfit.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from scipy import optimize
from ..bases.model import SASModel

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
    :arg positiveBackground: *(optional)* Fixes the background values to positive only,
                     default: False

    :returns: (*sc*, *conval*): A tuple of an array containing the
              intensity scaling factor and background and the reduced
              chi-squared value.
    """
    _findBackground = None # True: find optimal background as well
    _positiveBackground = None # True: Fix background to positive values only

    def __init__(self, findBackground, positiveBackground, *args):
        self._findBackground = bool(findBackground)
        self._positiveBackground = bool(positiveBackground)

    def signBackground(self, b):
        """ Fixes the background to positive values if bool is set """
        if self._positiveBackground:
            return abs(b)
        else:
            return b

    @staticmethod
    def chi(sc, dataMeas, dataErr, dataCalc):
        """Chi calculation, difference of measured and calculated signal.
        """
        return (dataMeas - sc[0] * dataCalc - sc[1]) / dataErr
    
    @staticmethod
    def chiPosBg(sc, dataMeas, dataErr, dataCalc):
        """Chi calculation, difference of measured and calculated signal.
        """
        return (dataMeas - sc[0] * dataCalc - abs(sc[1])) / dataErr

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

    @staticmethod
    def aGoFsAlpha(dataMeas, dataErr, dataCalc):
        """The alternative Goodness-of-Fit value without alpha, i.e. multiplied
        by alpha, according to [Henn 2016]
        ( http://dx.doi.org/10.1107/S2053273316013206 )."""
        return sum((dataMeas - dataCalc)**2) / sum(dataErr**2)

    def dataScaled(self, data, sc):
        """Returns the input data scaled by the provided factor and background
        level applied if requested."""
        if self._findBackground:
            return (data * sc[0]) + self.signBackground(sc[1]) # apply background on request
        # else:
        return (data * sc[0])

    def fitLM(self, dataMeas, dataErr, dataCalc, sc):
        """Levenberg-Marquardt method"""
        func = self.chiNoBg
        if self._findBackground:
            func = self.chi
            if self._positiveBackground:
                func = self.chiPosBg
        sc, success = optimize.leastsq(func, sc, args = (dataMeas, dataErr, dataCalc),
                                       full_output = False)
        return sc

    def fitSimplex(self, dataMeas, dataErr, dataCalc, sc):
        """Downhill Simplex (aka Nelder-Mead) method"""
        def residual(xsc):
            return self.chiSqr(dataMeas, dataErr, self.dataScaled(dataCalc, xsc))
        sc = optimize.fmin(residual, sc, full_output = False, disp = 0)
        return sc

    def calc(self, data, modelData, sc, ver = 2):
        dataMeas = data.f.binnedData.flatten()
        dataErr = 1.
        if data.f.binnedDataU is not None:
            dataErr = data.f.binnedDataU.flatten()
            dataErr[dataErr == 0.0] = 1. # prevent division by zero
        dataCalc = modelData.chisqrInt
        if not len(dataMeas): # all data filtered
            return sc, 1., dataCalc, 1.

        # different data fit approaches: speed vs. stability (?)
        # find 2 values: scaling & background
        # least squares optimize them to match model data with measurement
        if ver == 2:
            sc = self.fitLM(dataMeas, dataErr, dataCalc, sc)
        else:
            sc = self.fitSimplex(dataMeas, dataErr, dataCalc, sc)

        if not self._findBackground:
            sc[1] = 0.0
        else:
            sc[1] = self.signBackground(sc[1])
        # calculate convergence value
        conval = self.chiSqr(dataMeas, dataErr, self.dataScaled(dataCalc, sc))
        aGoFs = self.aGoFsAlpha(dataMeas, dataErr, self.dataScaled(dataCalc, sc))
        # multiplied by the reciprocal of alpha
        aGoFs *= len(dataMeas) / (len(dataMeas) - modelData.numParams )
        return sc, conval, dataCalc, aGoFs

# vim: set ts=4 sts=4 sw=4 tw=0:
