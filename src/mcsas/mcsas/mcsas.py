# -*- coding: utf-8 -*-
# mcsas/mcsas.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

import numpy # For arrays
from numpy import (inf, array, reshape, shape, pi, diff, zeros,
                  size, sum, sqrt, log10,
                  isnan, newaxis)
from scipy.integrate import simps
# useful for debugging numpy RuntimeWarnings
# numpy.seterr(all = "raise", under = "ignore")
import time # Timekeeping and timing of objects
import copy
import logging
logging.basicConfig(level = logging.INFO)
import subprocess, os, shutil, tempfile

from ..utils import isList 
from ..bases.dataset import DataSet
from ..bases.algorithm import AlgorithmBase
from ..utils.parameter import isActiveFitParam
from ..utils.tests import isMac
from ..bases.model import ScatteringModel
from ..gui.utils import processEventLoop
from .backgroundscalingfit import BackgroundScalingFit
from ..log import timestampFormatted

from . import PlotResults
from . import McSASParameters
from ..dataobj import SASData

class ConvBuffer(object):
    """Implements kind of a ring buffer to store a limited number
    of (x,y) value pairs and calculate statistics of it continuously."""
    _minConv = None
    _testConvVariance = None
    _y = None
    _x = None
    _size = 20 # number of values for variance or slope calc.
    _huge = 9.9e199 # 1e200 is max in the GUI, but need _size buffer too

    def __init__(self, minConv, testConvVariance = False):
        self._minConv = minConv
        self._testConvVariance = testConvVariance
        self.clear()

    def clear(self):
        self._x = []
        self._y = []

    def __len__(self):
        return len(self._x)

    def set(self, i, value):
        self._x.append(i)
        self._y.append(value)
        if len(self) > self._size: # limit the number of values to remember
            del self._x[0] # remove the oldest, front
            del self._y[0]

    @property
    def key(self):
#        return numpy.abs(self._slope()) # negative slope of chisqr usually
        if len(self._y) > 1:
            return numpy.array(self._y).var() # variance, need at least 1 val
        else:
            return self._huge

    def _slope(self):
        # https://stackoverflow.com/a/9538936
        # ((X*Y).mean(axis=1) - X.mean()*Y.mean(axis=1)) / ((X**2).mean() - (X.mean())**2)
        x = numpy.array(self._x)
        y = numpy.array(self._y)
        # calculate slope x vs. y values
        d = ((x**2).mean() - x.mean()**2)
        if d == 0.:
            return self._huge
        return ((x*y).mean() - x.mean()*y.mean()) / d

    def reached(self, conval):
        if self._testConvVariance:
            # all initial values have to be replaced and
            # the variance of the last conv. values has to be very small
            return (len(self) == self._size and self.key < self._minConv)
        else:
            return conval < self._minConv

class DistribRecorder(object):
    """Records the probability distribution function during the optimization.
    """
    _pdfDir = "test_pdf" # directory created in OS TMP dir
    _enabled = False
    _xrange = None
    _idx = None
    _rejected = None # list of rejected candidate values
    _dataMeas = None
    _dataVec = None # positions at which the measured signal was captured

    def reject(self, value):
        if self._enabled:
            self._rejected.append(value)

    def capture(self, pValues, modelData, numIter, numMoves):
        if not self._enabled:
            return
        def getPDF(values):
            """Calculates the PDF function from a set of values."""
            nBins = 30
            if not len(values):
                zeros = numpy.zeros(nBins)
                return zeros, zeros
            hist, binEdges = numpy.histogram(values, bins = nBins,
                                             density = True)
            binCenter = binEdges[:-1] + numpy.diff(binEdges)*.5
            return binCenter, hist
        if self._idx is None:
            self._idx = 0
        # set up output file names
        basename = "r{:03d}_{:05d}".format(self._numRep, self._idx)
        fnData = os.path.join(self._pdfDir, basename+".txt")
        fnImg  = os.path.join(self._pdfDir, basename+".png")
        # get the PDF function
        binsAcc, accepted = getPDF(self._param.unit().toDisplay(
                                   pValues[:,0].flatten()))
        binsRej, rejected = getPDF(self._param.unit().toDisplay(
                                   numpy.array(self._rejected)))
        pdfData = numpy.vstack((binsAcc, accepted, binsRej, rejected)).T
        # get the measured and model signal
        signalData = numpy.vstack(
                (self._dataVec, self._dataMeas, modelData)).T
        res = numpy.absolute(signalData[:,1]-signalData[:,2])
        # fill up the distribution data to match the signal length
        fillCount = signalData.shape[0] - pdfData.shape[0]
        fillArr = numpy.full((fillCount, pdfData.shape[1]), numpy.nan)
        pdfData = numpy.concatenate((pdfData, fillArr))
        combined = numpy.concatenate((pdfData, signalData), axis = 1)
        # store the distribution and the signal in a common file
        numpy.savetxt(fnData, combined)
        # plot the data from file, create one frame
        # FIXME: set xrange in gnuplot to param min/max
        cmd = ['[ -z "$(which gnuplot)" ] || \
                [ -f \"{fnData}\" ] && \
                gnuplot -e \'\
                set terminal pngcairo size 1280,720 dashed enhanced;\
                set output \"{fnImg}\";\
                set grid;\
                set xrange [{xr}];\
                set yrange [0:*];\
                set multiplot layout 2, 1 \
                    title \"move {move}/{it}, rep {rep}/{numReps}\" noenhanced;\
                plot \"{fnData}\" using 3:4 with lines lc \"red\" lw 2 \
                                    title \"radii rejected\",\
                     \"{fnData}\" using 1:2 with lines lc \"blue\" lw 2 \
                                    title \"radii accepted\";\
                set logscale x; set format xy \"%g\"; set y2tics;\
                unset xrange; set yrange [-0.05:1.05];\
                plot \"{fnData}\" using 5:6 with points pt 7 ps 1\
                        lc \"black\" title \"measured\",\
                     \"{fnData}\" using 5:7 with lines lw 2\
                        lc \"red\" title \"model\",\
                     \"{fnData}\" using 5:(abs($6-$7)) axes x1y2\
                        with lines lw 2\
                        lc \"blue\" title \"residual\",\
                     \"{fnData}\" using 5:({resMean}) axes x1y2 \
                        with line dt \"-\" lw 2\
                        lc \"blue\" title \"mean(residual) = {resMean:.2g}\",\
                     \"{fnData}\" using 5:({resVar}) axes x1y2 \
                        with line dt \"-.\" lw 2\
                        lc \"blue\" title \"var(residual) = {resVar:.2g}\",\
                     \"{fnData}\" using 5:({resArea}) axes x1y2 \
                        with line dt \".\" lw 2\
                        lc \"blue\" title \"area(residual) = {resArea:.2g}\";\
                \' && echo rm -f \"{fnData}\"'
                .format(fnData = fnData, fnImg = fnImg,
                        rep = self._numRep, numReps = self._numReps,
                        move = numMoves, it = numIter,
                        xr = self._xrange,
                        resMean = res.mean(), resVar = res.var(), resArea =
                        simps(res, x = signalData[:,0], even='first'))
                   .replace("\n","")]
#            print("cmd:", cmd)
#            return
        process = subprocess.Popen(cmd, shell = True,
                                   stdout = subprocess.PIPE)
        self._idx += 1
        # rename: num=0; for fn in *.png; do ln -s $fn "$(printf "l%05d.png" $num)"; num=$((num+1)); done
        # video: cat *.png | avconv -framerate 10 -i l%05d.png -threads 4 -c:v libvpx test.webm
        # cleanup remaining files first

    def finalize(self):
        if not os.path.exists(self._pdfDir):
            return
        cmd = ['[ -z "$(which avconv)" ] || (cd "{wd}" && \
                avconv -r 10 -i r{rep:03d}_%05d.png -threads 4 \
                       -c:v libvpx -crf 10 -b:v 2M \
                       {outfn} > /dev/null 2>&1 &)'
                .format(wd = self._pdfDir, rep = self._numRep,
                        outfn = "distrib+fit.webm")]
        subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)

    def clear(self):
        if not os.path.exists(self._pdfDir):
            return
        shutil.rmtree(self._pdfDir)
#        subprocess.Popen("rm -f test_pdf_*", shell = True, stdout = subprocess.PIPE)

    def setup(self):
        if not self._enabled:
            return
        self.clear()
        os.makedirs(self._pdfDir)
        self._xrange = ":".join(("{:g}".format(v)
                                for v in self._param.unit().toDisplay(
                                self._param.activeRange())))
        self._rejected = []

    def __init__(self, params, data, numRep, numReps, enabled = False):
        self._enabled = enabled
        self._numRep  = numRep
        self._numReps = numReps
        self._pdfDir = os.path.join(tempfile.gettempdir(),
                    timestampFormatted()+"_progress_r{:03d}".format(numRep))
        if not len(params):
            self._enabled = False
        if self._enabled:
            self._param = params[0]
            # get measured data for comparison, see BackgroundScalingFit.calc()
            self._dataMeas = data.f.binnedData.flatten()
            self._dataVec  = data.x0.binnedData.flatten()
            self.setup()

class McSAS(AlgorithmBase):
    r"""
    Main class containing all functions required to do Monte Carlo fitting.

    **Required:**

        - *data*: The dataset to fit. 
                   Has to be an instance of :py:class:SASData
        - *model*: The scattering model object to assume.
                   It has to be an instance of :py:class:`ScatteringModel`.

        For more settings, see mcsas/mcsasparameters.json

    **Returns:**

    A McSAS object with the following Results stored in the *result* member
    attribute. These can be extracted using
    McSAS.result[<parameterIndexNumber>]['<Keyword>']
    where the *parameterIndexNumber* indicates which shape parameter 
    information is requested.
    E.g. an ellipsoid has 3: width, height and orientation.
    (Some information is only stored in *parameterIndexNumber = 0* (default)).

    **Keyword** may be one of the following:

        *fitMeasValMean*: 1D array (*common result*)
            The fitted measVal, given as the mean of all numReps Results.
        *fitX0*: 1D array (*common result*)
            Corresponding q values
            (may be different than the input q if *X0Bounds* was used).
        *fitMeasValStd*: array (*common result*)
            Standard deviation of the fitted I(q), calculated as the standard 
            deviation of all numReps results.
        *contribs*: size array (numContribs x numReps) (*common result*)
            Collection of numContribs contributions fitted to best represent 
            the provided I(q) data. Contains the Results of each of 
            *numReps* iterations. This can be used for rebinning without 
            having to re-optimize.
        *scalingFactors*: size array (2 x numReps) (*common result*)
            Scaling and background values for each repetition.
            Used to display background level in data and fit plot.
        *histogramXLowerEdge*: array
            histogram bin left edge position (x-axis in histogram).
        *histogramXMean*: array
            Center positions for the size histogram bins
            (x-axis in histogram, used for errorbars).
        *histogramXWidth*: array
            histogram bin width
            (x-axis in histogram, defines bar plot bar widths).
        *volumeHistogramYMean*: array
            Volume-weighted particle size distribution values for
            all numReps Results (y-axis bar height).
        *numberHistogramYMean*: array
            Number-weighted analogue of the above *volumeHistogramYMean*.
        *volumeHistogramRepetitionsY*: size array (self.histogramBins x numReps)
            Volume-weighted particle size distribution bin values for
            each fit repetition (the mean of which is *volumeHistogramYMean*, 
            and the sample standard deviation is *volumeHistogramYStd*).
        *numberHistogramRepetitionsY*: size array (self.histogramBins x numReps)
            Number-weighted particle size distribution bin values for
            each MC fit repetition.
        *volumeHistogramYStd*: array
            Standard deviations of the corresponding volume-weighted size
            distribution bins, calculated from *numReps* repetitions of the
            model fitting function.
        *numberHistogramYStd*: array
            Standard deviation for the number-weigthed distribution.
        *volumeFraction*: size array (numContribs x numReps)
            Volume fractions for each of numContribs contributions in each of
            *numReps* iterations.
        *numberFraction*: size array (numContribs x numReps)
            Number fraction for each contribution.
        *totalVolumeFraction*: size array (numReps)
            Total scatterer volume fraction for each of the *numReps* 
            iterations.
        *totalNumberFraction*: size array (numReps)
            Total number fraction.
        *minimumRequiredVolume*: size array (numContribs x numReps)
            Minimum required volume fraction for each contribution to become
            statistically significant.
        *minimumRequiredNumber*: size array (numContribs x numReps)
            Number-weighted analogue to *minimumRequiredVolume*.
        *volumeHistogramMinimumRequired*: size array (histogramXMean)
            Array with the minimum required volume fraction per bin to become
            statistically significant. Used to display minimum required level
            in histogram.
        *numberHistogramMinimumRequired*: size array (histogramXMean)
            Number-weighted analogue to *volumeHistogramMinimumRequired*.
        *scalingFactors*: size array (2 x numReps)
            Scaling and background values for each repetition. Used to display
            background level in data and fit plot.
        *totalVolumeFraction*: size array (numReps)
            Total scatterer volume fraction for each of the *numReps*
            iterations.
        *minimumRequiredVolume*: size array (numContribs x numReps)
            Minimum required volube fraction for each contribution to become
            statistically significant.
        *volumeHistogramMinimumRequired*: size array (histogramXMean)
            Array with the minimum required volume fraction per bin to become
            statistically significant. Used to display minimum required level
            in histogram.

    """

    data = None    # user provided data to work with
    model = None
    result = None
    _convBuffer = None # remembers the last few goodness of fit values

    # there are several ways to accomplish this depending on where/when
    # McSASParameters() should be called: when creating an instance or
    # for creating different types with different settings ...
    @classmethod
    def factory(cls):
        # calling factory of parent class which is AlgorithmBase currently
        return super(McSAS, cls).factory("McSAS",
                                         *McSASParameters().parameters)
    def __init__(self, *args, **kwargs):
        super(McSAS, self).__init__(*args, **kwargs)
        # assist the user by setting the minimum required goodness of fit
        # variance to a useful value
        def setDefaultVarianceLimit():
            if self.testConvVariance() and self.convergenceCriterion() > .01:
                self.convergenceCriterion.setValue(1e-14)
        self.testConvVariance.setOnValueUpdate(setDefaultVarianceLimit)

    def calc(self, **kwargs):
        # initialize
        self.result = [] # TODO
        self.stop = False # TODO, move this into some simple result structure

        assert(self.data is not None)

        if self.model is None:
            # verify default model in case no model was provided
            if (McSASParameters.model is None or
                not isinstance(McSASParameters.model, ScatteringModel)):
                # load the model here, because it was reloaded by FindModels
                from models.sphere import Sphere
                McSASParameters.model = Sphere() # create instance
                logging.info("Default model not provided, setting to: {0}"
                        .format(str(McSASParameters.model.name())))
            self.model = McSASParameters.model
        logging.info("Using model: {0}".format(str(self.model.name())))
        if not self.model.paramCount():
            logging.warning("No parameters to analyse given! Breaking up.")
            return
        logging.info(
                "\n".join([u"Analysing parameters: "] +
                    [str(p) + u", active: " + str(isActiveFitParam(p))
                        for p in self.model.params()])
        )
        self.analyse()
        # continue if there are results only
        if not len(self.result):
            return
        self.histogram()

        ## Fix 2D mode
        # if self.data.f.values2d:
        #     # 2D mode, regenerate measVal
        #     # TODO: test 2D mode
        #     self.gen2DMeasVal()

    ######################################################################
    ####################### optimisation Functions #######################
    ######################################################################

    def analyse(self):
        """This function runs the Monte Carlo optimisation a multitude
        (*numReps*) of times. If convergence is not achieved, it will try 
        again for a maximum of *maxRetries* attempts.
        """
        # get settings
        numContribs = self.numContribs()
        numReps = self.numReps()
        if not any([isActiveFitParam(p) for p in self.model.params()]):
            numContribs, numReps = 1, 1
        # find out how many values a shape is defined by:
        contributions = zeros((numContribs, 
            self.model.activeParamCount(), 
            numReps))
        numIter = zeros(numReps)
        scalings = zeros(numReps)
        backgrounds = zeros(numReps)
        times = zeros(numReps)
        contribMeasVal = zeros([1, self.data.count, numReps])

        self._convBuffer = ConvBuffer(self.convergenceCriterion(),
                                      self.testConvVariance())
        # This is the loop that repeats the MC optimization numReps times,
        # after which we can calculate an uncertainty on the Results.
        for nr in range(numReps):
            elapsedStart = time.time() # for tracking elapsed time
            # keep track of how many failed attempts there have been
            nt = 0
            # do that MC thing! 
            convergence = inf
            self._convBuffer.clear()
            while not self._convBuffer.reached(convergence):
                if nt > self.maxRetries():
                    # this is not a coincidence.
                    # We have now tried maxRetries+2 times
                    logging.warning("Could not reach optimization criterion "
                                    "within {0} attempts, exiting..."
                                    .format(self.maxRetries() + 2))
                    if self.showIncomplete():
                        break
                    else:
                        return
                # retry in the case we were unlucky in reaching
                # convergence within MaximumIterations.
                (contributions[:, :, nr], contribMeasVal[:, :, nr],
                 convergence, details) = self.mcFit(
                                numContribs,
                                outputMeasVal = True, outputDetails = True,
                                nRun = nr)
                if any(array(contributions.shape) == 0):
                    break # nothing active, nothing to fit
                if self.stop:
                    logging.warning("Stop button pressed, exiting...")
                    if self.showIncomplete():
                        break
                    else:
                        return
                nt += 1
            # in minutes:
            # keep track of how many iterations were needed to reach converg.
            numIter[nr] = details.get('numIterations', 0)
            scalings[nr] = details.get('scaling', 1.0) 
            backgrounds[nr] = details.get('background', 0) 
            elapsedTime = (time.time() - elapsedStart)
            times[nr] = elapsedTime 

            tottime = times.sum() /60. # total elapsed time in minutes
            avetime = times[times > 0].mean() / 60. # average optimization time
            remtime = (avetime * numReps - tottime) # est. remaining time
            logging.info("finished optimization number {0} of {1}\n"
                    "  total elapsed time: {2} minutes\n"
                    "  average time per optimization {3} minutes\n"
                    "  total time remaining {4} minutes"
                    .format(nr+1, numReps, tottime, avetime, remtime))

        # store in output dict
        scalingsDDoF = 0
        if len(scalings) > 1: # prevent division by zero in numpy.std()
            scalingsDDoF = 1
        self.result.append(dict(
            contribs = contributions, # Rrep
            # what about modelDataMean? ...
            fitMeasValMean = contribMeasVal.mean(axis = 2),
            fitMeasValStd = contribMeasVal.std(axis = 2),
            fitX0 = self.data.x0.binnedData,
            # ... and dataMean
            dataX0 = self.data.x0.binnedData,
            dataMean = self.data.f.binnedData,
            dataStd = self.data.f.binnedDataU,
            # background details:
            scaling = (scalings.mean(),
                       scalings.std(ddof = scalingsDDoF)),
            background = (backgrounds.mean(),
                          backgrounds.std(ddof = scalingsDDoF)),
            times = times,
            # average number of iterations for all repetitions
            numIter = numIter.mean()))

    def mcFit(self, numContribs,
              outputMeasVal = False, outputDetails = False, nRun = None):
        """
        Object-oriented, shape-flexible core of the Monte Carlo procedure.
        Takes optional arguments:

        *outputMeasVal*:
            Returns the fitted measVal besides the Result

        *outputDetails*:
            details of the fitting procedure, number of iterations and so on

        *nRun*: "serial number" of run. Used to store results in parameters
            in the right place

        """
        data = self.data
        rset = numpy.zeros((numContribs, self.model.activeParamCount()))
        compensationExponent = self.compensationExponent()
        details = dict()
        # index of sphere to change. We'll sequentially change spheres,
        # which is perfectly random since they are in random order.
        
        if self.startFromMinimum():
            for idx, param in enumerate(self.model.activeParams()):
                mb = min(param.activeRange())
                if mb == 0: # FIXME: compare with EPS eventually?
                    mb = pi / (data.x0.limit[1])
                rset[:, idx] = numpy.ones(numContribs) * mb * .5
            logging.info("Starting from min {pName} = {pVal:g} {pUnit}."
                 .format(pName = param.name(),
                         pVal = param.unit().toDisplay(rset[:, idx].mean()),
                         pUnit = param.unit().displayMagnitudeName))
        else:
            rset = self.model.generateParameters(numContribs)

        modelData = self.model.calc(data, rset, compensationExponent)
        ft, vset, wset, sset = (modelData.cumInt, modelData.vset,
                                modelData.wset, modelData.sset)
        if any(array(rset.shape) == 0): # no active params, just return model intensity
            return [rset, modelData.chisqrInt, -1, {'scaling': 1, 'background': 0}]

        # Optimize the intensities and calculate convergence criterium
        # generate initial guess for scaling factor and background
        sc = numpy.array((1.0, data.f.limit[0]))
        if len(modelData.chisqrInt) and modelData.chisqrInt.max() != 0.0:
            # avoid numerical errors
            sc[0] = data.f.limit[1] / modelData.chisqrInt.max()
#        sc *= sum(wset)
        bgScalingFit = BackgroundScalingFit(self.findBackground(),
                                            self.fixed1stPoint())
        # for the least squares fit, normalize the intensity by the sum of
        # weights which is << 1 (for SAXS, usually it's the sum
        # of the scatterers volumes), though increasing ft and reducing the
        # scaling sc[0]; when histogramming, this gets reverted
        sc, conval, dummy, dummy2 = bgScalingFit.calc(data, modelData, sc,
                                                      ver = 1)
        # reoptimize with V2, there might be a slight discrepancy in the
        # residual definitions of V1 and V2 which would prevent optimization.
        sc, conval, dummy, dummy2 = bgScalingFit.calc(data, modelData, sc)
        logging.info("Initial Chi-squared value: {0}".format(conval))

        # start the MC procedure
        start = time.time()
        # progress tracking:
        numMoves, numIter, lastUpdate = 0, 0, 0
        # last successful chisqr
        self._convBuffer.clear()
        # running variable indicating which contribution to change
        ri = 0
        ftest = None

        pdfRecorder = DistribRecorder(self.model.activeParams(), data,
                                      nRun+1, self.numReps())
        # plot initial PDF
        pdfRecorder.capture(rset,
            bgScalingFit.dataScaled(modelData.chisqrInt, sc),
            numIter, numMoves)

        #NOTE: keep track of uncertainties in MC procedure through epsilon
        while (len(wset) > 1 and # see if there is a distribution at all
               not self._convBuffer.reached(conval) and
               numIter < self.maxIterations.value() and
               not self.stop):
            rt = self.model.generateParameters()
            # calculate contribution measVal:
            newModelData = self.model.calc(data, rt, compensationExponent)
            # Calculate new total measVal, subtract old measVal, add new:
            oldModelData = self.model.calc(data, rset[ri].reshape((1, -1)),
                                               compensationExponent)
            testModelData = self.model.getModelData(
                # is numerically stable (so far). Can calculate final uncertainty
                # based on number of valid "moves" and sys.float_info.epsilon
                ft - oldModelData.cumInt + newModelData.cumInt,
                vset,
                # not as intended but sufficient for now
                wset.sum() - wset[ri] + newModelData.wset,
                sset) # surface from testModelData is not used
#            ftest = (ft - oldModelData.cumInt + newModelData.cumInt)
#            wtest = wset.sum() - wset[ri] + newModelData.wset
            # optimize measVal and calculate convergence criterium
            # using version two here for a >10 times speed improvement
            sct, convalt, dummy, aGoFs = bgScalingFit.calc(
                                                    data, testModelData, sc)
            # test if the radius change is an improvement:
            if convalt < conval: # it's better
                # replace current settings with better ones
                rset[ri], sc, conval = rt, sct, convalt
                ft, wset[ri] = testModelData.cumInt, newModelData.wset
                # updating unused data for completeness as well
                vset[ri], sset[ri] = newModelData.vset, newModelData.sset
                self._convBuffer.set(numMoves, conval)
                # output current distribution
                pdfRecorder.capture(rset,
                    bgScalingFit.dataScaled(testModelData.chisqrInt, sc),
                    numIter, numMoves)
                text = "/{0:.3g}".format(self.convergenceCriterion())
                if self.testConvVariance():
                    text = ", var= {0:.3g}".format(self._convBuffer.key) + text
                logging.info("rep {rep}/{reps}, good iter {it}: "
                             "Chisqr= {cs:.3g}{t}, aGoFs= {opt:.3g}\r"
                             .format(it = numIter, cs = conval,
                                     t = text, rep = nRun+1,
                                     reps = self.numReps(), opt = aGoFs))
                numMoves += 1
            else:
                pdfRecorder.reject(rt) # store the rejected rset candidate

            if time.time() - lastUpdate > 0.25:
                # update twice a sec max -> speedup for fast models
                # because output takes much time especially in GUI
                # TODO: don't need this:
                # if we calc in separate processes (multiprocessing) the gui 
                # would have its own thread and will not be blocked by calc

                # process events, check for user input
                processEventLoop()
                lastUpdate = time.time()
            # move to next contribution in list, loop if last contribution
            ri = (ri + 1) % numContribs
            numIter += 1 # add one to the iteration number

        #print # for progress print in the loop
        if numIter >= self.maxIterations.value():
            logging.warning("Exited due to max. number of iterations ({0}) "
                            "reached".format(numIter))
        else:
            logging.info("normal exit")

        # Post-MC operations:
        # the +0.001 prevents a divide by zero error on some Windows systems.
        elapsed = time.time() - start + 1e-3
        logging.info("Number of iterations per second: {0}".format(
                        numIter/elapsed))
        logging.info("Number of valid moves: {0} out of {1} iterations."
                     .format(numMoves, numIter))
        chisqrInfo = ("Final Chi-squared value: {0}".format(conval))
        if self.testConvVariance():
            chisqrInfo += ", var= {0:.3g}".format(self._convBuffer.key)
        logging.info(chisqrInfo)
        details.update({'numIterations': numIter,
            'numMoves': numMoves,
            'elapsed': elapsed})

        modelData = self.model.getModelData(ft, vset, wset, sset)
        sc, conval, ifinal, dummy = bgScalingFit.calc(data, modelData, sc)
        details.update({'scaling': sc[0], 'background': sc[1]})

        result = [rset]
        if outputMeasVal:
            result.append((ifinal * sc[0] + sc[1]))
        result.append(conval)
        if outputDetails:
            result.append(details)
        if nRun is not None:
            #store results in parameter
            for idx, param in enumerate(self.model.activeParams()):
                param.setActiveVal(rset[:, idx], index = nRun) 
        pdfRecorder.finalize()
        # returning <rset, measVal, conval, details>
        return result

    #####################################################################
    #################### Post-optimisation Functions ####################
    #####################################################################

    def histogram(self, contribs = None):
        """
        Takes the *contribs* result from the :py:meth:`McSAS.analyse` function
        and calculates the corresponding volume- and number fractions for each
        contribution as well as the minimum observability limits. It will
        subsequently bin the Result across the range for histogramming 
        purposes.

        While the volume-weighted distribution will be in absolute units
        (providing volume fractions of material within a given size range),
        the number distributions have been normalized to 1.
        
        Output a list of dictionaries with one dictionary per shape parameter:

            *histogramXLowerEdge*: array
                histogram bin left edge position (x-axis in histogram)
            *histogramXMean*: array
                Center positions for the size histogram bins
                (x-axis in histogram, used for errorbars)
            *histogramXWidth*: array
                histogram bin width (x-axis in histogram,
                defines bar plot bar widths)
            *volumeHistogramYMean*: array
                Volume-weighted particle size distribution values for
                all *numReps* Results (y-axis bar height)
            *numberHistogramYMean*: array
                Number-weighted analogue of the above *volumeHistogramYMean*
            *volumeHistogramRepetitionsY*: size (histogramBins x numReps) 
                array Volume-weighted particle size distribution bin values for 
                each MC fit repetition (whose mean is *volumeHistogramYMean*, 
                and whose sample standard deviation is *volumeHistogramYStd*)
            *numberHistogramRepetitionsY*: size (histogramBins x numReps) 
                array Number-weighted particle size distribution bin values
                for each MC fit repetition
            *volumeHistogramYStd*: array
                Standard deviations of the corresponding volume-weighted size
                distribution bins, calculated from *numReps* repetitions of
                the model fitting function
            *numberHistogramYStd*: array
                Standard deviation for the number-weigthed distribution
            *volumeFraction*: size (numContribs x numReps) array
                Volume fractions for each of numContribs contributions 
                in each of numReps iterations
            *numberFraction*: size (numContribs x numReps) array
                Number fraction for each contribution
            *totalVolumeFraction*: size (numReps) array
                Total scatterer volume fraction for each of the *numReps*
                iterations
            *totalNumberFraction*: size (numReps) array
                Total number fraction 
            *minimumRequiredVolume*: size (numContribs x numReps) array
                minimum required volume fraction for each contribution to
                become statistically significant.
            *minimumRequiredNumber*: size (numContribs x numReps) array
                number-weighted analogue to *minimumRequiredVolume*
            *volumeHistogramMinimumRequired*: size (histogramXMean) array 
                array with the minimum required volume fraction per bin to
                become statistically significant. Used to display minimum
                required level in histogram.
            *numberHistogramMinimumRequired*: size (histogramXMean) array
                number-weighted analogue to *volumeHistogramMinimumRequired*
            *scalingFactors*: size (2 x numReps) array
                Scaling and background values for each repetition. Used to
                display background level in data and fit plot.
        """
        if not isList(self.result) or not len(self.result):
            logging.info("There are no results to histogram, breaking up.")
            return
        if contribs is None:
            contribs = self.result[0]['contribs']
        if not all(numpy.array(contribs.shape, dtype = bool)):
            # testing for any active parameters (contribs[1])
            logging.info("Nothing to histogram, giving up.")
            return
        numContribs, dummy, numReps = contribs.shape

        # volume fraction for each contribution
        volumeFraction = zeros((numContribs, numReps))
        # number fraction for each contribution
        numberFraction = zeros((numContribs, numReps))
        volSqrFraction = zeros((numContribs, numReps)) # aka intensity
        surfaceFraction = zeros((numContribs, numReps))
        # volume frac. for each histogram bin
        minReqVol = zeros((numContribs, numReps)) 
        # number frac. for each histogram bin
        minReqNum = zeros((numContribs, numReps))
        minReqVolSqr = zeros((numContribs, numReps))
        minReqSurface = zeros((numContribs, numReps))
        totalVolumeFraction = zeros((numReps))
        totalNumberFraction = zeros((numReps))
        totalVolSqrFraction = zeros((numReps))
        totalSurfaceFraction = zeros((numReps))
        # MeasVal scaling factors for matching to the experimental
        # scattering pattern (Amplitude A and flat background term b,
        # defined in the paper)
        scalingFactors = zeros((2, numReps))

        # data, store it in result too, enables to postprocess later
        # store the model instance too
        data = self.data
        bgScalingFit = BackgroundScalingFit(self.findBackground(),
                                            self.fixed1stPoint())
        # calc vol/num fraction and scaling factors for each repetition
        for ri in range(numReps):
            rset = contribs[:, :, ri] # single set of R for this calculation
            # compensated volume for each sphere vset:
            modelData = self.model.calc(data, rset, self.compensationExponent())
            if not len(modelData.cumInt):
                continue
            ## TODO: same code than in mcfit pre-loop around line 1225 ff.
            # initial guess for the scaling factor.
            sc = numpy.array([data.f.limit[1] / modelData.chisqrInt.max(), data.f.limit[0]])
            # optimize scaling and background for this repetition
            sc, conval, dummy, dummy2 = bgScalingFit.calc(data, modelData, sc)
            scalingFactors[:, ri] = sc # scaling and bgnd for this repetition.
            # calculate individual volume fractions:
            # here, the weight reverts intensity normalization effecting the
            # scaling sc[0] during optimization, it does not influence
            # the resulting volFrac
            volumeFraction[:, ri] = modelData.volumeFraction(sc[0])
            totalVolumeFraction[ri] = sum(volumeFraction[:, ri])
            numberFraction[:, ri] = volumeFraction[:, ri]/modelData.vset.flatten()
            totalNumberFraction[ri] = sum(numberFraction[:, ri])
            volSqrFraction[:, ri] = volumeFraction[:, ri]*modelData.vset.flatten()
            totalVolSqrFraction[ri] = sum(volSqrFraction[:, ri])
            surfaceFraction[:, ri] = numberFraction[:, ri]*modelData.sset.flatten()
            totalSurfaceFraction[ri] = sum(surfaceFraction[:, ri])

            # calc observability for each sphere/contribution
            for c in range(numContribs):
                # observability: the maximum contribution for
                # that sphere to the total scattering pattern
                # NOTE: no need to compensate for p_c here, we work with
                # volume fraction later which is compensated by default.
                # additionally, we actually do not use this value.
                # again, partial intensities for this size only required
                partialModelData = self.model.calc(data, rset[c].reshape((1, -1)),
                                                   self.compensationExponent())
                # dividing by zero tends to go towards infinity,
                # when chosing the minimum those can be ignored
                weightedInt = data.f.binnedDataU * volumeFraction[c, ri]
                partialCumIntScaled = sc[0] * partialModelData.chisqrInt
                indices = (partialCumIntScaled != 0.)
                minReqVol[c, ri] = (
                    weightedInt[indices] / partialCumIntScaled[indices]).min()
                minReqNum[c, ri] = minReqVol[c, ri] / modelData.vset[c]
                minReqVolSqr[c, ri] = (minReqNum[c, ri]
                        * minReqVol[c, ri] * minReqVol[c, ri])
                minReqSurface[c, ri] = (minReqNum[c, ri] * modelData.sset[c])

            if 0 != totalNumberFraction[ri]:
                numberFraction[:, ri] /= totalNumberFraction[ri]
                minReqNum[:, ri]      /= totalNumberFraction[ri]
            if 0 != totalVolSqrFraction[ri]:
                volSqrFraction[:, ri] /= totalVolSqrFraction[ri]
                minReqVolSqr[:, ri]   /= totalVolSqrFraction[ri]
            if 0 != totalSurfaceFraction[ri]:
                surfaceFraction[:, ri] /= totalSurfaceFraction[ri]
                minReqSurface[:, ri]   /= totalSurfaceFraction[ri]

        fractions = dict(vol = (volumeFraction, minReqVol),
                         num = (numberFraction, minReqNum),
                         int = (volSqrFraction, minReqVolSqr),
                         surf = (surfaceFraction, minReqSurface))

        # now we histogram over each variable
        # for each variable parameter we define,
        # we need to histogram separately.
        for paramIndex, param in enumerate(self.model.activeParams()):
            param.histograms().calc(contribs, paramIndex, fractions) # new method

    def gen2DMeasVal(self):
        """
        This function is optionally run after the histogram procedure for
        anisotropic images, and will calculate the MC fit measVal in
        image form
        """
        contribs = self.result[0]['contribs']
        numContribs, dummy, numReps = contribs.shape

        # load original Dataset
        x0 = data.x0.binnedData
        # we need to recalculate the result in two dimensions
        kansas = shape(q) # we will return to this shape
        x0 = x0.sanitized.flatten()

        logging.info("Recalculating 2D measVal, please wait")
        # for each Result
        intAvg = zeros(shape(x0))
        # TODO: for which parameter?
        scalingFactors = self.result[0]['scalingFactors']
        for ri in range(numReps):
            logging.info('regenerating set {} of {}'.format(ri, numReps-1))
            rset = contribs[:, :, ri]
            # calculate their form factors
            # ft, vset, wset = self.model.calc(data, rset, compensationExponent)
            modelData = self.model.calc(data, rset, compensationExponent)
            # Optimize the intensities and calculate convergence criterium
            intAvg = (intAvg + modelData.chisqrInt*scalingFactors[0, ri]
                             + scalingFactors[1, ri])
        # print "Initial conval V1", Conval1
        intAvg /= numReps
        # mask (lifted from clipDataset)
        intAvg = intAvg[data.x0.validIndices]
        # shape back to imageform
        self.result[0]['measVal2d'] = reshape(intAvg, kansas)

    def plot(self, axisMargin = 0.3,
             outputFilename = None, autoClose = False):
        """Expects outputFilename to be of type gui.calc.OutputFilename."""
        # the interactive plot figure usually blocks the app until it is closed
        # -> no batch processing possible
        # -> we have to call matplotlib plot in another thread
        # on linux it does not block, can show multiple figures
        # passing the model only, keeping the object hierarchy which is pickled
        # as small as possible
        # TODO: fix pickle support for Models, pickling Sphere() raises the
        # memoize() AssertionsError on Windows
        # It does not happen without the special __init__() constructor of Sphere
        # It does not like circular references (param in hist and hists in param),
        # reminder: histograms contain Parameter references too
        # http://bytes.com/topic/python/answers/37656-assertionerror-pickles-memoize-function#post141263
        # https://stackoverflow.com/questions/23706736/assertion-error-madness-with-qt-and-pickle
        # "autoClose" closes figure after saving for runs of many files. 

        # remove circular parameter references for pickling/forwarding
        # need a list of histograms only, with params for meta info
        # soon to be replaced by HDF parsing (TODO)
        histograms = []
        for p in self.model.activeParams():
            # calls FitParameter.__init__() which sets histogram.param
            newParam = p.copy()
            newParam.setActive(False) # remove the old histograms
            for h in p.histograms():
                newHist = copy.copy(h) # new hist containing old param
                newHist.param = newParam
                histograms.append(newHist)
        # arguments for plotting process below
        modelData = dict(activeParamCount = self.model.activeParamCount(),
                         histograms = histograms)
        pargs = [self.result, self.data]
        pkwargs = dict(axisMargin = axisMargin, outputFilename = outputFilename,
                       modelData = modelData,
                       autoClose = autoClose, logToFile = False)
        if isMac():
            PlotResults(*pargs, **pkwargs)
        else:
            from multiprocessing import Process, Queue
            # multithreaded plotting also logs to file
            pkwargs["logToFile"] = True
            q = Queue() # allow communication between processes
            pkwargs["queue"] = q
            proc = Process(target = PlotResults,
                           args = pargs, kwargs = pkwargs)
            proc.start() # memoize() AssertionsError on Windows raised here
            if not autoClose:
                return # keeps the plot window open
            # wait for the plot window to finish drawing, then close it
            # FIXME: use QTimer for polling in background
            interval = 0.1 # time steps to check if plotting is done
            while q.empty(): # usually needs 1.5sec on my box [ingo]
                time.sleep(interval)
            proc.terminate()

# vim: set ts=4 sts=4 sw=4 tw=0:
