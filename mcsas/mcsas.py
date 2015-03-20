# -*- coding: utf-8 -*-
# mcsas/mcsas.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

from __future__ import absolute_import # PEP328
import numpy # For arrays
from numpy import (inf, array, reshape, shape, pi, diff, zeros,
                  size, sum, sqrt, log10,
                  isnan, newaxis)
from scipy import optimize
from itertools import izip
import time # Timekeeping and timing of objects
import copy
import logging
logging.basicConfig(level = logging.INFO)

from utils import isList 
from bases.dataset import DataSet
from bases.algorithm import AlgorithmBase
from utils.parameter import isActiveParam
from utils.tests import isMac
from models.scatteringmodel import ScatteringModel
from models.sphere import Sphere
from gui.utils import processEventLoop

from . import PlotResults
from . import McSASParameters
from sasdata import SASData

class McSAS(AlgorithmBase):
    r"""
    Main class containing all functions required to do Monte Carlo fitting.

    **Required input Parameters:**

        - *Q*: 1D or 2D array of q-values
        - *I*: corresponding intensity values of the same shape
        - *IError*: corresponding intensity uncertainties of the same shape
        - *model*: The scattering model object to assume.
                   It has to be an instance of :py:class:`ScatteringModel`.

    **Optional input Parameters:**

        - *Psi*: 2D array
            Detector angle values, only required for 2D pattern fitting.
        - *compensationExponent*: float, default: :math:`1.5 \over 3`
            Parameter used to compensate the :math:`volume^2` scaling of each
            sphere contribution to the simulated I(q).

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

        *fitIntensityMean*: 1D array (*common result*)
            The fitted intensity, given as the mean of all numReps Results.
        *fitQ*: 1D array (*common result*)
            Corresponding q values
            (may be different than the input q if *QBounds* was used).
        *fitIntensityStd*: array (*common result*)
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

    **Internal Variables**
    
    :py:attr:`self.Dataset`
        Where Q, Psi, I and IError is stored, original Dataset.
    :py:attr:`self.FitData`
        May be populated with a subset of the aforementioned Dataset, limited
        to q-limits or psi limits and to positive I values alone.
    :py:attr:`self.Parameters`
        Where the fitting and binning settings are stored.
    :py:attr:`self.Result`
        Where all the analysis Results are stored. I do not think this needs
        separation after all into Results of analysis and Results of
        interpretation. However, this is a list of dicts, one per variable
        (as the method, especially for 2D analysis, can deal with more than
        one random values. analysis Results are stored along with the
        histogrammed Results of the first variable with index [0]:
    :py:attr:`self.Functions`
        Where the used functions are defined, this is where shape changes,
        smearing, and various forms of integrations are placed.

    """

    # user provided data to work with
    dataOriginal = None
    model = None
    result = None

    # there are several ways to accomplish this depending on where/when
    # McSASParameters() should be called: when creating an instance or
    # for creating different types with different settings ...
    @classmethod
    def factory(cls):
        # calling factory of parent class which is AlgorithmBase currently
        return super(McSAS, cls).factory("McSAS",
                                         *McSASParameters().parameters)

    def calc(self, **kwargs):
        # initialize
        self.result = [] # TODO
        self.stop = False # TODO, move this into some simple result structure
        # # set data values
        # self.setData(kwargs)
        # # set supplied kwargs and passing on
        # self.setParameter(kwargs)
        # apply q and psi limits and populate self.FitData
        if self.dataOriginal is not None:
            self.dataOriginal.qMin = self.qMin()
            self.dataOriginal.qMax = self.qMax()
            self.dataOriginal.pMin = self.psiMin()
            self.dataOriginal.pMax = self.psiMax()
            self.dataOriginal.eMin = self.eMin()
            self.dataOriginal.maskZeroInt = self.maskZeroInt()
            self.dataOriginal.maskNegativeInt = self.maskNegativeInt()
        if (McSASParameters.model is None or
            not isinstance(McSASParameters.model, ScatteringModel)):
            McSASParameters.model = Sphere() # create instance
            logging.info("Default model not provided, setting to: {0}"
                    .format(str(McSASParameters.model.name())))
        if self.model is None:
            self.model = McSASParameters.model
        logging.info("Using model: {0}".format(str(self.model.name())))
        if not self.model.paramCount():
            logging.warning("No parameters to analyse given! Breaking up.")
            return
        logging.info(
                "\n".join([u"Analysing parameters: "] +
                    [unicode(p) + u", active: " + unicode(isActiveParam(p))
                        for p in self.model.params()])
        )
        self.analyse()
        # continue if there are results only
        if not len(self.result):
            return
        self.histogram()

        if self.dataOriginal.is2d:
            # 2D mode, regenerate intensity
            # TODO: test 2D mode
            self.gen2DIntensity()

    ######################################################################
    ##################### Pre-optimisation Functions #####################
    ######################################################################

    def optimScalingAndBackground(self, intObs, intCalc, intError, sc, ver = 2):
        """
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
                  10x faster than version 1
        :arg outputIntensity: *(optional)* Return the scaled intensity as
                              third output argument, default: False
        :arg background: *(optional)* Enables a flat background contribution,
                         default: True

        :returns: (*sc*, *conval*): A tuple of an array containing the
                  intensity scaling factor and background and the reduced
                  chi-squared value.
        """
        def csqr(sc, intObs, intCalc, intError):
            """Least-squares intError for use with scipy.optimize.leastsq"""
            return (intObs - sc[0]*intCalc - sc[1]) / intError
        
        def csqr_nobg(sc, intObs, intCalc, intError):
            """Least-squares intError for use with scipy.optimize.leastsq,
            without background """
            return (intObs - sc[0]*intCalc) / intError

        def csqr_v1(intObs, intCalc, intError):
            """Least-squares for data with known intError,
            size of parameter-space not taken into account."""
            return sum(((intObs - intCalc)/intError)**2) / size(intObs)

        intObs = intObs.flatten()
        intCalc = intCalc.flatten()
        intError = intError.flatten()
        # we need an McSAS instance anyway to call this method
        background = self.findBackground.value()
        if ver == 2:
            # uses scipy.optimize.leastsqr
            if background:
                sc, dummySuccess = optimize.leastsq(
                        csqr, sc, args = (intObs, intCalc, intError),
                        full_output = False)
                conval = csqr_v1(intObs, sc[0]*intCalc + sc[1], intError)
            else:
                sc, dummySuccess = optimize.leastsq(
                        csqr_nobg, sc, args = (intObs, intCalc, intError),
                        full_output = False)
                sc[1] = 0.0
                conval = csqr_v1(intObs, sc[0]*intCalc, intError)
        else:
            # using scipy.optimize.fmin
            # Background can be set to False to just find the scaling factor.
            if background:
                sc = optimize.fmin(
                    lambda sc: csqr_v1(intObs, sc[0]*intCalc + sc[1], intError),
                    sc, full_output = False, disp = 0)
                conval = csqr_v1(intObs, sc[0]*intCalc + sc[1], intError)
            else:
                sc = optimize.fmin(
                    lambda sc: csqr_v1(intObs, sc[0]*intCalc, intError),
                    sc, full_output = False, disp = 0)
                sc[1] = 0.0
                conval = csqr_v1(intObs, sc[0]*intCalc, intError)

        return sc, conval

    def _passthrough(self,In):
        """A passthrough mechanism returning the input unchanged"""
        return In


    ######################################################################
    ####################### optimisation Functions #######################
    ######################################################################

    def analyse(self):
        """This function runs the Monte Carlo optimisation a multitude
        (*numReps*) of times. If convergence is not achieved, it will try 
        again for a maximum of *maxRetries* attempts.
        """
        data = self.dataOriginal
        # get settings
        priors = McSASParameters.priors
        prior = McSASParameters.prior
        numContribs = self.numContribs()
        numReps = self.numReps()
        minConvergence = self.convergenceCriterion()
        if not any([isActiveParam(p) for p in self.model.params()]):
            numContribs, numReps = 1, 1
        # find out how many values a shape is defined by:
        contributions = zeros((self.numContribs(), 
            self.model.activeParamCount(), 
            self.numReps()))
        numIter = zeros(numReps)
        backgrounds = zeros(numReps)
        times = zeros(numReps)
        contribIntensity = zeros([1, len(data.q), numReps])

        priorsflag = False
        # This is the loop that repeats the MC optimization numReps times,
        # after which we can calculate an uncertainty on the Results.
        for nr in range(numReps):
            elapsedStart = time.time() # for tracking elapsed time
            if (len(prior) <= 0 and len(priors) > 0) or priorsflag:
                # this flag needs to be set as prior will be set after
                # the first pass
                priorsflag = True
                McSASParameters.prior = priors[:, :, nr%size(priors, 2)]
            # keep track of how many failed attempts there have been
            nt = 0
            # do that MC thing! 
            convergence = inf
            while convergence > minConvergence:
                # retry in the case we were unlucky in reaching
                # convergence within MaximumIterations.
                nt += 1
                (contributions[:, :, nr], contribIntensity[:, :, nr],
                 convergence, details) = self.mcFit(
                                numContribs, minConvergence,
                                outputIntensity = True, outputDetails = True,
                                nRun = nr)
                if len(contributions) <= 1:
                    break # nothing to fit
                if self.stop:
                    logging.warning("Stop button pressed, exiting...")
                    return
                if nt > self.maxRetries():
                    # this is not a coincidence.
                    # We have now tried maxRetries+2 times
                    logging.warning("Could not reach optimization criterion "
                                    "within {0} attempts, exiting..."
                                    .format(self.maxRetries() + 2))
                    return
            # in minutes:
            # keep track of how many iterations were needed to reach converg.
            numIter[nr] = details.get('numIterations', 0)
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
        self.result.append(dict(
            contribs = contributions, # Rrep
            fitIntensityMean = contribIntensity.mean(axis = 2),
            fitIntensityStd = contribIntensity.std(axis = 2),
            fitQ = data.q,
            # background details:
            background = (backgrounds.mean(), backgrounds.std(ddof = 1)),
            times = times,
            # average number of iterations for all repetitions
            numIter = numIter.mean()))

    def mcFit(self, numContribs, minConvergence,
              outputIntensity = False, outputDetails = False, nRun = None):
        """
        Object-oriented, shape-flexible core of the Monte Carlo procedure.
        Takes optional arguments:

        *outputIntensity*:
            Returns the fitted intensity besides the Result

        *outputDetails*:
            details of the fitting procedure, number of iterations and so on

        *nRun*: "serial number" of run. Used to store results in parameters
            in the right place

        """
        data = self.dataOriginal
        prior = McSASParameters.prior
        rset = numpy.zeros((numContribs, self.model.activeParamCount()))
        compensationExponent = self.compensationExponent()
        details = dict()
        # index of sphere to change. We'll sequentially change spheres,
        # which is perfectly random since they are in random order.
        
        q = data.q 
        # generate initial set of spheres
        if size(prior) == 0:
            if self.startFromMinimum():
                for idx, param in enumerate(self.model.activeParams()):
                    mb = min(param.activeRange())
                    if mb == 0: # FIXME: compare with EPS eventually?
                        mb = pi / (q.max())
                    rset[:, idx] = numpy.ones(numContribs) * mb * .5
            else:
                rset = self.model.generateParameters(numContribs)
        # A prior set can be provided, f.ex. to continue an interrupted run.
        elif (prior.shape[0] != 0) and (numContribs == 0):
            # In this case, rset will assume whatever size prior happens to be
            numContribs = prior.shape[0]
            rset = prior
        elif prior.shape[0] == numContribs:
            # same size prior as requested
            rset = prior
        elif prior.shape[0] < numContribs:
            # prior is smaller in size than requested
            logging.info("size of prior is smaller than numContribs. "
                         "duplicating random prior values")
            randomIndices = numpy.random.randint(prior.shape[0],
                            size = numContribs - prior.shape[0])
            rset = numpy.concatenate((prior, prior[randomIndices, :]))
            logging.info("size now: {}".format(rset.shape))
        elif prior.shape[0] > numContribs:
            # prior is larger in size than requested.
            logging.info("Size of prior is larger than numContribs. "
                         "removing random prior values")
            # remaining choices
            randomIndices = numpy.random.randint(prior.shape[0],
                                                 size = numContribs)
            rset = prior[randomIndices, :]
            logging.info("size now: {}".format(rset.shape))

        # NOTE: put prior into each Parameter, initially
        it, vset = self.calcModel(data, rset, compensationExponent)
        vst = sum(vset**2) # total volume squared

        # Optimize the intensities and calculate convergence criterium
        # SMEAR function goes here
        it = self.model.smear(it)
        # generate initial guess for scaling factor and background
        sci = data.i.max() / it.max() # init. guess for the scaling factor
        bgi = data.i.min()
        sc, conval = self.optimScalingAndBackground(
                data.i, it/vst, data.u, numpy.array([sci, bgi]), ver = 1)
        # reoptimize with V2, there might be a slight discrepancy in the
        # residual definitions of V1 and V2 which would prevent optimization.
        sc, conval = self.optimScalingAndBackground(
                data.i, it/vst, data.u, sc)
        logging.info("Initial Chi-squared value: {0}".format(conval))

        # start the MC procedure
        start = time.time()
        # progress tracking:
        numMoves, numIter, lastUpdate = 0, 0, 0
        # running variable indicating which contribution to change
        ri = 0
        #NOTE: keep track of uncertainties in MC procedure through epsilon
        while (len(vset) > 1 and # see if there is a distribution at all
               conval > minConvergence and
               numIter < self.maxIterations.value() and
               not self.stop):
            rt = self.model.generateParameters()
            # calculate contribution intensity:
            itt, vtt = self.calcModel(data, rt, compensationExponent)
            # Calculate new total intensity, subtract old intensity, add new:
            itest = None
            io, dummy = self.calcModel(data, rset[ri].reshape((1, -1)), 
                    compensationExponent)
            itest = (it.flatten() - io + itt) # is this numerically stable?
            # is numerically stable (so far). Can calculate final uncertainty
            # based on number of valid "moves" and sys.float_info.epsilon

            # SMEAR function goes here
            itest = self.model.smear(itest)
            vstest = (sqrt(vst) - vset[ri])**2 + vtt**2
            # optimize intensity and calculate convergence criterium
            # using version two here for a >10 times speed improvement
            sct, convalt = self.optimScalingAndBackground(
                                    data.i, itest/vstest, data.u, sc)
            # test if the radius change is an improvement:
            if convalt < conval: # it's better
                # replace current settings with better ones
                rset[ri], sc, conval = rt, sct, convalt
                it, vset[ri], vst = itest, vtt, vstest
                logging.info("Improvement in iteration number {0}, "
                             "Chi-squared value {1:f} of {2:f}\r"
                             .format(numIter, conval, minConvergence))
                numMoves += 1

            if time.time() - lastUpdate > 0.25:
                # update twice a sec max -> speedup for fast models
                # because output takes much time especially in GUI
                # process events, check for user input
                # TODO: don't need this:
                # if we calc in separate processes (multiprocessing) the gui 
                # would have its own thread and will not be blocked by calc
                processEventLoop()
                lastUpdate = time.time()
            # move to next sphere in list, loop if last sphere
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
        logging.info("Number of valid moves: {0}".format(numMoves))
        logging.info("Final Chi-squared value: {0}".format(conval))
        details.update({'numIterations': numIter,
            'numMoves': numMoves,
            'elapsed': elapsed})

        ifinal = it / sum(vset**2)
        ifinal = self.model.smear(ifinal)
        sc, conval = self.optimScalingAndBackground(
                            data.i, ifinal, data.u, sc)
        details.update({'scaling': sc[0], 'background': sc[1]})

        result = [rset]
        if outputIntensity:
            result.append((ifinal * sc[0] + sc[1]))
        result.append(conval)
        if outputDetails:
            result.append(details)
        if nRun is not None:
            #store results in parameter
            for idx, param in enumerate(self.model.activeParams()):
                param.setActiveVal(rset[:, idx], index = nRun) 
        # returning <rset, intensity, conval, details>
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
        numContribs, dummy, numReps = contribs.shape

        # volume fraction for each contribution
        volumeFraction = zeros((numContribs, numReps))
        # number fraction for each contribution
        numberFraction = zeros((numContribs, numReps))
        # volume fraction for each contribution
        qm = zeros((numContribs, numReps))
        # volume frac. for each histogram bin
        minReqVol = zeros((numContribs, numReps)) 
        # number frac. for each histogram bin
        minReqNum = zeros((numContribs, numReps))
        totalVolumeFraction = zeros((numReps))
        totalNumberFraction = zeros((numReps))
        # Intensity scaling factors for matching to the experimental
        # scattering pattern (Amplitude A and flat background term b,
        # defined in the paper)
        scalingFactors = zeros((2, numReps))

        # data, store it in result too, enables to postprocess later
        # store the model instance too
        data = self.dataOriginal
        q = data.q
        intensity = data.i
        intError = data.u

        # calc vol/num fraction and scaling factors for each repetition
        for ri in range(numReps):
            rset = contribs[:, :, ri] # single set of R for this calculation
            # compensated volume for each sphere in the set
            it, vset = self.calcModel(data, rset, self.compensationExponent())
            vst = sum(vset**2) # total compensated volume squared 
            it = self.model.smear(it)
            
            # Now for each sphere, calculate its volume fraction
            # (p_c compensated):
            # compensated volume for each sphere in
            # the set Vsa = 4./3*pi*Rset**(3*PowerCompensationFactor)
            # Vsa = VOLfunc(Rset, PowerCompensationFactor)
            vsa = vset # vset did not change here
            # And the real particle volume:
            # compensated volume for each sphere in
            # the set Vsa = 4./3*pi*Rset**(3*PowerCompensationFactor)
            # Vpa = VOLfunc(Rset, PowerCompensationFactor = 1.)
            # dummy, vpa = self.calcModel(data, rset, compensationExponent = 1.0)
            dummy, vpa = self.calcModel(data, rset, 
                    compensationExponent = 1.0, useSLD = True)
            ## TODO: same code than in mcfit pre-loop around line 1225 ff.
            # initial guess for the scaling factor.
            sci = intensity.max() / it.max()
            bgi = intensity.min()
            # optimize scaling and background for this repetition
            sc, conval = self.optimScalingAndBackground(
                    intensity, it, intError, (sci, bgi))
            scalingFactors[:, ri] = sc # scaling and bgnd for this repetition.
            # a set of volume fractions
            # volumeFraction[:, ri] = (
            #         sc[0] * vsa**2/(vpa * self.deltaRhoSquared())
            #         ).flatten()
            volumeFraction[:, ri] = (
                    sc[0] * vsa**2/(vpa)
                    ).flatten()
            totalVolumeFraction[ri] = sum(volumeFraction[:, ri])
            numberFraction[:, ri] = volumeFraction[:, ri]/vpa.flatten()
            totalNumberFraction[ri] = sum(numberFraction[:, ri])

            # calc observability for each sphere/contribution
            for c in range(numContribs):
                # observability: the maximum contribution for
                # that sphere to the total scattering pattern
                # NOTE: no need to compensate for p_c here, we work with
                # volume fraction later which is compensated by default.
                # additionally, we actually do not use this value.
                # again, partial intensities for this size only required
                ir, dummy = self.calcModel(data, rset[c].reshape((1, -1)), 
                        self.compensationExponent())
                # determine where this maximum observability is
                # of contribution c (index)
                qmi = numpy.argmax(ir.flatten()/it.flatten())
                # point where the contribution of c is maximum
                qm[c, ri] = q[qmi]
                minReqVol[c, ri] = (
                        intError * volumeFraction[c, ri]
                                / (sc[0] * ir)).min()
                minReqNum[c, ri] = minReqVol[c, ri] / vpa[c]

            numberFraction[:, ri] /= totalNumberFraction[ri]
            minReqNum[:, ri] /= totalNumberFraction[ri]

        fractions = dict(vol = (volumeFraction, minReqVol),
                         num = (numberFraction, minReqNum))

        # now we histogram over each variable
        # for each variable parameter we define,
        # we need to histogram separately.
        for paramIndex, param in enumerate(self.model.activeParams()):
            param.histograms().calc(contribs, paramIndex, fractions) # new method

    def gen2DIntensity(self):
        """
        This function is optionally run after the histogram procedure for
        anisotropic images, and will calculate the MC fit intensity in
        image form
        """
        contribs = self.result[0]['contribs']
        numContribs, dummy, numReps = contribs.shape

        # load original Dataset
        q = data.q(clipData = False)
        # we need to recalculate the result in two dimensions
        kansas = shape(q) # we will return to this shape
        q = q.flatten()

        logging.info("Recalculating 2D intensity, please wait")
        # for each Result
        intAvg = zeros(shape(q))
        # TODO: for which parameter?
        scalingFactors = self.result[0]['scalingFactors']
        for ri in range(numReps):
            logging.info('regenerating set {} of {}'.format(ri, numReps-1))
            rset = contribs[:, :, ri]
            # calculate their form factors
            it, vset = self.calcModel(data, rset, compensationExponent)
            vst = sum(vset**2) # total volume squared
            # Optimize the intensities and calculate convergence criterium
            it = self.model.smear(it)
            intAvg = intAvg + it*scalingFactors[0, ri] + scalingFactors[1, ri]
        # print "Initial conval V1", Conval1
        intAvg /= numReps
        # mask (lifted from clipDataset)
        validIndices = data.validIndices
        intAvg = intAvg[validIndices]
        # shape back to imageform
        self.result[0]['intensity2d'] = reshape(intAvg, kansas)

    def plot(self, axisMargin = 0.3, parameterIdx = None,
             outputFilename = None):
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
        # get histograms to pass them to plotting
        histograms = [parHist for histograms in
                (p.histograms() for p in self.model.activeParams())
                    for parHist in histograms]
        # modify copies only below,
        # otherwise internal object structure gets inconsistent
        histograms = [copy.copy(h) for h in histograms]
        # remove circular references first for pickling
        for i in range(len(histograms)):
            newParam = histograms[i].param.copy()
            newParam.setActive(False)
            histograms[i].param = newParam
        # arguments for plotting process below
        modelData = dict(activeParamCount = self.model.activeParamCount(),
                         histograms = histograms
                         )
        plotArgs = [self.result, self.dataOriginal, axisMargin,
                    parameterIdx, outputFilename, modelData]
        if isMac():
            plotArgs.append(False) # logToFile, for multithreaded plotting below only
            PlotResults(*plotArgs)
        else:
            from multiprocessing import Process
            plotArgs.append(True) # logToFile
            proc = Process(target = PlotResults, args = plotArgs)
            proc.start()

    # move this to ScatteringModel eventually?
    # which additional output might by useful/required?
    def calcModel(self, data, rset, compensationExponent = None, useSLD = False):
        """Calculates the total intensity and scatterer volume contributions
        using the current model.
        *rset* number columns equals the number of active parameters."""
        # remember parameter values
        params = self.model.activeParams()
        oldValues = [p() for p in params]
        it = 0
        vset = zeros(rset.shape[0])
        # call the model for each parameter value explicitly
        # otherwise the model gets complex for multiple params incl. fitting
        for i in numpy.arange(rset.shape[0]):
            for p, v in izip(params, rset[i]):
                p.setValue(v)
            vset[i] = self.model.vol(compensationExponent, useSLD = useSLD)
            # calculate their form factors
            ffset = self.model.ff(data)
            # a set of intensities
            it += ffset**2 * vset[i]**2
        # restore previous parameter values
        for p, v in izip(params, oldValues):
            p.setValue(v)
        return it.flatten(), vset

# vim: set ts=4 sts=4 sw=4 tw=0:
