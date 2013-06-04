# -*- coding: utf-8 -*-
# McSAS.py
# Find the reST syntax at http://sphinx-doc.org/rest.html

r"""
Overview
========
A class and supplementary functions for Monte-Carlo fitting of SAXS patterns.
It is released under a `Creative Commons CC-BY-SA license
<http://creativecommons.org/licenses/by-sa/3.0/>`_.
Please cite as::

    Brian R. Pauw et al., J. Appl. Cryst. 46, (2013), pp. 365--371
        doi: http://dx.doi.org/10.1107/S0021889813001295

Classes and Functions Defined in This File
------------------------------------------

 - :py:class:`McSAS() <McSAS.McSAS>`:
   A class containing all the Functions required to perform a
   Monte Carlo analysis on small-angle scattering data.
    
Made possible with help from (amongst others)
---------------------------------------------

 - | Samuel Tardif
   | Derivations (mostly observability) and checking of mathematics
 - | Jan Skov Pedersen
   | checking of mathematics
 - | Pawel Kwasniewski <kwasniew@esrf.fr>
   | Code cleanup and documentation
 - | Ingo Bressler <ingo.bressler@bam.de>
   | Code cleanup, modification and documentation

A Note on Units
---------------

Internally, all length units are in meters, all angle units in degrees
clockwise from top. *Intensity* is in
:math:`\left[ 1 \over {m \cdot sr} \right]`,
*q* in :math:`\left[ 1 \over m \right]`.
The electron density contrast squared,
*DeltaRhoSquared* is in :math:`\left[ m^{-4} \right]`.
Other units may be used, but if absolute units are supplied and absolute
volume fractions required, meters are necessitated.

Example Usage
-------------

*For detailed usage, please see the* :doc:`quickstart`

Fitting a single Dataset using all automatic and default Parameters
(may go wrong on poorly conditioned input, needs sensibly-spaced datapoints
and good uncertainty estimates).
The Dataset is considered to consist of three variables *Q*, *I* and *IError*::

 McSAS(Q = Q, I = I, IError = IError, plot = True)

Optional Parameters can be supplied in parameter-value pairs to finetune
optimisation behaviour::

 A = McSAS(Q = Q, I = I, IError = numpy.maximum(0.01 * I, E),
           numContribs = 200, convergenceCriterion = 1,
           contribParamBounds = array([0.5e-9, 35e-9]),
           maxIterations = 1e5, histogramXScale = 'log',
           deltaRhoSquared = 1e30, numReps = 100, plot = True)

Module Documentation
====================
"""

import numpy # For arrays
from numpy import (inf, array, isfinite, reshape, prod, shape, pi, diff, zeros,
                  arange, size, sin, cos, sum, sqrt, log10,
                  isnan, ndim, newaxis)
from scipy import optimize
import os # Miscellaneous operating system interfaces
import time # Timekeeping and timing of objects
import sys # For printing of slightly more advanced messages to stdout
from abc import ABCMeta, abstractmethod
import inspect
import logging
logging.basicConfig(level = logging.INFO)

from dataset import DataSet
from utils import isList
from utils.algorithmbase import AlgorithmBase
from utils.parameter import ParameterFloat, ParameterNumerical, Parameter
from utils.numbergenerator import RandomUniform, RandomExponential

class PropertyNames(object):
    _cache = None

    @classmethod
    def properties(cls):
        """Returns all attributes configured in this class."""
        nameList = dir(cls)
        hashValue = hash(repr(nameList))
        if not cls._cache or cls._cache[0] != hashValue:
            result = [(name, getattr(cls, name)) for name in nameList
                      if not name.startswith("_") and
                      not inspect.ismethod(getattr(cls, name))]
            cls._cache = hashValue, result
        return cls._cache[1]

    @classmethod
    def propNames(cls):
        return [propName for propName, dummy in cls.properties()]

class ScatteringModel(AlgorithmBase, PropertyNames):
    __metaclass__ = ABCMeta
    compensationExponent = 0.5 # default

    def updateParamBounds(self, bounds):
        if not isList(bounds):
            bounds = [bounds,]
        if not isinstance(bounds, list):
            bounds = list(bounds)
        return bounds

    # it doesn't belong to the model?
    # should be instrumentation geometry ...
    def smear(self, arg):
        return arg

    @abstractmethod
    def vol(self, paramValues, compensationExponent = None):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters."""
        if len(self) == 0 and paramValues is None:
            return True
        return paramValues.shape[1] == len(self)

    @abstractmethod
    def ff(self, dataset, paramValues):
        """Calculates the Rayleigh function of this model."""
        if len(self) == 0 and paramValues is None:
            return True
        return paramValues.shape[1] == len(self)

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = numpy.zeros((count, len(self)))
        for idx, param in self:
            # generate numbers in different range for each parameter
            lst[:, idx] = param.generate(count = count)
        # output count-by-nParameters array
        return lst

class Sphere(ScatteringModel):
    """Form factor of a sphere"""
    shortName = "Sphere"
    parameters = (ParameterFloat.make("radius", 1.0,
                    valueRange = (0., numpy.inf),
                    generator = RandomUniform,
                    suffix = "nm", decimals = 1), )
    parameters[0].isActive = True

    def __init__(self):
        ScatteringModel.__init__(self)
        self.radius.setValueRange((1.0, 1e4))

    def updateParamBounds(self, bounds):
        bounds = ScatteringModel.updateParamBounds(self, bounds)
        if len(bounds) < 1:
            return
        if len(bounds) == 1:
            logging.warning("Only one bound provided, "
                            "assuming it denotes the maximum.")
            bounds.insert(0, self.radius.valueRange(0))
        elif len(bounds) > 2:
            bounds = bounds[0:2]
        logging.info("Updating lower and upper contribution parameter bounds "
                     "to: ({0}, {1}).".format(bounds[0], bounds[1]))
        self.radius.setValueRange((min(bounds), max(bounds)))

    def vol(self, paramValues, compensationExponent = None):
        assert ScatteringModel.vol(self, paramValues)
        if compensationExponent is None:
            compensationExponent = self.compensationExponent
        result = (pi*4./3.) * paramValues**(3. * compensationExponent)
        return result

    def ff(self, dataset, paramValues):
        assert ScatteringModel.ff(self, dataset, paramValues)
        r = paramValues.flatten()
        q = dataset[:, 0]
        qr = numpy.outer(q, r)
        result = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.)
        return result

class McSASParameters(PropertyNames):
    model = None
    contribParamBounds = ()
    numContribs = 200
    maxIterations = 1e5
    numReps = 100
    qBounds = None
    psiBounds = None
    priors = () # of shape Rrep, to be used as initial guess for
                # analyse(). It will pass on a Prior to MCFit.
    prior = ()  # of shape Rset, to be used as initial guess for
                # MCFit function
    histogramBins = 50
    histogramXScale = 'log'
    histogramWeighting = 'volume' # can be "volume" or "number"
    deltaRhoSquared = 1.0
    convergenceCriterion = 1.0
    startFromMinimum = False
    maxRetries = 5
    maskNegativeInt = False
    maskZeroInt = False
    lowMemoryFootprint = False
    doPlot = False

class SASData(DataSet):
    """Represents one set of data from a unique source (a file, for example).
    """
    _prepared = None
    _sizeBounds = None

    @classmethod
    def load(cls, filename):
        """Factory method for creating SASData objects from file."""
        if not os.path.isfile(filename):
            logging.warning("File '{0}' does not exist!".format(filename))
            return
        logging.info("Loading '{0}' ...".format(filename))
        sasFile = PDHFile(filename)
        sasData = cls(sasFile.name, sasFile.data)
        return sasData

    @property
    def is2d(self):
        """Returns true if this dataset contains two-dimensional data with
        psi information available."""
        return self.origin.shape[1] > 3 # psi column is present

    @property
    def prepared(self):
        return self._prepared

    @prepared.setter
    def prepared(self, data):
        self._prepared = data

    def __init__(self, title, data):
        DataSet.__init__(self, title, data)

    def setOrigin(self, data):
        DataSet.setOrigin(self, data)
        # determining sizeBounds from q vector
        q = data[:, 0]
        self._sizeBounds = pi / array([q.max(),
                                       min(abs(q.min()),
                                           abs(diff(q).min()))])

    def clip(self, *args, **kwargs):
        self.prepared = self.origin[
                self.clipMask(self.origin, *args, **kwargs)]

    @staticmethod
    def clipMask(data, qBounds = None, psiBounds = None,
                        maskNegativeInt = False, maskZeroInt = False):
        """If q and/or psi limits are supplied in self.Parameters,
        clips the Dataset to within the supplied limits. Copies data to
        :py:attr:`self.FitData` if no limits are set.
        """
        # some shortcut functions, not performance critical as this function
        # is called rarely
        def q(indices):
            return data[indices, 0]
        def intensity(indices):
            return data[indices, 1]
        def psi(indices):
            return data[indices, 3]

        # init indices: index array is more flexible than boolean masks
        validIndices = numpy.where(numpy.isfinite(data[:, 0]))[0]
        def cutIndices(mask):
            validIndices = validIndices[mask]

        # Optional masking of negative intensity
        if maskZeroInt:
            # FIXME: compare with machine precision (EPS)?
            cutIndices(intensity(validIndices) == 0.0)
        if maskNegativeInt:
            cutIndices(intensity(validIndices) > 0.0)
        if isList(qBounds):
            # excluding the lower q limit may prevent q = 0 from appearing
            cutIndices(q(validIndices) > min(qBounds))
            cutIndices(q(validIndices) <= max(qBounds))
        if isList(psiBounds) and data.shape[1] > 3: # psi in last column
            # excluding the lower q limit may prevent q = 0 from appearing
            cutIndices(psi(validIndices) > min(psiBounds))
            cutIndices(psi(validIndices) <= max(psiBounds))

        return validIndices

    @property
    def sizeBounds(self):
        return self._sizeBounds

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
        - *contribParamBounds*: list
            Two-element vector or list indicating upper and lower size
            bounds of the particle radii used in the fitting procedure. If
            not provided, these will be estimated as:
            :math:`R_{max} = {pi \over q_{min}}` and
            :math:`R_{min} = {pi \over q_{max}}`. Units in meter.
        - *numContribs*: int, default: 200
            Number of spheres used for the MC simulation
        - *maxIterations*: int, default: 1e5
            Maximum number of iterations for the :py:func:`MCFit` function
        - *compensationExponent*: float, default: :math:`1.5 \over 3`
            Parameter used to compensate the :math:`volume^2` scaling of each
            sphere contribution to the simulated I(q).
        - *numReps*: int, default: 100
            Number of repetitions of the MC fit for determination of final
            histogram uncertainty.
        - qBounds*: list, default: [0, inf]
            Limits on the fitting range in q.
            Units in :math:`m^{-1}`
        - *histogramBins*: int, default: 50
            Number of bins used for the histogramming procedure.
        - *histogramXScale*: string, default: 'log'
            Can be set to 'log' for histogramming on a logarithmic size scale,
            recommended for q- and/or size-ranges spanning more than a decade.
        - *histogramWeighting*: string, default: 'volume'
            Can be set to 'number' to force plotting of number-weighted
            distributions
        - *deltaRhoSquared*: float, default: 1
            Scattering contrast - when known it will be used to calculate the
            absolute volume fraction of each contribution.
            Units in :math:`m^{-4}`
        - *convergenceCriterion*: float, default: 1
            Convergence criterion for the least-squares fit. The fit converges
            once the :math:`normalized \chi^2 < convergenceCriterion`. If
            convergence is reached with `convergenceCriterion == 1`, the model
            describes the data (on average) to within the uncertainty, and thus
            all information has been extracted from the scattering pattern.
        - *startFromMinimum*: bool, default: False
            If set to False, the starting configuration is a set of spheres
            with radii uniformly sampled between the given or estimated
            bounds. If set to True, the starting configuration is a set of
            spheres with radii set to the lower given or estimated Bound
            (if not zero). Practically, this makes little difference and this
            feature might be depreciated.
        - *maxRetries*: int, default: 5
            If a single MC optimization fails to reach convergence within
            *maxIterations*, it may just be due to bad luck. The procedure
            will try to redo that MC optimization for a maximum of
            *maxRetries* tries before concluding that it is not bad luck
            but bad input.
        - *doPlot*: Bool, default: False
            If set to True, will generate a plot showing the data and fit, as
            well as the Resulting size histogram.
        - *lowMemoryFootprint*: Bool, default: False
            For 2D pattern fitting, or for fitting patterns with a very large
            number of datapoints or contributions, it may make sense to turn
            this option on in order for intensity generating functions not to
            take up much memory. The cost for this is perhaps a 20-ish percent
            reduction in speed.

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
        *volumeHistogramRepetitionsY*: size array (McSASParameters.histogramBins x numReps)
            Volume-weighted particle size distribution bin values for
            each fit repetition (the mean of which is *volumeHistogramYMean*, 
            and the sample standard deviation is *volumeHistogramYStd*).
        *numberHistogramRepetitionsY*: size array (McSASParameters.histogramBins x numReps)
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

    dataset = None # user provided data to work with
    model = None
    result = None
    shortName = "McSAS"
    parameters = (ParameterNumerical.make("numContribs", 200,
                    displayName = "number of contributions",
                    valueRange = (1, 1e6)),
                  ParameterNumerical.make("numReps", 100,
                    displayName = "number of repetitions",
                    valueRange = (1, 1e6)),
                  ParameterNumerical.make("maxIterations", 1e5,
                    displayName = "maximum iterations",
                    valueRange = (1, 1e100)),
                  ParameterNumerical.make("histogramBins", 50,
                    displayName = "number of histogram bins",
                    valueRange = (1, 1e6)),
                  ParameterFloat.make("convergenceCriterion", 1.0,
                    displayName = "convergence criterion",
                    valueRange = (0., numpy.inf)),
    )


    def __init__(self, **kwargs):
        """
        The constructor, takes keyword-value input Parameters. They can be
        one of the aforementioned parameter keyword-value pairs.
        This does the following:

            1. Initialises the variables to the right type
            2. Parses the input
            3. Stores the supplied data twice, one original and one for fitting
                (the latter of which may be clipped or otherwise adjusted)
            4. Applies Q- and optional Psi- limits to the data
            5. Reshapes FitData to 1-by-n dimensions
            6. Sets the function references
            7. Calculates the shape parameter bounds if not supplied
            8. Peforms simple checks on validity of input
            9. Runs the analyse() function which applies the MC fit multiple
               times
            10. Runs the histogram() procedure which processes the MC result
            11. Optionally recalculates the resulting intensity in the same
                shape as the original (for 2D Datasets)
            12. Optionally displays the results graphically.

        .. document private Functions
        .. automethod:: optimScalingAndBackground
        """
        AlgorithmBase.__init__(self)
        # initialize
        self.result = [] # TODO

    def calc(self, **kwargs):
        # set data values
        self.setData(kwargs)
        # set supplied kwargs and passing on
        self.setParameter(kwargs)
        # apply q and psi limits and populate self.FitData
        self.dataset.clip(McSASParameters.qBounds,
                          McSASParameters.psiBounds,
                          McSASParameters.maskNegativeInt,
                          McSASParameters.maskZeroInt)
        if (McSASParameters.model is None or
            not isinstance(McSASParameters.model, ScatteringModel)):
            McSASParameters.model = Sphere() # create instance
        if self.model is None:
            self.model = McSASParameters.model
        logging.info("Using model: {0}".format(str(self.model.name)))
        if not len(self.model):
            logging.warning("No parameters to analyse given! Breaking up.")
            return
        logging.info(
                "\n".join(["Analysing parameters: "]+
                    [str(p) for i, p in self.model])
        )
        self.checkParameters()
        self.analyse()
        self.histogram()

        if self.dataset.is2d:
            # 2D mode, regenerate intensity
            # TODO: test 2D mode
            self.gen2DIntensity()

        if McSASParameters.doPlot:
            self.plot()

    def setData(self, kwargs):
        """Sets the supplied data in the proper location. Optional argument
        *Dataset* can be set to ``fit`` or ``original`` to define which
        Dataset is set. Default is ``original``.
        """
        isOriginal = True
        try:
            kind = kwargs['Dataset'].lower()
            if kind not in ('fit', 'original'):
                raise ValueError
            if kind == 'fit':
                isOriginal = False
        except:
            pass
        # expecting flat arrays, TODO: check this earlier
        data = tuple((kwargs.pop(n, None) for n in ('Q', 'I', 'IError', 'Psi')))
        if data[0] is None:
            raise ValueError("No q values provided!")
        if data[1] is None:
            raise ValueError("No intensity values provided!")
        if data[2] is None:
            # ValueError instead? Is it mandatory to have IError?
            logging.warning("No intensity uncertainties provided!")
            # TODO: generate some
        # data[3]: PSI is optional, only for 2D required
        # TODO: is psi mandatory in 2D? Should ierror be mandatory?
        # can psi be present without ierror?

        # make single array: one row per intensity and its associated values
        # selection of intensity is shorter this way: dataset[validIndices]
        # enforce a certain shape here (removed ReshapeFitdata)
        data = numpy.vstack([d for d in data if d is not None]).T
        if isOriginal:
            self.dataset = SASData("SAS data provided", data)
        else:
            self.dataset = SASData("SAS data provided", None)
            self.dataset.prepared = data
        McSASParameters.contribParamBounds = list(self.dataset.sizeBounds)

    def setParameter(self, kwargs):
        """Sets the supplied Parameters given in keyword-value pairs for known
        setting keywords (unknown key-value pairs are skipped).
        If a supplied parameter is one of the function names, it is stored in
        the self.Functions dict.
        """
        for key in kwargs.keys():
            found = False
            for cls in McSASParameters, ScatteringModel:
                if key in cls.propNames():
                    value = kwargs.pop(key)
                    setattr(cls, key, value)
                    param = getattr(self, key, None)
                    if isinstance(param, Parameter):
                        param.setValue(value)
                    found = True
                    break
            if not found:
                logging.warning("Unknown McSAS parameter specified: '{0}'"
                                .format(key))

    ######################################################################
    ##################### Pre-optimisation Functions #####################
    ######################################################################

    def SetFunction(self, kwargs):
        """Defines Functions. In particular the following are specified:

        - The parameter bounds estimation function *BOUNDS*. Should be able
          to take input argument ContributionParameterBounds to update, i
          should set the parameter bounds in
          ``self.parameter['ContributionParameterBounds']``

        - The random number generation function *RAND* This must take its
          Parameters from self, and have an optional input argument specifying
          the number of sets to return (for MC initialization). It should
          return a set of Nsets-by-nvalues to be used directly in *FF*. This
          may be depreciated soon as is can be generated from within.

        - The Form-factor function *FF*. If called, this should get the
          required information from self and a supplied Nsets-by-nvalues
          shape-specifying parameter array. It should return an Nsets-by-q
          array. Orientational averaging should be part of the form-factor
          function (as it is most efficiently calculated there), so several
          form factor Functions can exist for non-spherical objects.

        - The shape volume calculation function *VOL*, which must be able to
          deal with input argument *PowerCompensationFactor*, ranging from 
          0 to 1. Should accept an Nsets-by-nvalues array returning an Nsets 
          number of (PowerCompensationFactor-compensated)-volumes. 

        - The smearing function *SMEAR*. Should take information from self
          and an input Icalc, to output an Ismear of the same length.

        This function will actually use the supplied function name as function
        reference.
        """
        for kw in kwargs.keys():
            if kw in self.Functions.keys():
                if callable(kwargs[kw]):
                    self.Functions[kw] = kwargs[kw]
                else:
                    # Make it into a function handle/pointer.
                    self.Functions[kw] = getattr(self, kwargs[kw])

    def checkParameters(self):
        """Checks for the Parameters, for example to make sure
        histbins is defined for all, or to check if all Parameters fall
        within their limits.
        For now, all I need is a check that McSASParameters.histogramBins is a 1D vector
        with n values, where n is the number of Parameters specifying
        a shape.
        """
        def fixLength(valueList):
            if not isList(valueList):
                valueList = [valueList for dummy in range(len(self.model))]
            elif len(valueList) > len(self.model):
                del valueList[len(self.model):]
            elif len(valueList) < len(self.model):
                nMissing = len(self.model) - len(valueList)
                valueList.extend([valueList[0] for dummy in range(nMissing)])
            return valueList

        McSASParameters.histogramBins = fixLength(self.histogramBins.value())
        McSASParameters.histogramXScale = fixLength(
                                            McSASParameters.histogramXScale)
        self.model.updateParamBounds(McSASParameters.contribParamBounds)

    def optimScalingAndBackground(self, intObs, intCalc, intError, sc, ver = 2,
            outputIntensity = False, background = True):
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
        if ver == 2:
            """uses scipy.optimize.leastsqr"""
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
            """using scipy.optimize.fmin"""
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

        if outputIntensity:
            return sc, conval, sc[0]*intCalc + sc[1]
        else:
            return sc, conval

    ######################## Shape Functions ######################
    def EllContributionParameterBounds_2D(self):
        """This function will take the q and psi input bounds and outputs
        properly formatted two-element size bounds for ellipsoids. Ellipsoids
        are defined by their equatorial radius, meridional radius and axis
        misalignment (default -45 to 45 degrees in Psi).
        """
        ContributionParameterBounds = \
                self.GetParameter('ContributionParameterBounds')
        q = self.GetData('Q')
        # reasonable, but not necessarily correct, Parameters
        QBounds = array([pi / numpy.max(q),
                         pi / numpy.min((abs(numpy.min(q)),
                                       abs(numpy.min(diff(q)))))])
        if len(ContributionParameterBounds) == 0:
            print "ContributionParameterBounds not provided, so set related "\
                    "to minimum q or minimum q step and maximum q. Lower and "\
                    "upper bounds are {} and {}".format(QBounds[0], QBounds[1])
            ContributionParameterBounds = numpy.array([QBounds[0], QBounds[1],
                                  QBounds[0], QBounds[1],
                                  -45, 45])
        elif len(ContributionParameterBounds) == 6:
            pass
        else:
            print "Wrong number of ContributionParameterBounds provided, "\
                    "defaulting to {} and {} for radii, -45, 45 for "\
                    "misalignment".format(QBounds[0], QBounds[1])
            ContributionParameterBounds = numpy.array([QBounds[0], QBounds[1],
                                  QBounds[0], QBounds[1],
                                  -45, 45])
        ContributionParameterBounds = \
                numpy.array([numpy.min(ContributionParameterBounds[0:2]), 
                    numpy.max(ContributionParameterBounds[0:2]),
                    numpy.min(ContributionParameterBounds[2:4]), 
                    numpy.max(ContributionParameterBounds[2:4]), 
                    numpy.min(ContributionParameterBounds[4:6]), 
                    numpy.max(ContributionParameterBounds[4:6])])

        self.SetParameter(ContributionParameterBounds = 
                ContributionParameterBounds)

    def random_uniform_ell(self, Nell = 1):
        """Random number generator for generating uniformly-sampled
        size- and orientation Parameters for ellipsoids.
        """
        # get Parameters from self
        ContributionParameterBounds = \
                self.GetParameter('ContributionParameterBounds') 
        # generate Nsph random numbers
        Rset = zeros((Nell, 3))
        Rset[:, 0] = numpy.random.uniform(
                numpy.min(ContributionParameterBounds[0]),
                numpy.max(ContributionParameterBounds[1]), Nell)
        Rset[:, 1] = numpy.random.uniform(
                numpy.min(ContributionParameterBounds[2]),
                numpy.max(ContributionParameterBounds[3]), Nell)
        Rset[:, 2] = numpy.random.uniform(
                numpy.min(ContributionParameterBounds[4]),
                numpy.max(ContributionParameterBounds[5]), Nell)
        # output Nsph-by-3 array
        return Rset

    def random_logR_ell(self, Nell = 1):
        """Random number generator which behaves like its uniform counterpart,
        but with a higher likelihood of sampling smaller sizes.
        May speed up some fitting procedures.
        """
        #get Parameters from self
        ContributionParameterBounds = \
                self.GetParameter('ContributionParameterBounds')
        #generate Nsph random numbers
        Rset = zeros((Nell, 3))
        Rset[:, 0] = 10**(numpy.random.uniform(
            log10(numpy.min(ContributionParameterBounds[0])),
            log10(numpy.max(ContributionParameterBounds[1])), Nell))
        Rset[:, 1] = 10**(numpy.random.uniform(
            log10(numpy.min(ContributionParameterBounds[2])), 
            log10(numpy.max(ContributionParameterBounds[3])), Nell))
        Rset[:, 2] = numpy.random.uniform(
                numpy.min(ContributionParameterBounds[4]),
                numpy.max(ContributionParameterBounds[5]), Nell)
        # output Nsph-by-3 array
        return Rset

    def random_logR_oblate_ell(self, Nell = 1):
        """Random number generator which behaves like its uniform counterpart,
        but with a higher likelihood of sampling smaller sizes.
        May speed up some fitting procedures.
        """
        # get Parameters from self
        ContributionParameterBounds = \
                self.GetParameter('ContributionParameterBounds')
        # generate Nsph random numbers
        Rset = zeros((Nell, 3))
        Rset[:, 0] = 10**(numpy.random.uniform(
            log10(numpy.min(ContributionParameterBounds[0])),
            log10(numpy.max(ContributionParameterBounds[1])), Nell))
        for Ni in range(Nell):
            Rset[Ni, 1] = 10**(numpy.random.uniform(
                log10(numpy.min(ContributionParameterBounds[2])),
                log10(numpy.minimum(ContributionParameterBounds[3],
                    Rset[Ni,0])), 1))
        Rset[:,2]=numpy.random.uniform(
                numpy.min(ContributionParameterBounds[4]),
                numpy.max(ContributionParameterBounds[5]), Nell)
        # output Nsph-by-3 array
        return Rset

    def random_logR_prolate_ell(self, Nell = 1):
        """Random number generator which behaves like its uniform counterpart,
        but with a higher likelihood of sampling smaller sizes.
        May speed up some fitting procedures.
        """
        # get Parameters from self
        ContributionParameterBounds = \
                self.GetParameter('ContributionParameterBounds')
        # generate Nsph random numbers
        Rset = zeros((Nell, 3))
        Rset[:, 0] = 10**(numpy.random.uniform(
                            log10(numpy.min(ContributionParameterBounds[0])),
                            log10(numpy.max(ContributionParameterBounds[1])), 
                            Nell))
        for Ni in range(Nell):
            Rset[Ni, 1] = \
                    10**(numpy.random.uniform(
                        log10(numpy.maximum(Rset[Ni, 0], 
                            ContributionParameterBounds[2])), 
                        log10(ContributionParameterBounds[3]), 1))
        Rset[:, 2] = \
                numpy.random.uniform(numpy.min(ContributionParameterBounds[4]),
                        numpy.max(ContributionParameterBounds[5]), Nell)
        # output Nsph-by-3 array
        return Rset

    def vol_ell(self, Rset, PowerCompensationFactor = []):
        """Calculates the volume of an ellipsoid, taking 
        PowerCompensationFactor from input or preset Parameters.
        """
        if PowerCompensationFactor == []:
            PowerCompensationFactor = \
                    self.GetParameter('PowerCompensationFactor')
        return ((4.0/3*pi) * Rset[:, 0]**(2*PowerCompensationFactor) * 
                Rset[:, 1]**(PowerCompensationFactor))[:, newaxis]

    def FF_ell_2D(self, Rset = [], Q = [], Psi = []):
        """Calculates the form factor for oriented ellipsoids,
        normalized to 1 at Q = 0.

        :arg Rset: is n-by-3::

                R1 = Rset[:, 0]
                R2 = Rset[:, 1]
                R3 = Rset[:, 2]

            **R1 < R2**:
                prolate ellipsoid (cigar-shaped)
            **R1 > R2**:
                oblate ellipsoid (disk-shaped)
        
        Rotation is offset from perfect orientation (psi-rot)

        **Note**: All 2D Functions should be able to potentially take
        externally supplied Q and Psi vectors.
        """
        # degrees to radians, forget the dot and get yourself into a
        # non-floating point mess, even though pi is floating point ...
        d_to_r = 1. / 360 * 2 * pi
        if Q == []:
            q = self.GetData('Q')     # 1-by-N
            psi = self.GetData('Psi') # 1-by-N
        else:
            # externally supplied data
            q = Q
            psi = Psi
        R1, R2, rot = Rset[:, 0], Rset[:, 1], Rset[:, 2]
        NR = prod(shape(R1))
        if NR == 1:
            # option 1:
            sda = sin((psi-rot) * d_to_r)
            cda = cos((psi-rot) * d_to_r)
            r = sqrt(R1**2 * sda**2 + R2**2 * cda**2)
            qr = q*r
            Fell = 3*(sin(qr) - qr*cos(qr)) / (qr**3)
            ##quicker? no, 20% slower:
            #Fell=3*(
            #        sin(q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2))
            #        -q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2)
            #        *cos(q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #            +R2**2*cos((psi-rot)*d_to_r)**2)))/
            #               ((q*sqrt(R1**2*sin((psi-rot)*d_to_r)**2
            #                +R2**2*cos((psi-rot)*d_to_r)**2))**3)
        else: # calculate a series
            Fell = zeros([NR, prod(shape(q))])
            for Ri in range(size(R1)):
                sda = sin((psi-rot[Ri]) * d_to_r)
                cda = cos((psi-rot[Ri]) * d_to_r)
                r = sqrt(R1[Ri]**2 * sda**2 + R2[Ri]**2 * cda**2)
                qr = q*r
                Fell[Ri, :] = 3*(sin(qr) - qr*cos(qr)) / (qr**3)

        return Fell # this will be n-by-len(q) array

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
        data = self.dataset.prepared
        # get settings
        priors = McSASParameters.priors
        prior = McSASParameters.prior
        numContribs = self.numContribs.value()
        numReps = self.numReps.value()
        maxRetries = McSASParameters.maxRetries
        # find out how many values a shape is defined by:
        contributions = zeros((numContribs, len(self.model), numReps))
        numIter = zeros(numReps)
        contribIntensity = zeros([1, len(data), numReps])
        start = time.time() # for time estimation and reporting

        # This is the loop that repeats the MC optimization numReps times,
        # after which we can calculate an uncertainty on the Results.
        priorsflag = False
        for nr in range(numReps):
            if (len(prior) <= 0 and len(priors) > 0) or priorsflag:
                # this flag needs to be set as prior will be set after
                # the first pass
                priorsflag = True
                McSASParameters.prior = priors[:, :, nr%size(priors, 2)]
            # keep track of how many failed attempts there have been
            nt = 0
            # do that MC thing! 
            convergence = inf
            while convergence > self.convergenceCriterion.value():
                # retry in the case we were unlucky in reaching
                # convergence within MaximumIterations.
                nt += 1
                (contributions[:, :, nr], contribIntensity[:, :, nr],
                 convergence, details) = self.mcFit(outputIntensity = True,
                                                    outputDetails = True)
                if nt > maxRetries:
                    # this is not a coincidence.
                    # We have now tried maxRetries+2 times
                    logging.warning("Could not reach optimization criterion "
                                    "within {0} attempts, exiting..."
                                    .format(maxRetries+2))
                    return
            # keep track of how many iterations were needed to reach converg.
            numIter[nr] = details['numIterations']

            # in minutes:
            tottime = (time.time() - start)/60. # total elapsed time
            avetime = (tottime / (nr+1)) # average time per MC optimization
            remtime = (avetime*numReps - tottime) # est. remaining time
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
            fitQ = data[:, 0],
            # average number of iterations for all repetitions
            numIter = numIter.mean()))

    def mcFit(self, outputIntensity = False,
                    outputDetails = False, outputIterations = False):
        """
        Object-oriented, shape-flexible core of the Monte Carlo procedure.
        Takes optional arguments:

        *outputIntensity*:
            Returns the fitted intensity besides the Result

        *outputDetails*:
            details of the fitting procedure, number of iterations and so on

        *outputIterations*:
            Returns the Result on every successful iteration step, useful for
            visualising the entire Monte Carlo optimisation procedure for
            presentations.
        """
        data = self.dataset.prepared
        numContribs = self.numContribs.value()
        prior = McSASParameters.prior
        rset = numpy.zeros((numContribs, len(self.model)))
        details = dict()
        # index of sphere to change. We'll sequentially change spheres,
        # which is perfectly random since they are in random order.
        
        q = data[:, 0]
        # generate initial set of spheres
        if size(prior) == 0:
            if McSASParameters.startFromMinimum:
                for idx, param in self.model:
                    mb = min(param.valueRange())
                    if mb == 0: # FIXME: compare with EPS eventually?
                        mb = pi / q.max()
                    rset[:, idx] = numpy.ones(numContribs) * mb * .5
            else:
                rset = self.model.generateParameters(numContribs)
        elif prior.shape[0] != 0: #? and size(numContribs) == 0:
                                  # (didnt understand this part)
            numContribs = prior.shape[0]
            rset = prior
        elif prior.shape[0] == numContribs:
            rset = prior
        elif prior.shape[0] < numContribs:
            print "size of prior is smaller than numContribs. "\
                  "duplicating random prior values"
            randomIndices = numpy.random.randint(prior.shape[0],
                            size = numContribs - prior.shape[0])
            rset = numpy.concatenate((prior, prior[randomIndices, :]))
            print "size now:", rset.shape
        elif prior.shape[0] > numContribs:
            print "Size of prior is larger than numContribs. "\
                  "removing random prior values"
            # remaining choices
            randomIndices = numpy.random.randint(prior.shape[0],
                                                 size = numContribs)
            rset = prior[randomIndices, :]
            print "size now:", rset.shape
        
        vset = self.model.vol(rset)
        if not McSASParameters.lowMemoryFootprint:
            # calculate their form factors
            ffset = self.model.ff(data, rset)
            # calculate the intensities
            iset = ffset**2 * numpy.outer(numpy.ones(ffset.shape[0]), vset**2)
            # the total intensity - eq. (1)
            # intensities for each q in a _row_
            it = iset.sum(axis = 1)
        else:
            it = 0
            for i in numpy.arange(rset.shape[0]):
                # calculate their form factors
                ffset = self.model.ff(data, rset[i].reshape((1, -1)))
                # a set of intensities
                it += ffset**2 * vset[i]**2
        vst = sum(vset**2) # total volume squared

        # Optimize the intensities and calculate convergence criterium
        # SMEAR function goes here
        it = self.model.smear(it)
        intensity = data[:, 1]
        intError = data[:, 2]
        sci = intensity.max() / it.max() # init. guess for the scaling factor
        bgi = intensity.min()
        sc, conval = self.optimScalingAndBackground(
                intensity, it/vst, intError, numpy.array([sci, bgi]), ver = 1)
        # reoptimize with V2, there might be a slight discrepancy in the
        # residual definitions of V1 and V2 which would prevent optimization.
        sc, conval = self.optimScalingAndBackground(
                intensity, it/vst, intError, sc)
        logging.info("Initial Chi-squared value: {0}".format(conval))

        if outputIterations:
            # Output each iteration, starting with number 0. Iterations will
            # be stored in details['paramDistrib'], details['intensityFitted'],
            # details['convergenceValue'], details['scalingFactor'] and
            # details['priorUnaccepted'] listing the unaccepted number of
            # moves before the recorded accepted move.

            # new iterations will (have to) be appended to this, cannot be
            # zero-padded due to size constraints
            details['paramDistrib'] = rset[:, newaxis]
            details['intensityFitted'] = (it/vst*sc[0] + sc[1])[:, newaxis]
            details['convergenceValue'] = conval[newaxis]
            details['scalingFactor'] = sc[:, newaxis]
            details['priorUnaccepted'] = numpy.array(0)[newaxis]

        # start the MC procedure
        intObs = data[:, 1]
        intError = data[:, 2]
        start = time.time()
        numMoves = 0 # tracking the number of moves
        numNotAccepted = 0
        numIter = 0
        ri = 0
        while (conval > self.convergenceCriterion.value() and
               numIter < self.maxIterations.value()):
            rt = self.model.generateParameters()
            ft = self.model.ff(data, rt)
            vtt = self.model.vol(rt)
            itt = (ft**2 * vtt**2).flatten()
            # Calculate new total intensity
            itest = None
            if not McSASParameters.lowMemoryFootprint:
                # we do subtractions and additions, which give us another
                # factor 2 improvement in speed over summation and is much
                # more scalable
                itest = (it - iset[:, ri] + itt)
            else:
                fo = self.model.ff(data, rset[ri].reshape((1, -1)))
                io = (fo**2 * vset[ri]**2).flatten()
                itest = (it.flatten() - io + itt)

            # SMEAR function goes here
            itest = self.model.smear(itest)
            vstest = (sqrt(vst) - vset[ri])**2 + vtt**2
            # optimize intensity and calculate convergence criterium
            # using version two here for a >10 times speed improvement
            sct, convalt = self.optimScalingAndBackground(
                                    intObs, itest/vstest, intError, sc)
            # test if the radius change is an improvement:
            if convalt < conval: # it's better
                # replace current settings with better ones
                rset[ri], sc, conval = rt, sct, convalt
                it, vset[ri], vst = itest, vtt, vstest
                if not McSASParameters.lowMemoryFootprint:
                    iset[:, ri] = itt
                print ("Improvement in iteration number {0}, "
                             "Chi-squared value {1:f} of {2:f}\r"
                             .format(numIter, conval,
                                 self.convergenceCriterion.value())),
                numMoves += 1
                if outputIterations:
                    # output each iteration, starting with number 0. 
                    # Iterations will be stored in details['paramDistrib'],
                    # details['intensityFitted'], details['convergenceValue'],
                    # details['scalingFactor'] and details['priorUnaccepted']
                    # listing the unaccepted number of moves before the
                    # recorded accepted move.

                    # new iterations will (have to) be appended to this,
                    # cannot be zero-padded due to size constraints
                    details['paramDistrib'] = numpy.dstack(
                            (details['paramDistrib'], rset[:, :, newaxis]))
                    details['intensityFitted'] = numpy.hstack(
                            (details['intensityFitted'],
                             (itest/vstest*sct[0] + sct[1]).T))
                    details['convergenceValue'] = numpy.concatenate(
                            (details['convergenceValue'], convalt[newaxis]))
                    details['scalingFactor'] = numpy.hstack(
                            (details['scalingFactor'], sct[:, newaxis]))
                    details['priorUnaccepted'] = numpy.concatenate(
                            (details['priorUnaccepted'],
                             numpy.array((numNotAccepted, ))))
                numNotAccepted = 0
            else:
                # number of non-accepted moves,
                # resets to zero after on accepted move
                numNotAccepted += 1
            # move to next sphere in list, loop if last sphere
            ri = (ri + 1) % (numContribs)
            numIter += 1 # add one to the iteration number

        print # for progress print in the loop
        if numIter >= self.maxIterations.value():
            logging.warning("Exited due to max. number of iterations ({0}) "
                            "reached".format(numIter))
        else:
            logging.info("normal exit")
        # the +0.001 seems necessary to prevent a divide by zero error
        # on some Windows systems.
        elapsed = time.time() - start + 1e-3
        logging.info("Number of iterations per second: {0}".format(
                        numIter/elapsed))
        logging.info("Number of valid moves: {0}".format(numMoves))
        logging.info("Final Chi-squared value: {0}".format(conval))
        details['numIterations'] = numIter
        details['numMoves'] = numMoves
        details['elapsed'] = elapsed

        ifinal = it / sum(vset**2)
        ifinal = self.model.smear(ifinal)
        sc, conval = self.optimScalingAndBackground(
                            intObs, ifinal, intError, sc)

        result = [rset]
        if outputIntensity:
            result.append((ifinal * sc[0] + sc[1]))
        result.append(conval)
        if outputDetails:
            result.append(details)
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
        data = self.dataset.prepared
        q = data[:, 0]
        intensity = data[:, 1]
        intError = data[:, 2]

        # loop over each repetition
        for ri in range(numReps):
            rset = contribs[:, :, ri] # single set of R for this calculation
            # compensated volume for each sphere in the set
            vset = self.model.vol(rset)
            ## TODO: same code than in mcfit pre-loop around line 1225 ff.
            if not McSASParameters.lowMemoryFootprint:
                # Form factors, all normalized to 1 at q=0.
                ffset = self.model.ff(data, rset)
                # Calculate the intensities
                # Intensity for each contribution as used in the MC calculation
                iset = ffset**2 * numpy.outer(
                                        numpy.ones(ffset.shape[0]),
                                        vset**2)
                # total intensity of the scattering pattern
                it = iset.sum(axis = 1)
            else:
                it = 0
                for i in numpy.arange(rset.shape[0]):
                    # calculate their form factors
                    ffset = self.model.ff(data, rset[i].reshape((1, -1)))
                    # a set of intensities
                    it += ffset**2 * vset[i]**2

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
            vpa = self.model.vol(rset, compensationExponent = 1.0) 
            ## TODO: same code than in mcfit pre-loop around line 1225 ff.
            # initial guess for the scaling factor.
            sci = intensity.max() / it.max()
            bgi = intensity.min()
            # optimize scaling and background for this repetition
            sc, conval = self.optimScalingAndBackground(
                    intensity, it, intError, (sci, bgi))
            scalingFactors[:, ri] = sc # scaling and bgnd for this repetition.
            # a set of volume fractions
            volumeFraction[:, ri] = (
                    sc[0] * vsa**2/(vpa * McSASParameters.deltaRhoSquared)
                    ).flatten()
            totalVolumeFraction[ri] = sum(volumeFraction[:, ri])
            numberFraction[:, ri] = volumeFraction[:, ri]/vpa.flatten()
            totalNumberFraction[ri] = sum(numberFraction[:, ri])

            for c in range(numContribs): # for each sphere
                # calculate the observability (the maximum contribution for
                # that sphere to the total scattering pattern)
                # NOTE: no need to compensate for p_c here, we work with
                # volume fraction later which is compensated by default.
                # additionally, we actually do not use this value.
                if not McSASParameters.lowMemoryFootprint:
                    # determine where this maximum observability is
                    # of contribution c (index)
                    qmi = numpy.argmax(iset[:, c]/it)
                    # point where the contribution of c is maximum
                    qm[c, ri] = q[qmi]
                    minReqVol[c, ri] = (
                            intError * volumeFraction[c, ri]
                                    / (sc[0] * iset[:, c])).min()
                else:
                    ffset = self.model.ff(data, rset[c].reshape((1, -1)))
                    ir = (ffset**2 * vset[c]**2).flatten()
                    # determine where this maximum observability is
                    # of contribution c (index)
                    qmi = numpy.argmax(ir.flatten()/it)
                    # point where the contribution of c is maximum
                    qm[c, ri] = q[qmi]
                    minReqVol[c, ri] = (
                            intError * volumeFraction[c, ri]
                                    / (sc[0] * ir)).min() / vpa[c]

                minReqNum[c, ri] = minReqVol[c, ri] / vpa[c]

            numberFraction[:, ri] /= totalNumberFraction[ri]
            minReqNum[:, ri] /= totalNumberFraction[ri]

        # now we histogram over each variable
        # for each variable parameter we define,
        # we need to histogram separately.
        for paramIndex, param in self.model:
            # Now bin whilst keeping track of which contribution ends up in
            # which bin: set bin edge locations
            if McSASParameters.histogramXScale[paramIndex] == 'lin':
                # histogramXLowerEdge contains #histogramBins+1 bin edges,
                # or class limits.
                histogramXLowerEdge = numpy.linspace(
                        min(param.valueRange()),
                        max(param.valueRange()),
                        McSASParameters.histogramBins[paramIndex] + 1)
            else:
                histogramXLowerEdge = 10**numpy.linspace(
                        log10(min(param.valueRange())),
                        log10(max(param.valueRange())),
                        McSASParameters.histogramBins[paramIndex] + 1)

            def initHist(reps = 0):
                """Helper for histogram array initialization"""
                shp = McSASParameters.histogramBins[paramIndex]
                if reps > 0:
                    shp = (McSASParameters.histogramBins[paramIndex], reps)
                return numpy.zeros(shp)

            # total volume fraction contribution in a bin
            volHistRepY = initHist(numReps)
            # total number fraction contribution in a bin
            numHistRepY = initHist(numReps)
            # minimum required number of contributions /in a bin/ to make
            # a measurable impact
            minReqVolBin = initHist(numReps)
            minReqNumBin = initHist(numReps)
            histogramXMean = initHist()
            volHistMinReq = initHist()
            numHistMinReq = initHist()

            for ri in range(numReps):
                # single set of R for this calculation
                rset = contribs[:, paramIndex, ri]
                for bini in range(McSASParameters.histogramBins[paramIndex]):
                    # indexing which contributions fall into the radius bin
                    binMask = (  (rset >= histogramXLowerEdge[bini])
                               * (rset <  histogramXLowerEdge[bini + 1]))
                    # y contains the volume fraction for that radius bin
                    volHistRepY[bini, ri] = sum(volumeFraction[binMask, ri])
                    numHistRepY[bini, ri] = sum(numberFraction[binMask, ri])
                    if not any(binMask):
                        minReqVolBin[bini, ri] = 0
                        minReqNumBin[bini, ri] = 0
                    else:
                        # why? ignored anyway
                        # minReqVolBin[bini, ri] = minReqVol[binMask, ri].max()
                        minReqVolBin[bini, ri] = minReqVol[binMask, ri].mean()
                        # ignored anyway
                        # minReqNumBin[bini, ri] = minReqNum[binMask, ri].max()
                        minReqNumBin[bini, ri] = minReqNum[binMask, ri].mean()
                    if isnan(volHistRepY[bini, ri]):
                        volHistRepY[bini, ri] = 0.
                        numHistRepY[bini, ri] = 0.
            for bini in range(McSASParameters.histogramBins[paramIndex]):
                histogramXMean[bini] = histogramXLowerEdge[bini:bini+2].mean()
                vb = minReqVolBin[bini, :]
                volHistMinReq[bini] = vb[vb < inf].max()
                nb = minReqNumBin[bini, :]
                numHistMinReq[bini] = nb[nb < inf].max()
            volHistYMean = volHistRepY.mean(axis = 1)
            numHistYMean = numHistRepY.mean(axis = 1)
            volHistYStd = volHistRepY.std(axis = 1)
            numHistYStd = numHistRepY.std(axis = 1)

            # store the results, we'll fix this later by a proper structure
            if paramIndex >= len(self.result):
                self.result.append(dict())
            self.result[paramIndex].update(dict(
                histogramXLowerEdge = histogramXLowerEdge,
                histogramXMean = histogramXMean,
                histogramXWidth = diff(histogramXLowerEdge),
                volumeHistogramRepetitionsY = volHistRepY,
                numberHistogramRepetitionsY = numHistRepY,
                volumeHistogramYMean = volHistYMean,
                volumeHistogramYStd = volHistYStd,
                numberHistogramYMean = numHistYMean,
                numberHistogramYStd = numHistYStd,
                volumeHistogramMinimumRequired = volHistMinReq,
                minimumRequiredVolume = minReqVol,
                volumeFraction = volumeFraction,
                totalVolumeFraction = totalVolumeFraction,
                numberHistogramMinimumRequired = numHistMinReq,
                minimumRequiredNumber = minReqNum,
                numberFraction = numberFraction,
                totalNumberFraction = totalNumberFraction,
                scalingFactors = scalingFactors))

    def gen2DIntensity(self):
        """
        This function is optionally run after the histogram procedure for
        anisotropic images, and will calculate the MC fit intensity in
        image form
        """
        print self.result[0].keys()
        contribs = self.result[0]['contribs']
        numContribs, dummy, numReps = contribs.shape

        # load original Dataset
        data = self.dataset.origin
        q = data[:, 0]
        # we need to recalculate the result in two dimensions
        kansas = shape(q) # we will return to this shape
        q = q.flatten()

        print "Recalculating final result, this may take some time"
        # for each Result
        intAvg = zeros(shape(q))
        # TODO: for which parameter?
        scalingFactors = self.result[0]['scalingFactors']
        for ri in range(numReps):
            print 'regenerating set {} of {}'.format(ri, numReps-1)
            rset = contribs[:, :, ri]
            # calculate their form factors
            vset = self.model.vol(rset)
            # calculate the intensities, same code as in mcfit(), see above
            if not McSASParameters.lowMemoryFootprint:
                # Form factors, all normalized to 1 at q=0.
                ffset = self.model.ff(data, rset)
                # Calculate the intensities
                # Intensity for each contribution as used in the MC calculation
                iset = ffset**2 * numpy.outer(numpy.ones(ffset.shape[0]),
                                              vset**2)
                # the total intensity of the scattering pattern
                it = iset.sum(axis = 1)
            else:
                it = 0
                for i in numpy.arange(rset.shape[0]):
                    # calculate their form factors
                    ffset = self.model.ff(data, rset[i].reshape((1, -1)))
                    # a set of intensities
                    it += ffset**2 * vset[i]**2
            vst = sum(vset**2) # total volume squared
            # Optimize the intensities and calculate convergence criterium
            it = self.model.smear(it)
            intAvg = intAvg + it*scalingFactors[0, ri] + scalingFactors[1, ri]
        # print "Initial conval V1", Conval1
        intAvg /= numReps
        # mask (lifted from clipDataset)
        validIndices = SASData.clipMask(data,
                                        McSASParameters.qBounds,
                                        McSASParameters.psiBounds,
                                        McSASParameters.maskNegativeInt,
                                        McSASParameters.maskZeroInt)
        intAvg = intAvg[validIndices]
        # shape back to imageform
        self.result[0]['intensity2d'] = reshape(intAvg, kansas)

    def exportCSV(self, filename, *args, **kwargs):
        """
        This function writes a semicolon-separated csv file to [filename]
        containing an arbitrary number of output variables *\*args*.
        In case of variable length columns, empty fields will contain ''.

        Optional last argument is a keyword-value argument:
        :py:arg:paramIndex = [integer], indicating which shape parameter it
        is intended to draw upon. VariableNumber can also be a list or array
        of integers, of the same length as the number of output variables
        *\*args* in which case each output variable is matched with a shape
        parameter index. Default is zero.

        Input arguments should be names of fields in *self.result*.
        For example::

            McSAS.exportCSV('hist.csv', 'histogramXLowerEdge',
                            'histogramXWidth', 'volumeHistogramYMean',
                            'volumeHistogramYStd', paramIndex = 0)

        I.e. just stick on as many columns as you'd like. They will be
        flattened by default. A header with the result keyword names will be
        added.
        
        Existing files with the same filename will be overwritten by default.
        """
        vna = numpy.zeros(len(args), dtype = int)
        vni = kwargs.pop('paramIndex', None)
        if vni is not None:
            if isList(vni):
                if len(vni) != len(args):
                    print("Error in exportCSV, supplied list of "
                          "variablenumbers does not have the length of 1 or"
                          "the same length as the list of output variables.")
                    return
                for vi in range(len(args)):
                    vna[vi] = vni[vi]
            else:
                # single integer value
                vna += vni
                
        # uses sprintf rather than csv for flexibility
        ncol = len(args)
        # make format string used for every line, don't need this
        # linestr=''
        # for coli in range(ncol):
        #    linestr = linestr+'{'+'};'
        # strip the last semicolon, add a newline
        # linestr = linestr[0:-1]+'\n'

        inlist = list()
        for argi in range(len(args)):
            inlist.append(self.result[vna[argi]][args[argi]].flatten())
        # find out the longest row
        nrow = 0
        for argi in range(len(args)):
            nrow = numpy.max((nrow, len(inlist[argi])))
        # now we can open the file:
        fh = open(filename, 'w')
        emptyfields = 0
        # write header:
        linestr = ''
        for coli in range(ncol):
            linestr = linestr + '{};'.format(args[coli])
        linestr = linestr[0:-1] + '\n'
        fh.write(linestr)
        for rowi in range(nrow):
            linestr = ''
            for coli in range(ncol):
                # print 'rowi {} coli {} len(args[coli]) {}'
                # .format(rowi,coli,len(args[coli]))
                # we ran out of numbers for this arg
                if len(inlist[coli]) <= rowi:
                    linestr = linestr + ';' # add empty field
                    emptyfields += 1
                else:
                    linestr = linestr + '{};'.format(inlist[coli][rowi])
            linestr = linestr[0:-1] + '\n'

            fh.write(linestr)

        fh.close()
        print "{} lines written with {} columns per line, "\
              "and {} empty fields".format(rowi,ncol,emptyfields)

    def plot(self, axisMargin = 0.3):
        """
        This function plots the output of the Monte-Carlo procedure in two
        windows, with the left window the measured signal versus the fitted
        intensity (on double-log scale), and the righthand window the size
        distribution.
        """
        import matplotlib.font_manager as fm
        from matplotlib.pyplot import figure, xticks, yticks, errorbar
        from matplotlib.pyplot import bar, plot, grid, legend, title, xlim
        fontFamilyArial = "Arial", "Bitstream Vera Sans", "sans-serif"
        fontFamilyTimes = "Times", "DejaVu Serif", "serif"
        def setAxis(ah):
            """Sets the axes Parameters. axtyp can be one of 'q' or 'R'"""
            plotfont = fm.FontProperties(
                        # this only works for macs, doesn't it?
                        # family = 'Courier New Bold',
                        # fname = '/Library/Fonts/Courier New Bold.ttf')
                        family = fontFamilyArial)
            textfont = fm.FontProperties(
                        # Baskerville.ttc does not work when saving to eps
                        # family = 'Times New Roman',
                        # fname = '/Library/Fonts/Times New Roman.ttf')
                        family = fontFamilyTimes)
            # setAxis font and ticks
            ah.set_yticklabels(ah.get_yticks(), fontproperties = plotfont,
                               size = 'large')
            ah.set_xticklabels(ah.get_xticks(), fontproperties = plotfont,
                               size = 'large')
            ah.set_xlabel(ah.get_xlabel(), fontproperties = textfont,
                          size = 'x-large')
            ah.set_ylabel(ah.get_ylabel(), fontproperties = textfont,
                          size = 'x-large')
            # qAxis.set_yticklabels(qAxis.get_yticks(),
            #                      fontproperties = plotfont)
            # qAxis.set_xticklabels(qAxis.get_xticks(),
            #                      fontproperties = plotfont)
            # sizeAxis.spines['bottom'].set_color('black')
            ah.spines['bottom'].set_lw(2)
            ah.spines['top'].set_lw(2)
            ah.spines['left'].set_lw(2)
            ah.spines['right'].set_lw(2)
            ah.tick_params(axis = 'both', colors = 'black', width = 2,
                           which = 'major', direction = 'in', length = 6)
            ah.tick_params(axis = 'x', colors = 'black', width = 2,
                           which = 'minor', direction = 'in', length = 3)
            ah.tick_params(axis = 'y', colors = 'black', width = 2,
                           which = 'minor', direction = 'in', length = 3)
            # qAxis.spines['bottom'].set_lw(2)
            # qAxis.spines['top'].set_lw(2)
            # qAxis.spines['left'].set_lw(2)
            # qAxis.spines['right'].set_lw(2)
            # qAxis.tick_params(axis = 'both', colors='black',width=2,
            #                  which='major',direction='in',length=6)
            # qAxis.tick_params(axis = 'x', colors='black',width=2,
            #                  which='minor',direction='in',length=3)
            # qAxis.tick_params(axis = 'y', colors='black',width=2,
            #                  which='minor',direction='in',length=3)
            locs, labels = xticks()
            xticks(locs, map(lambda x: "%g" % x, locs))
            locs, labels = yticks()
            yticks(locs, map(lambda x: "%g" % x, locs))
            return ah

        result = self.result[0]
        # check how many result plots we need to generate: maximum three.
        nhists = len(McSASParameters.histogramXScale)

        # set plot font
        plotfont = fm.FontProperties(
                    size = 'large',
                    family = fontFamilyArial)
        textfont = fm.FontProperties(
                    # Baskerville.ttc does not work when saving to eps
                    size = 'large',
                    family = fontFamilyTimes)
        # initialize figure and axes
        fig = figure(figsize = (7*(nhists+1), 7), dpi = 80,
                     facecolor = 'w', edgecolor = 'k')
        # load original Dataset
        data = self.dataset.origin
        q = data[:, 0]
        intensity = data[:, 1]
        intError = data[:, 2]
        if self.dataset.is2d:
            # 2D data
            psi = data[:, 3]
            # we need to recalculate the result in two dimensions
            # done by gen2DIntensity function
            intensity2d = self.result['intensity2d']
            intShow = intensity.copy()
            # quadrant 1 and 4 are simulated data, 2 and 3 are measured data
            intShow[(psi >   0) * (psi <=  90)] = intensity2d[
                    (psi >   0) * (psi <=  90)]
            intShow[(psi > 180) * (psi <= 270)] = intensity2d[
                    (psi > 180) * (psi <= 270)]
            # xalimits=(-numpy.min(q[:,0]),numpy.max(q[:,-1]))
            # yalimits=(-numpy.min(q[0,:]),numpy.max(q[-1,:]))
            xmidi = int(round(size(q, 1)/2))
            ymidi = int(round(size(q, 0)/2))
            QX = numpy.array([-q[ymidi, 0], q[ymidi, -1]])
            QY = numpy.array([-q[0, xmidi], q[-1, xmidi]])
            extent = (QX[0], QX[1], QY[0], QY[1])

            qAxis = fig.add_subplot(1, (nhists+1), 1, axisbg = (.95, .95, .95),
                                   xlim = QX, ylim = QY, xlabel = 'q_x, 1/m',
                                   ylabel = 'q_y, 1_m')
            imshow(log10(intShow), extent = extent, origin = 'lower')
            qAxis = setAxis(qAxis)
            colorbar()
        else:
            qAxis = fig.add_subplot(
                        1, (nhists+1), 1,
                        axisbg = (.95, .95, .95),
                        xlim = (q.min() * (1 - axisMargin),
                                q.max() * (1 + axisMargin)),
                        ylim = (intensity[intensity != 0].min()
                                    * (1 - axisMargin),
                                intensity.max()
                                    * (1 + axisMargin)
                        ),
                        xscale = 'log', yscale = 'log',
                        xlabel = 'q, 1/m', ylabel = 'intensity, 1/(m sr)'
                    )
            qAxis = setAxis(qAxis)
            errorbar(q, intensity, intError, zorder = 2, fmt = 'k.',
                     ecolor = 'k', elinewidth = 2, capsize = 4, ms = 5,
                     label = 'Measured intensity', lw = 2,
                     solid_capstyle = 'round', solid_joinstyle = 'miter')
            grid(lw = 2, color = 'black', alpha = .5, dashes = [1, 6],
                 dash_capstyle = 'round', zorder = -1)
            aq = numpy.sort(result['fitQ'])
            aI = result['fitIntensityMean'][0, 
                    numpy.argsort(result['fitQ'])]
            plot(aq, aI, 'r-', lw = 3, label = 'MC Fit intensity', zorder = 4)
            plot(aq, numpy.mean(result['scalingFactors'][1, :]) + 0*aq,
                 'g-', linewidth = 3,
                 label = 'MC Background level:\n\t ({0:03.3g})'
                         .format(numpy.mean(result['scalingFactors'][1, :])),
                 zorder = 3)
            leg = legend(loc = 1, fancybox = True, prop = textfont)
        title('Measured vs. Fitted intensity',
              fontproperties = textfont, size = 'x-large')
        sizeAxis = list()
        for histi in range(len(self.result)):
            # get data:
            res = self.result[histi]
            histXLowerEdge = res['histogramXLowerEdge']
            histXMean = res['histogramXMean']
            histXWidth = res['histogramXWidth']
            if McSASParameters.histogramWeighting == 'volume':
                volHistYMean = res['volumeHistogramYMean']
                volHistMinReq = res['volumeHistogramMinimumRequired']
                volHistYStd = res['volumeHistogramYStd']
            elif McSASParameters.histogramWeighting == 'number':
                volHistYMean = res['numberHistogramYMean']
                volHistMinReq = res['numberHistogramMinimumRequired']
                volHistYStd = res['numberHistogramYStd']
            else: 
                print "Incorrect value for histWeighting: "\
                      "should be either 'volume' or 'number'"

            # prep axes
            if McSASParameters.histogramXScale[histi] == 'log':
                # quick fix with the [0] reference. Needs fixing, this
                # plotting function should be rewritten to support multiple
                # variables.
                sizeAxis.append(fig.add_subplot(
                                    1, (nhists + 1), histi + 2,
                                    axisbg = (.95, .95, .95),
                                    xlim = (histXLowerEdge.min()
                                                * (1 - axisMargin),
                                            histXLowerEdge.max()
                                                * (1 + axisMargin)
                                    ),
                                    ylim = (0, volHistYMean.max()
                                                * (1 + axisMargin)
                                    ),
                                    xlabel = 'Radius, m',
                                    ylabel = '[Rel.] Volume Fraction',
                                    xscale = 'log')
                            )
            else:
                sizeAxis.append(fig.add_subplot(
                                    1, (nhists + 1), histi + 2,
                                    axisbg = (.95, .95, .95),
                                    xlim = (histXLowerEdge.min()
                                            - (1 - axisMargin)
                                            * histXLowerEdge.min(),
                                            histXLowerEdge.max()
                                            * (1 + axisMargin)
                                    ),
                                    ylim = (0, volHistYMean.max()
                                                * (1 + axisMargin)
                                    ),
                                    xlabel = 'Radius, m',
                                    ylabel = '[Rel.] Volume Fraction'))

            sizeAxis[histi] = setAxis(sizeAxis[histi])
            # fill axes
            bar(histXLowerEdge[0:-1], volHistYMean, 
                    width = histXWidth, color = 'orange',
                    edgecolor = 'black', linewidth = 1, zorder = 2,
                    label = 'MC size histogram')
            plot(histXMean, volHistMinReq, 'ro', 
                    ms = 5, markeredgecolor = 'r',
                    label = 'Minimum visibility limit', zorder = 3)
            errorbar(histXMean, volHistYMean, volHistYStd,
                    zorder = 4, fmt = 'k.', ecolor = 'k',
                    elinewidth = 2, capsize = 4, ms = 0, lw = 2,
                    solid_capstyle = 'round', solid_joinstyle = 'miter')
            legend(loc = 1, fancybox = True, prop = textfont)
            title('Radius size hist', fontproperties = textfont,
                  size = 'x-large')
            # reapply limits in x
            xlim((histXLowerEdge.min() * (1 - axisMargin),
                  histXLowerEdge.max() * (1 + axisMargin)))

        fig.subplots_adjust(left = 0.1, bottom = 0.11,
                            right = 0.96, top = 0.95,
                            wspace = 0.23, hspace = 0.13)
        # trigger plot window popup
        import pylab
        pylab.show()
        
    def rangeInfo(self, valueRange = [0, inf], paramIndex = 0):
        """Calculates the total volume or number fraction of the MC result
        within a given range, and returns the total numer or volume fraction
        and its standard deviation over all nreps as well as the first four
        distribution moments: mean, variance, skewness and kurtosis
        (Pearson's definition).
        Will use the *histogramWeighting* parameter for determining whether to 
        return the volume or number-weighted values.

        Input arguments are:

            *valueRange*
              The radius range in which the moments are to be calculated
            *paramIndex*
              Which shape parameter the moments are to be calculated for
              (e.g. 0 = width, 1 = length, 2 = orientation)

        Returns a 4-by-2 array, with the values and their sample standard
        deviations over all *numRepetitions*.
        """
        contribs = self.result[0]['contribs']
        numContribs, dummy, numReps = contribs.shape

        volumeFraction = self.result[paramIndex]['volumeFraction']
        numberFraction = self.result[paramIndex]['numberFraction']
        totalVolumeFraction = self.result[paramIndex]['totalVolumeFraction']
        totalNumberFraction = self.result[paramIndex]['totalNumberFraction']
        # Intensity scaling factors for matching to the experimental
        # scattering pattern (Amplitude A and flat background term b,
        # defined in the paper)
        scalingFactors = self.result[paramIndex]['scalingFactors']

        val = numpy.zeros(numReps) # total value
        mu  = numpy.zeros(numReps) # moments..
        var = numpy.zeros(numReps) # moments..
        skw = numpy.zeros(numReps) # moments..
        krt = numpy.zeros(numReps) # moments..

        # loop over each repetition
        for ri in range(numReps):
            # the single set of R for this calculation
            rset = contribs[:, paramIndex, ri]
            validRange = (  (rset > min(valueRange))
                          * (rset < max(valueRange)))
            rset = rset[validRange]
            # compensated volume for each sphere in the set
            vset = volumeFraction[validRange, ri]
            nset = numberFraction[validRange, ri]

            if McSASParameters.histogramWeighting == 'volume':
                val[ri] = sum(vset)
                mu[ri]  = sum(rset * vset)/sum(vset)
                var[ri] = sum( (rset - mu[ri])**2 * vset )/sum(vset)
                sigma   = numpy.sqrt(abs(var[ri]))
                skw[ri] = (  sum( (rset-mu[ri])**3 * vset )
                           / (sum(vset) * sigma**3))
                krt[ri] = ( sum( (rset-mu[ri])**4 * vset )
                           / (sum(vset) * sigma**3))
            elif McSASParameters.histogramWeighting == 'number':
                val[ri] = sum(nset)
                mu[ri]  = sum(rset * nset)/sum(nset)
                var[ri] = sum( (rset-mu[ri])**2 * nset )/sum(nset)
                sigma   = numpy.sqrt(abs(var[ri]))
                skw[ri] = ( sum( (rset-mu[ri])**3 * nset )
                           / (sum(nset) * sigma**3))
                krt[ri] = ( sum( (rset-mu[ri])**4 * nset )
                           / (sum(nset) * sigma**4))
            else:
                logging.error("Moment calculation: "
                              "unrecognised histogramWeighting value!")
                return None

        #now we can calculate the intensity contribution by the subset of
        #spheres highlighted by the range:
        logging.info("Calculating partial intensity contribution of range")
        # loop over each repetition
        partialIntensities = numpy.zeros(
                (numReps, self.result[0]['fitQ'].shape[0]))

        data = self.dataset.prepared
        for ri in range(numReps):
            rset = contribs[:, paramIndex, ri]
            validRange = (  (rset > min(valueRange))
                          * (rset < max(valueRange)))
            rset = rset[validRange][:, newaxis]
            # compensated volume for each sphere in the set
            # TODO: same code as in gen2dintensity and other places
            vset = self.model.vol(rset)
            if not McSASParameters.lowMemoryFootprint:
                # Form factors, all normalized to 1 at q=0.
                ffset = self.model.ff(data, rset)
                # Calculate the intensities
                # Intensity for each contribution as used in the MC calculation
                iset = ffset**2 * numpy.outer(numpy.ones(ffset.shape[0]),
                                              vset**2)
                # the total intensity of the scattering pattern
                it = iset.sum(axis = 1)
            else:
                it = 0
                for i in numpy.arange(rset.shape[0]):
                    # calculate their form factors
                    ffset = self.model.ff(data, rset[i].reshape((1, -1)))
                    # a set of intensities
                    it += ffset**2 * vset[i]**2
            it = self.model.smear(it)
            sc = scalingFactors[:, ri] # scaling and background for this repetition.
            # a set of volume fractions
            partialIntensities[ri, :] = it.flatten() * sc[0]

        rangeInfoResult = dict(
                partialIntensityMean = partialIntensities.mean(axis = 0),
                partialIntensityStd  = partialIntensities.std(axis = 0),
                partialQ =             self.result[0]['fitQ'].flatten(),
                totalValue =           [val.mean(), val.std(ddof = 1)],
                mean =                 [ mu.mean(),  mu.std(ddof = 1)],
                variance =             [var.mean(), var.std(ddof = 1)],
                skew =                 [skw.mean(), skw.std(ddof = 1)],
                kurtosis =             [krt.mean(), krt.std(ddof = 1)]
        )
        self.result.update(rangeInfoResult)
        return rangeInfoResult

# vim: set ts=4 sts=4 sw=4 tw=0:
