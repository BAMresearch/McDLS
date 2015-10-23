# -*- coding: utf-8 -*-
# models/scatteringmodel.py

import os.path
import logging
import inspect
from abc import ABCMeta, abstractmethod
from itertools import izip
from numpy import arange, zeros, argmax, hstack
from utils import isList, isNumber, mixedmethod, testfor
from bases.algorithm import AlgorithmBase
from utils.propertynames import PropertyNames
from utils.parameter import isActiveParam
from dataobj import SASData

class ScatteringModel(AlgorithmBase, PropertyNames):
    __metaclass__ = ABCMeta
    # compensationExponent = 1./2 # default, overridden with that from JSON dict

    # it doesn't belong to the model?
    # should be instrumentation geometry ...
    def smear(self, arg):
        return arg

    @abstractmethod
    def volume(self):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters.
        Reimplement this for new models."""
        raise NotImplemented

    def vol(self, compensationExponent = None, useSLD = False):
        self.compensationExponent = compensationExponent
        # calling user provided custom model
        if useSLD and hasattr(self, "absVolume"):
            # FIXME: perhaps, adding a default ScatteringModel.absVolume() forwarding volume()?
            v = self.absVolume()
        else:
            v = self.volume()
        # volume always returns a single value
        assert isNumber(v)
        return v

    @abstractmethod
    def formfactor(self, dataset):
        """Calculates the Rayleigh function of this model.
        Reimplement this for new models."""
        raise NotImplemented

    def ff(self, dataset):
        # calling user provided custom model
        i = self.formfactor(dataset)
        # there has to be one intensity value for each q-vector
        assert i.size == dataset.q.size
        return i

    @abstractmethod
    def calcIntensity(self, data, compensationExponent = None, useSLD = False):
        """Calculates the model intensity which is later compared to the data.
        Returns a tuple containing an array of the calculated intensities for
        the grid provided with the data and the volume of a single particle
        based on the model parameters.
        """
        raise NotImplemented

    def calc(self, data, pset, compensationExponent = None, useSLD = False):
        """Calculates the total intensity and scatterer volume contributions
        using the current model.
        *pset* number columns equals the number of active parameters.
        """
        # remember parameter values
        params = self.activeParams()
        oldValues = [p() for p in params] # this sucks. But we dont want to loose the user provided value
        cumInt = zeros(data.i.shape) # cumulated intensities
        vset = zeros(pset.shape[0])
        # call the model for each parameter set explicitly
        # otherwise the model gets complex for multiple params incl. fitting
        for i in arange(pset.shape[0]): # for each contribution
            for p, v in izip(params, pset[i]): # for each fit param within
                p.setValue(v)
            # result squared or not is model type dependent
            it, vset[i] = self.calcIntensity(data,
                    compensationExponent = compensationExponent,
                    useSLD = useSLD)
            # a set of intensities
            cumInt += it
        # restore previous parameter values
        for p, v in izip(params, oldValues):
            p.setValue(v)
        return cumInt.flatten(), vset

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = zeros((count, self.activeParamCount()))
        for idx, param in enumerate(self.activeParams()):
            # generate numbers in different range for each parameter
            #only for active parameters, otherwise it may try to generate
            #random values for a boolean-type parameter.
            if isActiveParam(param):
                lst[:, idx] = param.generate(count = count)
        # output count-by-nParameters array
        return lst

    def updateParamBounds(self, bounds):
        if not isList(bounds):
            bounds = [bounds,]
        if not isinstance(bounds, list):
            bounds = list(bounds)
        return bounds

    @mixedmethod
    def activeParams(setforcls):
        """returns all "active" parameters of this algorithm"""
        aPars = [par for par in setforcls.params() if isActiveParam(par)]
        return tuple(aPars)

    @mixedmethod
    def activeParamCount(setforcls):
        return len(setforcls.activeParams())

    @mixedmethod
    def activeParamNames(setforcls):
        namelist = list()
        for param in setforcls.activeParams():
            namelist.append(param.displayName())
        return namelist
        
    # helpers for model testing below

    @mixedmethod
    def update(selforcls, paramDict):
        """Update parameter values based on provided dict with parameter
        names as keys."""
        selforcls.fixTestParams(paramDict)
        for key, value in paramDict.iteritems():
            try:
                p = getattr(selforcls, key)
            except: pass
            else:
                p.setValue(p.dtype(value))

    @mixedmethod
    def fixTestParams(selforcls, params):
        """Eventually convert test parameters to the current model.
        Override this method in your model if parameters differ in order or
        meaning from those provided with the test data.
        *params*: A dict with mapping of parameter name to value.
        """
        return params

    @classmethod
    def getParametersFromFilename(cls, filename):
        """Derives model parameters for testing from reference data file."""
        pnames = os.path.splitext(os.path.basename(filename))[0]
        pnames = pnames.split("-")
        errorMsg = ("Could not infer {model} parameters from '{fn}'!"
                    .format(model = cls.name(), fn = filename))
        testfor(len(pnames) > cls.paramCount(), NameError, errorMsg)
        # exclude the base name at front
        pnames = tuple(enumerate(pnames[1:]))
        # store all values by index for reference in fixTestParams()
        result = dict(pnames)
        # map values to parameter names beginning at front
        for i, valueStr in pnames:
            try:
                p = cls.param(i) # raises IndexError eventually
                result[p.name()] = p.dtype(valueStr)
            except IndexError: continue
        return result

    @classmethod
    def test(cls, filename):
        """Regression test of a scattering model. File names are expected
        to contain the parameter values which produce the provided intensity.
        Otherwise implement fixTestParams() for the particular model.

        - *filename*: Name of the file in cls.testDataDir to test against.
        - *cls.testRelErr*: Acceptable mean of relative error against reference
                            intensity. Default: 1e-5
        - *cls.testVolExp*: Volume compensation exponent, sets the amount of
                            volume contribution the intensity is scaled by.
        - *cls.testDataDir*: Directory of test data relative to program dir.
                             Default: "testdata"
        """
        relerr = getattr(cls, "testRelErr", 1e-5)
        datadir = getattr(cls, "testDataDir", "testdata")
        volumeExponent = getattr(cls, "testVolExp", 1.0)
        filename = os.path.abspath(os.path.join(datadir, filename))
        dataset = SASData.load(filename)
        if dataset is None:
            return
        model = cls()
        testParams = model.getParametersFromFilename(dataset.filename)
        model.update(testParams)
        # intensity how it's calculated in SASfit
        intensity = (model.vol(None,
                               compensationExponent = volumeExponent)
                     * model.ff(dataset, None))**2.
        # computing the relative error to reference data
        delta = abs((dataset.i - intensity) / dataset.i)
        dmax = argmax(delta)
        testfor(delta.mean() < relerr, AssertionError,
                "Could not verify {model} intensity against\n'{fn}',"
                "\nmean(delta): {mean} >= relErr: {relerr}"
                "\nQ, Int_ref, Int_calc, Delta around dmax({dmax}):"
                "\n{data}"
                .format(model = cls.name(), fn = filename,
                        mean = delta.mean(), relerr = relerr,
                        dmax = dmax, data = hstack((
                            dataset.q.reshape(-1, 1),
                            dataset.i.reshape(-1, 1),
                            intensity.reshape(-1, 1),
                            delta.reshape(-1, 1)))[max(0, dmax-4):dmax+5]
                        )
                )

class SASModel(ScatteringModel):
    __metaclass__ = ABCMeta

    def slitSmear(self, data, slitShape = "trapezoid", shapeParam = [0., 0.] nSmearSteps = 25):

        def squareSlit(q, shapeParam, n):

            slitWidth = shapeParam[0]
            # prepare integration steps dU:
            dU = np.logspace(np.log10(q.min() / 10.),
                    np.log10(slitWidth / 2.), num = n - 1)
            dU = np.concatenate(([0,], dU)) [np.newaxis, :]
            return dU, np.ones((dU.size,)) / slitWidth

        def trapzSlit(q, shapeParam, n):
            """ defines integration over trapezoidal slit. Top of trapezoid 
            has width xt, bottom of trapezoid has width xb. Note that xb > xt"""
            xt, xb = shapeParam[0], shapeParam[1]

            # ensure things are what they are supposed to be
            assert (xt > 0.)
            if xb < xt:
                xb = xt # should use square profile in this case.

            # prepare integration steps dU:
            dU = np.logspace(np.log10(q.min() / 10.),
                    np.log10(xb / 2.), num = n)
            dU = np.concatenate(([0,], dU)) [np.newaxis, :]

            if xb == xt: 
                y = 1. - (dU * 0.)
            else:
                y = 1. - (dU - xt) / (xb - xt)

            y = np.clip(y, 0., 1.)
            y[dU < xt] = 1.
            Area = (xt + 0.5 * (xb - xt))
            return dU, y / Area
        
        # now we do the actual smearing
        assert isinstance(data.q, numpy.ndarray)
        assert (data.q.ndim == 1)

        # define smearing profile
        if slitShape == "trapezoid":
            slitFcn = trapzSlit
        else:
            slitFcn = squareSlit
        dU, weightFunc = slitFcn(data.q, shapeParam, nSmearSteps)

        # calculate the intensities at sqrt(q**2 + dU **2)
        locs = np.sqrt((q[:,np.newaxis] + 0 * dU)**2 + (0 * q[:, np.newaxis] + dU)**2)
        if fparams is None:
            fVal = fhandle(locs) * (0 * q[:, np.newaxis] + 2 * weightFunc)
        else:
            fVal = fhandle(locs, **fparams) * (0 * q[:, np.newaxis] + 2 * weightFunc)
        sI = np.sqrt(np.trapz(fVal**2, x = dU, axis = 1)) # * 2.
        return sI


    def calcIntensity(self, data, compensationExponent = None, useSLD = False):
        v = self.vol(compensationExponent = compensationExponent,
                     useSLD = useSLD)
        # calculate their form factors
        ff = self.ff(data)
        # a set of intensities
        it = ff**2 * v**2
        return it, v

        # # this section should go in calcIntensity
        # if self.slitWidthTrapzTop() > 0.:
        #     # smearing
        #     return smear(dataset.q, baseCalc, 
        #             fparams = {"radius": self.radius()}, 
        #             slitWidthTrapzTop = self.slitWidthTrapzTop(), 
        #             slitWidthTrapzBottom = self.slitWidthTrapzBottom(), 
        #             nIntSteps = self.smearingSteps())
        # else:
        #     # no smearing
        #     return baseCalc(dataset.q, self.radius())

class DLSModel(ScatteringModel):
    __metaclass__ = ABCMeta

    def calcIntensity(self, data, compensationExponent = None, useSLD = False):
        v = self.vol(compensationExponent = compensationExponent,
                     useSLD = useSLD)
        # calculate their form factors
        ff = self.ff(data)
        # a set of intensities
        it = v**2 * ff
        return it, v

# vim: set ts=4 sts=4 sw=4 tw=0:

        def smear(q, fhandle, fparams = None, slitWidthTrapzTop = None, slitWidthTrapzBottom = None, nIntSteps = 50):
            """ clean smearing program for obtaining slit-smeared scattering patterns. 
            usage: 
            *q*: A one-dimensional scattering vector
            *fhandle*: A function handle to a scattering pattern calculator. Must accept q
            *fparams*: A dictionary of keyword-parameter sets passed on to the calculator
            *slitwidth*: The width of the slit in units of q
            """
        
