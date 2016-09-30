# -*- coding: utf-8 -*-
# models/scatteringmodel.py

import os.path
import logging
import inspect
from math import sqrt
import numpy as np
from abc import ABCMeta, abstractmethod
from itertools import izip
from numpy import arange, zeros, argmax, hstack
from utils import isList, isNumber, mixedmethod, testfor, classname
from bases.algorithm import AlgorithmBase
from utils.parameter import isActiveParam

class ModelData(object):
    _cumInt = None
    _vset = None
    _wset = None
    _sset = None

    def hdfWrite(self, hdf):
        hdf.writeMembers(self, "cumInt", "vset", "wset", "volumeFraction")

    @property
    def cumInt(self):
        """Returns the cumulated model intensity or signal."""
        return self._cumInt

    @property
    def chisqrInt(self):
        """Make the model intensity comparable to the measured intensity. The
        difference of both will be calculated in BackgroundScalingFit in order
        to perform the chi-square test."""
        return self.cumInt

    @property
    def vset(self):
        """Returns the associated set of volumes."""
        return self._vset

    @property
    def wset(self):
        """Returns the associated set of weights."""
        return self._wset

    @property
    def sset(self):
        """Returns the associated set of surfaces."""
        return self._sset

    def __init__(self, cumInt, vset, wset, sset):
        assert cumInt is not None
        assert vset is not None
        assert wset is not None
        assert sset is not None
        self._cumInt = cumInt.flatten()
        self._vset = vset.flatten()
        self._wset = wset.flatten()
        self._sset = sset.flatten()

    def volumeFraction(self, scaling):
        """Returns the volume fraction based on the provided scaling factor to
        match this model data to the measured data. Assumes that the weights
        'self.wset' contain the scatterer volume squared."""
        return (self.wset * scaling / self.vset).flatten()

class SASModelData(ModelData):
    pass

class DLSModelData(ModelData):

    @property
    def chisqrInt(self):
        """Normalize and square the cumulated model intensities before passing
        them to the chi-square test. This is g1(tau)², the first-order
        correlation function squared."""
        return (self.cumInt / sum(self.wset))**2

    def volumeFraction(self, scaling):
        """Using the square root of the scaling factor to determine the volume
        fraction because the model intensities are squared after cumulation and
        normalization during post-processing."""
        return super(DLSModelData, self).volumeFraction(sqrt(scaling))

class ScatteringModel(AlgorithmBase):
    __metaclass__ = ABCMeta

    @abstractmethod
    def volume(self):
        """Calculates the volume of this model, taking compensationExponent
        into account from input or preset parameters.
        Reimplement this for new models."""
        raise NotImplementedError

    def absVolume(self):
        """Forwarding to usual volume() by default.
        Can be overridden to include SLD."""
        return self.volume()

    def _volume(self, compensationExponent = None):
        """Wrapper around the user-defined function."""
        self.compensationExponent = compensationExponent
        # calling user provided custom model
        return self.absVolume()

    @abstractmethod
    def weight(self):
        """A weighting function for the form factor.
        With SAXS, it is usually the volume squared."""
        raise NotImplementedError

    def _weight(self, compensationExponent = None):
        """Wrapper around the user-defined function."""
        self.compensationExponent = compensationExponent
        return self.weight()

    @abstractmethod
    def formfactor(self, dataset):
        """Calculates the Rayleigh function of this model.
        Reimplement this for new models."""
        raise NotImplementedError

    def _formfactor(self, dataset):
        """Wrapper around the user-defined function."""
        # calling user provided custom model
        i = self.formfactor(dataset)
        return i

    def surface(self):
        """Returns the surface area of a single scatterer. Used for the
        surface weighted distribution histogram. Returns 0 by default.
        Reimplement this for a model."""
        return 0

    def hdfWrite(self, hdf):
        hdf.writeAttributes(modelName = classname(self))
        for p in self.params():
            logging.debug("Writing model parameter: {} to HDF5".format(p.name()))
            hdf.writeMember(self, p.name())

    @abstractmethod
    def calcIntensity(self, data, compensationExponent = None):
        """Calculates the model intensity which is later compared to the data.
        Returns a tuple containing an array of the calculated intensities for
        the grid provided with the data and the volume of a single particle
        based on the model parameters.
        Has to be implemented in derived classes specific to a certain type of
        measurement.
        """
        raise NotImplementedError

    def calc(self, data, pset, compensationExponent = None):
        """Calculates the total intensity and scatterer volume contributions
        using the current model.
        *pset* number columns equals the number of active parameters.
        Returns a ModelData object for a certain type of measurement.
        """
        # remember parameter values
        params = self.activeParams()
        oldValues = [p() for p in params] # this sucks. But we dont want to lose the user provided value
        cumInt = zeros(data.f.binnedData.shape) # cumulated intensities
        vset = zeros(pset.shape[0])
        wset = zeros(pset.shape[0])
        sset = zeros(pset.shape[0])
        # call the model for each parameter set explicitly
        # otherwise the model gets complex for multiple params incl. fitting
        for i in arange(pset.shape[0]): # for each contribution
            for p, v in izip(params, pset[i]): # for each fit param within
                p.setValue(v)
            # result squared or not is model type dependent
            it, vset[i], wset[i], sset[i] = self.calcIntensity(data,
                          compensationExponent = compensationExponent)
            # a set of intensities
            cumInt += it
        # restore previous parameter values
        for p, v in izip(params, oldValues):
            p.setValue(v)
        return self.modelDataType()(cumInt.flatten(), vset, wset, sset)

    @abstractmethod
    def modelDataType(self):
        """Returns the appropriate ModelData class for this type of model.
        """
        raise NotImplementedError

    def generateParameters(self, count = 1):
        """Generates a set of parameters for this model using the predefined
        Parameter.generator. Allows for different random number distributions.
        """
        lst = zeros((count, self.activeParamCount()))
        for idx, param in enumerate(self.activeParams()):
            # generate numbers in different range for each active parameter
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
        return # this does not work anymore
        relerr = getattr(cls, "testRelErr", 1e-5)
        datadir = getattr(cls, "testDataDir", "testdata")
        volumeExponent = getattr(cls, "testVolExp", 1.0)
        filename = os.path.abspath(os.path.join(datadir, filename))
        #dataset = SASData.load(filename)
        if dataset is None:
            return
        model = cls()
        testParams = model.getParametersFromFilename(dataset.filename)
        model.update(testParams)
        # intensity how it's calculated in SASfit
        intensity = (model.vol(None,
                               compensationExponent = volumeExponent)
                     * model._formfactor(dataset, None))**2.
        # computing the relative error to reference data
        delta = abs((dataset.f.binnedData - intensity) / dataset.f.binnedData)
        dmax = argmax(delta)
        testfor(delta.mean() < relerr, AssertionError,
                "Could not verify {model} intensity against\n'{fn}',"
                "\nmean(delta): {mean} >= relErr: {relerr}"
                "\nQ, Int_ref, Int_calc, Delta around dmax({dmax}):"
                "\n{data}"
                .format(model = cls.name(), fn = filename,
                        mean = delta.mean(), relerr = relerr,
                        dmax = dmax, data = hstack((
                            dataset.x0.binnedData.reshape(-1, 1),
                            dataset.f.binnedData.reshape(-1, 1),
                            intensity.reshape(-1, 1),
                            delta.reshape(-1, 1)))[max(0, dmax-4):dmax+5]
                        )
                )

class SASModel(ScatteringModel):
    __metaclass__ = ABCMeta
    canSmear = False # Indicates a model function which supports smearing...

    def modelDataType(self):
        return SASModelData

    def __init__(self):
        # just checking:
        super(SASModel, self).__init__()
        logging.debug("SASData init method called")

    def getQ(self, dataset):
        """ This is a function that returns Q. In case of smearing, dataset itself
        is a 2D matrix of Q-values. When smearing is not enabled, dataset.q contains
        a 1D vector of q.

        I do realize that this is not a good way of doing things. This should be
        replaced at a given point in time by a better solution within sasdata.
        """

        if isinstance(dataset, np.ndarray):
            q = dataset
        else:
            q = dataset.q
        return q

    def weight(self):
        r"""Calculates an intensity weighting used during fitting. It is based
        on the scatterers volume. It can be modified by a user-defined
        compensation exponent *c*. The default value is :math:`c={2 \over 3}`

        :math:`w(r) = v(r)^{2c}`
        """
        return self.volume()**(2 * self.compensationExponent)

    def calcIntensity(self, data, compensationExponent = None):
        r"""Returns the intensity *I*, the volume :math:`v_{abs}` and the
        intensity weights *w* for a single parameter contribution over all *q*:

        :math:`I(q,r) = F^2(q,r) \cdot w(r)`
        """
        v = self._volume(compensationExponent = compensationExponent)
        w = self._weight(compensationExponent = compensationExponent)
        s = self.surface()

        if ((data.config.smearing is not None) and
                self.canSmear and
                data.config.smearing.doSmear() and # serves same purpose as first
                data.config.smearing.inputValid()):
            # inputValid can be removed once more appropriate limits are set in GUI

            # TODO: fix after change from x0Fit to x0:
            locs = data.locs # [data.x0.validIndices] # apply xlimits
            # the ff functions might only accept one-dimensional q arrays
            # kansas = locs.shape
            # locs = locs.reshape((locs.size))
            ff = self._formfactor(locs) # .reshape(kansas)
            qOffset, weightFunc = data.config.smearing.prepared
#            import sys
#            print >>sys.__stderr__, "prepared"
#            print >>sys.__stderr__, unicode(data.config.smearing)
            it = 2 * np.trapz(ff**2 * w * weightFunc,
                    x = qOffset, axis = 1)
        else:
            # calculate their form factors
            ff = self._formfactor(data)
            # a set of intensities
            it = ff**2 * w
        return it, v, w, s

class DLSModel(ScatteringModel):
    __metaclass__ = ABCMeta
    _scatteringVector = None

    def modelDataType(self):
        return DLSModelData

    @property
    def scatteringVector(self):
        return self._scatteringVector

    def calcIntensity(self, data, compensationExponent = None):
        self._scatteringVector = data.scatteringVector
        v = self._volume(compensationExponent = compensationExponent)
        w = self._weight(compensationExponent = compensationExponent)
        # calculate their form factors
        ff = self._formfactor(data)
        # a set of intensities
        it = ff * w
        return it, v, w

# vim: set ts=4 sts=4 sw=4 tw=0:
