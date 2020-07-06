# -*- coding: utf-8 -*-
# bases/model/modeldata.py

class ModelData(object):
    _cumInt = None
    _vset = None
    _wset = None
    _sset = None
    _numParams = 0

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

    @property
    def numParams(self):
        """Returns the number of active (fitted) parameters."""
        return self._numParams

    def __init__(self, cumInt, vset, wset, sset, numParams):
        assert cumInt is not None
        assert vset is not None
        assert wset is not None
        assert sset is not None
        self._cumInt = cumInt.flatten()
        self._vset = vset.flatten()
        self._wset = wset.flatten()
        self._sset = sset.flatten()
        self._numParams = abs(numParams)

    def volumeFraction(self, scaling):
        """Returns the volume fraction based on the provided scaling factor to
        match this model data to the measured data. Assumes that the weights
        'self.wset' contain the scatterer volume squared."""
        return (self.wset * scaling / self.vset).flatten()

class SASModelData(ModelData):
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
