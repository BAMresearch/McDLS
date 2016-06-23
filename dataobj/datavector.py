# -*- coding: utf-8 -*-
# dataobj/datavector.py

"""
A class describing a vector with limits, units, mask and uncertainties
"""

from __future__ import absolute_import # PEP328
import numpy as np
import h5py

from utils.units import Unit, NoUnit
from utils.hdf5base import h5w, HDF5Mixin

class DataVector(HDF5Mixin):
    """ a class for combining aspects of a particular vector of data.
    This is intended only as a storage container without additional functionality.
    """
    _name = None # descriptor for axes and labels, unicode string at init.
    _rawData = None # Relevant data directly from input file, non-SI units, unsanitized
    _rawDataU = None # raw uncertainties, unsanitized
    _siData = None # copy raw data in si units, sanitized, often accessed
    _siDataU = None # uncertainties in si units, sanitized, with eMin taken into account
    _binnedData = None # binned, sanitized data
    _binnedDataU = None # binned, sanitized data
    _unit = None # instance of unit
    _limit = None # two-element vector with min-max
    _validIndices = None # valid indices.
    # specify which values are to be stored in a HDF5 file.
    _h5Datasets = ["rawData", "rawDataU",
                   "siData", "siDataU",
                   "binnedData", "binnedDataU",
                   "validIndices", "limit"]
    _h5Callers = ["unit"] # writeHDF will be called.
    _h5Attrs = ["name"]

    def hdfWrite(self, hdf):
        hdf.writeMembers(self, "rawData", "rawDataU", "siData", "siDataU",
                "binnedData", "binnedDataU", "validIndices", "limit", "unit")

    def __init__(self, name, raw, rawU = None, unit = None):
        self._name = name
        self._rawData = raw
        self._rawDataU = rawU
        self.unit = unit
        self.validIndices = np.arange(self.rawData.size) # sets limits as well

    @property
    def name(self):
        return unicode(self._name)

    @property
    def validIndices(self):
        return self._validIndices

    @validIndices.setter
    def validIndices(self, indices):
        assert indices.min() >= 0
        assert indices.max() <= self.siData.size
        self._validIndices = indices
        self._limit = [self.sanitized.min(), self.sanitized.max()]

    @property
    def sanitized(self):
        return self.siData.copy()[self.validIndices]

    @sanitized.setter
    def sanitized(self, val):
        assert(val.size == self.validIndices.size)
        self.siData[self.validIndices] = val

    @property
    def sanitizedU(self):
        if self.siDataU is None:
            return None
        return self.siDataU.copy()[self.validIndices]

    @sanitizedU.setter
    def sanitizedU(self, val):
        assert(val.size == self.validIndices.size)
        self.siDataU[self.validIndices] = val

    # siData
    @property
    def siData(self):
        return self._siData

    @siData.setter
    def siData(self, vec):
        self._siData = vec

    # siDataU, uncertainties on siData
    @property
    def siDataU(self):
        return self._siDataU

    @siDataU.setter
    def siDataU(self, vec):
        self._siDataU = vec

    # binnedDataU, uncertainties on siData
    @property
    def binnedData(self):
        if self._binnedData is not None:
            return self._binnedData
        else:
            return self.sanitized

    @binnedData.setter
    def binnedData(self, vec):
        self._binnedData = vec

    # binnedDataU, uncertainties on siData
    @property
    def binnedDataU(self):
        if self._binnedDataU is not None:
            return self._binnedDataU
        else:
            return self.sanitizedU

    @binnedDataU.setter
    def binnedDataU(self, vec):
        self._binnedDataU = vec

    # raw data values
    @property
    def rawData(self):
        return self._rawData

    # raw uncertainties
    @property
    def rawDataU(self):
        return self._rawDataU

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, newUnit):
        if not isinstance(newUnit, Unit):
            self._unit = NoUnit()
            self.siData = self.rawData.copy()
            if self.rawDataU is not None:
                self.siDataU = self.rawDataU.copy()
        else:
            self._unit = newUnit
            self.siData = self.unit.toSi(self.rawData)
            if self.rawDataU is not None:
                self.siDataU = self.unit.toSi(self.rawDataU)

    # TODO: define min/max properties for convenience?
    @property
    def limit(self):
        return self._limit

    @property
    def limsString(self):
        return u"{0:.3g} ≤ {valName} ({magnitudeName}) ≤ {1:.3g}".format(
                self.unit.toDisplay(self.limit[0]),
                self.unit.toDisplay(self.limit[1]),
                magnitudeName = self.unit.displayMagnitudeName,
                valName = self.name)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
