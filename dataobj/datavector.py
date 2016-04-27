# -*- coding: utf-8 -*-
# scattering/datavector.py

"""
A class describing a vector with limits, units, mask and uncertainties
"""

from __future__ import absolute_import # PEP328
import numpy as np
import h5py

from utils.units import Unit, NoUnit

def h5w(wloc, field, hDat, hType = "dataset"):
    """ 
        writes dataset *hDat* to HDF5 location *wloc*, deleting if exists 
        htype can be "dataset" or "attribute"
    """
    # remove old field, only removes link, does not reclaim!
    # ideally, new h5py file should be generated on end:
    # http://stackoverflow.com/questions/11194927/deleting-information-from-an-hdf5-file
    if field in wloc:
        del wloc[field]
    if "dataset" in hType:
        wloc.create_dataset(field, data = hDat, compression = "gzip")
    else:
        wloc[field] = hDat

class DataVector(object):
    """ a class for combining aspects of a particular vector of data.
    This is intended only as a storage container without additional functionality.
    """
    _name = None # descriptor for axes and labels, unicode string at init.
    _raw = None # Relevant data directly from input file, non-SI units, unsanitized
    _rawU = None # raw uncertainties, unsanitized
    _siData = None # copy raw data in si units, often accessed
    _siDataU = None # uncertainties in si units, sanitized, with eMin taken into account
    _unit = None # instance of unit
    _limit = None # two-element vector with min-max
    _validIndices = None # valid indices. 
    # specify which values are to be stored in a HDF5 file. 
    _h5Fields = ["raw", "rawU", "siData", "siDataU", "validIndices", "limit"]
    _h5Callers = ["unit"] # writeHDF will be called. 
    
    def writeHDF(self, filename, loc):
        """ 
        Writes the vector to an HDF5 output file *filename*, at location *loc*. 
        This location should be e.g. "/mcentry01/[ sas | dls ]data01/". The 
        name of the vector (self.name) is appended to this path in this method.
        """
        with h5py.File(filename) as h5f:
            wloc = h5f.require_group(loc + self.name) # unicode's no problem
            for field in self._h5Fields:
                hDat = getattr(self, field, None)
                if hDat is not None:
                    h5w(wloc, field, hDat, hType = "dataset")

                    # the require_dataset solution does not work for lists
            # write unit:
            h5w(wloc, "unit", self.unit.name(), hType = "attribute")
            # for attribute in self._h5Attrs:
            #     hAttr = getattr(self, attribute, None)
            #     if hAttr is not None:
            #         wloc[attribute] = hAttr

    def __init__(self, name, raw, rawU = None, unit = None):
        self._name = name
        self._raw = raw
        self._rawU = rawU
        self.unit = unit
        self.validIndices = np.arange(self.raw.size) # sets limits as well

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

    # raw data values
    @property
    def raw(self):
        return self._raw

    # raw uncertainties
    @property
    def rawU(self):
        return self._rawU

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, newUnit):
        if not isinstance(newUnit, Unit):
            self._unit = NoUnit()
            self.siData = self.raw.copy()
            if self.rawU is not None:
                self.siDataU = self.rawU.copy()
        else:
            self._unit = newUnit
            self.siData = self.unit.toSi(self.raw)
            if self.rawU is not None:
                self.siDataU = self.unit.toSi(self.rawU)

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
