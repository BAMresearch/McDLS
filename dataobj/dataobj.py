# -*- coding: utf-8 -*-
# scattering/dataobj.py

"""
Represents input data associated with a measurement.
"""

from __future__ import absolute_import # PEP328
import os # Miscellaneous operating system interfaces
from numpy import all as np_all
import numpy as np

from abc import ABCMeta, abstractproperty, abstractmethod
from bases.dataset import DataSet, DisplayMixin
from utils.units import Unit, NoUnit
from utils import classproperty

class DataVector(object):
    """ a class for combining aspects of a particular vector of data.
    This is intended only as a storage container without additional functionality.
    """
    _name = None # descriptor for axes and labels, unicode string at init.
    _raw = None # Relevant data directly from input file, non-SI units, unsanitized
    _siData = None # copy raw data in si units, often accessed
    _unit = None # instance of unit
    _limit = None # two-element vector with min-max
    _validIndices = None # valid indices. 
    _editable = False # whether raw can be written or not
    
    def __init__(self, name, raw, unit = None, limit = None, editable = False):
        self._name = name
        self._raw = raw
        self.unit = unit
        self.limit = limit
        assert(isinstance(editable, bool))
        self._editable = editable

    @property
    def name(self):
        return unicode(self._name)

    @property
    def validIndices(self):
        if self._validIndices is None:
            return np.arange(self.raw.size)
        else:
            return self._validIndices

    @validIndices.setter
    def validIndices(self, value):
        # assert (value.max() <= self.raw.size)
        self._validIndices = value

    @property
    def sanitized(self):
        return self.siData.copy()[self.validIndices]

    @sanitized.setter
    def sanitized(self, val):
        assert(val.size == self.validIndices.size)
        self.siData[self.validIndices] = val

    @property
    def siData(self):
        return self._siData

    @siData.setter
    def siData(self, vec):
        self._siData = vec

    @property
    def raw(self):
        return self._raw

    @property
    def editable(self):
        return self._editable

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, newUnit):
        if not isinstance(newUnit, Unit):
            self._unit = NoUnit()
            self._siData = self.raw.copy()
        else:
            self._unit = newUnit
            self._siData = self.unit.toSi(self.raw)

    @property
    def limit(self):
        return self._limit

    @limit.setter
    def limit(self, value):
        self.setLimit(value)

    def setLimit(self, newLimit): # override&referencable
#        print('Limit value: {}'.format(value))
        if newLimit is None:
            self._limit = [self.siData.min(), self.siData.max()]
        else:
            self._limit = [np.maximum(np.min(newLimit), self.siData.min()),
                           np.minimum(np.max(newLimit), self.siData.max())]

    @property
    def limsString(self):
        return u"{0:.3g} ≤ {valName} ({magnitudeName}) ≤ {1:.3g}".format(
                self.unit.toDisplay(self.limit[0]),
                self.unit.toDisplay(self.limit[1]),
                magnitudeName = self.unit.displayMagnitudeName,
                valName = self.name)

class DataObj(DataSet, DisplayMixin):
    """General container for data loaded from file. It offers specialised
    methods to derive information from the provided data.
    """
    __metaclass__ = ABCMeta
    _filename = None
    _config = None

    # interface for basic DataVectors

    # These are to be set by the particular application dataset: 
    # i.e.: x = q, y = psi, f = I for SAS, x = tau, f = (G1 - 1) for DLS
    @property
    def x0(self):
        """First sampling vector."""
        raise NotImplementedError

    @property
    def x1(self):
        """Second sampling vector."""
        raise NotImplementedError

    @property
    def x2(self):
        """Third sampling vector."""
        raise NotImplementedError

    @property
    def f(self):
        """The measurement vector."""
        raise NotImplementedError

    @property
    def fu(self):
        """The measurement uncertainty regarding the measurement vector *f*.
        """
        raise NotImplementedError

    # other common meta data

    @classproperty
    @classmethod
    def sourceName(cls):
        """Returns the name of the measurement method."""
        raise NotImplementedError

    @property
    def filename(self):
        return self._filename

    def setFilename(self, fn):
        """Stores the absolute path to this data file.
        Should be reviewed when data sets can be created from several files."""
        if fn is None or not os.path.isfile(fn):
            return
        self._filename = os.path.abspath(fn)

    @abstractproperty
    def count(self):
        raise NotImplementedError

    @property
    def is2d(self):
        return False

    def accumulate(self, others):
        return None

    @property
    def config(self):
        return self._config

    def setConfig(self, config = None):
        """Set the extended configuration for this data and returns true if
        the configuration was different and an update was necessary."""
        if not isinstance(config, self.configType):
            return False
        # always replacing the config if it's valid, it originates from here
        self._config = config
        self.updateConfigMeta()
        return True

    def updateConfigMeta(self):
        """Updates general meta data of the config object
        based on this data set."""
        descr = self.config.xLow.displayName().format(x = self.x0.name)
        self.config.xLow.setDisplayName(descr)
        descr = self.config.xHigh.displayName().format(x = self.x0.name)
        self.config.xHigh.setDisplayName(descr)
        descr = self.config.fMaskZero.displayName().format(f = self.f.name)
        self.config.fMaskZero.setDisplayName(descr)
        descr = self.config.fMaskNeg.displayName().format(f = self.f.name)
        self.config.fMaskNeg.setDisplayName(descr)

    @abstractproperty
    def configType(self):
        """Returns a compatible DataConfig type."""
        raise NotImplementedError

    @abstractproperty
    def modelType(self):
        """Returns a compatible ScatteringModel type."""
        raise NotImplementedError

    def __init__(self, **kwargs):
        super(DataObj, self).__init__(**kwargs)

    def __eq__(self, other):
        return (np_all(self.rawArray == other.rawArray)
                and self.title == other.title
                and self.filename == other.filename)

    def __neq__(self, other):
        return not self.__eq__(other)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
