# -*- coding: utf-8 -*-
# scattering/dataobj.py

"""
Represents input data associated with a measurement.
"""

from __future__ import absolute_import # PEP328
import os # Miscellaneous operating system interfaces
from numpy import all as np_all
import numpy as np

# related to the class below
from abc import ABCMeta, abstractproperty, abstractmethod
from bases.dataset import DataSet, DisplayMixin
from utils.units import Unit, NoUnit

import sys

class DataVector(object):
    """ a class for combining aspects of a particular vector of data.
    This is intended only as a storage container without additional functionality.
    """
    _name = None # descriptor for axes and labels, unicode string at init.
    _raw = None # 
    _unit = None # instance of unit
    _limit = None # two-element vector with min-max
    _validIndices = None # valid indices. 
    _editable = False # whether raw can be written or not
    
    def __init__(self, name, raw, unit = None, limit = None, editable = False):
        self._name = name
        self._raw = raw
        self.limit = limit
        self.unit = unit
        self.editable = editable

    @property
    def name(self):
        return unicode(self._name)

    @property
    def validIndices(self):
        if self._validIndices is None:
            return range(self.raw.size)
        else:
            return self._validIndices
    @validIndices.setter
    def validIndices(self, value):
        # assert (value.max() <= self.raw.size)
        self._validIndices = value

    @property
    def value(self):
        return self.origin.copy()[self.validIndices]
    @value.setter
    def value(self, val):
        assert(val.size == np.size(self.validIndices))
        self.unit.toDisplay(val)
        self.raw = val
        self.raw[self.validIndices] = self.unit.toDisplay(val)

    @property
    def origin(self):
        return self.unit.toSi(self.raw)

    @property
    def raw(self):
        return self._raw
    @raw.setter
    def raw(self, value):
        assert(self.editable)
        assert(value.size == self.raw.size)
        self._raw = value
    
    @property
    def editable(self):
        return self._editable
    @editable.setter
    def editable(self, value):
        assert(isinstance(value, bool))
        self._editable = value

    @property
    def unit(self):
        return self._unit
    @unit.setter
    def unit(self, value):
        if value is None:
            self._unit = NoUnit
        else:
            assert(isinstance(value, Unit))
            self._unit = value

    @property
    def limit(self):
        return self._limit
    @limit.setter
    def limit(self, value):
        print('Limit value: {}'.format(value))
        if value is None:
            self._limit = [self.raw.min(), self.raw.max()]
        else:
            self._limit = [np.maximum(np.min(value), self.origin.min()), 
                    np.minimum(np.max(value), self.origin.max())]

    @property
    def limsString(self):
        return u"{0:.3g} ≤ {valName} ({magnitudeName}) ≤ {1:.3g}".format(
                self.unit.toDisplay(self.limit[0]),
                self.unit.toDisplay(self.limit[1]),
                magnitudeName = self.unit.displayMagnitudeName,
                valName = self.name)

# formerly known as 'ScatteringData', better? also for the module?
class DataObj(DataSet, DisplayMixin):
    """General container for data loaded from file. It offers specialised
    methods to derive information from the provided data.
    """
    __metaclass__ = ABCMeta
    _filename = None
    _config = None
    _x0 = None
    _x1 = None
    _x2= None
    _f = None
    _fu = None

    # These are to be set by the particular application dataset: 
    # i.e.: x = q, y = psi, f = I for SAS, x = tau, f = g1 for DLS
    @abstractproperty
    def x0(self): # sampling vector 1
        return self._x0

    # @abstractproperty # sampling vector 2
    # def x1(self):
    #     return self._x1

    # @abstractproperty # sampling vector 3
    # def x2(self):
    #     return self._x2

    @abstractproperty # measurement vector
    def f(self):
        return self._f

    @abstractproperty # measurement vector uncertainty
    def fu(self):
        return self._fu

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
        if config is None:
            return False
        if self.config is None:
            self._config = config
            return True
        if self.config == config:
            return False
        self._config = config.copy()
        return True

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
