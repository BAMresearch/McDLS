# -*- coding: utf-8 -*-
# dataobj/dlsdata.py

"""
Represents data from dynamic light scattering (DLS) measurement.

"""

from __future__ import absolute_import # PEP328

import logging
import numpy as np # For arrays
from utils import classproperty
from utils.units import Length, ScatteringVector, ScatteringIntensity, Angle

from dataobj.dataobj import DataObj

class DLSData(DataObj):
    """Represents one data set.
    """

    @property
    def dataContent(self):
        """Shows the content of the loaded data: DLS?"""
        return ""

    @classproperty
    @classmethod
    def displayDataDescr(cls):
        return () # not used atm

    @classproperty
    @classmethod
    def displayData(cls):
        return ("title", None, "dataContent", None, None, None)

    def __init__(self, **kwargs):
        super(DLSData, self).__init__(**kwargs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
