# -*- coding: utf-8 -*-
# dataobj/__init__.py

from __future__ import absolute_import # PEP328

__all__ = ["DataObj", "SASData", "DLSData"]

from dataobj.dataobj import DataObj
from dataobj.sasdata import SASData
from dataobj.dlsdata import DLSData

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
