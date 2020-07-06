# -*- coding: utf-8 -*-
# dataobj/__init__.py


__all__ = ["DataObj", "SASData", "DLSData", "DataConfig"]

from dataobj.dataobj import DataObj
from dataobj.datavector import DataVector
from dataobj.dataconfig import DataConfig
from dataobj.sasconfig import SASConfig
from dataobj.sasdata import SASData
from dataobj.dlsdata import DLSData

if __name__ == "__main__":
    import doctest
    doctest.testmod()

# vim: set ts=4 sts=4 sw=4 tw=0:
