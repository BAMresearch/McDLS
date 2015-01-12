# -*- coding: utf-8 -*-
# gui/qt.py

from __future__ import absolute_import # PEP328

#QTLIB = "PyQt4"
QTLIB = "PySide"
try:
    __import__(QTLIB)
except ImportError:
    #maybe we have PySide as fallback option
    QTLIB="PySide"
    print("Failed to import '{}'!".format(QTLIB))

import sys
import os
import importlib

thismodule = sys.modules[__name__]

# make qt modules available to be imported from this file/module
def provideQModules():
    for libname in ["QtGui", "QtCore", "QtSvg", "QtXml"]:
        mod = importlib.import_module(".{0}".format(libname), QTLIB)
        sys.modules[libname] = mod
        setattr(thismodule, libname, mod)

try:
    provideQModules()
except ImportError: # works around sphinx docs error
    pass

def pluginDirs():
    libpath = sys.modules[QTLIB].__path__ # PySide.__path__
    for pdir in [os.path.join(p, "plugins") for p in libpath]:
        yield pdir

# vim: set ts=4 sts=4 sw=4 tw=0:
