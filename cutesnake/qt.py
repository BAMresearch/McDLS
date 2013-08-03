# -*- coding: utf-8 -*-
# qt.py

__all__ = ["QtGui", "QtCore", "QtSvg"]
PYQT = "PyQt4" #what if I only have PySide? is there a nice way for checking?
try:
    __import__(PYQT)
except ImportError:
    #maybe we have PySide as fallback option
    PYQT="PySide"

import sys
import importlib

thismodule = sys.modules[__name__]

# make qt modules available to be imported from this file/module
for libname in __all__:
    mod = importlib.import_module(".{0}".format(libname), PYQT)
    sys.modules[libname] = mod
    setattr(thismodule, libname, mod)

# vim: set ts=4 sts=4 sw=4 tw=0:
