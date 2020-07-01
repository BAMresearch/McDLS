# -*- coding: utf-8 -*-
# gui/qt.py

from __future__ import absolute_import # PEP328

__all__ = []

import sys
import os
import importlib

import PySide2
from PySide2 import QtWidgets, QtGui, QtCore, QtSvg, QtXml
from PySide2 import __path__ as PySidePath

try:
    sys.modules["QtWidgets"] = QtWidgets
    sys.modules["QtGui"] = QtGui
    sys.modules["QtCore"] = QtCore
    sys.modules["QtSvg"] = QtSvg
    sys.modules["QtXml"] = QtXml
except ImportError: # works around sphinx docs error
    raise
    pass

def pluginDirs():
    libpath = PySidePath
    for pdir in [os.path.join(p, "plugins") for p in libpath]:
        yield pdir

# vim: set ts=4 sts=4 sw=4 tw=0:
