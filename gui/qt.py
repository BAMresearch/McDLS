# -*- coding: utf-8 -*-
# gui/qt.py

from __future__ import absolute_import # PEP328

__all__ = []

import sys
import os
import importlib

try:
    import PySide
    from PySide import QtGui, QtCore, QtSvg, QtXml
    from PySide import __path__ as PySidePath
except ModuleNotFoundError: # PySide is outdated with python3
    import PySide2
    from PySide2 import QtGui, QtCore, QtWidgets, QtSvg, QtXml
    from PySide2 import __path__ as PySidePath

try:
    sys.modules["QtGui"] = QtGui
    if 'QtWidgets' in locals():
        sys.modules["QtGui"] = QtWidgets
    sys.modules["QtGui5"] = QtGui
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
