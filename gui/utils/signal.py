# -*- coding: utf-8 -*-
# gui/utils/signal.py

from __future__ import absolute_import # PEP328
from gui.qt import QtCore

if not hasattr(QtCore, 'Signal'):  ## for pyside compatibility
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.Slot = QtCore.pyqtSlot

Signal = QtCore.Signal

# vim: set ts=4 sts=4 sw=4 tw=0:
