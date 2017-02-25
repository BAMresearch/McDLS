# -*- coding: utf-8 -*-
# gui/utils/translate.py

from __future__ import absolute_import # PEP328
from gui.qt import QtGui, QtCore
from QtGui import QApplication

def tr(s):
    return QApplication.translate(None, s)

fromUtf8 = str

# vim: set ts=4 sts=4 sw=4 tw=0:
