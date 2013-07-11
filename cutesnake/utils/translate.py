# -*- coding: utf-8 -*-
# translate.py

from cutesnake.qt import QtGui, QtCore
from QtCore import QString
from QtGui import QApplication

def tr(s):
    return QApplication.translate(None, s)

try:
    fromUtf8 = QString.fromUtf8
except AttributeError:
    fromUtf8 = lambda s: s

# vim: set ts=4 sts=4 sw=4 tw=0:
