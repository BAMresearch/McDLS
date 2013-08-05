# -*- coding: utf-8 -*-
# translate.py

from cutesnake.qt import QtGui, QtCore
from QtGui import QApplication

def tr(s):
    return QApplication.translate(None, s)

fromUtf8 = unicode

# vim: set ts=4 sts=4 sw=4 tw=0:
