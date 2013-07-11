# -*- coding: utf-8 -*-
# __init__.py

from cutesnake.qt import QtGui
from QtGui import QApplication
from time import sleep

def processEventLoop():
    QApplication.processEvents()
#    QApplication.sendPostedEvents()
#    QApplication.flush()

# vim: set ts=4 sw=4 sts=4 tw=0:
