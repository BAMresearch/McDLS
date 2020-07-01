# -*- coding: utf-8 -*-
# gui/utils/__init__.py

from __future__ import absolute_import # PEP328
from QtWidgets import QApplication

def processEventLoop():
    QApplication.processEvents()
#    QApplication.sendPostedEvents()
#    QApplication.flush()

# vim: set ts=4 sw=4 sts=4 tw=0:
