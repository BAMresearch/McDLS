# -*- coding: utf-8 -*-
# gui/datawidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.scientrybox import SciEntryBox
from gui.settingswidget import SettingsWidget

class DataWidget(SettingsWidget):

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Data settings")
        # create a new layout
        layout = QGridLayout(self)
        layout.setObjectName("dataLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)
        self._widgets = [] # containers for all inputs of one param
        entries = []
        # create inputs for a subset of calculator parameters
        # allowed parameters could be configurable from file too
        for p in ( "qMin", 
                    "qMax", 
                    "doSmear", 
                    "slitUmbra",
                    "slitPenumbra"):
            p = getattr(self.algorithm, p, None)
            if p is None: continue
            container = self.makeSetting(entries, p)
            self._widgets.append(container)

    @property
    def algorithm(self):
        return self.calculator.algo

# vim: set ts=4 sts=4 sw=4 tw=0:
