# -*- coding: utf-8 -*-
# gui/datawidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout, QLabel)
from utils import isList
from gui.utils.signal import Signal
from gui.bases.mixins.titlehandler import TitleHandler
from gui.algorithmwidget import AlgorithmWidget, SettingsGridWidget 
from dataobj import DataObj, SASConfig

import sys

class DataWidget(QWidget):
    sigConfig = Signal(SASConfig)

    def __init__(self, parent):
        super(DataWidget, self).__init__(parent)
        self.title = TitleHandler.setup(self, "Data Settings")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)

        self.buildUi(SASConfig())

    def buildUi(self, config):
        AlgorithmWidget.removeWidgets(self)
        
        # descriptive top row label
        lbl = QLabel(self)
        self.layout().addWidget(lbl)
        self.headLabel = lbl

        self._widgets = self.makeConfigUi(config)
        self.layout().addStretch()

    def makeConfigUi(self, config):
        if config is None: # not isinstance(config, DataConfig)
            return []
        # create a new layout
        w = SettingsGridWidget(self, algorithm = config)
        w.sigBackendUpdated.connect(self.onBackendUpdate)
        self.layout().addWidget(w)
        lst = [w]
        if getattr(config, "smearing", None) is not None:
            lst += self.makeConfigUi(config.smearing)
        return lst

    def onBackendUpdate(self):
        if not isList(self._widgets) or not len(self._widgets):
            return
        [w.onBackendUpdate() for w in self._widgets]
        # emit the first which contains the other
        self.sigConfig.emit(self._widgets[0].algorithm)

    def onDataSelected(self, dataobj):
        if not isinstance(dataobj, DataObj):
            return
        self.headLabel.setText("Configure all data sets measured by {}:"
                               .format(dataobj.sourceName))
        self.buildUi(dataobj.config)


# vim: set ts=4 sts=4 sw=4 tw=0:
