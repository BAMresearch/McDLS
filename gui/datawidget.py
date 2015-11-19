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
from gui.scientrybox import SciEntryBox
from gui.algorithmwidget import AlgorithmWidget, rearrangeWidgets
from dataobj import DataObj, SASConfig

import sys
from bases.algorithm.algorithmbase import AlgorithmBase

class ConfigWidget(AlgorithmWidget):

    def __init__(self, parent, algorithm, showParams = None):
        """Additional arguments: *showParams* is a list of parameter names to
        show in this widget. If not specified it shows all available parameters
        by default."""
        super(ConfigWidget, self).__init__(parent, algorithm)
        self.title = TitleHandler.setup(self, self.algorithm.name())
        layout = QGridLayout(self)
        layout.setObjectName("configLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout = layout

        if not isList(showParams) or not len(showParams):
            showParams = [p.name() for p in self.algorithm.params()]
        self._widgets = tuple(self.makeWidgets(*showParams))

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        rearrangeWidgets(self.gridLayout, self._widgets, targetWidth)

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
        # create a new layout
        w = ConfigWidget(self, algorithm = config)
        w.sigBackendUpdated.connect(self.onBackendUpdate)
        self.layout().addWidget(w)
        lst = [w]
        if getattr(config, "smearing", None) is not None:
            lst += self.makeConfigUi(config.smearing)
        return lst

    def onBackendUpdate(self):
        [w.onBackendUpdate() for w in self._widgets]
        # emit the first which contains the other
        self.sigConfig.emit(self._widgets[0].algorithm)

    def onDataSelected(self, dataobj):
        if not isinstance(dataobj, DataObj):
            return
        if dataobj.config is not None:
            print >>sys.__stderr__, "onDataSelected", dataobj.config.name()
        print >>sys.__stderr__, "onDataSelected", repr(dataobj)
        self.headLabel.setText("Configure all data sets measured by {}:"
                               .format(dataobj.sourceName))

# vim: set ts=4 sts=4 sw=4 tw=0:
