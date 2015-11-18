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
from gui.settingswidget import SettingsWidget, rearrangeWidgets
from dataobj import DataObj, SASConfig

import sys
from bases.algorithm.algorithmbase import AlgorithmBase

class ConfigWidget(SettingsWidget):
    _algo = None

    @property
    def algorithm(self):
        return self._algo

    @algorithm.setter
    def algorithm(self, algo):
        assert isinstance(algo, AlgorithmBase)
        self._algo = algo

    def __init__(self, *args, **kwargs):
        """Additional arguments: *algorithm* and ??"""
        super(ConfigWidget, self).__init__(*args, **kwargs)
        self.algorithm = kwargs.pop("algorithm", None)
        self.title = TitleHandler.setup(self, self.algorithm.name())
        layout = QGridLayout(self)
        layout.setObjectName("configLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout = layout

        pNames = kwargs.pop("showParams", None)
        if not isList(pNames) or not len(pNames):
            pNames = [p.name() for p in self.algorithm.params()]
        print >>sys.__stderr__, "ConfigWidget", pNames
        self._widgets = tuple(self.makeWidgets(*pNames))

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        rearrangeWidgets(self.gridLayout, self._widgets, targetWidth)

class DataWidget(QWidget):
    sigConfig = Signal(SASConfig)

    def __init__(self, parent, calculator, *args, **kwargs):
        super(DataWidget, self).__init__(parent, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Data Settings")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)

        args = (self, calculator) + args
        self.buildUi(SASConfig(), args)

    def buildUi(self, config, args):
        SettingsWidget.removeWidgets(self)
        
        # descriptive top row label
        lbl = QLabel(self)
        self.layout().addWidget(lbl)
        self.headLabel = lbl

        self._widgets = self.makeConfigUi(config, args)
        self.layout().addStretch()

    def makeConfigUi(self, config, args):
        # create a new layout
        w = ConfigWidget(*args, algorithm = config)
        w.sigBackendUpdated.connect(self.onBackendUpdate)
        self.layout().addWidget(w)
        lst = [w]
        if getattr(config, "smearing", None) is not None:
            lst += self.makeConfigUi(config.smearing, args)
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
