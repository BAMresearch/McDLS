# -*- coding: utf-8 -*-
# gui/datawidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout, QLabel)
from utils import isList
from bases.algorithm import AlgorithmBase
from gui.utils.signal import Signal
from gui.bases.mixins.titlehandler import TitleHandler
from gui.algorithmwidget import AlgorithmWidget, SettingsGridWidget 
from dataobj import DataObj, DataConfig

class DataWidget(QWidget):
    sigConfig = Signal((object,))
    _widgets = None

    def __init__(self, parent):
        super(DataWidget, self).__init__(parent)
        self.title = TitleHandler.setup(self, "Data Settings")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)

    def buildUi(self, dataobj):
        AlgorithmWidget.removeWidgets(self)
        if not isinstance(dataobj, DataObj):
            return
        
        # descriptive top row label
        lbl = QLabel(self)
        self.layout().addWidget(lbl)
        self.headLabel = lbl
        self.headLabel.setText("Configure all data sets measured by {}:"
                               .format(dataobj.sourceName))
        self._widgets = self.makeConfigUi(dataobj.config)
        self.layout().addStretch()

    onDataSelected = buildUi

    def makeConfigUi(self, config):
        if not isinstance(config, AlgorithmBase):
            return []
        # create a new layout
        w = SettingsGridWidget(self, algorithm = config)
        # first, disable the call to AlgorithmWidget.onBackendUpdate
        w.sigBackendUpdated.disconnect()
        # use this onBackendUpdate() instead which updates all subwidgets
        w.sigBackendUpdated.connect(self.onBackendUpdate)
        self.layout().addWidget(w)
        lst = [w]
        if hasattr(config, "smearing"):
            lst += self.makeConfigUi(config.smearing)
        return lst

    def onBackendUpdate(self):
        if not isList(self._widgets) or not len(self._widgets):
            return
        [w.onBackendUpdate() for w in self._widgets]
        # emit the first which contains the other(s)
        # this must not trigger onDataSelected,
        # fileWidget.sigUpdatedData disabled in MainWindow
        self.sigConfig.emit(self._widgets[0].algorithm)

# vim: set ts=4 sts=4 sw=4 tw=0:
