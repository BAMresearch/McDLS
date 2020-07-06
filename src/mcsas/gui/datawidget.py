# -*- coding: utf-8 -*-
# gui/datawidget.py

import logging
from .qt import QtCore, QtGui
from QtCore import Qt
from QtWidgets import (QWidget, QGridLayout, QVBoxLayout, QLabel)

from ..utils import isList
from ..bases.algorithm import AlgorithmBase
from .utils.signal import Signal
from .bases.mixins import TitleHandler, AppSettings
from .algorithmwidget import AlgorithmWidget, SettingsGridWidget 
from ..dataobj import DataObj, DataConfig

class DataWidget(QWidget, AppSettings):
    sigConfig = Signal((object,))
    # contains child widgets built from given DataConfig instances
    _widgets = None
    # remembers which algorithm UI settings were restored already
    # restore algorithm settings only once during program lifetime
    _restored = None

    def __init__(self, parent, appSettings):
        super(DataWidget, self).__init__(parent)
        self.title = TitleHandler.setup(self, "Data Settings")
        self.appSettings = appSettings # forwarded to SettingsGridWidget below
        self._widgets = []
        self._restored = set()
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)

    def clearUi(self):
        # remember the latest values before a widget is removed
        self.storeSession()
        AlgorithmWidget.removeWidgets(self)
        # prevents storeSession() calls on already removed widgets
        self._widgets = []

    def buildUi(self, dataobj):
        self.clearUi()
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
        self.restoreSession()

    onDataSelected = buildUi

    def makeConfigUi(self, config):
        if not isinstance(config, AlgorithmBase):
            return []
        # create a new layout
        w = SettingsGridWidget(self, algorithm = config,
                               appSettings = self.appSettings)
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
        if not len(self._widgets):
            return
        [w.onBackendUpdate() for w in self._widgets]
        # emit the first DataConfig which contains the other(s) as well
        # this must not trigger onDataSelected,
        # fileWidget.sigUpdatedData disabled in MainWindow
        self.sigConfig.emit(self._widgets[0].algorithm)

    def onEmptyDataList(self):
        """Forgets which data settings were already restored."""
        self._restored = set()

    def storeSession(self):
        if not len(self._widgets):
            return
        self.setRootGroup()
        for w in self._widgets:
            w.storeSession()

    def restoreSession(self):
        if not len(self._widgets):
            return
        self.setRootGroup()
        for w in self._widgets:
            if type(w.algorithm) in self._restored:
                continue # restore settings only once
            w.restoreSession()
            self._restored.add(type(w.algorithm))

# vim: set ts=4 sts=4 sw=4 tw=0:
