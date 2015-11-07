# -*- coding: utf-8 -*-
# gui/datawidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.scientrybox import SciEntryBox
from gui.settingswidget import SettingsWidget, rearrangeWidgets
from dataobj import SASConfig

import sys

class SmearingWidget(SettingsWidget):
    _smearingConfig = None

    @property
    def algorithm(self):
        return self._smearingConfig

    def __init__(self, *args, **kwargs):
        super(SmearingWidget, self).__init__(*args, **kwargs)
        self.title = TitleHandler.setup(self, "Smearing Settings")
        self._smearingConfig = kwargs.pop("smearingConfig", None)
        layout = QGridLayout(self)
        layout.setObjectName("configLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self._widgets = tuple(self.makeWidgets("umbra", "penumbra"))

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        rearrangeWidgets(self.layout(), self._widgets, targetWidth)

    def onUpdate(self, paramTuple):
        param = paramTuple[0]
        print >>sys.__stderr__, "onUpdate", param.name()
        print >>sys.__stderr__, unicode(self._smearingConfig)
        self.updateUi()

class DataWidget(SettingsWidget):
    _dataConfig = None

    @property
    def algorithm(self):
        return self._dataConfig

    def __init__(self, *args, **kwargs):
        super(DataWidget, self).__init__(*args, **kwargs)
        self.title = TitleHandler.setup(self, "Data Settings")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(hlayout)

        self._dataConfig = SASConfig()
        # create a new layout
        configWidget = QWidget()
        layout = QGridLayout(configWidget)
        layout.setObjectName("configLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.configLayout = layout
        self._widgets = tuple(self.makeWidgets("qLow", "qHigh"))
        self.sigBackendUpdated.connect(self.onUpdate)
        hlayout.addWidget(configWidget)

        kwargs["smearingConfig"] = self._dataConfig.smearing
        self.smearingWidget = SmearingWidget(*args, **kwargs)
        hlayout.addWidget(self.smearingWidget)

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        rearrangeWidgets(self.configLayout, self._widgets, targetWidth)

    def onUpdate(self):
        print >>sys.__stderr__, unicode(self._dataConfig)
        self.smearingWidget.updateUi()

# vim: set ts=4 sts=4 sw=4 tw=0:
