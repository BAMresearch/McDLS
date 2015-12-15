# -*- coding: utf-8 -*-
# gui/optimizationwidget.py

from __future__ import absolute_import # PEP328

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout, QGroupBox)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.algorithmwidget import AlgorithmWidget, rearrangeWidgets
from gui.settingsgroup import DefaultSettings, AdvancedSettings


class OptimizationWidget(AlgorithmWidget):

    @property
    def uiWidgets(self):
        if hasattr(self, "advanced"):
            return (self.advanced,)
        else:
            return ()

    def __init__(self, *args):
        super(OptimizationWidget, self).__init__(*args)
        self.title = TitleHandler.setup(self, "Optimization")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(hlayout)

        self.defaults = DefaultSettings(self, widgets = tuple(self.makeWidgets(
            "convergenceCriterion", "numReps", "autoClose", "seriesStats")))
        # creating an ui entry with settings persistence via store/restore
        self.advanced = self._makeEntry("showAdvanced", bool, False,
            widgetType = AdvancedSettings, widgets = tuple(self.makeWidgets(
                "numContribs", "compensationExponent", 
                "findBackground", "maxIterations", "showIncomplete")))
        hlayout.addWidget(self.defaults)
        hlayout.addWidget(self.advanced)
        self.sigValueChanged.connect(self.advanced.updateWidgets)

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        self.defaults.rearrangeWidgets(targetWidth)
        self.advanced.rearrangeWidgets(targetWidth)
        # add empty spacer at the bottom
        self.layout().addStretch()

# vim: set ts=4 sts=4 sw=4 tw=0:
