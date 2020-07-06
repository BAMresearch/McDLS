# -*- coding: utf-8 -*-
# gui/optimizationwidget.py

from QtCore import Qt
from QtWidgets import (QWidget, QGridLayout, QVBoxLayout, QGroupBox)

from .bases.mixins.titlehandler import TitleHandler
from .algorithmwidget import AlgorithmWidget, rearrangeWidgets
from .settingsgroup import DefaultSettings, AdvancedSettings
from ..dataobj import DataObj, SASConfig

class OptimizationWidget(AlgorithmWidget):
    _tempCompExp = None

    @property
    def uiWidgets(self):
        if hasattr(self, "advanced"):
            return [self.advanced,]
        else:
            return []

    def __init__(self, *args):
        super(OptimizationWidget, self).__init__(*args)
        self.title = TitleHandler.setup(self, "Optimization")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(hlayout)

        self.defaults = DefaultSettings(self, widgets = tuple(self.makeWidgets(
            "convergenceCriterion", "numReps", "autoClose")))
        # creating an ui entry with settings persistence via store/restore
        self.advanced = self._makeEntry("showAdvanced", bool, False,
            widgetType = AdvancedSettings, widgets = tuple(self.makeWidgets(
                "numContribs", "compensationExponent", 
                "findBackground", "positiveBackground", "maxIterations", 
                "showIncomplete", "seriesStats")))
        hlayout.addWidget(self.defaults)
        hlayout.addWidget(self.advanced)
        self.sigValueChanged.connect(self.advanced.updateWidgets)

    def resizeWidgets(self, targetWidth):
        """Creates a new layout with appropriate row/column count."""
        self.defaults.rearrangeWidgets(targetWidth)
        self.advanced.rearrangeWidgets(targetWidth)
        # add empty spacer at the bottom
        self.layout().addStretch()

    def onDataSelected(self, dataobj):
        """Sets defaults for certain types of DataConfig selected,
        respectively fixes some values."""
        if not isinstance(dataobj, DataObj):
            return
        if isinstance(dataobj.config, SASConfig):
            if self._tempCompExp is not None: # restore the original value
                self.set("compensationExponent", self._tempCompExp)
        elif self.get("compensationExponent", 1.0) != 1.0:
            # remember the original value
            self._tempCompExp = self.get("compensationExponent", None)
            self.set("compensationExponent", 1.0) # fix its value

    def storeSession(self, *args, **kwargs):
        self.setRootGroup()
        super(OptimizationWidget, self).storeSession(*args, **kwargs)

    def restoreSession(self, *args, **kwargs):
        super(OptimizationWidget, self).restoreSession(*args, **kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
