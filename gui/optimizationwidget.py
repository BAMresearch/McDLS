# -*- coding: utf-8 -*-
# gui/optimizationwidget.py

from __future__ import absolute_import # PEP328

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout, QGroupBox)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.settingswidget import SettingsWidget, rearrangeWidgets

class SettingsGroup(object):
    _widgets = None # preconfigured widgets to show in this group

    def __init__(self, *args, **kwargs):
        self._widgets = kwargs.pop("widgets", None)
        super(SettingsGroup, self).__init__(*args, **kwargs)
        layout = QGridLayout(self)
        layout.setObjectName(self.objectName() + "Layout")
        topMargin = 10
        layout.setContentsMargins(0, topMargin, 0, 0)

    def rearrangeWidgets(self, targetWidth):
        rearrangeWidgets(self.layout(), self._widgets, targetWidth)

class DefaultSettings(SettingsGroup, QWidget):
    def __init__(self, *args, **kwargs):
        super(DefaultSettings, self).__init__(*args, **kwargs)

class AdvancedSettings(SettingsGroup, QGroupBox):
    def __init__(self, *args, **kwargs):
        super(AdvancedSettings, self).__init__(*args, **kwargs)
        self.setTitle("Advanced Settings")
#        set by SettingsWidget._makeEntry()
#        self.setCheckable(True)
#        self.setChecked(False)
        self.clicked[bool].connect(self.showAdvanced)
        self.updateWidgets()

    def updateWidgets(self, widget = None):
        if widget is None or id(widget) == id(self):
            self.showAdvanced(self.isChecked())

    def showAdvanced(self, show = False):
        func = QWidget.hide
        if show:
            func = QWidget.show
        for w in self._widgets:
            func(w)

class OptimizationWidget(SettingsWidget):

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
            "convergenceCriterion", "numReps", "autoClose")))
        # creating an ui entry with settings persistence via store/restore
        self.advanced = self._makeEntry("showAdvanced", bool, False,
            widgetType = AdvancedSettings, widgets = tuple(self.makeWidgets(
                "numContribs", "compensationExponent", "findBackground", "maxIterations")))
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
