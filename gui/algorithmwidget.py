# -*- coding: utf-8 -*-
# gui/algorithmwidget.py

from __future__ import absolute_import # PEP328

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout, QVBoxLayout, QGroupBox)
from gui.bases.mixins.titlehandler import TitleHandler
from gui.settingswidget import SettingsWidget

def rearrangeWidgets(layout, widgets, targetWidth):
    def getNumCols():
        width = 0
        for i, w in enumerate(widgets):
            width += w.sizeHint().width()
            if width > targetWidth:
                return i
        return len(widgets)

    SettingsWidget.clearLayout(layout)
    numCols = max(1, getNumCols())
    # add them again with new column count
    for i, w in enumerate(widgets):
        layout.addWidget(w, i / numCols, i % numCols, Qt.AlignTop)
    # add empty spacer at the bottom
    layout.addWidget(QWidget(), layout.rowCount(), 0)
    layout.setRowStretch(layout.rowCount() - 1, 1)

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

class AlgorithmWidget(SettingsWidget):

    @property
    def algorithm(self):
        return self.calculator.algo

    @property
    def uiWidgets(self):
        if hasattr(self, "advanced"):
            return (self.advanced,)
        else:
            return ()

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Algorithm")
        # basic row oriented layout
        hlayout = QVBoxLayout(self)
        hlayout.setObjectName("baseLayout")
        hlayout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(hlayout)

        self.defaults = DefaultSettings(self, widgets = tuple(self._widgets(
            "convergenceCriterion", "numReps", "autoClose")))
        # creating an ui entry with settings persistence via store/restore
        self.advanced = self._makeEntry("showAdvanced", bool, False,
            widgetType = AdvancedSettings, widgets = tuple(self._widgets(
                "numContribs", "compensationExponent", "findBackground")))
        hlayout.addWidget(self.defaults)
        hlayout.addWidget(self.advanced)
        self.sigValueChanged.connect(self.advanced.updateWidgets)

    # create inputs for a subset of calculator parameters
    # allowed parameters could be configurable from file too
    def _widgets(self, *args):
        for p in args:
            p = getattr(self.algorithm, p, None)
            if p is None: continue
            yield self.makeSetting(p)

    def resizeEvent(self, resizeEvent):
        """Creates a new layout with appropriate row/column count."""
        # basically, reacts to the size change by spawning scroll bar
        targetWidth = resizeEvent.size().width()
        self.defaults.rearrangeWidgets(targetWidth)
        self.advanced.rearrangeWidgets(targetWidth)

# vim: set ts=4 sts=4 sw=4 tw=0:
