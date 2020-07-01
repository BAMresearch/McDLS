# -*- coding: utf-8 -*-
# gui/settingsgroup.py

from __future__ import absolute_import # PEP328

from builtins import object
from QtCore import Qt
from QtWidgets import (QWidget, QGridLayout, QGroupBox)
from gui.algorithmwidget import rearrangeWidgets

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
        self.setStyleSheet("QGroupBox { border: 0px; padding-top: 10px; }")
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

# vim: set ts=4 sts=4 sw=4 tw=0:
