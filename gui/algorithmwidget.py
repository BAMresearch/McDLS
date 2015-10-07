# -*- coding: utf-8 -*-
# gui/algorithmwidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import (QWidget, QGridLayout)
from gui.bases.mixins.titlehandler import TitleHandler

from gui.scientrybox import SciEntryBox

from gui.settingswidget import SettingsWidget

class AlgorithmWidget(SettingsWidget):

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Algorithm")
        # create a new layout
        layout = QGridLayout(self)
        layout.setObjectName("algorithmLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)
        self._widgets = [] # containers for all inputs of one param
        entries = []
        # create inputs for a subset of calculator parameters
        # allowed parameters could be configurable from file too
        #for i, p in enumerate(("qUnits", "iUnits", "convergenceCriterion", 
        for i, p in enumerate(( "convergenceCriterion", 
                    "numReps", 
                    "numContribs", 
                    "findBackground",
                    "autoClose")):
            p = getattr(self.algorithm, p, None)
            if p is None: continue
            container = self.makeSetting(entries, p)
            self._widgets.append(container)

    @property
    def algorithm(self):
        return self.calculator.algo

    def resizeEvent(self, resizeEvent):
        """Creates a new layout with appropriate row/column count."""
        # basically, reacts to the size change by spawning scroll bar
        scrollArea = self.parent().parent()
        targetWidth = resizeEvent.size().width()
        def getNumCols():
            width = 0
            for i, w in enumerate(self._widgets):
                width += w.sizeHint().width()
                if width > targetWidth:
                    return i
            return len(self._widgets)

        layout = self.layout()
        numCols = max(1, getNumCols())
        if layout.count() <= 0 or layout.columnCount() != numCols:
            self.clearLayout(layout)
        else:
            return
        # add them again with new column count
        for i, w in enumerate(self._widgets):
            layout.addWidget(w, i / numCols, i % numCols, Qt.AlignTop)
        # add empty spacer at the bottom
        layout.addWidget(QWidget(), layout.rowCount(), 0)
        layout.setRowStretch(layout.rowCount() - 1, 1)

# vim: set ts=4 sts=4 sw=4 tw=0:
