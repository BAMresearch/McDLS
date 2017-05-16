# -*- coding: utf-8 -*-
# gui/modelwidget.py

from __future__ import absolute_import # PEP328
from builtins import str
from builtins import range
import sys
import logging
import imp
import os

from gui.qt import QtCore, QtGui
from gui.utils.signal import Signal, tryDisconnect
from QtGui import (QWidget, QVBoxLayout, QComboBox)
from gui.bases.mixins import TitleHandler, AppSettings
from utils import isString
from utils.findmodels import FindModels
from bases.model import ScatteringModel
from gui.scientrybox import SciEntryBox
from gui.algorithmwidget import AlgorithmWidget
from gui.utils.displayexception import DisplayException
from dataobj import DataObj

FIXEDWIDTH = 240

class ModelWidget(AlgorithmWidget):
    _calculator = None
    _statsWidget = None # RangeList for (re-)storing histogram settings
    _models = None

    def __init__(self, parent, calculator, *args):
        super(ModelWidget, self).__init__(parent, None, *args)
        self._calculator = calculator
        self.title = TitleHandler.setup(self, "Model")
        # get all loadable ScatteringModels from model directory
        self._models = FindModels()

        layout = QVBoxLayout(self)
        layout.setObjectName("modelLayout")
        self.setLayout(layout)

        self.modelBox = QComboBox(self)
        self.modelBox.setFixedWidth(FIXEDWIDTH)
        layout.addWidget(self.modelBox)
        self.modelWidget = QWidget(self)
        paramLayout = QVBoxLayout(self.modelWidget)
        self.modelWidget.setLayout(paramLayout)
        layout.addWidget(self.modelWidget)

    def setStatsWidget(self, statsWidget):
        """Sets the statistics widget to use for updating ranges."""
        assert(isinstance(statsWidget, AppSettings))
        self._statsWidget = statsWidget

    def onDataSelected(self, dataobj):
        """Gets the data which is currently selected in the UI and rebuilds
        the model selection box based on compatible models."""
        if not isinstance(dataobj, DataObj):
            return
        try:
            self.modelBox.currentIndexChanged[int].disconnect()
        except:
            pass
        self.modelBox.clear()
        # build selection list of available models
        for modelid in self._models:
            cls, dummy = self._models[modelid]
            if cls is None or not issubclass(cls, dataobj.modelType):
                continue
            category = modelid.split('.')[0:-2]
            if category[0] == self._models.rootName:
                del category[0]
            # set up the display name of the model, may contain category
            displayName = " / ".join(category + [cls.name(),])
            # store the unique model identifier alongside its title
            self.modelBox.addItem(displayName, modelid)
        self.modelBox.setCurrentIndex(-1) # select none first
        self.modelBox.currentIndexChanged[int].connect(self._selectModelSlot)
        # trigger signal by switching from none -> 0
        self.modelBox.setCurrentIndex(0)

    def storeSession(self, section = None):
        if self.appSettings is None or self.model is None:
            return
        model = self.model.name()
        self.setRootGroup()
        self.appSettings.beginGroup(self.objectName())
        self.appSettings.setValue("model", model)
        super(ModelWidget, self).storeSession(model)
        self.appSettings.endGroup()
        self._statsWidget.storeSession()

    def restoreSession(self, model = None):
        """Load last known user settings from persistent app settings."""
        if self.appSettings is None:
            return
        if model is None:
            # get the last model used and select it
            self.appSettings.beginGroup(self.objectName())
            model = self.appSettings.value("model")
            self.appSettings.endGroup()
            # calls restoreSession(model) and storeSession()
            # mind the QSettings.group()
            if not isString(model): # default model if none set
                model = "Sphere"
            self.selectModel(model)
        else:
            self.appSettings.beginGroup(self.objectName())
            super(ModelWidget, self).restoreSession(model)
            self.appSettings.endGroup()
            self._statsWidget.restoreSession()

    def _selectModelSlot(self, key = None):
        # get the model by its unique name in the items data field
        modelid = self.modelBox.itemData(key)
        if modelid not in self._models:
            return
        model, dummy = self._models[modelid]
        if model is None or not issubclass(model, ScatteringModel):
            return
        # store current settings before changing the model
        self.storeSession()
        self._calculator.model = model() # instantiate the model class
        # remove parameter widgets from layout
        layout = self.modelWidget.layout()
        self.removeWidgets(self.modelWidget)
        # create new parameter widget based on current selection
        for p in self.model.params():
            try:
                widget = self.makeSetting(p, activeBtns = True)
                layout.addWidget(widget)
            except Exception as e:
                DisplayException(e, fmt = u"An error occurred on building a "
                    "UI for parameter '{p}' in model '{m}'."
                    "<p><nobr>{e}</nobr></p>"
                    "Please see the log for a traceback."
                    .format(p = p.name(), m = self.model.name(), e = "{e}"))
                self.removeWidgets(self.modelWidget)
                self.modelBox.setCurrentIndex(0)
                return
        layout.addStretch()
        # restore user settings for this model
        self.restoreSession(self.model.name())

    def selectModel(self, model):
        """*model*: string containing the name of the model to select.
        Calls _selectModelSlot() via signal."""
        if not isString(model):
            return
        index = 0
        # search the model with the provided name
        for i in range(0, self.modelBox.count()):
            if self.modelBox.itemText(i).lower() == model.lower().strip():
                index = i
                break
        # set the index found or the first one otherwise
        self.modelBox.setCurrentIndex(index)

    @AlgorithmWidget.algorithm.getter
    def algorithm(self):
        if self._calculator is None:
            return None
        return self._calculator.model

    model = algorithm

    def setSphericalSizeRange(self, minVal, maxVal):
        key = "radius"
        # get parameter display order of magnitude: 
        param = None
        for p in self.model.params():
            if key in p.name().lower():
                param = p
                break
        if param is None:
            logging.debug("No 'radius'-named parameter found, "
                          "not setting spherical size range!")
            return False # nothing to do
        keymin, keymax = key+"min", key+"max"
        if self.get(keymin) is not None and self.get(keymax) is not None:
            self.set(keymin, param.toDisplay(minVal))
            self.set(keymax, param.toDisplay(maxVal))
            return True
        return False

# vim: set ts=4 sts=4 sw=4 tw=0:
