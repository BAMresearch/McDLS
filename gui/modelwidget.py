# -*- coding: utf-8 -*-
# gui/modelwidget.py

from __future__ import absolute_import # PEP328
import sys
import logging
import imp
import os

from gui.qt import QtCore, QtGui
from gui.utils.signal import Signal
from QtGui import (QWidget, QVBoxLayout, QComboBox)
from gui.bases.mixins.titlehandler import TitleHandler
from utils import isString
from utils.findmodels import FindModels
from models.scatteringmodel import ScatteringModel
from gui.scientrybox import SciEntryBox
from gui.algorithmwidget import AlgorithmWidget
from dataobj import DataObj

FIXEDWIDTH = 240

from utils.devtools import DBG

class ModelWidget(AlgorithmWidget):
    sigModelChanged = Signal()
    _calculator = None
    _models = None

    def __init__(self, parent, calculator):
        super(ModelWidget, self).__init__(parent, None)
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

    def onDataSelected(self, dataobj):
        """Gets the data which is currently selected in the UI and rebuilds
        the model selection box based on compatible models."""
        if not isinstance(dataobj, DataObj):
            return
        try:
            self.modelBox.currentIndexChanged[str].disconnect()
        except:
            pass
        self.modelBox.clear()
        # build selection list of available models
        for modelid in self._models:
            DBG(modelid)
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
        self.appSettings.beginGroup(self.objectName())
        self.appSettings.setValue("model", model)
        super(ModelWidget, self).storeSession(model)
        self.appSettings.endGroup()

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

    def _selectModelSlot(self, key = None):
        # rebuild the UI without early signals
        try:
            self.sigRangeChanged.disconnect(self.sigModelChanged)
        except: sys.exc_clear()
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
            widget = self.makeSetting(p, activeBtns = True)
            layout.addWidget(widget)
        layout.addStretch()
        # restore user settings for this model
        self.restoreSession(self.model.name())
        self.sigRangeChanged.connect(self.sigModelChanged)
        self.sigModelChanged.emit()

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
