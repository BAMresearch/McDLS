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

from gui.scientrybox import SciEntryBox
from collections import OrderedDict

"""
from models.scatteringmodel import ScatteringModel
from models.sphere import Sphere
from models.kholodenko import Kholodenko
from models.gaussianchain import GaussianChain
from models.lmadensesphere import LMADenseSphere
from models.cylindersisotropic import CylindersIsotropic
from models.cylindersradiallyisotropic import CylindersRadiallyIsotropic
from models.ellipsoidsisotropic import EllipsoidsIsotropic
from models.ellipsoidalcoreshell import EllipsoidalCoreShell
from models.sphericalcoreshell import SphericalCoreShell

MODELS = OrderedDict((
    (Sphere.name(), Sphere),
    (CylindersIsotropic.name(), CylindersIsotropic),
    (CylindersRadiallyIsotropic.name(), CylindersRadiallyIsotropic),
    (EllipsoidsIsotropic.name(), EllipsoidsIsotropic),
    (EllipsoidalCoreShell.name(), EllipsoidalCoreShell),
    (SphericalCoreShell.name(), SphericalCoreShell),
    (GaussianChain.name(), GaussianChain),
    (LMADenseSphere.name(), LMADenseSphere),
    (Kholodenko.name(), Kholodenko),
))
"""
MODELS = OrderedDict()

FIXEDWIDTH = 120

from gui.algorithmwidget import AlgorithmWidget
from dataobj import DataObj

class ModelWidget(AlgorithmWidget):
    sigModelChanged = Signal()
    _calculator = None
    _modelNames = OrderedDict()
    _modelPaths = OrderedDict()

    def __init__(self, parent, calculator):
        super(ModelWidget, self).__init__(parent, None)
        FindModels()
        self._modelPaths = FindModels._orderedModelFiles
        self.importModels()
        self._calculator = calculator
        self.title = TitleHandler.setup(self, "Model")

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

    def importModels(self):
        # tries to import an ordered dict of models
        for mName in self._modelPaths.keys():
            print("trying to import mName: {} from path: {}".format(mName, self._modelPaths[mName]))
            # try importing
            mFunc = self.importModel(mName, self._modelPaths[mName])
            if mFunc is not None:
                print("importing mFunc: {}, using mName: {}".format(mFunc, mName))
                # MODELS.update({mFunc.name(): mFunc}) # complains it has no attribute "name"
                MODELS.update({mName: mFunc})
        

    def importModel(self, modelName, uri):
        path, fname = os.path.split(uri)

        fp, pathname, description = imp.find_module(modelName.lower(), [path])

        try:
            package = imp.load_module(modelName, fp, pathname, description)
            print("success trying to import model: {} from {}".format(modelName, fp))
        except:
            package = None
            print("failed trying to import model: {}".format(modelName))
        finally:
            if fp:
                fp.close()

        return package

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
        for name, cls in MODELS.iteritems():
            if cls is None or not issubclass(cls, dataobj.modelType):
                continue
            self.modelBox.addItem(name)
        self.modelBox.setCurrentIndex(-1) # select none first
        self.modelBox.currentIndexChanged[str].connect(self._selectModelSlot)
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
        model = MODELS.get(str(key), None)
        if model is None or not issubclass(model, ScatteringModel):
            return
        # store current settings before changing the model
        self.storeSession()
        self._calculator.model = model()
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
