# -*- coding: utf-8 -*-
# mainwindow.py

# sorry, this file is not well structured atm and is subject to change
# (the author)

import os.path
from cutesnake.qt import QtCore, QtGui
from cutesnake.utils.signal import Signal
from QtCore import Qt, QSettings, QString, QRegExp
from QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QPushButton,
                   QLabel, QCheckBox, QSizePolicy, QSpacerItem, QLayout,
                   QGroupBox, QComboBox)
from cutesnake.widgets.mainwindow import MainWindow as MainWindowBase
from cutesnake.widgets.logwidget import LogWidget
from cutesnake.utilsgui.filedialog import getOpenFiles
from cutesnake.widgets.settingswidget import SettingsWidget
from cutesnake.utils.lastpath import LastPath
from version import version
from calc import calc, SASData

INFOTEXT="""One or more selected files are read in and passed to Brian Pauws Monte-Carlo size distribution analysis program for 1D SAXS data.

The convergence criterion can be set by the user. If it is not reached no output is generated, just this log is saved to file. On success, the resulting size distribution and the data fit are stored to files with uncertainties.

Output files start with the base name of the input file. They have the current date+time appended to avoid overwriting existing results."""

CHANGESTEXT=QString(u"""

Changes in 0.0.9:
- added GUI to public McSAS repository

Changes in 0.0.8:
- new model: GaussianChain, verified against SASfit:
  http://sasfit.sf.net/manual/Gaussian_Chain#Gauss_2
- fixed volume function exponent in Kholodenko
  was v² instead of just v

Changes in 0.0.7:
- Using restructured McSAS
- building GUI for models and global settings dynamically
- new model: Kholodenkos worm-like structure, verified against SASfit

Changes in 0.0.6:
Fixed handling of negative values in PDH data files.

Changes in 0.0.5:
'Number-weighted distributions now come with correct-looking observability limits.'
 https://bitbucket.org/pkwasniew/mcsas/commits/81bbf84

""".replace('\n\n', '<hr />')).replace(
        QRegExp(r"(Changes in [0-9]\.[0-9]\.[0-9])"), r"<strong>\1</strong>")

from models.sphere import Sphere
from models.kholodenko import Kholodenko
from models.gaussianchain import GaussianChain
from models.lmadensesphere import LMADenseSphere

MODELS = {Sphere.name(): Sphere,
          Kholodenko.name(): Kholodenko,
          GaussianChain.name(): GaussianChain,
          LMADenseSphere.name(): LMADenseSphere}
FIXEDWIDTH = 120

class PropertyWidget(SettingsWidget):
    _optional = None
    _mcsasKeys = ("convergenceCriterion", "histogramBins", "numReps",
                  "numContribs")

    def keys(self):
        return self._mcsasKeys

    def selectModelSlot(self, key = None):
        model = MODELS.get(str(key), None)
        if key is not None and model is not None:
            SASData.mcsas.model = model()
        for child in self.modelParams.children():
            if not isinstance(child, QWidget):
                continue
            self.modelParams.layout().removeWidget(child)
            child.setParent(QWidget())
        for ptype in reversed(SASData.mcsas.model.params()):
            p = getattr(SASData.mcsas.model, ptype.name())
            container = self._makeSettingWidget(p, activeBtns = True)
            self.modelParams.layout().insertWidget(0, container)

    def updateModelParams(self):
        activeChanged = False
        for p in SASData.mcsas.model.params():
            # persistent due to changes to the class instead of instance
            key = p.name()
            newValue = self.get(key)
            if newValue is not None:
                p.setValue(newValue)
            minValue = self.get(key+"min")
            maxValue = self.get(key+"max")
            if minValue is not None and maxValue is not None:
                p.setValueRange((minValue, maxValue))
            newActive = self.get(key+"active")
            if isinstance(newActive, bool) and p.isActive != newActive:
                p.isActive = newActive
                activeChanged = True
        if activeChanged:
            self.selectModelSlot()

    def selectModel(self):
#        self.modelBox.setCurrentIndex(0)
        index = [i for i in range(0, self.modelBox.count())
                    if self.modelBox.itemText(i) == "Sphere"][0]
        self.modelBox.setCurrentIndex(index)

    def __init__(self, parent):
        SettingsWidget.__init__(self, parent)
        layout = QHBoxLayout(self)
        layout.setObjectName("layout")

        mcsasSettings = QGroupBox("McSAS settings")
        mcsasLayout = QVBoxLayout(mcsasSettings)
        mcsasLayout.setObjectName("mcsasLayout")
        for key in self._mcsasKeys:
            p = getattr(SASData.mcsas, key)
            container = self._makeSettingWidget(p)
            mcsasLayout.addWidget(container)
        mcsasSettings.setLayout(mcsasLayout)
        layout.addWidget(mcsasSettings)

        modelSettings = QGroupBox("Model settings")
        modelLayout = QVBoxLayout(modelSettings)
        modelLayout.setObjectName("modelLayout")
        self.modelBox = QComboBox(modelSettings)
        for name in MODELS.iterkeys():
            self.modelBox.addItem(name)
        self.modelBox.setFixedWidth(FIXEDWIDTH)
        modelLayout.addWidget(self.modelBox)
        self.modelParams = QWidget()
        paramLayout = QVBoxLayout(self.modelParams)
        self.modelParams.setLayout(paramLayout)
        self.modelParams.layout().addStretch()
        modelLayout.addWidget(self.modelParams)
        modelSettings.setLayout(modelLayout)
        layout.addWidget(modelSettings)
        self.modelBox.setCurrentIndex(-1)
        self.modelBox.currentIndexChanged[str].connect(self.selectModelSlot)
        self.setLayout(layout)
        self.sigValuesChanged.connect(self.updateModelParams)

    @staticmethod
    def _makeLabel(name):
        lbl = QLabel(name)
        lbl.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
        return lbl

    def _makeEntry(self, name, dtype, minvalue, maxvalue, value):
        ntry = self.getInputWidget(dtype)(self)
        ntry.setObjectName(name)
        ntry.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        ntry.setMaximum(maxvalue)
        ntry.setMinimum(minvalue)
        ntry.setValue(value)
        self.connectInputWidgets(ntry)
        return ntry

    def _makeSettingWidget(self, param, activeBtns = False):
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        lbl = self._makeLabel(param.displayName())
        lbl.setWordWrap(True)
        layout.addWidget(lbl)
        if param.isActive:
            ntryMin = self._makeEntry(
                    param.name()+"min", param.dtype,
                    type(param).min(), type(param).max(),
                    param.min())
            ntryMin.setPrefix("min: ")
            layout.addWidget(ntryMin)
            ntryMax = self._makeEntry(
                    param.name()+"max", param.dtype,
                    type(param).min(), type(param).max(),
                    param.max())
            ntryMax.setPrefix("max: ")
            layout.addWidget(ntryMax)
            ntryMin.setFixedWidth(FIXEDWIDTH)
            ntryMax.setFixedWidth(FIXEDWIDTH)
        else:
            ntry = self._makeEntry(param.name(), param.dtype,
                                   param.min(), param.max(),
                                   param.value())
            ntry.setFixedWidth(FIXEDWIDTH)
            layout.addWidget(ntry)
        if activeBtns:
            activeBtn = QPushButton("active", self)
            activeBtn.setObjectName(param.name()+"active")
            activeBtn.setCheckable(True)
            activeBtn.setChecked(param.isActive)
            activeBtn.setFixedWidth(FIXEDWIDTH*.5)
            layout.addWidget(activeBtn)
            activeBtn.clicked.connect(self.updateModelParams)
            self.connectInputWidgets(activeBtn)
        widget.setLayout(layout)
        return widget

class MainWindow(MainWindowBase):
    onCloseSignal = Signal()

    def __init__(self, parent = None):
        MainWindowBase.__init__(self, version, parent)

    def setupUi(self, *args):
        btnWidget = QWidget(self)
        btnLayout = QHBoxLayout()
        self.loadBtn = QPushButton("browse files ...")
        self.loadBtn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        self.propWidget = PropertyWidget(self)
        self.propWidget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        btnLayout.addWidget(self.loadBtn)
        btnLayout.addWidget(self.propWidget)
        btnWidget.setLayout(btnLayout)

        centralWidget = QWidget(self)
        centralLayout = QVBoxLayout()
        self.logWidget = LogWidget(self, appversion = version)
        self.onCloseSignal.connect(self.logWidget.onCloseSlot)
        self.logWidget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        centralLayout.addWidget(btnWidget)
        centralLayout.addWidget(self.logWidget)
        centralWidget.setLayout(centralLayout)

        self.setCentralWidget(centralWidget)

        self.loadBtn.pressed.connect(self.fileDialog)
        SASData.logWidget = self.logWidget
        SASData.settings = self.propWidget
        self.logWidget.append(INFOTEXT)
        if len(CHANGESTEXT):
            self.logWidget.append(CHANGESTEXT)
        self.logWidget.append("\n\r")

    def restoreSettings(self):
        MainWindowBase.restoreSettings(self)
        for name in ("convcrit", "nreps", "min", "max", "bins") + self.propWidget.keys():
            if self.propWidget.get(name) is None:
                continue
            self.propWidget.set(name, self.appSettings().value(name))
        value = self.appSettings().value("lastpath").toString()
        if os.path.isdir(value):
            LastPath.set(value)

    def storeSettings(self):
        MainWindowBase.storeSettings(self)
        for name in ("convcrit", "nreps", "min", "max", "bins") + self.propWidget.keys():
            self.appSettings().setValue(name, self.propWidget.get(name))
        self.appSettings().setValue("lastpath", LastPath.get())
        settings = self.appSettings()
        tempSettings = QSettings("/tmp/qsettings.test", QSettings.IniFormat)
        for key in settings.allKeys():
            if key in ('geometry', 'windowState', 'lastpath'):
                continue
            tempSettings.setValue(key, settings.value(key))
        tempSettings.sync()

    def fileDialog(self):
        filenames = getOpenFiles(self, "Load one or more data files",
                                 LastPath.get(), multiple = True)
        calc(filenames)

    def onStartup(self):
        self.propWidget.selectModel()
        calc(self.getCommandlineArguments())

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onCloseSignal.emit()

# vim: set ts=4 sts=4 sw=4 tw=0:
