# -*- coding: utf-8 -*-
# mainwindow.py

# sorry, this file is not well structured atm and is subject to change
# (the author)

import os.path
import re
import sys
import logging
from cutesnake.qt import QtCore, QtGui
from cutesnake.utils.signal import Signal
from QtCore import Qt, QSettings, QRegExp
from QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QPushButton,
                   QLabel, QCheckBox, QSizePolicy, QSpacerItem, QLayout,
                   QGroupBox, QComboBox, QApplication, QGridLayout,
                   QTreeWidgetItem, QTreeWidget, QToolBox)
from cutesnake.widgets.mainwindow import MainWindow as MainWindowBase
from cutesnake.widgets.logwidget import LogWidget
from cutesnake.widgets.datalist import DataList
from cutesnake.widgets.dockwidget import DockWidget
from cutesnake.widgets.mixins.titlehandler import TitleHandler
from cutesnake.utilsgui.filedialog import getOpenFiles
from cutesnake.widgets.settingswidget import SettingsWidget as SettingsWidgetBase
from cutesnake.utils.lastpath import LastPath
from cutesnake.utils import isList
from cutesnake.utilsgui.displayexception import DisplayException
from utils.parameter import ParameterNumerical
from version import version
from calc import Calculator
from sasdata import SASData
from mcsas.mcsas import McSAS

INFOTEXT = """One or more selected files are read in and passed to Brian Pauws Monte-Carlo size distribution analysis program for 1D SAXS data.

The convergence criterion can be set by the user. If it is not reached no output is generated, just this log is saved to file. On success, the resulting size distribution and the data fit are stored to files with uncertainties.

Output files start with the base name of the input file. They have the current date+time appended to avoid overwriting existing results."""

CHANGESTEXT = (u"""
Latest changes:
- More information shown in data list
- Correct handling when fitting multiple files
- Range estimate now uses minimum Q by itself, no longer considers Q-spacing
- LMA Dense Spheres, spherical and ellipsoidal core shell models work again
- Stability improvements and code cleanup
- Extended internal parameter functionality, using JSON defaults file
- Improvements towards implementing RangeInfo in the GUI

Changes in 0.0.11:
- distribution statistics log output and writing to stats file
- plain SAS evaluation (no fit) if no param is active
- all output files with extensions and storing in settings
- configuration file now *.cfg
- single start/stop button

Changes in 0.0.10:
- data file listing widget on top with short info
  - double-click uses sphere size est. for sphere model radius range
- start/stop buttons
  - allows to abort current calculation and restart
- column names in output files
- data file names stored with settings
- switch to enable/disable background level fit
- multiple plot figures on Windows supported

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

""".replace('\n\n', '<hr />'))
CHANGESTEXT = re.sub(r"([\s\w]*[cC]hanges.*\:)",
                     r"<strong>\1</strong>",
                     CHANGESTEXT)

from models.scatteringmodel import ScatteringModel
from models.sphere import Sphere
from models.kholodenko import Kholodenko
from models.gaussianchain import GaussianChain
from models.lmadensesphere import LMADenseSphere
from models.cylindersisotropic import CylindersIsotropic
from models.ellipsoidalcoreshell import EllipsoidalCoreShell
from models.sphericalcoreshell import SphericalCoreShell

MODELS = {Sphere.name(): Sphere,
          CylindersIsotropic.name(): CylindersIsotropic,
          EllipsoidalCoreShell.name(): EllipsoidalCoreShell,
          SphericalCoreShell.name(): SphericalCoreShell,
          GaussianChain.name(): GaussianChain,
          LMADenseSphere.name(): LMADenseSphere,
          Kholodenko.name(): Kholodenko
          }
FIXEDWIDTH = 120

def eventLoop(args):
    """Starts the UI event loop and get command line parser arguments."""
    app = QApplication(sys.argv)
    mw = MainWindow(args = args)
    mw.show()
    return app.exec_()

from cutesnake.dataset import DataSet, DisplayMixin
class ParameterRange(DataSet, DisplayMixin):
    """Represents a range tuple for a parameter of a model.
    Required for proper GUI construction."""
    _range = None

    @staticmethod
    def displayDataDescr():
        return ("lower", "upper")

    @property
    def displayData(self):
        return self.displayDataDescr()

    @property
    def lower(self):
        return self._range[0]

    @property
    def upper(self):
        return self._range[1]

    @classmethod
    def create(cls, *args):
        paramRange = cls(*args)
        return paramRange

    @staticmethod
    def sanitize(valueRange):
        return (min(valueRange), max(valueRange))

    def __init__(self, valueRange):
        self._range = self.sanitize(valueRange)
        DataSet.__init__(self, "({0}, {1})"
                               .format(*self._range), None)

    def __eq__(self, other):
        """Compares with another ParameterRange or tuple."""
        try:
            return self.lower == other.lower and self.upper == other.upper
        except AttributeError:
            return self._range == other

    def __neq__(self, other):
        return not self.__eq__(other)

from numpy import inf as numpy_inf
from QtGui import QDialog, QDoubleSpinBox
class RangeList(DataList):
    def addRange(self):
        """Creates a modal dialog window to ask the user for a range to be
        added. Returns a list of tuples (only one for now) or an empty list.
        """
        dialog = QDialog(self)
        dialog.setObjectName("AddRangeDialog")
        dialog.setWindowTitle("Add Range")
        dialog.setWindowModality(Qt.WindowModal)
        vlayout = QVBoxLayout(dialog)
        vlayout.setObjectName("vlayout")
        entryWidget = QWidget(self)
        btnWidget = QWidget(self)
        vlayout.addWidget(entryWidget)
        entryLayout = QHBoxLayout(entryWidget)
        lentry = QDoubleSpinBox(dialog)
        lentry.setPrefix("lower: ")
        lentry.setRange(-1e100, +1e100) # FIXME: float input format
        uentry = QDoubleSpinBox(dialog)
        uentry.setPrefix("upper: ")
        uentry.setRange(-1e100, +1e100)
        entryLayout.addWidget(lentry)
        entryLayout.addWidget(uentry)
        entryWidget.setLayout(entryLayout)
        vlayout.addWidget(btnWidget)
        btnLayout = QHBoxLayout(btnWidget)
        okBtn = QPushButton("add", self)
        okBtn.clicked.connect(dialog.accept)
        cancelBtn = QPushButton("cancel", self)
        cancelBtn.clicked.connect(dialog.reject)
        btnLayout.addWidget(okBtn)
        btnLayout.addWidget(cancelBtn)
        btnWidget.setLayout(btnLayout)
        dialog.setLayout(vlayout)
        lentry.selectAll() # select the first input by default
        if not dialog.exec_() or lentry.value() == uentry.value():
            return []
        return [(lentry.value(), uentry.value())]

    def loadData(self, ranges = None):
        """Overridden base class method for adding entries to the list."""
        # add only one item at a time into the list
        if ranges is None:
            ranges = self.addRange()
        # do not add duplicates
        ranges = [r for r in ranges if r not in self.data()]
        DataList.loadData(self, sourceList = ranges, showProgress = False,
                          processSourceFunc = ParameterRange.create)
        # align text to the right
        if len(self.topLevelItems()) > 0:
            item = self.topLevelItems()[-1]
            item.setTextAlignment(0, Qt.AlignRight)
            item.setTextAlignment(1, Qt.AlignRight)

    def isRemovableSelected(self):
        """Decides if selected items can be removed.
        Allow remove only if there is at least one item left."""
        return (len(self) - len(self.listWidget.selectedItems())) > 0

    def __init__(self, *args, **kwargs):
        DataList.__init__(self, *args, **kwargs)
        self.action("load").setText("add range") # fix default action name
        self.loadData([(0., numpy_inf)]) # default range
        self.listWidget.clearSelection()
        # note: derive the default range from parameters?
        # -> works only if statistics ranges are defined
        # per parameter individually, not for all as it is now

class SettingsWidget(SettingsWidgetBase):
    def __init__(self, parent):
        SettingsWidgetBase.__init__(self, parent)

    @staticmethod
    def _makeLabel(name):
        lbl = QLabel(name + ":")
        lbl.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
        return lbl

    def _makeEntry(self, name, dtype, value, minmax = None):
        ntry = self.getInputWidget(dtype)(self)
        ntry.setObjectName(name)
        if dtype is bool:
            ntry.setChecked(value)
        else:
            ntry.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
            if isList(minmax) and len(minmax):
                ntry.setMinimum(min(minmax))
                ntry.setMaximum(max(minmax))
            ntry.setValue(value)
        self.connectInputWidgets(ntry)
        return ntry

    def makeSetting(self, entries, param, activeBtns = False):
        """entries: Extended list of input widgets, for taborder elsewhere."""
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        lbl = self._makeLabel(param.displayName())
        lbl.setWordWrap(True)
        layout.addWidget(lbl)
        minmax = None
        if param.isActive() and isinstance(param, ParameterNumerical):
            # get default upper/lower bound from class
            minmax = type(param).min(), type(param).max()
            # user specified min/max within default upper/lower
            ntryMin = self._makeEntry(
                    param.name()+"min", param.dtype, param.min(), minmax)
            ntryMin.setPrefix("min: ")
            layout.addWidget(ntryMin)
            ntryMax = self._makeEntry(
                    param.name()+"max", param.dtype, param.max(), minmax)
            ntryMax.setPrefix("max: ")
            layout.addWidget(ntryMax)
            ntryMin.setFixedWidth(FIXEDWIDTH)
            ntryMax.setFixedWidth(FIXEDWIDTH)
            entries.extend((ntryMax, ntryMin)) # reversed order, see above
        else:
            if isinstance(param, ParameterNumerical):
                minmax = param.min(), param.max()
            ntry = self._makeEntry(param.name(), param.dtype, param.value(),
                                   minmax)
            ntry.setFixedWidth(FIXEDWIDTH)
            layout.addWidget(ntry)
            entries.append(ntry)
        if activeBtns:
            activeBtn = QPushButton("active", self)
            activeBtn.setObjectName(param.name()+"active")
            activeBtn.setCheckable(True)
            activeBtn.setChecked(param.isActive())
            activeBtn.setFixedWidth(FIXEDWIDTH*.5)
            layout.addWidget(activeBtn)
            activeBtn.clicked.connect(self._updateModelParams)
            self.connectInputWidgets(activeBtn)
        widget.setLayout(layout)
        return widget

    def removeLayout(self):
        if self.layout() is None:
            return
        for i in reversed(range(self.layout().count())):
            # reversed removal avoids renumbering eventually
            self.layout().takeAt(i)
        QWidget().setLayout(self.layout()) # removes layout

class AlgorithmWidget(SettingsWidget):
    def __init__(self, parent, calculator = None):
        SettingsWidget.__init__(self, parent)
        self.title = TitleHandler.setup(self, "Algorithm Settings")
        if not isinstance(calculator, Calculator):
            return
        self._widgets = []
        for i, p in enumerate(calculator.params()):
            container = self.makeSetting([], p)
            self._widgets.append(container)

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

        numCols = getNumCols()
        if self.layout() is not None:
            if numCols and self.layout().columnCount() != numCols:
                self.removeLayout()
            else:
                return
        # create a new layout
        layout = QGridLayout(self)
        layout.setObjectName("algorithmLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        # add them again with new column count
        for i, w in enumerate(self._widgets):
            layout.addWidget(w, i / numCols, i % numCols, Qt.AlignTop)
        # add empty spacer at the bottom
        layout.addWidget(QWidget(), layout.rowCount(), 0)
        layout.setRowStretch(layout.rowCount() - 1, 1)
        self.setLayout(layout)

class PropertyWidget(SettingsWidgetBase):
    _calculator = None

    def selectModel(self):
        index = [i for i in range(0, self.modelBox.count())
                    if self.modelBox.itemText(i) == "Sphere"][0]
        self.modelBox.setCurrentIndex(index)

    def keys(self):
        return self._calculator.paramNames()

    def calculator(self):
        """Returns a calculator object updated with current GUI settings."""
        self._updateModelParams()
        return self._calculator

    def setSphericalSizeRange(self, minVal, maxVal):
        if self.modelBox.currentText() != "Sphere":
            return
        key = "radius"
        keymin, keymax = key+"min", key+"max"
        if self.get(keymin) is not None and self.get(keymax) is not None:
            self.set(keymin, minVal)
            self.set(keymax, maxVal)

    def __init__(self, parent):
        SettingsWidgetBase.__init__(self, parent)
        self.title = TitleHandler.setup(self, "settings")
        self._calculator = Calculator()
        layout = QHBoxLayout(self)
        layout.setObjectName("layout")

        mcsasSettings = QGroupBox("McSAS settings")
        mcsasLayout = QVBoxLayout(mcsasSettings)
        mcsasLayout.setObjectName("mcsasLayout")
        # TODO: use signalmapper for update ...
        entries = []
        for p in self._calculator.params():
            container = self._makeSettingWidget(entries, p)
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
        self.modelBox.currentIndexChanged[str].connect(self._selectModelSlot)

        rangeStats = QGroupBox("Range Statistics")
        rangeLayout = QVBoxLayout(rangeStats)
        rangeLayout.setObjectName("rangeLayout")
        self.rangeWidget = RangeList(self, title = "ranges", withBtn = False)
        self.rangeWidget.setHeader(ParameterRange.displayDataDescr())
        rangeLayout.addWidget(self.rangeWidget)
        rangeLayout.addStretch()
        rangeStats.setLayout(rangeLayout)
        layout.addWidget(rangeStats)

        self.setLayout(layout)
        self.sigValuesChanged.connect(self._updateModelParams)
        self.sigValueChanged.connect(self._test)

    def _test(self, key):
        print >>sys.__stderr__, "changed:", key

    def _selectModelSlot(self, key = None):
        model = MODELS.get(str(key), None)
        if (key is not None and model is not None and
            issubclass(model, ScatteringModel)):
            self._calculator.model = model()
        for child in self.modelParams.children():
            if not isinstance(child, QWidget):
                continue
            self.modelParams.layout().removeWidget(child)
            child.setParent(QWidget())
        entries = []
        for p in reversed(self._calculator.modelParams()):
            container = self._makeSettingWidget(entries, p,
                                                activeBtns = True)
            self.modelParams.layout().insertWidget(0, container)
        entries.append(self.modelBox)
        for i in reversed(range(1, len(entries))):
            self.modelParams.setTabOrder(entries[i], entries[i-1])

    def _updateModelParams(self):
        activeChanged = False
        for p in self._calculator.modelParams():
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
            if isinstance(newActive, bool) and p.isActive() != newActive:
                p.setActive(newActive)
                activeChanged = True
            if p.isActive():
                # sync parameter ranges for statistics after calc
                p.histogram().resetRanges()
                for r in self.rangeWidget.data():
                    p.histogram().addRange(r.lower, r.upper)
        # update algo settings
        for p in self._calculator.params():
            value = self.get(p.name())
            if value is None:
                continue
            p.setValue(value)
        if activeChanged:
            self._selectModelSlot()

    @staticmethod
    def _makeLabel(name):
        lbl = QLabel(name)
        lbl.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
        return lbl

    def _makeEntry(self, name, dtype, value, minmax = None):
        ntry = self.getInputWidget(dtype)(self)
        ntry.setObjectName(name)
        if dtype is bool:
            ntry.setChecked(value)
        else:
            ntry.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
            if isList(minmax) and len(minmax):
                ntry.setMinimum(min(minmax))
                ntry.setMaximum(max(minmax))
            ntry.setValue(value)
        self.connectInputWidgets(ntry)
        return ntry

    def _makeSettingWidget(self, entries, param, activeBtns = False):
        """entries: Extended list of input widgets, for taborder elsewhere."""
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        lbl = self._makeLabel(param.displayName())
        lbl.setWordWrap(True)
        layout.addWidget(lbl)
        minmax = None
        if param.isActive() and isinstance(param, ParameterNumerical):
            # get default upper/lower bound from class
            minmax = type(param).min(), type(param).max()
            # user specified min/max within default upper/lower
            ntryMin = self._makeEntry(
                    param.name()+"min", param.dtype, param.min(), minmax)
            ntryMin.setPrefix("min: ")
            layout.addWidget(ntryMin)
            ntryMax = self._makeEntry(
                    param.name()+"max", param.dtype, param.max(), minmax)
            ntryMax.setPrefix("max: ")
            layout.addWidget(ntryMax)
            ntryMin.setFixedWidth(FIXEDWIDTH)
            ntryMax.setFixedWidth(FIXEDWIDTH)
            entries.extend((ntryMax, ntryMin)) # reversed order, see above
        else:
            if isinstance(param, ParameterNumerical):
                minmax = param.min(), param.max()
            ntry = self._makeEntry(param.name(), param.dtype, param.value(),
                                   minmax)
            ntry.setFixedWidth(FIXEDWIDTH)
            layout.addWidget(ntry)
            entries.append(ntry)
        if activeBtns:
            activeBtn = QPushButton("active", self)
            activeBtn.setObjectName(param.name()+"active")
            activeBtn.setCheckable(True)
            activeBtn.setChecked(param.isActive())
            activeBtn.setFixedWidth(FIXEDWIDTH*.5)
            layout.addWidget(activeBtn)
            activeBtn.clicked.connect(self._updateModelParams)
            self.connectInputWidgets(activeBtn)
        widget.setLayout(layout)
        return widget

class FileList(DataList):
    sigSphericalSizeRange = Signal((float, float))

    def loadData(self, fileList = None):
        if fileList is None or type(fileList) is bool:
            fileList = getOpenFiles(self, "Load one or more data files",
                                    LastPath.get(), multiple = True)
        # populates to data list widget with items based on the return of
        # processSourceFunc(filename)
        DataList.loadData(self, sourceList = fileList, showProgress = False,
                          processSourceFunc = SASData.load)

    def itemDoubleClicked(self, item, column):
        valueRange = item.data().sphericalSizeEst()
        self.sigSphericalSizeRange.emit(min(valueRange), max(valueRange))

class MainWindow(MainWindowBase):
    onCloseSignal = Signal()
    _args = None # python command line arguments parser
    _calculator = None # calculator calling algorithm on all data

    def __init__(self, parent = None, args = None):
        # calls setupUi() and restoreSettings()
        MainWindowBase.__init__(self, version, parent)
        self._args = args

    @property
    def calculator(self):
        """Returns a calculator object."""
        if self._calculator is None:
            self._calculator = Calculator()
        return self._calculator

    def setupUi(self, *args):
        # called in MainWindowBase.__init__()
        # file widget at the top
        self.toolbox = QToolBox(self)
        self._addToolboxItem(self._setupFileWidget())
        self._addToolboxItem(self._setupAlgoWidget())
        self._addToolboxItem(self._setupSettings())
        self.fileWidget.sigSphericalSizeRange.connect(
                        self.propWidget.setSphericalSizeRange)
        # put the log widget at the bottom
        self.addDockWidget(Qt.BottomDockWidgetArea, self._setupLogWidget())

        # set up central widget of the main window
        self.centralLayout = QVBoxLayout()
        # put buttons in central widget
        self.centralLayout.addWidget(self.toolbox)
        self.centralLayout.addWidget(self._setupButtons())
        centralWidget = QWidget(self)
        centralWidget.setLayout(self.centralLayout)
        centralWidget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        self.setCentralWidget(centralWidget)
        self.onStartupSignal.connect(self.initUI)

    def _addToolboxItem(self, widget):
        self.toolbox.addItem(widget, "{n}. {t}"
                                     .format(n = self.toolbox.count()+1,
                                             t = widget.title()))
        try:
            widget.layout().setContentsMargins(0, 0, 0, 0)
        except: pass

    def _setupSettings(self):
        """Set up property widget with settings."""
        propWidget = PropertyWidget(self)
        self.propWidget = propWidget
        return propWidget

    def _setupFileWidget(self):
        # set up file widget
        fileWidget = FileList(self, title = "data files", withBtn = False)
        fileWidget.setHeader(SASData.displayDataDescr())
        fileWidget.setToolTip(
                "Double click to use the estimated size for the model.")
        self.fileWidget = fileWidget
        return fileWidget

    def _setupAlgoWidget(self):
        """Set up property widget with settings."""
        self.algoWidget = AlgorithmWidget(self, self.calculator)
        return self.algoWidget

    def _setupLogWidget(self):
        """Set up widget for logging output."""
        logDock = DockWidget(self, LogWidget, appversion = version)
        logWidget = logDock.child
        self.onCloseSignal.connect(logWidget.onCloseSlot)
        logWidget.setSizePolicy(QSizePolicy.Preferred,
                                QSizePolicy.Expanding)
        logWidget.append(INFOTEXT)
        if len(CHANGESTEXT):
            logWidget.append(CHANGESTEXT)
        logWidget.append("\n\r")
        self.logWidget = logWidget
        return logDock

    def _setupButtons(self):
        """Set up buttons."""
        self.loadBtn = QPushButton("load files ...")
        self.loadBtn.pressed.connect(self.fileWidget.loadData)
        self.startStopBtn = QPushButton()
        self.startStopBtn.setCheckable(True)
        self.startStopBtn.clicked[bool].connect(self.onStartStopClick)
        btnLayout = QHBoxLayout()
        btnLayout.setContentsMargins(0, 0, 0, 0)
        for btn in self.loadBtn, self.startStopBtn:
            btn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
            btnLayout.addWidget(btn)
        btnWidget = QWidget(self)
        btnWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        btnWidget.setLayout(btnLayout)
        return btnWidget

    def restoreSettings(self):
        MainWindowBase.restoreSettings(self)
        settings = self.appSettings()
        for name in self.propWidget.keys():
            if self.propWidget.get(name) is None:
                continue # no UI element for this setting
            val = settings.value(name)
            self.propWidget.set(name, settings.value(name))
        try:
            value = unicode(settings.value("lastpath").toString())
        except AttributeError: # QVariant
            value = unicode(settings.value("lastpath"))
        if os.path.isdir(value):
            LastPath.set(value)

    def storeSettings(self):
        MainWindowBase.storeSettings(self)
        settings = self.appSettings()
        for name in self.propWidget.keys():
            value = self.propWidget.get(name)
            if value is None:
                continue # no UI element for this setting
            settings.setValue(name, value)
        settings.setValue("lastpath", LastPath.get())
        settings.sync()
        tempSettings = QSettings("/tmp/qsettings.test", QSettings.IniFormat)
        for key in settings.allKeys():
            if key in ('geometry', 'windowState', 'lastpath'):
                continue
            tempSettings.setValue(key, settings.value(key))
        tempSettings.sync()

    def fileDialog(self):
        filenames = getOpenFiles(self, "Load one or more data files",
                                 LastPath.get(), multiple = True)
        self.loadFiles(filenames)

    def initUI(self):
        self.propWidget.selectModel()
        self.fileWidget.loadData(getattr(self._args, "fnames", []))
        self.onStartStopClick(getattr(self._args, "start", False))
        self.logWidget.scrollToTop()

    def onStartStopClick(self, checked):
        if checked:
            self.startStopBtn.setText("stop")
            self.startStopBtn.setChecked(True)
            self.calc()
        # run this also for 'start' after calculation
        self.propWidget.calculator().stop()
        self.startStopBtn.setText("start")
        self.startStopBtn.setChecked(False)

    def calc(self):
        if len(self.fileWidget) <= 0:
            return
        self.logWidget.clear()
        calculator = self.propWidget.calculator()
        self.fileWidget.updateData(updateFunc = calculator,
                                   showProgress = False)

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onStartStopClick(False)
        self.onCloseSignal.emit()

# vim: set ts=4 sts=4 sw=4 tw=0:
