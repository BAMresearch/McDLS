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
from cutesnake.utils import isList, isString
from cutesnake.utils.tests import testfor
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
    _parameter = None #identifier to which parameter the range applies

    @staticmethod
    def displayDataDescr():
        return ("lower", "upper", "parameter")

    @property
    def displayData(self):
        return self.displayDataDescr()

    @property
    def lower(self):
        return self._range[0]

    @property
    def upper(self):
        return self._range[1]

    @property
    def parameter(self):
        return self._parameter

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

    def __str__(self):
        return str(self._range)

    __repr__ = __str__

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
        #need to identify to which parameter the range limits apply.
        pentry = QDoubleSpinBox(dialog) # FIXME: set to parameter pulldown
        pentry.setPrefix("parameter: ")
        pentry.setRange(0, 100)
        entryLayout.addWidget(lentry)
        entryLayout.addWidget(uentry)
        entryLayout.addWidget(pentry)
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
        return [(lentry.value(), uentry.value())], pentry.value() 

    def loadData(self, ranges = None):
        """Overridden base class method for adding entries to the list."""
        # add only one item at a time into the list
        if ranges is None:
            ranges, pval = self.addRange()
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

    def setupUi(self):
        self.clearSelection()
        self.listWidget.setRootIsDecorated(False)
        self.listWidget.setUniformRowHeights(True)
        self.listWidget.setItemsExpandable(False)
        self.listWidget.setAlternatingRowColors(True)
        self.action("load").setText("add range") # fix default action name
        self.loadData([(0., numpy_inf)]) # default range
        # note: derive the default range from parameters?
        # -> works only if statistics ranges are defined
        # per parameter individually, not for all as it is now

class SettingsWidget(SettingsWidgetBase):
    _calculator = None # calculator instance associated
    _appSettings = None

    def __init__(self, parent, calculator = None):
        SettingsWidgetBase.__init__(self, parent)
        self.sigValueChanged.connect(self._updateParam)
        assert isinstance(calculator, Calculator)
        self._calculator = calculator

    @property
    def appSettings(self):
        return self._appSettings

    @appSettings.setter
    def appSettings(self, settings):
        assert isinstance(settings, QSettings)
        self._appSettings = settings

    @property
    def calculator(self):
        return self._calculator

    @property
    def algorithm(self):
        """Retrieves AlgorithmBase object containing all parameters
        for this settings."""
        raise NotImplementedError

    @property
    def keys(self):
        """Returns all existing input names (for store/restore)."""
        if self.algorithm is None:
            return
        for p in self.algorithm.params():
            query = p.name()
            try:
                p.isActive() # fails for non-FitParameters
                query = QRegExp(p.name()+".*")
            except: pass
            for w in self.findChildren(QWidget, query):
                yield w.objectName()

    def storeSession(self, section = None):
        """Stores current UI configuration to persistent application settings.
        """
        if self.appSettings is None:
            return
        if section is None:
            section = self.objectName()
        self.appSettings.beginGroup(section)
        for key in self.keys:
            value = self.get(key)
            self.appSettings.setValue(key, value)
        self.appSettings.endGroup()

    def restoreSession(self, section = None):
        if self.appSettings is None:
            return
        if section is None:
            section = self.objectName()
        self.appSettings.beginGroup(section)
        for key in self.keys:
            value = self.appSettings.value(key)
            self.set(key, value)
        self.appSettings.endGroup()

    def _updateParam(self, widget):
        # get the parameter instance associated with this widget
        if self.algorithm is None:
            return
        p = None
        try:
            p = getattr(self.algorithm, widget.parameterName)
        except AttributeError:
            p = None
        if p is None:
            return
        self.updateParam(widget, p)

    def updateParam(self, widget, p):
        """Write UI settings back to the algorithm."""
        activeChanged = False
        # persistent name due to changes to the class instead of instance
        key = p.name()
        # get the parent of the updated widget and other input for this param
        parent = widget.parent()
        valueWidget = parent.findChild(QWidget, key)
        # get changed value range if any
        minValue = self.get(key+"min")
        maxValue = self.get(key+"max")
        if minValue is not None and maxValue is not None:
            # update value range for numerical parameters
            p.setValueRange((minValue, maxValue))
            # update bounds of the value input widget
            valueWidget.setMinimum(p.min())
            valueWidget.setMaximum(p.max())
        # update the value input widget itself
        newValue = self.get(key)
        if newValue is not None:
            p.setValue(newValue)
        # fit parameter related updates
        newActive = self.get(key+"active")
        minWidget = parent.findChild(QWidget, key+"min")
        maxWidget = parent.findChild(QWidget, key+"max")
        if isinstance(newActive, bool): # None for non-fit parameters
            # update active state for fit parameters
            p.setActive(newActive)
            if p.isActive():
                valueWidget.hide()
                minWidget.show()
                maxWidget.show()
                # sync parameter ranges for statistics after calc
                self._updateStatsRanges()
        if not isinstance(newActive, bool) or not newActive:
            valueWidget.show()
            try:
                minWidget.hide()
                maxWidget.hide()
            except: pass

    def setStatsWidget(self, statsWidget):
        """Sets the statistics widget to use for updating ranges."""
        assert(isinstance(statsWidget, DataList))
        self._statsWidget = statsWidget
        statsWidget.sigEditingFinished.connect(self._updateStatsRanges)

    def _updateStatsRanges(self, count = None, index = None):
        """Sets all statistics ranges to all parameters in the associated
        algorithm. Uses the previously configured statistics widget."""
        if self._statsWidget is None:
            return
        for p in self.algorithm.params():
            try:
                # works for active FitParameters only
                p.histogram().resetRanges()
            except:
                continue
            for r in self._statsWidget.data():
                p.histogram().addRange(r.lower, r.upper)

    @staticmethod
    def _makeLabel(name):
        lbl = QLabel(name + ":")
        lbl.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
        lbl.setWordWrap(True)
        return lbl

    def _makeEntry(self, name, dtype, value, minmax = None,
                   widgetType = None, parent = None):
        testfor(name not in self.keys, KeyError,
            "Input widget '{w}' exists already in '{s}'"
            .format(w = name, s = self.objectName()))
        if widgetType is None:
            widgetType = self.getInputWidget(dtype)
        if parent is None:
            parent = self
        widget = widgetType(parent)
        widget.setObjectName(name)
        if dtype is bool:
            widget.setCheckable(True)
            widget.setChecked(value)
        else:
            widget.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
            if isList(minmax) and len(minmax):
                widget.setMinimum(min(minmax))
                widget.setMaximum(max(minmax))
            widget.setValue(value)
        self.connectInputWidgets(widget)
        return widget

    def makeSetting(self, entries, param, activeBtns = False):
        """entries: Extended list of input widgets, for taborder elsewhere."""
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._makeLabel(param.displayName()))
        widget.setLayout(layout)
        minmaxValue, widgets = None, []
        # instead of create/remove widgets, show/hide them on active toggle
        if isinstance(param, ParameterNumerical):
            # create input boxes for user specified min/max
            # within default upper/lower from class definition
            minmaxValue = type(param).min(), type(param).max()
            for bound in "min", "max":
                w = self._makeEntry(param.name() + bound, param.dtype,
                                       getattr(param, bound)(),
                                       minmax = minmaxValue, parent = widget)
                w.setPrefix(bound + ": ")
                w.setFixedWidth(FIXEDWIDTH)
                widgets.append(w)
        # create scalar value input widget
        w = self._makeEntry(param.name(), param.dtype, param.value(),
                                      minmax = minmaxValue, parent = widget)
        w.setFixedWidth(FIXEDWIDTH)
        widgets.insert(len(widgets)/2, w)
        if activeBtns:
            w = self._makeEntry(param.name()+"active", bool,
                                param.isActive(),
                                widgetType = QPushButton,
                                parent = widget)
            w.setText("active")
            w.setFixedWidth(FIXEDWIDTH*.5)
            widgets.append(w)
        # add input widgets to the layout
        for w in widgets:
            layout.addWidget(w)
            entries.append(w)
            # store the parameter name
            w.parameterName = param.name()
        # configure UI accordingly (hide/show widgets)
        self.updateParam(widgets[-1], param)
        return widget

    @staticmethod
    def clearLayout(layout, newParent = None):
        """Removes all widgets from the given layout and reparents them if
        *newParent* is a sub class of QWidget"""
        if layout is None:
            return
        assert isinstance(layout, QLayout)
        if not isinstance(newParent, QWidget):
            newParent = None
        for i in reversed(range(layout.count())):
            # reversed removal avoids renumbering eventually
            item = layout.takeAt(i)
            if newParent is not None:
                try:
                    item.widget().setParent(newParent)
                except: pass

    @staticmethod
    def removeWidgets(widget):
        """Removes all widgets from the layout of the given widget."""
        SettingsWidget.clearLayout(widget.layout(), QWidget())

class AlgorithmWidget(SettingsWidget):

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Algorithm Settings")
        # create a new layout
        layout = QGridLayout(self)
        layout.setObjectName("algorithmLayout")
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)
        self._widgets = [] # containers for all inputs of one param
        entries = []
        # create inputs for a subset of calculator parameters
        # allowed parameters could be configurable from file too
        for i, p in enumerate(("convergenceCriterion", "histogramBins",
                              "numReps", "numContribs", "findBackground")):
            p = getattr(self.algorithm, p)
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

class ModelWidget(SettingsWidget):

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Model Settings")

        layout = QVBoxLayout(self)
        layout.setObjectName("modelLayout")
        self.setLayout(layout)

        self.modelBox = QComboBox(self)
        for name in MODELS.iterkeys():
            self.modelBox.addItem(name)
        self.modelBox.setFixedWidth(FIXEDWIDTH)
        layout.addWidget(self.modelBox)

        self.modelWidget = QWidget(self)
        paramLayout = QVBoxLayout(self.modelWidget)
        self.modelWidget.setLayout(paramLayout)
        layout.addWidget(self.modelWidget)

        self.modelBox.setCurrentIndex(-1)
        self.modelBox.currentIndexChanged[str].connect(self._selectModelSlot)

    def storeSession(self, section = None):
        if self.appSettings is None or self.algorithm is None:
            return
        model = self.algorithm.name()
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
        model = MODELS.get(str(key), None)
        if model is None or not issubclass(model, ScatteringModel):
            return
        # store current settings before changing the model
        self.storeSession()
        self.calculator.model = model()
        # remove parameter widgets from layout
        layout = self.modelWidget.layout()
        self.removeWidgets(self.modelWidget)
        # create new parameter widget based on current selection
        entries = [self.modelBox]
        for p in self.algorithm.params():
            widget = self.makeSetting(entries, p,
                                      activeBtns = True)
            layout.addWidget(widget)
        layout.addStretch()
        # restore user settings for this model
        self.restoreSession(self.calculator.model.name())

    def selectModel(self, model):
        """*model*: string containing the name of the model to select."""
        if not isString(model):
            return
        for i in range(0, self.modelBox.count()):
            if self.modelBox.itemText(i).lower() == model.lower().strip():
                self.modelBox.setCurrentIndex(i)
                break

    @property
    def algorithm(self):
        if self.calculator is None:
            return None
        return self.calculator.model

    def setSphericalSizeRange(self, minVal, maxVal):
        if self.modelBox.currentText() != "Sphere":
            return
        key = "radius"
        keymin, keymax = key+"min", key+"max"
        if self.get(keymin) is not None and self.get(keymax) is not None:
            self.set(keymin, minVal)
            self.set(keymax, maxVal)

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
        self._addToolboxItem(self._setupModelWidget())
        self._addToolboxItem(self._setupStatsWidget())
        # put the log widget at the bottom
        self.addDockWidget(Qt.BottomDockWidgetArea, self._setupLogWidget())

        # set up central widget of the main window
        self.centralLayout = QVBoxLayout()
        # put buttons in central widget
        self.centralLayout.addWidget(self.toolbox)
        self.centralLayout.addWidget(self._setupStartButton())
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

    def _setupFileWidget(self):
        # set up file widget
        fileWidget = FileList(self, title = "Data Files",
                              withBtn = False, nestedItems = False)
        fileWidget.setHeader(SASData.displayDataDescr())
        fileWidget.setToolTip(
                "Right-click to add datafiles.\n" +
                "Double click to use the estimated size for the model.")
        self.fileWidget = fileWidget
        return fileWidget

    def _setupAlgoWidget(self):
        """Set up property widget with settings."""
        self.algoWidget = AlgorithmWidget(self, self.calculator)
        return self.algoWidget

    def _setupModelWidget(self):
        """Set up property widget with settings."""
        self.modelWidget = ModelWidget(self, self.calculator)
        self.fileWidget.sigSphericalSizeRange.connect(
                self.modelWidget.setSphericalSizeRange)
        return self.modelWidget

    def _setupStatsWidget(self):
        """Set up property widget with settings."""
        # setup similar to the file widget
        self.statsWidget = RangeList(self,
                                     title = "Set up Statistics",
                                     withBtn = False, nestedItems = False)
        self.statsWidget.setToolTip(
                "Right-click to add additional ranges.\n" +
                "Keeping full range (0, inf) is highly recommended.")
        self.statsWidget.setHeader(ParameterRange.displayDataDescr())
        self.modelWidget.setStatsWidget(self.statsWidget)
        return self.statsWidget

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

    def _setupStartButton(self):
        """Set up "Start/Stop" - button."""
        self.startStopBtn = QPushButton()
        self.startStopBtn.setCheckable(True)
        self.startStopBtn.clicked[bool].connect(self.onStartStopClick)
        btnLayout = QHBoxLayout()
        btnLayout.setContentsMargins(0, 0, 0, 0)
        self.startStopBtn.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        btnLayout.addWidget(self.startStopBtn)
        btnWidget = QWidget(self)
        btnWidget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        btnWidget.setLayout(btnLayout)
        return btnWidget

    def restoreSettings(self):
        MainWindowBase.restoreSettings(self)
        settings = self.appSettings()
        for settingsWidget in self.algoWidget, self.modelWidget:
            settingsWidget.appSettings = self.appSettings()
            settingsWidget.restoreSession()
        try:
            value = unicode(settings.value("lastpath").toString())
        except AttributeError: # QVariant
            value = unicode(settings.value("lastpath"))
        if os.path.isdir(value):
            LastPath.set(value)

    def storeSettings(self):
        MainWindowBase.storeSettings(self)
        settings = self.appSettings()
        for settingsWidget in self.algoWidget, self.modelWidget:
            settingsWidget.storeSession()
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
        self.fileWidget.loadData(getattr(self._args, "fnames", []))
        self.onStartStopClick(getattr(self._args, "start", False))
        self.logWidget.scrollToTop()

    def onStartStopClick(self, checked):
        if checked:
            self.startStopBtn.setText("stop")
            self.startStopBtn.setChecked(True)
            self.calc()
        # run this also for 'start' after calculation
        self.calculator.stop()
        self.startStopBtn.setText("start")
        self.startStopBtn.setChecked(False)

    def calc(self):
        if len(self.fileWidget) <= 0:
            return
        self.logWidget.clear()
        self.fileWidget.updateData(updateFunc = self.calculator,
                                   showProgress = False)

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onStartStopClick(False)
        self.onCloseSignal.emit()

# vim: set ts=4 sts=4 sw=4 tw=0:
