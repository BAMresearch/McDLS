# -*- coding: utf-8 -*-
# gui/mainwindow.py

from __future__ import absolute_import # PEP328
import os.path
import re
import sys
import logging

from numpy import inf as numpy_inf
from gui.qt import QtCore, QtGui
from gui.utils.signal import Signal
from QtCore import Qt, QSettings, QRegExp, QFileInfo
from QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QPushButton,
                   QLabel, QCheckBox, QSizePolicy, QSpacerItem, QLayout,
                   QGroupBox, QComboBox, QApplication, QGridLayout,
                   QTreeWidgetItem, QTreeWidget, QToolBox, QPalette,
                   QDialog, QDoubleSpinBox, QSpinBox)
from bases.dataset import DataSet, DisplayMixin
from gui.bases.mainwindow import MainWindow as MainWindowBase
from gui.bases.logwidget import LogWidget
from gui.bases.datalist import DataList
from gui.bases.dockwidget import DockWidget
from gui.bases.mixins.titlehandler import TitleHandler
from gui.bases.settingswidget import SettingsWidget as SettingsWidgetBase
from gui.utils.filedialog import getOpenFiles
from utils.lastpath import LastPath
from bases.algorithm.parameter import ParameterFloat # instance for test
from utils import isList, isString, testfor
from gui.utils import processEventLoop
from gui.utils.displayexception import DisplayException
from utils.parameter import (ParameterNumerical, Histogram, FitParameterBase,
                                FitParameterNumerical)
from gui.version import version
from gui.calc import Calculator
from sasdata import SASData
from mcsas.mcsas import McSAS

from numpy import inf as numpy_inf
from QtGui import QDialog, QDoubleSpinBox, QSpinBox, QLineEdit, QDoubleValidator
from gui.scientrybox import SciEntryBox

INFOTEXT = """One or more selected files are read in and passed to Brian Pauws Monte-Carlo size distribution analysis program for 1D SAXS data.

The convergence criterion can be set by the user. If it is not reached no output is generated, just this log is saved to file. On success, the resulting size distribution and the data fit are stored to files with uncertainties.

Output files start with the base name of the input file. They have the current date+time appended to avoid overwriting existing results."""

CHANGESTEXT = (u"""
Latest changes:
- Number input boxes now allow scientific notation input
- Stability improvements and code cleanup
- Improved plotting routine
- Range statistics list functionality added
- Improved division of GUI items into vertical tabs
- More information shown in data list
- Correct handling when fitting multiple files
- Range estimate now uses minimum Q by itself, no longer considers Q-spacing
- LMA Dense Spheres, spherical and ellipsoidal core shell models work again
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
from models.cylindersradiallyisotropic import CylindersRadiallyIsotropic
from models.ellipsoidsisotropic import EllipsoidsIsotropic
from models.ellipsoidalcoreshell import EllipsoidalCoreShell
from models.sphericalcoreshell import SphericalCoreShell
from collections import OrderedDict

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
FIXEDWIDTH = 120

# required for svg graphics support
from gui.qt import QtSvg, QtXml, pluginDirs
from main import makeAbsolutePath

def eventLoop(args):
    """Starts the UI event loop and get command line parser arguments."""
    app = QApplication(sys.argv)
    for pluginDir in pluginDirs(): # required for svg graphics support
        app.addLibraryPath(pluginDir)
    mw = MainWindow(args = args)
    mw.show()
    return app.exec_()

def makeAlternatingRowColorsTransparent(widget):
    palette = widget.palette()
    color = palette.color(QPalette.AlternateBase)
    color.setAlphaF(0.4)
    palette.setColor(QPalette.AlternateBase, color)
    widget.setPalette(palette)

def setBackgroundStyleSheet(widget, imgpath):
    assert isinstance(widget, QWidget)
    makeAlternatingRowColorsTransparent(widget.listWidget)
    stylesheet = """
        #listWidget {{
            background-image:       url({path});
            background-repeat:      no-repeat;
            background-position:    center center;
            background-attachment:  fixed;
            background-color:       white;
        }}
    """
    # convert path to qt style formatting (separators, ...)
    imgpath = QFileInfo(makeAbsolutePath(imgpath)).absoluteFilePath()
    widget.setStyleSheet(stylesheet.format(path = imgpath))

class RangeDialog(QDialog):
    """Creates a modal dialog window to ask the user for a range to be
    added."""
    _model = None

    def __init__(self, parent = None, model = None):
        QDialog.__init__(self, parent)
        assert isinstance(model, ScatteringModel)
        self._model = model
        self.setObjectName("AddRangeDialog")
        self.setWindowTitle("Add Range")
        self.setWindowModality(Qt.WindowModal)
        vlayout = QVBoxLayout(self)
        vlayout.setObjectName("vlayout")
        vlayout.addWidget(self._createEntries())
        vlayout.addWidget(self._createButtons())
        self.setLayout(vlayout)

    def _createEntries(self):
        entryWidget = QWidget(self)
        entryLayout = QHBoxLayout(entryWidget)
        entryLayout.addWidget(self._createParamBox())
        entryLayout.addWidget(self._createLower())
        entryLayout.addWidget(self._createUpper())
        entryLayout.addWidget(self._createBins())
        entryLayout.addWidget(self._createXScale())
        entryLayout.addWidget(self._createYWeight())
        entryWidget.setLayout(entryLayout)
        self.pbox.setCurrentIndex(0)
        self.lentry.selectAll() # select the first input by default
        return entryWidget

    def _createLower(self):
        # add input for lower limit
        lentry = SciEntryBox(self)
        lentry.setPrefix("lower: ")
        self.lentry = lentry
        return lentry

    def _createUpper(self):
        # add input for upper limit
        uentry = SciEntryBox(self)
        uentry.setPrefix("upper: ")
        self.uentry = uentry
        return uentry

    def _createBins(self):
        # number of histogram bin input box
        bentry = QSpinBox(self)
        bentry.setPrefix("# bins: ")
        bentry.setRange(1, 200)
        bentry.setValue(50)
        bentry.setSingleStep(10)
        self.bentry = bentry
        return bentry

    def _createParamBox(self):
        # add a parameter choice list
        pbox = QComboBox(self)
        for p in self._model.params():
            if not hasattr(p, "isActive"):
                #not a fit parameter, a regular one
                continue
            if not p.isActive():
                continue
            # providing the internal param name as item data
            pbox.addItem(p.displayName(), p.name())
        # pbox.addItem("Test Dummy", "dummy")
        pbox.setCurrentIndex(-1)
        pbox.currentIndexChanged[int].connect(self._selectParam)
        self.pbox = pbox
        return pbox

    def _createXScale(self):
        # histogram scaling choice (X-axis only at this point)
        sentry = QComboBox(self)
        for name in Histogram.xscaling():
            sentry.addItem(name)
        self.sentry = sentry
        return sentry

    def _createYWeight(self):
        # histogram weighting input
        wentry = QComboBox(self)
        for name in Histogram.yweighting():
            wentry.addItem(name)
        self.wentry = wentry
        return wentry

    def _selectParam(self, index):
        """Gets the index within the selection box."""
        pname = self.pbox.itemData(index)
        p = getattr(self._model, pname)
        # perhaps, use testfor() for that:
        assert p is not None, "Could not find parameter from selection box" 
        if isinstance(p, ParameterFloat):
            # account for units conversion:
            # llim, ulim = type(p).displayValueRange()
            llim, ulim = p.displayValueRange()
            lval, uval = p.displayActiveRange()
        else:
            # llim, ulim = type(p).valueRange()
            llim, ulim = p.valueRange()
            lval, uval = p.activeRange()
        self.lentry.setRange(llim, ulim)
        self.uentry.setRange(llim, ulim)
        self.lentry.setValue(lval)
        self.uentry.setValue(uval)

    def _createButtons(self):
        btnWidget = QWidget(self)
        btnLayout = QHBoxLayout(btnWidget)
        okBtn = QPushButton("add", self)
        okBtn.clicked.connect(self.accept)
        cancelBtn = QPushButton("cancel", self)
        cancelBtn.clicked.connect(self.reject)
        btnLayout.addWidget(okBtn)
        btnLayout.addWidget(cancelBtn)
        btnWidget.setLayout(btnLayout)
        return btnWidget

    def output(self):
        if not self.exec_() or self.lentry.value() == self.uentry.value():
            return None
        p = None
        try:
            pname = self.pbox.itemData(self.pbox.currentIndex())
            p = getattr(self._model, pname)
        except:
            return None
        if isinstance(p, ParameterFloat):
            # take units into account,
            # convert from display units to SI units for internal use
            lval, uval = (p.toSi(self.lentry.value()),
                          p.toSi(self.uentry.value()))
        else: 
            lval, uval = (self.lentry.value(), self.uentry.value())

        output = Histogram(p, lval, uval,
                self.bentry.value(), self.sentry.currentText(),
                self.wentry.currentText())
        return output

class RangeList(DataList):
    _calculator = None

    def __init__(self, calculator = None, **kwargs):
        DataList.__init__(self, **kwargs)
        assert isinstance(calculator, Calculator)
        self._calculator = calculator
        self.sigRemovedData.connect(self.onRemoval)

    def onRemoval(self, removedHistograms):
        for hist in removedHistograms:
            hist.param.histograms().remove(hist)

    def updateHistograms(self):
        self.clear()
        lst = []
        for p in self._calculator.model.activeParams():
            lst += [h for h in p.histograms()]
        self.append(lst)

    def append(self, histList):
        if histList is None:
            return
        if not isList(histList):
            histList = [histList]
        DataList.loadData(self, sourceList = histList, showProgress = False,
                processSourceFunc = lambda x: x)
        self.fitColumnsToContents()

    def loadData(self, ranges = None):
        """Overridden base class method for adding entries to the list."""
        # add only one item at a time into the list
        # set up dialog window for "add range"
        dialog = RangeDialog(self, self._calculator.model)
        newHist = dialog.output()
        if not isinstance(newHist, Histogram):
            return
        # do not add duplicates
        if newHist in self.data():
            return
        # add it to the actual histogram list of the parameter
        newHist.param.histograms().append(newHist) # hmm, funny here
        # update the GUI based on that
        self.updateHistograms()

    def setupUi(self):
        setBackgroundStyleSheet(self, "./resources/background_ranges.svg")
        self.listWidget.setRootIsDecorated(False)
        self.listWidget.setUniformRowHeights(True)
        self.listWidget.setItemsExpandable(False)
        self.listWidget.setAlternatingRowColors(True)
        self.action("load").setText("add range") # fix default action name
# FIXME        self.loadData([(0., numpy_inf)]) # default range
        self.clearSelection()
        self.setHeader(Histogram.displayDataDescr())
        self.setToolTip(
            "Right-click to add additional ranges."
        )
        # self.addMenuEntry(name = "edit", text = "Edit selected", 
        #                   menuStates = "hasSelection",
        #                   callbacks = self.editEntry)
        self.addMenuEntry(name = "recalc", text = "recalc selected",
                          menuStates = "hasSelection",
                          callbacks = self.recalc)

    def editEntry(self):
        return
        # more on the TODO list... Will figure this out later.
    
    def recalc(self):
        return
        # does not work yet! missing arguments for Hist.calc()
        DataList.updateData(self, selectedOnly = True, showProgress = False,
                updateFunc = Histogram.calc,
                stopFunc = None)

class SettingsWidget(SettingsWidgetBase):
    _calculator = None # calculator instance associated
    _appSettings = None
    sigRangeChanged = Signal()

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
                query = QRegExp("^" + p.name() + ".*")
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
            try:
                self.set(key, value)
            except StandardError, e:
                logging.warn(e)
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
        def isNotNone(*args):
            return all((a is not None for a in args))
        # persistent name due to changes to the class instead of instance
        key = p.name()
        # get the parent of the updated widget and other input for this param
        parent = widget.parent()
        valueWidget = parent.findChild(QWidget, key)
        minWidget = parent.findChild(QWidget, key+"min")
        maxWidget = parent.findChild(QWidget, key+"max")
        # get changed value range if any
        minValue = self.get(key+"min")
        maxValue = self.get(key+"max")
        # update value range for numerical parameters
        if isinstance(p, FitParameterBase):
            if isinstance(p, ParameterFloat):
                if isNotNone(minValue, maxValue):
                    p.setDisplayActiveRange((minValue, maxValue))
                clippedVals = p.displayActiveRange() # get updated values
            else:
                if isNotNone(minValue, maxValue):
                    p.setActiveRange((minValue, maxValue))
                clippedVals = p.activeRange()
            # somehow move clippedVals back to widgets, does not update
            minWidget.setValue(min(clippedVals)) 
            maxWidget.setValue(max(clippedVals))
        # update the value input widget itself
        newValue = self.get(key)
        if newValue is not None:
            # clipping function in bases.algorithm.parameter def
            p.setDisplayValue(newValue)
            clippedVal = p.displayValue()
            try: 
                valueWidget.setValue(clippedVal)
            except AttributeError:
                pass
        # fit parameter related updates
        newActive = self.get(key+"active")
        if isinstance(newActive, bool): # None for non-fit parameters
            # update active state for fit parameters
            p.setActive(newActive)
            if p.isActive():
                valueWidget.hide()
                minWidget.show()
                maxWidget.show()
        if not isinstance(newActive, bool) or not newActive:
            valueWidget.show()
            try:
                minWidget.hide()
                maxWidget.hide()
            except: pass
        if isNotNone(minValue, maxValue):
            # the range was updated
            self.sigRangeChanged.emit()

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
        return # FIXME
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
            if dtype is float:
                widgetType = SciEntryBox
            else:
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
        if param is None:
            return None
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addStretch()

        try:
            suffix = param.suffix()
        except AttributeError:
            suffix = None

        if suffix is None:
            layout.addWidget(self._makeLabel(param.displayName()))
        else:
            layout.addWidget(self._makeLabel(u"{} ({})".format(
                param.displayName(), suffix)))

        widget.setLayout(layout)
        if isString(param.__doc__):
            #add description as tooltip if available for parameter
            widget.setToolTip(param.__doc__)

        # create scalar value input widget with min/max limits
        minmaxValue, widgets = None, []
        if isinstance(param, ParameterNumerical):
            minmaxValue = param.min(), param.max()
        if isinstance(param, ParameterFloat):
            minmaxValue = param.displayValueRange()

        w = self._makeEntry(param.name(), param.dtype, param.displayValue(),
                                      minmax = minmaxValue, parent = widget)
        try: # set special value text for lower bound, simple solution
            special = param.displayValues(w.minimum())
            w.setSpecialValueText(special)
        except: 
            pass
        w.setFixedWidth(FIXEDWIDTH)
        widgets.insert(len(widgets)/2, w)

        # Special widget settings for active fitting parameters:
        activeBtns = activeBtns and isinstance(param, FitParameterBase)
        if activeBtns:
            # create input boxes for user specified active fitting range
            # within default upper/lower from class definition
            # minmaxValue = type(param).min(), type(param).max()
            activeRange = {"min": param.displayActiveRange()[0],
                           "max": param.displayActiveRange()[1] }
            for bound in "min", "max":
                w = self._makeEntry(param.name() + bound, param.dtype,
                                    activeRange[bound],
                                    minmax = minmaxValue, parent = widget)
                w.setPrefix(bound + ": ")
                w.setFixedWidth(FIXEDWIDTH)
                widgets.append(w)
            # create *active* buttons for FitParameters only
            w = self._makeEntry(param.name()+"active", bool,
                                param.isActive(),
                                widgetType = QPushButton,
                                parent = widget)
            w.setText("Active")
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
                    "numContribs", "findBackground")):
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

class ModelWidget(SettingsWidget):
    sigModelChanged = Signal()

    def __init__(self, *args, **kwargs):
        SettingsWidget.__init__(self, *args, **kwargs)
        self.title = TitleHandler.setup(self, "Model")

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
        # rebuild the UI without early signals
        try:
            self.sigRangeChanged.disconnect(self.sigModelChanged)
        except: sys.exc_clear()
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
        self.sigRangeChanged.connect(self.sigModelChanged)
        self.sigModelChanged.emit()

    def selectModel(self, model):
        """*model*: string containing the name of the model to select."""
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

    @property
    def algorithm(self):
        # This should really be called "model", not "algorithm" to avoid 
        # confusion
        if self.calculator is None:
            return None
        return self.calculator.model

    def setSphericalSizeRange(self, minVal, maxVal):
        key = "radius"
        # get parameter display order of magnitude: 
        param = None
        for p in self.algorithm.params():
            if p.name() == key: # "is" does not work with unicode strings
                param = p
                break
        if param is None:
            logging.debug("No 'radius'-named parameter found, "
                          "not setting spherical size range!")
            return # nothing to do
        keymin, keymax = key+"min", key+"max"
        if self.get(keymin) is not None and self.get(keymax) is not None:
            self.set(keymin, param.toDisplay(minVal))
            self.set(keymax, param.toDisplay(maxVal))

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

    def setupUi(self):
        self.listWidget.setAlternatingRowColors(True)
        setBackgroundStyleSheet(self, "./resources/background_files.svg")

class ToolBox(QToolBox):
    """QToolBox containing the widgets for user settings.
    Used to propagate resize events to child widgets to enable responsive behaviour.
    On MacOS, fixes failed detection of size changes in child widget due to scroll area.
    """
    def resizeEvent(self, event):
        child = self.currentWidget()
        if isinstance(child, AlgorithmWidget):
            child.resizeEvent(event)

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
        # put the log widget at the bottom
        self.addDockWidget(Qt.BottomDockWidgetArea, self._setupLogWidget())
        # file widget at the top
        self.toolbox = ToolBox(self)
        self._addToolboxItem(self._setupFileWidget())
        self._addToolboxItem(self._setupAlgoWidget())
        self._addToolboxItem(self._setupModelWidget())
        self._addToolboxItem(self._setupStatsWidget())

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
        self.statsWidget = RangeList(parent = self,
                                     calculator = self.calculator,
                                     title = "Post-fit Analysis",
                                     withBtn = False, nestedItems = False)
        self.modelWidget.setStatsWidget(self.statsWidget)
        self.modelWidget.sigModelChanged.connect(self.statsWidget.updateHistograms)
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
        return
        # test for additionally storing settings to file
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
        processEventLoop()
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
                                   stopFunc = self.calculator.isStopped,
                                   showProgress = False)

    def closeEvent(self, closeEvent):
        super(MainWindow, self).closeEvent(closeEvent)
        self.onStartStopClick(False)
        self.onCloseSignal.emit()

    def keyPressEvent(self, keyEvent):
        if keyEvent.key() == Qt.Key_Escape and self.startStopBtn.isChecked():
            self.onStartStopClick(False) # hit 'stop'
            logging.info("Calculation aborted by user interupt!")

# vim: set ts=4 sts=4 sw=4 tw=0:
