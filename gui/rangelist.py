# -*- coding: utf-8 -*-
# gui/rangelist.py

from __future__ import absolute_import # PEP328
from builtins import zip
import logging

from numpy import inf as numpy_inf
from gui.qt import QtCore, QtGui
from QtCore import Qt, QFileInfo, QMargins
from QtGui import (QWidget, QHBoxLayout, QVBoxLayout, QPushButton,
                   QLabel, QComboBox, QPalette, QDialog, QSpinBox,
                   QCheckBox)
from gui.bases.datalist import DataList
from utils import isList, testfor
from utils.parameter import (ParameterNumerical, Histogram)
from gui.calc import Calculator

# do not remove, dialog will not work without this
from gui.scientrybox import SciEntryBox

from models.scatteringmodel import ScatteringModel

# required for svg graphics support
from gui.qt import QtSvg, QtXml, pluginDirs
from gui.liststyle import setBackgroundStyleSheet
from main import makeAbsolutePath

class RangeDialog(QDialog):
    """Creates a modal dialog window to ask the user for a range to be
    added."""
    _model = None
    _labels = None

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
        inputWidgets = (self._createParamBox(), self._createAutoRange(),
                        self._createLower(), self._createUpper(),
                        self._createBins(),
                        self._createXScale(), self._createYWeight())
        self._labels = dict()
        # assumes same ordering of entryWidgets above and Histogram.displayData
        for col, inputWidget in zip(Histogram.displayData, inputWidgets):
            fieldWidget = QWidget(self) # combines label + input
            fieldLayout = QVBoxLayout(fieldWidget)
            fieldLayout.setContentsMargins(QMargins())
            # create label, text is set in _selectParam()
            self._labels[col] = QLabel(self)
            self._labels[col].setAlignment(Qt.AlignHCenter)
            # stack label + input
            fieldLayout.addWidget(self._labels[col])
            fieldLayout.addWidget(inputWidget)
            fieldWidget.setLayout(fieldLayout)
            # add field to row of inputs
            entryLayout.addWidget(fieldWidget)
        entryWidget.setLayout(entryLayout)
        self.pbox.setCurrentIndex(0)
        self.lentry.selectAll() # select the first input by default
        return entryWidget

    def _createAutoRange(self):
        autoRange = QCheckBox(self)
        autoRange.stateChanged.connect(self._onAutoRangeChange)
        self.autoRange = autoRange
        return autoRange

    def _onAutoRangeChange(self, isChecked):
        isChecked = (isChecked == Qt.Checked) # no tristate
        enableLimits = not isChecked
        self.lentry.setEnabled(enableLimits)
        self.uentry.setEnabled(enableLimits)
        self._setLimits()

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

    def _getParam(self, index):
        pname = self.pbox.itemData(index)
        p = getattr(self._model, pname)
        # perhaps, use testfor() for that:
        assert p is not None, "Could not find parameter from selection box" 
        return p

    def _setLimits(self):
        p = self._getParam(self.pbox.currentIndex())
        lval, uval = p.activeRange()
        try:
            # account for units conversion
            lval, uval = p.toDisplay(lval), p.toDisplay(uval)
        except AttributeError:
            # raise # will break here, use it for debugging
            pass
        self.lentry.setValue(lval)
        self.uentry.setValue(uval)

    def _selectParam(self, index):
        """Gets the index within the selection box."""
        p = self._getParam(index)
        llim, ulim = p.valueRange()
        try:
            # account for units conversion
            llim, ulim = p.toDisplay(llim), p.toDisplay(ulim)
        except AttributeError:
            # raise # will break here, use it for debugging
            pass
        self.lentry.setRange(llim, ulim)
        self.uentry.setRange(llim, ulim)
        self.autoRange.setCheckState(Qt.Checked) # sets limits too
        # set labels
        for col, text in zip(Histogram.displayData, Histogram.displayDataDescr):
            if "lower" in col or "upper" in col:
                text = u"{t} ({u})".format(t = text, u = p.displayMagnitudeName())
            elif any([(l in text) for l in ("axis", "bins", "range")]):
                # break long descriptions to keep them short
                text = text.replace(" ", "\n", 1)
            self._labels[col].setText(text)

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
        lval, uval = (self.lentry.value(), self.uentry.value())
        try:
            # take units into account,
            # convert from display units to SI units for internal use
            lval, uval = (p.toSi(lval), p.toSi(uval))
        except AttributeError:
            pass

        hist = Histogram(p, lval, uval,
                self.bentry.value(), self.sentry.currentText(),
                self.wentry.currentText())
        hist.autoFollow = self.autoRange.isChecked()
        return hist

class RangeList(DataList):
    _calculator = None

    def __init__(self, calculator = None, **kwargs):
        DataList.__init__(self, **kwargs)
        assert isinstance(calculator, Calculator)
        self._calculator = calculator
        self.sigRemovedData.connect(self.onRemoval)
        self.listWidget.mouseDoubleClickEvent = self.mouseDoubleClickEvent

    def mouseDoubleClickEvent(self, event):
        """Shows RangeDialog on double click."""
        self.loadData()

    def itemUpdate(self, item, column):
        # handle autoFollow only, i.e. auto range update changes only
        if column != Histogram.displayData.index("autoFollow"):
            return
        # update histogram based on auto range update on/off
        hist = item.data()
        hist.autoFollow = (item.checkState(column) == Qt.Checked)
        if hist.autoFollow:
            hist.updateRange()
            item.update()

    def onRemoval(self, removedHistograms):
        for hist in removedHistograms:
            try:
                hist.param.histograms().remove(hist)
            except AttributeError:
                continue

    def updateHistograms(self):
        """Called after UI update by sigRangeChanged and/or sigModelChanged in
        MainWindow."""
        lst = []
        for p in self._calculator.model.activeParams():
            lst += [h for h in p.histograms()]
        if set(self.data()) == set(lst):
            return # no update if same hists exist already, order may differ
        self.clear()
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
        self.clearSelection()
        self.setHeader(Histogram.displayDataDescr)
        self.setToolTip(
            "Right-click to add additional ranges."
        )
        # self.addMenuEntry(name = "edit", text = "Edit selected", 
        #                   menuStates = "hasSelection",
        #                   callbacks = self.editEntry)
        # self.addMenuEntry(name = "recalc", text = "recalc histograms",
        #                   menuStates = "hasSelection",
        #                   callbacks = self.recalc)

    def editEntry(self):
        return
        # more on the TODO list... Will figure this out later.
    
    def recalc(self):
        return
        # this is not the _calculator instance I was looking for. 
        # self._calculator(recalc = True)
        # does not work yet! missing arguments for Hist.calc()
        DataList.updateData(self, selectedOnly = True, showProgress = False,
                updateFunc = Histogram.calc,
                stopFunc = None)

# vim: set ts=4 sts=4 sw=4 tw=0:
