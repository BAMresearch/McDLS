# -*- coding: utf-8 -*-
# gui/settingswidget.py

from __future__ import absolute_import # PEP328
import logging

from gui.qt import QtCore, QtGui
from gui.utils.signal import Signal
from QtCore import Qt, QSettings, QRegExp
from QtGui import (QWidget, QHBoxLayout, QPushButton,
                   QLabel, QLayout, QMessageBox)
from gui.bases.datalist import DataList
from gui.bases.settingswidget import SettingsWidget as SettingsWidgetBase
from bases.algorithm.parameter import ParameterFloat # instance for test
from utils import isList, isString, testfor
from utils.parameter import (ParameterNumerical, FitParameterBase)
from gui.calc import Calculator

from gui.scientrybox import SciEntryBox

FIXEDWIDTH = 120

def isNotNone(lst):
    if not isList(lst):
        return False
    return all((a is not None for a in lst))

def askForAutoUpdateRanges(parentWidget, userHists):
    mb = QMessageBox(parentWidget)
    text = ["{0}ser defined analysis range{1} appear to be set for this parameter:"]
    if len(userHists) == 1:
        text[0] = text[0].format("A u", "")
    else:
        text[0] = text[0].format("U", "s")
    for h in userHists:
        text.append("    " + str(h))
    text.append("Please confirm you want to switch from automatic to manual range setting")
    mb.setText("\n".join(text))
    mb.setWindowTitle("Automatic or manual analysis range?")
    okBtn = mb.addButton("Automatic", QMessageBox.AcceptRole)
    cnlBtn = mb.addButton("Manual", QMessageBox.RejectRole)
    mb.setDefaultButton(okBtn)
    mb.setEscapeButton(cnlBtn)
    return mb.exec_() == QMessageBox.AcceptRole

class SettingsWidget(SettingsWidgetBase):
    _calculator = None # calculator instance associated
    _appSettings = None
    sigRangeChanged = Signal()
    _autoUpdateHistograms = True # Remembers last user question

    def __init__(self, parent, calculator = None):
        SettingsWidgetBase.__init__(self, parent)
        self.sigValueChanged.connect(self.updateParam)
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
    def inputWidgets(self):
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
                yield w

    @property
    def keys(self):
        for w in self.inputWidgets:
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

    def _paramFromWidget(self, widget):
        if self.algorithm is None:
            return
        p = None
        try:
            p = getattr(self.algorithm, widget.parameterName)
        except AttributeError:
            p = None
        return p

    def updateParamRange(self, param, newRange):
        if not isNotNone(newRange): # any is None
            return None
        displayRange = None
        try:
            param.setDisplayActiveRange(newRange)
            displayRange = param.displayActiveRange() # get updated values
        except:
            param.setActiveRange(newRange)
            displayRange = param.activeRange()
        # check configured analysis ranges
        userHists = [h for h in param.histograms()
                        if h.xrange != param.activeRange()]
        # update histograms automatically?
        if self._statsWidget.userInput:
            self._autoUpdateHistograms = askForAutoUpdateRanges(self, userHists)
            self._statsWidget.userInput = False
        if self._autoUpdateHistograms:
            param.histograms().updateRanges(force = True)

        return displayRange

    def updateParam(self, widget):
        """Write UI settings back to the algorithm."""
        p = self._paramFromWidget(widget)
        if p is None:
            logging.error("updateParam({}) could not find associated parameter!"
                          .format(widget.objectName()))
            return
        # disable signals during ui updates
        self.sigValueChanged.disconnect(self.updateParam)
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
        newRange = self.updateParamRange(p, (minValue, maxValue))
        if isNotNone(newRange):
            minWidget.setValue(min(newRange))
            maxWidget.setValue(max(newRange))
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
        # enable signals again after ui updates
        self.sigValueChanged.connect(self.updateParam)
        if isNotNone(newRange):
            # the range was updated
            self.sigRangeChanged.emit()

    def updateAll(self):
        """Called in MainWindow on calculation start."""
        for w in self.inputWidgets:
            self.updateParam(w)

    def setStatsWidget(self, statsWidget):
        """Sets the statistics widget to use for updating ranges."""
        assert(isinstance(statsWidget, DataList))
        self._statsWidget = statsWidget

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
            # add description as tooltip if available for parameter
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
        self.updateParam(widgets[-1])
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

# vim: set ts=4 sts=4 sw=4 tw=0: