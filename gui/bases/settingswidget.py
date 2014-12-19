# -*- coding: utf-8 -*-
# gui/bases/settingswidget.py

from __future__ import absolute_import # PEP328
from gui.qt import QtCore, QtGui
from QtGui import (QWidget, QSpinBox, QDoubleSpinBox, QLineEdit, QCheckBox,
                   QAbstractButton, QAbstractSpinBox, QLineEdit)
from QtCore import QSignalMapper, QObject
from gui.utils.signal import Signal
from utils import isList, isString

class SettingsWidget(QWidget):
    """
    Provides access to user provided application settings.

    Call get('<objectname>') to get an input widget value.
    """
    sigValuesChanged = Signal()       # signals a change of any child widget
    sigValueChanged = Signal(QWidget) # signals a change of a specific widget
    _signalMapper = None
    _isUpdateRequired = None
    _inputWidget = { int:   QSpinBox,
                     float: QDoubleSpinBox,
                     str:   QLineEdit,
                     bool:  QCheckBox,
                   }

    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self._signalMapper = QSignalMapper()
        self._signalMapper.mapped[QObject].connect(self._editingFinishedSlot)
        try:
            self.setupUi(self)
        except AttributeError:
            pass
        self._isUpdateRequired = False
        self._connect(self)

    def getInputWidget(self, datatype):
        return self._inputWidget.get(datatype, QLineEdit)

    def get(self, key, defaultValue = None):
        child = self.findChild(QWidget, key)
        value = defaultValue
        if child is None:
            return value
        propValue = None
        for attr in "value", "checked", "text":
            try:
                valueGetter = getattr(child, attr)
            except:
                value = child.property(attr)
            else:
                try: # float conversion may fail
                    value = valueGetter()
                except ValueError:
                    value = None
                    # stop here and return None,
                    # otherwise another getter is tried which may return text
                    break
            if value is not None:
                break
        return value

    def set(self, key, value):
        if value is None:
            return
        child = self.findChild(QWidget, key)
        if child is None:
            return
        setter, dtype = None, None
        for attr in "checked", "value", "text":
            setterName = "set" + attr.title()
            try:
                setter = getattr(child, setterName)
                try:
                    oldValue = getattr(child, attr)()
                except:
                    oldValue = child.property(attr)
                dtype = type(oldValue)
            except AttributeError:
                continue
            else:
                if dtype is not None:
                    break
        try:
            if (issubclass(dtype, bool) and
                isString(value) and
                value.lower() == "false"):
                value = ""
            setter(dtype(value))
        except (TypeError, ValueError):
            raise IndexError("Could not set widget '{0}' to '{1}'!".
                             format(key, value))
        else:
            self.getEditingFinishedSignal(child).emit()
            child.update()

    def _connect(self, widget):
        self._connectInputWidget(widget)
        for child in widget.children():
            self._connect(child)

    def _connectInputWidget(self, widget):
        def isInputWidget(widget):
            return any([isinstance(widget, widgetType)
                        for widgetType in QAbstractButton,
                                          QAbstractSpinBox,
                                          QLineEdit
                        ])
        if not isInputWidget(widget):
            return
        self._signalMapper.setMapping(widget, widget)
        # connect appropriate signals, depends on widget type
        self.getEditingFinishedSignal(widget).connect(self._signalMapper.map)
        # value changed signals
        try: # QAbstractButton
            widget.clicked[bool].connect(self._valueChangedSlot)
        except AttributeError: pass
        try: # QCheckbox
            widget.stateChanged.connect(self._valueChangedSlot)
        except AttributeError: pass
        try: # QLineEdit
            widget.textChanged.connect(self._valueChangedSlot)
        except AttributeError: pass
        try: # QSpinbox
            widget.valueChanged.connect(self._valueChangedSlot)
        except AttributeError: pass

    connectInputWidgets = _connectInputWidget

    @staticmethod
    def getEditingFinishedSignal(widget):
        if not isinstance(widget, QWidget):
            return
        sig = None
        for sigName in "editingFinished", "clicked":
            sig = getattr(widget, sigName, None)
            if sig is not None:
                break
        return sig

    def _editingFinishedSlot(self, widget):
        """Called after input widget looses focus,
        usually only one at a time."""
        self._isUpdateRequired = False
        self.sigValuesChanged.emit()
        self.sigValueChanged.emit(widget)

    def _valueChangedSlot(self, dummy):
        self._isUpdateRequired = True

# vim: set sts=4 ts=4 sw=4 tw=0:
