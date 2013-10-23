# -*- coding: utf-8 -*-
# settingswidget.py

from cutesnake.qt import QtCore, QtGui
try: # remove if switching to PySide exclusively
    from QtCore import QVariant
except ImportError:
    class QVariant: # a dummy class for compatibility
        Int, Double, Bool, String = None, None, None, None
        toInt, toDouble, toBool, toString = None, None, None, None

from QtGui import QWidget, QSpinBox, QDoubleSpinBox, QLineEdit
from cutesnake.widgets.expdoublespinbox import ExpDoubleSpinBox
from cutesnake.utils.signal import Signal
from cutesnake.utils import isList

class SettingsWidget(QWidget):
    """
    Provides access to user provided application settings.

    Call get('<objectname>') to get an input widget value.
    """
    sigValuesChanged = Signal()
    _isUpdateRequired = None
    _inputWidget = { int:   QSpinBox,
                     float: QDoubleSpinBox,
                     'exp': ExpDoubleSpinBox,
                     str:   QLineEdit }
    _convertData = { QVariant.Int:    QVariant.toInt,
                     QVariant.Double: QVariant.toDouble,
                     QVariant.Bool:   QVariant.toBool,
                     QVariant.String:
                        (lambda variant: str(variant.toString())) }

    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
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
        if child == None:
            return defaultValue
        variant = None
        for key in "value", "checked", "text":
            variant = child.property(key)
            try:
                if variant.isValid():
                    break
            except AttributeError: # QVariant
                if variant is not None:
                    break
        if variant is None:
            return defaultValue
        if isinstance(variant, QVariant):
            if not variant.isValid():
                return defaultValue
            variant = self._convertData[variant.type()](variant)
        value = variant
        if value is None:
            return defaultValue
        # QVariant.toInt, QVariant.toDouble return tuple (value, ok)
        if isList(value) and len(value) == 2 and isinstance(value[1], bool):
            if not value[1]:
                return defaultValue
            value = value[0]
        return value

    def set(self, key, value):
        child = self.findChild(QWidget, key)
        for key in ("isChecked", "value", "text"):
            variant = child.property(key)
            try:
                if variant.isValid():
                    break
            except AttributeError: # QVariant
                if variant is not None:
                    break
        try:
            child.setProperty(key, value)
        except StandardError:
            raise IndexError("Could not set widget '{0}' to '{1}'!".
                             format(key, value))
        child.update()

    def _connect(self, widget):
        self._connectInputWidget(widget)
        for child in widget.children():
            self._connect(child)

    def _connectInputWidget(self, widget):
        def isInputWidget(widget):
            return any([isinstance(widget, widgetType)
                        for widgetType in self._inputWidget.values()])
        if isInputWidget(widget) and hasattr(widget, "isReadOnly"):
            widget.editingFinished.connect(self._editingFinishedSlot)
            if isinstance(widget, QLineEdit): # no valueChanged signal
                widget.textChanged.connect(self._valueChangedSlot)
            else:
                widget.valueChanged.connect(self._valueChangedSlot)

    connectInputWidgets = _connectInputWidget

    def _editingFinishedSlot(self):
        if not self._isUpdateRequired:
            return
        self._isUpdateRequired = False
        self.sigValuesChanged.emit()

    def _valueChangedSlot(self, dummy):
        self._isUpdateRequired = True

# vim: set sts=4 ts=4 sw=4 tw=0:
