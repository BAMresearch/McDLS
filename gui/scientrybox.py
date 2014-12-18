# -*- coding: utf-8 -*-
# mainwindow.py

from QtGui import QLineEdit, QDoubleValidator, QValidator
from QtCore import Qt
from utils import clip

class SciEntryValidator(QDoubleValidator):
    """Assumes the associated QLineEdit is provided as parent object
    in the constructor."""

    def __init__(self, parent):
        super(SciEntryValidator, self).__init__(parent)
        self.setNotation(QDoubleValidator.ScientificNotation)
        self.setRange(-1e200, 1e200)
        self.setDecimals(9)

    def validate(self, input, pos):
        state, value, pos = QDoubleValidator.validate(self, input, pos)
        if state is QValidator.State.Invalid:
            # do not accept any invalid text, usually not float compatible
            return state, value, pos
        elif state is QValidator.State.Acceptable:
            # verfiy the value range again,
            # QDoubleValidator.validate() fails sometimes on the range
            try:
                value = float(value)
                if value < self.bottom() or value > self.top():
                    raise ValueError
            except:
                state = QValidator.State.Intermediate
        # highlight the input on errorneous values
        self.parent().indicateCorrectness(state is QValidator.State.Acceptable)
        if self.parent().hasFocus(): # change nothing, by default
            return QValidator.State.Acceptable, input, pos
        # else: on focus left, fix it if required and validate again
        if state is not QValidator.State.Acceptable:
            input = self.fixup(input)
            state, value, pos = QDoubleValidator.validate(self, input, pos)
            if state is not QValidator.State.Acceptable: # keep focus
                self.parent().setFocus(Qt.OtherFocusReason)
        return state, unicode(value), pos

    def setTop(self, value):
        """Work around numerical issues for comparing numbers from user input
        with calculation results."""
        QDoubleValidator.setTop(self, float(str(value)))

    def setBottom(self, value):
        """Work around numerical issues for comparing numbers from user input
        with calculation results."""
        QDoubleValidator.setBottom(self, float(str(value)))

    def fixup(self, input):
        # restrict the value to the valid range
        try:
            input = str(clip(float(input), self.bottom(), self.top()))
        except: pass # do nothing if float conversion fails
        return input

class SciEntryBox(QLineEdit):
    def __init__(self, parent):
        QLineEdit.__init__(self, parent)

        val = SciEntryValidator(self)
        self.setValidator(val)
        self.setAlignment(Qt.AlignRight)
        self.setMaxLength(14)
        setattr(self, "setMinimum", val.setBottom)
        setattr(self, "setMaximum", val.setTop)
        setattr(self, "minimum", val.bottom)
        setattr(self, "maximum", val.top)
        setattr(self, "setDecimals", val.setDecimals)
        setattr(self, "decimals", val.decimals)

    def indicateCorrectness(self, isValid):
        if not isValid:
            self.setStyleSheet("background: rgb(255,160,160);")
        else:
            self.setStyleSheet("")

    def setRange(self, lo, hi):
        self.setMinimum(lo)
        self.setMaximum(hi)
        
    def setValue(self, value):
        fstr = "{:." + str(self.decimals()) + "g}" 
        self.setText(unicode(fstr.format(value)))

    def value(self):
        return float(self.text())

    def setPrefix(self, value):
        self.setPlaceholderText(value)

# vim: set ts=4 sts=4 sw=4 tw=0:
