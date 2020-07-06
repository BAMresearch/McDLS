# -*- coding: utf-8 -*-
# gui/mainwindow.py

from builtins import str
from QtGui import QDoubleValidator, QValidator
from QtWidgets import QLineEdit
from QtCore import Qt

from ..utils import clip

class SciEntryValidator(QDoubleValidator):
    """Assumes the associated QLineEdit is provided as parent object
    in the constructor."""

    def __init__(self, parent = None):
        # FIXME: On Windows, python crashes here, sometimes ...
        # everything looks ok, except that it 'stops working'
        super(SciEntryValidator, self).__init__(parent)
        self.setNotation(QDoubleValidator.ScientificNotation)
        # some defaults if nothing is set
        self.setRange(-1e200, 1e200)
        self.setDecimals(9)

    def validate(self, txt, pos):
        state, value, pos = QDoubleValidator.validate(self, txt, pos)
        if state is QValidator.State.Invalid:
            # do not accept any invalid text, usually not float compatible
            return state, value, pos
        # highlight the input on errorneous values, but accept them
        self.parent().indicateCorrectness(state is not QValidator.State.Invalid)
        if self.parent().hasFocus(): # change nothing, by default
            return QValidator.State.Acceptable, txt, pos
        # verify the value range in any case,
        # QDoubleValidator.validate() fails sometimes on the range
        txt = self.fixup(txt)
        state, txt, pos = QDoubleValidator.validate(self, txt, pos)
        if state is QValidator.State.Invalid: # keep focus
            self.parent().setFocus(Qt.OtherFocusReason)
        self.parent().indicateCorrectness(state is not QValidator.State.Invalid)
        if state is QValidator.State.Intermediate:
            # there is an issue with floating point format not being accepted
            # with ScientificNotation set, ignore this here since the value is clipped above
            # avoid endless loop because of this
            state = QValidator.State.Acceptable
        return state, txt, pos

    def fixup(self, txt):
        """Restricts the value to the valid range defined by setTop() and
        setBottom(). Limits the precision as well."""
        # restrict the value to the valid range
        try:
            txt = self.parent().fmt.format(
                clip(float(txt), self.bottom(), self.top()))
        except ValueError:
            pass # do nothing if float conversion fails
        return str(txt)

class SciEntryBox(QLineEdit):
    toolTipFmt = "A value between {lo} and {hi} (including)."

    def __init__(self, parent = None):
        super(SciEntryBox, self).__init__(parent)

        val = SciEntryValidator(self)
        self.setValidator(val)
        self.setAlignment(Qt.AlignRight)
        self.setMaxLength(14)
        setattr(self, "minimum", val.bottom)
        setattr(self, "maximum", val.top)
        setattr(self, "setDecimals", val.setDecimals)
        setattr(self, "decimals", val.decimals)

    @property
    def fmt(self):
        return self.numberFormat(self.decimals())

    @staticmethod
    def numberFormat(decimals = None):
        if decimals is None:
            return "{0:g}"
        # according to https://docs.python.org/2/library/string.html#formatspec
        return "{0:." + str(max(0, decimals - 3)) + "g}"

    @classmethod
    def updateToolTip(cls, widget, decimals = None):
        text = []
        try:
            parentText = widget.parent().toolTip()
            if len(parentText):
                text.append(parentText)
        except:
            pass
        fmt = cls.numberFormat(decimals)
        text.append(("<span>" + cls.toolTipFmt + "</span>").format(
                lo = fmt.format(widget.minimum()),
                hi = fmt.format(widget.maximum())))
        widget.setToolTip("<br />\n".join(text))

    def setMinimum(self, value):
        """Work around issues regarding round-off errors by the 'g' format type
        when comparing numbers from user input with calculation results."""
        self.validator().setBottom(float(self.fmt.format(value)))
        self.updateToolTip(self, self.decimals())

    def setMaximum(self, value):
        """Work around issues regarding round-off errors by the 'g' format type
        when comparing numbers from user input with calculation results."""
        self.validator().setTop(float(self.fmt.format(value)))
        self.updateToolTip(self, self.decimals())

    def indicateCorrectness(self, isValid):
        if not isValid:
            self.setStyleSheet("background: rgb(255,160,160);")
        else:
            self.setStyleSheet("")

    def setRange(self, lo, hi):
        self.setMinimum(lo)
        self.setMaximum(hi)
        
    def setValue(self, value):
        self.setText(str(self.fmt.format(value)))

    def value(self):
        return float(self.text())

    def setPrefix(self, value):
        self.setPlaceholderText(value)

# vim: set ts=4 sts=4 sw=4 tw=0:
