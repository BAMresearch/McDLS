# -*- coding: utf-8 -*-
# mainwindow.py

from QtGui import QLineEdit, QDoubleValidator
from QtCore import Qt
import sys

class MyValidator(QDoubleValidator):
    def validate(self, input, pos):
        result = QDoubleValidator.validate(self, input, pos)
        print >>sys.__stderr__, "validate '{}'".format(input), pos, result
        print >>sys.__stderr__, "validator:", self.bottom(), self.top(), self.decimals(), self.notation()
        return result

    def fixup(self, input):
        result = QDoubleValidator.fixup(self, input)
        print >>sys.__stderr__, "fixup '{}' -> '{}'".format(input, result)
        return result

class SciEntryBox(QLineEdit):
    def __init__(self, parent):
        QLineEdit.__init__(self, parent)

        lval = MyValidator()
        lval.setNotation(QDoubleValidator.ScientificNotation)
        lval.setRange(-1e100, 1e100)
        lval.setDecimals(9)
        self.setValidator(lval)
        self.setAlignment(Qt.AlignRight)
        self.setMaxLength(14)
        setattr(self, "setMinimum", lval.setBottom)
        setattr(self, "setMaximum", lval.setTop)
        setattr(self, "minimum", lval.bottom)
        setattr(self, "maximum", lval.top)
        setattr(self, "setDecimals", lval.setDecimals)
        setattr(self, "decimals", lval.decimals)

    def focusOutEvent(self, event):
        print >>sys.__stderr__, "focusOutEvent, loosing focus"
        QLineEdit.focusOutEvent(self, event)

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
