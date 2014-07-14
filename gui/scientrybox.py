# -*- coding: utf-8 -*-
# mainwindow.py

from QtGui import QLineEdit, QDoubleValidator
from QtCore import Qt

class SciEntryBox(QLineEdit):
    def __init__(self, parent):
        QLineEdit.__init__(self, parent)

        lval = QDoubleValidator()
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
        
    def setValue(self, value):
        fstr = "{:." + str(self.decimals()) + "g}" 
        self.setText(unicode(fstr.format(value)))

    def value(self):
        return float(self.text())

    def setPrefix(self, value):
        self.setPlaceholderText(value)

# vim: set ts=4 sts=4 sw=4 tw=0:
