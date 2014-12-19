# -*- coding: utf-8 -*-
# gui/bases/expdoublespinbox.py

from __future__ import absolute_import # PEP328
from gui.qt import QtGui
from QtGui import QDoubleSpinBox

class ExpDoubleSpinBox(QDoubleSpinBox):
    """Double spin box with exponentially increasing value."""
    
    def stepBy(self, steps):
        value = self.value() * 10.**steps
        if value < self.minimum():
            value = self.minimum()
        if value > self.maximum():
            value = self.maximum()
        if self.value() <= 0.0:
            value = 10.**-self.decimals()
        self.setValue(value)

# vim: set ts=4 sts=4 sw=4 tw=0:
