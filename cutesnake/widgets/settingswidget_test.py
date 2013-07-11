# -*- coding: utf-8 -*-
# settingswidget.py

from cutesnake.qt import QtGui
from QtGui import QSpinBox, QDoubleSpinBox, QLineEdit, QApplication
from cutesnake.widgets.expdoublespinbox import ExpDoubleSpinBox
from settingswidget import SettingsWidget
from cutesnake.utils import EPS

class TestSettings(SettingsWidget):
    def setupUi(self, dummy):
        pass

app = QApplication([])
w = TestSettings()
 
def testIntegerInputBox():
    sp = QSpinBox(w)
    sp.setObjectName("sp")
    sp.setValue(42)
    assert w.get("sp") == 42

def testTextualInputBox():
    le = QLineEdit(w)
    le.setObjectName("le")
    text = "c633296c74bb1e0268cf0df7f2cc50c26589bc3e"
    le.setText(text)
    assert w.get("le") == text

def testFloatingPointInputBox():
    dsp = QDoubleSpinBox(w)
    dsp.setObjectName("dsp")
    testvalue = 23.424242
    dsp.setValue(testvalue)
    assert abs(w.get("dsp") - round(testvalue, dsp.decimals())) < EPS

def testExponentialFloatingPointBox():
    edsp = ExpDoubleSpinBox(w)
    edsp.setObjectName("edsp")
    testvalue = 42.23235
    edsp.setValue(testvalue)
    assert abs(w.get("edsp") - round(testvalue, edsp.decimals())) < EPS

# vim: set sts=4 ts=4 sw=4 tw=0:
