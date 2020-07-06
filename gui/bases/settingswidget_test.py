# -*- coding: utf-8 -*-
# gui/bases/settingswidget_test.py

from QtWidgets import QSpinBox, QDoubleSpinBox, QLineEdit, QApplication
from gui.bases.settingswidget import SettingsWidget
from utils import EPS

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

# vim: set sts=4 ts=4 sw=4 tw=0:
