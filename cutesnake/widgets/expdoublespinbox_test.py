# -*- coding: utf-8 -*-
# expdoublespinbox_test.py

from expdoublespinbox import ExpDoubleSpinBox
from cutesnake.utilsgui import DialogInteraction

def test():
    edsp = DialogInteraction.instance(ExpDoubleSpinBox)
    edsp.setValue(0)
    edsp.stepBy(0)
    oldValue = edsp.value()
    assert oldValue > 0

    edsp.stepBy(1)
    oldDiff = abs(oldValue - edsp.value())
    assert oldValue < edsp.value()

    edsp.stepBy(1)
    diff = abs(oldValue - edsp.value())
    assert oldDiff < diff
    assert oldValue < edsp.value()

    edsp.stepBy(-2)
    assert oldValue == edsp.value()

# vim: set ts=4 sts=4 sw=4 tw=0:
