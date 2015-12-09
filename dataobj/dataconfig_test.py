# -*- coding: utf-8 -*-
# dataobj/dataconfig_test.py

from __future__ import absolute_import # PEP328
from functools import partial
from nose.tools import assert_raises
from dataobj.dataconfig import DataConfig

def testDefaults():
    dc = DataConfig()
    assert dc.x0Low()  == 0.0
    assert dc.x0High() >= 1e200
    assert dc.x1Low()  == 0.0
    assert dc.x1High() >= 1e200
    assert dc.fMaskZero() == False
    assert dc.fMaskNeg()  == False

def testSetter():
    dc = DataConfig()
    dc.x0Low.setValue(-1.5)
    dc.x0High.setValue(12.)
    assert dc.x0Low()  == 0
    assert dc.x0High() == 12

def testLimits():
    dc = DataConfig()
    dc.x0Low.setValue(-2.5)
    dc.x0High.setValue(1.234e2)
    dc.setX0ValueRange((0., 10.))
    assert dc.x0Low()  == 0.
    assert dc.x0High() == 10.

def testCallbacks():
    class X0CallbackRun(StandardError):
        pass
    class X1CallbackRun(StandardError):
        pass
    class FMasksCallbackRun(StandardError):
        pass
    def test(err, *args):
        raise err
    dc = DataConfig()
    dc.register("x0limits", partial(test, X0CallbackRun))
    dc.register("x1limits", partial(test, X1CallbackRun))
    dc.register("fMasks", partial(test, FMasksCallbackRun))
    dc.x0Low.setValue(-2.5) # outside value range, not updated
    assert dc.x0Low()  == 0.
    assert_raises(X0CallbackRun, dc.x0High.setValue, 2.5)
    assert dc.x0High() == 2.5
    dc.x1Low.setValue(-1.5) # outside value range, not updated
    assert dc.x1Low()  == 0.
    assert_raises(X1CallbackRun, dc.x1High.setValue, 1.5)
    assert dc.x1High() == 1.5
    assert_raises(FMasksCallbackRun, dc.fMaskZero.setValue, True)
    assert dc.fMaskZero()
    assert_raises(FMasksCallbackRun, dc.fMaskNeg.setValue, True)
    assert dc.fMaskNeg()

def testCopy():
    dc = DataConfig()
    dc.x0Low.setValue(-2.5)
    dc.x0High.setValue(1.234e2)
    dc2 = dc.copy()
    assert dc == dc2
    # test if they behave individually now
    dc.x1Low.setValue(.5)
    assert dc.x1Low() == .5
    assert dc2.x1Low() == 0
    assert dc2.x1Low() != dc.x1Low()

def testSerialize():
    def dummyFunc(*args):
        pass
    dc = DataConfig()
    dc.x0Low.setValueRange((1, 12))
    dc.x0Low.setValue(.5)
    dc.x0High.setValue(1.234e2)
    assert dc.x0Low() == 1
    assert dc.x0High() == 1.234e2
    dc.register("x0limits", dummyFunc)
    dc.register("x1limits", dummyFunc)
    dc.register("fMasks",   dummyFunc)
    import pickle
    data = pickle.dumps(dc)
    dc2 = pickle.loads(data)
    assert dc == dc2
    for key in dc.__dict__:
        if "callback" in key:
            continue
        v0, v1 = getattr(dc, key), getattr(dc2, key)
        assert v0 == v1, \
                "{key} does not match!".format(key = key)
        if v0 is None: # None is a singleton: always same id()
            continue
        assert id(v0) != id(v1), \
                "Same value Ids for {key}!".format(key = key)

# vim: set ts=4 sts=4 sw=4 tw=0:
