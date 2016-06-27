# -*- coding: utf-8 -*-
# dataobj/sasconfig_test.py

from __future__ import absolute_import # PEP328
from functools import partial
from nose.tools import assert_raises
from dataobj.sasconfig import SASConfig

def testDefaults(): # from dataconfig_test
    sc = SASConfig()
    assert sc.x0Low()  == 0.0
    assert sc.x0High() >= 1e200
    assert sc.x1Low()  == 0.0
    assert sc.x1High() >= 1e200
    assert sc.fMaskZero() == False
    assert sc.fMaskNeg()  == False

def testSerialize():
    def dummyFunc(*args):
        pass
    sc = SASConfig()
    sc.x0Low.setValueRange((1, 12))
    sc.x0Low.setValue(.5)
    sc.x0High.setValue(1.234e2)
    assert sc.x0Low() == 1
    assert sc.x0High() == 1.234e2
    sc.register("x0limits", dummyFunc)
    sc.register("x1limits", dummyFunc)
    sc.register("fMasks",   dummyFunc)
    import pickle
    data = pickle.dumps(sc)
    sc2 = pickle.loads(data)
    assert sc == sc2
    for key in sc.__dict__:
        if "callback" in key:
            continue
        v0, v1 = getattr(sc, key), getattr(sc2, key)
        assert v0 == v1, \
                "{key} does not match!".format(key = key)
        if v0 is None: # None is a singleton: always same id()
            continue
        assert id(v0) != id(v1), \
                "Same value Ids for {key}!".format(key = key)

# vim: set ts=4 sts=4 sw=4 tw=0:
