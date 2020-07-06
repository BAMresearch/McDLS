# -*- coding: utf-8 -*-
# bases/model/scatteringmodel_test.py

from . import SASModel
from ...utils.parameter import FitParameter

class DummyModel(SASModel):
    shortName = "DummyModelName"
    parameters = (
        FitParameter(name = "testp0", value = 3.4, valueRange = (.5, 6.7)),
        FitParameter(name = "testp1", value = 4, valueRange = (1, 5))
    )
    def calcIntensity(*args): pass
    def volume(*args): pass
    def formfactor(*args): pass

# TODO: does not work for active parameters, see McSAS.plot()
#       -> serialization support for histograms/fitparameters has to be added
def testSerialization():
    import pickle
    a = DummyModel.factory()()
    data = pickle.dumps(a)
    a2 = pickle.loads(data)
    print a
    print a2
    assert a == a2

# vim: set ts=4 sts=4 sw=4 tw=0:
