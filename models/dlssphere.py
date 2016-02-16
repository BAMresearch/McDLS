# -*- coding: utf-8 -*-
# models/dlssphere.py

import numpy
from numpy import pi, exp, sqrt, sin, cos
from bases.algorithm import RandomUniform
from utils.parameter import FitParameter
from scatteringmodel import DLSModel
from utils.units import NM

class DLSSphere(DLSModel):
    """Sphere Form factor for DLS"""
    shortName = "Sphere (DLS)"
    parameters = (FitParameter("radius",
                    NM.toSi(100.), unit = NM,
                    displayName = "Hydrodynamic Radius",
                    valueRange = (0., numpy.inf),
                    activeRange = NM.toSi((1., 1000.)),
                    generator = RandomUniform,
                    decimals = 1), 
                  )

    def __init__(self):
        super(DLSSphere, self).__init__()
        self.radius.setActive(True)

    def volume(self):
        return (pi*4./3.) * self.radius()**3

    def _ffSphere(self):
        qr = self.angles * self.radius()
        ff = 3. * (sin(qr) - qr * cos(qr)) / (qr**3.) # usual sphere ff
        #return 1.0
        return ff

    def weight(self):
        # a compExp. < 1 reduces the volume contribution to the amplitude
        # compExp. << 1 (e.g. 1e-4): equal volume contrib. for all radii
        return self._ffSphere() * self.volume()**self.compensationExponent

    def formfactor(self, data):
        res = exp( data.tauGamma.sanitized / self.radius() )
        # in case of multi-angles, helpers from MultiDataVector
        # res = data.tauGamma.unflatten(res) * ff
        # res = data.tauGamma.flatten(res)
        return res

DLSSphere.factory()

def test():
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
