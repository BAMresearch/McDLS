# -*- coding: utf-8 -*-
# __init__.py

from algorithmbase import AlgorithmBase
from numbergenerator import NumberGenerator, RandomUniform, RandomExponential
from parameter import (ParameterBase, ParameterFloat, ParameterNumerical,
                       ParameterBoolean, ParameterLog, ParameterNameError)
from parameter import factory as makeParameterType
from cutesnake.utils.tests import assertName

def Parameter(name, value, **kwargs):
    """User interface to parameter.factory()"""
    kwargs.update(dict(name = name, value = value))
    assertName(name, ParameterNameError)
    return makeParameterType(**kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
