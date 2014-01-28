# -*- coding: utf-8 -*-
# utils/parameter.py

from cutesnake.utils import mixedmethod
from cutesnake.algorithm import (ParameterBase, ParameterFloat,
                                 ParameterNumerical, ParameterBoolean,
                                 ParameterLog, ParameterString)
from cutesnake.algorithm import Parameter as ParameterFactory

class FitParameterBase(ParameterBase):
    ParameterBase.setAttributes(locals(), isActive = False)

    @mixedmethod
    def attributes(selforcls):
        res = super(FitParameterBase, selforcls).attributes()
        return res

    @classmethod
    def factory(cls, **kwargs):
        cls = super(FitParameterBase, cls).factory(**kwargs)
        return cls

class FitParameterString(FitParameterBase, ParameterString):
    pass

class FitParameterBoolean(FitParameterBase, ParameterBoolean):
    pass

class FitParameterNumerical(FitParameterBase, ParameterNumerical):
    pass

class FitParameterFloat(FitParameterBase, ParameterFloat):
    pass

class FitParameterLog(FitParameterBase, ParameterLog):
    pass

def Parameter(*args, **kwargs):
    return ParameterFactory(*args, paramTypes = (
                FitParameterBoolean,
                FitParameterFloat,
                FitParameterNumerical,
                FitParameterString,
                FitParameterBase
            ), **kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
