# -*- coding: utf-8 -*-
# algorithmbase.py

from classproperty import classproperty
from utils import isString, isList
from parameter import Parameter, ParameterError

class AlgorithmBase(object):
    """Base class for all data filtering algorithms."""
    shortName = None # name to display in GUI
    parameters = None # list of parameters

    @classproperty
    @classmethod
    def name(cls):
        return cls.shortName

    @classmethod
    def makeDefault(cls):
        return cls()

    @classmethod
    def assertions(cls):
        for attrname in ('shortName', ):
            assert (isString(getattr(cls, attrname)) and
                    len(getattr(cls, attrname)) > 0), \
                "{0} attribute '{1}' is not set!".format(cls.name, attrname)
        assert AlgorithmBase.parameters is None, \
                "Do not add parameters to the base class!"
        assert isList(cls.parameters), \
                ("Please set {0} attribute 'parameters' to an empty tuple"
                 " == () if no parameters needed!".format(cls.name))
        for p in cls.parameters: # verify parameters
            assert isinstance(p, Parameter)

    def copy(self):
        f = self.makeDefault()
        for p in self.parameters:
            setattr(f, p.name, getattr(self, p.name))
        return f

    def __init__(self, settings = None):
        self.assertions()
        # init with default values
        for p in self.parameters:
            assert not hasattr(self, p.name), \
                    "Parameter '{0}' already set!".format(p.name)
            setattr(self, p.name, p.defaultValue)

    def __str__(self):
        text = [self.name]
        for p in self.parameters:
            text.append("{0}: '{1}'".format(p.name, getattr(self, p.name)))
        return "; ".join(text)

    def __eq__(self, other):
        if (self.name != other.name or
            self.parameters != other.parameters):
            return False
        for p in self.parameters:
            if getattr(self, p.name) != getattr(other, p.name):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

# vim: set ts=4 sts=4 sw=4 tw=0:
