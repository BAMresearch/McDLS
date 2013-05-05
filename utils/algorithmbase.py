# -*- coding: utf-8 -*-
# algorithmbase.py

from classproperty import classproperty
from utils import isString, isList
from parameter import Parameter, ParameterError, testfor

class AlgorithmError(StandardError):
    pass

class AlgorithmNameError(AlgorithmError):
    pass

class AlgorithmParameterError(AlgorithmError):
    pass

class AlgorithmBase(object):
    """Base class for all data filtering algorithms."""
    shortName = None # name to display in GUI
    parameters = None # list of parameter classes

    @classproperty
    @classmethod
    def name(cls):
        return cls.shortName

    @classmethod
    def makeDefault(cls):
        return cls()

    @classmethod
    def __len__(cls):
        if isList(cls.parameters):
            return len(cls.parameters)
        return 0

    def copy(self):
        f = self.makeDefault()
        for p in self.parameters:
            setattr(f, p.name, getattr(self, p.name))
        return f

    def __init__(self):
        testfor(isString(self.shortName) and len(self.shortName) > 0,
                AlgorithmNameError,
                "{0} name not set!".format(type(self).__name__))
        testfor(AlgorithmBase.parameters is None, AlgorithmParameterError,
                "Do not add parameters to the base class!")
        # fix parameters attributes if set to single class without list
        if self.parameters is None:
            setattr(type(self), 'parameters', ())
        elif not isList(self.parameters):
            setattr(type(self), 'parameters', (self.parameters, ))
        # verify parameter types
        for p in self.parameters:
            testfor(isinstance(p, type) and issubclass(p, Parameter),
                    AlgorithmParameterError,
                    "Expected a subclass of Parameter, got '{}'!"
                    .format(type(p)))
        # init parameter instances
        for p in self.parameters:
            testfor(not hasattr(self, p.name), AlgorithmParameterError,
                    "Parameter '{0}' already set!".format(p.name))
            setattr(self, p.name, p())

    def __str__(self):
        text = [ self.name ]
        for p in self.parameters:
            text.append("  " + str(getattr(self, p.name)))
        return "\n".join(text)

    @property
    def paramIter(self):
        """Iterator over all Parameter instances."""
        for p in self.parameters:
            yield getattr(self, p.name)

# vim: set ts=4 sts=4 sw=4 tw=0:
