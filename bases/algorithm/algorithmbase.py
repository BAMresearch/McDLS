# -*- coding: utf-8 -*-
# bases/algorithm/algorithmbase.py

from __future__ import absolute_import # PEP328
from utils import isString, isList, testfor, assertName
from utils.mixedmethod import mixedmethod
from utils.classproperty import classproperty
from bases.algorithm.parameter import ParameterBase, ParameterError

class AlgorithmError(StandardError):
    pass

class AlgorithmNameError(AlgorithmError):
    pass

class AlgorithmParameterError(AlgorithmError):
    pass

# a) parameters should be set by providing a ParameterBase class
# b) parameters should be access by their name in both instance AND class
# Sphere.radius.setValue()
# as well as
# sph = Sphere()
# sph.radius.setValue()
# How to define an Algorithm? By python-style class definition?
# (b) requires the algorithm to be constructed programmatically, doing sth
# after class definition was loaded

# How to let the user define his model easily?
# -> needs to supply some functions, in the current case: ff() and vol()
#    -> in the filter case: input() and output(), different use
# -> the user implements these function, by accessing Spere.radius
#    -> a staticfunction or similar would be counterintuitive because
#       parameters would not be available there by logic (from implementers POV)
# -> therefore, we need some kind of subclass, does multiple inheritance help?
#    inheriting from AlgorithmBase and the user class?

# Does an AlgorithmBase class contain Parameter types or instances?
# Changing a Parameter type changes all future Parameters of that type.
# -> allows to define e.g. RadiusParameter which can be reused for several
#    models.

class AlgorithmBase(object):
    """Base class for all data filtering algorithms."""
    _name = None # name to display in GUI
    _parameters = None # list of parameters types or instances

    @classmethod
    def setName(cls, name):
        assertName(name, AlgorithmNameError)
        cls._name = name

    @classmethod
    def name(cls):
        return cls._name

    @classmethod
    def setParams(cls, *parameters):
        """Expects a list of ParameterBase classes/types and sets them as
        class attributes to this class. They will become instances later,
        please see __init__()"""
        for i, p in enumerate(parameters):
            testfor(isinstance(p, type) and issubclass(p, ParameterBase),
                    AlgorithmParameterError, "{name}: Expected a "
                    "ParameterBase for parameter {index}, got {type}!"
                    .format(name = cls.__name__, index = i, type = type(p)))
            cls.setParam(p)

    @mixedmethod
    def setParam(selforcls, p):
        """Sets the given parameter as an attribute of this object. Use this
        method instead of setattr()."""
        # sanitize the parameter name
        name = unicode(p.name())
        replacements = dict([(ord(char), None) for char in u' \t\n\r'])
        attrname = name.translate(replacements)
        # use the clean name
        if not isList(selforcls._parameters):
            selforcls._parameters = []
        old = getattr(selforcls, attrname, None)
        oldIdx = -1
        if old is not None:
            try:
                p.setOnValueUpdate(old.onValueUpdate())
            except AttributeError:
                pass
            try:
                oldIdx = selforcls._parameters.index(old)
            except ValueError:
                pass
        setattr(selforcls, attrname, p)
        # maintain a list of all parameters, replace the old by the new to
        # keep the ordering or append the new if it didn't exist already
        if oldIdx >= 0:
            selforcls._parameters[oldIdx] = p
        else:
            selforcls._parameters.append(p)

    @mixedmethod
    def params(selforcls):
        """All parameters of this algorithm."""
        return selforcls._parameters

    @mixedmethod
    def param(selforcls, index):
        return selforcls._parameters[index]

    @mixedmethod
    def paramCount(selforcls):
        return len(selforcls._parameters)

    @classmethod
    def factory(cls, name = None, *parameters):
        """This sets the algorithm up and adds the defined parameters as class
        attributes. They become instance attributes automatically at
        initialization."""
        if name is None:
            if hasattr(cls, "shortName"):
                # for backwards compatibility
                name = cls.shortName
            else:
                name = cls.__name__
        cls.setName(name)
        if not len(parameters) and hasattr(cls, "parameters"):
            # for backwards compatibility
            parameters = []
            for baseCls in reversed(cls.__mro__):
                # get parameters from parent classes as well
                try:
                    # add parameters to final list without duplicates
                    parameters.extend([p for p in baseCls.parameters
                                          if p not in parameters])
                except:
                    continue
        cls.setParams(*parameters)
        return cls

    def __init__(self):
        """Creates instances from defined parameters and replaces the class
        attributes accordingly."""
        testfor(self._name is not None, AlgorithmNameError,
                "No name provided! "
                "Set up {name} by calling factory() first."
                .format(name = type(self).__name__))
        testfor(self._parameters is not None, AlgorithmParameterError,
                "Parameters not configured! "
                "Set up {name} by calling factory() first."
                .format(name = type(self).__name__))
        paramTypes = self.params()
        self._parameters = None
        for ptype in paramTypes:
            self.setParam(ptype())

    def __str__(self):
        text = [ self.name() ]
        for i, p in enumerate(self.params()):
            text.append(u"  {0}: {1}".format(i, unicode(getattr(self, p.name()))))
        return "\n".join(text)

    @classmethod
    def makeDefault(cls):
        return cls()

    def copy(self):
        other = self.makeDefault()
        for p in self.params():
            other.setParam(p.copy())
        return other

    def __eq__(self, other):
        if (self.name() != other.name() or
            self.paramCount() != other.paramCount()):
            return False
        for p in self.params():
            if getattr(self, p.name()) != getattr(other, p.name()):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __reduce__(self):
        return (_unpickleAlgo, (type(self), self.name(), self.params(),))

# not used atm
def _unpickleAlgo(cls, name, params):
    algo = cls.makeDefault()
    for p in params:
        algo.setParam(p)
    return algo

if __name__ == "__main__":
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
