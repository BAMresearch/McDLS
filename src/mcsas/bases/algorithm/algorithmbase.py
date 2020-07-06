# -*- coding: utf-8 -*-
# bases/algorithm/algorithmbase.py

from builtins import str
from builtins import object
from ...utils import isString, isList, testfor, assertName
from ...utils.mixedmethod import mixedmethod
from ...utils import classproperty, classname
from .parameter import ParameterBase, ParameterError

class AlgorithmError(Exception):
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
        cls._parameters = []
        for i, p in enumerate(parameters):
            testfor(isinstance(p, type) and issubclass(p, ParameterBase),
                    AlgorithmParameterError, "{name}: Expected a "
                    "ParameterBase type for parameter {index}, got {type}!"
                    .format(name = cls.__name__, index = i, type = p))
            cls.setParam(p)

    @mixedmethod
    def setParam(selforcls, p):
        """Sets the given parameter as an attribute of this object. Use this
        method instead of setattr()."""
        attrname = p.name()
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
        if selforcls._parameters is None:
            return []
        return selforcls._parameters

    @mixedmethod
    def param(selforcls, index):
        return selforcls._parameters[index]

    @mixedmethod
    def paramCount(selforcls):
        return len(selforcls._parameters)

    @property
    def showParams(self):
        """A list of parameter names which defines the parameters and their
        ordering shown in a UI. To be overridden in sub classes."""
        return [p.name() for p in self.params()]

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
            for baseCls in reversed(cls.mro()):
                # get parameters from parent classes as well
                params = []
                try:
                    params = baseCls.params()
                except AttributeError:
                    pass
                try:
                    params += [p for p in baseCls.parameters
                                 if p not in params]
                except AttributeError:
                    pass
                # add parameters to final list without duplicates
                parameters.extend([p for p in params if p not in parameters])
            # cls.parameters = None # enforce usage of params()
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
            text.append(u"  {0}: {1}".format(i, str(getattr(self, p.name()))))
        return "\n".join(text)

    @classmethod
    def makeDefault(cls):
        return cls()

    def __eq__(self, other):
        if (not isinstance(other, type(self)) or
            self.name() != other.name() or
            self.paramCount() != other.paramCount()):
            return False
        for p in self.params():
            if getattr(self, p.name()) != getattr(other, p.name()):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getstate__(self):
        """Prepares the internal state for pickle by removing the attribute
        for each parameter which is usually set by __init__() and during
        unpickle in __setstate__()."""
        state = self.__dict__.copy()
        keys = set(state.keys())
        for cls in type(self).mro():
            if (issubclass(cls, AlgorithmBase) # prevents recursion loop
                or not hasattr(cls, "__getstate__")): # excludes other mixins
                continue
            parentState = cls.__getstate__(self)
            keys = keys & set(parentState.keys()) # sync common keys only
            state.update([(key, parentState[key]) for key in keys])
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        # self.params() contains parameter instances already at this point
        # AlgorithmBase.__init__ expects types to be initialized
        parameters = self.params()
        self._parameters = [type(p) for p in parameters]
        self.__init__() # call custom constructor code
        for p in parameters:
            newParam = getattr(self, p.name(), None)
            if not isinstance(newParam, ParameterBase):
                continue
            newParam.setAttributes(**p.attributes())

    def hdfWrite(self, hdf):
        hdf.writeAttributes(name = self.name(), cls = classname(self))
        hdf.writeMembers(self, *[p.name() for p in self.params()])

    def update(self, other):
        """Copy parameter values from another algorithm of the same type."""
        if not isinstance(other, type(self)):
            return
        for p in self.params():
            otherP = getattr(other, p.name())
            p.setValue(otherP.value()) # possibly runs value callback

if __name__ == "__main__":
    pass

# vim: set ts=4 sts=4 sw=4 tw=0:
