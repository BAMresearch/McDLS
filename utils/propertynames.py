# -*- coding: utf-8 -*-
# utils/propertynames.py

import inspect

class PropertyNames(object):
    _cache = None

    @classmethod
    def properties(cls):
        """Returns all attributes configured in this class."""
        nameList = dir(cls)
        hashValue = hash(repr(nameList))
        if not cls._cache or cls._cache[0] != hashValue:
            result = [(name, getattr(cls, name)) for name in nameList
                      if not name.startswith("_") and
                      not inspect.ismethod(getattr(cls, name))]
            cls._cache = hashValue, result
        return cls._cache[1]

    @classmethod
    def propNames(cls):
        return [propName for propName, dummy in cls.properties()]

# vim: set ts=4 sts=4 sw=4 tw=0:
