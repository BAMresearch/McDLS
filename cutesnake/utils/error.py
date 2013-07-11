# -*- coding: utf-8 -*-
# error.py

"""
Some Error classes.
"""

class AppError(StandardError):
    msg = ""
    """
    General error with descriptive message to be forwarded and shown to
    the user in a message box.
    """
    def __init__(self, msg = ""):
        StandardError.__init__(self, self.getMessage(msg))

    @classmethod
    def getMessage(cls, msg = ""):
        if len(msg) == 0:
            return cls.msg
        return msg

class InitError(AppError):
    pass

class EmptySelection(AppError):
    pass

class FileError(AppError):
    def __init__(self, msg, fn):
        AppError.__init__(self,
            "{0}\n\n'{1}'".format(msg, fn))

class LoadError(AppError):
    pass

class VerboseError(AppError):
    """Verbose Exception containing class and method where it was raised.
    A cache makes sure VerboseErrors with the same name can be matched."""
    _cache = None

    def __init__(self):
        AppError.__init__(self, self.prefix(self.msg))

    @classmethod
    def prefix(cls, msg):
        """Retrieves some meta info to decribe the error further."""
        frame = inspect.stack()[2]
        inClass = frame[0].f_locals.get('self', None)
        if inClass is None:
            inClass = frame[0].f_locals.get('cls', type(None))
        else:
            inClass = type(inClass)
        className = inClass.__name__
        if className == "NoneType":
            className = ""
        methodName = frame[3]
        return "{0}.{1}: {2}".format(className, methodName, msg)

    @classmethod
    def new(cls, name, msg = ""):
        """Creates an exception dynamically.
        Remembers previous ones for testing"""
        if cls._cache is None:
            cls._cache = dict()
        if name not in cls._cache:
            cls._cache[name] = type(name, (cls,), {'msg': msg})
        return cls._cache[name]

# vim: set ts=4 sts=4 sw=4 tw=0: 
