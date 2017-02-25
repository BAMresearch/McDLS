# -*- coding: utf-8 -*-
# utils/pickleinstancemethods.py

from future import standard_library
standard_library.install_aliases()
from types import MethodType
import copyreg

# https://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods#post2155350

def _pickle_method(method):
    func_name = method.__func__.__name__
    obj = method.__self__
    cls = method.__self__.__class__
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

copyreg.pickle(MethodType, _pickle_method, _unpickle_method)

# vim: set ts=4 sts=4 sw=4 tw=0:
