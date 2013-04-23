# -*- coding: utf-8 -*-
# origin.py

from abc import ABCMeta
from numpy import ndarray as np_ndarray

class Origin(object):
    """
    Memorizes the original data before processing it eventually.
    """
    __metaclass__ = ABCMeta
    _origin = None

    def __init__(self, data = None):
        if data is not None:
            self.setOrigin(data)

    @property
    def origin(self):
        return self._origin

    def setOrigin(self, data):
        assert self.isValidInput(data)
        self._origin = data
        self._origin.setflags(write = False)

    @property
    def valid(self):
        return self.isValidInput(self.origin)

    @classmethod
    def isValidInput(cls, data):
        return (data is not None and
                isinstance(data, np_ndarray))

# vim: set ts=4 sts=4 sw=4 tw=0:
