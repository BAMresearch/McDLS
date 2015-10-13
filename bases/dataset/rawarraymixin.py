# -*- coding: utf-8 -*-
# bases/dataset/rawarraymixin.py

from numpy import ndarray as np_ndarray

class RawArrayMixin(object):
    """
    Memorizes the original data before processing it eventually.
    """
    _rawArray = None

    def __init__(self, rawArray = None, **kwargs):
        super(RawArrayMixin, self).__init__(**kwargs)
        if rawArray is not None:
            self.setRawArray(rawArray)

    @property
    def rawArray(self):
        return self._rawArray

    def setRawArray(self, rawArray):
        assert self.isValidInput(rawArray)
        self._rawArray = rawArray
        self._rawArray.setflags(write = False)

    @property
    def valid(self):
        return self.isValidInput(self.rawArray)

    @classmethod
    def isValidInput(cls, rawArray):
        return (rawArray is not None and
                isinstance(rawArray, np_ndarray))

# vim: set ts=4 sts=4 sw=4 tw=0:
