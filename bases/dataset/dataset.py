# -*- coding: utf-8 -*-
# bases/dataset/dataset.py

from builtins import object
from abc import ABCMeta, abstractproperty
from .titlemixin import TitleMixin
from .rawarraymixin import RawArrayMixin
from utils import isString, classproperty
from future.utils import with_metaclass
from functools import reduce

class DisplayMixin(with_metaclass(ABCMeta, object)):
    """Provides additional data to display in a list or tree view."""

    @abstractproperty
    @classproperty
    @classmethod
    def displayDataDescr(cls):
        """Head row title for the properties to be shown when listed in a
        DataList UI element. See also displayData()."""
        return ("data title", )

    @abstractproperty
    @classproperty
    @classmethod
    def displayData(cls):
        """Properties to be shown when listed in a DataList UI element.
        See also displayDataDescr()."""
        return ("title", )

    @property
    def isRemovable(self):
        """Returns if this data set may be removed
        (e.g. from data lists in a GUI)"""
        return True

    def __init__(self, **kwargs):
        super(DisplayMixin, self).__init__(**kwargs)

class ResultMixin(with_metaclass(ABCMeta, object)):
    @property
    def result(self):
        """Supposed to return a list of *Result* types."""
        return []

    def __init__(self, **kwargs):
        super(ResultMixin, self).__init__(**kwargs)

class DataSet(TitleMixin, RawArrayMixin):
    """Container base class for all kinds of data to be passed around in the
    UI. Knows its originally loaded (raw) data array and its title to be shown
    in a UI."""

    def __init__(self, **kwargs):
        super(DataSet, self).__init__(**kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
