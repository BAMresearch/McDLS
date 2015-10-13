# -*- coding: utf-8 -*-
# bases/dataset/dataset.py

from abc import ABCMeta, abstractproperty
from titlemixin import TitleMixin
from rawarraymixin import RawArrayMixin
from hierarchicalmixin import HierarchicalMixin
from utils import isString, classproperty

class DisplayMixin(object):
    """Provides additional data to display in a list or tree view."""
    __metaclass__ = ABCMeta

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

class ResultMixin(object):
    __metaclass__ = ABCMeta

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

class HierarchicalDataSet(DataSet, HierarchicalMixin):

    def addChild(self, dataset, position = None):
        assert isinstance(dataset, DataSet)
        dataset.setParent(self)
        if position is None:
            self._children.append(dataset)
        else:
            self._children.insert(position, dataset)
        return dataset

    def fullTitle(self, titleFormat = None):
        """TitleFormat defines how parent titles will be concatenated with
        children titles"""
        if not isString(titleFormat):
            titleFormat = "{0} > {1}"
        def parents():
            dataset = self
            while dataset is not None:
                yield dataset
                dataset = dataset.parent
        return reduce(lambda a, b: titleFormat.format(b, a),
                      [ds.title for ds in parents()])

    def __init__(self, *args, **kwargs):
        super(HierarchicalDataSet, self).__init__(**kwargs)

# vim: set ts=4 sts=4 sw=4 tw=0:
