# -*- coding: utf-8 -*-
# dataset.py

from abc import ABCMeta, abstractproperty
from titlemixin import TitleMixin
from originmixin import OriginMixin
from hierarchicalmixin import HierarchicalMixin
from cutesnake.utils import isString

class DisplayMixin(object):
    """Provides additional data to display in a list or tree view."""
    __metaclass__ = ABCMeta

    @property
    def displayData(self):
        return ("title", )

    @property
    def isRemovable(self):
        """Returns if this data set may be removed
        (e.g. from data lists in a GUI)"""
        return True

class ResultMixin(object):
    __metaclass__ = ABCMeta

    @property
    def result(self):
        """Supposed to return a list of *Result* types."""
        return []

class DataSet(TitleMixin, OriginMixin):

    def __init__(self, title, data, parent = None):
        TitleMixin.__init__(self, title)
        OriginMixin.__init__(self, data)

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

class HierachicalDataSet(DataSet, HierarchicalMixin):

    def __init__(self, *args, **kwargs):
        DataSet.__init__(self, *args)
        HierarchicalMixin.__init__(self, **kwargs)

    def addChild(self, dataset, position = None):
        assert isinstance(dataset, DataSet)
        dataset.setParent(self)
        if position is None:
            self._children.append(dataset)
        else:
            self._children.insert(position, dataset)
        return dataset

# vim: set ts=4 sts=4 sw=4 tw=0:
