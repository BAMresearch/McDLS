# -*- coding: utf-8 -*-
# dataset.py

from title import Title
from origin import Origin
from hierarchical import Hierarchical
from utils import isString

class DataSet(Title, Origin, Hierarchical):

    def __init__(self, title, data, parent = None):
        Title.__init__(self, title)
        Origin.__init__(self, data)
        Hierarchical.__init__(self, parent)

    @property
    def isRemovable(self):
        """Returns if this data set may be removed
        (e.g. from data lists in the GUI)"""
        return True

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

# vim: set ts=4 sts=4 sw=4 tw=0:
