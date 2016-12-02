# -*- coding: utf-8 -*-
# bases/dataset/hierarchicalmixin.py

from builtins import object
from abc import ABCMeta
from future.utils import with_metaclass

class HierarchicalMixin(with_metaclass(ABCMeta, object)):
    """
    Manages a list of objects of the same type which makes up a tree
    structure.
    """
    _children = None
    _parent = None

    def __init__(self, parent = None, **kwargs):
        super(HierarchicalMixin, self).__init__(**kwargs)
        self._children = []
        self.setParent(parent)

    @property
    def parent(self):
        return self._parent

    def setParent(self, parent):
        self._parent = parent

    @property
    def topLevelParent(self):
        p = self
        if p.parent is None:
            return None
        # traverse hierarchy
        while p.parent is not None:
            p = p.parent
        return p

    def __len__(self):
        return len(self._children)

    def __iter__(self):
        for c in self._children:
            yield c

    def __getitem__(self, position):
        return self._children[position]

    def __delitem__(self, position):
        del self._children[position]

# vim: set ts=4 sts=4 sw=4 tw=0:
