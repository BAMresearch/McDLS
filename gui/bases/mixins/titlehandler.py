# -*- coding: utf-8 -*-
# gui/bases/mixins/titlehandler.py

from QtCore import QObject
from QtWidgets import QWidget
from bases.dataset import TitleMixin
from utils import isString, isCallable

class TitleHandler(TitleMixin):
    _updateFunc = None

    @classmethod
    def setup(cls, parent, title):
        """Gets a title and the widget this title belongs to."""
        th = cls(title)
        th.registerUpdateFunc(parent.setWindowTitle)
        th.update(title)
        objectName = th.title.title().replace(" ","")
        parent.setObjectName(objectName)
        return th

    def __init__(self, title):
        TitleMixin.__init__(self, title)
        self._updateFunc = []

    def __call__(self, obj = None):
        if obj is not None:
            self.update(obj)
        return self.title

    def update(self, obj):
        title = obj
        if not isString(title):
            return
        TitleMixin.title.fset(self, title)
        for func in self._updateFunc:
            func(title)

    def registerUpdateFunc(self, func):
        assert isCallable(func)
        self._updateFunc.append(func)

# vim: set ts=4 sts=4 sw=4 tw=0:
