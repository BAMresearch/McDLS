# -*- coding: utf-8 -*-
# titlehandler.py

from cutesnake.qt import QtCore, QtGui
from QtCore import QObject, QString
from QtGui import QWidget
from cutesnake.dataset import DataSet, TitleMixin
from cutesnake.utils import isString

class TitleHandler(TitleMixin):
    _updateFunc = None

    @classmethod
    def setup(cls, parent, title):
        """Gets a title and the widget this title belongs to."""
        th = cls(title)
        th.registerUpdateFunc(parent.setWindowTitle)
        th.update(title)
        parent.setObjectName(th.title)
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
        if isinstance(obj, DataSet):
            title = obj.fullTitle()
        if not isString(title):
            return
        TitleMixin.title.fset(self, title)
        for func in self._updateFunc:
            func(title)

    def registerUpdateFunc(self, func):
        assert callable(func)
        self._updateFunc.append(func)

# vim: set ts=4 sts=4 sw=4 tw=0:
