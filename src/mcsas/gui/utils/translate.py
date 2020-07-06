# -*- coding: utf-8 -*-
# gui/utils/translate.py

from QtWidgets import QApplication

def tr(s):
    return QApplication.translate(None, s)

fromUtf8 = str

# vim: set ts=4 sts=4 sw=4 tw=0:
