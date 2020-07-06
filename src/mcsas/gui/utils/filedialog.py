# -*- coding: utf-8 -*-
# gui/utils/filedialog.py

"""
File dialogs and convenience functions.

"""

import os
from QtWidgets import QFileDialog, QDialog

from ...utils import isList

def fileDialogType():
    if sys.platform.startswith('win'):
        return QDialog
    return QFileDialog

def makeFilter(filterList):
    if filterList is None:
        return ""
    return ";;".join(filterList)

def fileDialog(parent, labeltext, path, directory = False,
               readOnly = True):
    """
    Opens a dialog to select one or more files for reading.

    Alternative to native file dialogs.
    """
    if directory and not os.path.isdir(path):
        path = os.path.dirname(path)
    dialog = QFileDialog(parent, labeltext, path)
    options = dialog.options() | QFileDialog.ReadOnly
    mode = dialog.fileMode()
    if directory:
        mode = QFileDialog.Directory
        options |= QFileDialog.ShowDirsOnly
    if not readOnly:
        options &= ~QFileDialog.ReadOnly
    dialog.setFileMode(mode)
    dialog.setOptions(options)
    dialog.setViewMode(QFileDialog.Detail)
    return dialog

def getOpenFiles(parent, labeltext, path,
                 filefilter = None, multiple = True):
    kwargs = {'options': QFileDialog.ReadOnly,
              'caption': labeltext,
              'dir': path, # why does pyside have different kwarg keys?
              'filter': makeFilter(filefilter)}
    if multiple:
        res = QFileDialog.getOpenFileNames(parent, **kwargs)
        if isList(res) and isList(res[0]):
            res = res[0]
        elif not isList(res):
            res = list(res)
        return res
    else:
        return [QFileDialog.getOpenFileName(parent, **kwargs)]

def getSaveFile(parent, labeltext, path, filefilter):
    fnfilter = makeFilter(filefilter)
    fileList = QFileDialog.getSaveFileName(
        parent,
        caption = labeltext,
        dir = path,
        filter = fnfilter)
    if isList(fileList):
        fileList = [fn for fn in fileList if fnfilter not in fn]
        fileList = fileList[0] if len(fileList) else None
    return fileList

def getSaveDirectory(parent, labeltext, path):
    if not os.path.isdir(path):
        path = os.path.dirname(path)
    dirList = QFileDialog.getExistingDirectory(
        parent,
        caption = labeltext,
        dir = path)
    return dirList

# vim: set ts=4 sts=4 sw=4 tw=0:

