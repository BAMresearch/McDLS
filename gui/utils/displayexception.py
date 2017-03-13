# -*- coding: utf-8 -*-
# gui/utils/displayexception.py

from __future__ import absolute_import # PEP328
from builtins import str
import logging
import traceback
import inspect
from gui.utils.translate import tr
from gui.qt import QtGui
from QtGui import QMessageBox, QApplication

class DisplayException(QMessageBox):
    """
    Displays a message box with the text of a provided exception.
    *level*: one of 'question', 'info', 'warning', 'critical'

    >>> from utilsgui import DisplayException, DialogInteraction
    >>> from utils import AppError
    >>> de = DialogInteraction.queryInstance(DisplayException,
    ...                                      AppError("test"),
    ...                                      button = 'default')
    >>> type(de).__name__
    'DisplayException'
    >>> de.result() == 0
    True
    """
    _icons =   {'question': QMessageBox.Question,
                'info':     QMessageBox.Information,
                'warning':  QMessageBox.Warning,
                'critical': QMessageBox.Critical }
    _titles =  {'question': tr("Question"),
                'info':     tr("Information"),
                'warning':  tr("Warning"),
                'critical': tr("An Error Occurred") }
    _logging = {'question': logging.info,
                'info':     logging.info,
                'warning':  logging.warning,
                'critical': logging.error }

    def __init__(self, exception, level = 'critical', fmt = None):
        QMessageBox.__init__(self, parent = QApplication.activeWindow())
        icon = self._icons.get(level, QMessageBox.NoIcon)
        self.setIcon(icon)
        title = self._titles.get(level, "")
        self.setWindowTitle(title)
        okBtn = self.addButton("ok", QMessageBox.AcceptRole)
        self.setDefaultButton(okBtn)
        if fmt is None:
            fmt = u"<nobr>{e}</nobr>"
        className, methodName = self.classAndMethodName()
        excName, excMsg = type(exception).__name__, str(exception)
        text = u"{3}\n[ {2} in {0}.{1} ]".format(className, methodName, excName, excMsg)
        self.setText(fmt.format(e = text.replace("\n", "<br />")))
        logfunc = self._logging.get(level, logging.error)
        logfunc(traceback.format_exc()) # log exception traceback
        self.exec_()

    @classmethod
    def classAndMethodName(cls):
        """Retrieves some meta info to decribe the error further."""
        trace = inspect.trace()[-1]
        frame = trace[0]
        inClass = frame.f_locals.get('self', None)
        if inClass is None:
            inClass = frame.f_locals.get('cls', type(None))
        else:
            inClass = type(inClass)
        className = ""
        if hasattr(inClass, "__name__"):
            className = inClass.__name__
        elif isinstance(inClass, super):
            className = inClass.__self_class__.__name__
        if className == "NoneType":
            className = ""
        return className, trace[3]+"()"

# vim: set ts=4 sts=4 sw=4 tw=0:
