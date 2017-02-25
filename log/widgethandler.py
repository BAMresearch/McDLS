# -*- coding: utf-8 -*-
# widgethandler.py

import logging

class WidgetHandler(logging.Handler):
    """
    A logging.Handler which appends messages to some widget.append().

    >>> import logging
    >>> from PyQt4.QtGui import QPlainTextEdit
    >>> from utilsgui import DialogInteraction
    >>> from logdock import LogHandler
    >>> pte = DialogInteraction.instance(QPlainTextEdit)
    >>> handler = LogHandler(pte)
    >>> msg = "testest"
    >>> handler.emit(logging.LogRecord("", 0, "", 0, msg, (), None))
    >>> str(pte.toPlainText()) == msg
    True
    """
    _widget = None # (Qt) widget to write log data to

    def __init__(self, widget):
        logging.Handler.__init__(self)
        self._widget = widget
        assert hasattr(widget, 'append'), \
                "append() method required for given widget!"

    @property
    def widget(self):
        return self._widget

    def emit(self, record):
        if self._widget is None:
            return
        msg = self.format(record)
        self._widget.append(msg)

# vim: set sts=4 ts=4 sw=4 tw=0: 
