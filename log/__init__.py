# -*- coding: utf-8 -*-
# log/__init__.py

import logging
from .log import (timestampFormat, timestamp, timestampFormatted, formatter,
                  replaceStdOutErr, replaceHandler, addHandler, removeHandler)
from .widgethandler import WidgetHandler

log = logging.getLogger()

def getWidgetHandlers():
    """Returns all active WidgetHandlers for logging."""
    return list((h for h in log.handlers if isinstance(h, WidgetHandler)))

def disableFormatter():
    for h in log.handlers:
        h.setFormatter(logging.Formatter())

def enableFormatter():
    for h in log.handlers:
        h.setFormatter(formatter())

# vim: set ts=4 sw=4 sts=4 tw=0:
