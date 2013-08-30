# -*- coding: utf-8 -*-
# __init__.py

import logging
from log import (timestampFormat, timestamp, replaceStdOutErr,
                 replaceHandler, addHandler, removeHandler)
from widgethandler import WidgetHandler

log = logging.getLogger()

# vim: set ts=4 sw=4 sts=4 tw=0:
