# -*- coding: utf-8 -*-
# sink.py

"""
Interface and convenience methods for general logging.
"""

import logging
from cStringIO import StringIO

class Sink(object):
    buf = None

    def process(self, msg, func):
        if msg is None or len(msg) <= 0:
            return
        if self.buf is None:
            self.buf = StringIO()
        if msg[0] == '\t':
            self.buf.write('\n')
        self.buf.write(msg)
        if msg[-1] in ('\n', '\r'):
            # call the logging function, removing duplicate newlines
            func(self.buf.getvalue().rstrip("\n"))
            self.buf = StringIO()

    def flush(self):
        pass

class StdOutSink(Sink):
    def write(self, msg):
        self.process(msg, logging.info)

class StdErrSink(Sink):
    def write(self, msg):
        self.process(msg, logging.error)

# vim: set ts=4 sts=4 sw=4 tw=0:
