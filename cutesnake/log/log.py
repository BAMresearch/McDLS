# -*- coding: utf-8 -*-
# log.py

"""
Interface and convenience methods for general logging.
"""

import sys
import time
import logging
import cStringIO

class Sink(object):
    buf = None

    def process(self, msg, func):
        if msg is None or len(msg) <= 0:
            return
        if self.buf is None:
            self.buf = cStringIO.StringIO()
        if msg[0] == '\t':
            self.buf.write('\n')
        self.buf.write(msg)
        if msg[-1] in ('\n', '\r'):
            func(self.buf.getvalue())
            if self.buf is not None:
                self.buf.close()
                self.buf = None

    def flush(self):
        pass

class StdOutSink(Sink):
    def write(self, msg):
        self.process(msg, logging.info)

class StdErrSink(Sink):
    def write(self, msg):
        self.process(msg, logging.error)

class Log(object):
    _formatter = None

    @classmethod
    def formatter(cls):
        """
        >>> from utils import Log
        >>> formatter = Log.formatter()

        #>>> type(formatter)
        #<class 'logging.Formatter'>
        >>> formatter._fmt
        '%(asctime)s %(levelname)-8s %(message)s'
        >>> formatter.datefmt
        '%Y-%m-%d %H:%M:%S'
        """
        if cls._formatter is None:
            cls._formatter = logging.Formatter(
                fmt='%(asctime)s %(levelname)-8s %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S')
        return cls._formatter

    @classmethod
    def timestampFormat(cls):
        return "%Y-%m-%d_%H-%M-%S"

    @classmethod
    def timestamp(cls):
        return time.strftime(cls.timestampFormat())

    @classmethod
    def hasHandler(cls):
        return len(logging.getLogger().handlers) > 0

    @classmethod
    def setConsoleHandler(cls):
        handler = logging.StreamHandler(sys.stderr)
        cls.setHandler(handler)

    @classmethod
    def setHandler(cls, handler):
        if handler is None:
            return
        cls.clearHandler()
        handler.setFormatter(cls.formatter())
        logging.getLogger().addHandler(handler)
        logging.getLogger().setLevel(logging.NOTSET)

    @classmethod
    def clearHandler(cls):
        for h in logging.getLogger().handlers:
            logging.getLogger().removeHandler(h)

    @classmethod
    def replaceStdOutErr(cls, sout = None, serr = None):
        if sout is None:
            sout = StdOutSink()
        if serr is None:
            serr = StdErrSink()
        sys.stdout = sout
        sys.stderr = serr

# vim: set ts=4 sts=4 sw=4 tw=0:
