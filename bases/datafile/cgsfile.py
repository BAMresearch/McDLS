# -*- coding: utf-8 -*-
# bases/datafile/pdhfile.py

import logging
import datetime
from utils import isInteger
#from utils.classproperty import classproperty
from bases.datafile import AsciiFile
from utils.classproperty import classproperty

import sys

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

class CGSFile(AsciiFile):
    """Parsing *ALV-7004 CGS-8F Data, Single Run Data*
    Read-only at the moment."""
    # FIXME? unordered?
    _knownProps = set(("date", "time", "timestamp", "sampleName", "sampleMemo",))

    @classproperty
    @classmethod
    def knownProperties(cls):
        return cls._knownProps

    @classmethod
    def setPropertyGetters(cls):
        for name in cls.knownProperties:
            intName = "_" + name
            setattr(cls, intName, None) # init value
            setattr(cls, name, _makeProperty(intName))

    @classproperty
    @classmethod
    def extensions(cls):
        return ("asc", )

    @classmethod
    def signatures(cls):
        return ("ALV-7004 CGS-8F Data, Single Run Data",)

    @classmethod
    def verifySignatures(cls, fileSignature):
        if fileSignature in cls.signatures():
            return True
        msg = ("Found unknown file signature:",
               "    '{0}'! ".format(fileSignature),
               "Please report this new file type to "
               "the program authors.")
#        [logging.warning(m) for m in msg]
        raise IOError("\n".join(msg))
        return False

    def parseLines(self, asciiLines, **kwargs):
        self.verifySignatures(asciiLines[0].strip())
        print >>sys.__stderr__,dir(self)
        for line in asciiLines[1:10]:
            key, value = line.strip().split(':', 1)
            key, value = key.strip(), value.strip().strip('"')
            key = key[0].lower() + key[1:]
            print >>sys.__stderr__, "k: '{}', v: '{}'".format(key, value)
            if key == "samplename":
                key = "sampleName" # fix this key from file
            if key not in self.knownProperties:
                continue
            parseFunc = getattr(self, "parse" + key[0].upper() + key[1:], None)
            if parseFunc is None:
                continue
            parseFunc(value)
            print >>sys.__stderr__, "  --> p: '{}'".format(getattr(self, key, None))

        for prop in self.knownProperties:
            print >>sys.__stderr__, "{}: \t{}".format(prop, getattr(self, prop))
        sys.exit(0)

    def parseDate(self, text):
        values = [int(v) for v in text.split('/')]
        values.reverse()
        self._date = datetime.date(*values)
        self.setTimestamp()

    def parseTime(self, text):
        values = [int(v) for v in text.split(':')]
        self._time = datetime.time(*values)
        self.setTimestamp()

    def setTimestamp(self):
        if not (isinstance(self.date, datetime.date) and
                isinstance(self.time, datetime.time)):
            return
        self._timestamp = datetime.datetime.combine(self.date, self.time)

    def parseSampleName(self, text):
        self._sampleName = text

CGSFile.setPropertyGetters()

# vim: set ts=4 sw=4 sts=4 tw=0:
