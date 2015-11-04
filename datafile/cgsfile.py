# -*- coding: utf-8 -*-
# datafile/pdhfile.py

from __future__ import absolute_import # PEP328
import logging
import datetime
import re
from numpy import array as np_array
from collections import OrderedDict
from utils import isList
from utils.units import Angle, Temperature, DynamicViscosity, Length, MSec
from datafile import AsciiFile
from utils import classproperty
from dataobj import DLSData

import sys

def text2num(valueText, dtype = float):
    value = None
    try:
        value = dtype(valueText)
    except:
        return None
    return value

def _makeProperty(varName):
    def getter(selforcls):
        return getattr(selforcls, varName)
    return property(getter)

class CGSFile(AsciiFile):
    """Parsing *ALV-7004 CGS-8F Data, Single Run Data*
    Read-only at the moment."""
    # a set of attributes found in the data file:
    # pairs of the name in the file for lookup and the associated member name
    # if both are the same (besides of case) only one name is required
    _knownProps = OrderedDict((
                    ("Date", "date"),
                    ("Time", "time"),
                    ("timestamp", None),
                    ("Samplename", "sampleName"),
                    ("SampMemo", "sampleMemo"),
                    ("Temperature", "temperature"),
                    ("Viscosity", "viscosity"),
                    ("Refractive Index", "refractiveIndex"),
                    ("Wavelength", "wavelength"),
                    ("Angle", "angles"),
                    ("Duration", "duration"),
                    ("Runs", "runs"),
                    ("Mode", "mode"),
                    ("MeanCR", "meanCR"),
                    ("Atten", "attenuation"),
                    ("RelSens", "relSens"),
                    ("DC", "dc"),
                    ("Correlation", "correlation"),
                    ("Count Rate", "countRate"),
                    ("Monitor Diode", "monitorDiode"),
                    # dict containing attribute specific units from file if found
                    ("units", None),
    ))
    # regexp that matches the keys including index number and unit text
    # match group consists of: <key>, (<idx>), <idx>, [<unit>], unit
    # whereby <key> must not consist of digits,
    # ending digits within keys are interpreted as indices
    _keyPattern = re.compile("(?P<key>[^\(\[0-9]+)"             # key
                             "([\(\[]?(?P<idx>[0-9])[\)\]]?)?"  # index
                             "(\s*\[(?P<unit>[^0-9\[\]]+)\])?"  # unit
                             )
    _tauUnit = MSec # default for known file signatures

    @classproperty
    @classmethod
    def knownProperties(cls):
        return cls._knownProps

    @classmethod
    def setPropertyGetters(cls):
        for inName, memName in cls.knownProperties.iteritems():
            if memName is None or not len(memName):
                memName = inName
            intName = "_" + memName
            setattr(cls, intName, None) # init value
            setattr(cls, memName, _makeProperty(intName))

    @classproperty
    @classmethod
    def fileFilter(cls):
        return (("DLS ALV Data", "asc"),)

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
        i = 0
        self.verifySignatures(asciiLines[i].strip())
        while i < len(asciiLines) - 1:
            i += 1
            line = asciiLines[i].strip()
#            print >>sys.__stderr__, "line", i, line
            if not len(line):
                continue
            if line.strip('"') in self.knownProperties:
                # parse float array
                i, rawArray = self.readArray(asciiLines, startLine = i+1)
                key, idx, unit = self.processKey(line.strip('"'))
                setattr(self, "_"+key, rawArray)
                continue
            key, value = self.splitKeyValue(line)
#            print >>sys.__stderr__, u"k: '{}', v: '{}'".format(key, value)
            key, idx, unit = self.processKey(key)
#            print >>sys.__stderr__, "processKey out", key, idx, unit
            if key is None:
                continue # attribute not supported
            funcName = "parse" + key[0].upper() + key[1:]
            parseFunc = getattr(self, funcName, None)
#            print >>sys.__stderr__, "parseFunc", parseFunc
            if parseFunc is None:
                logging.warning("Parse function '{}' not found!".format(funcName))
                continue
            parseFunc(value, text2num(idx, dtype = int))
            self.setUnit(key, unit)
#            print >>sys.__stderr__, "  --> p: '{}'".format(getattr(self, key, None))

        for dataName, memName in self.knownProperties.iteritems():
            if memName is None or not len(memName): memName = dataName
            value = getattr(self, memName)
#            print >>sys.__stderr__, "{} ({}): \t{}".format(memName, type(value), value)
#            try:
#                print >>sys.__stderr__, value.shape
#            except: pass

    def splitKeyValue(self, line):
        if not len(line):
            return None, None
        key, value = None, None
        try:
            key, value = line.split(':', 1)
        except:
            try:
                lst = line.split()
                key, value = " ".join(lst[:-1]), lst[-1]
            except:
                pass
        key, value = key.strip(), value.strip().strip('"')
        return key, value

    def processKey(self, inKey):
#        print >>sys.__stderr__, "processKey  in", inKey
        outKey = self.knownProperties.get(inKey, None)
        if outKey is not None:
            return outKey, None, None
        result = self._keyPattern.match(inKey)
        if result is None or len(result.groups()) < 5:
            return None, None, None
#        print >>sys.__stderr__, "match:", result.groups()
        inKey = result.group("key").strip()
        outKey = self.knownProperties.get(inKey, None)
        if outKey is None or not len(outKey):
            return None, None, None
        return outKey, result.group("idx"), result.group("unit")

    def setUnit(self, key, unit):
        if key is None or not len(key) or unit is None or not len(unit):
            return
        if self._units is None:
            self._units = dict()
        self._units[key] = unit

    @property
    def units(self):
        return self._units

    def parseDate(self, text, *args):
        values = [int(v) for v in text.split('/')]
        values.reverse()
        self._date = datetime.date(*values)
        self.setTimestamp()

    def parseTime(self, text, *args):
        values = [int(v) for v in text.split(':')]
        self._time = datetime.time(*values)
        self.setTimestamp()

    def setTimestamp(self):
        if not (isinstance(self.date, datetime.date) and
                isinstance(self.time, datetime.time)):
            return
        self._timestamp = datetime.datetime.combine(self.date, self.time)

    def parseSampleName(self, text, *args):
        self._sampleName = text

    def setNumberValue(self, propName, valueText, dtype = float):
        setattr(self, "_"+propName, text2num(valueText, dtype))

    def setListValue(self, propName, value, idx):
        propName = "_"+propName
        if idx is None:
            setattr(self, propName, value)
        else:
            if not isList(getattr(self, propName)):
                setattr(self, propName, list())
            lst = getattr(self, propName)
            idx = min(idx, len(lst))
            lst.insert(idx, value)

    def setArray(self, propName, rawArray):
        if propName not in self.knownProperties:
            return False
        setattr(self, propName, rawArray)
        return True

    def parseSampleMemo(self, *args):
        self.setListValue("sampleMemo", *args)

    def parseTemperature(self, text, *args):
        self.setNumberValue("temperature", text)

    def parseViscosity(self, text, *args):
        self.setNumberValue("viscosity", text)

    def parseRefractiveIndex(self, text, *args):
        self.setNumberValue("refractiveIndex", text)

    def parseWavelength(self, text, *args):
        self.setNumberValue("wavelength", text)

    def parseAngles(self, text, *args):
        self.setListValue("angles", text2num(text), *args)

    def parseDuration(self, text, *args):
        self.setNumberValue("duration", text)

    def parseRuns(self, text, *args):
        self.setNumberValue("runs", text, dtype = int)

    def parseMode(self, text, *args):
        self._mode = text

    def parseMeanCR(self, text, *args):
        self.setListValue("meanCR", text2num(text), *args)

    def parseAttenuation(self, text, *args):
        self.setListValue("attenuation", text2num(text), *args)

    def parseRelSens(self, text, *args):
        self.setListValue("relSens", text2num(text), *args)

    def parseDc(self, text, *args):
        self.setListValue("dc", text2num(text), *args)

    def parseMonitorDiode(self, text, *args):
        self.setNumberValue("monitorDiode", text)

    def getDataObj(self):
        dlsData = DLSData(title = self.name)
        dlsData.setFilename(self.filename)
        dlsData.setSampleName(self.sampleName)
        dlsData.setSampleDescription(
                " ".join([s for s in self.sampleMemo if len(s)]))
        dlsData.setTemperature(Temperature(self.units['temperature'])
                    .toSi(self.temperature))
        dlsData.setViscosity(DynamicViscosity(self.units['viscosity'])
                    .toSi(self.viscosity))
        dlsData.setRefractiveIndex(self.refractiveIndex) # dimensionless
        dlsData.setWavelength(Length(self.units['wavelength'])
                .toSi(self.wavelength))
        # scattering angles
        unit = self.units['angles']
        angles = np_array(self.angles)
        dlsData.setAngles(Angle(unit).toSi(angles))
        # correlation array
        dlsData.setTau(self._tauUnit.toSi(self.correlation[:, 0]))
        dlsData.setCorrelation(self.correlation[:, 1:])
        return dlsData

CGSFile.setPropertyGetters()

# vim: set ts=4 sw=4 sts=4 tw=0:
