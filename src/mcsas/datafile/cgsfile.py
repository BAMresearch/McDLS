# -*- coding: utf-8 -*-
# datafile/pdhfile.py

import logging
import datetime
import re
import os.path
from numpy import array as np_array
from collections import OrderedDict

from ..utils.units import Angle, Temperature, DynamicViscosity, Length, MSec, Sec
from . import AsciiFile
from ..utils import classproperty, isList
from ..dataobj import DLSData

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
                    ("measIndex", None),
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
        for inName, memName in cls.knownProperties.items():
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
        return ("ALV-7004 CGS-8F Data, Single Run Data",
                "Simulation Data")

    @classmethod
    def verifySignatures(cls, fileSignature):
        if any((sig in fileSignature) for sig in cls.signatures()):
            return True
        msg = ("Found unknown file signature:",
               "    '{0}'! ".format(fileSignature),
               "Please report this new file type to "
               "the program authors.")
#        [logging.warning(m) for m in msg]
        raise IOError("\n".join(msg))
        return False

    def parseLines(self, asciiLines, **kwargs):
        """Will be called by the base class' readFile() method."""
        i = 0
        # use the file signature as fallback replacement for the sample name
        self._sampleName = asciiLines[i].strip()
        self.verifySignatures(self.sampleName)
        while i < len(asciiLines) - 1:
            i += 1
            line = asciiLines[i].strip()
            if not len(line):
                continue
            if line.strip('"') in self.knownProperties:
                # parse float array
                i, arr = self.readArray(asciiLines, startLine = i+1)
                key, idx, unit = self.processKey(line.strip('"'))
                setattr(self, "_"+key, arr)
                continue
            key, value = self.splitKeyValue(line)
            key, idx, unit = self.processKey(key)
            if key is None:
                continue # attribute not supported
            funcName = "parse" + key[0].upper() + key[1:]
            parseFunc = getattr(self, funcName, None)
            if parseFunc is None:
                logging.warning("Parse function '{}' not found!".format(funcName))
                continue
            parseFunc(value, text2num(idx, dtype = int))
            self.setUnit(key, unit)

        for dataName, memName in self.knownProperties.items():
            if memName is None or not len(memName): memName = dataName
            value = getattr(self, memName)

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
        outKey = self.knownProperties.get(inKey, None)
        if outKey is not None:
            return outKey, None, None
        result = self._keyPattern.match(inKey)
        if result is None or len(result.groups()) < 5:
            return None, None, None
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
        try:
            values = [int(v) for v in text.split('/')]
        except ValueError: # int(v) fails
            return
        values.reverse()
        self._date = datetime.date(*values)
        self.setTimestamp()

    def parseTime(self, text, *args):
        try:
            values = [int(v) for v in text.split(':')]
        except ValueError: # int(v) fails
            return
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
            if not isList(value):
                value = [value,]
            setattr(self, propName, value)
        else:
            if not isList(getattr(self, propName)):
                setattr(self, propName, list())
            lst = getattr(self, propName)
            idx = min(idx, len(lst))
            lst.insert(idx, value)

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

    def setFilename(self, fn):
        super(CGSFile, self).setFilename(fn)
        self._setMeasIndexFromFilename(self.filename)

    def _setMeasIndexFromFilename(self, fn):
        fn = os.path.basename(fn)
        fn, dummy = os.path.splitext(fn)
        result = re.match(".*(?P<mgroup>[0-9]{4})_(?P<meas>[0-9]{4})$", fn)
        groupIndex, measIndex = None, None
        try:
            groupIndex = int(result.group("mgroup"))
        except:
            pass
        try:
            measIndex = int(result.group("meas"))
        except:
            pass
        self._measIndex = groupIndex, measIndex

    def getDataObj(self):
        dlsData = DLSData(title = self.sampleName)
        dlsData.setFilename(self.filename)
        if any([idx is not None for idx in self.measIndex]):
            dlsData.setMeasIndices((self.measIndex,))
        else:
            dlsData.setMeasIndices(())
        dlsData.setSampleName(self.sampleName)
        if isList(self.sampleMemo):
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
        aunit = Angle(self.units['angles'])
        if self.angles is not None:
            dlsData.setAngles(aunit, aunit.toSi(np_array(self.angles)))
        # correlation array
        dlsData.setTau(self._tauUnit, self.correlation[:, 0])
        dlsData.setCorrelation(self.correlation[:, 1:])
        if isList(self.countRate):
            dlsData.setCapTime(Sec, self.countRate[:, 0])
            dlsData.setCountRate(self.countRate[:, 1:])
        dlsData.initConfig() # has to be called after setup
        return dlsData

CGSFile.setPropertyGetters()

# vim: set ts=4 sw=4 sts=4 tw=0:
