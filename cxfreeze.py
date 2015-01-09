# -*- coding: utf-8 -*-
# cxfreeze.py
#
# Usage on target platform:
#  python cxfreeze.py build_exe
# tested on Windows with MinGW/MSYS

import sys
import gui.version

OLDVERSIONNUM = gui.version.version.number()
if len(sys.argv) > 2:
    alternateVersion = "-".join(sys.argv[2:])
    del sys.argv[2:]
    print "Using an alternate version: '{0}'".format(alternateVersion)
    gui.version.version.updateFile(gui.version, alternateVersion)
    reload(gui.version)

import re
import os.path
import subprocess
import hashlib
import platform
from cx_Freeze import setup, Executable
from gui.version import version
from utils import isWindows, isLinux, isMac, testfor

def sanitizeVersionNumber(number):
    """Removes non-digits to be compatible with pywin32"""
    if not isWindows():
        return number
    match = re.match(r'[\d.]+', number)
    if match is None:
        return ""
    number = match.group(0)
    try:
        bits = [int(i) for i in number.split(".")]
        assert(len(bits))
    except (AssertionError, IndexError, TypeError, ValueError):
        number = "0" # do not set version number in win32 binary
    return number

class Archiver(object):
    _name = None # to be defined by sub classes
    _path = None

    def __init__(self):
        self._path = self._getPath()
        testfor(os.path.isfile(self._path), OSError,
                "{name}: '{path}' not found!".format(
                    name = self._name,
                    path = self._path))

    @property
    def execName(self):
        return os.path.basename(self._path)

    def getLogFilename(self):
        return os.path.splitext(self.execName)[0] + ".log" # 7z.log

    def archive(self, targetPath):
        raise NotImplementedError

class Archiver7z(Archiver):
    _name = "7-Zip"

    @staticmethod
    def _getPath():
        """Retrieves 7-Zip install path and binary from Windows Registry.
        This way, we do not rely on PATH containing 7-Zip.
        """
        path = None
        if isWindows():
            import win32api, win32con
            for key in (win32con.HKEY_LOCAL_MACHINE, win32con.HKEY_CURRENT_USER):
                try:
                    key = win32api.RegConnectRegistry(None, key)
                    key = win32api.RegOpenKeyEx(key, r"SOFTWARE\7-Zip")
                    path, dtype = win32api.RegQueryValueEx(key, "Path")
                    path = os.path.abspath(os.path.join(path, "7z.exe"))
                except:
                    pass
        else:
            p = subprocess.Popen(["which", "7z"],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE)
            out, err = p.communicate()
            path = out.strip()
        return path

    def archive(self, targetPath):
        """Expects an absolute target directory path"""
        fnPackage = targetPath + ".7z"
        fnLog = self.getLogFilename()
        with open(fnLog, 'w') as fd:
            retcode = subprocess.call([self._path, "a", "-t7z", "-mx=9",
                                       fnPackage, targetPath],
                                       stdout = fd,
                                       stderr = fd)
        if not os.path.exists(fnPackage):
            fnPackage = None
        return fnPackage

class ArchiverZip(Archiver):
    _name = "Zip"

    @staticmethod
    def _getPath():
        path = None
        if isWindows():
            pass # none yet
        else:
            p = subprocess.Popen(["which", "zip"],
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE)
            out, err = p.communicate()
            path = out.strip()
        return path

    def archive(self, targetPath):
        """Expects an absolute target directory path"""
        fnPackage = targetPath + ".zip"
        fnLog = self.getLogFilename()
        with open(fnLog, 'w') as fd:
            retcode = subprocess.call([self._path, "-r9",
                                       fnPackage, targetPath],
                                       stdout = fd,
                                       stderr = fd)
        if not os.path.exists(fnPackage):
            fnPackage = None
        return fnPackage

archiver = None
if isMac():
    archiver = ArchiverZip()
else:
    archiver = Archiver7z()

# target (temp) dir for mcsas package
TARGETDIR = "{name}-{ver}_{plat}".format(
                name = version.name(),
                ver = version.number(),
                plat = platform.system().lower())

BASE = None
EXEC_SUFFIX = ""
if sys.platform in "win32":
    BASE = "Win32GUI"
    EXEC_SUFFIX = ".exe"

TARGETNAME = version.name() + EXEC_SUFFIX

# (source, target) pairs
# without a target the file is placed in the top level directory of the package
INCLUDEFILES = [
        "mcsas/mcsasparameters.json",
        ("resources/background_files.svg", "resources/background_files.svg"),
        ("resources/background_ranges.svg", "resources/background_ranges.svg"),
        "dejavuserif.ttf",
]
if isLinux():
    INCLUDEFILES += [
        "/usr/lib/liblapack.so.3",
        "/usr/lib/libblas.so.3",
        "/usr/lib/x86_64-linux-gnu/libgfortran.so.3",
        "/usr/lib/x86_64-linux-gnu/libquadmath.so.0",
        "/usr/lib/x86_64-linux-gnu/libpyside-python2.7.so.1.2",
        "/usr/lib/x86_64-linux-gnu/libshiboken-python2.7.so.1.2",
    ]
if isWindows():
    INCLUDEFILES += [
        'Microsoft.VC90.CRT',
    ]

BUILDOPTIONS = dict(
    compressed = False,
    include_files = INCLUDEFILES,
    packages = [],
    includes = ["PySide", "PySide.QtCore", "PySide.QtGui",
                "PySide.QtSvg", "PySide.QtXml",
                "multiprocessing",
                "scipy.sparse.csgraph._validation",
                "scipy.sparse.linalg.dsolve.umfpack", "scipy.integrate.vode",
                "scipy.integrate.lsoda", "matplotlib",
                "matplotlib.backends.backend_qt4agg",
                # savefig dependencies?
                "matplotlib.backends.backend_tkagg", "Tkinter", "FileDialog",
                ],
#    excludes = ["Tkinter"],
#    icon = "res/img/brianpauw.ico",
    path = [os.getcwd()] + sys.path,
    build_exe = TARGETDIR,
    silent = False,
    copy_dependent_files = True,
)

# OSX bundle building
MACOPTIONS = dict(bundle_name = TARGETDIR)
DMGOPTIONS = dict(volume_label = TARGETDIR)
if isMac():
    BUILDOPTIONS.pop("build_exe") # bdist_mac expects plain 'build' directory
    BUILDOPTIONS["includes"] = [ # order in which they were requested
        "PySide", "PySide.QtCore", "PySide.QtGui", "PySide.QtSvg", "PySide.QtXml",
        "scipy.sparse.csgraph._validation",
        "scipy.sparse.linalg.dsolve.umfpack",
        "matplotlib.backends.backend_qt4agg",
        "scipy.integrate.vode",
        "scipy.integrate.lsoda",
    ]
    BUILDOPTIONS["bin_includes"] = [
        "libpyside-python2.7.1.1.dylib",
        "libshiboken-python2.7.1.1.dylib",
    ]
    if False:
        BUILDOPTIONS["replace_paths"] = [
        ]
    #BUILDOPTIONS["copy_dependent_files"] = False

setup(
    name = version.name(),
    version = sanitizeVersionNumber(version.number()),
    description = version.name(),
    long_description = ("GUI for Monte-Carlo size distribution analysis"),
#    url = "http://sourceforge.net/",
    license = "Creative Commons CC-BY-SA",
    author = u"Brian R. Pauw",
    author_email = "brian@stack.nl",
    contact = u"Brian R. Pauw",
    contact_email = "brian@stack.nl",
    maintainer = u"Ingo Breßler",
    maintainer_email = "ingo.bressler@bam.de",
    # additional metadata for modified version of cx_Freeze
    # displayed in copyright info
    download_url = "2014, https://bitbucket.org/pkwasniew/mcsas",
    # company
    classifiers = u"NIMS\r\n"
                  u"National Institute for Materials Science, \r\n\r\n"
                  u"1-2-1 Sengen, 305-0047, "
                  u"305-0047, Tsukuba, Japan",
    options = dict(build_exe = BUILDOPTIONS,
                   bdist_mac = MACOPTIONS,
                   bdist_dmg = DMGOPTIONS),
    executables = [Executable("main.py", base = BASE,
                              targetName = TARGETNAME)])

# zip:
# zip -r9 MCSAS-1.0.zip MCSAS-1.0/
# package the freezed program into an 7z archive
PACKAGEFN = archiver.archive(TARGETDIR)

# calc a checksum of the package
def hashFile(filename, hasher, blocksize = 65536):
    with open(filename, 'rb') as fd:
        buf = fd.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = fd.read(blocksize)
    return hasher.hexdigest(), os.path.basename(filename)
hashValue = hashFile(PACKAGEFN, hashlib.sha256())

# write the checksum to file
with open(version.name() + ".sha", 'w') as fd:
    fd.write(" *".join(hashValue))

# restore initially modified version
if version.number() != OLDVERSIONNUM:
    version.updateFile(gui.version, OLDVERSIONNUM)

# vim: set ts=4 sw=4 sts=4 tw=0:
