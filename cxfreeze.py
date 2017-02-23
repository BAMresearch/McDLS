# -*- coding: utf-8 -*-
# cxfreeze.py

"""
Overview
========

Creates a standalone program package for a particular platform to be run by
restricted users without installing any additional packages.

This script is executable and has to be run on the platform for which a
package shall be created. Please follow the instructions below for each
particular platform.

Common Package Dependencies:

    - `Python 2.7 <https://www.python.org/download/releases/2.7/>`_
    - `Qt 4.8 <http://qt-project.org/doc/qt-4.8/qt4-8-intro.html>`_ + `PySide <http://qt-project.org/wiki/Category:LanguageBindings::PySide::Downloads>`_
    - `NumPy and SciPy <http://www.scipy.org/scipylib/download.html>`_
        In order to work around freeze failures with newer versions it is
        recommended to stick with Numpy 1.7 and SciPy 1.12 which was tested
        successfully.
    - `matplotlib <http://matplotlib.org/downloads.html>`_

In addition to the dependencies of the MCSAS package listed above the
`cx_Freeze package <http://cx-freeze.readthedocs.org/en/latest/>`_
is used for freezing the python source code structure into a standalone
package.

Working with Source Code Repositories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to download the latest source code repositories of packages such as
MCSAS or cx_Freeze a client to `Git`_ and `Mercurial`_ is required. There are
several available, for both Mac OS X and Windows the `SourceTree`_ program
is recommended.

.. _Git: https://en.wikipedia.org/wiki/Git_%28software%29
.. _Mercurial: https://en.wikipedia.org/wiki/Mercurial
.. _SourceTree: http://www.sourcetreeapp.com

Windows
-------
A self-contained archive consisting of ``MCSAS.exe`` and all necessary
libraries and files is created by the following command executed in the
MCSAS folder::

    > python cxfreeze.py build_exe

Requirements
^^^^^^^^^^^^
On a fresh installation of Windows 7 the following packages are required:

    - `Python 2.7.9 <https://www.python.org/ftp/python/2.7.9/python-2.7.9.msi>`_

    - `PySide 1.2.1 <https://download.qt.io/official_releases/pyside/PySide-1.2.1.win32-py2.7.exe>`_

    - `NumPy 1.7.1 <http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/numpy-1.7.1-win32-superpack-python2.7.exe>`_

    - `SciPy 0.12.0 <http://sourceforge.net/projects/scipy/files/scipy/0.12.0/scipy-0.12.0-win32-superpack-python2.7.exe>`_

    - `matplotlib 1.4.2 <https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.4.2/windows/matplotlib-1.4.2.win32-py2.7.exe>`_ and its requirements:

        - `Six 1.9.0 <https://pypi.python.org/packages/3.3/s/six/six-1.9.0-py2.py3-none-any.whl>`_
            Install it on the command line by::

                pip install six-1.9.0-py2.py3-none-any.whl

        - `dateutil 2.4.0 <https://pypi.python.org/packages/py2.py3/p/python-dateutil/python_dateutil-2.4.0-py2.py3-none-any.whl>`_

        - `pyparsing 2.0.3 <http://sourceforge.net/projects/pyparsing/files/pyparsing/pyparsing-2.0.3/pyparsing-2.0.3.win32-py2.7.exe>`_

    - `cx_Freeze 4.3.4 <https://pypi.python.org/packages/2.7/c/cx_Freeze/cx_Freeze-4.3.4.win32-py2.7.exe>`_

    - `pywin32 219 <http://sourceforge.net/projects/pywin32/files/pywin32/Build%20219/pywin32-219.win32-py2.7.exe>`_

    - h5py HDF5 support, install one of the precompiled Windows packages, such as
      `h5py-2.4.0.win32-py2.7.exe <https://pypi.python.org/packages/df/6d/b9463b64fa8ff6bd71ae7c5a5b3591dbc1c8536af0926f0bfb712a4356fe/h5py-2.4.0.win32-py2.7.exe#md5=032d7d10a9938aefbcc278e2fe5e2864>`_

Mac OS X
--------
After installing the required packages below a disk image file (.dmg)
consisting of the application bundle is created by::

    $ /usr/local/bin/python2 cxfreeze.py bdist_dmg

Alternatively, for testing purposes the bundle can be created without
packaging into a disk image by::

    $ /usr/local/bin/python2 cxfreeze.py bdist_mac

Requirements
^^^^^^^^^^^^
On a fresh installation of OS X 10.8 the following packages are required:

    - Xcode command line tools: for build essentials such as a compiler
        ( `xcode461_cltools_10_86938245a.dmg <https://developer.apple.com/downloads/download.action?path=Developer_Tools/command_line_tools_os_x_mountain_lion_for_xcode__march_2013/xcode461_cltools_10_86938245a.dmg>`_ )

    - `Python 2.7.9 <https://www.python.org/ftp/python/2.7.9/python-2.7.9-macosx10.6.pkg>`_

    - `Qt 4.8.6 <http://download.qt.io/official_releases/qt/4.8/4.8.6/qt-opensource-mac-4.8.6-1.dmg>`_

    - `PySide 1.2.1 / Qt 4.8 <http://pyside.markus-ullmann.de/pyside-1.2.1-qt4.8.5-py27apple-developer-signed.pkg>`_

    - `NumPy 1.7.1 <http://downloads.sourceforge.net/project/numpy/NumPy/1.7.1/numpy-1.7.1-py2.7-python.org-macosx10.6.dmg>`_

    - `SciPy 0.12.0 <http://sourceforge.net/projects/scipy/files/scipy/0.12.0/scipy-0.12.0-py2.7-python.org-macosx10.6.dmg>`_

    - `matplotlib 1.4.2 <https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.4.2/mac/matplotlib-1.4.2-cp27-none-macosx_10_6_intel.macosx_10_9_intel.macosx_10_9_x86_64.whl>`_
        Install it on the command line by::

                $ /usr/local/bin/pip install matplotlib-1.4.2-*.whl

    - h5py HDF5 support, install HDF5 from source first::

            $ cd hdf5-src
            $ ./configure --prefix=/usr/local
            $ make && sudo make install
            $
            $ pip2 install h5py

    - `a modified cx_Freeze 4.3.4 <https://bitbucket.org/ibressler/cx_freeze>`_
      with local modifications for successful app freezing on OS X

        Download the source and install it on the command line by::

                $ hg clone https://bitbucket.org/ibressler/cx_freeze
                $ cd cx_freeze
                $ hg co 4.x
                $ /usr/local/bin/python2 setup.py install

Ubuntu/Linux
------------
Similar to the procedure on Windows a self-contained archive containing all
necessary libraries and files is created by::

    $ python cxfreeze.py build_exe

Requirements
^^^^^^^^^^^^
On a fresh installation of Ubuntu Linux 14.04 LTS the following packages
are required:

    - ``apt-get install git build-essential python-setuptools python-dev liblapack3-dev libfreetype6-dev tk-dev``

    - PySide 1.2.1

    - NumPy 1.7.1

    - SciPy 0.12.0

    - matplotlib 1.4.2

    - cx_Freeze 4.3.4

Internals
=========
"""

import sys
import gui.version

OLDVERSIONNUM = gui.version.version.number()
if len(sys.argv) > 2 and __name__ == "__main__":
    alternateVersion = "-".join(sys.argv[2:])
    del sys.argv[2:]
    print "Using an alternate version: '{0}'".format(alternateVersion)
    gui.version.version.updateFile(gui.version, alternateVersion)
    reload(gui.version)

import re
import os
import subprocess
import hashlib
import platform
import tempfile
from cx_Freeze import setup, Executable
from gui.version import version
from utils import isWindows, isLinux, isMac, testfor, mcopen
from utils.findmodels import FindModels

def sanitizeVersionNumber(number):
    """Removes non-digits to be compatible with pywin32"""
    if not isWindows():
        return number
    match = re.match(r'(\d+)[^\d](\d+)[^\d](\d+)', number)
    if match is None:
        return "0"
    number = ".".join(match.groups())
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
        testfor(self._path is not None and os.path.isfile(self._path),
                OSError,
                "{name}: '{path}' not found!".format(
                    name = self._name,
                    path = self._path))

    @property
    def execName(self):
        return os.path.basename(self._path)

    def getLogFilename(self):
        return os.path.splitext(self.execName)[0] + ".log" # 7z.log

    def archive(self, targetPath):
        """Creates an archive from the given absolute target directory path.
        The archive file will have the base name of the last directory of the
        given path.
        """
        raise NotImplementedError

class Archiver7z(Archiver):
    _name = "7-Zip"
    _types = ("7z", "zip")
    _ext = None
    _type = None

    def __init__(self, filetype = "7z"):
        super(Archiver7z, self).__init__()
        if filetype not in self._types:
            filetype = self._types[0]
        self._type = filetype
        self._ext = filetype

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
        targetPath = os.path.abspath(targetPath)
        if not os.path.isdir(targetPath):
            return None
        fnPackage = targetPath + "." + self._ext
        fnLog = os.path.abspath(self.getLogFilename())
        with mcopen(fnLog, 'w') as fd:
            retcode = subprocess.call(
                [self._path, "a", "-t" + self._type, "-mx=9",
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
        targetPath = os.path.abspath(targetPath)
        if not os.path.isdir(targetPath):
            return None
        fnPackage = targetPath + ".zip"
        fnLog = os.path.abspath(self.getLogFilename())
        with mcopen(fnLog, 'w') as fd:
            retcode = subprocess.call([self._path, "-r9",
                                       fnPackage, targetPath],
                                       stdout = fd,
                                       stderr = fd)
        if not os.path.exists(fnPackage):
            fnPackage = None
        return fnPackage

def includeModels(includeFilesLst):
    modelsDir = FindModels.rootName
    for path, relFile in FindModels.candidateFiles():
        src = os.path.join(path, relFile)
        dst = os.path.join(modelsDir, relFile)
        includeFilesLst.append((src, dst))

if __name__ == "__main__":

    # using zip by default, its preinstalled everywhere
    archiver = None
    # TODO: perhaps switch to pythons builtin zip?
    if isWindows():
        archiver = Archiver7z(filetype = "zip")
    else:
        archiver = ArchiverZip()

    PACKAGENAME = "{name} {ver}".format(
                    name = version.name(),
                    ver = version.number())
    # target (temp) dir for mcsas package
    TARGETDIR = "{pckg} for {plat}".format(
                    pckg = PACKAGENAME,
                    plat = platform.system().title())

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
            ("resources/icon/mcsas.ico", "resources/icon/mcsas.ico"),
            "dejavuserif.ttf", "dejavumono.ttf",
    ]
    # python3 compatibility fix
    import lib2to3
    lib23_path = os.path.dirname(lib2to3.__file__)
    INCLUDEFILES.append(lib23_path)
    # testing
    includeModels(INCLUDEFILES)

    if isLinux():
        INCLUDEFILES += [
            "/usr/lib/liblapack.so.3",
            "/usr/lib/libblas.so.3",
        ]
        for lib in ("libgfortran.so.3", "libquadmath.so.0",
            "libpyside-python2.7.so.1.2", "libshiboken-python2.7.so.1.2",
            "libQtGui.so.4", "libQtCore.so.4", "libQtSvg.so.4", "libQtXml.so.4",
            "libaudio.so.2"):
            INCLUDEFILES.append(os.path.join("/usr/lib/x86_64-linux-gnu", lib))
    if isWindows():
        INCLUDEFILES += [
            "Microsoft.VC90.CRT",
        ]
    if isMac():
        INCLUDEFILES += [
            ("resources/icon/mcsas.icns", "resources/icon/mcsas.icns"),
            "/usr/lib/system/libdnsinfo.dylib"
        ]

    BUILDOPTIONS = dict(
        compressed = False,
        include_files = INCLUDEFILES,
        packages = [],
        excludes = ["lib2to3", # python3 compatibility fix
                    "models.sphere"], # source file added for dyn. loading
        includes = ["PySide", "PySide.QtCore", "PySide.QtGui",
                    "PySide.QtSvg", "PySide.QtXml",
                    "multiprocessing",
                    "scipy.sparse.csgraph._validation",
                    "scipy.special._ufuncs_cxx",
                    "scipy.sparse.linalg.dsolve.umfpack", "scipy.integrate.vode",
                    "scipy.integrate.lsoda", "matplotlib",
                    "matplotlib.backends.backend_qt4agg",
                    # savefig dependencies?
                    "matplotlib.backends.backend_tkagg", "Tkinter", "FileDialog",
                    "h5py",
                    "UserList", "UserString",
                    ],
        # Icons for Windows using this approach:
        # https://stackoverflow.com/a/10819673
        icon = "resources/icon/mcsas.ico",
        path = [os.getcwd()] + sys.path,
        build_exe = TARGETDIR,
        silent = False,
        copy_dependent_files = True,
    )

    # OSX bundle building
    MACOPTIONS = dict(
        bundle_name = PACKAGENAME,
        # creating icns for OSX: https://stackoverflow.com/a/20703594
        iconfile = "resources/icon/mcsas.icns"
    )
    DMGOPTIONS = dict(volume_label = PACKAGENAME)
    if isMac():
        BUILDOPTIONS.pop("build_exe") # bdist_mac expects plain 'build' directory
        BUILDOPTIONS["includes"] = [ # order in which they were requested
            "PySide", "PySide.QtCore", "PySide.QtGui", "PySide.QtSvg", "PySide.QtXml",
            "scipy.sparse.csgraph._validation",
            "scipy.sparse.linalg.dsolve.umfpack",
            "matplotlib.backends.backend_qt4agg",
            "scipy.integrate.vode",
            "scipy.integrate.lsoda",
            "h5py", "UserList", "UserString",
        ]
        BUILDOPTIONS.pop("icon")
        # tcl/tk is installed by default
        BUILDOPTIONS["bin_excludes"] = ["Tcl", "Tk", "libgcc_s.1.dylib"]
        BUILDOPTIONS["excludes"] += ["Tkinter"]
        os.environ["DYLD_FRAMEWORK_PATH"] = ":".join((
                "/Library/Frameworks", "/System/Library/Frameworks"))
        os.environ["DYLD_LIBRARY_PATH"] = ":".join((
                "/usr/lib", "/usr/local/lib"))

    setup(
        name = version.name(),
        version = sanitizeVersionNumber(version.number()),
        description = version.name(),
        long_description = ("GUI for Monte-Carlo size distribution analysis"),
        license = "Creative Commons CC-BY-SA",
        author = u"Brian R. Pauw",
        author_email = "brian@stack.nl",
        contact = u"Brian R. Pauw",
        contact_email = "brian@stack.nl",
        maintainer = u"Ingo Breßler",
        maintainer_email = "ingo.bressler@bam.de",
        # additional metadata for modified version of cx_Freeze
        # displayed in copyright info
        download_url = "2015, http://www.mcsas.net",
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

    # package the freezed program into an archive file
    PACKAGEFN = archiver.archive(TARGETDIR)
    if PACKAGEFN is None or not os.path.isfile(PACKAGEFN):
        sys.exit(0) # nothing to do

    # calc a checksum of the package
    def hashFile(filename, hasher, blocksize = 65536):
        os.path.abspath(filename)
        with mcopen(filename, 'rb') as fd:
            buf = fd.read(blocksize)
            while len(buf) > 0:
                hasher.update(buf)
                buf = fd.read(blocksize)
        return hasher.hexdigest(), os.path.basename(filename)
    hashValue = hashFile(PACKAGEFN, hashlib.sha256())

    # write the checksum to file
    with mcopen(os.path.abspath(version.name() + ".sha"), 'w') as fd:
        fd.write(" *".join(hashValue))

    # restore initially modified version
    if version.number() != OLDVERSIONNUM:
        version.updateFile(gui.version, OLDVERSIONNUM)

# vim: set ts=4 sw=4 sts=4 tw=0:
