# -*- coding: utf-8 -*-
# __init__.py

import sys
from cutesnake.qt import QtGui
from QtGui import QApplication
from cutesnake.log import Log
from mainwindow import MainWindow

def eventLoop():
    app = QApplication(sys.argv)
    mw = MainWindow()
    # run initialization stuff
    Log.replaceStdOutErr()
    if not Log.hasHandler():
        Log.setConsoleHandler()
    mw.show()
    return app.exec_()

# vim: set ts=4 sts=4 sw=4 tw=0:
