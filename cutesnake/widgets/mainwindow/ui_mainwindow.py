# -*- coding: utf-8 -*-

from cutesnake.qt import QtCore, QtGui
from cutesnake.utils.translate import fromUtf8

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(fromUtf8("MainWindow"))
        MainWindow.resize(800, 600)
        MainWindow.setStyleSheet(fromUtf8("\n"
"   QMainWindow::separator:hover {\n"
"    background: darkgrey;\n"
"   }\n"
"  "))
        MainWindow.setDocumentMode(True)
        MainWindow.setDockOptions(QtGui.QMainWindow.AnimatedDocks|QtGui.QMainWindow.ForceTabbedDocks)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(fromUtf8("centralwidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setObjectName(fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        pass

import mainwindow_rc
