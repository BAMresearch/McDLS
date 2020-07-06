# -*- coding: utf-8 -*-
# gui/utils/progressdialog.py

from gui.utils.translate import tr
from gui.utils import processEventLoop
from QtCore import Qt
from QtWidgets import QProgressDialog, QApplication

class ProgressDialog(QProgressDialog):
    """
    A progress dialog to visualize calculation status of the program.

    Tests the ProgressDialog with 10 steps, cancels in between.
    >>> import time
    >>> from utilsgui import ProgressDialog, DialogInteraction
    >>> pd = DialogInteraction.queryInstance(ProgressDialog,
    ...                                      slot = 'cancel',
    ...                                      count = 11)
    >>> for i in range(0, pd.maximum()):
    ...      dummy = pd.update()
    ...      time.sleep(0.2)
    >>> pd.wasCanceled()
    True

    Tests again, without cancel this time.
    >>> pd.reset()
    >>> for i in range(0, pd.maximum()):
    ...      dummy = pd.update()
    ...      time.sleep(0.2)
    >>> pd.wasCanceled()
    False
    """

    def __init__(self, parent = None,
                 title = tr("Please wait ..."), count = 0):
        QProgressDialog.__init__(self, parent)
        self.setWindowTitle(title)
        self.setWindowModality(Qt.ApplicationModal)
        self.setRange(0, count)
        self.setValue(0)
        self.setMinimumDuration(300)

    def update(self):
        """Updates progress status and returns True if canceled,
        False otherwise."""
        self.setValue(self.value()+1)
        processEventLoop() # process canceled signal eventually
        return self.wasCanceled()

# vim: set ts=4 sts=4 sw=4 tw=0:
