# -*- coding: utf-8 -*-
# gui/utils/signal.py

from ..qt import QtCore

if not hasattr(QtCore, 'Signal'):  ## for pyside compatibility
    QtCore.Signal = QtCore.pyqtSignal
    QtCore.Slot = QtCore.pyqtSlot

Signal = QtCore.Signal

def tryDisconnect(sig, slot):
    """Tries to disconnect signal *sig* from slot *slot*."""
    try:
        sig.disconnect(slot)
    except RuntimeError: # Failed to disconnect signal
        pass

# vim: set ts=4 sts=4 sw=4 tw=0:
