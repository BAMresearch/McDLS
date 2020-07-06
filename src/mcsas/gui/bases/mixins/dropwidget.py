# -*- coding: utf-8 -*-
# gui/bases/mixins/dropwidget.py

class DropWidget(object):
    """
    Drag&Drop support for widgets which inherit from this.
    """
    def dragEnterEvent(self, ev):
        if ev.mimeData().hasFormat('text/uri-list'):
            ev.acceptProposedAction()

    def dropEvent(self, ev):
        lst = [str(url.toLocalFile()) for url in ev.mimeData().urls()]
        self.sigReceivedUrls.emit(lst)
        ev.acceptProposedAction()

# vim: set ts=4 sts=4 sw=4 tw=0: 
