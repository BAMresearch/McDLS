# -*- coding: utf-8 -*-
# gui/utils/dialoginteraction.py

from __future__ import print_function
from QtCore import (Qt, QThread, QMetaObject)
from QtWidgets import (QApplication, QMessageBox)

class DialogInteraction(QThread):
    """Simulates user interaction on blocking widgets.

    Calls a specified slot on all widgets of a given type after a
    given delay. Please see tests in other classes for examples.
    """
    _application = None
    _properties = None
    _propertyNames = ('instance', 'queryType', 'slot', 'button', 'delay')
    _defaultDelay = 1.0
    _success = None

    @classmethod
    def application(cls):
        if cls._application is None:
            cls._application = QApplication([])
        return cls._application

    @classmethod
    def instance(cls, instanceType, *args, **kwargs):
        if instanceType is None:
            return None
        cls.application()
        instance = instanceType(*args, **kwargs)
        return instance

    @classmethod
    def query(cls, queryType, blockingFunc, *args, **kwargs):
        if queryType is None:
            return
        slot = kwargs.pop('slot', None)
        button = kwargs.pop('button', None)
        delay = kwargs.pop('delay', cls._defaultDelay)
        if blockingFunc is None or (slot is None and button is None):
            return
        di = cls({'queryType': queryType,
                  'slot':      slot,
                  'button':    button,
                  'delay':     delay})
        result = None
        di.start()
        try:
            result = blockingFunc(*args, **kwargs)
        except Exception:
            raise
        finally:
            di.quit()
            di.wait()
        return result

    @classmethod
    def queryInstance(cls, instanceType, *args, **kwargs):
        instance = cls.query(instanceType,
                             cls.instance, instanceType,
                             *args, **kwargs)
        return instance

    def __init__(self, config):
        """
        *config*: A *dict* containing argument/value pairs. For valid
        arguments, please see _propertyNames.
        """
        parent = config.pop('parent', None)
        QThread.__init__(self, parent)
        self._properties = dict(config)
        for key in self._propertyNames:
            if key not in self._properties:
                self._properties[key] = None
        # fix values
        delay = self._get('delay', self._defaultDelay)
        self._properties['delay'] = int(delay*1e6)
        if self._get('queryType') is None:
            self._properties['queryType'] = self._get('instance')
        self._success = False # check success of query

    def _get(self, key, default = None):
        return self._properties.get(key, default)

    def _callSlot(self, widget):
        if self._get('slot') is None:
            return
        QMetaObject.invokeMethod(widget, self._get('slot'),
                                 Qt.QueuedConnection)
        self._success = True

    def _clickButton(self, widget):
        button = None
        if self._get('button') == 'escape':
            button = widget.escapeButton()
        elif self._get('button') == 'accept':
            for btn in widget.buttons():
                if widget.buttonRole(btn) == QMessageBox.AcceptRole:
                    button = btn
                    break
        elif self._get('button') == 'default':
            button = widget.defaultButton()
        # call click() on selected button
        QMetaObject.invokeMethod(button, 'click',
                                 Qt.QueuedConnection)
        self._success = True

    def _matchedWidgets(self, queryType):
        for widget in QApplication.allWidgets():
            if type(widget) is queryType:
                yield widget

    def run(self):
        if self._get('queryType') is None:
            return
        self.usleep(self._get('delay'))
        for widget in self._matchedWidgets(self._get('queryType')):
            if self._get('button') is not None:
                self._clickButton(widget)
            elif self._get('slot') is not None:
                self._callSlot(widget)
        if not self._success:
            print ("Could not find widget of type {0}"
                   .format(self._get('queryType')))

# vim: set sts=4 ts=4 sw=4 tw=0:
