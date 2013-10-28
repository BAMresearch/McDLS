# -*- coding: utf-8 -*-
# contextmenuwidget.py

from collections import OrderedDict
from cutesnake.utils import isString, isList, isMap
from cutesnake.utils.translate import tr
from cutesnake.widgets.mixins.titlehandler import TitleHandler
from cutesnake.qt import QtCore, QtGui
from QtCore import Qt
from QtGui import QAction

def escapeAmp(text):
    return text.replace("&", "&&")

class ContextMenuWidget(object):
    """
    Menu states are mutual exclusive states of available context menu actions.

    Menu states map to boolean methods of this object. The names of a state
    equal names of methods of this object otherwise the menu state is useless.
    updateMenu() tests the states, i.e. boolean methods for a true return
    value. The first one specifies the menu to be build. If none matches, the
    'default' state will be used.

    Menu states are defined by updateMenuState(). The provided list of entry
    names has to be a partial set of those returned by menuEntries(). Entries
    may be defined by addMenuEntry().

    TODO: Instead of multiple inheritance, convert it to a standalone object.
    Problems: need to store bound methods as callbacks, to get state on update()
    The widget this menu applies to needs to be stored only once at init.
    """
    allMenuStates = "*"
    _actions = None
    _states = None

    def __init__(self):
        self.setContextMenuPolicy(Qt.ActionsContextMenu)
        self._actions = dict()
        self._states = OrderedDict() # preserves order

    def _updateTitle(self, text):
        """Updates a menu title with the widget title."""
        text = unicode(text)
        if '%' in text:
            title = self.title
            if isinstance(title, TitleHandler):
                title = title()
            text = text.replace("%1", "{0}")
            return text.format(title.lower())
        return text

    def addMenuEntry(self, name = None, text = None, toolTip = None,
                     shortCut = None, checkable = False, checked = False,
                     callbacks = None, menuStates = None):
        """
        Argument 'callbacks' is supposed to be a list of methods
        which will be connected to the QAction.triggered signal.
        """
        a = QAction(self)
        assert isString(name)
        a.setObjectName("action{0}".format(name.title()))
        assert text is not None
        a.setText(self._updateTitle(escapeAmp(text)))
        if toolTip is not None:
            a.setToolTip(self._updateTitle(toolTip))
        if shortCut is not None:
            a.setShortcut(shortCut)
            a.setShortcutContext(Qt.WidgetShortcut)
        a.setCheckable(bool(checkable))
        a.setChecked(bool(checked))
        if not isList(callbacks):
            callbacks = (callbacks, )
        for callback in callbacks:
            if callback is None:
                continue
            assert callable(callback), tr("Specified callbacks have to be "
                                          "a method or a function!")
            a.triggered.connect(callback)
        self.addMenuEntryAction(name, a, menuStates)
        return a

    def addMenuEntryAction(self, name, action, menuStates = None):
        """
        Common menu states: default and all. 'all' includes the default state.
        Everything else expects a method of the same name which evaluates to
        True or False. It is used to show or hide the respective context menu
        action.
        Actions bound to the same state appear grouped in the menu. The
        underlying ordered dict preserves the order of states.
        """
        action.setParent(self)
        action.setEnabled(True)
        if name in self._actions:
            name = "{0}_{1}".format(name, hash(action))
        self._actions[name] = action
        if menuStates is None:
            menuStates = self.allMenuStates
        if not isList(menuStates):
            menuStates = (menuStates,)
        for stateName in menuStates:
            self.addToMenuState(stateName, name)

    def addToMenuState(self, stateName, *entryNames):
        assert isString(stateName)
        if stateName != self.allMenuStates: # default state
            assert hasattr(self, stateName)
        assert isList(entryNames)
        entryList = self._states.get(stateName, [])
        for name in entryNames:
            if name in self._actions:
                entryList.append(name)
        self._states[stateName] = entryList

    def addMenuSeparator(self, menuStates = None):
        a = QAction(self)
        a.setSeparator(True)
        self.addMenuEntryAction("sep", a, menuStates)

    def removeMenuEntries(self, stateName):
        for name in self.menuEntries(stateName):
            del self._actions[name]
        try:
            del self._states[stateName]
        except KeyError:
            pass
        self.updateMenu()

    def menuEntries(self, stateName):
        return self._states.get(stateName, [])[:] # copy!

    def updateMenu(self, widget = None):
        if widget is None:
            widget = self
        # remove all actions first
        for action in widget.actions():
            widget.removeAction(action)
        # default actions
        availableActions = self.menuEntries(self.allMenuStates) # copy!
        # try to find a set of actions depending on the current state
        for attrname, actions in self._states.iteritems():
            if (not hasattr(self, attrname) or
                not getattr(self, attrname)()):
                continue
            availableActions += actions
        # set actions for current state
        prevAction = None
        for actionName in availableActions:
            action = self._actions.get(actionName)
            if (isinstance(prevAction, QAction) and
                prevAction.isSeparator() and
                action.isSeparator()):
                # skip adjacent separators
                continue
            assert isinstance(action, QAction)
            widget.addAction(action)

# vim: set ts=4 sts=4 sw=4 tw=0:
