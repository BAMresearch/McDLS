# -*- coding: utf-8 -*-
# classproperty.py

class classproperty(property):
    """
    Subclass property to make classmethod properties possible.

    Use it like this:

    @classproperty
    @classmethod
    def var(cls):
        <code>

    Getters only, see
    http://stackoverflow.com/questions/128573/using-property-on-classmethods
    """
    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

# vim: set ts=4 sts=4 sw=4 tw=0:
