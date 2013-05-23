# -*-coding: utf8-*-
# mixedmethod.py

"""This module implements a mixedmethod() decorator for class definitions.
Basic idea found here:
http://www.daniweb.com/software-development/python/code/406393/mixedmethod-mixes-class-and-instance-method-into-one
"""

from functools import partial

class mixedmethod(object):
    """This decorator mutates a function defined in a class into a 'mixed'
    class and instance method. Usage::
    
        class Spam:
            @mixedmethod
            def egg(selforcls, *args, **kwargs):
                # selforcls is the instance when called on an instance
                # selforcls is the class when called on the class
                pass
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, instance, cls):
        if instance is None:
            return partial(self.func, cls)
        else:
            return partial(self.func, instance)

if __name__ == '__main__':
    
    class Spam(object):
        bla = None
        
        @mixedmethod
        def ham(selforcls, val):
            selforcls.bla = val

    egg = Spam()
    egg.ham(5)
    Spam.ham(4)

    print "Spam: ", Spam.bla, "egg:", egg.bla

# vim: set ts=4 sts=4 sw=4 tw=0:
