.. Find the reST syntax at http://sphinx-doc.org/rest.html

*********************************
How to generate the documentation
*********************************

Requirements
============

* `Python <http://python.org/>`_, of course
* `Sphinx <http://sphinx-doc.org/>`_ package
* For Latex/PDF generation, there should be a latex environment installed

Generate a PDF document
=======================
::

    cd <mcsas>/doc
    make latexpdf

The resulting ``McSAS.pdf`` can be found in ``<mcsas>/doc/_build/latex/``.

Generate HTML pages
===================
::

    cd <mcsas>/doc
    make html

The entry point ``index.html`` can be found in ``<mcsas>/doc/_build/html/``.

Update Source Code Documentation
================================
::
    $ sphinx-apidoc --separate --force --doc-project=MCSAS --output-dir=doc/source/code .

This command automatically generates sphinx documentation files for all
source code files in the directory. It assumes the current working
directory is the MCSAS root directory.

.. vim: set ts=4 sts=4 sw=4 tw=0:
