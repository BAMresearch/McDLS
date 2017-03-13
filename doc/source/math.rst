.. Find the reST syntax at http://sphinx-doc.org/rest.html
.. TODO: download mathjax in conf.py and put it into the static dir
.. http://www.sphinx-doc.org/en/stable/ext/math.html

*********************************
The Math Behind
*********************************

====
SAXS
====

In :class:`models.sphere` the form factor is defined as:

.. math::

    \def\v#1{{\bf #1}}
    ff_{sph}(q, r) &= { 3 ~ sin(qr) - qr ~ cos(qr) \over (qr)^3 } \\
    v_{sph}(r) &= {\tiny {4 \over 3}} ~ \pi ~ r^3 \\
    v_{sph,abs}(r, \Delta\rho) &= \Delta\rho^2 ~ v_{sph}(r)

Where *q* is
    the scattering vector loaded from the data file and possibly preprocessed,
    respectively filtered by defining min/max *q* or masking invalid values
    equal or below zero.

*r* denotes
    the radius of the sphere set in the user interface (UI) or varied during
    optimization.

:math:`\Delta\rho` denotes
    the scattering length density difference constant of the model against the
    solution which is defined in the UI.

.. related questions: http://stackoverflow.com/q/7825263

.. automethod:: models.sphere.Sphere.formfactor
    :noindex:

.. automethod:: models.sphere.Sphere.volume
    :noindex:

.. automethod:: models.sphere.Sphere.absVolume
    :noindex:

.. automethod:: bases.model.SASModel.weight
    :noindex:

.. automethod:: bases.model.SASModel.calcIntensity
    :noindex:

.. vim: set ts=4 sts=4 sw=4 tw=0:
